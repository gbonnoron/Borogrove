from math import log, sqrt, erfc, ceil, floor

# Loading security estimator script for LWE

load("lwe-estimator/estimator.py")


# ##############    Utility functions    #################

def safelog2(x):
    if x > 0:
        return log(x)/log(2)
    else:
        return -float("inf")


def sum_square(x):
    # 1 + 2^2 + 3^2 + ... + x^2
    t = x - 1
    return t * (t+1) * (2*t + 1) / 6


def variance_uniform_modB(B, balanced=True):
    if balanced:
        # Variance of a random variable mod B, represented in [1-B/2 , ... B-1/2]
        b1 = int(floor((B-1.) / 2))
        b2 = int(ceil((B-1.) / 2))
        v = (1. / B) * (sum_square(b1) + sum_square(b2))
    else:
        # Variance of a random variable mod B, represented in [0, ... B-1]
        v = (1. / B) * (sum_square(B-1))
    return v


def GaussianOutsideRangeProba(sigma, t):
    return erfc(t / (sqrt(2) * sigma))


def DecryptionErrorProbability(c):
    return GaussianOutsideRangeProba(c.err_stddev, c.Q/(2.*c.t))


# ##############    Data types    #################


class CLWE_CT:
    def __init__(self, d, n, Q, t, err, secret_density=None):
        self.d = d
        self.n = n
        self.Q = Q
        self.t = t
        self.err_stddev = err
        self.secret_density = secret_density

    def __str__(self):
        res = "CLWE Ciphertext:\t d:%d, n:%d, \tQ:2^%.2f, \tt:%d, \te:2^%.2f" % (
            self.d, self.n, log(self.Q)/log(2), self.t, log(self.err_stddev)/log(2))
        if self.secret_density is not None:
            res += "\tsd:%.3f"%self.secret_density
        return res

    def security(self):
        raise ValueError("In the current scheme, CLWE_CT is always \
                          generated from something else, so you shouldn't \
                          need to estimate their hardness.")
        # dim = self.n * self.d
        # h = int(floor(self.secret_density * dim))
        # alpha = RR(sqrt(2*pi) * self.err_stddev / self.Q)
        # return estimate_lwe(dim, alpha, self.Q, secret_distribution=secret_distribution)


class CGSW_CT:
    def __init__(self, d, n, Q, t, B, err, secret_density=None):
        self.d = d
        self.n = n
        self.Q = Q
        self.t = t
        self.err_stddev = err
        self.B = B
        self.K = int(ceil(log(Q)/log(B)))
        self.secret_density = secret_density

    def __str__(self):
        res = "CGSW Ciphertext:\t d:%d, n:%d, \tQ:2^%.2f, \tt:%d, \tB:2^%.2f, \tK:%d, \te:2^%.2f" % (
            self.d, self.n, log(self.Q)/log(2), self.t, log(self.B)/log(2), self.K, log(self.err_stddev)/log(2))
        if self.secret_density is not None:
            res += "\tsd:%.3f"%self.secret_density
        return res

    def security(self):
        dim = self.n * self.d
        h = int(floor(self.secret_density * dim))
        alpha = RR(sqrt(2*pi) * self.err_stddev / self.Q)
        return estimate_lwe(dim, alpha, self.Q, secret_distribution=((-1,1),h))


class KS_Key:
    def __init__(self, d, n, nn, Q, B, err, secret_density=None):
        self.d = d
        self.n = n
        self.nn = nn
        self.Q = Q
        self.err_stddev = err
        self.B = B
        self.K = int(ceil(log(Q)/log(B)))
        self.secret_density = secret_density

    def __str__(self):
        res = "KeySwitching Key:\t d:%d, n->nn:%d->%d, \tQ:2^%.2f, \tB:2^%.2f, \tK:%d, \te:2^%.2f" % (
            self.d, self.n, self.nn, log(self.Q)/log(2), log(self.B)/log(2), self.K, log(self.err_stddev)/log(2))
        if self.secret_density is not None:
            res += "\tsd:%.3f"%self.secret_density
        return res

    def security(self):
        dim = self.nn * self.d

        if self.secret_density is None:
            secret_distribution=True
        else:
            h = int(floor(self.secret_density * dim))
            secret_distribution=((-1, 1), h)
        
        alpha = RR(sqrt(2*pi) * self.err_stddev / self.Q)
        return estimate_lwe(dim, alpha, self.Q, secret_distribution=secret_distribution)


# class KS_CRT_Key:
#     def __init__(self, d, n, nn, Q, err, secret_density=None):
#         self.d = d
#         self.n = n
#         self.nn = nn
#         #self.Qis = Qis
#         self.Q = Q
#         # for Qi in Qis:
#         #     self.Q *= Qi

#         self.err_stddev = err
#         self.K = len(Qis)

#         self.secret_density = secret_density

#     def __str__(self):
#         res = "KeySwitchingCRT Key:\t d:%d, n->nn:%d->%d, \tQ:2^%.2f, \tK:%d, \te:2^%.2f" % (
#             self.d, self.n, self.nn, log(self.Q)/log(2), self.K, log(self.err_stddev)/log(2))
#         if self.secret_density is not None:
#             res += "\tsd:%.3f"%self.secret_density
#         return res

    def security(self):
        dim = self.nn * self.d
        h = int(floor(self.secret_density * dim))
        alpha = RR(sqrt(2*pi) * self.err_stddev / self.Q)
        return estimate_lwe(dim, alpha, self.Q, secret_distribution=((-1,1),h))


# ##############    Operations    #################
# The function below reflect the bounds presented in section ``7.3 Heuristic error propagation''


def SumPow2(c, k):
    # Implicitly assume all inputs are independent. In our context, there is no reason
    # to duplicate an input bit.
    v = sqrt(sum([2**(2*i) for i in range(k)]))
    cc = CLWE_CT(c.d, c.n, c.Q, c.t, v * c.err_stddev, c.secret_density)
    return cc


def ModSwitch(QQ, c):
    if not isinstance(c, CLWE_CT):
        raise ValueError("Modulus Switching only applies to CLWE Ciphertexts")

    if c.secret_density is None:
        raise ValueError("Density of the secret undefined")

    rounding_variance = ((c.d * c.n * c.secret_density) + 1.) * 1./12

    print "(ModSwitch adding std. dev. %.2f)" % (safelog2(rounding_variance)/2)
    final_error = sqrt((c.err_stddev * QQ * 1. / c.Q)**2 + rounding_variance)
    cc = CLWE_CT(c.d, c.n, QQ, c.t, final_error, secret_density=c.secret_density)
    return cc


def KeySwitch(ks, c):
    if not isinstance(ks, KS_Key):
        raise ValueError("ks is not a Key Switching Key")
    if not isinstance(c, CLWE_CT):
        raise ValueError("c is not a CLWE Ciphertext")
    if not (c.n == ks.n and c.d == ks.d and c.Q == ks.Q):
        raise ValueError("Dimension/Modulus Mismatch in key Switching")

    # this is a average value, different from the worst-case bound given in the paper as Remark 8
    # the d^2 is replaced by d. This is fairplay because a lot of such operation are used and
    # one can argue they will average out.

    added_variance = (ks.err_stddev**2) * variance_uniform_modB(ks.B) * ks.n * ks.d * ks.K
    print "(KS adding std. dev. %.2f)" % (safelog2(added_variance)/2)
    final_error = sqrt(c.err_stddev**2 + added_variance)
    cc = CLWE_CT(c.d, ks.nn, c.Q, c.t, final_error, secret_density = ks.secret_density)
    return cc


def ExtMult(C, c):
    raise ValueError("Should not be used sirectly in the current scheme")
    # if not isinstance(C, CGSW_CT):
    #     raise ValueError("C is not a CGSW Ciphertext")
    # if not isinstance(c, CLWE_CT):
    #     raise ValueError("c is not a CLWE Ciphertext")
    # if not (c.n == C.n and c.d == C.d and c.Q == C.Q and c.t == C.t and c.secret_density == C.secret_density):
    #     raise ValueError("Dimension/Modulus Mismatch in key Switching")

    # added_variance = (C.err_stddev**2) * variance_uniform_modB(C.B) * C.n * C.d * C.K
    # final_error = sqrt(c.err_stddev**2 + added_variance)
    # cc = CLWE_CT(c.d, c.n, c.Q, c.t, final_error, secret_density=c.secret_density)
    # return cc


def Galois(c):
    if not isinstance(c, CLWE_CT):
        raise ValueError("c is not a CLWE Ciphertext")
    return c


def ExtExpInner(C, ks, l):
    if not isinstance(C, CGSW_CT):
        raise ValueError("C is not a CGSW Ciphertext")
    if not isinstance(ks, KS_Key):
        raise ValueError("ks is not a Key Switching Key")
    if not (ks.d == C.d) and (ks.n == C.n) and (ks.Q == C.Q) and (ks.B == C.B) and (ks.secret_density==C.secret_density):
        raise ValueError("ks and C parameters are not matching")
    if not (ks.n == 1):
    	raise ValueError("This is not used in our scheme")

    variance = C.d * C.K * l * variance_uniform_modB(C.B) * (C.err_stddev**2 + 2 * ks.err_stddev**2)
    cc = CLWE_CT(C.d, C.n, C.Q, C.t, sqrt(variance), secret_density=C.secret_density)

    return cc


def ExpCRT(c1, c2):
    if not (c1.n == 1 and c2.n == 1):
        raise ValueError("ExpCRT not implemented for n!=1")
    if not (c1.Q == c2.Q and c1.t == c2.t):
        raise ValueError("Dimension/Modulus Mismatch in key Switching")

    final_error = sqrt((c1.t * c1.err_stddev * c2.err_stddev) ** 2 + c1.err_stddev**2 + c2.err_stddev**2)
    cc = CLWE_CT(c1.d * c2.d, 3, c1.Q, c1.t, final_error, secret_density=c1.secret_density * c2.secret_density)
    return cc




# This is the code for the version optimized with TensorGadget trick.
#
# def FunExpExtract(ks, c, p, q):
#     assert c.d == p*q
#     f = c.d / 2
#     sumV_Qi = sum([variance_uniform_modB(Qi)**2 for Qi in ks.Qis])
#     added_variance = 3 * ks.err_stddev**2 * p * q * sumV_Qi
#     final_error = sqrt(c.err_stddev**2 + added_variance)
#     cc = CLWE_CT(1, p, c.Q, c.t, final_error * sqrt(f), secret_density = ks.secret_density)
#     return cc

def FunExpExtract(ks, c, p, q):
    assert c.d == p*q
    f = c.d / 2
    V_B = variance_uniform_modB(ks.B)
    added_variance = 3 * ks.err_stddev**2 * 2 * ks.K * p * q * V_B
    print "(KS adding std. dev. %.2f)" % (safelog2(added_variance)/2)    
    final_error = sqrt(c.err_stddev**2 + added_variance)
    cc = CLWE_CT(1, p, c.Q, c.t, final_error * sqrt(f), secret_density = ks.secret_density)
    return cc



# ##############    Parameters    #################


p = 1439
q = 1447
l = 600


e1 = 2**(1)
Q1 = 2**56
B1 = 2**8
sd1 = .66

Q2 = 2**38

e3 = 2.**(1)
sd3 = .66
Q3 = 2**56
B3 = 2**8

Q4 = Q3
e4 = 2.**(33)
B4 = 2**4
sd4 = .66

Q5 = p * q

t = 2**6
k = 6

###############    Scheme    #################

print " \n === ExtExpInner"

BSp = CGSW_CT(p, 1, Q1, t, B1, e1, secret_density=sd1)
KSp = KS_Key(p, 1, 1, Q1, B1, e1, secret_density=sd1)
print "BSp :", BSp
print "KSp :", KSp
BSq = CGSW_CT(q, 1, Q1, t, B1, e1, secret_density=sd1)
KSq = KS_Key(q, 1, 1, Q1, B1, e1, secret_density=sd1)
print "BSq :", BSq
print "KSq :", KSq

print ""
cp = ExtExpInner(BSp, KSp, l)
cq = ExtExpInner(BSq, KSq, l)

print "cq :", cq
print "cp :", cp

print " \n === ModSwitch"
c2p = ModSwitch(Q2, cp)
print "c2p :", c2p
c2q = ModSwitch(Q2, cq)
print "c2q :", c2q

print " \n === ExpCRT"
cpqt = ExpCRT(c2p, c2q)
print "cpqt :", cpqt

print " \n === ModSwitch"
cpq = ModSwitch(Q3, cpqt)
print "cpq :", cpq


print " \n === FunExpExtract"
KS = KS_Key(p, 3*q, 1, Q3, B3, e3, secret_density=sd3)
print "KS :", KS
c = FunExpExtract(KS, cpq, p, q)
print "c :", c

print " \n === SumPow2"
c_sum = SumPow2(c, k)
print "c_sum :", c_sum

print " \n === Final KeySwitch"
KSfinal = KS_Key(1, p, l, Q4, B4, e4, secret_density=sd4)
print "KSl :", KSfinal

c_pre_final = KeySwitch(KSfinal, c=c_sum)
print "c_pre_final :", c_pre_final


print " \n === Final ModSwitch"
c_final = ModSwitch(Q5, c_pre_final)
print "c_final :", c_final





print
print
ep = DecryptionErrorProbability(c_final)
print " Final failure probability %.4e = 2^%.2f" % (ep, safelog2(ep))
# exit(1)


###############    Security    #################

print
print "==================="
print "     SECURITY      "
print "==================="
print
print "BSp :", BSp
#BSp.security()
print "KSp :", KSp
#KSp.security()

print "BSq :", BSq
#BSq.security()
print "KSq :", KSq
print
KSq.security()

print
print
print "KS :", KS
print
KS.security()

print
print
print "KSl :", KSfinal
print
KSfinal.security()



# print deltaLWE(p, Q, err, 2**(-16))
# print deltaLWE(p, QQ, ee, 2**(-16))
# print deltaLWE(l, p*q, le, 2**(-16))
