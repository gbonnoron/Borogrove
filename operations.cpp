/**
  * @file operations.cpp
  * @date Spring 2017
  * @author Guillaume Bonnoron
  * @brief Definition of some elementary operations of the scheme
  */

#include <ctime>
#include <chrono>
#include "operations.h"
#include "circulant_ring.h"
#include "rgsw.h"
#include "ksw_key.h"
#include <cmath>

void print_param()
{
    std::cout << "T=" << T << ", ";
    std::cout << "Q=" << Q << ", ";
    std::cout << "Qp=" << Qp << ", ";
    std::cout << "P1=" << P1 << ", ";
    std::cout << "P2=" << P2 << ", ";
    std::cout << "N=" << N << ", ";
    std::cout << "B_Q=" << B_Q << ", ";
    std::cout << "K_Q=" << K_Q << ", ";
    std::cout << "B_Qp=" << B_Qp << ", ";
    std::cout << "K_Qp=" << K_Qp << std::endl;
}


template<class R, uint64_t n>
R dot_product(const R* const v1, const R* const v2)
{
    R res;
    res = 0;
    R *v2_rw = new R[n];
    for (size_t i = 0 ; i < n; ++i)
    {
        v2_rw[i] = v2[i];
        v2_rw[i].decomp_fft();
    }
    for (size_t i = 0 ; i < n; ++i)
    {
        res += v1[i] * v2_rw[i];
        //std::cout << "res " << res << " " << v1[i] << " * " << v2[i] << std::endl;
    }
    /*
    std::cout << "v1=(";
    for (size_t i = 0 ; i < n-1; ++i)
        std::cout << v1[i] << ",";
    std::cout << v1[n-1] << ")" << std::endl;
    std::cout << "v2=(";
    for (size_t i = 0 ; i < n-1; ++i)
        std::cout << v2[i] << ",";
    std::cout << v2[n-1] << ")" << std::endl;
    std::cout << "res=" << res << std::endl;
    */
    delete[] v2_rw;
    return res;
}

void exp_crt(RLWE<Rp12_crt, 3> &cpq, const RLWE<Rp1_crt, 1> &cp, const RLWE<Rp2_crt, 1> &cq)
{
    //const RLWE<Rp1_crt, 1> c_p = cp.galois(inv_mod_P1[P2 % P1]);
    const RLWE<Rp1_crt, 1> c_p = cp.galois(P2inv_modP1);
    const Rp1_crt a_p = c_p.get_a()[0];
    const Rp1_crt b_p = c_p.get_b();

    //const RLWE<Rp2_crt, 1> c_q = cq.galois(inv_mod_P2[P1 % P2]);
    const RLWE<Rp2_crt, 1> c_q = cq.galois(P1inv_modP2);
    const Rp2_crt a_q = c_q.get_a()[0];
    const Rp2_crt b_q = c_q.get_b();

    Rp12_crt a[3];
    a[0] = QpoverT_inv*a_p.tensor_product<P2, FFT_DIM1, FFT_DIM2>(a_q);
    a[1] = QpoverT_inv*a_p.tensor_product<P2, FFT_DIM1, FFT_DIM2>(b_q);
    a[2] = QpoverT_inv*b_p.tensor_product<P2, FFT_DIM1, FFT_DIM2>(a_q);
    Rp12_crt b = QpoverT_inv*b_p.tensor_product<P2, FFT_DIM1, FFT_DIM2>(b_q);

    cpq = RLWE<Rp12_crt, 3>(a, b);
}

void crt_key(Rp12 *spq, const Rp1 &sp, const Rp2 &sq)

{
    //const Rp1_p s_p = sp.galois(inv_mod_P1[P2 % P1]);
    //const Rp2_p s_q = sq.galois(inv_mod_P2[P1 % P2]);
    const Rp1_p s_p = sp.galois(P2inv_modP1);
    const Rp2_p s_q = sq.galois(P1inv_modP2);

    spq[0] = - s_p.tensor_product<P2, FFT_DIM1, FFT_DIM2>(s_q);
    spq[1] = s_p.tensor_product<P2, FFT_DIM1, FFT_DIM2>(Rp2_p(1));
    spq[2] = Rp1_p(1).tensor_product<P2, FFT_DIM1, FFT_DIM2>(s_q);
}

template<class R, class Zp, uint64_t p>
RLWE1<R> ext_exp_inner(const uint64_t q, const uint64_t t, const size_t l, const Zp *y, const RGSW<R> *C, const KSWKey<R, 1, 1, B_Q, K_Q> *S)
{
    R zeros[1];
    zeros[0] = 0;
    R b((int64_t)(q/t));
    RLWE1<R> result(zeros, b);

    Zp current;
    size_t one[l], other[l];
    size_t count_one = 0, count_other = 0;

    for (size_t i = 0 ; i < l ; ++i)
    {
        if (y[i] == 0)  // nothing to ext_exp_mult_add
            continue;
        if (y[i] == 1)  // no galois+ksw required
            one[count_one++] = i;
        else
            other[count_other++] = i;
    }

    if (count_other > 0)
    {
        for (size_t i = 0 ; i < count_other-1 ; ++i)
        {
            current = y[other[i]] * y[other[i+1]].inv();
            result.ext_mult_inplace(C[other[i]]);
            if (current != 1)
            {
                result.galois_inplace((size_t) current);
                result.key_switch_inplace(S[(size_t) current]);
            }
        }
        current = y[other[count_other - 1]];
        result.ext_mult_inplace(C[other[count_other - 1]]);
        result.galois_inplace((size_t) current);
        result.key_switch_inplace(S[(size_t) current]);
    }
    for (size_t i = 0 ; i < count_one ; ++i)
        result.ext_mult_inplace(C[one[i]]);

    return result;
}

LWE fun_exp_extract(const fftw_complex *f, const Rp12LWE &c_Zm, const KSWKeyRp12 &S)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> ksw_time, mult_time, trace_time;

    //Rp12 f(F);
    //std::cout << "f: " << f << std::endl;

    start = std::chrono::system_clock::now();
    RLWE<Rp12 , 1> c_1 = c_Zm.key_switch(S);
    end = std::chrono::system_clock::now();
    ksw_time = end - start;

    start = std::chrono::system_clock::now();
    c_1.mult(f);
    end = std::chrono::system_clock::now();
    mult_time = end - start;

    start = std::chrono::system_clock::now();
    RLWE<Rp1_p, 1> c_2 = c_1.template trace<Rp1_p, P1, FFT_DIM1>();
    Zqp *a = c_2.get_a()->get_data();
    Rz c[P1 + 1];
    c[0] = a[0];
    for (size_t i = 1 ; i < P1 ; ++i)
        c[P1 - i] = a[i];
    c[P1] = c_2.get_b().get_data()[0];
    end = std::chrono::system_clock::now();
    trace_time = end - start;
    std::cerr << "In N_f -> KeySwitch: " << ksw_time.count() << "s, Mult: " << mult_time.count() << "s, Tr*:" << trace_time.count() << "s" << std::endl;
    return LWE(c);
}

template <class R, class Z, class R_crt, uint64_t p>
RLWE<R_crt, 1> accumulation(const RLWE<CirculantRing<Zp12, 1, 1>, N> &ab, const RGSW<R> BS[N], const KSWKey<R> KS[p])
{
    /// 1. Prepare the -a_p
    Z *ab_p = new Z[N+1];
    ab.mod<Z>(ab_p);
    for (size_t i = 0 ; i < N ; ++i)
        ab_p[i] = -ab_p[i];

    /// 2. prepare T^{b_p}
    Zq *Tb_coefs = new Zq[p];
    std::fill(Tb_coefs, Tb_coefs + p, 0);
    Tb_coefs[(size_t)ab_p[N]] = 1;
    R Tb(Tb_coefs, false);

    /// 3. Perform the Inner product in the exponent and the mod-switch
    RLWE<R, 1> acc_p = ext_exp_inner<R, Z, p>(Q, T, N, ab_p, BS, KS);
    acc_p.mult(Tb);
    delete[] ab_p;

    return acc_p.template mod_switch<R_crt, Qcrt, Q>();
}

Rp12LWE preparation(const LWE &Em, const KSWKeyLWE &S_lwe, const Rp1GSW Xsi[N], const KSWKeyRp1 KSp1[P1], const Rp2GSW Ysi[N], const KSWKeyRp2 KSp2[P2])
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> ksw_time, msw_time, acc1_time, acc2_time, crt_time;

    start = std::chrono::system_clock::now();
    RLWE<Rz, N> Em_short = Em.key_switch(S_lwe);
    end = std::chrono::system_clock::now();
    ksw_time = end - start;

    start = std::chrono::system_clock::now();
    RLWE<CirculantRing<Zp12, 1, 1>, N> ab_short = Em_short.mod_switch<CirculantRing<Zp12, 1, 1>, P1*P2, Qp>();
    end = std::chrono::system_clock::now();
    msw_time = end - start;

    start = std::chrono::system_clock::now();
    RLWE<Rp1_crt, 1> crt_in_p = accumulation<Rp1, Zp1, Rp1_crt, P1>(ab_short, Xsi, KSp1);
    end = std::chrono::system_clock::now();
    acc1_time = end - start;

    start = std::chrono::system_clock::now();
    RLWE<Rp2_crt, 1> crt_in_q = accumulation<Rp2, Zp2, Rp2_crt, P2>(ab_short, Ysi, KSp2);
    end = std::chrono::system_clock::now();
    acc2_time = end - start;

    RLWE<Rp12_crt, 3> result;
    start = std::chrono::system_clock::now();
    exp_crt(result, crt_in_p, crt_in_q);
    end = std::chrono::system_clock::now();
    crt_time = end - start;
    std::cerr << "In L_c -> KeySwitchLWE: " << ksw_time.count() << "s, ModSwitch: " << msw_time.count() << "s, ExtExpInner 1: " << acc1_time.count() << "s, ExtExpInner 2: " << acc2_time.count() << "s, ExpCRT: " << crt_time.count() << "s" << std::endl;

    return result.template mod_switch<Rp12, Qp, Qcrt>();
}

LWE combination(const size_t k, const int64_t *coefs, const LWE *c_i)
{
    LWE result = c_i[0] * coefs[0];
    for (size_t i = 1 ; i < k ; ++i)
        result += c_i[i] * coefs[i];
    return result;
}

LWE gate(const size_t k, const int64_t *coefs, const LWE *c_i, const fftw_complex *f, const KSWKeyLWE &S_lwe, const Rp1GSW Xsi[N], const KSWKeyRp1 KSp1[P1], const Rp2GSW Ysi[N], const KSWKeyRp2 KSp2[P2], const KSWKeyRp12 &S)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> comb_time, prep_time, fun_time;

    start = std::chrono::system_clock::now();
    LWE Em = combination(k, coefs, c_i);
    end = std::chrono::system_clock::now();
    comb_time = end - start;

    start = std::chrono::system_clock::now();
    Rp12LWE c_Zm = preparation(Em, S_lwe, Xsi, KSp1, Ysi, KSp2);
    end = std::chrono::system_clock::now();
    prep_time = end - start;

    start = std::chrono::system_clock::now();
    LWE Fm = fun_exp_extract(f, c_Zm, S);
    end = std::chrono::system_clock::now();
    fun_time = end - start;

    std::cerr << "Combination: " << comb_time.count() << "s, L_c: " << prep_time.count() << "s, N_f: " << fun_time.count() << "s" << std::endl;
    return Fm;
}

template <class R, uint64_t n, int64_t p>
void gen_bootstrapping_keys(const uint64_t q, const uint64_t t, RGSW<R> Tsi[n], const Rz s[n], const R &s_p, const double variance)
{
    Zt *coefs = new Zt[p];
    for (size_t j = 0 ; j < p ; ++j)
        coefs[j] = 0;
    for (size_t i = 0 ; i < n ; ++i)
    {
        const size_t position = (p + (int64_t)(s[i].get_data()[0])) % p;
        coefs[position] = 1;
        CirculantRing<Zt, p> T_si(coefs); // =T^{s_i}
        Tsi[i].encrypt(&s_p, T_si, variance);
        coefs[position] = 0;
    }
    delete[] coefs;
}

template <class R, uint64_t p>
void gen_keyswitching_keys(KSWKey<R, 1, 1, B_Q, K_Q> keys[p], const R &s_p, const double variance)
{
    R s_alpha[1];
    for (size_t i = 1 ; i < p ; ++i)
    {
        s_alpha[0] = s_p.galois(i);
        keys[i].init(s_alpha, &s_p, variance);
    }
}

void gen_funexpextract_key(KSWKeyRp12 *S, const Rp12 s_pq[3], const Rz s[P1], const double variance)
{
    Zqp *coefs = new Zqp[P1*P2];
    for (size_t i = 0 ; i < P1*P2 ; ++i)
        coefs[i] = 0;
    for (size_t i = 0 ; i < P1 ; ++i)
        coefs[i*P2] = s[i].get_data()[0];
    Rp12 s2(coefs, false);
    S->init(s_pq, &s2, variance);
}

template void gen_bootstrapping_keys<Rp1, N, P1>(const uint64_t q, const uint64_t t, Rp1GSW Tsi[N], const Rz s[N], const Rp1 &s_p, const double variance);
template void gen_bootstrapping_keys<Rp2, N, P2>(const uint64_t q, const uint64_t t, Rp2GSW Tsi[N], const Rz s[N], const Rp2 &s_p, const double variance);

template void gen_keyswitching_keys<Rp1, P1>(KSWKeyRp1 keys[P1], const Rp1 &s_p, const double variance);
template void gen_keyswitching_keys<Rp2, P2>(KSWKeyRp2 keys[P2], const Rp2 &s_p, const double variance);

template RLWE1<Rp1> ext_exp_inner<Rp1, Zp1, P1>(const uint64_t q, const uint64_t t, const size_t l, const Zp1 *y, const Rp1GSW *C, const KSWKeyRp1 *S);
template RLWE1<Rp2> ext_exp_inner<Rp2, Zp2, P2>(const uint64_t q, const uint64_t t, const size_t l, const Zp2 *y, const Rp2GSW *C, const KSWKeyRp2 *S);

template RLWE<Rp1_crt, 1> accumulation<Rp1, Zp1, Rp1_crt, P1>(const RLWE<CirculantRing<Zp12, 1, 1>, N> &ab, const Rp1GSW BS[N], const KSWKeyRp1 KS[P1]);
template RLWE<Rp2_crt, 1> accumulation<Rp2, Zp2, Rp2_crt, P2>(const RLWE<CirculantRing<Zp12, 1, 1>, N> &ab, const Rp2GSW BS[N], const KSWKeyRp2 KS[P2]);

template Rz       dot_product<Rz,      P1>(const Rz * const v1, const Rz * const v2);
template Rz       dot_product<Rz,       N>(const Rz * const v1, const Rz * const v2);
template Rp1      dot_product<Rp1,      1>(const Rp1 * const v1, const Rp1 * const v2);
template Rp2      dot_product<Rp2,      1>(const Rp2 * const v1, const Rp2 * const v2);
template Rp12     dot_product<Rp12,     1>(const Rp12 * const v1, const Rp12 * const v2);
template Rp12     dot_product<Rp12,     3>(const Rp12 * const v1, const Rp12 * const v2);
template Rp1_crt  dot_product<Rp1_crt,  1>(const Rp1_crt * const v1, const Rp1_crt * const v2);
template Rp2_crt  dot_product<Rp2_crt,  1>(const Rp2_crt * const v1, const Rp2_crt * const v2);
template Rp12_crt dot_product<Rp12_crt, 3>(const Rp12_crt * const v1, const Rp12_crt * const v2);
template CirculantRing<Zp12, 1, 1> dot_product<CirculantRing<Zp12, 1, 1>, N>(const CirculantRing<Zp12, 1, 1> * const v1, const CirculantRing<Zp12, 1, 1> * const v2);
