#ifndef HE8_OPERATIONS_H
#define HE8_OPERATIONS_H

/**
  * @file operations.h
  * @date Spring 2017
  * @author Guillaume Bonnoron
  * @brief Declaration of some elementary operations of the scheme
  *
  * Some operations are not directly included in RLWE or RGSW classes. They are regrouped here.
  */

//#include "param.h"
#include "rlwe.h"
#include "predicate.h"
#include "ksw_key.h"

void print_param();

/**
 * @brief Compute the dot product of two vectors
 * @param v1 the first vector
 * @param v2 the second vector
 * @tparam R the class representing vector elements
 * @tparam n the size of the vectors
 * @return the dot product <v1.v2> in R
 */
template<class R, uint64_t n>
R dot_product(const R* const v1, const R* const v2);

/**
 * @brief ExpCRT operation
 * @param cp the RLWE encryption in Rp
 * @param cq the RLWE encryption in Rq
 * @param cpq the computed RLWE encryption in Rpq
 */
void exp_crt(RLWE<Rp12_crt, 3> &cpq, const RLWE<Rp1_crt, 1> &cp, const RLWE<Rp2_crt, 1> &cq);

/**
 * @brief Computes the key for ciphertexts from ExpCRT
 * @param sp the key of the RLWE encryption in Rp
 * @param sq the key of the RLWE encryption in Rq
 * @param spq the resulting key for the RLWE encryption in Rpq
 */
void crt_key(Rp12 *spq, const Rp1 &sp, const Rp2 &sq);

/**
 * @brief ExtExpInner operation (optimized)
 * @param q the modulus of the encryptions
 * @param t the plaintext modulus
 * @param l the length of the inner product
 * @param y the l plaintext elements
 * @param C the l RGSW encryptions of T^{x_i}
 * @param S the array of key-switching keys
 * @tparam R the class of the ring elements
 * @tparam Zp the class of elements of y
 * @tparam p the modulus of the input elements
 * @tparam B the integer basis for the decomposition
 * @tparam K the size of the decomposition
 * @return an RLWE encryption of T^{<x, y>}
 */
template<class R, class Zp, uint64_t p>
RLWE1<R> ext_exp_inner(const uint64_t q, const uint64_t t, const size_t l, const Zp *y, const RGSW<R> *C, const KSWKey<R, 1, 1, B_Q, K_Q> *S);

/**
 * @brief fun_exp_extract
 * @param F the function to evaluate on m
 * @param c_Zm a RLWE encryption in Rpq of Z^m
 * @param S a Key-Switching key from s_pq in Rpq to s
 * @return an LWE encryption of dimension n of F(m)
 */
LWE fun_exp_extract(const fftw_complex *f, const Rp12LWE &c_Zm, const KSWKeyRp12 &S);

/**
 * @brief Accumulation, part of the preparation phase
 * @param ab encryption of m mod pq
 * @param BS Bootstrapping keys (RGSW encryptions of s_i)
 * @param KS KeySwitching Keys phi_alpha(s)-> s for alpha in Zp1
 * @tparam R the ring of the key elements
 * @tparam Z the class of the intermediary elements
 * @tparam R_p the ring of the resulting elements
 * @tparam p the modulus
 * @return an RLWE encryption of Z^{-<a.s> mod p}
 */
template <class R, class Z, class R_p, uint64_t p>
RLWE<R_p, 1> accumulation(const RLWE<CirculantRing<Zp12, 1, 1>, N> &ab, const RGSW<R> BS[N], const KSWKey<R> KS[p]);

/**
 * @brief Preparation phase, from LWE(m) to Rp12LWE(Z^m)
 * @param Em LWE encryption of m
 * @param Xsi Rp1GSW encryptions of s_i
 * @param KSp1 KeySwitching Keys phi_alpha(s)-> s for alpha in Zp1
 * @param Ysi Rp2GSW encryptions of s_i
 * @param KSp2 KeySwitching Keys phi_alpha(s)-> s for alpha in Zp2
 * @return an Rp12LWE encryption of Z^m
 */
Rp12LWE preparation(const LWE &Em, const KSWKeyLWE &S_lwe, const Rp1GSW Xsi[N], const KSWKeyRp1 KSp1[P1], const Rp2GSW Ysi[N], const KSWKeyRp2 KSp2[P2]);

/**
 * @brief Combination phase, from LWE(m_k) to LWE(m) where m = coef_0*m_0 + coef_1*m_1 + coef_{k-1}*m_{k-1}
 * @param k the number of element to combine
 * @param coefs the k coefs for the combination
 * @param c_i the elements
 * @return the combination of the inputs
 */
LWE combination(const size_t k, const int64_t *coefs, const LWE *c_i);

/**
 * @brief Compute the function F on the input data
 * @param k the number of inputs
 * @param coefs the coefficients for the combination
 * @param c_i the inputs (LWE encryptions)
 * @param F the function to compute
 * @param Xsi the bootstrapping keys in the first ring
 * @param KSp1 the Key-switching keys in the first ring
 * @param Ysi the bootstrapping keys in the second ring
 * @param KSp2 the Key-switching keys in the second ring
 * @param S the key-switching key for the function extraction
 * @return an LWE encryption of F(m_1, ..., m_k)
 */
LWE gate(const size_t k, const int64_t *coefs, const LWE *c_i, const fftw_complex *f, const KSWKeyLWE &S_lwe, const Rp1GSW Xsi[N], const KSWKeyRp1 KSp1[P1], const Rp2GSW Ysi[N], const KSWKeyRp2 KSp2[P2], const KSWKeyRp12 &S);
/**
 * @brief Generate the bootstrapping keys for ExtExpInner
 * @param q ciphertext modulus
 * @param t plaintext modulus
 * @param Tsi the RGSW encryptions of T^s_i
 * @param s the LWE key
 * @param s_p the key for the RGSW encryptions
 * @tparam R the class for the ring elements
 * @tparam n the length of the LWE key
 * @tparam p the dimension of the key
 */
template <class R, uint64_t n, int64_t p>
void gen_bootstrapping_keys(const uint64_t q, const uint64_t t, RGSW<R> Tsi[n], const Rz s[n], const R &s_p, const double variance);

/**
 * @brief Generate the key-switching keys for ExtExpInner
 * @param keys the resulting keys
 * @param s_p the secret key for which to compute the KS keys
 * @tparam R the class for the ring elements
 * @tparam p the dimension of the key
 * @tparam B the integer basis for the decomposition
 * @tparam K the size of the decomposition
 */
template <class R, uint64_t p>
void gen_keyswitching_keys(KSWKey<R, 1, 1, B_Q, K_Q> keys[p], const R &s_p, const double variance);
/**
 * @brief Generate the key-switching key for FunExpExtract
 * @param S the key-switching key
 * @param s_pq the key in R_pq after ExpCRT
 * @param s the LWE key
 */
void gen_funexpextract_key(KSWKeyRp12 *S, const Rp12 s_pq[3], const Rz s[P1], const double variance);

#endif // HE8_OPERATIONS_H
