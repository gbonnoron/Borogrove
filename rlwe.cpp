/**
  * @file rlwe.cpp
  * @date Spring 2017
  * @author Guillaume Bonnoron
  * @brief Definition of RLWE encryption and other operations
  */

#include <iostream>
#include <memory>
#include <cmath>
#include "rlwe.h"
#include "param.h"
#include "gadget.h"
#include "operations.h"
#include "ksw_key.h"
#include "circulant_ring.h"
//#include "rgsw.h"

template<class R, uint64_t n>
RLWE<R, n>::RLWE()
{
    a = new R[n];
}

template<class R, uint64_t n>
RLWE<R, n>::RLWE(const RLWE &src)
{
    a = new R[n];
    std::copy(src.a, src.a + n, a);
    b = src.b;
}

template<class R, uint64_t n>
RLWE<R, n>::RLWE(R *new_a, R new_b)
{
    a = new R[n];
    if (new_a != nullptr)
        std::copy(new_a, new_a + n, a);
    b = new_b;
}

template<class R, uint64_t n>
RLWE<R, n>::RLWE(R *new_c, bool copy)
{
    b = new_c[n];
    if (copy)
    {
        a = new R[n];
        std::copy(new_c, new_c + n, a);
    }
    else
    {
        a = std::move(new_c);
    }
}

template<class R, uint64_t n>
RLWE<R, n>::~RLWE()
{
    delete[] a;
}

template<class R, uint64_t n>
template <class Zp>
void RLWE<R, n>:: mod(Zp* ab) const
{
    for (size_t i = 0 ; i < n ; ++i)
        ab[i] = Zp((int64_t)a[i].get_data()[0]);
    ab[n] = Zp((int64_t)b.get_data()[0]);
}

template<class R, uint64_t n>
template<class Rt>
void RLWE<R, n>::encrypt(const R* const s, const Rt &m, const uint64_t q, const uint64_t t, const double variance)
{
    const R e = R::sample_e(variance);
    //R e;
    //e = 0;

    for (size_t i = 0 ; i < n ; ++i)
        a[i] = R::sample_a();

    /*
    std::cout << "m\t" << m << std::endl;
    std::cout << "as\t" << dot_product<R, n>(a, s) << std::endl;
    std::cout << "e\t" << e << std::endl;
    std::cout << "q/t\t" << R(q/t) << std::endl;
    std::cout << "mq/t\t" << R(m)*R(q/t) << std::endl;
    std::cout << "as+e+mq/t\t" <<  dot_product<R, n>(a, s) + e + R(m)*R(q/t) << std::endl;
    */
    //b = m;
    //b *= (q/t);
    size_t size = Rt::get_d();
    int64_t *b_data = new int64_t[size];
    for (size_t i = 0 ; i < size ; ++i)
        b_data[i] = floor((double)(m.get_data()[i]) * (double) q / (double) t + 0.5);
    b = R(b_data);
    delete[] b_data;

    //b = llround((double)m * (double) q / (double) t);
    b += dot_product<R, n>(s, a) + e;
    return;
}

template<class R, uint64_t n>
template<class Rt>
void RLWE<R, n>::decrypt(const R* const s, Rt &dec, const uint64_t q, const uint64_t t) const
{
    /*
    std::cout << "b\t\t" << b << std::endl;
    std::cout << "as\t\t" <<dot_product<R, n>(a, s) << std::endl;
    std::cout << "b-as\t\t" <<b - dot_product<R, n>(a, s) << std::endl;
    */
    dec = (b - dot_product<R, n>(s, a)).exact_rounding(q, t);
    //std::cout << "dec\t\t" << dec << std::endl;
}

template<class R, uint64_t n>
double RLWE<R, n>::noise(const R* const s, const uint64_t q, const uint64_t t) const
{
    R res = b - dot_product<R, n>(s, a);
    R signed_noise = res - (q/t) * res.exact_rounding(q, t);
    //std::cout << "noise = " << signed_noise << std::endl;
    return signed_noise.get_norm();
}

template<class R, uint64_t n>
template<class Rp, uint64_t p, uint64_t d_fft>
RLWE<Rp, n> RLWE<R, n>::trace() const
{
    Rp new_c[n + 1];
    for (size_t i = 0 ; i < n ; ++i)
        new_c[i] = a[i].template trace<p, d_fft>();
    new_c[n] = b.template trace<p, d_fft>();
    return RLWE<Rp, n>(new_c);
}

template<class R, uint64_t n>
template <class R2, uint64_t q2, uint64_t q>
RLWE<R2, n> RLWE<R, n>::mod_switch() const
{
    R2 a_dst[n], b_dst;
    for (size_t i = 0 ; i < n ; ++i)
        a_dst[i] = a[i].template rounding<R2>(q, q2);
    b_dst = b.template rounding<R2>(q, q2);
    return RLWE<R2, n>(a_dst, b_dst);
}

template<class R, uint64_t n>
template <uint64_t np, uint64_t B, uint64_t K>
RLWE<R, np> RLWE<R, n>::key_switch(const KSWKey<R, n, np, B, K> &ksw) const
{
    R ai_decomp[n * K];
    for (size_t i = 0 ; i < n ; ++i)
        Gadget<R, B, K>::g_invT(ai_decomp + i * K, a[i]);
    R* res = R::keymult(ai_decomp, ksw);
    for (size_t i = 0 ; i < np ; ++i)
        res[i] = - res[i];
    res[np] = b - res[np];

    return RLWE<R, np>(std::move(res), false);
}

template<class R, uint64_t n>
RLWE<R, n> RLWE<R, n>::ext_mult(const RGSW<R> &rgsw) const
{
    RLWE<R, n> result(*this);
    result.ext_mult_inplace(rgsw);
    return result;
}

template<class R, uint64_t n>
void RLWE<R, n>::ext_mult_inplace(const RGSW<R> &rgsw)
{
    R v[(n + 1) * K_Q];
    for (size_t i = 0 ; i < n ; ++i)
        Gadget<R, B_Q, K_Q>::g_invT(v + i * K_Q, U_inv * a[i]);
    Gadget<R, B_Q, K_Q>::g_invT(v + n * K_Q, U_inv * b);
    R* res = v*rgsw;

    std::copy(res, res+n, a);
    b = res[n];
    delete[] res;
}

template<class R, uint64_t n>
void RLWE<R, n>::galois_inplace(const uint64_t alpha)
{
    for (size_t i = 0 ; i < n ; ++i)
        a[i].galois_inplace(alpha);
    b.galois_inplace(alpha);
}

template<class R, uint64_t n>
RLWE<R, n> RLWE<R, n>::galois(const uint64_t alpha) const
{
    RLWE<R, n> result(*this);
    result.galois_inplace(alpha);
    return result;
}

template<class R, uint64_t n>
RLWE<R, n>& RLWE<R, n>::operator=(const RLWE& other)
{
    if (this != &other)
    {
        std::copy(other.a, other.a + n, a);
        b = other.b;
    }
    return *this;
}

template<class R, uint64_t n>
RLWE<R, n>& RLWE<R, n>::operator=(RLWE&& other) noexcept
{
    if (this != &other)
    {
        delete[] a;
        //C++14 : a = std::exchange(other.a, nullptr);
        R *old_value = std::move(other.a);
        other.a = std::forward<R *>(nullptr);
        a = old_value;

        b = other.b;
    }
    return *this;
}

template<class R>
void RLWE1<R>::ext_exp_mult_add(const RGSW<R> &rgsw, const uint64_t alpha, const KSWKey<R, 1, 1, B_Q, K_Q>& S_alpha, const uint64_t beta, const KSWKey<R, 1, 1, B_Q, K_Q>& S_beta)
{
    if (beta != 1)
    {
        galois_inplace(beta);
        key_switch_inplace<B_Q, K_Q>(S_beta);
    }
    ext_mult_inplace(rgsw);
    if (alpha != 1)
    {
        galois_inplace(alpha);
        key_switch_inplace<B_Q, K_Q>(S_alpha);
    }
}
template<class R>
template <uint64_t B, uint64_t K>
void RLWE1<R>::key_switch_inplace(const KSWKey<R, 1, 1, B, K> &ksw)
{
    R ai_decomp[1 * K];
    for (size_t i = 0 ; i < 1 ; ++i)
        Gadget<R, B, K>::g_invT(ai_decomp + i * K, a[i]);
    R* res = R::keymult(ai_decomp, ksw);

    a[0] = -res[0];
    b -= res[1];
    delete[] res;
}

template class RLWE<Rz      , P1>;
template class RLWE<Rz      , N>;
template class RLWE<CirculantRing<Zp12, 1, 1>, N>;
template class RLWE<Rp1     , 1>;
template class RLWE<Rp2     , 1>;
template class RLWE<Rp1_crt , 1>;
template class RLWE<Rp2_crt , 1>;
template class RLWE<Rp12_crt, 3>;
template class RLWE<Rp12    , 3>;
template class RLWE<Rp12    , 1>;

// LWE
template RLWE<Rz, N> RLWE<Rz, P1>::key_switch<N, B_Qp_2, K_Qp_2>(const KSWKeyLWE &ksw) const;
template RLWE<CirculantRing<Zp12, 1, 1>, N> RLWE<Rz, N>::mod_switch<CirculantRing<Zp12, 1, 1>, P1*P2, Qp>() const;

// ExtExpInner
template void RLWE<CirculantRing<Zp12, 1, 1>, N>::mod<Zp1>(Zp1* ab) const;
template void RLWE<CirculantRing<Zp12, 1, 1>, N>::mod<Zp2>(Zp2* ab) const;
template void RLWE1<Rp1>::key_switch_inplace<B_Q, K_Q>(const KSWKeyRp1 &ksw);
template void RLWE1<Rp2>::key_switch_inplace<B_Q, K_Q>(const KSWKeyRp2 &ksw);
template void Rp1LWE::ext_exp_mult_add(const Rp1GSW &rgsw, const uint64_t alpha, const KSWKeyRp1 &S_alpha, const uint64_t beta,  const KSWKeyRp1 &S_beta);
template void Rp2LWE::ext_exp_mult_add(const Rp2GSW &rgsw, const uint64_t alpha, const KSWKeyRp2 &S_alpha, const uint64_t beta,  const KSWKeyRp2 &S_beta);

// ExpCRT
template RLWE<Rp1_crt, 1> RLWE<Rp1, 1>::mod_switch<Rp1_crt, Qcrt, Q>() const;
template RLWE<Rp2_crt, 1> RLWE<Rp2, 1>::mod_switch<Rp2_crt, Qcrt, Q>() const;
template Rp12LWE RLWE<Rp12_crt, 3>::mod_switch<Rp12, Qp, Qcrt>() const;

// FunExpExtract
template RLWE<Rp12, 1> Rp12LWE::key_switch<1, B_Qp, K_Qp>(const KSWKeyRp12 &ksw) const;
template RLWE<Rp1_p, 1> RLWE<Rp12 , 1>::trace<Rp1_p, P1, FFT_DIM1>() const;

// input/output
template void LWE::encrypt<CirculantRing<Zt, 1, 1> >(const Rz *s, const CirculantRing<Zt, 1, 1> &m, const uint64_t q, const uint64_t t, const double variance);
template void LWE::decrypt<CirculantRing<Zt, 1, 1> >(const Rz *s, CirculantRing<Zt, 1, 1> &dec, const uint64_t q, const uint64_t t) const;

// Key-switching keys
template void RLWE<Rz, N>::encrypt<Rz>(const Rz *s, const Rz &m, const uint64_t q, const uint64_t t, const double variance);
template void RLWE<Rp1, 1>::encrypt<Rp1>(const Rp1 *s, const Rp1 &m, const uint64_t q, const uint64_t t, const double variance);
template void RLWE<Rp2, 1>::encrypt<Rp2>(const Rp2 *s, const Rp2 &m, const uint64_t q, const uint64_t t, const double variance);
template void RLWE<Rp12, 1>::encrypt<Rp12>(const Rp12 *s, const Rp12&m, const uint64_t q, const uint64_t t, const double variance);


// For tests only
template void RLWE<Rp1, 1>::encrypt<CirculantRing<Zt, P1> >(const Rp1 *s, const CirculantRing<Zt, P1> &m, const uint64_t q, const uint64_t t, const double variance);
template void RLWE<Rp1, 1>::decrypt<CirculantRing<Zt, P1> >(const Rp1 *s, CirculantRing<Zt, P1> &dec, const uint64_t q, const uint64_t t) const;
template void Rp12LWE::encrypt<CirculantRing<Zt, P1*P2, FFT_DIM2> >(const Rp12 *s, const CirculantRing<Zt, P1*P2, FFT_DIM2> &m, const uint64_t q, const uint64_t t, const double variance);
template void Rp12LWE::decrypt<CirculantRing<Zt, P1*P2, FFT_DIM2> >(const Rp12 *s, CirculantRing<Zt, P1*P2, FFT_DIM2> &dec, const uint64_t q, const uint64_t t) const;

template void RLWE<Rz, N>::encrypt<CirculantRing<Zt, 1, 1> >(const Rz *s, const CirculantRing<Zt, 1, 1>  &m, const uint64_t q, const uint64_t t, const double variance);
template void RLWE<Rp2, 1>::encrypt<CirculantRing<Zt, P2> >(const Rp2 *s, const CirculantRing<Zt, P2> &m, const uint64_t q, const uint64_t t, const double variance);
template void RLWE<CirculantRing<Zp12, 1, 1>, N>::decrypt<CirculantRing<Zt, 1, 1> >(const CirculantRing<Zp12, 1, 1> *s, CirculantRing<Zt, 1, 1> &dec, const uint64_t q, const uint64_t t) const;
template void RLWE<Rp12, 1>::decrypt<CirculantRing<Zt, P1*P2, FFT_DIM2> >(const Rp12 *s, CirculantRing<Zt, P1*P2, FFT_DIM2> &dec, const uint64_t q, const uint64_t t) const;

// For stats only
template void RLWE<Rp1_crt, 1>::encrypt<CirculantRing<Zt, P1> >(const Rp1_crt *s, const CirculantRing<Zt, P1> &m, const uint64_t q, const uint64_t t, const double variance);
template void RLWE<Rp2_crt, 1>::encrypt<CirculantRing<Zt, P2> >(const Rp2_crt *s, const CirculantRing<Zt, P2> &m, const uint64_t q, const uint64_t t, const double variance);
template void RLWE<Rp12_crt, 3>::encrypt<CirculantRing<Zt, P1*P2, FFT_DIM2> >(const Rp12_crt *s, const CirculantRing<Zt, P1*P2, FFT_DIM2> &m, const uint64_t q, const uint64_t t, const double variance);
