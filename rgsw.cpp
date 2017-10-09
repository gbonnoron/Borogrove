/**
  * @file rgsw.cpp
  * @date Spring 2017
  * @author Guillaume Bonnoron
  * @brief Definition of RGSW encryption and other operations
  */

#include <iostream>
#include "operations.h"
#include "rgsw.h"

template <class R>
template<class Rt>
void RGSW<R>::encrypt(const R* const s, const Rt &m, const double variance)
{
    for (size_t i = 0 ; i < 2*K_Q ; ++i)
    {
        for (size_t j = 0 ; j < 2-1 ; ++j)
            A[i][j] = R::sample_a();

        const R e = R::sample_e(variance);
        A[i][2-1] = dot_product<R, 2-1>(s, A[i]) + e;
        //std::cout << "a.s\t" << dot_product<R, 2-1>(A[i], s) << std::endl;
        //std::cout << "e\t" << e << std::endl;
        //std::cout << "a.s+e\t" << A[i][2-1] << std::endl;
    }
    R mu;
    mu = U * R(m);
    //std::cout << "mu " << mu << std::endl;
    for (size_t i = 0 ; i < K_Q ; ++i)
    {
        const R tmp = G[i][0] * mu;
        //std::cout << "g "<< G[i][0] << std::endl;
        //std::cout << "mu*g "<< tmp << std::endl;
        for (size_t j = 0 ; j < 2 ; ++j)
            A[i + K_Q * j][j] += tmp;
    }
    for (size_t i = 0 ; i < 2*K_Q ; ++i)
        for (size_t j = 0 ; j < 2 ; ++j)
            A[i][j].decomp_fft(); //FFTf A
}
template <class R>
template <class Rt>
void RGSW<R>::decrypt(const R* const s, Rt &dec) const
{
    R cur_val = A[2 * K_Q - (K_Q-1) - 1][2-1];
    cur_val.recomp_fft(SPLIT_Q);
    for (size_t j = 0 ; j < 2-1 ; ++j)
    {
        R tmp = s[j];
        tmp *= A[2 * K_Q - (K_Q-1) - 1][j];  //FFTb before adding
        cur_val -= tmp;
    }
    dec = cur_val.exact_rounding(Q, T);
}

template class RGSW<Rp1>;
template class RGSW<Rp2>;

template void Rp1GSW::encrypt<CirculantRing<Zt, P1> >(const Rp1* const s, const CirculantRing<Zt, P1> &m, const double variance);
template void Rp2GSW::encrypt<CirculantRing<Zt, P2> >(const Rp2* const s, const CirculantRing<Zt, P2> &m, const double variance);

//For tests only
template void Rp1GSW::decrypt<CirculantRing<Zt, P1> >(const Rp1* const s, CirculantRing<Zt, P1> &dec) const;
