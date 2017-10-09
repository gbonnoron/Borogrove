/**
  * @file ksw_key.h
  * @date Spring 2017
  * @author Guillaume Bonnoron
  * @brief Definition of Key-Switching key class
  */

#include <assert.h>
#include "ksw_key.h"
#include "rlwe.h"

template<class R, uint64_t n, uint64_t np, uint64_t B, uint64_t K>
KSWKey<R, n, np, B, K>::KSWKey() : init_ok(false) {}

template<class R, uint64_t n, uint64_t np, uint64_t B, uint64_t K>
KSWKey<R, n, np, B, K>::KSWKey(const R s[n], const R sp[np], const double variance)
{
    //std::cout << "KSWKey in basis " << B << " of size " << K << std::endl;
    for (size_t i = 0 ; i < n ; ++i)
    {
        //std::cout << i << " " << std::flush;
        int64_t pow = 1;
        //std::cout << s[i] << std::endl;
        for (size_t j = 0 ; j < K ; ++j)
        {
            S[i][j] = RLWE<R, np>();
            const R tmp = s[i] * pow;
            //std::cout << "-> " << tmp << std::endl;
            S[i][j].encrypt(sp, tmp, 1, 1, variance);
            S[i][j].decomp_fft();//FFTf S[i][j]
            pow *= B;
        }
    }
    init_ok = true;
}

template<class R, uint64_t n, uint64_t np, uint64_t B, uint64_t K>
void KSWKey<R, n, np, B, K>::init(const R s[n], const R sp[np], const double variance)
{
    for (size_t i = 0 ; i < n ; ++i)
    {
        int64_t pow = 1;
        //std::cout << s[i] << std::endl;
        for (size_t j = 0 ; j < K ; ++j)
        {
            S[i][j] = RLWE<R, np>();
            R tmp = s[i] * pow;
            //tmp.decomp_fft();
            //std::cout << "-> " << tmp << std::endl;
            S[i][j].encrypt(sp, tmp, 1, 1, variance);
            S[i][j].decomp_fft();//FFTf S[i][j]
            pow *= B;
        }
    }
    init_ok = true;
}

template class KSWKey<Rp1>;
template class KSWKey<Rp2>;
template class KSWKey<Rp12, 3, 1, B_Qp, K_Qp>;
template class KSWKey<Rz, P1, N, B_Qp_2, K_Qp_2>;
