#ifndef HE8_GADGET_H
#define HE8_GADGET_H

/**
  * @file gadget.h
  * @date Spring 2017
  * @author Guillaume Bonnoron
  * @brief Definition of Gadget tool
  *
  * The Gadget allows to compute decomposition of integer elements into a basis B.
  * It has no internal data, just class methods.
  * The class is templated on the ring it operates on, the basis that is used for the decomposition and the size of these decompositions.
  */

#include <cmath>
#include <cstring>

/**
 * @class Gadget
 * @brief The tool for decomposition of elements in R smaller than B^K-1 into basis B
 * @tparam R the class representing elements
 * @tparam B the integer basis for the decomposition
 * @tparam K the size of the decomposition
 */
template <class R, uint64_t B, uint64_t K>
class Gadget
{
public:
    /**
     * @brief Computes the gadget matrix of size n: G_n = I_{n+1} (x) g
     * @tparam n the size of the gadget
     * @return the gadget matrix
     */
    template<uint64_t n>
    static void gadget(int64_t G[n*K][n])
    {
        for (size_t i = 0 ; i < n * K ; ++i)
            for (size_t j = 0 ; j < n ; ++j)
                G[i][j] = 0;
        int64_t Bexpo = 1;
        for (size_t i = 0 ; i < K ; ++i)
        {
            for (size_t j = 0 ; j < n ; ++j)
                G[i + K*j][j] = Bexpo;
            Bexpo *= B;
        }
    }
    /**
     * @brief Compute the decomposition
     * @param decomp the returned decomposition
     * @param val the element to decompose
     */
    static void g_invT(R *decomp, const R &val)
    {
        //std::cout << "Signed decomp in basis " << B << " of size " << K << std::endl;
        //std::cout << val << std::endl;
        R tmp(val);
        for (size_t expo = 0 ; expo < K ; ++expo)
        {
            decomp[expo] = tmp % B;
            decomp[expo].balance(B);
            //std::cout << tmp << std::endl;
            //std::cout << "=> " << Bexpo << " *\t" <<decomp[expo] << std::endl;
            tmp -= decomp[expo];
            tmp = tmp / B;
        }
    }
    /**
     * @brief Compute the decomposition and outputs the results in FFT
     * @param dec_fft the returned decomposition in FFT
     * @param val the element to decompose
     */
    template <class Ring>
    static void g_invT_fft(fftw_complex *dec_fft[K], const Ring& val)
    {
        Ring decomp[K];
        g_invT(decomp, val);
        for (size_t i = 0 ; i < K ; ++i)
            decomp[i].compute_fft(dec_fft[i]);
    }
    /**
     * @brief Compute the decomposition of a pure tensor
     * @param decomp the returned decomposition, in FFT
     * @param val1, val2 the elements whose tensor is to be decomposed
     * /
    template<class R2>
    static void g_invT_tensor(fftw_complex *decomp[K], const R &val1, const R2&val2)
    {
        /// Compute val1 decomposition in FFT
        fftw_complex *dec1_fft[K];
        for (size_t i = 0 ; i < K ; ++i)
            dec1_fft[i] = fftw_alloc_complex(P1);
        g_invT_fft<R>(dec1_fft, val1);

        /// Compute val2 decomposition in FFT
        fftw_complex *dec2_fft[K];
        for (size_t i = 0 ; i < K ; ++i)
            dec2_fft[i] = fftw_alloc_complex(P2);
        g_invT_fft<R2>(dec2_fft, val2);

        /// Initialize the FFT coefficients of the result
        for (size_t j = 0 ; j < K ; ++j)
            std::memset(decomp[j], 0, (P1*P2/2+1) * 2 * sizeof(double));

        fftw_complex *fft_bigcache = fftw_alloc_complex(P1*P2);
        for (size_t j = 0 ; j < K ; ++j)
        {
            for (size_t i = 0 ; i <= j ; ++i)
            {
                tensor_product(fft_bigcache, dec1_fft[i], dec2_fft[j-i]);
                for (size_t k = 0 ; k < P1*P2 ; ++k)
                {
                    decomp[j][k][0] += fft_bigcache[k][0];
                    decomp[j][k][1] += fft_bigcache[k][1];
                }
            }
        }
        fftw_free(fft_bigcache);

        for (size_t i = 0 ; i < K ; ++i)
        {
            fftw_free(dec1_fft[i]);
            fftw_free(dec2_fft[i]);
        }
    }
    */
    /**
     * @brief Accessor to the basis for the decomposition
     * @return B
     */
    static inline uint64_t get_basis()
    {
        return B;
    }
    /**
     * @brief Accessor to the size of the decomposition
     * @return K
     */
    static inline uint64_t get_size()
    {
        return K;
    }
};

#endif // HE8_GADGET_H
