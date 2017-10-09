#ifndef HE8_FFT_H
#define HE8_FFT_H

/**
  * @file fft.h
  * @date Spring 2017
  * @author Guillaume Bonnoron
  * @brief Declaration of FFT utility
  */

#include <fftw3.h>

/**
 * @class FFT
 * @brief Self-contained FFT calculation object
 * @tparam n the FFT dimension
 *
 * The input data should be provided to the in
 */
template<size_t n>
class FFT
{
private:
    double *in;                 ///< array for input data
    fftw_complex *out;          ///< array for output data
    fftw_plan forward_plan;     ///< Forward plan for FFTW
    fftw_plan backward_plan;    ///< Backward plan for FFTW
public:
    /**
     * @brief Default constructor that does the memory allocation and the plans setup.
     */
    FFT();
    /**
      * @brief Default destructor
      */
    ~FFT();
    /**
     * @brief Computes a forward FFT
     * @param res the resulting FFT representation
     */
    void forward(fftw_complex *res);
    /**
     * @brief Computes a backward FFT
     * @param data the FFT representation to invert
     */
    void backward(const fftw_complex *data);
    /**
     * @brief Accessor to the input array for initialization
     * @return the pointer to the internal input array
     */
    inline double *get_in() const
    {
        return in;
    }
};

namespace fft {

/**
 * @brief Compute the component-wise product of the fft representation
 * @param res_fft resulting coefficients
 * @param lhs_fft lefthand side coefficients
 * @param rhs_fft righthand side coefficients
 */
template <size_t n>
void product(fftw_complex *res_fft, const fftw_complex * lhs_fft, const fftw_complex * rhs_fft);
/**
 * @brief Compute the FFT of the tensor product x (x) y
 * @param res_fft resulting coefficients
 * @param lhs_fft lefthand side coefficients
 * @param rhs_fft righthand side coefficients
 */
template< size_t p, size_t q>
void tensor_product(fftw_complex *res_fft, const fftw_complex * lhs_fft, const fftw_complex * rhs_fft);

}

#endif // HE8_FFT_H
