#ifndef HE8_RING_H
#define HE8_RING_H

/**
  * @file circulant_ring.h
  * @date Spring 2017
  * @author Guillaume Bonnoron
  * @brief Declaration of the CirculantRing class
  *
  * The CirculantRing class provides everything needed to represent and operate with elements in the circulant ring of degree d, eg. R = Z[T]/(T^d-1)
  * The degree is template for better performance. Several constructors and operators are defined for ease of use.
  * For each degree an FFT class is associated for efficient multiplication.
  */

#include <iostream>
#include <utility>
#include <memory>
#include <cstring>
#include <assert.h>
#include "param.h"
#include "fft.h"
#include "integer_mod.h"
#include "predicate.h"

template <class R>
class RGSW;
template<class R, uint64_t n, uint64_t np, uint64_t B, uint64_t K>
class KSWKey;

/**
 * @class CirculantRing
 * @brief Class representing elements in the circulant ring of degree d.
 * @tparam d the degree of the ring
 */
template <class Int, uint64_t d, uint64_t d_fft = FFT_DIM1>
class CirculantRing
{
private:
    static FFT<d_fft> fft;   ///< FFT object associated to the CirculantRing class, for each degree
    Int *x;                 ///< Internal representation of the ring element
    fftw_complex *msb_fft;
    fftw_complex *lsb_fft;
    bool fft_status;
    static Int tmp[d];
    static fftw_complex fft_cache1[d_fft/2+1] __attribute__ ((aligned (32)));
    static fftw_complex fft_cache2[d_fft/2+1] __attribute__ ((aligned (32)));
    static fftw_complex fft_cache3[d_fft/2+1] __attribute__ ((aligned (32)));
public:
    /**
     * @brief Default constructor, doing only memory allocation
     */
    CirculantRing();
    /**
     * @brief Copy constructor
     * @param src ring element to duplicate
     */
    CirculantRing (const CirculantRing &src);
    /**
     * @brief Copy constructor where src has another element representation
     * @param src ring element to reinterpret and duplicate
     */
    template <class Int2>
    CirculantRing (const CirculantRing<Int2, d, d_fft> &src);
    /**
     * @brief Copy constructor where src has another element representation
     * @param src ring element to reinterpret and duplicate
     */
    template <uint64_t d_fft2>
    CirculantRing (const CirculantRing<Int, d, d_fft2> &src);
    /**
     * @brief Constructs a ring element from a representation
     * @param values an array corresponding to the internal representation
     * @param copy indicates if the representation is copied to a new internal array or just referenced
     */
    CirculantRing (Int *values, bool copy=true);
    /**
     * @brief Constructs a ring element from a representation
     * @param msb_fft the FFT representation for the most signficant part
     * @param lsb_fft the FFT representation for the least signficant part
     * @param copy indicates if the representation is copied to a new internal array or just referenced
     */
    CirculantRing (fftw_complex *msb_fft_values, fftw_complex *lsb_fft_values, bool copy=true);
    /**
     * @brief Constructs a ring element from a constant term
     * @param val the constant term in the internal representation (Int)
     */
    CirculantRing (const Int val);
    CirculantRing (const int64_t *values);
    CirculantRing (const Fun &F);

    /**
     * @brief Default destructor
     */
    ~CirculantRing();
    /**
     * @brief Accessor to the FFT dimension
     * @return the FFT dimension
     */
    static inline size_t get_d_fft()
    {
        return d_fft;
    }

    static inline size_t get_d()
    {
        return d;
    }
    /**
     * @brief Accessor to the internal representation
     * @return a pointer to the internal representation
     */
    inline Int *get_data() const
    {
        assert(fft_status == false);
        return x;
    }
    /**
     * @brief Accessor to the msb FFT internal representation
     * @return a pointer to the msb FFT internal representation
     */
    inline fftw_complex *get_msb_fft() const
    {
        assert(fft_status);
        return msb_fft;
    }
    /**
     * @brief Accessor to the lsb FFT internal representation
     * @return a pointer to the lsb FFT internal representation
     */
    inline fftw_complex *get_lsb_fft() const
    {
        assert(fft_status);
        return lsb_fft;
    }
    /**
     * @brief Accessor to the internal flag fft_status
     * @return its value
     */
    inline bool get_fft_status() const
    {
        return fft_status;
    }
    /**
     * @brief Compute the 2-norm of the coefficient representation
     * @return the norm
     */
    inline double get_norm() const
    {
        double res = 0;
        for (size_t i = 0 ; i < d ; ++i)
        {
            //std::cout << (double)x[i] << std::endl;
            res += (double)x[i] * (double)x[i];
            //std::cout << res << std::endl;
        }
        return res;
    }
    inline bool is_zero() const
    {
        bool result = true;
        for (size_t i = 0 ; i < d && result ; ++i)
            result = result && (x[i] == 0);
        return result;

    }
    static inline FFT<d_fft> *get_fft_obj()
    {
        return &fft;
    }
    /**
     * @brief Class method for Gaussian sampling of ring element
     * @param stddev the standard deviation
     * @return a ring element whose coefficients are drawn from a Gaussian distribution of standard deviation stddev
     */
    static CirculantRing gaussian_sample(const double stddev);
    /**
     * @brief Class method uniformly sampling a ring element
     * @return a ring element whose coefficients are drawn uniformly
     */
    static CirculantRing uniform_sample();
    /**
     * @brief Class method sampling a binary ring element
     * @return a ring element with coefficients in {0, 1}
     */
    static CirculantRing binary_sample();
    /**
     * @brief Class method uniformly sampling a ring element a with contraint a(1) = 0
     * @return a ring element whose d-1 coefficients are drawn uniformly and the last meets the contraint
     */
    static CirculantRing sample_a();
    /**
     * @brief Class method sampling a ring elements with d/3 1, d/3 -1 and 0 elsewhere with operator norm < 0.8 sqrt(d ln d)
     * @return a ring element satisfying the contraints
     */
    static CirculantRing sample_s(const double density = 1/3);
    /**
     * @brief Class method sampling a ring elements serving as error term in encryptions
     * @return a ring element satisfying the contraints for errors
     */
    static CirculantRing sample_e(const double variance = 8);
    /**
     * @brief Compute the fft of the current
     * @param x_fft the FFT coefficients
     */
    void compute_fft(fftw_complex *x_fft) const;
    /**
     * @brief Update x from x_fft
     */
    inline void compute_fft_inv(fftw_complex *x_fft)
    {
        compute_fft_inv(x_fft, x);
    }
    /**
     * @brief Compute the inverse FFT of x_fft and store to values
     * @param x_fft the FFT coefficients
     * @param values the result of the inverse FFT
     */
    void compute_fft_inv(fftw_complex *x_fft, Int *values) const;
    /**
     * @brief Split the polynomial into the most and least significant parts and computes the FFT of each
     */
    void decomp_fft();
    /**
     * @brief Split the polynomial into the most and least significant parts and computes the FFT of each
     */
    void recomp_fft(const int64_t B);
    /**
     * @brief Split the ring element into two halves, one with the most significant part and the other with the least significant part
     * @param msb the most significant part
     * @param lsb the least signigicant part
     * @return the threshold value between the part
     */
    int64_t split(CirculantRing &msb, CirculantRing &lsb) const;
    /**
     * @brief Recompose the ring element from two halves, one with the most significant part and the other with the least significant part
     * @param msb the most significant part
     * @param lsb the least signigicant part
     * @return the threshold value between the part
     */
    void unsplit(CirculantRing &msb, CirculantRing &lsb);
    /**
     * @brief Shift the values by a factor B from least to most significant
     * @param B the shifting factor
     * @return the shifted values
     */
    inline void shift(const int64_t B)
    {
        for (size_t i = 0 ; i < d ; ++i)
            x[i] = ((int64_t)x[i] % B) * B;
    }
    /**
     * @brief Scales the ring element by a factor new_modulus/old_modulus and round (out-of-place)
     * @param old_modulus denominator of the ratio
     * @param new_modulus numerator of the ratio
     * @return the scaled ring element
     */
    CirculantRing exact_rounding(const uint64_t old_modulus, const uint64_t new_modulus) const;
    /**
     * @brief Scales the ring element by a factor new_modulus/old_modulus and probabilistically round it (out-of-place)
     * @param old_modulus denominator of the ratio
     * @param new_modulus numerator of the ratio
     * @return the scaled ring element
     */
    template <class R2>
    R2 rounding(const uint64_t old_modulus, const uint64_t new_modulus) const;

    /**
     * @brief Computes the Trace of the ring element over another circulant ring of degree p (templated)
     * @tparam p the degree of the other circulant ring
     * @tparam d2_fft the fft dimension of the result ring element
     * @return a element in CirculantRing<Int, p>
     */
    template <uint64_t p, uint64_t d2_fft>
    CirculantRing<Int, p, d2_fft> trace() const;

    /**
     * @brief Applies the Galois morphism to the ring element (out-of-place)
     * @param alpha the morphism is psi_alpha : T -> T^alpha
     * @return the resulting ring element
     */
    CirculantRing galois(const uint64_t alpha) const;
    /**
     * @brief Applies the Galois morphism to the ring element (in-place)
     * @param alpha the morphism is psi_alpha : T -> T^alpha
     */
    void galois_inplace(const uint64_t alpha);

    /**
     * @brief Computes the tensor product between the ring element of degree d and another of degree d2
     * @tparam d2 the degree of the other ring
     * @tparam d2_fft the fft dimension of the other ring
     * @tparam d3_fft the fft dimension of the result ring
     * @param rhs a ring element in a circulant ring of degree d2
     * @return the tensor product in the circulant ring of degree d x d2
     */
    template<uint64_t d2, uint64_t d2_fft, uint64_t d3_fft>
    CirculantRing<Int, d*d2, d3_fft> tensor_product(const CirculantRing<Int, d2, d2_fft>& rhs) const;

    /**
     * @brief Balances the coefficients between -B/2 and B/2
     * @param B the range
     */
    inline void balance(const int64_t B)
    {
        for (size_t i = 0 ; i < d ; ++i)
            x[i].balance(B);
    }
    /**
     * @brief Fallback generic multiplication function for cases where ring elements are in one FFT
     * @param rhs the ring element to multiply with
     */
    void mult(const fftw_complex *f);
    /**
     * @brief Fallback generic multiplication function for cases where ring elements are small nor in FFT
     * @param rhs the ring element to multiply with
     */
    void mult(const CirculantRing &rhs);

    /*
    void *operator new[](size_t size)
    {
        size_t space = size * sizeof(Int) + ALIGNMENT;
        void *tmp = new char[space];
        std::align(ALIGNMENT, space-ALIGNMENT, tmp, space);
        return reinterpret_cast<CirculantRing *>(tmp);
    }
    */

    inline CirculantRing& operator+=(const CirculantRing& rhs)
    {
        for (size_t i = 0 ; i < d ; ++i)
            x[i] += rhs.x[i];
        return *this;
    }
    inline CirculantRing& operator-=(const CirculantRing& rhs)
    {
        for (size_t i = 0 ; i < d ; ++i)
            x[i] -= rhs.x[i];
        return *this;
    }
    inline CirculantRing& operator*=(const int64_t& rhs)
    {
        const Int tmp(rhs);
        for (size_t i = 0 ; i < d ; ++i)
            x[i] *= tmp;
        return *this;
    }
    CirculantRing& operator=(const CirculantRing& other);
    CirculantRing& operator=(CirculantRing&& other) noexcept;
    CirculantRing& operator=(const int& other);
    //CirculantRing& operator*=(CirculantRing& rhs);
    CirculantRing& operator*=(const CirculantRing& rhs);
    friend inline CirculantRing operator+(CirculantRing lhs, const CirculantRing& rhs) { lhs += rhs; return lhs; }
    friend inline CirculantRing operator-(CirculantRing lhs, const CirculantRing& rhs) { lhs -= rhs; return lhs; }
    friend inline CirculantRing operator*(CirculantRing lhs, const CirculantRing& rhs) { lhs *= rhs; return lhs; }
    //friend inline CirculantRing operator*(CirculantRing lhs, CirculantRing& rhs) { lhs *= rhs; return lhs; }
    friend inline bool operator!=(const CirculantRing& lhs, const CirculantRing& rhs) { return !(lhs == rhs); }
    friend inline CirculantRing operator-(CirculantRing lhs) { return Int(0) - lhs; }
    friend inline CirculantRing operator*(const int64_t &lhs, const CirculantRing &rhs) { return rhs * lhs; }
    friend CirculantRing operator*(const CirculantRing &lhs, const int64_t& rhs)
    {
        Int *res = (Int *)aligned_alloc(ALIGNMENT, (d * sizeof(Int) / ALIGNMENT + 1) * ALIGNMENT);
        for (size_t i = 0 ; i < d ; ++i)
            res[i] = lhs.x[i] * rhs;
        return CirculantRing(res, false);
    }

    friend CirculantRing operator/(const CirculantRing &lhs, const int64_t& rhs)
    {
        Int *res = (Int *)aligned_alloc(ALIGNMENT, (d * sizeof(Int) / ALIGNMENT + 1) * ALIGNMENT);
        for (size_t i = 0 ; i < d ; ++i)
            res[i] = lhs.x[i] / rhs;
        return CirculantRing(res, false);
    }
    friend inline CirculantRing operator%(CirculantRing lhs, const int64_t& rhs)
    {
        for (size_t i = 0 ; i < d ; ++i)
            //lhs.x[i] = ((int64_t)lhs.x[i]) % rhs;  //TODO better?
            lhs.x[i] = lhs.x[i] % rhs;  //TODO better?
        return lhs;
    }
    friend inline bool operator==(const CirculantRing& lhs, const CirculantRing& rhs)
    {
        bool res = true;
        for (size_t i = 0 ; i < d  && res ; ++i)
            res = res && (lhs.x[i] == rhs.x[i]);
        return res;
    }
    friend std::ostream& operator<<(std::ostream& os, const CirculantRing<Int, d, d_fft>& obj)
    {
        /*
        for (size_t i = d-1 ; i >= 1 ; --i)
            os << obj.x[i] << "*X^" << i << " + ";
        os << obj.x[0] << std::flush;
        */
        os << "(";
        for (size_t i = 0 ; i < d-1 ; ++i)
            os << obj.x[i] << ",";
        os << obj.x[d-1] << ")" << std::flush;
        return os;
    }
    /**
     * @brief Vector-Matrix Multiplication
     * @param val the row vector to multiply to the ciphertext matrix
     * @param rhs the ciphertext
     * @return a vector row vector equal to val * rhs
     */
    friend CirculantRing* operator*(CirculantRing* val, const RGSW<CirculantRing> &rhs)
    {
        static const int64_t B = Int::get_split_threshold();
        CirculantRing *res = new CirculantRing[2];
        CirculantRing res_lsb[2];
        fftw_complex val_fft[2* K_Q][d_fft/2+1];
        for (size_t j = 0 ; j < 2*K_Q ; ++j)
            val[j].compute_fft(val_fft[j]);

        for (size_t i = 0 ; i < 2 ; ++i)
        {
            std::memset(fft_cache1, 0, (d_fft/2+1) * 2 * sizeof(double));
            std::memset(fft_cache2, 0, (d_fft/2+1) * 2 * sizeof(double));
            for (size_t j = 0 ; j < 2*K_Q ; ++j)
            {
                fft::product<d_fft>(fft_cache3, val_fft[j], rhs(j, i).get_msb_fft());
                for (size_t k = 0 ; k < d_fft/2+1 ; ++k)
                {
                    fft_cache1[k][0] += fft_cache3[k][0];
                    fft_cache1[k][1] += fft_cache3[k][1];
                }
                fft::product<d_fft>(fft_cache3, val_fft[j], rhs(j, i).get_lsb_fft());
                for (size_t k = 0 ; k < d_fft/2+1 ; ++k)
                {
                    fft_cache2[k][0] += fft_cache3[k][0];
                    fft_cache2[k][1] += fft_cache3[k][1];
                }
                //res[i] += val[j] * rhs(j, i);
            }
            res[i].compute_fft_inv(fft_cache1);
            res[i].shift(B);

            res_lsb[i].compute_fft_inv(fft_cache2);
            res[i] += res_lsb[i];
        }
        return res;
    }

    /**
     * @brief Vector-Matrix Multiplication
     * @param val the row vector to multiply to the ciphertext matrix
     * @param rhs the ciphertext
     * @return a vector row vector equal to val * rhs
     */
    template <uint64_t n, uint64_t np, uint64_t B, uint64_t K>
    static CirculantRing* keymult(CirculantRing* val, const KSWKey<CirculantRing, n, np, B, K> &rhs)
    {
        static const int64_t thres = Int::get_split_threshold();
        CirculantRing *res = new CirculantRing[np + 1];
        CirculantRing res_lsb;
        fftw_complex *val_fft[n * K];
        for (size_t j = 0 ; j < n*K ; ++j)
        {
            val_fft[j] = (fftw_complex *) aligned_alloc(ALIGNMENT, ((d_fft/2+1) * sizeof(fftw_complex) / ALIGNMENT + 1) * ALIGNMENT);
            val[j].compute_fft(val_fft[j]);
        }

        for (size_t idx = 0 ; idx < np ; ++idx)
        {
            std::memset(fft_cache1, 0, (d_fft/2+1) * 2 * sizeof(double));
            std::memset(fft_cache2, 0, (d_fft/2+1) * 2 * sizeof(double));
            for (size_t i = 0 ; i < n ; ++i)
            {
                for (size_t j = 0 ; j < K ; ++j)
                {
                    fft::product<d_fft>(fft_cache3, val_fft[i * K + j], rhs(i, j).get_a()[idx].get_msb_fft());
                    for (size_t k = 0 ; k < d_fft/2+1 ; ++k)
                    {
                        fft_cache1[k][0] += fft_cache3[k][0];
                        fft_cache1[k][1] += fft_cache3[k][1];
                    }
                    fft::product<d_fft>(fft_cache3, val_fft[i * K + j], rhs(i, j).get_a()[idx].get_lsb_fft());
                    for (size_t k = 0 ; k < d_fft/2+1 ; ++k)
                    {
                        fft_cache2[k][0] += fft_cache3[k][0];
                        fft_cache2[k][1] += fft_cache3[k][1];
                    }
                }
            }
            res[idx].compute_fft_inv(fft_cache1);
            res[idx].shift(thres);

            res_lsb.compute_fft_inv(fft_cache2);
            res[idx] += res_lsb;
            //res[idx] -= ai_decomp[j] * ksw(i, j).get_a()[idx];
        }
        {
            std::memset(fft_cache1, 0, (d_fft/2+1) * 2 * sizeof(double));
            std::memset(fft_cache2, 0, (d_fft/2+1) * 2 * sizeof(double));
            for (size_t i = 0 ; i < n ; ++i)
            {
                for (size_t j = 0 ; j < K ; ++j)
                {
                    fft::product<d_fft>(fft_cache3, val_fft[i * K + j], rhs(i, j).get_b().get_msb_fft());
                    for (size_t k = 0 ; k < d_fft/2+1 ; ++k)
                    {
                        fft_cache1[k][0] += fft_cache3[k][0];
                        fft_cache1[k][1] += fft_cache3[k][1];
                    }
                    fft::product<d_fft>(fft_cache3, val_fft[i * K + j], rhs(i, j).get_b().get_lsb_fft());
                    for (size_t k = 0 ; k < d_fft/2+1 ; ++k)
                    {
                        fft_cache2[k][0] += fft_cache3[k][0];
                        fft_cache2[k][1] += fft_cache3[k][1];
                    }
                }
            }
            res[np].compute_fft_inv(fft_cache1);
            res[np].shift(thres);

            res_lsb.compute_fft_inv(fft_cache2);
            res[np] += res_lsb;
            //res[np] -= ai_decomp[j] * ksw(i, j).get_b();
        }
        for (size_t j = 0 ; j < n*K ; ++j)
            free(val_fft[j]);
        return res;
    }
};

template<class Int, uint64_t d, uint64_t d_fft>
FFT<d_fft> CirculantRing<Int, d, d_fft>::fft;
template<class Int, uint64_t d, uint64_t d_fft>
Int CirculantRing<Int, d, d_fft>::tmp[d];
template<class Int, uint64_t d, uint64_t d_fft>
fftw_complex CirculantRing<Int, d, d_fft>::fft_cache1[d_fft/2+1] __attribute__ ((aligned (32)));
template<class Int, uint64_t d, uint64_t d_fft>
fftw_complex CirculantRing<Int, d, d_fft>::fft_cache2[d_fft/2+1] __attribute__ ((aligned (32)));
template<class Int, uint64_t d, uint64_t d_fft>
fftw_complex CirculantRing<Int, d, d_fft>::fft_cache3[d_fft/2+1] __attribute__ ((aligned (32)));

typedef CirculantRing<Zqp,     1, 1> Rz;
typedef CirculantRing<Zq ,    P1> Rp1;
typedef CirculantRing<Zq ,    P2> Rp2;
typedef CirculantRing<Zqp,    P1> Rp1_p;
typedef CirculantRing<Zqp,    P2> Rp2_p;
typedef CirculantRing<Zqcrt,  P1> Rp1_crt;
typedef CirculantRing<Zqcrt,  P2> Rp2_crt;
typedef CirculantRing<Zqp,   P1*P2, FFT_DIM2> Rp12;
typedef CirculantRing<Zqcrt, P1*P2, FFT_DIM2> Rp12_crt;

#endif  // HE8_RING_H

