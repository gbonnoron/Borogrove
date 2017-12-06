/**
  * @file circulant_ring.cpp
  * @date Spring 2017
  * @author Guillaume Bonnoron
  * @brief Definition of the CirculantRing class
  */

#include <fftw3.h>
#include <utility>
#include <cmath>
#include <random>
#include <limits>
#include <cstring>
#include <algorithm>
#include "circulant_ring.h"
#include "operations.h"
#include "predicate.h"

template<class Int, uint64_t d, uint64_t d_fft>
CirculantRing<Int, d, d_fft>::CirculantRing() : msb_fft(nullptr), lsb_fft(nullptr), fft_status(false)
{
    x = align_alloc<Int>(ALIGNMENT, d);
}
template<class Int, uint64_t d, uint64_t d_fft>
CirculantRing<Int, d, d_fft>::CirculantRing (const CirculantRing &src)
{
    fft_status = src.fft_status;
    if (fft_status)
    {
        x = nullptr;
        msb_fft = align_alloc<fftw_complex>(ALIGNMENT, d_fft/2+1);
        lsb_fft = align_alloc<fftw_complex>(ALIGNMENT, d_fft/2+1);
        if (d < 2000)
        {
            std::memcpy(msb_fft, src.msb_fft, (d_fft/2+1) * 2 * sizeof(double));
            std::memcpy(lsb_fft, src.lsb_fft, (d_fft/2+1) * 2 * sizeof(double));
        }
        else
        {
            for (size_t i = 0 ; i < d_fft/2 + 1 ; ++i)
            {
                msb_fft[i][0] = src.msb_fft[i][0];
                msb_fft[i][1] = src.msb_fft[i][1];
                lsb_fft[i][0] = src.lsb_fft[i][0];
                lsb_fft[i][1] = src.lsb_fft[i][1];
            }
        }
    }
    else
    {
        x = align_alloc<Int>(ALIGNMENT, d);
        msb_fft = nullptr;
        lsb_fft = nullptr;
        std::copy(src.x, src.x + d, x);
    }
}
template<class Int, uint64_t d, uint64_t d_fft>
template <class Int2>
CirculantRing<Int, d, d_fft>::CirculantRing (const CirculantRing<Int2, d, d_fft> &src)
{
    assert(! src.get_fft_status());
    x = align_alloc<Int>(ALIGNMENT, d);
    Int2 *src_x = src.get_data();
    for (size_t i = 0 ; i < d ; ++i)
        x[i] = src_x[i];
    fft_status = false;
    msb_fft = nullptr;
    lsb_fft = nullptr;
}

template<class Int, uint64_t d, uint64_t d_fft>
template <uint64_t d_fft2>
CirculantRing<Int, d, d_fft>::CirculantRing (const CirculantRing<Int, d, d_fft2> &src)
{
    assert(! src.get_fft_status());
    const Int *src_x = src.get_data();
    x = align_alloc<Int>(ALIGNMENT, d);
    std::copy(src_x, src_x + d, x);
    fft_status = false;
    msb_fft = nullptr;
    lsb_fft = nullptr;
}

template<class Int, uint64_t d, uint64_t d_fft>
CirculantRing<Int, d, d_fft>::CirculantRing (Int *values, bool copy)
{
    if (copy)
    {
        x = align_alloc<Int>(ALIGNMENT, d);
        std::copy(values, values + d, x);
    }
    else
        x = std::move(values);
    fft_status = false;
    msb_fft = nullptr;
    lsb_fft = nullptr;
}

template<class Int, uint64_t d, uint64_t d_fft>
CirculantRing<Int, d, d_fft>::CirculantRing (fftw_complex *msb_fft_values, fftw_complex *lsb_fft_values, bool copy)
{
    if (copy)
    {
        msb_fft = align_alloc<fftw_complex>(ALIGNMENT, d_fft/2+1);
        lsb_fft = align_alloc<fftw_complex>(ALIGNMENT, d_fft/2+1);
        std::memcpy(msb_fft, msb_fft_values, (d_fft/2+1) * 2 * sizeof(double));
        std::memcpy(lsb_fft, lsb_fft_values, (d_fft/2+1) * 2 * sizeof(double));
    }
    else
    {
        msb_fft = std::move(msb_fft_values);
        lsb_fft = std::move(lsb_fft_values);
    }
    fft_status = true;
    x = nullptr;
}

template<class Int, uint64_t d, uint64_t d_fft>
CirculantRing<Int, d, d_fft>::CirculantRing (const Int val)
{
    x = align_alloc<Int>(ALIGNMENT, d);
    x[0] = val;
    std::fill(x + 1, x + d, 0);
    fft_status = false;
    msb_fft = nullptr;
    lsb_fft = nullptr;
}
template<class Int, uint64_t d, uint64_t d_fft>
CirculantRing<Int, d, d_fft>::CirculantRing (const int64_t *values)
{
    x = align_alloc<Int>(ALIGNMENT, d);
    for (size_t i = 0 ; i < d ; ++i)
        x[i] = values[i];
    fft_status = false;
    msb_fft = nullptr;
    lsb_fft = nullptr;
}
template<class Int, uint64_t d, uint64_t d_fft>
CirculantRing<Int, d, d_fft>::CirculantRing (const Fun &F)
{
    x = align_alloc<Int>(ALIGNMENT, d);
    x[0] = F(0);
    for (size_t i = 1 ; i < d ; ++i)
        x[d - i] = F(i);
    fft_status = false;
    msb_fft = nullptr;
    lsb_fft = nullptr;
}
template<class Int, uint64_t d, uint64_t d_fft>
CirculantRing<Int, d, d_fft>::~CirculantRing()
{
    if (fft_status)
    {
        free(msb_fft);
        free(lsb_fft);
        msb_fft = nullptr;
        lsb_fft = nullptr;
    }
    else
    {
        free(x);
        x = nullptr;
    }
}

template<class Int, uint64_t d, uint64_t d_fft>
CirculantRing<Int, d, d_fft>& CirculantRing<Int, d, d_fft>::operator=(const int &other)
{
    assert(! fft_status);
    x[0] = Int(other);
    std::fill(x + 1, x + d, 0);
    fft_status = false;
    msb_fft = nullptr;
    lsb_fft = nullptr;
    return *this;
}

template<class Int, uint64_t d, uint64_t d_fft>
CirculantRing<Int, d, d_fft>& CirculantRing<Int, d, d_fft>::operator=(const CirculantRing& other)
{
    if (this != &other)
    {
        if (other.fft_status)
        {
            if (! fft_status)
            {
                free(x);
                msb_fft = align_alloc<fftw_complex>(ALIGNMENT, d_fft/2+1);
                lsb_fft = align_alloc<fftw_complex>(ALIGNMENT, d_fft/2+1);
                fft_status = true;
                x = nullptr;
            }
            std::memcpy(msb_fft, other.msb_fft, (d_fft/2+1) * 2 * sizeof(double));
            std::memcpy(lsb_fft, other.lsb_fft, (d_fft/2+1) * 2 * sizeof(double));
        }
        else
        {
            if (fft_status)
            {
                free(msb_fft);
                free(lsb_fft);
                x = align_alloc<Int>(ALIGNMENT, d);
                fft_status = false;
                msb_fft = nullptr;
                lsb_fft = nullptr;
            }
            std::copy(other.x, other.x + d, x);
        }
    }
    return *this;
}
template<class Int, uint64_t d, uint64_t d_fft>
CirculantRing<Int, d, d_fft>& CirculantRing<Int, d, d_fft>::operator=(CirculantRing&& other) noexcept
{
    if (this != &other)
    {
        if (fft_status)
        {
            free(msb_fft);
            free(lsb_fft);
            msb_fft = nullptr;
            lsb_fft = nullptr;
        }
        else
        {
            free(x);
            x = nullptr;
        }

        if (other.fft_status)
        {
            //C++14 : msb_fft = std::exchange(other.msb_fft, nullptr);
            fftw_complex *old_msb = std::move(other.msb_fft);
            other.msb_fft = std::forward<fftw_complex *>(nullptr);
            msb_fft = old_msb;
            //C++14 : lsb_fft = std::exchange(other.lsb_fft, nullptr);
            fftw_complex *old_lsb = std::move(other.lsb_fft);
            other.lsb_fft = std::forward<fftw_complex *>(nullptr);
            lsb_fft = old_lsb;
        }
        else
        {
            //C++14 : x = std::exchange(other.x, nullptr);
            Int *old_value = std::move(other.x);
            other.x = std::forward<Int *>(nullptr);
            x = old_value;
        }
        fft_status = other.fft_status;
    }
    return *this;
}

template<class Int, uint64_t d, uint64_t d_fft>
void CirculantRing<Int, d, d_fft>::mult(const fftw_complex *f)
{
    const int64_t B = Int::get_split_threshold();
    CirculantRing lsb, msb;
    split(msb, lsb);

    //a1*b0
    msb.compute_fft(fft_cache1);
    fft::product<d_fft>(fft_cache2, fft_cache1, f);

    compute_fft_inv(fft_cache2);
    for (size_t i = 0 ; i < d ; ++i)
        x[i] = ((int64_t)x[i] % B) * B;

    //a0*b0
    lsb.compute_fft(fft_cache1);
    fft::product<d_fft>(fft_cache2, fft_cache1, f);
    compute_fft_inv(fft_cache2, tmp);
    for (size_t i = 0 ; i < d ; ++i)
        x[i] += tmp[i];
}

template<class Int, uint64_t d, uint64_t d_fft>
void CirculantRing<Int, d, d_fft>::mult(const CirculantRing& rhs)
{
    rhs.compute_fft(fft_cache3);

    const int64_t B = Int::get_split_threshold();
    CirculantRing lsb, msb;
    split(msb, lsb);

    //a1*b0
    msb.compute_fft(fft_cache1);
    fft::product<d_fft>(fft_cache2, fft_cache1, fft_cache3);

    compute_fft_inv(fft_cache2);
    shift(B);

    //a0*b0
    lsb.compute_fft(fft_cache1);
    fft::product<d_fft>(fft_cache2, fft_cache1, fft_cache3);
    compute_fft_inv(fft_cache2, tmp);
    for (size_t i = 0 ; i < d ; ++i)
        x[i] += tmp[i];
}

template<class Int, uint64_t d, uint64_t d_fft>
CirculantRing<Int, d, d_fft>& CirculantRing<Int, d, d_fft>::operator*=(const CirculantRing& rhs)
{
    assert (rhs.fft_status);
    const int64_t B = Int::get_split_threshold();
    compute_fft(fft_cache1);

    fft::product<d_fft>(fft_cache2, fft_cache1, rhs.get_msb_fft());

    compute_fft_inv(fft_cache2);
    shift(B);

    fft::product<d_fft>(fft_cache2, fft_cache1, rhs.get_lsb_fft());
    compute_fft_inv(fft_cache2, tmp);
    for (size_t i = 0 ; i < d ; ++i)
        x[i] += tmp[i];

    return *this;
}

template<class Int, uint64_t d, uint64_t d_fft>
void CirculantRing<Int, d, d_fft>::compute_fft(fftw_complex *x_fft) const
{
    double *in = fft.get_in();
    for (size_t i = 0 ; i < d ; ++i)
        in[i] = (double) x[i];
    std::fill(in + d, in + d_fft, 0.0);
    fft.forward(x_fft);
}

template<class Int, uint64_t d, uint64_t d_fft>
void CirculantRing<Int, d, d_fft>::compute_fft_inv(fftw_complex *x_fft, Int *values) const
{
    //std::cout.precision(std::numeric_limits<double>::max_digits10);
    double *in = fft.get_in();
    fft.backward(x_fft);

    if (d > 1)
    {
        for (size_t i = 0 ; i < d ; ++i)
        {
            values[i] = (int64_t)std::floor(in[i] + 0.5) + (int64_t)std::floor(in[i+d]+0.5);
            //values[i] = (int64_t)(in[i] + (double)0.5) + (int64_t)(in[i + d] + (double)0.5); //TODO faster
            //std::cout << std::fixed << in[i] << " + " << std::fixed << in[i+d] << " = " << x[i] << std::endl;
        }
    }
    else
    {
        values[0] = (int64_t)std::floor(in[0] + 0.5);
    }
}

template<class Int, uint64_t d, uint64_t d_fft>
void CirculantRing<Int, d, d_fft>::decomp_fft()
{
    assert(! fft_status);
    CirculantRing msb, lsb;
    split(msb, lsb);

    msb_fft = align_alloc<fftw_complex>(ALIGNMENT, d_fft/2+1);
    lsb_fft = align_alloc<fftw_complex>(ALIGNMENT, d_fft/2+1);
    msb.compute_fft(msb_fft);
    lsb.compute_fft(lsb_fft);
    free(x);
    x = nullptr;
    fft_status = true;
}

template<class Int, uint64_t d, uint64_t d_fft>
void CirculantRing<Int, d, d_fft>::recomp_fft(const int64_t B)
{
    assert(fft_status);
    Int *values = align_alloc<Int>(ALIGNMENT, d);
    compute_fft_inv(msb_fft, values);

    x = align_alloc<Int>(ALIGNMENT, d);
    for (size_t i = 0 ; i < d ; ++i)
        x[i] = (values[i] % B) * B;
    //std::cout << "recomp msb " << *this << std::endl;

    compute_fft_inv(lsb_fft, values);
    for (size_t i = 0 ; i < d ; ++i)
        x[i] += values[i];
    //std::cout << "recomp lsb " << *this << std::endl;

    free(msb_fft);
    free(lsb_fft);
    free(values);
    fft_status = false;
}

template<class Int, uint64_t d, uint64_t d_fft>
int64_t CirculantRing<Int, d, d_fft>::split(CirculantRing &msb, CirculantRing &lsb) const
{
    Int *most_significant = align_alloc<Int>(ALIGNMENT, d);
    Int *least_significant = align_alloc<Int>(ALIGNMENT, d);
    int64_t used_threshold = x[0].split(most_significant[0], least_significant[0]);
    for (size_t i = 1 ; i < d ; ++i)
        x[i].split(most_significant[i], least_significant[i]);
    msb = CirculantRing(most_significant, false);
    lsb = CirculantRing(least_significant, false);

    return used_threshold;
}

template<class Int, uint64_t d, uint64_t d_fft>
void CirculantRing<Int, d, d_fft>::unsplit(CirculantRing &msb, CirculantRing &lsb)
{
    for (size_t i = 0 ; i < d ; ++i)
        x[i].unsplit(msb.x[i], lsb.x[i]);
}

/*
inline static int64_t bernouilli(const long double p)
{
    return (rand() % 1000 < (p - std::floor(p)) * 1000) ? 1 : 0;
}
*/
template<class Int, uint64_t d, uint64_t d_fft>
template<class R2>
R2 CirculantRing<Int, d, d_fft>::rounding(const uint64_t old_modulus, const uint64_t new_modulus) const
{
    int64_t *res = align_alloc<int64_t>(ALIGNMENT, d);
    for (size_t i = 0 ; i < d ; ++i)
    {
        res[i] = (int64_t)std::floor((double)x[i] * (double)new_modulus / (double) old_modulus + 0.5);
        //std::cout << x[i] << " * " << new_modulus << " / " << old_modulus << " = " << quotient << " ";
        //std::cout << res[i] << " " << (int64_t)((int64_t)(quotient) + bernouilli(quotient)) << std::endl;
    }
    R2 result(res);
    //std::cout << result << std::endl;
    free(res);
    return result;
}
template<class Int, uint64_t d, uint64_t d_fft>
CirculantRing<Int, d, d_fft> CirculantRing<Int, d, d_fft>::exact_rounding(const uint64_t old_modulus, const uint64_t new_modulus) const
{
    Int *res = align_alloc<Int>(ALIGNMENT, d);
    for (size_t i = 0 ; i < d ; ++i)
        res[i] = (int64_t)std::floor((double)x[i] * (double)new_modulus / (double) old_modulus + 0.5);
    return CirculantRing<Int, d, d_fft>(res, false);
}


template<class Int>
static Int sample_long(const double variance)
{
    double r = 0;
    for (size_t i = 0 ; i < 12 ; ++i)
        r += static_cast<double> (rand()) / static_cast<double> (RAND_MAX) - 0.5;
    return (int64_t)floor(r * std::sqrt(variance) + 0.5);
}

template<class Int, uint64_t d, uint64_t d_fft>
CirculantRing<Int, d, d_fft> CirculantRing<Int, d, d_fft>::gaussian_sample(const double variance)
{
    Int *res = align_alloc<Int>(ALIGNMENT, d);
    for (size_t i = 0 ; i < d ; ++i)
        res[i] = sample_long<Int>(variance);
    return CirculantRing<Int, d, d_fft>(res, false);
}


template<class Int, uint64_t d, uint64_t d_fft>
CirculantRing<Int, d, d_fft> CirculantRing<Int, d, d_fft>::uniform_sample()
{
    Int *res = align_alloc<Int>(ALIGNMENT, d);
    for (size_t i = 0 ; i < d ; ++i)
        //res[i] = rand();
        res[i] = (((int64_t)rand())<<31) + (int64_t)rand();
    return CirculantRing<Int, d, d_fft>(res, false);
}

template<class Int, uint64_t d, uint64_t d_fft>
CirculantRing<Int, d, d_fft> CirculantRing<Int, d, d_fft>::binary_sample()
{
    Int *res = align_alloc<Int>(ALIGNMENT, d);
    for (size_t i = 0 ; i < d ; ++i)
        res[i] = rand() % 2;
    return CirculantRing<Int, d, d_fft>(res, false);
}

template<class Int, uint64_t d, uint64_t d_fft>
CirculantRing<Int, d, d_fft> CirculantRing<Int, d, d_fft>::sample_a()
{
    Int *res = align_alloc<Int>(ALIGNMENT, d);
    if (d == 1)
    {
        res[0] = (((int64_t)rand())<<31) + (int64_t)rand();
    }
    else
    {
        Int sum;
        sum = 0;
        for (size_t i = 0 ; i < d-1 ; ++i)
        {
            res[i] = (((int64_t)rand())<<31) + (int64_t)rand();
            sum += res[i];
        }
        res[d-1] = -sum;
    }
    return CirculantRing<Int, d, d_fft>(res, false);
}

template<class Int, uint64_t d, uint64_t d_fft>
CirculantRing<Int, d, d_fft> CirculantRing<Int, d, d_fft>::sample_s(const double density)
{
    Int *res = align_alloc<Int>(ALIGNMENT, d);
    if (d == 1)
    {
        res[0] = rand() % 3 - 1;
    }
    else
    {
        std::fill(res, res + d, 0);

        size_t position, defined;
        defined = 0;
        while (defined < d*density)
        {
            position = rand() % d;
            if (res[position] == 0)
            {
                res[position] = 1;
                ++defined;
            }
        }
        defined = 0;
        while (defined < d*density)
        {
            position = rand() % d;
            if (res[position] == 0)
            {
                res[position] = -1;
                ++defined;
            }
        }
    }
    /*
    for (size_t i = 0 ; i < d ; ++i)
        std::cout << res[i] << " ";
    std::cout << std::endl;
    */
    return CirculantRing<Int, d, d_fft>(res, false);
}

template<class Int, uint64_t d, uint64_t d_fft>
CirculantRing<Int, d, d_fft> CirculantRing<Int, d, d_fft>::sample_e(const double variance)
{
    if (variance > 1e6)
        return gaussian_sample(variance);
    else
    {
        Int *res = align_alloc<Int>(ALIGNMENT, d);
        if (d == 1)
        {
            res[0] = rand() % 3 - 1;
        }
        else
        {
            std::fill(res, res + d, 0);

            for (size_t i = 0 ; i < variance*d/2 ; ++i)
            {
                res[rand() % d] += 1;
                res[rand() % d] -= 1;
            }
        }
        return CirculantRing<Int, d, d_fft>(res, false);
    }
}

template<class Int, uint64_t d, uint64_t d_fft>
CirculantRing<Int, d, d_fft> CirculantRing<Int, d, d_fft>::galois(const uint64_t alpha) const
{
    CirculantRing<Int, d, d_fft> result(*this);
    result.galois_inplace(alpha);
    return result;
}

template<class Int, uint64_t d, uint64_t d_fft>
void CirculantRing<Int, d, d_fft>::galois_inplace(const uint64_t alpha)
{
    Int *tmp = align_alloc<Int>(ALIGNMENT, d);
    size_t current = 0;
    for (size_t i = 0 ; i < d ; ++i)
    {
        tmp[current] = x[i];
        current += alpha;
        current %= d;
    }
    std::copy(tmp, tmp + d, x);
    free(tmp);
}

template<class Int, uint64_t d, uint64_t d_fft>
template<uint64_t d2, uint64_t d2_fft, uint64_t d3_fft>
CirculantRing<Int, d*d2, d3_fft> CirculantRing<Int, d, d_fft>::tensor_product(const CirculantRing<Int, d2, d2_fft>& rhs) const
{
    const uint64_t new_dim = d*d2;
    Int *result = align_alloc<Int>(ALIGNMENT, d*d2);
    Int *rhs_x = rhs.get_data();
    size_t current = 0;
    size_t current_i = 0;
    for (size_t i = 0 ; i < d ; ++i)
    {
        current = current_i;
        for (size_t j = 0 ; j < d2 ; ++j)
        {
            result[current] = x[i] * rhs_x[j];
            current += d;
            current %= new_dim;
        }
        current_i += d2;
        current_i %= new_dim;
    }
    return CirculantRing<Int, d*d2, d3_fft>(result, false);
}

template<class Int, uint64_t pq, uint64_t d_fft>
template <uint64_t p, uint64_t d2_fft>
CirculantRing<Int, p, d2_fft> CirculantRing<Int, pq, d_fft>::trace() const
{
    const uint64_t q = pq/p;
    Int *x_res = align_alloc<Int>(ALIGNMENT, p);

    size_t idx = 0;
    for (size_t i = 0 ; i < p ; ++i)
    {
        x_res[i] = x[idx];
        idx += q;
    }
    return CirculantRing<Int, p, d2_fft>(x_res, false);
}

template class CirculantRing<Zt   ,  1, 1>;
template class CirculantRing<Zt   , P1>;
template class CirculantRing<Zt   , P2>;
template class CirculantRing<Zqp  ,  1, 1>;
template class CirculantRing<Zp12 ,  1, 1>;
template class CirculantRing<Zq   , P1>;
template class CirculantRing<Zq   , P2>;
template class CirculantRing<Zqcrt, P1>;
template class CirculantRing<Zqcrt, P2>;
template class CirculantRing<Zqcrt, P1*P2, FFT_DIM2>;
template class CirculantRing<Zqp  , P1*P2, FFT_DIM2>;

template Rp1::CirculantRing<Zt>(const CirculantRing<Zt, P1> &src);
template Rp2::CirculantRing<Zt>(const CirculantRing<Zt, P2> &src);

template CirculantRing<Zt, P1>::CirculantRing<Zq>(const Rp1 &src);
template CirculantRing<Zt, P2>::CirculantRing<Zq>(const Rp2 &src);
template CirculantRing<Zt, 1, 1>::CirculantRing<Zqp>(const Rz &src);

template Rp12_crt Rp1_crt::tensor_product<P2, FFT_DIM1, FFT_DIM2>(const Rp2_crt &rhs) const;
template Rp12 Rp1_p::tensor_product<P2, FFT_DIM1, FFT_DIM2>(const Rp2_p &rhs) const;

template CirculantRing<Zp12, 1, 1> Rz::rounding<CirculantRing<Zp12, 1, 1> >(const uint64_t old_modulus, const uint64_t new_modulus) const;
template Rp1_crt Rp1::rounding<Rp1_crt>(const uint64_t old_modulus, const uint64_t new_modulus) const;
template Rp2_crt Rp2::rounding<Rp2_crt>(const uint64_t old_modulus, const uint64_t new_modulus) const;
template Rp12 Rp12_crt::rounding<Rp12>(const uint64_t old_modulus, const uint64_t new_modulus) const;

template Rp1_p Rp12::trace<P1, FFT_DIM1>() const;

// For tests only
template class CirculantRing<Zt, P1*P2, FFT_DIM2>;
template            CirculantRing<Zt, 1, 1>::CirculantRing<Zp12>(const CirculantRing<Zp12, 1, 1> &src);
template CirculantRing<Zt, P1*P2, FFT_DIM2>::CirculantRing<Zqp>(const Rp12 &src);
template          CirculantRing<Zp12, 1, 1>::CirculantRing<Zqp>(const Rz &src);

// For stats only
template Rp1_crt::CirculantRing<Zq>(const Rp1 &src);
template Rp12_crt::CirculantRing<Zqp>(const Rp12 &src);
template Rp1::CirculantRing<Zqcrt>(const Rp1_crt &src);
template Rp2::CirculantRing<Zqcrt>(const Rp2_crt &src);
template Rp12::CirculantRing<Zqcrt>(const Rp12_crt &src);
