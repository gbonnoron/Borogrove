#ifndef HE8_RLWE_H
#define HE8_RLWE_H

/**
  * @file rlwe.h
  * @date Spring 2017
  * @author Guillaume Bonnoron
  * @brief Declaration of RLWE encryption and other operations
  *
  * This class defines how to encrypt and decrypt following the symmetric LWE encryption scheme over ring.
  * It can be used for LWE (ring of degree 1), Ring-LWE (ring of arbitrary degree) and n-dimensionnal variant of Ring-LWE (array of ring elements)
  * The object holds the ciphertext c=(a, b) as an array of ring elements for a (of size 1 for LWE/Ring-LWE), plus a ring element for b.
  * It is templated with the ring, the ciphetext modulus and the dimension
  */

#include <iostream>
#include "circulant_ring.h"
#include "rgsw.h"

template<class R, uint64_t n, uint64_t np, uint64_t B, uint64_t K>
class KSWKey;

/**
 * @class RLWE
 * @brief Class for encryption/decryption and other associated operations
 * @tparam R the class for the elements (=Z for LWE)
 * @tparam n the ciphertext dimension
 */
template<class R, uint64_t n>
class RLWE
{
protected:
    R *a;    ///< array of ring elements of size n
    R b;     ///< one ring element
public:
    /**
     * @brief Default constructor, doing only memory allocation
     */
    RLWE();
    /**
     * @brief Copy constructor
     * @param src RLWE object to duplicate
     */
    RLWE(const RLWE &src);
    /**
     * @brief Constructs an RLWE encryption from a representation
     * @param a an array of size n
     * @param b a ring element
     */
    RLWE(R *new_a, R new_b);
    /**
     * @brief Constructs an RLWE encryption from a representation
     * @param c an array of size n+1, b is in the last coefficient
     */
    RLWE(R *new_c, bool copy=true);
    /**
     * @brief Default destructor
     */
    ~RLWE();
    /**
     * @brief Accessor to the internal representation
     * @return a pointer to the internal array
     */
    inline R* get_a() const
    {
        return a;
    }
    /**
     * @brief Accessor to the internal representation
     * @return the ring element b
     */
    inline R get_b() const
    {
        return b;
    }

    /**
     * @brief Encryption function. It stores the ciphertext in the object
     * @param s the secret to use for the encryption
     * @param m the plaintext data to encrypt
     * @param q ciphertext modulus
     * @param t plaintext modulus
     * @param variance variance for the noise
     * @tparam Rt the plaintext set
     */
    template <class Rt>
    void encrypt(const R* const s, const Rt &m, const uint64_t q, const uint64_t t, const double variance);
    /**
     * @brief Decryption function. It decrypts the internal ciphertext
     * @param s the secret to use for the decryption
     * @param dec receive the decryption of the internal ciphertext
     * @param q ciphertext modulus
     * @param t plaintext modulus
     * @tparam Rt the plaintext set
     */
    template <class Rt>
    void decrypt(const R* const s, Rt &dec, const uint64_t q, const uint64_t t) const;

    /**
     * @brief Return the noise in the ciphertext
     * @param s the secret key
     * @param q the ciphertext modulus
     * @param t the plaintext modulus
     * @return the noise magnitude
     */
    double noise(const R* const s, const uint64_t q, const uint64_t t) const;

    /**
     * @brief Computes the Trace of each internal ring elements over another ring of degree p (templated)
     * @tparam Rp the other ring
     * @tparam p the degree of the other ring
     * @return a RLWE encryption on the ring Rp
     */
    template<class Rp, uint64_t p, uint64_t d_fft>
    RLWE<Rp, n> trace() const;

    /**
     * @brief Compute the modulo on each internal elements, makes sense for LWE only
     * @param an array of size n+1 with ring elements modulo rhs
     * @tparam Zp the class for the resulting elements
     */
    template <class Zp>
    void mod(Zp* ab) const;

    /**
     * @brief ModSwitch operation (out-of-place)
     * @tparam q the old ciphertext modulus
     * @tparam q2 the new ciphertext modulus
     * @tparam R2 the new ciphertext set
     * @return an RLWE encryption of the same plaintext under the new modulus
     */
    template<class R2, uint64_t q2, uint64_t q>
    RLWE<R2, n> mod_switch() const;

    /**
     * @brief KeySwitch operation (out-of-place)
     * @param ksw the Key-Switching key
     * @tparam B the basis for the decomposition
     * @tparam K the size of the decompositions
     * @return an RLWE encryption of the same plaintext under the new key
     */
    template <uint64_t np, uint64_t B, uint64_t K>
    RLWE<R, np> key_switch(const KSWKey <R, n, np, B, K> &ksw) const;

    /**
     * @brief ExtMult operation (out-of-place)
     * @param rgsw the RGSW ciphertext to multiply with
     * @return the RLWE encryption result of the external multiplication
     */
    RLWE ext_mult(const RGSW<R> &rgsw) const;
    /**
     * @brief ExtMult operation (in-place)
     * @param rgsw the RGSW ciphertext to multiply with
     */
    void ext_mult_inplace(const RGSW<R> &rgsw);

    /**
     * @brief Applies the Galois morphism to each internal ring elements (out-of-place)
     * @param alpha the morphism is psi_alpha : T -> T^alpha
     * @return the resulting RLWE encryption
     */
    RLWE galois(const uint64_t alpha) const;
    /**
     * @brief Applies the Galois morphism to each internal ring elements (in-place)
     * @param alpha the morphism is psi_alpha : T -> T^alpha
     */
    void galois_inplace(const uint64_t alpha);

    /**
     * @brief Fallback generic multiplication function for cases where ring elements are small nor in FFT
     * @param rhs the ring element to multiply with
     */
    inline void mult(const fftw_complex *rhs)
    {
        for (size_t i = 0 ; i < n ; ++i)
            a[i].mult(rhs);
        b.mult(rhs);
    }
    /**
     * @brief Fallback generic multiplication function for cases where ring elements are small nor in FFT
     * @param rhs the ring element to multiply with
     */
    inline void mult(const R &rhs)
    {
        for (size_t i = 0 ; i < n ; ++i)
            a[i].mult(rhs);
        b.mult(rhs);
    }

    /**
     * @brief Compute FFT from data for each element in internal representation
     */
    inline void decomp_fft()
    {
        for (size_t i = 0 ; i < n ; ++i)
            a[i].decomp_fft();
        b.decomp_fft();
    }
    /**
     * @brief Compute data from fft for each element in internal representation
     */
    inline void recomp_fft()
    {
        for (size_t i = 0 ; i < n ; ++i)
            a[i].recomp_fft(SPLIT_Q);
        b.recomp_fft(SPLIT_Q);
    }

    inline RLWE& operator+=(const RLWE& rhs)
    {
        for (size_t i = 0 ; i < n ; ++i)
            a[i] += rhs.a[i];
        b += rhs.b;
        return *this;
    }
    inline RLWE& operator-=(const RLWE& rhs)
    {
        for (size_t i = 0 ; i < n ; ++i)
            a[i] -= rhs.a[i];
        b -= rhs.b;
        return *this;
    }
    inline RLWE& operator*=(const R& rhs)
    {
        for (size_t i = 0 ; i < n ; ++i)
            a[i] *= rhs;
        b *= rhs;
        return *this;
    }
    inline RLWE& operator*=(const int64_t& rhs)
    {
        for (size_t i = 0 ; i < n ; ++i)
            a[i] *= rhs;
        b *= rhs;
        return *this;
    }
    RLWE& operator=(const RLWE& other);
    RLWE& operator=(RLWE&& other) noexcept;
    friend inline RLWE operator+(RLWE lhs, const RLWE& rhs) { lhs += rhs; return lhs; }
    friend inline RLWE operator-(RLWE lhs, const RLWE& rhs) { lhs -= rhs; return lhs; }
    friend inline RLWE operator*(RLWE lhs, const int64_t &rhs) { lhs *= rhs; return lhs; }
    friend inline RLWE operator*(RLWE lhs, const R &rhs) { lhs *= rhs; return lhs; }
    friend inline RLWE operator*(const R &lhs, RLWE rhs) { rhs *= lhs; return rhs; }

    friend std::ostream& operator<<(std::ostream& os, const RLWE& obj)
    {
        os << "(";
        for (size_t i = 0 ; i < n ; ++i)
            os << obj.a[i] << ", ";
        os << obj.b << ")";
        return os;
    }
};

template <class R>
class RLWE1 : public RLWE<R, 1>
{
private:
    using RLWE<R, 1>::a;
    using RLWE<R, 1>::b;
public:
    using RLWE<R, 1>::ext_mult_inplace;
    using RLWE<R, 1>::galois_inplace;
    using RLWE<R, 1>::RLWE;
    /**
     * @brief KeySwitch operation (in-place). The new and old keys have the same dimension
     * @param ksw the Key-Switching key
     * @tparam B the basis for the decomposition
     * @tparam K the size of the decompositions
     */
    template <uint64_t B, uint64_t K>
    void key_switch_inplace(const KSWKey <R, 1, 1, B, K> &ksw);
    /**
     * @brief ExtExpMultAdd operation (in-place)
     * @param rgsw the RGSW encryption to multiply and add
     * @param alpha the coefficient for the multiplication
     * @param S_alpha a Key-Switching key from psi_alpha(s) to s
     * @param beta inverse of alpha
     * @param S_beta a Key-Switching key from psi_beta(s) to s
     * @tparam B the basis for the decomposition in the Key-switching keys
     * @tparam K the size of the decompositions
     */
    void ext_exp_mult_add(const RGSW<R> &rgsw, const uint64_t alpha, const KSWKey<R, 1, 1, B_Q, K_Q>& S_alpha, const uint64_t beta, const KSWKey<R, 1, 1, B_Q, K_Q>& S_beta);
};

typedef RLWE<Rz  , P1>     LWE;
typedef RLWE1<Rp1>  Rp1LWE;
typedef RLWE1<Rp2>  Rp2LWE;
typedef RLWE<Rp12, 3> Rp12LWE;


#endif  // HE8_RLWE_H

