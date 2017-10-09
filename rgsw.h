#ifndef HE8_RGSW_H
#define HE8_RGSW_H

/**
  * @file rgsw.h
  * @date Spring 2017
  * @author Guillaume Bonnoron
  * @brief Declaration of RGSW encryption
  *
  * This class defines how to encrypt and decrypt following the symmetric GSW encryption scheme over ring.
  * The object holds the ciphertext as a matrix.
  * For each template setting, the class also includes the gadget matrix
  */

#include "param.h"
#include "gadget.h"
#include <iostream>

template<class R, uint64_t n>
class RLWE;

/**
 * @class RGSW
 * @brief Store Ring-GSW ciphertexts and provide encryption/decryption routines
 * @tparam R the class for representing its elements
 */
template <class R>
class RGSW
{
private:
    static int64_t G[2*K_Q][2];     ///< the Gadget matrix for fixed class
    R A[2*K_Q][2];            ///< the ciphertext
public:
    /**
     * @brief Initialize the gadget matrix for the class in use. To be called once per setting by the user
     */
    static void G_init()
    {
        Gadget<R, B_Q, K_Q>::template gadget<2>(G);
        /*
        for (size_t j=0 ; j < 2 ; ++j)
        {
            for (size_t i=0 ; i < 2*K_Q ; ++i) std::cout << G[i][j] << std::endl;
            std::cout << std::endl;
        }
        */
    }

    /**
     * @brief Encryption function. It stores the ciphertext in the object
     * @param s the secret to use for the encryption
     * @param m the plaintext data to encrypt
     * @tparam Rt the plaintext set
     */
    template<class Rt>
    void encrypt(const R* const s, const Rt &m, const double variance);

    /**
     * @brief Decryption function. It decrypts the internal ciphertext
     * @param s the secret to use for the decryption
     * @param dec receives the decryption of the internal ciphertext
     * @tparam Rt the plaintext set
     */
    template <class Rt>
    void decrypt(const R* const s, Rt &dec) const;

    inline const R &operator() (const size_t i, const size_t j) const
    {
        return A[i][j];
    }

    friend std::ostream& operator<<(std::ostream& os, const RGSW& obj)
    {
        for (size_t i = 0 ; i < 2*K_Q ; ++i)
        {
            for (size_t j = 0 ; j < 2 ; ++j)
                os << obj.A[i][j] << " ";
            os << std::endl;
        }
        return os;
    }
};
template <class R>
int64_t RGSW<R>::G[2*K_Q][2];

typedef RGSW<Rp1> Rp1GSW;
typedef RGSW<Rp2> Rp2GSW;

#endif  // HE8_RGSW_H

