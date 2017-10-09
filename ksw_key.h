#ifndef HE8_KSW_KEY_H
#define HE8_KSW_KEY_H

/**
  * @file ksw_key.h
  * @date Spring 2017
  * @author Guillaume Bonnoron
  * @brief Definition of Key-Switching key class
  *
  * The KSWkey object is a convenient container for a Key-Switching for secret key s to secret key s' (sp)
  * It contains only basic initialization functions and accessors and is extensively templated for performance.
  */

#include <assert.h>
#include "rlwe.h"

/**
 * @class KSWKey
 * @brief The key-Switching key from s to sp, aka RLWE encryptions under sp of scaled s.
 * @tparam R the class representing ring elements
 * @tparam n the dimension of the source key s
 * @tparam B the integer basis for the decomposition
 * @tparam K the size of the decomposition
 */
template<class R, uint64_t n = 1, uint64_t np = 1, uint64_t B = B_Q, uint64_t K = K_Q>
class KSWKey
{
private:
    /**
     * @brief S contains the decompositions of the RLWE encryptions of each element in s under sp
     *
     * There are as many lines as the dimension of s (n)
     * Each line contains K RLWE encryptions, one for each s*B^i, 0 <= i < K
     * The cells are RLWE encryptions under key sp. They have dimension 1.
     */
    RLWE<R, np> S[n][K];
    bool init_ok;       ///< true if the initialization is done, false otherwise
public:
    /**
     * @brief Default constructor. It does not initialize the structure
     */
    KSWKey();
    /**
     * @brief Constructor that computes all the required encryptions. It initializes the structure
     * @param s the source key
     * @param sp the destination key
     */
    KSWKey(const R s[n], const R sp[np], const double variance);
    /**
     * @brief Initialization function that computes all the required encryptions.
     * @param s the source key
     * @param sp the destination key
     */
    void init(const R s[n], const R sp[np], const double variance);
    /**
     * @brief Returns the encryption of S[idx_n]*B^idx_K
     * @param idx_n the source key element
     * @param idx_K the idx_K-th power
     * @return RLWE under sp of s[idx_n] * B^idx_K
     */
    inline const RLWE<R, np> &operator() (const size_t idx_n, const size_t idx_K) const
    {
        assert(init_ok);
        return S[idx_n][idx_K];
    }
    /**
     * @brief Accessor to the basis for the decomposition
     * @return B
     */
    inline uint64_t get_basis() const
    {
        return B;
    }
    /**
     * @brief Accessor to the size of the decomposition
     * @return K
     */
    inline uint64_t get_size() const
    {
        return K;
    }
};

typedef KSWKey<Rp1>                    KSWKeyRp1;
typedef KSWKey<Rp2>                    KSWKeyRp2;
typedef KSWKey<Rp12, 3, 1, B_Qp, K_Qp> KSWKeyRp12;
typedef KSWKey<Rz, P1, N, B_Qp_2, K_Qp_2>  KSWKeyLWE;

#endif  // HE8_KSW_KEY_H

