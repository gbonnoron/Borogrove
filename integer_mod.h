#ifndef HE8_INTEGER_MOD_H
#define HE8_INTEGER_MOD_H

/**
  * @file integer_mod.h
  * @date Spring 2017
  * @author Guillaume Bonnoron
  * @brief Declaration of the IntegerMod class
  *
  * The class IntegerMod purpose is to represent element in Z/qZ where q is a template parameter. It offers the usual operators in order for arithmetic operations to be valid.
  * The internal representation is on unsigned type in the range [0, q-1].
  * Attention should be raise concerning the operations *, / and % which "lift" the element in Z/aZ back to Z to do the computation. The elements in ]q/2, q-] behave as their negative representation.
  * For example, IntegerMod<5>(4) * 2 = -1 * 2 = -2. It can easily be fed back to IntegerMod<5> and the result is what you expect it to be if the computation had been done in Z/5Z all along.
  */


#include <cmath>
#include "assert.h"
#include "param.h"

/**
 * @class IntegerMod
 * @brief Class to handle elements in Z/qZ
 * @tparam q the modulus
 */
template <int64_t q>
class IntegerMod
{
protected:
    int64_t value;
public:
    /**
     * @brief Default constructor
     */
    IntegerMod() {  }
    /**
     * @brief Copy constructor for the same modulus
     * @param src the element to duplicate
     */
    IntegerMod(const IntegerMod &src)
    {
        value = src.value;
    }
    /**
     * @brief Copy constructor for a different modulus
     * @param src the element to reinterpret and duplicate
     * @tparam q2 the modulus of the src element
     */
    template<int64_t q2>
    IntegerMod(const IntegerMod<q2> &src)
    {
        value = src.get_value() % q;
        value -= (value / (q/2+1)) * q;
    }
    /**
     * @brief Constructor from unsigned integer type
     * @param src the element to interpret and store
     */
    IntegerMod(const uint64_t &src)
    {
        value = src % q;
        value -= (value / (q/2+1)) * q;
    }
    /**
     * @brief Constructor from signed integer type
     * @param src the element to interpret and store
     */
    IntegerMod(const int64_t &src)
    {
        value = src % q;
        value -= (value / (q/2+1)) * q;
    }
    /**
     * @brief Another constructor from signed integer type
     * @param src the element to interpret and store
     */
    IntegerMod(const int &src)
    {
        value = src % q;
        value -= (value / (q/2+1)) * q;
    }
    /**
     * @brief Accessor to the internal representation
     * @return the internal value
     */
    inline int64_t get_value() const
    {
        return value;
    }
    /**
     * @brief Accessor to the modulus
     * @return the modulus
     */
    inline static int64_t get_mod()
    {
        return q;
    }
    /**
     * @brief Computes the absolute value of the element
     * @return the absolute value
     */
    inline int64_t abs() const
    {
        if (value < 0)
            return -value;
        else
            return value;
    }
    /**
     * @brief Accessor to the split threshold
     * @return the threshold
     */
    inline static int64_t get_split_threshold()
    {
        return (1<<24);
    }
    /**
     * @brief Split the element into two halves, one with the most significant part and the other with the least significant part
     * @param msb the most significant part
     * @param lsb the least signigicant part
     * @return the threshold value between the part
     */
    inline int64_t split(IntegerMod &msb, IntegerMod &lsb) const
    {
        msb = 0;
        lsb = value;
        //std::cout << value << " = " << msb << " * " << B << " + " << lsb << std::endl;
        return get_split_threshold();
    }
    /**
     * @brief Recompose the element from two halves
     * @param msb the most significant part
     * @param lsb the least signigicant part
     * @return the threshold value between the part
     */
    inline void unsplit(IntegerMod &msb, IntegerMod &lsb)
    {
        assert(msb.get_value() == 0);
        value = lsb.value;
        //std::cout << value << " = " << msb << " * " << B << " + " << lsb << std::endl;
    }

    /**
     * @brief Balances the coefficients between -B/2 and B/2
     * @param B the range
     */
    inline void balance(const int64_t B)
    {
        value -= (value / (B/2)) * B;
        value %= B;
    }
    inline IntegerMod& operator=(const IntegerMod& other)
    {
        value = other.value;
        return *this;
    }
    inline IntegerMod& operator+=(const IntegerMod& rhs)
    {
        value += rhs.value;
        value -= (value / (q/2+1)) * q;
        return *this;
    }
    inline IntegerMod& operator++()
    {
        *this += 1;
        return *this;
    }
    inline IntegerMod& operator--()
    {
        *this -= 1;
        return *this;
    }
    inline IntegerMod& operator-=(const IntegerMod& rhs)
    {
        *this += -rhs;
        return *this;
    }
    IntegerMod& operator*=(const IntegerMod& rhs)
    {
        assert(q != Q && q != Qp && q != Qcrt);
        value = (value * rhs.value) % q;

        value -= (value / (q/2+1)) * q;

        return *this;
    }
    inline int64_t operator%(const int64_t &rhs)
    {
        return value % rhs;
    }

    inline int64_t operator*(const int64_t &rhs)
    {
        return value * rhs;
    }
    friend inline int64_t operator/(const IntegerMod &lhs, const int64_t& rhs)
    {
        return lhs.value / rhs;
    }
    friend inline IntegerMod operator+(IntegerMod lhs, const IntegerMod& rhs) { lhs += rhs; return lhs; }
    friend inline IntegerMod operator-(IntegerMod lhs, const IntegerMod& rhs) { lhs -= rhs; return lhs; }
    friend inline IntegerMod operator*(IntegerMod lhs, const IntegerMod& rhs) { lhs *= rhs; return lhs; }
    //friend uint64_t operator*(const uint64_t &lhs, const IntegerMod& rhs) { return rhs * lhs; }
    friend inline bool operator!=(const IntegerMod& lhs, const IntegerMod& rhs) { return !(lhs == rhs); }
    friend inline IntegerMod operator-(IntegerMod lhs) { return -lhs.value; }

    friend inline bool operator==(const IntegerMod& lhs, const IntegerMod& rhs)
    {
        return lhs.value == rhs.value;
    }
    friend std::ostream& operator<<(std::ostream& os, const IntegerMod<q>& obj)
    {
        os << obj.value << std::flush;
        return os;
    }
    inline explicit operator double() const
    {
        return (double)value;
    }
    inline explicit operator int64_t() const
    {
        return value;
    }
    inline explicit operator size_t() const
    {
        if (value < 0)
            return value + q;
        else
            return value;
    }
};// __attribute__ ((__aligned__(ALIGNMENT)));

template <int64_t q>
class IntegerPrimeMod : public IntegerMod<q>
{
private:
    using IntegerMod<q>::value;
    using IntegerMod<q>::IntegerMod;
public:
    /**
     * @brief Default constructor
     */
    IntegerPrimeMod() {  }
    IntegerPrimeMod(const IntegerMod<q> &src)
    {
        value = (int64_t)src;
    }
    /**
     * @brief Computes the inverse, very stupidly of the element
     * @return the inverse
     */
    inline IntegerPrimeMod inv() const
    {
        if (q == P1)
            return inv_mod_P1[(size_t)(*this)];
        if (q == P2)
            return inv_mod_P2[(size_t)(*this)];
        //TODO problem here
        return 0;
    }
};// __attribute__ ((__aligned__(ALIGNMENT)));

template <int64_t q>
class IntegerBigMod : public IntegerMod<q>
{
protected:
    static const int64_t B = (q == Q) ? SPLIT_Q : (q == Qp) ? SPLIT_Qp : SPLIT_Qcrt;
    using IntegerMod<q>::value;
    using IntegerMod<q>::IntegerMod;
public:
    using IntegerMod<q>::operator*;
    /**
     * @brief Default constructor
     */
    IntegerBigMod() {  }
    IntegerBigMod(const IntegerMod<q> &src)
    {
        value = (int64_t)src;
    }
    /**
     * @brief Accessor to the split threshold
     * @return the threshold
     */
    inline static int64_t get_split_threshold()
    {
        assert(q == Qcrt);
        return B;
    }
    /**
     * @brief Split the element into two halves, one with the most significant part and the other with the least significant part
     * @param msb the most significant part
     * @param lsb the least signigicant part
     * @return the threshold value between the part
     */
    inline int64_t split(IntegerBigMod &msb, IntegerBigMod &lsb) const
    {
        assert(q == Qcrt);
        msb = value / B;
        lsb = value % B;
        //std::cout << value << " = " << msb << " * " << B << " + " << lsb << std::endl;
        return B;
    }
    /**
     * @brief Recompose the element from two halves
     * @param msb the most significant part
     * @param lsb the least signigicant part
     * @return the threshold value between the part
     */
    inline void unsplit(IntegerBigMod &msb, IntegerBigMod &lsb)
    {
        assert(q == Qcrt);
        value = msb.get_value() * B + lsb.get_value();
        //std::cout << value << " = " << msb << " * " << B << " + " << lsb << std::endl;
    }

    IntegerBigMod& operator*=(const IntegerBigMod& rhs)
    {
        assert(q == Qcrt);
        const int64_t lhs_m = value / B;
        const int64_t lhs_l = value % B;
        const int64_t rhs_m = rhs.value / B;
        const int64_t rhs_l = rhs.value % B;

        value = ((lhs_l * rhs_m + lhs_m * rhs_l) % B) * B;
        value += lhs_l * rhs_l;

        value -= (value / (q/2+1)) * q;

        return *this;
    }
    friend inline IntegerBigMod operator*(IntegerBigMod lhs, const IntegerBigMod& rhs) { lhs *= rhs; return lhs; }
};// __attribute__ ((__aligned__(ALIGNMENT)));

template <int64_t q>
class IntegerBigPow2Mod : public IntegerBigMod<q>
{
private:
    static const int64_t B = (q == Q) ? SPLIT_Q : SPLIT_Qp;
    static const int64_t split_shift = 28;
    static const int64_t split_mask = (1 << split_shift) - 1;
    static const int64_t shift = 56;
    static const int64_t mask = (1LL << shift) - 1;
    using IntegerBigMod<q>::value;
    using IntegerBigMod<q>::IntegerBigMod;
public:
    using IntegerMod<q>::operator*;
    /**
     * @brief Default constructor
     */
    IntegerBigPow2Mod() {  }
    /**
     * @brief Copy constructor for a different modulus
     * @param src the element to reinterpret and duplicate
     * @tparam q2 the modulus of the src element
     */
    template<int64_t q2>
    IntegerBigPow2Mod(const IntegerBigPow2Mod<q2> &src)
    {
        value = src.get_value() & mask;
    }
    /**
     * @brief Constructor from signed integer type
     * @param src the element to interpret and store
     */
    IntegerBigPow2Mod(const int64_t &src)
    {
        value = src & mask;
        value -= (value / (q/2+1)) * q;
    }
    /**
     * @brief Accessor to the split threshold
     * @return the threshold
     */
    inline static int64_t get_split_threshold()
    {
        return B;
    }
    /**
     * @brief Split the element into two halves, one with the most significant part and the other with the least significant part
     * @param msb the most significant part
     * @param lsb the least signigicant part
     * @return the threshold value between the part
     */
    inline int64_t split(IntegerBigPow2Mod &msb, IntegerBigPow2Mod &lsb) const
    {
        msb = value / B;
        lsb = value % B;
        //msb = value >> split_shift;
        //lsb = value & split_mask;
        //std::cout << value << " = " << msb << " * " << B << " + " << lsb << std::endl;
        return B;
    }
    /**
     * @brief Recompose the element from two halves
     * @param msb the most significant part
     * @param lsb the least signigicant part
     * @return the threshold value between the part
     */
    inline void unsplit(IntegerBigPow2Mod &msb, IntegerBigPow2Mod &lsb)
    {
        value = (msb.get_value() << split_shift) + lsb.get_value();
        //value = (msb.get_value() * B) + lsb.get_value();
        //std::cout << value << " = " << msb << " * " << B << " + " << lsb << std::endl;
    }
    IntegerBigPow2Mod& operator*=(const IntegerBigPow2Mod& rhs)
    {
        const int64_t lhs_m = value >> split_shift;
        const int64_t lhs_l = value & split_mask;
        const int64_t rhs_m = rhs.value >> split_shift;
        const int64_t rhs_l = rhs.value & split_mask;

        value = ((lhs_l * rhs_m + lhs_m * rhs_l) & split_mask) << split_shift;
        value += lhs_l * rhs_l;

        value -= (value / (q/2+1)) * q;

        return *this;
    }
    friend inline IntegerBigPow2Mod operator-(IntegerBigPow2Mod lhs) { return -lhs.value; }
    friend inline IntegerBigPow2Mod operator*(IntegerBigPow2Mod lhs, const IntegerBigPow2Mod& rhs) { lhs *= rhs; return lhs; }
    inline IntegerBigPow2Mod& operator+=(const IntegerBigPow2Mod& rhs)
    {
        value += rhs.value;
        value -= (value / (q/2+1)) * q;
        return *this;
    }
    inline IntegerBigPow2Mod& operator-=(const IntegerBigPow2Mod& rhs)
    {
        *this += -rhs;
        return *this;
    }
};// __attribute__ ((__aligned__(ALIGNMENT)));

template class IntegerMod<T>;
template class IntegerMod<P1*P2>;
template class IntegerBigPow2Mod<Q>;
template class IntegerBigMod<Qcrt>;
template class IntegerPrimeMod<P1>;
template class IntegerPrimeMod<P2>;

typedef IntegerMod<T>         Zt;
typedef IntegerMod<P1*P2>     Zp12;
typedef IntegerBigPow2Mod<Q>  Zq;
typedef IntegerBigPow2Mod<Qp> Zqp;
typedef IntegerBigMod<Qcrt>   Zqcrt;
typedef IntegerPrimeMod<P1>   Zp1;
typedef IntegerPrimeMod<P2>   Zp2;

#endif // HE8_INTEGER_MOD_H
