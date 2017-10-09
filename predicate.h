#ifndef HE8_PREDICATE_H
#define HE8_PREDICATE_H

/**
  * @file predicate.h
  * @date Spring 2017
  * @author Guillaume Bonnoron
  * @brief Declaration of the class Fun
  *
  * The Fun object stores everything needed to compute the function. The most important aspect is its operator ()
  */

#include <iostream>
#include "param.h"
#include "integer_mod.h"

/**
 * @class Fun
 * @brief The container with its call operator
 */
class Fun
{
public:
    /**
     * @brief Allows to evaluate the function represented by the object
     * @param i the point where to evaluate the function
     * @return F(i)
     */
    inline Zt operator()(const Zp12 i) const
    {
        return eval((int64_t)std::floor((double)i * (double)T / (double)(P1*P2) + 0.5));
    }
    inline Zt eval(const Zt i) const
    {
#if INPUT_BIT > 10
        /** 63 bits threshold **/
        int64_t input = (int64_t)i;
        int64_t result = (input < 0) ? 1 : 0;
        //std::cout << i << " > " << input << " => " << __builtin_popcount(input) << " ==> " << result << std::endl;
#else
        /** 6 bits parity **/
        int64_t input = (int64_t)i;
        input = input & 0x3F;
        int64_t result = __builtin_popcount(input) % 2;
#endif
        Zt res(result);
        return res;
    }
};

#endif  // HE8_PREDICATE_H

