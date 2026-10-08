// See LICENCE file at project root
#ifndef FMATH_HPP
#define FMATH_HPP

#include <cmath>

/**
 * @author Berenger Bramas (berenger.bramas@inria.fr)
 * Please read the license
 *
 * indirections to std math.
 * some specialized templates for FReal == float | double
 */

/** namespace to grab altogether indirections to STL math */
namespace FMath
{
    /** To get absolute value */
    template <class NumType>
    static NumType Abs(const NumType inV){ return (inV < 0 ? -inV : inV); }

    /** To get max between 2 values */
    template <class NumType>
    static NumType Max(const NumType inV1, const NumType inV2) { return (inV1 > inV2 ? inV1 : inV2); }

    /** To get min between 2 values */
    template <class NumType>
    static NumType Min(const NumType inV1, const NumType inV2) { return (inV1 < inV2 ? inV1 : inV2); }

    /** To get pow of 2 */
    static int pow2(const int power) { return (1 << power); }

    /** To know if a value is between two others */
    template <class NumType>
    static bool Between(const NumType inValue, const NumType inMin, const NumType inMax) 
        { return ( inMin <= inValue && inValue < inMax ); }

    /** sqrt (double version) */
    static double Sqrt(const double inValue) { return sqrt(inValue); }

    /** atan2 (double version), return value is given in radians and is in [-pi,pi], inclusive. */
    static double Atan2(const double inValue1,const double inValue2) { return atan2(inValue1,inValue2); }

} // end namespace

#endif //FMATH_HPP

