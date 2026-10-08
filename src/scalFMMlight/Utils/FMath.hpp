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
    /** To get pow of 2 */
    static int pow2(const int power) { return (1 << power); }

} // end namespace

#endif //FMATH_HPP

