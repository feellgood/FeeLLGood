#ifndef FG_PYREPR_H
#define FG_PYREPR_H

/** \file pyRepr.h
 * Formatting of floating point numbers as Python's str() and repr() do, so that the C++ tools
 * write the same text as the former Python ones.
 */

#include <string>

namespace fg
    {
/** Shortest decimal representation of x that reads back to the same double, formatted as
 * Python's repr(): fixed notation with at least one decimal ("1.0", "0.0001") when the
 * decimal exponent is in [-4, 16), scientific notation otherwise ("1e-05", "1.5e+16"), and
 * "nan", "inf", "-inf". */
std::string pyRepr(double x);
    }  // namespace fg

#endif
