#include "pyRepr.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdlib>

namespace fg
    {
std::string pyRepr(double x)
    {
    if (std::isnan(x))
        {
        return "nan";
        }
    if (std::isinf(x))
        {
        return x > 0 ? "inf" : "-inf";
        }

    // shortest round trip digits, in the form [-]d[.ddd]e(+|-)XX
    char buf[64];
    const auto res = std::to_chars(buf, buf + sizeof(buf), x, std::chars_format::scientific);
    const std::string sci(buf, res.ptr);
    const std::size_t ePos = sci.find('e');
    std::string mantissa = sci.substr(0, ePos);
    const int exponent = std::atoi(sci.c_str() + ePos + 1);

    std::string sign;
    if (mantissa[0] == '-')
        {
        sign = "-";
        mantissa.erase(0, 1);
        }
    std::string digits;
    for (char c : mantissa)
        {
        if (c != '.')
            {
            digits += c;
            }
        }
    const int n = static_cast<int>(digits.size());

    if (exponent >= -4 && exponent < 16)
        {
        std::string s;
        if (exponent >= 0)
            {
            std::string intPart = digits.substr(0, std::min(n, exponent + 1));
            intPart.append(std::max(0, exponent + 1 - n), '0');
            std::string fracPart = (n > exponent + 1) ? digits.substr(exponent + 1) : "0";
            s = intPart + "." + fracPart;
            }
        else
            {
            s = "0." + std::string(-exponent - 1, '0') + digits;
            }
        return sign + s;
        }

    std::string s = digits.substr(0, 1);
    if (n > 1)
        {
        s += "." + digits.substr(1);
        }
    const int absExp = std::abs(exponent);
    s += (exponent < 0 ? "e-" : "e+");
    if (absExp < 10)
        {
        s += "0";
        }
    s += std::to_string(absExp);
    return sign + s;
    }
    }  // namespace fg
