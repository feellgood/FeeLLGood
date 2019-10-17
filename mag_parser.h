/*
 * MagnetizationParser: parser for the magnetization expressions.
 */

#ifndef mag_parser_h
#define mag_parser_h

#include "pt3D.h"
#include <string>
#include <memory>

/**
 * \brief Parse and evaluate the expressions that give the components of
 * the initial magnetization as a function of the position coordinates
 * (x, y, z).
 *
 * A single parser object is used for the three components of the
 * magnetization because the independent variables (node coordinates)
 * are shared by the three expressions.
 *
 * This class **cannot be copied**, as it holds a `unique_ptr`.
 */
class MagnetizationParser {
public:
    MagnetizationParser();
    ~MagnetizationParser();

    /**
     * Compile the expressions for the magnetization components.
     */
    void set_expressions(
            const std::string &Mx,
            const std::string &My,
            const std::string &Mz);

    /**
     * Evaluate the magnetization at point `p`. Returns the _normalized_
     * magnetization. This should only called _after_ defining the
     * expressions with set_expressions().
     */
    Pt::pt3D get_magnetization(const Pt::pt3D &p);

private:
    class Impl;
    std::unique_ptr<Impl> pimpl; /**< Pointer to the internal implementation. */
};

#endif /* mag_parser_h */
