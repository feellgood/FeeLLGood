/*
 * MagnetizationParser: parser for the magnetization expressions.
 */

#ifndef mag_parser_h
#define mag_parser_h

#include "pt3D.h"
#include <memory>
#include <string>

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
class MagnetizationParser
    {
public:
    MagnetizationParser();
    ~MagnetizationParser();

    /**
     * Compile the expressions for the magnetization components.
     */
    void set_expressions(const std::string &Mx, const std::string &My, const std::string &Mz);

    /**
     * Evaluate the magnetization at point `p`. Returns the _normalized_
     * magnetization. This should only called _after_ defining the
     * expressions with set_expressions().
     */
    Pt::pt3D get_magnetization(const Pt::pt3D &p) const;

private:
    class Impl3Dprm;
    std::unique_ptr<Impl3Dprm> pimpl3Dprm; /**< Pointer to the internal implementation. */
    };

/**
 * \brief Parse and evaluate the expressions that give the components of
 * the time dependant applied field as a function of time t.
 *
 * This class **cannot be copied**, as it holds a `unique_ptr`.
 */

class TimeDepFieldParser
    {
public:
    TimeDepFieldParser();
    ~TimeDepFieldParser();

    /** setter for the expressions of applied field B */
    void set_expressions(const std::string &Bx, const std::string &By, const std::string &Bz);

    /** getter for time dependant field value */
    Pt::pt3D get_timeDepField(const double t_val) const;

private:
    class Impl1Dprm;
    std::unique_ptr<Impl1Dprm> pimpl1Dprm; /**< internal implementation pointer */
    };

#endif /* mag_parser_h */
