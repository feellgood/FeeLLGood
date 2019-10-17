/*
 * Implementation of MagnetizationParser.
 */

#include "mag_parser.h"
#include <exprtk.hpp>
#include <iostream>

/**
 * \brief Internal implementation of MagnetizationParser.
 *
 * This implementation is based on the
 * [exprtk](http://www.partow.net/programming/exprtk/index.html)
 * expression parser. That parser is implemented as a single huge,
 * template-rich header file, that takes very long to compile. For this
 * reason, we hide this implementation from the MagnetizationParser
 * interface using the pimpl (Pointer to IMPLementation) idiom.
 */
struct MagnetizationParser::Impl {
    Impl() {
        s_table.add_variable("x", x);
        s_table.add_variable("y", y);
        s_table.add_variable("z", z);
        s_table.add_constants();
        Mx.register_symbol_table(s_table);
        My.register_symbol_table(s_table);
        Mz.register_symbol_table(s_table);
    }

    /**
     * Internal implementation of
     * MagnetizationParser::set_expressions().
     */
    void set_expressions(
            const std::string &str_Mx,
            const std::string &str_My,
            const std::string &str_Mz) {
        exprtk::parser<double> parser;
        parser.compile(str_Mx, Mx);
        parser.compile(str_My, My);
        parser.compile(str_Mz, Mz);
    }

    /**
     * Internal implementation of
     * MagnetizationParser::get_magnetization().
     */
    Pt::pt3D get_magnetization(const Pt::pt3D &p) {
        x = p.x();
        y = p.y();
        z = p.z();
        Pt::pt3D mag = Pt::pt3D(Mx.value(), My.value(), Mz.value());
        mag.normalize();
        return mag;
    }

    double x; /**< working variable x for the exprtk parser */
    double y; /**< working variable y for the exprtk parser */
    double z; /**< working variable z for the exprtk parser */
    exprtk::symbol_table<double> s_table; /**< symbol table for the exprtk parser  */
    exprtk::expression<double> Mx; /**< Mx expression */
    exprtk::expression<double> My; /**< My expression */
    exprtk::expression<double> Mz; /**< Mz expression */
};

/*
 * MagnetizationParser is just a proxy for it's internal implementation.
 */

MagnetizationParser::MagnetizationParser() : pimpl(new Impl) {}

MagnetizationParser::~MagnetizationParser() = default;

void MagnetizationParser::set_expressions(
            const std::string &Mx,
            const std::string &My,
            const std::string &Mz)
{
    pimpl->set_expressions(Mx, My, Mz);
}

Pt::pt3D MagnetizationParser::get_magnetization(const Pt::pt3D &p)
{
    return pimpl->get_magnetization(p);
}
