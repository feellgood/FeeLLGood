/*
 * Parsers/evaluators for the analytic expressions given by the user:
 *  - MagnetizationParser: parser for the initial magnetization M(x, y, z)
 *  - TimeDepFieldParser: parser for the applied field B(t)
 */

#ifndef expression_parser_h
#define expression_parser_h

#include "pt3D.h"
#include <duktape.h>
#include <string>

/**
 * Generic parser and evaluator. This is the parent class of the other two parsers, and holds the
 * code shared between them. It handles a JavaScript function that computes the three components of
 * a vector, which depends on either a single parameter t or on (x, y, z).
 */
class VectorParser
    {
protected:
    VectorParser();

    /**
     * Compile the expressions of the x, y and z components, which depend on the specified
     * parameters, into a JavaScript function returning a 3-element array. Leave the compiled
     * function as the sole value on the Duktape stack.
     */
    void set_expressions(const std::string &parameters, const std::string &expr_x,
                         const std::string &expr_y, const std::string &expr_z);

    /**
     * Compute a vector from the given scalar argument.
     */
    Pt::pt3D get_vector(double arg) const;

    /**
     * Compute a vector from the given vector argument.
     */
    Pt::pt3D get_vector(const Pt::pt3D &arg) const;

private:
    /**
     * Abort with an suitable error message if `err` is an actual error, in which case the top of
     * the stack is assumed to hold the corresponding Error object.
     */
    void die_if_error(duk_int_t err) const;

    /**
     * Push the provided function expression to the stack.
     */
    void push_function(const std::string &js_function) const;

    /**
     * Get the given component from the array at the top of the stack. Preserve the stack state.
     */
    double get_vector_component(int idx) const;

    /**
     * Compute a vector. This must be called after duk_dup(ctx, -1) and one or more calls to
     * duk_push_number(). `argument_count` should match the number of arguments pushed.
     */
    Pt::pt3D compute_vector(int argument_count) const;

    /**
     * Ducktape context holding the internal state of the interpreter.
     */
    duk_context *ctx;
    };

/**
 * Parse and evaluate the expressions that give the components of the initial magnetization as a
 * function of the position coordinates (x, y, z).
 */
class MagnetizationParser : VectorParser
    {
public:
    /**
     * Compile the expressions for the magnetization components.
     */
    void set_expressions(const std::string &Mx, const std::string &My, const std::string &Mz)
        {
        VectorParser::set_expressions("x,y,z", Mx, My, Mz);
        }

    /**
     * Evaluate the magnetization at point `p`. Returns the _normalized_ magnetization. This should
     * only be called _after_ defining the expressions with set_expressions().
     */
    Pt::pt3D get_magnetization(const Pt::pt3D &p) const { return get_vector(p).normalize(); }
    };

/**
 * Parse and evaluate the expressions that give the components of the time-dependent applied field
 * as a function of time t.
 */
class TimeDepFieldParser : public VectorParser
    {
public:
    /**
     * Compile the expressions for the field components.
     */
    void set_expressions(const std::string &Bx, const std::string &By, const std::string &Bz)
        {
        VectorParser::set_expressions("t", Bx, By, Bz);
        }

    /**
     * Evaluate the applied field at time `t`. This should only be called _after_ defining the
     * expressions with set_expressions().
     */
    Pt::pt3D get_timeDepField(const double t) const { return get_vector(t); }
    };

#endif /* expression_parser_h */
