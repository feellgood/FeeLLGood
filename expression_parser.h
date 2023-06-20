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
 * code shared between them. It handles three expressions that compute the three components of a
 * vector, which depend on unspecified parameters.
 */
class VectorParser
    {
protected:
    VectorParser();

    /**
     * Compile the expressions of the x, y and z components, which depend on the specified
     * parameters, into JavaScript functions. Leave the compiled functions on the Duktape stack at
     * positions (0, 1, 2).
     */
    void set_expressions(const std::string &parameters, const std::string &expr_x,
                         const std::string &expr_y, const std::string &expr_z);

    /**
     * Prepare the Duktape interpreter to call the function at the given stack position.
     * Positions (0, 1, 2) correspond to components (x, y, z) of the computed vector.
     */
    void prepare_call(int position) const;

    /**
     * Push a function argument to the Duktape stack.
     */
    void push_argument(double value) const;

    /**
     * Compute a component of the vector. This must be called after prepare_call() and one or more
     * calls to push_argument(). `argument_count` should match the number of arguments pushed.
     */
    double get_component(int argument_count) const;

private:
    /**
     * Abort with an suitable error message if `err` is an actual error, in which case the top of
     * the stack is assumed to hold the corresponding Error object.
     */
    void die_if_error(duk_int_t err) const;

    /**
     * Build a function that takes the given parameters and evaluates the given expression.
     * Push it to the stack.
     */
    void push_function(const std::string &parameters, const std::string &expression);

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
    Pt::pt3D get_magnetization(const Pt::pt3D &p) const;

private:
    /**
     * Get a single vector component. `component_index` should be either 0, 1, or 2.
     */
    double get_vector_component(const Pt::pt3D &p, int component_index) const;
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
    Pt::pt3D get_timeDepField(const double t) const;

private:
    /**
     * Get a single vector component. `component_index` should be either 0, 1, or 2.
     */
    double get_vector_component(double t, int component_index) const;
    };

#endif /* expression_parser_h */
