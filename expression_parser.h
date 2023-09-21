/*
 * Parser/evaluator for the analytic expressions given by the user:
 * the initial magnetization M(x, y, z) and the applied field B(t).
 */

#ifndef expression_parser_h
#define expression_parser_h

#include "pt3D.h"
#include <duktape.h>
#include <string>

/**
 * This class handles a JavaScript function that computes the three components of a vector, which
 * depends on either a single parameter t or on (x, y, z).
 */
class VectorParser
    {
public:
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

#endif /* expression_parser_h */
