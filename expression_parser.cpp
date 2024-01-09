/*
 * Implementation of VectorParser.
 */

#include "expression_parser.h"
#include <iostream>

/* Seed the JavaScript interpreter with some useful math functions. */
static const char js_library[] = R"--(
/* Copy all the contents of the Math object (sqrt, sin, log, PI...) to the global object. */
Object.getOwnPropertyNames(Math).forEach(
    function(name) { globalThis[name] = Math[name]; }
);

/* Shifted exp() and log(). */
function expm1(x) { return abs(x)<2e-8 ? x : exp(x) - 1; }
function log1p(x) { return abs(x)<2e-8 ? x : log(1 + x); }

/* Hyperbolic functions. */
function cosh(x) { return exp(x - LN2) + exp(-x - LN2); }
function sinh(x) { return abs(x)<1e-5 ? x : exp(x - LN2) - exp(-x - LN2); }
function tanh(x) { return abs(x)>20 ? sign(x) : expm1(2*x)/(exp(2*x) + 1); }

/* Inverse hyperbolic functions. */
function acosh(x) { return x>1e8 ? log(x) + LN2 : log(x + sqrt(x*x - 1)); }
function asinh(x) { return abs(x)<1e-5 ? x : sign(x) * (log(abs(x)/2 + hypot(1, x)/2) + LN2); }
function atanh(x) { return abs(x)<1e-5 ? x : log((1 + x)/(1 - x)) / 2; }
)--";

VectorParser::VectorParser()
    {
    ctx = duk_create_heap_default();
    duk_int_t err = duk_peval_string(ctx, js_library);
    die_if_error(err);
    duk_pop(ctx);  // drop the result of the evaluation
    }

VectorParser::~VectorParser() { duk_destroy_heap(ctx); }

void VectorParser::die_if_error(duk_int_t err) const
    {
    if (err == DUK_ERR_NONE) return;
    std::cerr << "Script " << duk_safe_to_string(ctx, -1) << '\n';
    exit(EXIT_FAILURE);
    }

// Compile a function expression and leave it as the sole value on the Duktape stack.
void VectorParser::set_function(const std::string &js_function) const
    {
    duk_set_top(ctx, 0);  // clear the stack
    duk_int_t err = duk_pcompile_string(ctx, DUK_COMPILE_FUNCTION, js_function.c_str());
    die_if_error(err);
    if (duk_get_top(ctx) != 1)
        {
        std::cerr << "Compilation of JavaScript function damaged the stack.\n";
        exit(EXIT_FAILURE);
        }
    }

void VectorParser::set_expressions(const std::string &parameters, const std::string &expr_x,
                                   const std::string &expr_y, const std::string &expr_z)
    {
    set_function("function(" + parameters + ") { return [(" + expr_x + "), (" + expr_y + "), ("
                 + expr_z + ")]; }");
    }

double VectorParser::get_vector_component(int idx) const
    {
    duk_push_int(ctx, idx);                      // [ ... v ] -> [ ... v idx ]
    duk_bool_t success = duk_get_prop(ctx, -2);  //           -> [ ... v v[idx] ]
    if (!success)
        {
        std::cerr << "JavaScript function did not return a 3-vector.\n";
        exit(EXIT_FAILURE);
        }
    double val = duk_get_number(ctx, -1);
    duk_pop(ctx);  // -> [ ... v ]
    return val;
    }

Eigen::Vector3d VectorParser::compute_vector(int argument_count) const
    {
    duk_int_t err = duk_pcall(ctx, argument_count);  // [ f f args... ] -> [ f v ]
    die_if_error(err);
    double x = get_vector_component(0);
    double y = get_vector_component(1);
    double z = get_vector_component(2);
    duk_pop(ctx);  // -> [ f ]
    return Eigen::Vector3d(x, y, z);
    }

Eigen::Vector3d VectorParser::get_vector(double arg) const
    {
    duk_dup(ctx, -1);           // -> [ f f ]
    duk_push_number(ctx, arg);  // -> [ f f arg ]
    return compute_vector(1);
    }

Eigen::Vector3d VectorParser::get_vector(const Eigen::Ref<Eigen::Vector3d> arg) const
    {
    duk_dup(ctx, -1);               // -> [ f f ]
    duk_push_number(ctx, arg.x());  // -> [ f f arg.x ]
    duk_push_number(ctx, arg.y());  // -> [ f f arg.x arg.y ]
    duk_push_number(ctx, arg.z());  // -> [ f f arg.x arg.y arg.z ]
    return compute_vector(3);
    }
