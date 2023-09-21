/*
 * Implementation of VectorParser.
 */

#include "expression_parser.h"
#include <iostream>

/* Copy all the contents of the Math object (sqrt, sin, log, PI...) to the global object. */
static const char js_library[] = "Object.getOwnPropertyNames(Math).forEach("
                                 "    function(name) { globalThis[name] = Math[name]; }"
                                 ");";

VectorParser::VectorParser()
    {
    ctx = duk_create_heap_default();
    duk_int_t err = duk_peval_string(ctx, js_library);
    die_if_error(err);
    duk_pop(ctx);  // drop the result of the evaluation
    }

void VectorParser::die_if_error(duk_int_t err) const
    {
    if (err == DUK_ERR_NONE) return;
    std::cerr << "Script " << duk_safe_to_string(ctx, -1) << '\n';
    exit(EXIT_FAILURE);
    }

void VectorParser::push_function(const std::string &js_function) const
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
    push_function("function(" + parameters + ") { return [(" + expr_x + "), (" + expr_y + "), ("
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

Pt::pt3D VectorParser::compute_vector(int argument_count) const
    {
    duk_int_t err = duk_pcall(ctx, argument_count);  // [ f f args... ] -> [ f v ]
    die_if_error(err);
    double x = get_vector_component(0);
    double y = get_vector_component(1);
    double z = get_vector_component(2);
    duk_pop(ctx);  // -> [ f ]
    return Pt::pt3D(x, y, z);
    }

Pt::pt3D VectorParser::get_vector(double arg) const
    {
    duk_dup(ctx, -1);           // -> [ f f ]
    duk_push_number(ctx, arg);  // -> [ f f arg ]
    return compute_vector(1);
    }

Pt::pt3D VectorParser::get_vector(const Pt::pt3D &arg) const
    {
    duk_dup(ctx, -1);               // -> [ f f ]
    duk_push_number(ctx, arg.x());  // -> [ f f arg.x ]
    duk_push_number(ctx, arg.y());  // -> [ f f arg.x arg.y ]
    duk_push_number(ctx, arg.z());  // -> [ f f arg.x arg.y arg.z ]
    return compute_vector(3);
    }
