/*
 * Implementation of MagnetizationParser and TimeDepFieldParser.
 */

#include "expression_parser.h"
#include <iostream>

/* Copy all the contents of the Math object (sqrt, sin, log, PI...) to the global object. */
static const char js_library[] = "Object.getOwnPropertyNames(Math).forEach("
                                 "    function(name) { globalThis[name] = Math[name]; }"
                                 ");";

/*******************************************************************************
 * VectorParser implementation.
 */

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

void VectorParser::push_function(const std::string &parameters, const std::string &expression)
    {
    std::string js_function = "function(" + parameters + ") { return (" + expression + "); }";
    int stack_size_before = duk_get_top(ctx);
    duk_int_t err = duk_pcompile_string(ctx, DUK_COMPILE_FUNCTION, js_function.c_str());
    die_if_error(err);
    int stack_size_after = duk_get_top(ctx);
    if (stack_size_after != stack_size_before + 1)
        {
        std::cerr << "Compilation of " << expression << " damaged the stack.\n";
        exit(EXIT_FAILURE);
        }
    }

void VectorParser::set_expressions(const std::string &parameters, const std::string &expr_x,
                                   const std::string &expr_y, const std::string &expr_z)
    {
    duk_set_top(ctx, 0);  // clear the stack
    push_function(parameters, expr_x);
    push_function(parameters, expr_y);
    push_function(parameters, expr_z);
    }

void VectorParser::prepare_call(int position) const
    {
    duk_dup(ctx, position);  // copy the function to the top of the stack
    }

void VectorParser::push_argument(double value) const { duk_push_number(ctx, value); }

double VectorParser::get_component(int argument_count) const
    {
    duk_int_t err = duk_pcall(ctx, argument_count);
    die_if_error(err);
    double val = duk_get_number(ctx, -1);
    duk_pop(ctx);
    return val;
    }

/*******************************************************************************
 * MagnetizationParser implementation.
 */

double MagnetizationParser::get_vector_component(const Pt::pt3D &p, int component_index) const
    {
    prepare_call(component_index);
    push_argument(p.x());
    push_argument(p.y());
    push_argument(p.z());
    return get_component(3);
    }

Pt::pt3D MagnetizationParser::get_magnetization(const Pt::pt3D &p) const
    {
    double Mx = get_vector_component(p, 0);
    double My = get_vector_component(p, 1);
    double Mz = get_vector_component(p, 2);
    Pt::pt3D M(Mx, My, Mz);
    M.normalize();
    return M;
    }

/*******************************************************************************
 * TimeDepFieldParser implementation.
 */

double TimeDepFieldParser::get_vector_component(double t, int component_index) const
    {
    prepare_call(component_index);
    push_argument(t);
    return get_component(1);
    }

Pt::pt3D TimeDepFieldParser::get_timeDepField(const double t) const
    {
    double Bx = get_vector_component(t, 0);
    double By = get_vector_component(t, 1);
    double Bz = get_vector_component(t, 2);
    Pt::pt3D B(Bx, By, Bz);
    return B;
    }
