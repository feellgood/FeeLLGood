/*
 * Implementation of MagnetizationParser and TimeDepFieldParser.
 */

#include "mag_parser.h"
#include <duktape.h>
#include <iostream>

/***********************************************************************
 * Code shared between both parsers.
 */

/* Copy all the contents of the Math object (sqrt, sin, log, PI...) to the global object. */
static const char js_library[] = "Object.getOwnPropertyNames(Math).forEach("
                                 "    function(name) { globalThis[name] = Math[name]; }"
                                 ");";

/**
 * Generic parser and evaluator, parametrized by the type of the parameter:
 *  - T = double for the field, which is a function of time
 *  - T = Pt::pt3D for the initial magnetization, which is a function of position.
 */
template<typename T>
struct VectorParser
    {
    VectorParser()
        {
        ctx = duk_create_heap_default();
        duk_int_t err = duk_peval_string(ctx, js_library);
        die_if_error(err);
        duk_pop(ctx);  // drop the result of the evaluation
        }

    /**
     * Abort with an suitable error message if `err` is an actual error,
     * in which case the top of the stack is assumed to hold the corresponding Error object.
     */
    void die_if_error(duk_int_t err)
        {
        if (err == DUK_ERR_NONE) return;
        std::cerr << "Script " << duk_safe_to_string(ctx, -1) << '\n';
        exit(EXIT_FAILURE);
        }

    /**
     * Build a function that takes the given parameters and evaluates the given expression.
     * Push it to the stack.
     */
    void push_function(std::string expression)
        {
        std::string js_function = "function(" + js_params + ") { return (" + expression + "); }";
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

    /**
     * Compile the expressions into JavaScript functions.
     * Leave them on the stack at positions (0, 1, 2).
     */
    void set_expressions(const std::string &expr_x, const std::string &expr_y,
                         const std::string &expr_z)
        {
        duk_set_top(ctx, 0);  // clear the stack
        push_function(expr_x);
        push_function(expr_y);
        push_function(expr_z);
        }

    /**
     * Dummy method to push the parameter to the stack, will be specialized later.
     */
    void push(const T &value) {}

    /**
     * Get a single vector component. `component_index` is either 0, 1, or 2.
     */
    double get_vector_component(const T &param, int component_index)
        {
        duk_dup(ctx, component_index);  // copy the function to the top of the stack
        push(param);
        duk_int_t err = duk_pcall(ctx, dim);
        die_if_error(err);
        double val = duk_get_number(ctx, -1);
        duk_pop(ctx);
        return val;
        }

    /**
     * Get the vector corresponding to the provided parameter.
     */
    Pt::pt3D get_vector(const T &param)
        {
        double x = get_vector_component(param, 0);
        double y = get_vector_component(param, 1);
        double z = get_vector_component(param, 2);
        return Pt::pt3D(x, y, z);
        }

    duk_context *ctx;             /**< Ducktape context holding the functions */
    static const int dim;         /**< dimension of the parameter (1 or 3) */
    static std::string js_params; /**< JavaScript parameter list: either "t" or "x,y,z" */
    };

/*
 * VectorParser template specializations for `double` and `Pt::pt3D`.
 * Hide them from Doxygen, as it does not seem to understand them:
 * it says "warning: no matching class member found for [...]"
 */
//! @cond Doxygen_Suppress

template<>
const int VectorParser<double>::dim = 1;

template<>
std::string VectorParser<double>::js_params = "t";

template<>
void VectorParser<double>::push(const double &value)
    {
    duk_push_number(ctx, value);
    }

template<>
const int VectorParser<Pt::pt3D>::dim = 3;

template<>
std::string VectorParser<Pt::pt3D>::js_params = "x,y,z";

template<>
void VectorParser<Pt::pt3D>::push(const Pt::pt3D &value)
    {
    duk_push_number(ctx, value.x());
    duk_push_number(ctx, value.y());
    duk_push_number(ctx, value.z());
    }

//! @endcond

/** ********************************************************************
 * \brief Internal implementation of MagnetizationParser.
 */
struct MagnetizationParser::Impl3Dprm : VectorParser<Pt::pt3D>
    {
    };

/*
 * MagnetizationParser is just a proxy for it's internal implementation.
 */

MagnetizationParser::MagnetizationParser() : pimpl3Dprm(new Impl3Dprm) {}

MagnetizationParser::~MagnetizationParser() = default;

void MagnetizationParser::set_expressions(const std::string &Mx, const std::string &My,
                                          const std::string &Mz)
    {
    pimpl3Dprm->set_expressions(Mx, My, Mz);
    }

Pt::pt3D MagnetizationParser::get_magnetization(const Pt::pt3D &p) const
    {
    Pt::pt3D M = pimpl3Dprm->get_vector(p);
    M.normalize();
    return M;
    }

/** ********************************************************************
 * \brief Internal implementation of TimeDepFieldParser.
 */
struct TimeDepFieldParser::Impl1Dprm : VectorParser<double>
    {
    };

/*
 * TimeDepFieldParser is a proxy for it's internal implementation.
 */

TimeDepFieldParser::TimeDepFieldParser() : pimpl1Dprm(new Impl1Dprm) {}

TimeDepFieldParser::~TimeDepFieldParser() = default;

void TimeDepFieldParser::set_expressions(const std::string &Bx, const std::string &By,
                                         const std::string &Bz)
    {
    pimpl1Dprm->set_expressions(Bx, By, Bz);
    }

Pt::pt3D TimeDepFieldParser::get_timeDepField(const double t_val) const
    {
    return pimpl1Dprm->get_vector(t_val);
    }
