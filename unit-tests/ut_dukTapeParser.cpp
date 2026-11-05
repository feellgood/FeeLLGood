#define BOOST_TEST_MODULE dukTapeParserTest

#include <boost/test/unit_test.hpp>

#include <cstdlib>
#include <iostream>
#include <string>

#include "expression_parser.h"
#include "ut_config.h"  // for tolerance UT_TOL macro

BOOST_AUTO_TEST_SUITE(ut_dukTapeParser)

/*-----------------------------------------------------*/
/* test ExpressionParser with function "function(t) {return 0.1*t*t;}" */
/*-----------------------------------------------------*/

BOOST_AUTO_TEST_CASE(dukTape_ScalarParser, *boost::unit_test::tolerance(UT_TOL))
    {
    std::string str_expr("function(t) {return 0.1*t*t;}");
    ExpressionParser func_parser;
    func_parser.set_function(str_expr);
    double result = func_parser.get_scalar(1.0);
    BOOST_TEST( result == 0.1);
    result = func_parser.get_scalar(2.0);
    BOOST_TEST( result == 0.4);
    }

/*-----------------------------------------------------*/
/* test ExpressionParser with function "function(x,y,z) {return [(x*y*z),(x+y),(y-z)];}" */
/*-----------------------------------------------------*/

BOOST_AUTO_TEST_CASE(dukTape_VectorParser, *boost::unit_test::tolerance(UT_TOL))
    {
    std::string str_expr("function(x,y,z) {return [(x*y*z),(x+y),(y-z)];}");
    ExpressionParser func_parser;
    func_parser.set_function(str_expr);
    Eigen::Vector3d val(1.0,2.0,3.0);
    Eigen::Vector3d result = func_parser.get_vector(val);
    BOOST_TEST( result[0] == 6.0);
    BOOST_TEST( result[1] == 3.0);
    BOOST_TEST( result[2] == -1.0);
    val = Eigen::Vector3d(0.0,2.0,3.0);
    result = func_parser.get_vector( val );
    BOOST_TEST( result[0] == 0.0);
    BOOST_TEST( result[1] == 2.0);
    BOOST_TEST( result[2] == -1.0);
    }

/*-----------------------------------------------------*/
/* test ExpressionParser calling a global function jn_multiplier(t), defined before and */
/* after the function that calls it (as Bext calling the multiplier of the current densities) */
/*-----------------------------------------------------*/

BOOST_AUTO_TEST_CASE(dukTape_GlobalFunction, *boost::unit_test::tolerance(UT_TOL))
    {
    std::string str_jn("function(t) {return t < 1 ? t : 1;}");
    ExpressionParser time_parser;
    time_parser.define_global_function("jn_multiplier", str_jn);
    time_parser.set_function("function(t) {return 2*jn_multiplier(t);}");
    BOOST_TEST( time_parser.get_scalar(0.25) == 0.5);
    BOOST_TEST( time_parser.get_scalar(3.0) == 2.0);

    ExpressionParser field_parser;
    field_parser.set_expressions("t", "jn_multiplier(t)", "0", "-jn_multiplier(t)");
    field_parser.define_global_function("jn_multiplier", str_jn);
    Eigen::Vector3d result = field_parser.get_vector(0.5);
    BOOST_TEST( result[0] == 0.5);
    BOOST_TEST( result[1] == 0.0);
    BOOST_TEST( result[2] == -0.5);
    }

/*-----------------------------------------------------*/
/* test ExpressionParser::is_set(), true once a function is given, after evaluations too */
/*-----------------------------------------------------*/

BOOST_AUTO_TEST_CASE(dukTape_IsSet, *boost::unit_test::tolerance(UT_TOL))
    {
    ExpressionParser parser;
    BOOST_TEST( !parser.is_set() );
    parser.define_global_function("jn_multiplier", "function(t) {return 2*t;}");
    BOOST_TEST( !parser.is_set() );
    parser.set_function("function(t) {return jn_multiplier(t) + 1;}");
    BOOST_TEST( parser.is_set() );
    BOOST_TEST( parser.get_scalar(1.5) == 4.0 );
    BOOST_TEST( parser.is_set() );
    }

BOOST_AUTO_TEST_SUITE_END()
