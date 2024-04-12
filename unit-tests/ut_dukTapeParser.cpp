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

BOOST_AUTO_TEST_SUITE_END()
