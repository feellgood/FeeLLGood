#define BOOST_TEST_MODULE algebraTest

#include <boost/test/unit_test.hpp>

#include <iostream>

#include "ut_config.h"  // for tolerance UT_TOL macro
#include "../algebra/algebra.h"

using namespace algebra;

BOOST_AUTO_TEST_SUITE(ut_algebra)

/**
test on scaled, dot, norm, sub functions
*/
BOOST_AUTO_TEST_CASE(test_scaled, *boost::unit_test::tolerance(UT_TOL))
    {
    std::vector<double> x {1,4,-2,sqrt(42)};
    double alpha(2.0);
    std::vector<double> y(4);
    y.assign(x.begin(),x.end()); // y <- x
    scaled(alpha,y); // y *= alpha
    
    double result = dot(y,y); // result = scalar product y with himself
    std::cout << "dot(y,y)= " << result << std::endl;
    BOOST_CHECK( result == 252.0 );
    
    scaled(1.0/alpha,y); // y *= 1.0/alpha
    
    sub(x,y); // y -= x 
    BOOST_CHECK( norm(y) == 0.0 );
    }

/** test on p_direct */
BOOST_AUTO_TEST_CASE(test_p_direct, *boost::unit_test::tolerance(UT_TOL))
    {
    std::vector<double> x {1,4,-2,sqrt(42)};
    std::vector<double> y {1,-3,2,0.5};
    std::vector<double> z(4);
    p_direct(x,y,z);

    BOOST_CHECK( z[0] == 1.0 );
    BOOST_CHECK( z[1] == -12.0 );
    BOOST_CHECK( z[2] == -4.0 );
    BOOST_CHECK( z[3] == sqrt(42)/2.0 );
    }

/** test on add */
BOOST_AUTO_TEST_CASE(test_add, *boost::unit_test::tolerance(UT_TOL))
    {
    std::vector<double> x {1,4,-2,sqrt(42)};
    std::vector<double> y {1,-3,2,0.5};
    add(x,y); // y += x

    BOOST_CHECK( y[0] == 2.0 );
    BOOST_CHECK( y[1] == 1.0 );
    BOOST_CHECK( y[2] == 0.0 );
    BOOST_CHECK( y[3] == (0.5 + sqrt(42)) );
    }

/** test on sub */
BOOST_AUTO_TEST_CASE(test_sub, *boost::unit_test::tolerance(UT_TOL))
    {
    std::vector<double> x {1,4,-2,sqrt(42)};
    std::vector<double> y {1,-3,2,0.5};
    sub(x,y); // y -= x

    BOOST_CHECK( y[0] == 0.0 );
    BOOST_CHECK( y[1] == -7.0 );
    BOOST_CHECK( y[2] == 4.0 );
    BOOST_CHECK( y[3] == (0.5 - sqrt(42)) );
    }

BOOST_AUTO_TEST_SUITE_END()

