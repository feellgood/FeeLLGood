#define BOOST_TEST_MODULE algTest

#include <boost/test/unit_test.hpp>

#include <iostream>

#include "ut_config.h"  // for tolerance UT_TOL macro
#include "../alg/alg.h"

BOOST_AUTO_TEST_SUITE(ut_alg)

/**
test on scaled, dot, norm, sub functions
*/
BOOST_AUTO_TEST_CASE(test_scaled, *boost::unit_test::tolerance(UT_TOL))
    {
    std::vector<double> x {1,4,-2,sqrt(42)};
    double alpha(2.0);
    std::vector<double> y(4);
    alg::scaled(x,alpha,y); // y = alpha*x
    
    double result = alg::dot(y,y); // result = scalar product y with himself
    std::cout << "dot(y,y)= " << result << std::endl;
    BOOST_CHECK( result == 252.0 );
    
    alg::scaled(1.0/alpha,y); // y *= 1.0/alpha
    
    alg::sub(x,y); // y -= x 
    BOOST_CHECK( alg::norm(y) == 0.0 );
    }

BOOST_AUTO_TEST_SUITE_END()

