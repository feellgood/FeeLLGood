#define BOOST_TEST_MODULE algebraTest

#include <boost/test/unit_test.hpp>

#include <iostream>

#include "ut_config.h"  // for tolerance UT_TOL macro
#include "../algebra/algebra.h"
#include "../algebra/cg.h"

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

/** test on scaled_add */
BOOST_AUTO_TEST_CASE(test_scaled_add, *boost::unit_test::tolerance(UT_TOL))
    {
    std::vector<double> x {1,4,-2,sqrt(42)};
    std::vector<double> y {1,-3,2,0.5};
    double alpha(2.0);
    scaled_add(x,alpha,y); // y += alpha*x

    BOOST_CHECK( y[0] == 3.0 );
    BOOST_CHECK( y[1] == 5.0 );
    BOOST_CHECK( y[2] == -2.0 );
    BOOST_CHECK( y[3] == (0.5 + 2.0*sqrt(42)) );
    }

/** test on applyMask */
BOOST_AUTO_TEST_CASE(test_applyMask)
    {
    std::vector<double> x {1,4,-2,sqrt(42)};
    std::vector<int> mask {1,3};
    applyMask(mask,x);
    
    BOOST_TEST( x[0] == 1.0 );
    BOOST_TEST( x[1] == 0.0 );
    BOOST_TEST( x[2] == -2.0 );
    BOOST_TEST( x[3] == 0.0 );
    }

/** tests on r_sparseMat built from w_sparseMat */
BOOST_AUTO_TEST_CASE(test_w_sparseMat, *boost::unit_test::tolerance(UT_TOL))
    {
    const int N=4;
    w_sparseMat m(N);
    m.insert( 1,1,3.14 );
    m.insert( 0,0,1.0 );
    m.insert( 2,2,5.0 );
    m.insert( 3,3,42.0 );
    m.insert( 1,3,-10.0 );
    m.insert( 1,3,10.0 );
    m.insert( 0,3,0.5 );
    r_sparseMat bob(m);
    std::vector<double> x {1.0,1.0,1.0,1.0};
    std::vector<double> y(N);
    m.insert( 2,2,5.0 );// m modification must not affect bob
    mult(bob,x,y); //y = m*x
    BOOST_CHECK( y[0] == 1.5 );
    BOOST_CHECK( y[1] == 3.14 );
    BOOST_CHECK( y[2] == 5.0 );
    BOOST_CHECK( y[3] == 42.0 );
    }

/** tests on r_sparseMat built from MatrixShape */
BOOST_AUTO_TEST_CASE(test_matrix_shape, *boost::unit_test::tolerance(UT_TOL))
    {
    const int N=4;
    MatrixShape shape(N);
    shape[1].insert(1);
    shape[0].insert(0);
    shape[2].insert(2);
    shape[3].insert(3);
    shape[1].insert(3);
    shape[0].insert(3);
    r_sparseMat m(shape);
    m.add(1, 1, 3.14);
    m.add(0, 0, 1);
    m.add(2, 2, 5);
    m.add(3, 3, 42);
    m.add(1, 3, -10);
    m.add(1, 3, 10);
    m.add(0, 3, 0.5);
    std::vector<double> x {1, 1, 1, 1};
    std::vector<double> y(N);
    mult(m,x,y); //y = m*x
    BOOST_CHECK( y[0] == 1.5 );
    BOOST_CHECK( y[1] == 3.14 );
    BOOST_CHECK( y[2] == 5.0 );
    BOOST_CHECK( y[3] == 42.0 );

    m.clear();
    m.add(0, 3, 1);
    m.add(1, 1, 2);
    m.add(1, 3, 3);
    m.add(2, 2, 4);
    m.add(3, 3, 6);
    mult(m, x, y);
    BOOST_CHECK( y[0] == 1 );
    BOOST_CHECK( y[1] == 5 );
    BOOST_CHECK( y[2] == 4 );
    BOOST_CHECK( y[3] == 6 );
    }

/** test on gradient conjugate algorithm */
BOOST_AUTO_TEST_CASE(test_cg, *boost::unit_test::tolerance(10.0*UT_TOL))
    {
    const int N=4;
    w_sparseMat m(N);
    m.insert( 1,1,3.14 );
    m.insert( 0,0,1.0 );
    m.insert( 2,2,5.0 );
    m.insert( 3,3,42.0 );
    m.insert( 1,3,-10.0 );
    m.insert( 1,3,10.0 );
    m.insert( 0,3,0.5 );
    r_sparseMat bob(m);
    std::vector<double> b {1.0,1.0,1.0,1.0};
    std::vector<double> x(N);
    iteration<double> algo_it("cg",1e-6,false,700);
    cg(algo_it,bob,x,b);
    double res = algo_it.get_res();
    std::cout << "CG test:\nresidu= " << res << std::endl;
    BOOST_CHECK(res < 1e-6);// default iteration tol is 1e-8

    std::vector<double> y(N);
    mult(bob,x,y);
    for(int i=0;i<N;i++)
        {
        std::cout << y[i] << "; ref val= " << b[i] << std::endl;
        double result = sq(y[i] - b[i]);
        std::cout << "result(should be zero)= " << result << std::endl;
        BOOST_TEST( result == 0.0 );
        }
    }

BOOST_AUTO_TEST_SUITE_END()

