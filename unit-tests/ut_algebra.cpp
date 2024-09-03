#define BOOST_TEST_MODULE algebraTest

#include <boost/test/unit_test.hpp>

#include <iostream>

#include "ut_config.h"  // for tolerance UT_TOL macro
#include "../algebra/algebra.h"
#include "../algebra/cg.h"
#include "../algebra/bicg.h"

using namespace algebra;
inline double sq(const double x) { return x * x; } // same as Nodes::sq

BOOST_AUTO_TEST_SUITE(ut_algebra)

/**
test on scaled, dot, norm, sub functions
*/
BOOST_AUTO_TEST_CASE(test_scaled, *boost::unit_test::tolerance(UT_TOL))
    {
    std::vector<double> x {1,4,-2,sqrt(42)};
    double alpha(2.0);
    std::vector<double> y(4);
    scaled(x,alpha,y); // y = alpha*x
    
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

/** test on v_coeff constructor, getVal and operators == < */
BOOST_AUTO_TEST_CASE(test_v_coeff, *boost::unit_test::tolerance(UT_TOL))
    {
    v_coeff bob(2,sqrt(2));
    
    BOOST_TEST( bob.getVal() == sqrt(2) );
    BOOST_TEST( bob._i == 2 );

    bob.add(3.14);
    BOOST_CHECK( bob.getVal() == (3.14 + sqrt(2)) );

    v_coeff jeff(2,4.56);
    BOOST_TEST( bob == jeff );
    BOOST_TEST( (jeff < bob) == false );
    jeff._i = 1;
    BOOST_TEST( jeff < bob );
    }

/** elementary tests on w_sparseVect: constructor, insert, exist, getVal methods */
BOOST_AUTO_TEST_CASE(test_w_sparseVect, *boost::unit_test::tolerance(UT_TOL))
    {
    w_sparseVect bob;
    v_coeff jeff(2,4.56);
    bob.insert(jeff);
    BOOST_CHECK( bob.exist(2) );
    BOOST_CHECK( bob.exist(3) == false );
    double result = bob.getVal(2);
    BOOST_CHECK( result == 4.56 );
    BOOST_CHECK( bob.getVal(42) == 0.0 );
    }

/** more tests on w_sparseVect */
BOOST_AUTO_TEST_CASE(advanced_test_w_sparseVect, *boost::unit_test::tolerance(UT_TOL))
    {
    w_sparseVect x;
    v_coeff bob(2,sqrt(2));
    x.insert(bob);
    v_coeff jeff(2,4.56);
    x.insert(jeff);
    v_coeff cat(1,1.23);
    x.insert(cat);
    BOOST_CHECK( x.getVal(2) == (4.56 + sqrt(2)) );
    BOOST_TEST ( x.getVal(42) == 0 );

    x.insert(v_coeff(10,7.89));
    BOOST_TEST( x.exist(10) );
    BOOST_TEST( !(x.exist(11)) );

    x.insert(v_coeff(15,0));
    BOOST_TEST( x.exist(15) ); /// this is weird ...

    x.insert(v_coeff(-1,3.14));
    BOOST_TEST( x.exist(-1) ); /// this might be weird too ...
    }

/** tests on r_sparseVect */
BOOST_AUTO_TEST_CASE(test_r_sparseVect, *boost::unit_test::tolerance(UT_TOL))
    {
    w_sparseVect x;
    v_coeff bob(2,sqrt(2));
    x.insert(bob);
    v_coeff jeff(2,4.56);
    x.insert(jeff);
    v_coeff cat(1,1.23);
    x.insert(cat);
    x.insert(v_coeff(10,7.89));
    x.insert(v_coeff(15,0));
    r_sparseVect y(x);

    BOOST_TEST( !(y.exist(11)) );
    BOOST_TEST( y.exist(10) );
    BOOST_TEST( y.exist(15) ); //the coeff val is zero but not filtered out by w_sparseVect.insert method
    BOOST_CHECK( y.getVal(2) == (4.56 + sqrt(2)) );
    }

/** tests on r_sparseVect.dot */
BOOST_AUTO_TEST_CASE(test_r_sparseVect_dot, *boost::unit_test::tolerance(UT_TOL))
    {
    w_sparseVect x;
    v_coeff bob(2,sqrt(2));
    x.insert(bob);
    v_coeff jeff(2,4.56);
    x.insert(jeff);
    v_coeff cat(1,1.23);
    x.insert(cat);
    x.insert(v_coeff(10,7.89));
    x.insert(v_coeff(15,0));
    r_sparseVect y(x);

    std::vector<double> z {0,1,0.5};
    BOOST_CHECK( y.dot(z) == (1.23 + 0.5*(4.56 + sqrt(2))) );
    }

/** tests on (w|r)_sparseMat */
BOOST_AUTO_TEST_CASE(test_w_sparseMat, *boost::unit_test::tolerance(UT_TOL))
    {
    const int N=4;
    w_sparseMat m(N);
    BOOST_TEST( m.getDim() == N );
    m.insert( m_coeff(1,1,3.14) );
    m.insert( m_coeff(0,0,1.0) );
    m.insert( m_coeff(2,2,5.0) );
    m.insert( m_coeff(3,3,42.0) );
    m.insert( m_coeff(1,3,-10.0) );
    m.insert( m_coeff(1,3,10.0) );
    m.insert( m_coeff(0,3,0.5) );
    r_sparseMat bob(m);
    std::vector<double> x {1.0,1.0,1.0,1.0};
    std::vector<double> y(N);
    m.insert( m_coeff(2,2,5.0) );// m modification must not affect bob
    mult(bob,x,y); //y = m*x
    BOOST_CHECK( y[0] == 1.5 );
    BOOST_CHECK( y[1] == 3.14 );
    BOOST_CHECK( y[2] == 5.0 );
    BOOST_CHECK( y[3] == 42.0 );
    }

/** test on gradient conjugate algorithm */
BOOST_AUTO_TEST_CASE(test_cg, *boost::unit_test::tolerance(UT_TOL))
    {
    const int N=4;
    w_sparseMat m(N);
    BOOST_TEST( m.getDim() == N );
    m.insert( m_coeff(1,1,3.14) );
    m.insert( m_coeff(0,0,1.0) );
    m.insert( m_coeff(2,2,5.0) );
    m.insert( m_coeff(3,3,42.0) );
    m.insert( m_coeff(1,3,-10.0) );
    m.insert( m_coeff(1,3,10.0) );
    m.insert( m_coeff(0,3,0.5) );
    r_sparseMat bob(m);
    std::vector<double> b {1.0,1.0,1.0,1.0};
    std::vector<double> x(N);
    iteration algo_it;
    double res = cg(bob,x,b,algo_it);
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

/** test on stabilized bi-gradient conjugate algorithm */
BOOST_AUTO_TEST_CASE(test_bicg, *boost::unit_test::tolerance(UT_TOL))
    {
    const int N=4;
    w_sparseMat m(N);
    BOOST_TEST( m.getDim() == N );
    m.insert( m_coeff(1,1,3.14) );
    m.insert( m_coeff(0,0,1.0) );
    m.insert( m_coeff(2,2,5.0) );
    m.insert( m_coeff(3,3,42.0) );
    m.insert( m_coeff(1,3,-10.0) );
    m.insert( m_coeff(1,3,10.0) );
    m.insert( m_coeff(0,3,0.5) );
    r_sparseMat bob(m);
    std::vector<double> b {1.0,1.0,1.0,1.0};
    std::vector<double> x(N);
    iteration algo_it;
    double res = bicg(bob,x,b,algo_it);
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

