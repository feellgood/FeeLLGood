#define BOOST_TEST_MODULE tinyTest

#include <boost/test/unit_test.hpp>

#include "tiny.h"

BOOST_AUTO_TEST_SUITE(ut_tiny)

/*---------------------------------------*/
/* minus one test: check boost is fine   */
/*---------------------------------------*/

BOOST_AUTO_TEST_CASE(Stupid)
{
float x = 1.0;
BOOST_CHECK(x != 0.0f);
}

/*-----------------------------------------------------*/
/* zero lvl tests : direct elementary member functions */ 
/*-----------------------------------------------------*/

BOOST_AUTO_TEST_CASE(tiny_sub)
{
const int M = 3;
const int N = 2;
double A[M][N] = {{1,2},{3,4},{5,6}};
double B[M][N] = {{4,31},{-5,0},{0,1}};
double C[M][N] = {{0}};
const double val_ref[M][N] = {{-3,-29},{8,4},{5,5}};
tiny::sub<double,M,N>(A,B,C);

bool result = (C[0][0] == val_ref[0][0]);
result &= (C[0][1] == val_ref[0][1]);
result &= (C[1][0] == val_ref[1][0]);
result &= (C[1][1] == val_ref[1][1]);
result &= (C[2][0] == val_ref[2][0]);
result &= (C[2][1] == val_ref[2][1]);

std::cout << "basic test substraction" << std::endl;
BOOST_CHECK( result == true );
}

/*---------------------------------------*/
/* first lvl tests : nested calculus,... */
/*---------------------------------------*/


/*---------------------------------------*/
/* second lvl tests : pure mathematics   */
/*---------------------------------------*/

BOOST_AUTO_TEST_CASE(tiny_dist, * boost::unit_test::tolerance(1e-15))
{
const int M = 3;
const int N = 2;
double A[M][N] = {{1,2},{3,4},{5,6}};
double B[M][N] = {{4,31},{-5,0},{1e-17,1}};
double C[M][N] = {{0}};
const double val_ref[M][N] = {{-3,-29},{8,4},{5,5}};
tiny::sub<double,M,N>(A,B,C);

double result = tiny::dist<double,M,N>(C,val_ref);

std::cout << "test dist = frobenius_norm(A-B)" << std::endl;
BOOST_TEST( result == 0.0 );
}


BOOST_AUTO_TEST_SUITE_END()
