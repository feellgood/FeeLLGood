#define BOOST_TEST_MODULE geometryTest

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <ctime>

#include "pt3D.h"

BOOST_AUTO_TEST_SUITE(ut_pt3D_arithmetic)

/* minus one test: check boost is fine */

BOOST_AUTO_TEST_CASE(Stupid)
{
float x = 1.0;
BOOST_CHECK(x != 0.0f);
}

/* zero lvl tests : direct elementary member functions */ 

BOOST_AUTO_TEST_CASE(pt3D_getter)
{
Pt::pt3D X(1,2,3);

double result = X(Pt::IDX_X)+X(Pt::IDX_Y)*X(Pt::IDX_Z);
double x = 7.0;
BOOST_CHECK(x == result);
}

BOOST_AUTO_TEST_CASE(pt3D_op_add)
{
Pt::pt3D X(1,2,3);
Pt::pt3D Y(3,2,-4);

Pt::pt3D result = X+Y;
double x = 4.0;
double y = 4.0;
double z = -1.0;
BOOST_CHECK( (x == result(Pt::IDX_X)) && (y == result(Pt::IDX_Y)) && (z == result(Pt::IDX_Z)) );
}

BOOST_AUTO_TEST_CASE(pt3D_op_sub)
{
Pt::pt3D X(1,2,3);
Pt::pt3D Y(3,2,-4);

Pt::pt3D result = X-Y;
double x = -2.0;
double y = 0.0;
double z = 7.0;
BOOST_CHECK( (x == result(Pt::IDX_X)) && (y == result(Pt::IDX_Y)) && (z == result(Pt::IDX_Z)) );
}

BOOST_AUTO_TEST_CASE(pt3D_op_mult)
{
Pt::pt3D X(1,0,0);
Pt::pt3D Y(0,1,0);

Pt::pt3D result = X*Y;
double x = 0.0;
double y = 0.0;
double z = 1.0;
BOOST_CHECK( (x == result(Pt::IDX_X)) && (y == result(Pt::IDX_Y)) && (z == result(Pt::IDX_Z)) );
}

BOOST_AUTO_TEST_CASE(pt3D_op_pScal)
{
Pt::pt3D X(1,2,3);
Pt::pt3D Y(3,2,-4);

double result = Pt::pScal(X,Y);
double x = -5.0;
BOOST_CHECK( x == result );
}

BOOST_AUTO_TEST_CASE(pt3D_op_pDirect)
{
Pt::pt3D X(1,2,3);
Pt::pt3D Y(3,2,-4);

Pt::pt3D result = Pt::pDirect(X,Y);
double x = 3.0;
double y = 4.0;
double z = -12.0;
BOOST_CHECK( (x == result(Pt::IDX_X))&&(y == result(Pt::IDX_Y))&&(z == result(Pt::IDX_Z)) );
}

/* first lvl tests : nested calculus,... */

BOOST_AUTO_TEST_CASE(pt3D_product_left_and_right)
{
Pt::pt3D X(1,3,5);
Pt::pt3D Y(-2,4,1);

Pt::pt3D r1 = 3*X;
double x1 = (r1 - X*3).norm();
Pt::pt3D r2 = 0.5*Y;

double x2 = 2.0;
BOOST_CHECK( (x1 == 0.0) && (x2 == r2.maxLength())  );
}

/* second lvl tests : pure mathematics */

BOOST_AUTO_TEST_CASE(pt3D_triple_product)
{
srand(time(0));
const Pt::pt3D X(((double)rand())/RAND_MAX,((double)rand())/RAND_MAX,((double)rand())/RAND_MAX);
const Pt::pt3D Y(((double)rand())/RAND_MAX,((double)rand())/RAND_MAX,((double)rand())/RAND_MAX);
const Pt::pt3D Z((double)(rand())/RAND_MAX,((double)rand())/RAND_MAX,((double)rand())/RAND_MAX);

double x1 = Pt::pTriple(X,Y,Z);
double x2 = Pt::pTriple(Y,Z,X);
double x3 = Pt::pTriple(Z,X,Y);

std::cout << "this nasty test check that triple product of 3D vectors is invariant by circular permutation, whatever is the numerical error" << std::endl;
std::cout << "X=" << X << std::endl;
std::cout << "Y=" << Y << std::endl;
std::cout << "Z=" << Z << std::endl;

std::cout << "pTriple(X,Y,Z)" << x1 << std::endl;
std::cout << "pTriple(Y,Z,X)" << x2 << std::endl;
std::cout << "pTriple(Z,X,Y)" << x3 << std::endl;

std::cout << "x1-x2=" << x1-x2 << std::endl;
std::cout << "x2-x3=" << x2-x3 << std::endl;
std::cout << "x3-x1=" << x3-x1 << std::endl;

std::cout << "__DBL_EPSILON__ =" << __DBL_EPSILON__ << std::endl;

BOOST_CHECK( (fabs(x1-x2)  < __DBL_EPSILON__) && (fabs(x2-x3)  < __DBL_EPSILON__) && ( fabs(x3-x1)  < __DBL_EPSILON__)  );
}




BOOST_AUTO_TEST_SUITE_END()
