#define BOOST_TEST_MODULE geometryTest

#include <boost/test/unit_test.hpp>

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

BOOST_AUTO_TEST_CASE(pt3D_triple_product)
{
Pt::pt3D X(1,3,5);
Pt::pt3D Y(-2,4,1);
Pt::pt3D Z();

Pt::pt3D r1 = 3*X;
double x1 = (r1 - X*3).norm();
Pt::pt3D r2 = 0.5*Y;

double x2 = 2.0;
BOOST_CHECK( (x1 == 0.0) && (x2 == r2.maxLength())  );
}




BOOST_AUTO_TEST_SUITE_END()
