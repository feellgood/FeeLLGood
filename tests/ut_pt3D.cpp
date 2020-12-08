#define BOOST_TEST_MODULE geometryTest

#include <boost/test/unit_test.hpp>

#include "pt3D.h"

BOOST_AUTO_TEST_SUITE(ut_pt3D_arithmetic)

BOOST_AUTO_TEST_CASE(Stupid)
{
float x = 1.0;
BOOST_CHECK(x != 0.0f);
}

BOOST_AUTO_TEST_CASE(pt3D_getter)
{
Pt::pt3D X(1,2,3);

double result = X(Pt::IDX_X)+X(Pt::IDX_Y)*X(Pt::IDX_Z);
double x = 7.0;
BOOST_CHECK(x == result);
}

BOOST_AUTO_TEST_SUITE_END()
