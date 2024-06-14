#define BOOST_TEST_MODULE regionTest

#include <boost/test/unit_test.hpp>
#include<string>
#include <iostream>

#include "region.h"
#include "ut_config.h"  // for tolerance UT_TOL macro

#include "tetra.h"

BOOST_AUTO_TEST_SUITE(ut_region)

BOOST_AUTO_TEST_CASE(region_getters)
    {
    const int nb(42);
    std::string my_name("my_region");
    region<Tetra::Tet,Tetra::prm> bob(my_name,nb,LLG|ELECTROSTAT);
    bob.infos();
    BOOST_CHECK(bob.getName() == my_name);
    BOOST_CHECK(bob.getIndex() == nb);
    }

BOOST_AUTO_TEST_CASE(region_lalala, *boost::unit_test::tolerance(UT_TOL))
    {
    BOOST_TEST( (1.0+0.5*UT_TOL) == 1.0);
    }

BOOST_AUTO_TEST_SUITE_END()
