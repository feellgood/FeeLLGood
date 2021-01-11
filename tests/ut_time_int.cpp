#define BOOST_TEST_MODULE time_integrationTest

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <random>

#include "config.h" // for tolerance UT_TOL macro

BOOST_AUTO_TEST_SUITE(ut_time_int)

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

/*---------------------------------------*/
/* first lvl tests : nested calculus,... */
/*---------------------------------------*/

/*---------------------------------------*/
/* second lvl tests : pure mathematics   */
/*---------------------------------------*/


BOOST_AUTO_TEST_SUITE_END()
