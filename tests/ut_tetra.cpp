#define BOOST_TEST_MODULE tetraTest

#include <boost/test/unit_test.hpp>

#include "tetra.h"

BOOST_AUTO_TEST_SUITE(ut_tetra)

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
BOOST_AUTO_TEST_CASE(Tet_constructor)
{
Tetra::Tet tet(nullptr,0,0,0,0,0,0,0);

BOOST_CHECK( (tet.getN() == Tetra::N) && (tet.getNPI() == Tetra::NPI) );
}


/*---------------------------------------*/
/* first lvl tests : nested calculus,... */
/*---------------------------------------*/


/*---------------------------------------*/
/* second lvl tests : pure mathematics   */
/*---------------------------------------*/


BOOST_AUTO_TEST_SUITE_END()
