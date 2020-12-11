#define BOOST_TEST_MODULE facetteTest

#include <boost/test/unit_test.hpp>

#include "facette.h"
#include "node.h"

BOOST_AUTO_TEST_SUITE(ut_facette)

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

BOOST_AUTO_TEST_CASE(Fac_constructor1)
{
std::vector<Nodes::Node> v_nodes;  

Facette::Fac f(0);

BOOST_CHECK( (f.getN() == Facette::N) && (f.getNPI() == Facette::NPI) );
}

BOOST_AUTO_TEST_CASE(Fac_constructor2)
{
std::vector<Nodes::Node> v_nodes;

Facette::Fac f(&v_nodes,0,0,0,0,0);

BOOST_CHECK( (f.getN() == Facette::N) && (f.getNPI() == Facette::NPI) );
}


/*---------------------------------------*/
/* first lvl tests : nested calculus,... */
/*---------------------------------------*/


/*---------------------------------------*/
/* second lvl tests : pure mathematics   */
/*---------------------------------------*/


BOOST_AUTO_TEST_SUITE_END()
