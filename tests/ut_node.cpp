#define BOOST_TEST_MODULE nodeTest

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <random>

#include "node.h"

BOOST_AUTO_TEST_SUITE(ut_node)

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

BOOST_AUTO_TEST_CASE(node_e_p)
{
std::random_device rd;

std::mt19937 gen(rd());// random number generator: standard Mersenne twister initialized with seed rd()
std::uniform_real_distribution<> distrib(0.0,M_PI);
    
Nodes::Node n;

n.theta_sph = distrib(gen);
n.phi_sph = 2*distrib(gen);
n.u0 = Pt::pt3D(1,3,5); // u0 should be a unit vector but the tests here check more
Pt::pt3D X = n.calc_ep();

BOOST_CHECK( (fabs(X.norm() - 1.0) < __DBL_EPSILON__)&&( fabs(Pt::pScal(n.u0,X)) < __DBL_EPSILON__  ) );
}


BOOST_AUTO_TEST_SUITE_END()
