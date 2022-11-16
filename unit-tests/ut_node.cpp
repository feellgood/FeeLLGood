#define BOOST_TEST_MODULE nodeTest

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <random>

#include "ut_config.h" // for tolerance UT_TOL macro
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

BOOST_AUTO_TEST_CASE(node_get_p_lvl0)
{
Nodes::Node n;
Pt::pt3D pPos(1.0,0.0,0.0);

n.p = pPos;

BOOST_CHECK(Nodes::get_p(n).x() == 1.0);
BOOST_CHECK(Nodes::get_p(n).y() == 0.0);
BOOST_CHECK(Nodes::get_p(n).z() == 0.0);
}


/*---------------------------------------*/
/* first lvl tests : nested calculus,... */
/*---------------------------------------*/

BOOST_AUTO_TEST_CASE(node_get_p_lvl1, * boost::unit_test::tolerance(UT_TOL))
{
Nodes::Node n;
Pt::pt3D pPos(1.0,3.0,5.0);

n.p = pPos;
n.p.normalize();

BOOST_TEST(Nodes::get_p(n).norm() == 1.0);
}

/*---------------------------------------*/
/* second lvl tests : pure mathematics   */
/*---------------------------------------*/

// Build a unit vector from the cylindrical coordinates (theta, z).
// If theta and z are uniformly distributed in [-pi, pi] and [-1, 1]
// respectively, the resulting vector is isotropically distributed.
static Pt::pt3D unit_vector(double theta, double z)
{
double r = sqrt(1 - z*z);
return Pt::pt3D(r * cos(theta), r * sin(theta), z);
}

BOOST_AUTO_TEST_CASE(node_e_p, * boost::unit_test::tolerance(10.0*UT_TOL))
{
unsigned sd = my_seed();
std::mt19937 gen(sd);
std::uniform_real_distribution<> distrib(-1.0,1.0);
    
Nodes::Node n;

// test the orthonormality of the basis (u0, ep, eq)
n.u0 = unit_vector(M_PI * distrib(gen), distrib(gen));
n.setBasis(M_PI * distrib(gen));

if (!DET_UT) std::cout << "seed =" << sd << std::endl;
BOOST_TEST(n.u0.norm() == 1.0);
BOOST_TEST(n.ep.norm() == 1.0);
BOOST_TEST(n.eq.norm() == 1.0);
BOOST_TEST(Pt::pScal(n.u0, n.ep) == 0.0);
BOOST_TEST(Pt::pScal(n.ep, n.eq) == 0.0);
BOOST_TEST(Pt::pScal(n.eq, n.u0) == 0.0);
}

BOOST_AUTO_TEST_CASE(node_evol, * boost::unit_test::tolerance(1e3*UT_TOL))
{
unsigned sd = my_seed();
std::mt19937 gen(sd);
std::uniform_real_distribution<> distrib(-1.0,1.0);
    
Nodes::Node n;

double vp = distrib(gen);
double vq = distrib(gen);
double dt = distrib(gen) + 1.0;

n.u0 = unit_vector(M_PI * distrib(gen), distrib(gen));
n.setBasis(M_PI * distrib(gen));

Pt::pt3D v = vp*n.ep + vq*n.eq;

n.make_evol(vp,vq,dt);

if (!DET_UT) std::cout << "seed =" << sd << std::endl;
std::cout << "v = " << v << std::endl;
std::cout << "simplify[v] = " << n.v << std::endl;

BOOST_TEST( Pt::dist(n.v,v)  == 0.0 );
}


BOOST_AUTO_TEST_SUITE_END()
