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
Facette::Fac f(0);

BOOST_CHECK( (f.getN() == Facette::N) && (f.getNPI() == Facette::NPI) );
}

BOOST_AUTO_TEST_CASE(Fac_full_constructor)
{
Facette::Fac f(nullptr,0,0,0,0,0,0);

BOOST_CHECK( (f.getN() == Facette::N) && (f.getNPI() == Facette::NPI) );
}



BOOST_AUTO_TEST_CASE(Fac_constructor3)
{
Facette::Fac f(0,0,0,0,0);

BOOST_CHECK( (f.getN() == Facette::N) && (f.getNPI() == Facette::NPI) );
}

/*---------------------------------------*/
/* first lvl tests : nested calculus,... */
/*---------------------------------------*/


/*---------------------------------------*/
/* second lvl tests : pure mathematics   */
/*---------------------------------------*/

BOOST_AUTO_TEST_CASE(Fac_calc_surf, * boost::unit_test::tolerance(UT_TOL))
{
int nbNod = 3;
std::shared_ptr<Nodes::Node[]> node = std::shared_ptr<Nodes::Node[]>(new Nodes::Node[nbNod],std::default_delete<Nodes::Node[]>() ); 

Pt::pt3D p1(1,0,0),p2(0,1,0),p3(1,1,0),u0(0,0,0),v0(0,0,0),u(0,0,0),v(0,0,0);
double theta_sph(0),phi_sph(0),phi0(0),phi(0),phiv0(0),phiv(0);

Nodes::Node n1 = {p1,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv};
Nodes::Node n2 = {p2,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv};
Nodes::Node n3 = {p3,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv};

node[0] = n1;
node[1] = n2;
node[2] = n3;

Facette::Fac f(node,nbNod,0,0,0,1,2);
double s = 0.5;
BOOST_CHECK( f.calc_surf() == s );
}


BOOST_AUTO_TEST_SUITE_END()
