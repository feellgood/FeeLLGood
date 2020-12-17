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
BOOST_AUTO_TEST_CASE(Tet_calc_vol, * boost::unit_test::tolerance(UT_TOL))
{
int nbNod = 4;
std::shared_ptr<Nodes::Node[]> node = std::shared_ptr<Nodes::Node[]>(new Nodes::Node[nbNod],std::default_delete<Nodes::Node[]>() ); 

Pt::pt3D p0(0,0,0),p1(1,0,0),p2(0,1,0),p3(0,0,1),u0(0,0,0),v0(0,0,0),u(0,0,0),v(0,0,0);
double theta_sph(0),phi_sph(0),phi0(0),phi(0),phiv0(0),phiv(0);

Nodes::Node n0 = {p0,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv};
Nodes::Node n1 = {p1,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv};
Nodes::Node n2 = {p2,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv};
Nodes::Node n3 = {p3,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv};

node[0] = n0;
node[1] = n1;
node[2] = n2;
node[3] = n3;

Tetra::Tet t(node,nbNod,0,0,1,2,3,4);//carefull with indices (starting from 1)
double result = 1/6.0;
double vol = t.calc_vol();
std::cout << "vol(tetra) =" << vol << std::endl;
BOOST_TEST( vol == result );
}

BOOST_AUTO_TEST_SUITE_END()
