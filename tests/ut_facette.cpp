#define BOOST_TEST_MODULE facetteTest

#include <boost/test/unit_test.hpp>

#include <random>

#include "config.h"
#include "facette.h"
#include "node.h"
#include "tiny.h"

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
std::cout << "1 param constructor" << std::endl;
Facette::Fac f(0);

std::cout << "indices:" << f.ind[0] << ";" << f.ind[1] <<";" << f.ind[2] << std::endl;
BOOST_CHECK( (f.getN() == Facette::N) && (f.getNPI() == Facette::NPI) );
}

BOOST_AUTO_TEST_CASE(Fac_full_constructor)
{
std::cout << "7 param constructor" << std::endl;
Facette::Fac f(nullptr,0,0,0,0,0,0);

std::cout << "indices:" << f.ind[0] << ";" << f.ind[1] <<";" << f.ind[2] << std::endl;
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
std::cout << "calc_surf test" << std::endl;
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

Facette::Fac f(node,nbNod,0,0,1,2,3);// carefull with the index shift

std::cout << "indices:" << f.ind[0] << ";" << f.ind[1] <<";" << f.ind[2] << std::endl;

double s = 0.5;
BOOST_TEST( f.calc_surf() == s );
}

BOOST_AUTO_TEST_CASE(Fac_interpolation, * boost::unit_test::tolerance(UT_TOL))
{
std::cout << "surf interpolation test" << std::endl;
int nbNod = 3;
std::shared_ptr<Nodes::Node[]> node = std::shared_ptr<Nodes::Node[]>(new Nodes::Node[nbNod],std::default_delete<Nodes::Node[]>() ); 

std::random_device rd;
std::mt19937 gen(rd());// random number generator: standard Mersenne twister initialized with seed rd()
std::uniform_real_distribution<> distrib(0.0,1.0);


Pt::pt3D p1(1,0,0),p2(0,1,0),p3(1,1,0),u0(0,0,0),v0(0,0,0),u(0,0,0),v(0,0,0);
double theta_sph(0),phi_sph(0),phi0(0),phi(0),phiv0(0),phiv(0);

Nodes::Node n1 = {p1,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv};
Nodes::Node n2 = {p2,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv};
Nodes::Node n3 = {p3,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv};

node[0] = n1;
node[1] = n2;
node[2] = n3;

node[0].u0 = Pt::pt3D(M_PI*distrib(gen),2*M_PI*distrib(gen));
node[1].u0 = Pt::pt3D(M_PI*distrib(gen),2*M_PI*distrib(gen));
node[2].u0 = Pt::pt3D(M_PI*distrib(gen),2*M_PI*distrib(gen));

Facette::Fac f(node,nbNod,0,0,1,2,3);// carefull with the index shift

Pt::pt3D _u[Facette::NPI];
f.interpolation<Pt::pt3D>(Nodes::get_u0,_u);

double vec_nod[Pt::DIM][Facette::N];
for (int i=0;i<Pt::DIM;i++)
    for (int j=0;j<Facette::N;j++)
        { vec_nod[i][j] = node[j].u0(i); } // hidden transposition here
double result[Pt::DIM][Facette::NPI];
// a[N][NPI]
tiny::mult<double, Pt::DIM, Facette::N, Facette::NPI> (vec_nod, Facette::a, result);

double diff_r = 0;

for (int i=0;i<Pt::DIM;i++)
    for (int j=0;j<Facette::NPI;j++)
        diff_r += Pt::sq(result[i][j] - _u[j](i));

std::cout << "raw difference result =" << diff_r << std::endl;
BOOST_TEST( sqrt(diff_r) == 0.0 );
}


BOOST_AUTO_TEST_SUITE_END()
