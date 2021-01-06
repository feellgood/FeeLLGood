#define BOOST_TEST_MODULE tetraTest

#include <boost/test/unit_test.hpp>

#include <random>

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

BOOST_AUTO_TEST_CASE(Tet_nod_interpolation, * boost::unit_test::tolerance(UT_TOL))
{
int nbNod = 4;
std::shared_ptr<Nodes::Node[]> node = std::shared_ptr<Nodes::Node[]>(new Nodes::Node[nbNod],std::default_delete<Nodes::Node[]>() ); 

std::random_device rd;
std::mt19937 gen(rd());// random number generator: standard Mersenne twister initialized with seed rd()
std::uniform_real_distribution<> distrib(0.0,1.0);

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

for (int i=0;i<nbNod;i++) { node[i].u0 = Pt::pt3D(M_PI*distrib(gen),2*M_PI*distrib(gen)); }

Tetra::Tet t(node,nbNod,0,0,1,2,3,4);//carefull with indices (starting from 1)

// ref code (with minimal adaptations of integrales method in file MuMag_integrales.cc of src_Tube_scalfmm_thiaville_ec_mu_oersted_thiele_dyn20180903.tgz )
double _u_nod[3][Tetra::N], _u[3][Tetra::NPI];

for (int ie=0; ie<Tetra::N; ie++)
    {
    int i= t.ind[ie];
    Nodes::Node &nod = node[i];
    for (int d=0; d<3; d++) { _u_nod[d][ie]   = nod.u0(d); }
    }
tiny::mult<double, 3, Tetra::N, Tetra::NPI> (_u_nod, Tetra::a, _u);
// end ref code


// code to check
Pt::pt3D dUdx[Tetra::NPI], dUdy[Tetra::NPI], dUdz[Tetra::NPI];
Pt::pt3D U[Tetra::NPI];

t.interpolation(Nodes::get_u0,U,dUdx,dUdy,dUdz);
// end code to check

double distance = 0.0;
double n_u = 0.0;
double n_U = 0.0;
for(int i=0;i<Tetra::N;i++)
    for(int j=0;j<Pt::DIM;j++)
        {
        distance += Pt::sq( _u[j][i] - U[i](j) );
        n_u += Pt::sq( _u[j][i] );
        n_U += Pt::sq( U[i](j) );
        }
    
std::cout << "distance^2 =" << distance << std::endl;
BOOST_TEST( sqrt(distance) == 0.0 );
std::cout << "sum of the square of the matrix components of ref code =" << n_u << std::endl;
std::cout << "sum of the square of the matrix components of code to test =" << n_U << std::endl;
BOOST_TEST( n_u == n_U ); //let's be paranoid
}

BOOST_AUTO_TEST_SUITE_END()
