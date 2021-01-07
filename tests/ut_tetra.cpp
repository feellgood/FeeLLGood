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
BOOST_AUTO_TEST_CASE(Tet_inner_tables, * boost::unit_test::tolerance(UT_TOL))
{
// this test is dedicated to  check dadx,dady,dadz and weight tables, those values are initilized once by Tet constructor
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

// ref code (with minimal adaptations of dad(x|y|z) in file Mesh_hat.cc of src_Tube_scalfmm_thiaville_ec_mu_oersted_thiele_dyn20180903.tgz )
double _dadx[Tetra::N][Tetra::NPI];
double _dady[Tetra::N][Tetra::NPI];
double _dadz[Tetra::N][Tetra::NPI];
double da[Tetra::N][Pt::DIM];
double J[Pt::DIM][Pt::DIM];
double nod[Pt::DIM][Tetra::N];
double weight[Tetra::NPI];

for (int ie=0; ie<Tetra::N; ie++)
    {
    int i= t.ind[ie];
    nod[0][ie] = node[i].p(0);
    nod[1][ie] = node[i].p(1);
    nod[2][ie] = node[i].p(2);
    }

tiny::mult<double,Pt::DIM,Tetra::N,Pt::DIM>(nod, Tetra::dadu, J);
double detJ = Pt::det(J);//lu_det(J);
Pt::inverse(J,detJ);
tiny::mult<double,Tetra::N,Pt::DIM,Pt::DIM>(Tetra::dadu, J, da);

for (int npi=0; npi<Tetra::NPI; npi++)
    {
    for (int ie=0; ie<Tetra::N; ie++)
        {
        _dadx[ie][npi]= da[ie][0];
        _dady[ie][npi]= da[ie][1];
        _dadz[ie][npi]= da[ie][2];
        }
    weight[npi]= detJ*Tetra::pds[npi];
    }

// end ref code

// the Tet constructor is computing dad(x|y|z)

double result_dadx(0.0),result_dady(0.0),result_dadz(0.0),result_w(0.0);
for (int npi=0; npi<Tetra::NPI; npi++)
    {
    for (int ie=0; ie<Tetra::N; ie++)
        {
        result_dadx += Pt::sq( _dadx[ie][npi] - t.dadx[ie][npi] );
        result_dady += Pt::sq( _dady[ie][npi] - t.dady[ie][npi] );
        result_dadz += Pt::sq( _dadz[ie][npi] - t.dadz[ie][npi] );
        }
    result_w += Pt::sq( weight[npi] - t.weight[npi] );
    }

std::cout << "sq_frob norm diff dadx,dady,dadz= " << result_dadx << " ; " << result_dady << " ; " << result_dadz << std::endl;
BOOST_TEST( sqrt(result_dadx) == 0.0 );
BOOST_TEST( sqrt(result_dady) == 0.0 );
BOOST_TEST( sqrt(result_dadz) == 0.0 );

std::cout<< "frob norm diff weight= "<< result_w <<std::endl;
BOOST_TEST( sqrt(result_w) == 0.0 );
}

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

double sq_dist(double _x[Pt::DIM][Tetra::NPI],Pt::pt3D X[Tetra::NPI])
{
double val(0.0);

for(int i=0;i<Tetra::N;i++)
    for(int j=0;j<Pt::DIM;j++)
        {val += Pt::sq( _x[j][i] - X[i](j) );}

return val;
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

for (int i=0;i<nbNod;i++)
    {
    node[i].u0 = Pt::pt3D(M_PI*distrib(gen),2*M_PI*distrib(gen));
    node[i].v0 = Pt::pt3D(M_PI*distrib(gen),2*M_PI*distrib(gen));
    }

Tetra::Tet t(node,nbNod,0,0,1,2,3,4);//carefull with indices (starting from 1)

// ref code (with minimal adaptations of integrales method in file MuMag_integrales.cc of src_Tube_scalfmm_thiaville_ec_mu_oersted_thiele_dyn20180903.tgz )
double _u_nod[3][Tetra::N], _u[3][Tetra::NPI];
double dudx[3][Tetra::NPI], dudy[3][Tetra::NPI], dudz[3][Tetra::NPI];

double _v_nod[3][Tetra::N], _v[3][Tetra::NPI];
double dvdx[3][Tetra::NPI], dvdy[3][Tetra::NPI], dvdz[3][Tetra::NPI];

for (int ie=0; ie<Tetra::N; ie++)
    {
    int i= t.ind[ie];
    Nodes::Node &nod = node[i];
    for (int d=0; d<3; d++)
        { 
        _u_nod[d][ie]   = nod.u0(d);
        _v_nod[d][ie]   = nod.v0(d);
        }
    }
tiny::mult<double, 3, Tetra::N, Tetra::NPI> (_u_nod, Tetra::a, _u);
tiny::mult<double, 3, Tetra::N, Tetra::NPI> (_u_nod, t.dadx, dudx);
tiny::mult<double, 3, Tetra::N, Tetra::NPI> (_u_nod, t.dady, dudy);
tiny::mult<double, 3, Tetra::N, Tetra::NPI> (_u_nod, t.dadz, dudz);

tiny::mult<double, 3, Tetra::N, Tetra::NPI> (_v_nod, Tetra::a, _v);
tiny::mult<double, 3, Tetra::N, Tetra::NPI> (_v_nod, t.dadx, dvdx);
tiny::mult<double, 3, Tetra::N, Tetra::NPI> (_v_nod, t.dady, dvdy);
tiny::mult<double, 3, Tetra::N, Tetra::NPI> (_v_nod, t.dadz, dvdz);
// end ref code


// code to check
Pt::pt3D dUdx[Tetra::NPI], dUdy[Tetra::NPI], dUdz[Tetra::NPI];
Pt::pt3D dVdx[Tetra::NPI], dVdy[Tetra::NPI], dVdz[Tetra::NPI];
Pt::pt3D U[Tetra::NPI],V[Tetra::NPI];

t.interpolation(Nodes::get_u0,U,dUdx,dUdy,dUdz);
t.interpolation(Nodes::get_v0,V,dVdx,dVdy,dVdz);
// end code to check

double n_u = tiny::frob_norm<double,3,Tetra::NPI>(_u);
double n_dudx = tiny::frob_norm<double,3,Tetra::NPI>(dudx);
double n_dudy = tiny::frob_norm<double,3,Tetra::NPI>(dudy);
double n_dudz = tiny::frob_norm<double,3,Tetra::NPI>(dudz);

double n_v = tiny::frob_norm<double,3,Tetra::NPI>(_v);
double n_dvdx = tiny::frob_norm<double,3,Tetra::NPI>(dvdx);
double n_dvdy = tiny::frob_norm<double,3,Tetra::NPI>(dvdy);
double n_dvdz = tiny::frob_norm<double,3,Tetra::NPI>(dvdz);

double dist_uU = sq_dist(_u,U);
double dist_dudx_dUdx = sq_dist(dudx,dUdx);
double dist_dudy_dUdy = sq_dist(dudy,dUdy);
double dist_dudz_dUdz = sq_dist(dudz,dUdz);

double dist_vV = sq_dist(_v,V);
double dist_dvdx_dVdx = sq_dist(dvdx,dVdx);
double dist_dvdy_dVdy = sq_dist(dvdy,dVdy);
double dist_dvdz_dVdz = sq_dist(dvdz,dVdz);

double n_U = sqrt(Pt::sq_frobenius_norm<Tetra::NPI>(U));
double n_dUdx = sqrt(Pt::sq_frobenius_norm<Tetra::NPI>(dUdx));
double n_dUdy = sqrt(Pt::sq_frobenius_norm<Tetra::NPI>(dUdy));
double n_dUdz = sqrt(Pt::sq_frobenius_norm<Tetra::NPI>(dUdz));

double n_V = sqrt(Pt::sq_frobenius_norm<Tetra::NPI>(V));
double n_dVdx = sqrt(Pt::sq_frobenius_norm<Tetra::NPI>(dVdx));
double n_dVdy = sqrt(Pt::sq_frobenius_norm<Tetra::NPI>(dVdy));
double n_dVdz = sqrt(Pt::sq_frobenius_norm<Tetra::NPI>(dVdz));

std::cout << "distance^2 (u,U) =" << dist_uU << std::endl;
BOOST_TEST( sqrt(dist_uU) == 0.0 );
BOOST_TEST( sqrt(dist_dudx_dUdx) == 0.0 );
BOOST_TEST( sqrt(dist_dudy_dUdy) == 0.0 );
BOOST_TEST( sqrt(dist_dudz_dUdz) == 0.0 );

std::cout << "distance^2 (v,V) =" << dist_vV << std::endl;
BOOST_TEST( sqrt(dist_vV) == 0.0 );
BOOST_TEST( sqrt(dist_dvdx_dVdx) == 0.0 );
BOOST_TEST( sqrt(dist_dvdy_dVdy) == 0.0 );
BOOST_TEST( sqrt(dist_dvdz_dVdz) == 0.0 );

// to avoid gag of comparing pure zeros we also check that matrices norm are equal
//let's be paranoid

std::cout << "frobenius norm of ref code u=" << n_u << std::endl;
std::cout << "frobenius norm of code to test U=" << n_U << std::endl;
BOOST_TEST( n_u == n_U );

std::cout << "frobenius norm of ref code dudx=" << n_dudx << std::endl;
std::cout << "frobenius norm of code to test dUdx=" << n_dUdx << std::endl;
BOOST_TEST( n_dudx == n_dUdx );

std::cout << "frobenius norm of ref code dudy=" << n_dudy << std::endl;
std::cout << "frobenius norm of code to test dUdy=" << n_dUdy << std::endl;
BOOST_TEST( n_dudy == n_dUdy );

std::cout << "frobenius norm of ref code dudz=" << n_dudz << std::endl;
std::cout << "frobenius norm of code to test dUdz=" << n_dUdz << std::endl;
BOOST_TEST( n_dudz == n_dUdz );

std::cout << "frobenius norm of ref code v=" << n_v << std::endl;
std::cout << "frobenius norm of code to test V=" << n_V << std::endl;
BOOST_TEST( n_v == n_V );

std::cout << "frobenius norm of ref code dvdx=" << n_dvdx << std::endl;
std::cout << "frobenius norm of code to test dVdx=" << n_dVdx << std::endl;
BOOST_TEST( n_dvdx == n_dVdx );

std::cout << "frobenius norm of ref code dvdy=" << n_dvdy << std::endl;
std::cout << "frobenius norm of code to test dVdy=" << n_dVdy << std::endl;
BOOST_TEST( n_dvdy == n_dVdy );

std::cout << "frobenius norm of ref code dvdz=" << n_dvdz << std::endl;
std::cout << "frobenius norm of code to test dVdz=" << n_dVdz << std::endl;
BOOST_TEST( n_dvdz == n_dVdz );
}

BOOST_AUTO_TEST_SUITE_END()
