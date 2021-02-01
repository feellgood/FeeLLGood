#define BOOST_TEST_MODULE anisotropyTest

#include <boost/test/unit_test.hpp>

#include <random>

#include "tetra.h"

BOOST_AUTO_TEST_SUITE(ut_anisotropy)

BOOST_AUTO_TEST_CASE(Tet_aniso_uniax, * boost::unit_test::tolerance(10.0*UT_TOL))
{
int nbNod = 4;
std::shared_ptr<Nodes::Node[]> node = std::shared_ptr<Nodes::Node[]>(new Nodes::Node[nbNod],std::default_delete<Nodes::Node[]>() ); 

std::random_device rd;
std::mt19937 gen(rd());// random number generator: standard Mersenne twister initialized with seed rd()
std::uniform_real_distribution<> distrib(0.0,1.0);

Pt::pt3D p0(0,0,0),p1(1,0,0),p2(0,1,0),p3(0,0,1),u0(0,0,0),v0(0,0,0);
double theta_sph(0),phi_sph(0),phi0(0),phi(0),phiv0(0),phiv(0);

Nodes::Node n0 = {p0,u0,v0,Pt::pt3D(0,0,0),Pt::pt3D(0,0,0),theta_sph,phi_sph,phi0,phi,phiv0,phiv,0};
Nodes::Node n1 = {p1,u0,v0,Pt::pt3D(0,0,0),Pt::pt3D(0,0,0),theta_sph,phi_sph,phi0,phi,phiv0,phiv,0};
Nodes::Node n2 = {p2,u0,v0,Pt::pt3D(0,0,0),Pt::pt3D(0,0,0),theta_sph,phi_sph,phi0,phi,phiv0,phiv,0};
Nodes::Node n3 = {p3,u0,v0,Pt::pt3D(0,0,0),Pt::pt3D(0,0,0),theta_sph,phi_sph,phi0,phi,phiv0,phiv,0};

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

double dt = distrib(gen);

Pt::pt3D uk = Pt::pt3D(M_PI*distrib(gen),2*M_PI*distrib(gen));

double uk00 = uk.x();
double uk01 = uk.y();
double uk02 = uk.z();
double K = distrib(gen);
double Js = 0.5 + distrib(gen);// add 0.5 to center Js around 1

double contrib_aniso(0);
double Kbis=2.0*K/Js;
std::cout << "uniaxial anisotropy test on a tetrahedron" << std::endl;
std::cout << "uk =" << uk << std::endl;

// ref code (with minimal adaptations of integrales method in file MuMag_integrales.cc of src_Tube_scalfmm_thiaville_ec_mu_oersted_thiele_dyn20180903.tgz )
double u_nod[3][Tetra::N];
double v_nod[3][Tetra::N];
double u[3][Tetra::NPI];
double v[3][Tetra::NPI];

for (int ie=0; ie<Tetra::N; ie++){
    int i= t.ind[ie];
    Nodes::Node &nod = node[i];
    for (int d=0; d<3; d++)
        {
        u_nod[d][ie]   = nod.u0(d);
        v_nod[d][ie]   = nod.v0(d);     
        }
    }

tiny::mult<double, 3, Tetra::N, Tetra::NPI> (u_nod, Tetra::a, u);
tiny::mult<double, 3, Tetra::N, Tetra::NPI> (v_nod, Tetra::a, v);

Pt::pt3D U[Tetra::NPI];
Pt::pt3D V[Tetra::NPI];
for(int npi=0;npi<Tetra::NPI;npi++)
    {
    U[npi]=Pt::pt3D(u[0][npi],u[1][npi],u[2][npi]);
    V[npi]=Pt::pt3D(v[0][npi],v[1][npi],v[2][npi]);
    }// deep copy instead of re-computing

for (int npi=0; npi<Tetra::NPI; npi++)
    {
    double uk0_u = uk00*u[0][npi] + uk01*u[1][npi] + uk02*u[2][npi]; 
    double uk0_v = uk00*v[0][npi] + uk01*v[1][npi] + uk02*v[2][npi]; 

    double uHau = 2*K/Js* uk0_u*uk0_u;

    double Ht[3];
    Ht[0]= 2*K/Js* uk0_v*uk00;
    Ht[1]= 2*K/Js* uk0_v*uk01;
    Ht[2]= 2*K/Js* uk0_v*uk02;

    double H[3];
    H[0]=(2*K/Js* uk0_u*uk00);
    H[1]=(2*K/Js* uk0_u*uk01);
    H[2]=(2*K/Js* uk0_u*uk02);

    Pt::pt3D H_aniso;
    contrib_aniso = t.calc_aniso_uniax(npi,uk,Kbis,THETA*dt,U,V,H_aniso); // code to test

    Pt::pt3D refHval = Pt::pt3D(H[0]+THETA*dt*Ht[0],H[1]+THETA*dt*Ht[1],H[2]+THETA*dt*Ht[2]);
    std::cout << "ref H value = " << refHval << "; H_aniso=" << H_aniso << std::endl;
    BOOST_TEST( Pt::dist(refHval, H_aniso) == 0.0 );
    
    BOOST_TEST( uHau == contrib_aniso );
    }
}

BOOST_AUTO_TEST_SUITE_END()
