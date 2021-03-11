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
double theta_sph(0),phi_sph(0),phi0(0),phi(0),phiv0(0),phiv(0),V(0);

Nodes::Node n1 = {p1,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv,V};
Nodes::Node n2 = {p2,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv,V};
Nodes::Node n3 = {p3,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv,V};

node[0] = n1;
node[1] = n2;
node[2] = n3;

Facette::Fac f(node,nbNod,0,0,1,2,3);// carefull with the index shift

std::cout << "indices:" << f.ind[0] << ";" << f.ind[1] <<";" << f.ind[2] << std::endl;

double s = 0.5;
BOOST_TEST( f.calc_surf() == s );
}

BOOST_AUTO_TEST_CASE(Fac_interpolation_pt3D, * boost::unit_test::tolerance(UT_TOL))
{
std::cout << "surf interpolation test" << std::endl;
int nbNod = 3;
std::shared_ptr<Nodes::Node[]> node = std::shared_ptr<Nodes::Node[]>(new Nodes::Node[nbNod],std::default_delete<Nodes::Node[]>() ); 

std::random_device rd;
std::mt19937 gen(rd());// random number generator: standard Mersenne twister initialized with seed rd()
std::uniform_real_distribution<> distrib(0.0,1.0);


Pt::pt3D p1(1,0,0),p2(0,1,0),p3(1,1,0),u0(0,0,0),v0(0,0,0),u(0,0,0),v(0,0,0);
double theta_sph(0),phi_sph(0),phi0(0),phi(0),phiv0(0),phiv(0),V(0);

Nodes::Node n1 = {p1,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv,V};
Nodes::Node n2 = {p2,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv,V};
Nodes::Node n3 = {p3,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv,V};

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

BOOST_AUTO_TEST_CASE(Fac_interpolation_double, * boost::unit_test::tolerance(UT_TOL))
{
std::cout << "surf interpolation test" << std::endl;
int nbNod = 3;
std::shared_ptr<Nodes::Node[]> node = std::shared_ptr<Nodes::Node[]>(new Nodes::Node[nbNod],std::default_delete<Nodes::Node[]>() ); 

std::random_device rd;
std::mt19937 gen(rd());// random number generator: standard Mersenne twister initialized with seed rd()
std::uniform_real_distribution<> distrib(0.0,1.0);


Pt::pt3D p1(1,0,0),p2(0,1,0),p3(1,1,0),u,v,u0,v0;
double theta_sph(0),phi_sph(0),phi0(0),phi(0),phiv0(0),phiv(0),V(0);

Nodes::Node n1 = {p1,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv,V};
Nodes::Node n2 = {p2,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv,V};
Nodes::Node n3 = {p3,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv,V};

node[0] = n1;
node[1] = n2;
node[2] = n3;

node[0].phi0 = distrib(gen);
node[1].phi0 = distrib(gen);
node[2].phi0 = distrib(gen);

Facette::Fac f(node,nbNod,0,0,1,2,3);// carefull with the index shift

double _p[Facette::NPI];
f.interpolation<double>(Nodes::get_phi0,_p);

double scal_nod[Facette::N];
    for (int j=0;j<Facette::N;j++)
        { scal_nod[j] = node[j].phi0; }

double result[Facette::NPI];
// a[N][NPI]
tiny::transposed_mult<double, Facette::N, Facette::NPI> (scal_nod, Facette::a, result);

double diff_r = 0;

for (int j=0;j<Facette::NPI;j++)
    diff_r += Pt::sq(result[j] - _p[j]);

std::cout << "raw difference result =" << diff_r << std::endl;
BOOST_TEST( sqrt(diff_r) == 0.0 );
}

BOOST_AUTO_TEST_CASE(Fac_potential_u, * boost::unit_test::tolerance(UT_TOL))
{
std::cout << "fac potential test on u" << std::endl;
int nbNod = 3;
std::shared_ptr<Nodes::Node[]> node = std::shared_ptr<Nodes::Node[]>(new Nodes::Node[nbNod],std::default_delete<Nodes::Node[]>() ); 

std::random_device rd;
std::mt19937 gen(rd());// random number generator: standard Mersenne twister initialized with seed rd()
std::uniform_real_distribution<> distrib(0.0,1.0);

Pt::pt3D p1(1,0,0),p2(0,1,0),p3(1,1,0);
Pt::pt3D u0(M_PI*distrib(gen),2*M_PI*distrib(gen)), u(M_PI*distrib(gen),2*M_PI*distrib(gen));
Pt::pt3D v0(2*distrib(gen)-1,2*distrib(gen)-1,2*distrib(gen)-1), v(2*distrib(gen)-1,2*distrib(gen)-1,2*distrib(gen)-1);

double theta_sph(0),phi_sph(0),phi0(0),phi(0),phiv0(0),phiv(0),V(0);

Nodes::Node n1 = {p1,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv,V};
Nodes::Node n2 = {p2,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv,V};
Nodes::Node n3 = {p3,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv,V};

node[0] = n1;
node[1] = n2;
node[2] = n3;

Facette::Fac f(node,nbNod,0,0,1,2,3);// carefull with the index shift
f.Ms = distrib(gen);

int i =0;
int Hv = 0;
// ref code (from MuMag_potential.cc with minimal adaptations of src_Tube_scalfmm_thiaville_ec_mu_oersted_thiele_dyn20180903.tgz )
double nx,ny,nz,Ms;

Pt::pt3D f_norm = f.calc_norm();
nx=f_norm.x();  ny=f_norm.y();  nz=f_norm.z(); Ms=f.Ms;

int ii  = (i+1)%3;
int iii = (i+2)%3;

int i_,ii_,iii_;
i_=f.ind[i];  ii_=f.ind[ii];  iii_=f.ind[iii];

double x1,x2,x3, y1,y2,y3, z1,z2,z3;

Nodes::Node &node1 = node[i_];
Nodes::Node &node2 = node[ii_];
Nodes::Node &node3 = node[iii_];

x1=node1.p.x();  x2=node2.p.x();  x3=node3.p.x();
y1=node1.p.y();  y2=node2.p.y();  y3=node3.p.y();
z1=node1.p.z();  z2=node2.p.z();  z3=node3.p.z();

double b,t,h;
b = sqrt( Pt::sq(x2-x1) + Pt::sq(y2-y1) + Pt::sq(z2-z1) );
t = (x2-x1)*(x3-x1) + (y2-y1)*(y3-y1) + (z2-z1)*(z3-z1);
h = 2.*f.surf;
t/=b;  h/=b;
double a = t/h;  double c = (t-b)/h;

double s1, s2, s3;
if (Hv) {
	s1 = node1.v(0)*nx + node1.v(1)*ny + node1.v(2)*nz;
	s2 = node2.v(0)*nx + node2.v(1)*ny + node2.v(2)*nz;
	s3 = node3.v(0)*nx + node3.v(1)*ny + node3.v(2)*nz;
   }
else {
	s1 = node1.u(0)*nx + node1.u(1)*ny + node1.u(2)*nz;
	s2 = node2.u(0)*nx + node2.u(1)*ny + node2.u(2)*nz;
	s3 = node3.u(0)*nx + node3.u(1)*ny + node3.u(2)*nz;
   }

double l = s1;
double j = (s2-s1)/b;
double k = t/b/h*(s1-s2) + (s3-s1)/h;

double cc1 = c*c+1;
double r = sqrt(h*h + (c*h+b)*(c*h+b));
double ll = log( (cc1*h + c*b + sqrt(cc1)*r) / (b*(c+sqrt(cc1))) );

double pot1, pot2, pot3, pot;
pot1 = b*b/pow(cc1,1.5)*ll + c*b*r/cc1 + h*r - c*b*b/cc1 - sqrt(a*a+1)*h*h;
pot1*= j/2.;

pot2 = -c*b*b/pow(cc1,1.5)*ll + b*r/cc1 - h*h/2. + h*h*log(c*h+b+r) - b*b/cc1;
pot2*= k/2.;

pot3 = h*log(c*h+b+r) - h + b/sqrt(cc1)*ll;
pot3*= l;

pot = pot1 + pot2 + pot3 + h*(k*h/2.+l)*(1-log(h*(a+sqrt(a*a+1)))) - k*h*h/4.;

double result_ref = Ms*pot;
// end ref code

double result_to_test = f.potential(Nodes::get_u,i);
std::cout << "raw difference result =" << result_to_test - result_ref << std::endl;
BOOST_TEST( result_to_test == result_ref );
}

BOOST_AUTO_TEST_CASE(Fac_potential_v, * boost::unit_test::tolerance(UT_TOL))
{
std::cout << "fac potential test on v" << std::endl;
int nbNod = 3;
std::shared_ptr<Nodes::Node[]> node = std::shared_ptr<Nodes::Node[]>(new Nodes::Node[nbNod],std::default_delete<Nodes::Node[]>() ); 

std::random_device rd;
std::mt19937 gen(rd());// random number generator: standard Mersenne twister initialized with seed rd()
std::uniform_real_distribution<> distrib(0.0,1.0);

Pt::pt3D p1(1,0,0),p2(0,1,0),p3(1,1,0);
Pt::pt3D u0(M_PI*distrib(gen),2*M_PI*distrib(gen)), u(M_PI*distrib(gen),2*M_PI*distrib(gen));
Pt::pt3D v0(2*distrib(gen)-1,2*distrib(gen)-1,2*distrib(gen)-1), v(2*distrib(gen)-1,2*distrib(gen)-1,2*distrib(gen)-1);

double theta_sph(0),phi_sph(0),phi0(0),phi(0),phiv0(0),phiv(0),V(0);

Nodes::Node n1 = {p1,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv,V};
Nodes::Node n2 = {p2,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv,V};
Nodes::Node n3 = {p3,u0,v0,u,v,theta_sph,phi_sph,phi0,phi,phiv0,phiv,V};

node[0] = n1;
node[1] = n2;
node[2] = n3;

Facette::Fac f(node,nbNod,0,0,1,2,3);// carefull with the index shift
f.Ms = distrib(gen);

int i =0;
int Hv = 1;
// ref code (from MuMag_potential.cc with minimal adaptations of src_Tube_scalfmm_thiaville_ec_mu_oersted_thiele_dyn20180903.tgz )
double nx,ny,nz,Ms;

Pt::pt3D f_norm = f.calc_norm();
nx=f_norm.x();  ny=f_norm.y();  nz=f_norm.z(); Ms=f.Ms;

int ii  = (i+1)%3;
int iii = (i+2)%3;

int i_,ii_,iii_;
i_=f.ind[i];  ii_=f.ind[ii];  iii_=f.ind[iii];

double x1,x2,x3, y1,y2,y3, z1,z2,z3;

Nodes::Node &node1 = node[i_];
Nodes::Node &node2 = node[ii_];
Nodes::Node &node3 = node[iii_];

x1=node1.p.x();  x2=node2.p.x();  x3=node3.p.x();
y1=node1.p.y();  y2=node2.p.y();  y3=node3.p.y();
z1=node1.p.z();  z2=node2.p.z();  z3=node3.p.z();

double b,t,h;
b = sqrt( Pt::sq(x2-x1) + Pt::sq(y2-y1) + Pt::sq(z2-z1) );
t = (x2-x1)*(x3-x1) + (y2-y1)*(y3-y1) + (z2-z1)*(z3-z1);
h = 2.*f.surf;
t/=b;  h/=b;
double a = t/h;  double c = (t-b)/h;

double s1, s2, s3;
if (Hv) {
	s1 = node1.v(0)*nx + node1.v(1)*ny + node1.v(2)*nz;
	s2 = node2.v(0)*nx + node2.v(1)*ny + node2.v(2)*nz;
	s3 = node3.v(0)*nx + node3.v(1)*ny + node3.v(2)*nz;
   }
else {
	s1 = node1.u(0)*nx + node1.u(1)*ny + node1.u(2)*nz;
	s2 = node2.u(0)*nx + node2.u(1)*ny + node2.u(2)*nz;
	s3 = node3.u(0)*nx + node3.u(1)*ny + node3.u(2)*nz;
   }

double l = s1;
double j = (s2-s1)/b;
double k = t/b/h*(s1-s2) + (s3-s1)/h;

double cc1 = c*c+1;
double r = sqrt(h*h + (c*h+b)*(c*h+b));
double ll = log( (cc1*h + c*b + sqrt(cc1)*r) / (b*(c+sqrt(cc1))) );

double pot1, pot2, pot3, pot;
pot1 = b*b/pow(cc1,1.5)*ll + c*b*r/cc1 + h*r - c*b*b/cc1 - sqrt(a*a+1)*h*h;
pot1*= j/2.;

pot2 = -c*b*b/pow(cc1,1.5)*ll + b*r/cc1 - h*h/2. + h*h*log(c*h+b+r) - b*b/cc1;
pot2*= k/2.;

pot3 = h*log(c*h+b+r) - h + b/sqrt(cc1)*ll;
pot3*= l;

pot = pot1 + pot2 + pot3 + h*(k*h/2.+l)*(1-log(h*(a+sqrt(a*a+1)))) - k*h*h/4.;

double result_ref = Ms*pot;
// end ref code

double result_to_test = f.potential(Nodes::get_v,i);
std::cout << "raw difference result =" << result_to_test - result_ref << std::endl;
BOOST_TEST( result_to_test == result_ref );
}


BOOST_AUTO_TEST_SUITE_END()
