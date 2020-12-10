#define BOOST_TEST_MODULE pt3DTest

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <random>

#include "pt3D.h"
#include "tiny.h"

BOOST_AUTO_TEST_SUITE(ut_pt3D_arithmetic)

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


BOOST_AUTO_TEST_CASE(pt3D_getter)
{
Pt::pt3D X(1,2,3);

double result = X(Pt::IDX_X)+X(Pt::IDX_Y)*X(Pt::IDX_Z);
double x = 7.0;
BOOST_CHECK(x == result);
}

BOOST_AUTO_TEST_CASE(pt3D_op_add)
{
Pt::pt3D X(1,2,3);
Pt::pt3D Y(3,2,-4);

Pt::pt3D result = X+Y;
double x = 4.0;
double y = 4.0;
double z = -1.0;
BOOST_CHECK( (x == result(Pt::IDX_X)) && (y == result(Pt::IDX_Y)) && (z == result(Pt::IDX_Z)) );
}

BOOST_AUTO_TEST_CASE(pt3D_op_sub)
{
Pt::pt3D X(1,2,3);
Pt::pt3D Y(3,2,-4);

Pt::pt3D result = X-Y;
double x = -2.0;
double y = 0.0;
double z = 7.0;
BOOST_CHECK( (x == result(Pt::IDX_X)) && (y == result(Pt::IDX_Y)) && (z == result(Pt::IDX_Z)) );
}

BOOST_AUTO_TEST_CASE(pt3D_op_mult)
{
Pt::pt3D X(1,0,0);
Pt::pt3D Y(0,1,0);

Pt::pt3D result = X*Y;
double x = 0.0;
double y = 0.0;
double z = 1.0;
BOOST_CHECK( (x == result(Pt::IDX_X)) && (y == result(Pt::IDX_Y)) && (z == result(Pt::IDX_Z)) );
}

BOOST_AUTO_TEST_CASE(pt3D_op_pScal)
{
Pt::pt3D X(1,2,3);
Pt::pt3D Y(3,2,-4);

double result = Pt::pScal(X,Y);
double x = -5.0;
BOOST_CHECK( x == result );
}

BOOST_AUTO_TEST_CASE(pt3D_op_pDirect)
{
Pt::pt3D X(1,2,3);
Pt::pt3D Y(3,2,-4);

Pt::pt3D result = Pt::pDirect(X,Y);
double x = 3.0;
double y = 4.0;
double z = -12.0;
BOOST_CHECK( (x == result(Pt::IDX_X))&&(y == result(Pt::IDX_Y))&&(z == result(Pt::IDX_Z)) );
}

BOOST_AUTO_TEST_CASE(pt3D_op_dist)
{
Pt::pt3D X(1,2,3);
Pt::pt3D Y(3,2,-4);

double x = sqrt(53.0);
BOOST_CHECK( fabs( x - Pt::dist(X,Y) ) < __DBL_EPSILON__ );
}

/*---------------------------------------*/
/* first lvl tests : nested calculus,... */
/*---------------------------------------*/


BOOST_AUTO_TEST_CASE(pt3D_product_left_and_right)
{
Pt::pt3D X(1,3,5);
Pt::pt3D Y(-2,4,1);

Pt::pt3D r1 = 3*X;
double x1 = (r1 - X*3).norm();
Pt::pt3D r2 = 0.5*Y;

double x2 = 2.0;
BOOST_CHECK( (x1 == 0.0) && (x2 == r2.maxLength())  );
}

/*---------------------------------------*/
/* second lvl tests : pure mathematics   */
/*---------------------------------------*/


BOOST_AUTO_TEST_CASE(pt3D_triple_product)
{
std::random_device rd;

std::mt19937 gen(rd());// random number generator: standard Mersenne twister initialized with seed rd()
std::uniform_real_distribution<> distrib(-1.0,1.0);

const Pt::pt3D X(distrib(gen),distrib(gen),distrib(gen));
const Pt::pt3D Y(distrib(gen),distrib(gen),distrib(gen));
const Pt::pt3D Z(distrib(gen),distrib(gen),distrib(gen));

double x1 = Pt::pTriple(X,Y,Z);
double x2 = Pt::pTriple(Y,Z,X);
double x3 = Pt::pTriple(Z,X,Y);

std::cout << "this test check that the triple product of three 3D random vectors is invariant by circular permutation, whatever is the numerical error" << std::endl;
std::cout << "X=" << X << std::endl;
std::cout << "Y=" << Y << std::endl;
std::cout << "Z=" << Z << std::endl;

std::cout << "pTriple(X,Y,Z)" << x1 << std::endl;
std::cout << "pTriple(Y,Z,X)" << x2 << std::endl;
std::cout << "pTriple(Z,X,Y)" << x3 << std::endl;

std::cout << "x1-x2=" << x1-x2 << std::endl;
std::cout << "x2-x3=" << x2-x3 << std::endl;
std::cout << "x3-x1=" << x3-x1 << std::endl;

std::cout << "__DBL_EPSILON__ =" << __DBL_EPSILON__ << std::endl;

BOOST_CHECK( (fabs(x1-x2)  < __DBL_EPSILON__) && (fabs(x2-x3)  < __DBL_EPSILON__) && ( fabs(x3-x1)  < __DBL_EPSILON__)  );
}

BOOST_AUTO_TEST_CASE(pt3D_det)
{
double M[Pt::DIM][Pt::DIM];
std::random_device rd;

std::mt19937 gen(rd());// random number generator: standard Mersenne twister initialized with seed rd()
std::uniform_real_distribution<> distrib(-1.0,1.0);

Pt::pt3D X(distrib(gen),distrib(gen),distrib(gen));

std::cout << "###################" << std::endl;
std::cout << "this test check that a random general rotation (filled using a random unit quaternions) respect sq(det-1) < __DBL_EPSILON__ " << std::endl;
std::cout << "X= " << X << std::endl;

while(Pt::norme2(X) > 1.0 )
    { X.x(distrib(gen)); X.y(distrib(gen)); X.z(distrib(gen)); }

double x = X.x();double y = X.y();double z = X.z();
double w = sqrt(1.0 - Pt::norme2(X) ); // w could be negative too

M[0][0] = 2.0*(x*x + w*w) - 1.0;
M[0][1] = 2.0*(x*y - z*w);
M[0][2] = 2.0*(x*z + y*w);
M[1][0] = 2.0*(x*y + z*w);
M[1][1] = 2.0*(y*y + w*w) - 1.0;
M[1][2] = 2.0*(y*z - x*w);
M[2][0] = 2.0*(x*z - y*w);
M[2][1] = 2.0*(y*z + x*w);
M[2][2] = 2.0*(z*z + w*w) - 1.0;

double result = Pt::det(M);

std::cout << "det(random_rot) -1=" << result-1.0 << std::endl;
BOOST_CHECK( (Pt::sq(result - 1.0)  < __DBL_EPSILON__ )  );
}

BOOST_AUTO_TEST_CASE(pt3D_inverse)
{
double M[Pt::DIM][Pt::DIM];
double inv_M[Pt::DIM][Pt::DIM];
std::random_device rd;

std::mt19937 gen(rd());// random number generator: standard Mersenne twister initialized with seed rd()
std::uniform_real_distribution<> distrib(-1.0,1.0);

Pt::pt3D X(distrib(gen),distrib(gen),distrib(gen));

std::cout << "###################" << std::endl;
std::cout << "this test check that a random general rotation multiplied by its inverse is ID : norm^2( M * M^-1 - M^-1 * M ) < __DBL_EPSILON__ " << std::endl;
std::cout << "* is performed using a tiny.h template , norm is the sqrt of the sum of the square of all components (Frobenius) " << std::endl;
std::cout << "X= " << X << std::endl;

while(Pt::norme2(X) > 1.0 )
    { X.x(distrib(gen)); X.y(distrib(gen)); X.z(distrib(gen)); }

double x = X.x();double y = X.y();double z = X.z();
double w = sqrt(1.0 - Pt::norme2(X) ); // w could be negative too

inv_M[0][0] = M[0][0] = 2.0*(x*x + w*w) - 1.0;
inv_M[0][1] = M[0][1] = 2.0*(x*y - z*w);
inv_M[0][2] = M[0][2] = 2.0*(x*z + y*w);
inv_M[1][0] = M[1][0] = 2.0*(x*y + z*w);
inv_M[1][1] = M[1][1] = 2.0*(y*y + w*w) - 1.0;
inv_M[1][2] = M[1][2] = 2.0*(y*z - x*w);
inv_M[2][0] = M[2][0] = 2.0*(x*z - y*w);
inv_M[2][1] = M[2][1] = 2.0*(y*z + x*w);
inv_M[2][2] = M[2][2] = 2.0*(z*z + w*w) - 1.0;

double det_M = Pt::det(M);

Pt::inverse(inv_M,det_M);

double M_inv_M[Pt::DIM][Pt::DIM];
double inv_M_M[Pt::DIM][Pt::DIM];

tiny::mult<double, Pt::DIM, Pt::DIM, Pt::DIM>(M,inv_M, M_inv_M);
tiny::mult<double, Pt::DIM, Pt::DIM, Pt::DIM>(inv_M, M, inv_M_M);

double result = 0.0;
for(int i=0;i<Pt::DIM;i++)
    for(int j=0;j<Pt::DIM;j++)
        result += Pt::sq(M_inv_M[i][j] - inv_M_M[i][j]);

std::cout << "norm^2(M M^-1 - M^-1 M) =" << result << std::endl;
BOOST_CHECK( result  < __DBL_EPSILON__ );
}

BOOST_AUTO_TEST_SUITE_END()
