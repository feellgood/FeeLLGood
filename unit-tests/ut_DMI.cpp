#define BOOST_TEST_MODULE DMI_Test

#include <boost/test/unit_test.hpp>
#include <random>

#include "ut_tools.h"
#include "ut_config.h"

BOOST_AUTO_TEST_SUITE(ut_DMI)

BOOST_AUTO_TEST_CASE(tet_DMI_bulk_BE, *boost::unit_test::tolerance(10.0*UT_TOL))
	{
	// in that test we check the structure of the contributions of DMI_bulk to RHS BE
	// and the values of the coeff, not the prefactor depending on material and other params
	using namespace Nodes;
    using namespace Tetra;

    const int nbNod = 4;
	std::vector<Nodes::Node> node;
	dummyNodes<nbNod>(node);	
	unsigned sd = my_seed();
	std::mt19937 gen(sd);
	std::uniform_real_distribution<> distrib(0.0, 1.0);
	// carefull with indices (starting from 1)
	Tetra::Tet t(node, 0, {1, 2, 3, 4});
	for (int i = 0; i < nbNod; i++)
	    { node[i].d[Nodes::CURRENT].u = rand_vec3d(M_PI * distrib(gen), 2 * M_PI * distrib(gen)); }	

    // begin code to test
    Eigen::Matrix<double, DIM, N> BE_to_check;
    BE_to_check.setZero();
    double cst = distrib(gen); // cst should be = DMI_bulk/(mu*Ms), replaced by random number
    Eigen::Matrix<double,DIM,NPI> U;
    Eigen::Matrix<double,DIM,NPI> dUdx;
    Eigen::Matrix<double,DIM,NPI> dUdy;
    Eigen::Matrix<double,DIM,NPI> dUdz;
    t.interpolation(Nodes::get_u<CURRENT>, U, dUdx, dUdy, dUdz);
    t.DMI_bulk_BE(cst,U,dUdx,dUdy,dUdz,BE_to_check);
    // end code to test

    //begin code ref
    double BE[DIM*N] = {0};

    for (int npi=0;npi<NPI;npi++)
        {
        const double w = t.weight[npi];
        for(int ie = 0;ie<N; ie++)
            {
            const double ai = a[ie][npi];

            // Dbulk/Js \omega \cdot (\nabla \times u)
            BE[    ie]  += cst* ai*(dUdy(2,npi)-dUdz(1,npi)) *w;
            BE[  N+ie]  += cst* ai*(dUdz(0,npi)-dUdx(2,npi)) *w;
            BE[2*N+ie]  += cst* ai*(dUdx(1,npi)-dUdy(0,npi)) *w;

            // Dbulk/Js m \cdot (\nabla \times \omega)
            const double dai_dx = t.da(ie,IDX_X);
            const double dai_dy = t.da(ie,IDX_Y);
            const double dai_dz = t.da(ie,IDX_Z);

            BE[    ie]  += cst* (dai_dz * U(1,npi) - dai_dy * U(2,npi)) *w;
            BE[  N+ie]  += cst* (dai_dx * U(2,npi) - dai_dz * U(0,npi)) *w;
            BE[2*N+ie]  += cst* (dai_dy * U(0,npi) - dai_dx * U(1,npi)) *w;
            }
        }
    //end code ref
			
    // this might be better to reshape BE than direct comparing BE's
    for (int i = 0; i < DIM; i++)
        for (int j = 0; j < N; j++)
            {
            std::cout << "BE=" << BE_to_check(i,j) << "\t; ref= " << BE[N*i+j] << std::endl;
            BOOST_TEST(BE_to_check(i,j) == BE[N*i+j]);
            }
    }

BOOST_AUTO_TEST_SUITE_END()
