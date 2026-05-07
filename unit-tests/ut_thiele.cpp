#define BOOST_TEST_MODULE thieleTest

#include <boost/test/unit_test.hpp>

#include <random>

#include "tiny.h"
#include "ut_tools.h"
#include "ut_config.h"

BOOST_AUTO_TEST_SUITE(ut_thiele)

BOOST_AUTO_TEST_CASE(thiele_on_single_tetra, *boost::unit_test::tolerance(UT_TOL))
    {
    using namespace Nodes; 
    const int nbNod = 4;
    std::vector<Nodes::Node> node;
    dummyNodes<nbNod>(node);

    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);

    for (int i = 0; i < nbNod; i++)
        { node[i].d[0].u = rand_vec3d(M_PI * distrib(gen), 2 * M_PI * distrib(gen)); }

    Tetra::Tet t(node, 0, {1, 2, 3, 4});// carefull with indices (starting from 1)

    // code to test: lambda in transform_reduce of thiele method in mesh class
    Eigen::Matrix<double,Nodes::DIM,Tetra::N> mag_nod;
    for (int i = 0; i< Tetra::N; i++)
        { mag_nod.col(i) = node[i].d[0].u; }//Mesh::getNode_u(te.ind[i]);

    Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> du_dz = mag_nod * (t.da.col(Nodes::IDX_Z)).replicate(1,Tetra::NPI);
    double valToTest(0);
    for (int npi=0;npi<Tetra::NPI;npi++)
        valToTest += du_dz.col(npi).dot(du_dz.col(npi)) * t.weight[npi];
    // end code to test

    // code ref
    double u_nod[3][Tetra::N];
    double dudz[3][Tetra::NPI];
    double dadz[Tetra::N][Tetra::NPI];

    for (int ie=0; ie<Tetra::N; ie++)
        {
        for (int d=0; d<3; d++)
            {u_nod[d][ie] = node[ie].d[0].u[d];}
        }

    for (int i=0; i<Tetra::N; i++)
        for (int j=0; j<Tetra::NPI; j++)
            { dadz[i][j] = t.da(i,Nodes::IDX_Z); } // yes this is silly ...

    tiny::mult<double, 3, Tetra::N, Tetra::NPI> (u_nod, dadz, dudz);
    double sum_gradu_sq = 0.;
    for (int npi=0; npi<Tetra::NPI; npi++)
        {
		sum_gradu_sq += (dudz[0][npi]*dudz[0][npi]+dudz[1][npi]*dudz[1][npi]+dudz[2][npi]*dudz[2][npi])* t.weight[npi];
		}
    // end code ref
    BOOST_TEST( sum_gradu_sq == valToTest );
    }

BOOST_AUTO_TEST_SUITE_END()

