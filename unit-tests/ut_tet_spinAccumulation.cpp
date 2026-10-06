#define BOOST_TEST_MODULE tetraSpinAccumulationTest

#include <boost/test/unit_test.hpp>

#include <random>

#include "tiny.h"
#include "ut_tools.h"
#include "ut_config.h"

BOOST_AUTO_TEST_SUITE(ut_tetra_spinAccumulation)

/* test of calc_Hst, return value Hst is a field, contribution of spin diffusion to LLG */
BOOST_AUTO_TEST_CASE(Tet_calc_Hst, *boost::unit_test::tolerance(1e-13))
    // 1e-13: the reference computes the prefactor with operations in a different order
    {
    using namespace Nodes;
    const int nbNod = 4;
    std::vector<double> V(nbNod); // scalar potential
    std::vector<Nodes::Node> node;
    dummyNodes<nbNod>(node);

    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);
    std::cout << "test on Hst with a tetrahedron\n V = ";
    for (int i = 0; i < nbNod; i++)
        {
        node[i].p[0] += distrib(gen)/10.0;
        node[i].p[1] += distrib(gen)/10.0;
        node[i].p[2] += distrib(gen)/10.0;
        V[i] = distrib(gen);
        std::cout << V[i] << '\t';
        }
    std::cout << std::endl;
    // carefull with indices (starting from 1)
    Tetra::Tet tet(node, 0, {1, 2, 3, 4});
    double Ms = distrib(gen);
    double sigma = distrib(gen);
    double N0 = distrib(gen);
    double l_sd = distrib(gen);
    double D0 = 2.0*sigma/(sq(CHARGE_ELECTRON)*N0);

    std::vector<Eigen::Vector3d> s(nbNod); // s is faking spin diffusion solution on the nodes
    for(int i=0;i<nbNod;i++)
        { s[i] = Eigen::Vector3d(distrib(gen),distrib(gen),distrib(gen)); }

    // code to test
    double prefactor = D0/(sq(l_sd)*gamma0*Ms);
    Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> Hst_to_test = Tetra::calc_Hst(tet,prefactor,s);
    // end code to test

    // ref version from June 2025 research version
    double Hst[DIM][Tetra::NPI] = {{0}};
    double nu0 = 1.0/mu0;
    double Js = mu0*Ms;
    const double Cst=D0/(sq(l_sd))/(gamma0*nu0*Js);
    double m_nod[3][Tetra::N];
    double m[3][Tetra::NPI];
    for(size_t ie=0;ie<Tetra::N;ie++)
        {
        size_t i = tet.ind[ie];
        for(size_t d=0;d<DIM;d++)
            { m_nod[d][ie]= s[i][d]; }
        }
    tiny::mult<double, DIM, Tetra::N, Tetra::NPI> (m_nod, Tetra::a, m);
    for(size_t npi=0;npi<Tetra::NPI;npi++)
        {
        Hst[IDX_X][npi] = Cst*m[IDX_X][npi];
        Hst[IDX_Y][npi] = Cst*m[IDX_Y][npi];
        Hst[IDX_Z][npi] = Cst*m[IDX_Z][npi];
        }
    // end ref version

    for(int i=0;i<DIM;i++)
        for(int j=0;j<Tetra::NPI;j++)
            {
            std::cout << Hst_to_test(i,j) << " should be " << Hst[i][j] << std::endl;
            BOOST_TEST( Hst_to_test(i,j) == Hst[i][j] );
            }
    }

BOOST_AUTO_TEST_CASE(Tet_Hst_order2_contribution,*boost::unit_test::tolerance(1e-6))
    {
    /* the second order term of the spin accumulation field must be the derivative of its tangent
     part Hst - (u.Hst) u in the direction v (v tangent to u), projected on the tangent plane:
     P_u (f(u + eps v) - f(u - eps v))/(2 eps) = -(u.Hst) v */
    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);
    if (!DET_UT) std::cout << "seed =" << sd << std::endl;

    std::vector<Nodes::Node> node;
    const int nbNod = 4;
    dummyNodes<nbNod>(node);
    Tetra::Tet t(node, 0, {1, 2, 3, 4});

    Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> Hst, U, V, H;
    for (int npi = 0; npi < Tetra::NPI; npi++)
        {
        U.col(npi) = rand_vec3d(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
        Eigen::Vector3d w = rand_vec3d(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
        V.col(npi) = w - w.dot(U.col(npi)) * U.col(npi);  // tangent to u
        Hst.col(npi) = (0.5 + distrib(gen)) * rand_vec3d(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
        }
    H = t.Hst_order2_contribution(Hst, U, V);

    const double eps = 1e-6;// eps is used to compute a directionnal derivative
    auto tangent = [](const Eigen::Vector3d &u, const Eigen::Vector3d &h)
        { return Eigen::Vector3d(h - u.dot(h) * u); };
    for (int npi = 0; npi < Tetra::NPI; npi++)
        {
        const Eigen::Vector3d u = U.col(npi), v = V.col(npi), h = Hst.col(npi);
        Eigen::Vector3d dH_fd = (tangent(u + eps*v, h) - tangent(u - eps*v, h))/(2*eps);
        dH_fd -= u.dot(dH_fd) * u;  // projection on the plane tangent to u
        double relative_err = (H.col(npi) - dH_fd).norm()/dH_fd.norm();
        std::cout << "npi= " << npi << " relative error on Hst'[v] = " << relative_err << std::endl;
        BOOST_TEST( relative_err == 0.0, "spin accumulation: second order term is not the derivative of Hst" );
        }
    }

BOOST_AUTO_TEST_SUITE_END()
