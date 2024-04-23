#define BOOST_TEST_MODULE tet_lumpingTest

#include <boost/test/unit_test.hpp>
#include <random>

#include "ut_tools.h"
#include "ut_config.h"

BOOST_AUTO_TEST_SUITE(ut_tet_lumping)

/*-----------------------------------------------------

class tetra Lumping is tested here

-------------------------------------------------------*/

BOOST_AUTO_TEST_CASE(tet_exchange_lumping, *boost::unit_test::tolerance(UT_TOL))
    {
    const int nbNod = 4;
    std::vector<Nodes::Node> node;
    dummyNodes<nbNod>(node);

    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);

    for (int i = 0; i < nbNod; i++)
        {
        node[i].d[0].u = rand_vec3d(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
        }

    // carefull with indices (starting from 1)
    Tetra::Tet t(node, 0, {1, 2, 3, 4});

    double prefactor = 1.0 + distrib(gen);
    Eigen::Matrix<double, 3*Tetra::N, 3*Tetra::N> AE;
    AE.setZero();
    // ref code
    Eigen::Matrix<double,Tetra::N,Tetra::NPI> dadx,dady,dadz;// matrices dad(x|y|z) are repetition of da parts
    for (int j = 0; j < Tetra::NPI; j++)
        {
        for (int i = 0; i < Tetra::N; i++)
            {
            dadx(i,j) = t.da(i,0);
            dady(i,j) = t.da(i,1);
            dadz(i,j) = t.da(i,2);
            }
        }
    
    for(int npi = 0; npi < Tetra::NPI; npi++)
        {
        const double w = t.weight[npi];

        for (int i = 0; i < Tetra::N; i++)
            {
            for (int j = 0; j < Tetra::N; j++)
                {
                double contrib = w * prefactor
                                 * (dadx(i,npi) * dadx(j,npi) + dady(i,npi) * dady(j,npi)
                                    + dadz(i,npi) * dadz(j,npi));

                AE(i,j) += contrib;
                AE(Tetra::N + i,Tetra::N + j) += contrib;
                AE(2*Tetra::N + i,2*Tetra::N + j) += contrib;
                }
            }
        }
    // end ref code

    Eigen::Matrix<double, 3*Tetra::N, 3*Tetra::N> AE_to_check;
    AE_to_check.setZero();
    t.exchange_lumping(prefactor,AE_to_check);
    for (int i = 0; i < Tetra::N; i++)
        for (int j = 0; j < Tetra::N; j++)
            {// AE is block diagonal, we just print the block
            std::cout << "ref val= " << AE_to_check(i,j) << "found = " << AE(i,j) << std::endl;
            }

    for (int i = 0; i < 3*Tetra::N; i++)
        for (int j = 0; j < 3*Tetra::N; j++)
            {
            BOOST_TEST(AE_to_check(i,j) == AE(i,j));
            }
    }

BOOST_AUTO_TEST_CASE(tet_lumping, *boost::unit_test::tolerance(UT_TOL))
    {
    Eigen::Matrix<double, 3*Tetra::N, 3*Tetra::N> AE_to_check;
    AE_to_check.setZero();
    double AE[3 * Tetra::N][3 * Tetra::N] = {{0}};

    const int nbNod = 4;
    std::vector<Nodes::Node> node;
    dummyNodes<nbNod>(node);

    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);
    
    for (int i = 0; i < nbNod; i++)
        {
        node[i].d[0].u = rand_vec3d(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
        }

    // carefull with indices (starting from 1)
    Tetra::Tet t(node, 0, {1, 2, 3, 4});
    
    double a = distrib(gen);
    double b = distrib(gen);

    timing prm_t = timing(1.0, std::min(a, b), std::max(a, b));

    double dt = prm_t.get_dt();
    double s_dt = THETA * dt;  // theta from theta scheme in config.h.in
    double alpha_LLG = 0.5;
    Eigen::Matrix<double,Tetra::NPI,1> uHeff;
    uHeff.setConstant(distrib(gen));
    Eigen::Matrix<double,Tetra::NPI,1> alfa = Tetra::calc_alpha_eff(dt, alpha_LLG, uHeff);
    double A = distrib(gen);         // Ae
    double Js = 0.5 + distrib(gen);  // 0.5 offset to center on 1 the Js value
    double Abis = 2.0 * A / Js;
    double TAUR = prm_t.TAUR;

    // ref code (with minimal adaptations of integrales method in file MuMag_integrales.cc of
    // src_Tube_scalfmm_thiaville_ec_mu_oersted_thiele_dyn20180903.tgz )
    double u_nod[3][Tetra::N];
    for (int ie = 0; ie < Tetra::N; ie++)
        {
        int i = t.ind[ie];
        Nodes::dataNode &d0 = node[i].d[0];
        for (int dim = 0; dim < Nodes::DIM; dim++)
            {
            u_nod[dim][ie] = d0.u(dim);  // aimantation
            }
        }

    for (int npi = 0; npi < Tetra::NPI; npi++)
        {
        double w, ai, dai_dx, dai_dy, dai_dz, daj_dx, daj_dy, daj_dz;
        double Dai_Daj;

        double R = dt / TAUR * abs(log(dt / TAUR));

        w = t.weight[npi];
        for (int ie = 0; ie < Tetra::N; ie++)
            {
            ai = Tetra::a[ie][npi];
            dai_dx = t.da(ie,0);
            dai_dy = t.da(ie,1);
            dai_dz = t.da(ie,2);

            AE[ie][ie] += alfa(npi) * ai * w;  // lumping
            AE[Tetra::N + ie][Tetra::N + ie] += alfa(npi) * ai * w;
            AE[2 * Tetra::N + ie][2 * Tetra::N + ie] += alfa(npi) * ai * w;

            AE[ie][2 * Tetra::N + ie] += +u_nod[1][ie] * ai * w;  // lumping
            AE[ie][Tetra::N + ie] += -u_nod[2][ie] * ai * w;
            AE[Tetra::N + ie][ie] += +u_nod[2][ie] * ai * w;
            AE[Tetra::N + ie][2 * Tetra::N + ie] += -u_nod[0][ie] * ai * w;
            AE[2 * Tetra::N + ie][Tetra::N + ie] += +u_nod[0][ie] * ai * w;
            AE[2 * Tetra::N + ie][ie] += -u_nod[1][ie] * ai * w;

            for (int je = 0; je < Tetra::N; je++)
                {
                daj_dx = t.da(je,0);
                daj_dy = t.da(je,1);
                daj_dz = t.da(je,2);
                Dai_Daj = dai_dx * daj_dx + dai_dy * daj_dy + dai_dz * daj_dz;

                AE[ie][je] += s_dt * (1. + R) * 2 * A / Js * Dai_Daj * w;
                AE[Tetra::N + ie][Tetra::N + je] += s_dt * (1. + R) * 2 * A / Js * Dai_Daj * w;
                AE[2 * Tetra::N + ie][2 * Tetra::N + je] +=
                        s_dt * (1. + R) * 2 * A / Js * Dai_Daj * w;
                }
            }
        }
    // end ref code

    // code to check
    t.lumping(alfa, prm_t.prefactor * s_dt * Abis, AE_to_check);
    // end code to check

    for (int i = 0; i < 3*Tetra::N; i++)
        for (int j = 0; j < 3*Tetra::N; j++)
            {
            BOOST_TEST(AE_to_check(i,j) == AE[i][j]);
            }
    if (!DET_UT) std::cout << "seed =" << sd << std::endl;
    }

BOOST_AUTO_TEST_SUITE_END()
