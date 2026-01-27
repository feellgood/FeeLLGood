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

    // begin code to test
    using namespace Tetra;
    Eigen::Matrix<double, 3*N, 3*N> AE_to_check;
    AE_to_check.setZero();
    Eigen::Matrix<double,N,N> exch_block = t.da*t.da.transpose();
    exch_block *= prefactor*t.weight.sum();
    AE_to_check.block<N,N>(0,0) += exch_block;
    AE_to_check.block<N,N>(N,N) += exch_block;
    AE_to_check.block<N,N>(2*N,2*N) += exch_block;
    // end code to test
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            {// AE is block diagonal, we just print the block
            std::cout << "ref val= " << AE_to_check(i,j) << "found = " << AE(i,j) << std::endl;
            }

    for (int i = 0; i < 3*N; i++)
        for (int j = 0; j < 3*N; j++)
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

BOOST_AUTO_TEST_CASE(tet_spin_diff_AE_filling, *boost::unit_test::tolerance(2.0*UT_TOL))
    {
    std::cout <<"\ttest on spin diffusion AE filling\n";
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

    const double lsf = distrib(gen);
    const double lsd = distrib(gen);
    const double D0 = distrib(gen);
    /* units: [D0/sq(lsf)] = s^-1 : it is 1/tau_sf */
    using algebra::sq;
    using namespace Tetra;
    using namespace Nodes;
    Eigen::Matrix<double,N,1> a_w = eigen_a*t.weight;
    Eigen::Matrix<double,N,1> diag = (D0/sq(lsf))*a_w;
    Eigen::Matrix<double,N,N> diagBlock = t.calcDiagBlock(D0,diag);
    AE_to_check.block<N,N>(0,0) += diagBlock;
    AE_to_check.block<N,N>(N,N) += diagBlock;
    AE_to_check.block<N,N>(2*N,2*N) += diagBlock;

    const double invTau_sd = D0/sq(lsd); //units: [D0/sq(lsd)] = s^-1 : it is 1/tau_sd
    diag = invTau_sd * a_w.cwiseProduct(t.calcOffDiagBlock(IDX_X));
    AE_to_check.block<N,N>(N,2*N).diagonal() += diag;
    AE_to_check.block<N,N>(2*N,N).diagonal() -= diag;
    diag = invTau_sd * a_w.cwiseProduct(t.calcOffDiagBlock(IDX_Y));
    AE_to_check.block<N,N>(0,2*N).diagonal() -= diag;
    AE_to_check.block<N,N>(2*N,0).diagonal() += diag;
    diag = invTau_sd * a_w.cwiseProduct(t.calcOffDiagBlock(IDX_Z));
    AE_to_check.block<N,N>(0,N).diagonal() += diag;
    AE_to_check.block<N,N>(N,0).diagonal() -= diag;

// start code ref (from june 2025)
double u_nod[3][Tetra::N];
    for (int ie = 0; ie < Tetra::N; ie++)
        {
        int i = t.ind[ie];
        Nodes::dataNode &d0 = node[i].d[0];
        for (int dim = 0; dim < Nodes::DIM; dim++)
            {
            u_nod[dim][ie] = d0.u(dim);  // magnetization
            }
        }
const double ilJ = 1.0/lsd;
for (size_t npi=0; npi<NPI; npi++)
    {
    double w, ai, dai_dx, dai_dy, dai_dz, aj, daj_dx, daj_dy, daj_dz, Dai_Daj;
    w = t.weight[npi];

    for (size_t ie=0; ie<N; ie++)
        {
        ai = Tetra::a[ie][npi];
        dai_dx = t.da(ie,0);
        dai_dy = t.da(ie,1);
        dai_dz = t.da(ie,2);
        AE[    ie][     ie] +=  D0/pow(lsf, 2.)* ai *w;
        AE[  N+ie][   N+ie] +=  D0/pow(lsf, 2.)* ai *w;
        AE[2*N+ie][ 2*N+ie] +=  D0/pow(lsf, 2.)* ai *w;

        AE[    ie][   N+ie] +=  D0*pow(ilJ, 2.)*u_nod[2][ie] * ai *w;
        AE[    ie][ 2*N+ie] += -D0*pow(ilJ, 2.)*u_nod[1][ie] * ai *w;

        AE[  N+ie][     ie] += -D0*pow(ilJ, 2.)*u_nod[2][ie] * ai *w;
        AE[  N+ie][ 2*N+ie] +=  D0*pow(ilJ, 2.)*u_nod[0][ie] * ai *w;

        AE[2*N+ie][     ie] +=  D0*pow(ilJ, 2.)*u_nod[1][ie] * ai *w;
        AE[2*N+ie][   N+ie] += -D0*pow(ilJ, 2.)*u_nod[0][ie] * ai *w;

        for (size_t je=0; je<N; je++)
            {
            aj =  Tetra::a[je][npi];
            daj_dx = t.da(je,0);
            daj_dy = t.da(je,1);
            daj_dz = t.da(je,2);
            Dai_Daj = dai_dx*daj_dx + dai_dy*daj_dy + dai_dz*daj_dz;
            AE[    ie][     je] += D0* Dai_Daj *w;
            AE[  N+ie][   N+je] += D0* Dai_Daj *w;
            AE[2*N+ie][ 2*N+je] += D0* Dai_Daj *w;
	        }
	    }
    }
// end code ref

    for (int i = 0; i < 3*Tetra::N; i++)
        for (int j = 0; j < 3*Tetra::N; j++)
            {
            //std::cout << "AE_to_check[" << i << "," << j << "]=" << AE_to_check(i,j) << "to compare to " << AE[i][j] << std::endl;
            BOOST_TEST(AE_to_check(i,j) == AE[i][j]);
            }
    }

BOOST_AUTO_TEST_CASE(tet_spin_diff_BE_filling, *boost::unit_test::tolerance(2.0*UT_TOL))
    {
    std::cout <<"\ttest on spin diffusion BE filling\n";
    const int nbNod = 4;
    std::vector<Nodes::Node> node;
    std::vector<double> V(nbNod);
    dummyNodes<nbNod>(node);

    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);

    for (int i = 0; i < nbNod; i++)
        {
        node[i].d[0].u = rand_vec3d(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
        V[i] = distrib(gen);
        }

    // carefull with indices (starting from 1)
    Tetra::Tet t(node, 0, {1, 2, 3, 4});
    double sigma = distrib(gen);
    double P = distrib(gen);

    using namespace Tetra;
    double beta = P;
    double C0 = sigma;
    std::vector<double> BE(3*N,0.0);
    double u_nod[3][N];
    Eigen::Matrix<double,N,1> V_nod;

    Eigen::Matrix<double,N,NPI> dadx,dady,dadz;// matrices dad(x|y|z) are repetition of da parts
    for (int j = 0; j < NPI; j++)
        {
        for (int i = 0; i < N; i++)
            {
            dadx(i,j) = t.da(i,0);
            dady(i,j) = t.da(i,1);
            dadz(i,j) = t.da(i,2);
            }
        }

        for (size_t ie=0; ie<N; ie++)
            {
            V_nod[ie] = V[t.ind[ie]];
            for(size_t d=0;d<3;d++) { u_nod[d][ie] = node[t.ind[ie]].d[0].u[d]; }
            }

    Eigen::Matrix<double,NPI,1> dVdx = V_nod.transpose() * dadx;//tiny::transposed_mult<double, N, NPI> (V_nod, dadx, dVdx);
    Eigen::Matrix<double,NPI,1> dVdy = V_nod.transpose() * dady;//tiny::transposed_mult<double, N, NPI> (V_nod, dady, dVdy);
    Eigen::Matrix<double,NPI,1> dVdz = V_nod.transpose() * dadz;//tiny::transposed_mult<double, N, NPI> (V_nod, dadz, dVdz);

    // start code to test
    std::vector<double> BE_to_check(3*N,0.0);
    const double cst0 = BOHRS_MUB*P*sigma/CHARGE_ELECTRON;
    // the three following messy lines mimic pre_compute() method
    std::vector< Eigen::Matrix<double,Nodes::DIM,NPI> > gradV;
    gradV.push_back(calc_gradV(t,V));
    Eigen::Matrix<double,Nodes::DIM,NPI> &_gradV = gradV[t.idx];

    for (size_t npi=0; npi<NPI; npi++)
        {
        const double cst0_w = cst0*t.weight[npi];

        for (size_t ie=0; ie<N; ie++)
            {
            const Eigen::Vector3d &m = node[t.ind[ie]].d[0].u;//magnetization
            Eigen::Vector3d grad_ai = t.da.row(ie);
            double tmp = cst0_w*grad_ai.dot( _gradV.col(npi) );
            BE_to_check[    ie] += tmp*m[0];
            BE_to_check[  N+ie] += tmp*m[1];
            BE_to_check[2*N+ie] += tmp*m[2];
            }
        }
    // end code to test

    // start code ref (from june 2025)
    for (size_t npi=0; npi<NPI; npi++)
        {
        double w = t.weight[npi];
        for (size_t ie=0; ie<N; ie++)
            {
            double dai_dx = t.da(ie,0);
            double dai_dy = t.da(ie,1);
            double dai_dz = t.da(ie,2);
            /* Changement de convention de signe pour le courant dans l'expression du ST tel que j = -C0 grad V */
            double Dai_DV = dai_dx * dVdx[npi] + dai_dy * dVdy[npi] + dai_dz * dVdz[npi];
            BE[    ie] += BOHRS_MUB*beta*C0/CHARGE_ELECTRON*u_nod[0][ie]* Dai_DV *w; // lumping
            BE[  N+ie] += BOHRS_MUB*beta*C0/CHARGE_ELECTRON*u_nod[1][ie]* Dai_DV *w; // lumping
            BE[2*N+ie] += BOHRS_MUB*beta*C0/CHARGE_ELECTRON*u_nod[2][ie]* Dai_DV *w; // lumping
	        }
        }
    // end code ref
    for(size_t i=0;i<BE.size();i++)
        {
        std::cout << BE_to_check[i] << " should be " << BE[i] << std::endl;
        BOOST_TEST( BE_to_check[i] == BE[i]);
        }
    }

BOOST_AUTO_TEST_SUITE_END()
