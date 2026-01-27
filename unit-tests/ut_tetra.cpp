#define BOOST_TEST_MODULE tetraTest

#include <boost/test/unit_test.hpp>

#include <random>

#include "tiny.h"
#include "ut_tools.h"
#include "ut_config.h"

BOOST_AUTO_TEST_SUITE(ut_tetra)

/*-----------------------------------------------------

 Tetra::Tet constructor is tested in ut_element.cpp
Lumping is tested in ut_tet_lumping.cpp

Here are tested constant table values, interpolations, various inner matrices and P matrix

---------------------------------------*/
BOOST_AUTO_TEST_CASE(Tet_inner_tables, *boost::unit_test::tolerance(UT_TOL))
    {
    using namespace Nodes;
    // this test is dedicated to  check dadx,dady,dadz and weight tables, those values are
    // initialized once by Tet constructor
    std::cout << "constructor test with 4 nodes in node vector\n";
    const int nbNod = 4;
    std::vector<Nodes::Node> node;
    dummyNodes<nbNod>(node);

    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);

    for (int i = 0; i < nbNod; i++)
        { node[i].d[0].u = rand_vec3d(M_PI * distrib(gen), 2 * M_PI * distrib(gen)); }

    Tetra::Tet t(node, 0, {1, 2, 3, 4});// carefull with indices (starting from 1)
    t.infos();

    // ref code (with minimal adaptations of dad(x|y|z) in file Mesh_hat.cc of
    // src_Tube_scalfmm_thiaville_ec_mu_oersted_thiele_dyn20180903.tgz )
    double _dadx[Tetra::N][Tetra::NPI];
    double _dady[Tetra::N][Tetra::NPI];
    double _dadz[Tetra::N][Tetra::NPI];
    double da[Tetra::N][DIM];
    double J[DIM][DIM];
    double nod[DIM][Tetra::N];
    double weight[Tetra::NPI];
    constexpr double dadu[Tetra::N][DIM] = {{-1., -1., -1.}, {1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

    for (int ie = 0; ie < Tetra::N; ie++)
        {
        int i = t.ind[ie];
        nod[0][ie] = node[i].p(0);
        nod[1][ie] = node[i].p(1);
        nod[2][ie] = node[i].p(2);
        }

    tiny::mult<double, DIM, Tetra::N, DIM>(nod, dadu, J);
    double detJ = det(J);// from ut_tools.h
    inverse(J, detJ);// from ut_tools.h
    tiny::mult<double, Tetra::N, DIM, DIM>(dadu, J, da);

    for (int npi = 0; npi < Tetra::NPI; npi++)
        {
        for (int ie = 0; ie < Tetra::N; ie++)
            {
            _dadx[ie][npi] = da[ie][0];
            _dady[ie][npi] = da[ie][1];
            _dadz[ie][npi] = da[ie][2];
            }
        weight[npi] = detJ * Tetra::pds[npi];
        }

    // end ref code

    // the Tet constructor is computing dad(x|y|z)

    double result_dadx(0.0), result_dady(0.0), result_dadz(0.0), result_w(0.0);
    for (int npi = 0; npi < Tetra::NPI; npi++)
        {
        for (int ie = 0; ie < Tetra::N; ie++)
            {
            result_dadx += sq(_dadx[ie][npi] - t.da(ie,0));
            result_dady += sq(_dady[ie][npi] - t.da(ie,1));
            result_dadz += sq(_dadz[ie][npi] - t.da(ie,2));
            }
        result_w += sq(weight[npi] - t.weight[npi]);
        }

    std::cout << "sq_frob norm diff dadx,dady,dadz= " << result_dadx << " ; " << result_dady
              << " ; " << result_dadz << std::endl;
    if (!DET_UT) std::cout << "seed =" << sd << std::endl;
    BOOST_TEST(sqrt(result_dadx) == 0.0);
    BOOST_TEST(sqrt(result_dady) == 0.0);
    BOOST_TEST(sqrt(result_dadz) == 0.0);

    std::cout << "frob norm diff weight= " << result_w << std::endl;
    BOOST_TEST(sqrt(result_w) == 0.0);
    }

BOOST_AUTO_TEST_CASE(Tet_calc_vol, *boost::unit_test::tolerance(UT_TOL))
    {
    const int nbNod = 4;
    std::vector<Nodes::Node> node;
    dummyNodes<nbNod>(node);

    // carefull with indices (starting from 1)
    Tetra::Tet t(node, 0, {1, 2, 3, 4});

    double result = 1 / 6.0;
    double vol = t.calc_vol();
    std::cout << "vol(tetra) =" << vol << std::endl;
    BOOST_TEST(vol == result);
    }

BOOST_AUTO_TEST_CASE(Tet_nod_interpolation, *boost::unit_test::tolerance(UT_TOL))
    {
    using namespace Nodes;
    const int nbNod = 4;
    std::vector<Nodes::Node> node;
    dummyNodes<nbNod>(node);
    
    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);

    for (int i = 0; i < nbNod; i++)
        {
        node[i].d[0].u = rand_vec3d(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
        node[i].d[0].v = rand_vec3d(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
        }

    // carefull with indices (starting from 1)
    Tetra::Tet t(node, 0, {1, 2, 3, 4});

    double t_dadx[Tetra::N][Tetra::NPI],t_dady[Tetra::N][Tetra::NPI],t_dadz[Tetra::N][Tetra::NPI];
    
    for(int i=0;i<Tetra::N;i++)
        for(int j=0;j<Tetra::NPI;j++)
            {
            t_dadx[i][j] = t.da(i,0);
            t_dady[i][j] = t.da(i,1);
            t_dadz[i][j] = t.da(i,2);
            }
    // ref code (with minimal adaptations of integrales method in file MuMag_integrales.cc of
    // src_Tube_scalfmm_thiaville_ec_mu_oersted_thiele_dyn20180903.tgz )
    double _u_nod[DIM][Tetra::N], _u[DIM][Tetra::NPI];
    double dudx[DIM][Tetra::NPI], dudy[DIM][Tetra::NPI], dudz[DIM][Tetra::NPI];

    double _v_nod[DIM][Tetra::N], _v[DIM][Tetra::NPI];
    double dvdx[DIM][Tetra::NPI], dvdy[DIM][Tetra::NPI], dvdz[DIM][Tetra::NPI];

    for (int ie = 0; ie < Tetra::N; ie++)
        {
        int i = t.ind[ie];
        Nodes::dataNode &d0 = node[i].d[0];
        for (int dim = 0; dim < DIM; dim++)
            {
            _u_nod[dim][ie] = d0.u(dim);
            _v_nod[dim][ie] = d0.v(dim);
            }
        }
    tiny::mult<double, DIM, Tetra::N, Tetra::NPI>(_u_nod, Tetra::a, _u);
    tiny::mult<double, DIM, Tetra::N, Tetra::NPI>(_u_nod, t_dadx, dudx);
    tiny::mult<double, DIM, Tetra::N, Tetra::NPI>(_u_nod, t_dady, dudy);
    tiny::mult<double, DIM, Tetra::N, Tetra::NPI>(_u_nod, t_dadz, dudz);

    tiny::mult<double, DIM, Tetra::N, Tetra::NPI>(_v_nod, Tetra::a, _v);
    tiny::mult<double, DIM, Tetra::N, Tetra::NPI>(_v_nod, t_dadx, dvdx);
    tiny::mult<double, DIM, Tetra::N, Tetra::NPI>(_v_nod, t_dady, dvdy);
    tiny::mult<double, DIM, Tetra::N, Tetra::NPI>(_v_nod, t_dadz, dvdz);
    // end ref code

    // code to check
    Eigen::Matrix<double,DIM,Tetra::NPI> dUdx;
    Eigen::Matrix<double,DIM,Tetra::NPI> dUdy;
    Eigen::Matrix<double,DIM,Tetra::NPI> dUdz;
    Eigen::Matrix<double,DIM,Tetra::NPI> dVdx;
    Eigen::Matrix<double,DIM,Tetra::NPI> dVdy;
    Eigen::Matrix<double,DIM,Tetra::NPI> dVdz;
    Eigen::Matrix<double,DIM,Tetra::NPI> U;
    Eigen::Matrix<double,DIM,Tetra::NPI> V;

    t.interpolation(Nodes::get_u<Nodes::CURRENT>, U, dUdx, dUdy, dUdz);
    t.interpolation(Nodes::get_v<Nodes::CURRENT>, V, dVdx, dVdy, dVdz);
    // end code to check

    double n_u = tiny::frob_norm<double, DIM, Tetra::NPI>(_u);
    double n_dudx = tiny::frob_norm<double, DIM, Tetra::NPI>(dudx);
    double n_dudy = tiny::frob_norm<double, DIM, Tetra::NPI>(dudy);
    double n_dudz = tiny::frob_norm<double, DIM, Tetra::NPI>(dudz);

    double n_v = tiny::frob_norm<double, DIM, Tetra::NPI>(_v);
    double n_dvdx = tiny::frob_norm<double, DIM, Tetra::NPI>(dvdx);
    double n_dvdy = tiny::frob_norm<double, DIM, Tetra::NPI>(dvdy);
    double n_dvdz = tiny::frob_norm<double, DIM, Tetra::NPI>(dvdz);

    double dist_uU = sq_dist(_u, U);
    double dist_dudx_dUdx = sq_dist(dudx, dUdx);
    double dist_dudy_dUdy = sq_dist(dudy, dUdy);
    double dist_dudz_dUdz = sq_dist(dudz, dUdz);

    double dist_vV = sq_dist(_v, V);
    double dist_dvdx_dVdx = sq_dist(dvdx, dVdx);
    double dist_dvdy_dVdy = sq_dist(dvdy, dVdy);
    double dist_dvdz_dVdz = sq_dist(dvdz, dVdz);

    double n_U = sqrt(sq_frobenius_norm<Tetra::NPI>(U));
    double n_dUdx = sqrt(sq_frobenius_norm<Tetra::NPI>(dUdx));
    double n_dUdy = sqrt(sq_frobenius_norm<Tetra::NPI>(dUdy));
    double n_dUdz = sqrt(sq_frobenius_norm<Tetra::NPI>(dUdz));

    double n_V = sqrt(sq_frobenius_norm<Tetra::NPI>(V));
    double n_dVdx = sqrt(sq_frobenius_norm<Tetra::NPI>(dVdx));
    double n_dVdy = sqrt(sq_frobenius_norm<Tetra::NPI>(dVdy));
    double n_dVdz = sqrt(sq_frobenius_norm<Tetra::NPI>(dVdz));

    if (!DET_UT) std::cout << "seed =" << sd << std::endl;
    std::cout << "distance^2 (u,U) =" << dist_uU << std::endl;
    BOOST_TEST(sqrt(dist_uU) == 0.0);
    BOOST_TEST(sqrt(dist_dudx_dUdx) == 0.0);
    BOOST_TEST(sqrt(dist_dudy_dUdy) == 0.0);
    BOOST_TEST(sqrt(dist_dudz_dUdz) == 0.0);

    std::cout << "distance^2 (v,V) =" << dist_vV << std::endl;
    BOOST_TEST(sqrt(dist_vV) == 0.0);
    BOOST_TEST(sqrt(dist_dvdx_dVdx) == 0.0);
    BOOST_TEST(sqrt(dist_dvdy_dVdy) == 0.0);
    BOOST_TEST(sqrt(dist_dvdz_dVdz) == 0.0);

    // to avoid gag of comparing pure zeros we also check that matrices norm are equal
    // let's be paranoid

    std::cout << "frobenius norm of ref code u=" << n_u << std::endl;
    std::cout << "frobenius norm of code to test U=" << n_U << std::endl;
    BOOST_TEST(n_u == n_U);

    std::cout << "frobenius norm of ref code dudx=" << n_dudx << std::endl;
    std::cout << "frobenius norm of code to test dUdx=" << n_dUdx << std::endl;
    BOOST_TEST(n_dudx == n_dUdx);

    std::cout << "frobenius norm of ref code dudy=" << n_dudy << std::endl;
    std::cout << "frobenius norm of code to test dUdy=" << n_dUdy << std::endl;
    BOOST_TEST(n_dudy == n_dUdy);

    std::cout << "frobenius norm of ref code dudz=" << n_dudz << std::endl;
    std::cout << "frobenius norm of code to test dUdz=" << n_dUdz << std::endl;
    BOOST_TEST(n_dudz == n_dUdz);

    std::cout << "frobenius norm of ref code v=" << n_v << std::endl;
    std::cout << "frobenius norm of code to test V=" << n_V << std::endl;
    BOOST_TEST(n_v == n_V);

    std::cout << "frobenius norm of ref code dvdx=" << n_dvdx << std::endl;
    std::cout << "frobenius norm of code to test dVdx=" << n_dVdx << std::endl;
    BOOST_TEST(n_dvdx == n_dVdx);

    std::cout << "frobenius norm of ref code dvdy=" << n_dvdy << std::endl;
    std::cout << "frobenius norm of code to test dVdy=" << n_dVdy << std::endl;
    BOOST_TEST(n_dvdy == n_dVdy);

    std::cout << "frobenius norm of ref code dvdz=" << n_dvdz << std::endl;
    std::cout << "frobenius norm of code to test dVdz=" << n_dVdz << std::endl;
    BOOST_TEST(n_dvdz == n_dVdz);
    }

BOOST_AUTO_TEST_CASE(Tet_nod_interpolation2, *boost::unit_test::tolerance(UT_TOL))
    {
    using namespace Nodes;
    const int nbNod = 4;
    std::vector<Nodes::Node> node;
    dummyNodes<nbNod>(node);

    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);

    for (int i = 0; i < nbNod; i++)
        {
        node[i].d[0].phi = distrib(gen);
        node[i].d[0].phiv = distrib(gen);
        }

    // carefull with indices (starting from 1)
    Tetra::Tet t(node, 0, {1, 2, 3, 4});

    double t_dadx[Tetra::N][Tetra::NPI],t_dady[Tetra::N][Tetra::NPI],t_dadz[Tetra::N][Tetra::NPI];
    
    for(int i=0;i<Tetra::N;i++)
        for(int j=0;j<Tetra::NPI;j++)
            {
            t_dadx[i][j] = t.da(i,0);
            t_dady[i][j] = t.da(i,1);
            t_dadz[i][j] = t.da(i,2);
            }

    // ref code (with minimal adaptations of integrales method in file MuMag_integrales.cc of
    // src_Tube_scalfmm_thiaville_ec_mu_oersted_thiele_dyn20180903.tgz )
    double negphi0_nod[Tetra::N], Hdx[Tetra::NPI], Hdy[Tetra::NPI], Hdz[Tetra::NPI];
    double negphiv0_nod[Tetra::N], Hvx[Tetra::NPI], Hvy[Tetra::NPI], Hvz[Tetra::NPI];

    for (int ie = 0; ie < Tetra::N; ie++)
        {
        int i = t.ind[ie];
        Nodes::dataNode &d0 = node[i].d[0];

        negphi0_nod[ie] = -d0.phi;
        negphiv0_nod[ie] = -d0.phiv;
        }

    tiny::transposed_mult<double, Tetra::N, Tetra::NPI>(negphi0_nod, t_dadx, Hdx);
    tiny::transposed_mult<double, Tetra::N, Tetra::NPI>(negphi0_nod, t_dady, Hdy);
    tiny::transposed_mult<double, Tetra::N, Tetra::NPI>(negphi0_nod, t_dadz, Hdz);

    tiny::transposed_mult<double, Tetra::N, Tetra::NPI>(negphiv0_nod, t_dadx, Hvx);
    tiny::transposed_mult<double, Tetra::N, Tetra::NPI>(negphiv0_nod, t_dady, Hvy);
    tiny::transposed_mult<double, Tetra::N, Tetra::NPI>(negphiv0_nod, t_dadz, Hvz);
    // end ref code

    // code to check
    Eigen::Matrix<double,DIM,Tetra::NPI> Hd, Hv;

    t.interpolation_field(Nodes::get_phi<Nodes::CURRENT>, Hd);
    t.interpolation_field(Nodes::get_phiv<Nodes::CURRENT>, Hv);
    // end code to check

    double n_Hdx = tiny::frob_norm<double, Tetra::NPI>(Hdx);
    double n_Hdy = tiny::frob_norm<double, Tetra::NPI>(Hdy);
    double n_Hdz = tiny::frob_norm<double, Tetra::NPI>(Hdz);
    double n_Hd_ref = sqrt(sq(n_Hdx) + sq(n_Hdy) + sq(n_Hdz));

    double n_Hvx = tiny::frob_norm<double, Tetra::NPI>(Hvx);
    double n_Hvy = tiny::frob_norm<double, Tetra::NPI>(Hvy);
    double n_Hvz = tiny::frob_norm<double, Tetra::NPI>(Hvz);
    double n_Hv_ref = sqrt(sq(n_Hvx) + sq(n_Hvy) + sq(n_Hvz));

    double dist_Hd = sq_dist(Hdx, Hdy, Hdz, Hd);
    double dist_Hv = sq_dist(Hvx, Hvy, Hvz, Hv);

    double n_Hd = sqrt(sq_frobenius_norm<Tetra::NPI>(Hd));
    double n_Hv = sqrt(sq_frobenius_norm<Tetra::NPI>(Hv));

    if (!DET_UT) std::cout << "seed =" << sd << std::endl;
    std::cout << "distance^2 Hd =" << dist_Hd << std::endl;
    BOOST_TEST(sqrt(dist_Hd) == 0.0);

    std::cout << "distance^2 Hv =" << dist_Hv << std::endl;
    BOOST_TEST(sqrt(dist_Hv) == 0.0);

    std::cout << "frobenius norm of ref code Hd=" << n_Hd_ref << std::endl;
    std::cout << "frobenius norm of code to test Hd=" << n_Hd << std::endl;
    BOOST_TEST(n_Hd_ref == n_Hd);

    std::cout << "frobenius norm of ref code Hv=" << n_Hv_ref << std::endl;
    std::cout << "frobenius norm of code to test Hv=" << n_Hv << std::endl;
    BOOST_TEST(n_Hv_ref == n_Hv);
    }

BOOST_AUTO_TEST_CASE(Tet_Pcoeff)
    {
    std::cout << "Tet test on Tetra::buildMatP()" << std::endl;
    const int N = Tetra::N;

    const int nbNod = 4;
    std::vector<Nodes::Node> node;
    dummyNodes<nbNod>(node);

    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);

    for (int i = 0; i < nbNod; i++)
        {
        node[i].d[0].u = rand_vec3d(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
        node[i].setBasis(2 * M_PI * distrib(gen));
        }

    Tetra::Tet t(node, 0, {1, 2, 3, 4});// carefull with indices (starting from 1)
    Eigen::Matrix<double,2*N,3*N> P;
    t.buildMatP(P);

    /* ref code */
    double Pref[2 * N][3 * N] = {{0}};  // P must be filled with zero

    for (int i = 0; i < N; i++)
        {
        const Eigen::Vector3d &ep = node[t.ind[i]].ep;
        Pref[i][i] = ep.x();
        Pref[i][N + i] = ep.y();
        Pref[i][2 * N + i] = ep.z();
        const Eigen::Vector3d &eq = node[t.ind[i]].eq;
        Pref[N + i][i] = eq.x();
        Pref[N + i][N + i] = eq.y();
        Pref[N + i][2 * N + i] = eq.z();
        }
    /* end ref code */

    for (int i = 0; i < 2 * N; i++)
        for (int j = 0; j < 3 * N; j++)
            {
            BOOST_CHECK(P(i,j) == Pref[i][j]);
            }
    }

/* test of two equivalent formulas for a gradient of V scalar potential on the nodes */
BOOST_AUTO_TEST_CASE(Tet_gradV, *boost::unit_test::tolerance(UT_TOL))
    {
    using namespace Nodes;
    const int nbNod = 4;
    std::vector<double> V(nbNod); // scalar potential
    std::vector<Nodes::Node> node;
    dummyNodes<nbNod>(node);


    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);
    std::cout << "test on gradV with a tetrahedron\n V = ";
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

    // code to test
    Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> gradV = Tetra::calc_gradV(tet,V);

    // ugly version
    Eigen::Matrix<double,Tetra::N,1> V_nod;
    for (size_t ie=0; ie<Tetra::N; ie++) { V_nod[ie] = V[ tet.ind[ie] ]; }
    Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> gradVbis;
    Eigen::Matrix<double,Tetra::N,Tetra::NPI> dadx;
    dadx.colwise() = tet.da.col(IDX_X); // colwise() means da.col(IDX_X) is repeated to build dadx
    Eigen::Matrix<double,Tetra::N,Tetra::NPI> dady;
    dady.colwise() = tet.da.col(IDX_Y);
    Eigen::Matrix<double,Tetra::N,Tetra::NPI> dadz;
    dadz.colwise() = tet.da.col(IDX_Z);
    // building explicitely dad(x|y|z) migth be avoided rewritting the following multiplications
    gradVbis.row(IDX_X) = V_nod.transpose() * dadx;// V_nod^T * dadx
    gradVbis.row(IDX_Y) = V_nod.transpose() * dady;// V_nod^T * dady
    gradVbis.row(IDX_Z) = V_nod.transpose() * dadz;// V_nod^T * dadz

    for(int i=0;i<Nodes::DIM;i++)
        for(int j=0;j<Tetra::NPI;j++)
            {
            std::cout << gradV(i,j) << " should be " << gradVbis(i,j) << std::endl;
            BOOST_TEST( gradV(i,j) == gradVbis(i,j) );
            }
    }

BOOST_AUTO_TEST_SUITE_END()
