#define BOOST_TEST_MODULE tetraTest

#include <boost/test/unit_test.hpp>

#include <random>

#include "tetra.h"
#include "ut_config.h"

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

// Tetra::Tet constructor is tested in ut_element.cpp

/*---------------------------------------*/
/* first lvl tests : nested calculus,... */
/*---------------------------------------*/

/*---------------------------------------*/
/* second lvl tests : pure mathematics   */
/*---------------------------------------*/
BOOST_AUTO_TEST_CASE(Tet_inner_tables, *boost::unit_test::tolerance(UT_TOL))
    {
    // this test is dedicated to  check dadx,dady,dadz and weight tables, those values are
    // initilized once by Tet constructor
    std::cout << "constructor test with 4 nodes in node vector\n";
    const int nbNod = 4;
    // std::shared_ptr<Nodes::Node[]> node(new Nodes::Node[nbNod]);
    std::vector<Nodes::Node> node(nbNod);

    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);

    Pt::pt3D p0(0, 0, 0), p1(1, 0, 0), p2(0, 1, 0), p3(0, 0, 1), u0(0, 0, 0), v0(0, 0, 0),
            u(0, 0, 0), v(0, 0, 0);
    double phi0(0), phi(0), phiv0(0), phiv(0);

    Nodes::Node n0 = {p0,   u0,  v0,    u,   v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0),
                      phi0, phi, phiv0, phiv};
    Nodes::Node n1 = {p1,   u0,  v0,    u,   v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0),
                      phi0, phi, phiv0, phiv};
    Nodes::Node n2 = {p2,   u0,  v0,    u,   v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0),
                      phi0, phi, phiv0, phiv};
    Nodes::Node n3 = {p3,   u0,  v0,    u,   v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0),
                      phi0, phi, phiv0, phiv};

    node[0] = n0;
    node[1] = n1;
    node[2] = n2;
    node[3] = n3;
    for (int i = 0; i < nbNod; i++)
        {
        node[i].u0 = Pt::pt3D(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
        }

    // carefull with indices (starting from 1)
    Tetra::Tet t(node, 0, {1, 2, 3, 4});
    t.infos();

    // ref code (with minimal adaptations of dad(x|y|z) in file Mesh_hat.cc of
    // src_Tube_scalfmm_thiaville_ec_mu_oersted_thiele_dyn20180903.tgz )
    double _dadx[Tetra::N][Tetra::NPI];
    double _dady[Tetra::N][Tetra::NPI];
    double _dadz[Tetra::N][Tetra::NPI];
    double da[Tetra::N][Pt::DIM];
    double J[Pt::DIM][Pt::DIM];
    double nod[Pt::DIM][Tetra::N];
    double weight[Tetra::NPI];

    for (int ie = 0; ie < Tetra::N; ie++)
        {
        int i = t.ind[ie];
        nod[0][ie] = node[i].p(0);
        nod[1][ie] = node[i].p(1);
        nod[2][ie] = node[i].p(2);
        }

    tiny::mult<double, Pt::DIM, Tetra::N, Pt::DIM>(nod, Tetra::dadu, J);
    double detJ = Pt::det(J);  // lu_det(J);
    Pt::inverse(J, detJ);
    tiny::mult<double, Tetra::N, Pt::DIM, Pt::DIM>(Tetra::dadu, J, da);

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
            result_dadx += Pt::sq(_dadx[ie][npi] - t.dadx[ie][npi]);
            result_dady += Pt::sq(_dady[ie][npi] - t.dady[ie][npi]);
            result_dadz += Pt::sq(_dadz[ie][npi] - t.dadz[ie][npi]);
            }
        result_w += Pt::sq(weight[npi] - t.weight[npi]);
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
    int nbNod = 4;
    std::vector<Nodes::Node> node(nbNod);

    Pt::pt3D p0(0, 0, 0), p1(1, 0, 0), p2(0, 1, 0), p3(0, 0, 1), u0(0, 0, 0), v0(0, 0, 0),
            u(0, 0, 0), v(0, 0, 0);
    double phi0(0), phi(0), phiv0(0), phiv(0);

    Nodes::Node n0 = {p0,   u0,  v0,    u,   v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0),
                      phi0, phi, phiv0, phiv};
    Nodes::Node n1 = {p1,   u0,  v0,    u,   v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0),
                      phi0, phi, phiv0, phiv};
    Nodes::Node n2 = {p2,   u0,  v0,    u,   v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0),
                      phi0, phi, phiv0, phiv};
    Nodes::Node n3 = {p3,   u0,  v0,    u,   v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0),
                      phi0, phi, phiv0, phiv};

    node[0] = n0;
    node[1] = n1;
    node[2] = n2;
    node[3] = n3;

    // carefull with indices (starting from 1)
    Tetra::Tet t(node, 0, {1, 2, 3, 4});

    double result = 1 / 6.0;
    double vol = t.calc_vol();
    std::cout << "vol(tetra) =" << vol << std::endl;
    BOOST_TEST(vol == result);
    }

double sq_dist(double _x[Pt::DIM][Tetra::NPI], Pt::pt3D X[Tetra::NPI])
    {
    double val(0.0);

    for (int i = 0; i < Tetra::N; i++)
        for (int j = 0; j < Pt::DIM; j++)
            {
            val += Pt::sq(_x[j][i] - X[i](j));
            }

    return val;
    }

double sq_dist(double _x[Tetra::NPI], double _y[Tetra::NPI], double _z[Tetra::NPI],
               Pt::pt3D X[Tetra::NPI])
    {
    double val(0.0);

    for (int i = 0; i < Tetra::N; i++)
        {
        val += Pt::sq(_x[i] - X[i].x());
        val += Pt::sq(_y[i] - X[i].y());
        val += Pt::sq(_z[i] - X[i].z());
        }

    return val;
    }

BOOST_AUTO_TEST_CASE(Tet_nod_interpolation, *boost::unit_test::tolerance(UT_TOL))
    {
    int nbNod = 4;
    std::vector<Nodes::Node> node(nbNod);

    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);

    Pt::pt3D p0(0, 0, 0), p1(1, 0, 0), p2(0, 1, 0), p3(0, 0, 1), u0(0, 0, 0), v0(0, 0, 0),
            u(0, 0, 0), v(0, 0, 0);
    double phi0(0), phi(0), phiv0(0), phiv(0);

    Nodes::Node n0 = {p0,   u0,  v0,    u,   v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0),
                      phi0, phi, phiv0, phiv};
    Nodes::Node n1 = {p1,   u0,  v0,    u,   v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0),
                      phi0, phi, phiv0, phiv};
    Nodes::Node n2 = {p2,   u0,  v0,    u,   v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0),
                      phi0, phi, phiv0, phiv};
    Nodes::Node n3 = {p3,   u0,  v0,    u,   v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0),
                      phi0, phi, phiv0, phiv};

    node[0] = n0;
    node[1] = n1;
    node[2] = n2;
    node[3] = n3;

    for (int i = 0; i < nbNod; i++)
        {
        node[i].u0 = Pt::pt3D(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
        node[i].v0 = Pt::pt3D(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
        }

    // carefull with indices (starting from 1)
    Tetra::Tet t(node, 0, {1, 2, 3, 4});

    // ref code (with minimal adaptations of integrales method in file MuMag_integrales.cc of
    // src_Tube_scalfmm_thiaville_ec_mu_oersted_thiele_dyn20180903.tgz )
    double _u_nod[3][Tetra::N], _u[3][Tetra::NPI];
    double dudx[3][Tetra::NPI], dudy[3][Tetra::NPI], dudz[3][Tetra::NPI];

    double _v_nod[3][Tetra::N], _v[3][Tetra::NPI];
    double dvdx[3][Tetra::NPI], dvdy[3][Tetra::NPI], dvdz[3][Tetra::NPI];

    for (int ie = 0; ie < Tetra::N; ie++)
        {
        int i = t.ind[ie];
        Nodes::Node &nod = node[i];
        for (int d = 0; d < 3; d++)
            {
            _u_nod[d][ie] = nod.u0(d);
            _v_nod[d][ie] = nod.v0(d);
            }
        }
    tiny::mult<double, 3, Tetra::N, Tetra::NPI>(_u_nod, Tetra::a, _u);
    tiny::mult<double, 3, Tetra::N, Tetra::NPI>(_u_nod, t.dadx, dudx);
    tiny::mult<double, 3, Tetra::N, Tetra::NPI>(_u_nod, t.dady, dudy);
    tiny::mult<double, 3, Tetra::N, Tetra::NPI>(_u_nod, t.dadz, dudz);

    tiny::mult<double, 3, Tetra::N, Tetra::NPI>(_v_nod, Tetra::a, _v);
    tiny::mult<double, 3, Tetra::N, Tetra::NPI>(_v_nod, t.dadx, dvdx);
    tiny::mult<double, 3, Tetra::N, Tetra::NPI>(_v_nod, t.dady, dvdy);
    tiny::mult<double, 3, Tetra::N, Tetra::NPI>(_v_nod, t.dadz, dvdz);
    // end ref code

    // code to check
    Pt::pt3D dUdx[Tetra::NPI], dUdy[Tetra::NPI], dUdz[Tetra::NPI];
    Pt::pt3D dVdx[Tetra::NPI], dVdy[Tetra::NPI], dVdz[Tetra::NPI];
    Pt::pt3D U[Tetra::NPI], V[Tetra::NPI];

    t.interpolation(Nodes::get_u0, U, dUdx, dUdy, dUdz);
    t.interpolation(Nodes::get_v0, V, dVdx, dVdy, dVdz);
    // end code to check

    double n_u = tiny::frob_norm<double, 3, Tetra::NPI>(_u);
    double n_dudx = tiny::frob_norm<double, 3, Tetra::NPI>(dudx);
    double n_dudy = tiny::frob_norm<double, 3, Tetra::NPI>(dudy);
    double n_dudz = tiny::frob_norm<double, 3, Tetra::NPI>(dudz);

    double n_v = tiny::frob_norm<double, 3, Tetra::NPI>(_v);
    double n_dvdx = tiny::frob_norm<double, 3, Tetra::NPI>(dvdx);
    double n_dvdy = tiny::frob_norm<double, 3, Tetra::NPI>(dvdy);
    double n_dvdz = tiny::frob_norm<double, 3, Tetra::NPI>(dvdz);

    double dist_uU = sq_dist(_u, U);
    double dist_dudx_dUdx = sq_dist(dudx, dUdx);
    double dist_dudy_dUdy = sq_dist(dudy, dUdy);
    double dist_dudz_dUdz = sq_dist(dudz, dUdz);

    double dist_vV = sq_dist(_v, V);
    double dist_dvdx_dVdx = sq_dist(dvdx, dVdx);
    double dist_dvdy_dVdy = sq_dist(dvdy, dVdy);
    double dist_dvdz_dVdz = sq_dist(dvdz, dVdz);

    double n_U = sqrt(Pt::sq_frobenius_norm<Tetra::NPI>(U));
    double n_dUdx = sqrt(Pt::sq_frobenius_norm<Tetra::NPI>(dUdx));
    double n_dUdy = sqrt(Pt::sq_frobenius_norm<Tetra::NPI>(dUdy));
    double n_dUdz = sqrt(Pt::sq_frobenius_norm<Tetra::NPI>(dUdz));

    double n_V = sqrt(Pt::sq_frobenius_norm<Tetra::NPI>(V));
    double n_dVdx = sqrt(Pt::sq_frobenius_norm<Tetra::NPI>(dVdx));
    double n_dVdy = sqrt(Pt::sq_frobenius_norm<Tetra::NPI>(dVdy));
    double n_dVdz = sqrt(Pt::sq_frobenius_norm<Tetra::NPI>(dVdz));

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
    int nbNod = 4;
    std::vector<Nodes::Node> node(nbNod);

    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);

    Pt::pt3D p0(0, 0, 0), p1(1, 0, 0), p2(0, 1, 0), p3(0, 0, 1), u0(0, 0, 0), v0(0, 0, 0),
            u(0, 0, 0), v(0, 0, 0);
    double phi0(0), phi(0), phiv0(0), phiv(0);

    Nodes::Node n0 = {p0,   u0,  v0,    u,   v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0),
                      phi0, phi, phiv0, phiv};
    Nodes::Node n1 = {p1,   u0,  v0,    u,   v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0),
                      phi0, phi, phiv0, phiv};
    Nodes::Node n2 = {p2,   u0,  v0,    u,   v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0),
                      phi0, phi, phiv0, phiv};
    Nodes::Node n3 = {p3,   u0,  v0,    u,   v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0),
                      phi0, phi, phiv0, phiv};

    node[0] = n0;
    node[1] = n1;
    node[2] = n2;
    node[3] = n3;

    for (int i = 0; i < nbNod; i++)
        {
        node[i].phi0 = distrib(gen);
        node[i].phiv0 = distrib(gen);
        }

    // carefull with indices (starting from 1)
    Tetra::Tet t(node, 0, {1, 2, 3, 4});

    // ref code (with minimal adaptations of integrales method in file MuMag_integrales.cc of
    // src_Tube_scalfmm_thiaville_ec_mu_oersted_thiele_dyn20180903.tgz )
    double negphi0_nod[Tetra::N], Hdx[Tetra::NPI], Hdy[Tetra::NPI], Hdz[Tetra::NPI];
    double negphiv0_nod[Tetra::N], Hvx[Tetra::NPI], Hvy[Tetra::NPI], Hvz[Tetra::NPI];

    for (int ie = 0; ie < Tetra::N; ie++)
        {
        int i = t.ind[ie];
        Nodes::Node &nod = node[i];

        negphi0_nod[ie] = -nod.phi0;
        negphiv0_nod[ie] = -nod.phiv0;
        }

    tiny::transposed_mult<double, Tetra::N, Tetra::NPI>(negphi0_nod, t.dadx, Hdx);
    tiny::transposed_mult<double, Tetra::N, Tetra::NPI>(negphi0_nod, t.dady, Hdy);
    tiny::transposed_mult<double, Tetra::N, Tetra::NPI>(negphi0_nod, t.dadz, Hdz);

    tiny::transposed_mult<double, Tetra::N, Tetra::NPI>(negphiv0_nod, t.dadx, Hvx);
    tiny::transposed_mult<double, Tetra::N, Tetra::NPI>(negphiv0_nod, t.dady, Hvy);
    tiny::transposed_mult<double, Tetra::N, Tetra::NPI>(negphiv0_nod, t.dadz, Hvz);
    // end ref code

    // code to check
    Pt::pt3D Hd[Tetra::NPI], Hv[Tetra::NPI];

    t.interpolation(Nodes::get_phi0, Hd);
    t.interpolation(Nodes::get_phiv0, Hv);
    // end code to check

    double n_Hdx = tiny::frob_norm<double, Tetra::NPI>(Hdx);
    double n_Hdy = tiny::frob_norm<double, Tetra::NPI>(Hdy);
    double n_Hdz = tiny::frob_norm<double, Tetra::NPI>(Hdz);
    double n_Hd_ref = sqrt(Pt::sq(n_Hdx) + Pt::sq(n_Hdy) + Pt::sq(n_Hdz));

    double n_Hvx = tiny::frob_norm<double, Tetra::NPI>(Hvx);
    double n_Hvy = tiny::frob_norm<double, Tetra::NPI>(Hvy);
    double n_Hvz = tiny::frob_norm<double, Tetra::NPI>(Hvz);
    double n_Hv_ref = sqrt(Pt::sq(n_Hvx) + Pt::sq(n_Hvy) + Pt::sq(n_Hvz));

    double dist_Hd = sq_dist(Hdx, Hdy, Hdz, Hd);
    double dist_Hv = sq_dist(Hvx, Hvy, Hvz, Hv);

    double n_Hd = sqrt(Pt::sq_frobenius_norm<Tetra::NPI>(Hd));
    double n_Hv = sqrt(Pt::sq_frobenius_norm<Tetra::NPI>(Hv));

    if (!DET_UT) std::cout << "seed =" << sd << std::endl;
    std::cout << "distance^2 Hd =" << dist_Hd << std::endl;
    BOOST_TEST(sqrt(dist_Hd) == 0.0);

    std::cout << "distance^2 Hv =" << dist_Hv << std::endl;
    BOOST_TEST(sqrt(dist_Hv) == 0.0);

    // to avoid gag of comparing pure zeros we also check that matrices norm are equal
    // let's be paranoid

    std::cout << "frobenius norm of ref code Hd=" << n_Hd_ref << std::endl;
    std::cout << "frobenius norm of code to test Hd=" << n_Hd << std::endl;
    BOOST_TEST(n_Hd_ref == n_Hd);

    std::cout << "frobenius norm of ref code Hv=" << n_Hv_ref << std::endl;
    std::cout << "frobenius norm of code to test Hv=" << n_Hv << std::endl;
    BOOST_TEST(n_Hv_ref == n_Hv);
    }

BOOST_AUTO_TEST_CASE(Tet_lumping, *boost::unit_test::tolerance(UT_TOL))
    {
    double AE[3 * Tetra::N][3 * Tetra::N] = {{0}};
    double AE_to_check[3 * Tetra::N][3 * Tetra::N] = {{0}};

    int nbNod = 4;
    std::vector<Nodes::Node> node(nbNod);

    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);

    Pt::pt3D p0(0, 0, 0), p1(1, 0, 0), p2(0, 1, 0), p3(0, 0, 1), u0(0, 0, 0), v0(0, 0, 0),
            u(0, 0, 0), v(0, 0, 0);
    double phi0(0), phi(0), phiv0(0), phiv(0);

    Nodes::Node n0 = {p0,   u0,  v0,    u,   v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0),
                      phi0, phi, phiv0, phiv};
    Nodes::Node n1 = {p1,   u0,  v0,    u,   v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0),
                      phi0, phi, phiv0, phiv};
    Nodes::Node n2 = {p2,   u0,  v0,    u,   v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0),
                      phi0, phi, phiv0, phiv};
    Nodes::Node n3 = {p3,   u0,  v0,    u,   v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0),
                      phi0, phi, phiv0, phiv};

    node[0] = n0;
    node[1] = n1;
    node[2] = n2;
    node[3] = n3;

    for (int i = 0; i < nbNod; i++)
        {
        node[i].u0 = Pt::pt3D(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
        }

    // carefull with indices (starting from 1)
    Tetra::Tet t(node, 0, {1, 2, 3, 4});
    
    double a = distrib(gen);
    double b = distrib(gen);

    timing prm_t = timing(1.0, std::min(a, b), std::max(a, b));

    double s_dt = THETA * prm_t.get_dt();  // theta from theta scheme in config.h.in

    double alpha_LLG = 0.5;
    double uHeff = distrib(gen);
    double alfa = prm_t.calc_alpha_eff(alpha_LLG, uHeff);

    double dt = prm_t.get_dt();
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
        Nodes::Node &nod = node[i];
        for (int d = 0; d < 3; d++)
            {
            u_nod[d][ie] = nod.u0(d);  // aimantation
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
            dai_dx = t.dadx[ie][npi];
            dai_dy = t.dady[ie][npi];
            dai_dz = t.dadz[ie][npi];

            AE[ie][ie] += alfa * ai * w;  // lumping
            AE[Tetra::N + ie][Tetra::N + ie] += alfa * ai * w;
            AE[2 * Tetra::N + ie][2 * Tetra::N + ie] += alfa * ai * w;

            AE[ie][2 * Tetra::N + ie] += +u_nod[1][ie] * ai * w;  // lumping
            AE[ie][Tetra::N + ie] += -u_nod[2][ie] * ai * w;
            AE[Tetra::N + ie][ie] += +u_nod[2][ie] * ai * w;
            AE[Tetra::N + ie][2 * Tetra::N + ie] += -u_nod[0][ie] * ai * w;
            AE[2 * Tetra::N + ie][Tetra::N + ie] += +u_nod[0][ie] * ai * w;
            AE[2 * Tetra::N + ie][ie] += -u_nod[1][ie] * ai * w;

            for (int je = 0; je < Tetra::N; je++)
                {
                daj_dx = t.dadx[je][npi];
                daj_dy = t.dady[je][npi];
                daj_dz = t.dadz[je][npi];
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
    for (int npi = 0; npi < Tetra::NPI; npi++)
        t.lumping(npi, alfa, prm_t.prefactor * s_dt * Abis, AE_to_check);
    // end code to check

    double val = tiny::dist<double, 3 * Tetra::N, 3 * Tetra::N>(AE, AE_to_check);
    if (!DET_UT) std::cout << "seed =" << sd << std::endl;
    std::cout << "distance = " << val << std::endl;
    BOOST_TEST(val == 0.0);
    }

BOOST_AUTO_TEST_CASE(Tet_Pcoeff)
    {
    std::cout << "Tet test on Nodes::Pcoeff template" << std::endl;
    const int N = Tetra::N;

    int nbNod = 4;
    std::vector<Nodes::Node> node(nbNod);

    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);

    Pt::pt3D p1(1, 0, 0), p2(0, 1, 0), p3(1, 1, 0), p4(0, 0, 1), u0(0, 0, 0), v0(0, 0, 0),
            u(0, 0, 0), v(0, 0, 0);
    double phi0(0), phi(0), phiv0(0), phiv(0);

    node[0] = {p1, u0, v0, u, v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0), phi0, phi, phiv0, phiv};
    node[1] = {p2, u0, v0, u, v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0), phi0, phi, phiv0, phiv};
    node[2] = {p3, u0, v0, u, v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0), phi0, phi, phiv0, phiv};
    node[3] = {p4, u0, v0, u, v, Pt::pt3D(0, 0, 0), Pt::pt3D(0, 0, 0), phi0, phi, phiv0, phiv};

    for (int i = 0; i < nbNod; i++)
        {
        node[i].u0 = Pt::pt3D(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
        node[i].setBasis(2 * M_PI * distrib(gen));
        }

    // carefull with indices (starting from 1)
    Tetra::Tet t(node, 0, {1, 2, 3, 4});

    double P[2 * N][3 * N] = {{0}};
    for (int i = 0; i < 2 * N; i++)
        for (int j = 0; j < 3 * N; j++)
            {
            P[i][j] = t.Pcoeff(i, j);
            }

    /* ref code */
    double Pref[2 * N][3 * N] = {{0}};  // P must be filled with zero

    for (int i = 0; i < N; i++)
        {
        const Nodes::Node &n = t.getNode(i);
        Pref[i][i] = n.ep.x();
        Pref[i][N + i] = n.ep.y();
        Pref[i][2 * N + i] = n.ep.z();
        Pref[N + i][i] = n.eq.x();
        Pref[N + i][N + i] = n.eq.y();
        Pref[N + i][2 * N + i] = n.eq.z();
        }
    /* end ref code */
    double normP = tiny::frob_norm<double, 2 * N, 3 * N>(P);
    if (!DET_UT) std::cout << "seed =" << sd << std::endl;

    std::cout << "frob norm(P) = " << normP
              << " ; frob norm(Pref) = " << tiny::frob_norm<double, 2 * N, 3 * N>(Pref)
              << std::endl;
    double result = tiny::dist<double, 2 * N, 3 * N>(P, Pref);
    std::cout << "dist(P,Pref) = " << result << std::endl;

    BOOST_CHECK(normP > ((double) 0));
    BOOST_CHECK(isfinite(result));
    for (int i = 0; i < 2 * N; i++)
        for (int j = 0; j < 3 * N; j++)
            {
            BOOST_CHECK(P[i][j] == Pref[i][j]);
            }
    }

BOOST_AUTO_TEST_SUITE_END()
