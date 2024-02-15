#define BOOST_TEST_MODULE anisotropyTest

#include <boost/test/unit_test.hpp>

#include <random>

#include <eigen3/Eigen/Dense>

#include "tetra.h"
#include "tiny.h"
#include "ut_tools.h"
#include "ut_config.h"

/**
 The purpose of these unit test is to check that the formulation of the anisotropies are formulated
 the same in reference code ref code (in file MuMag_integrales.cc of
 src_Tube_scalfmm_thiaville_ec_mu_oersted_thiele_dyn20180903.tgz ) and feellgood 'public'
 */

BOOST_AUTO_TEST_SUITE(ut_anisotropy)

BOOST_AUTO_TEST_CASE(anisotropy_uniax, *boost::unit_test::tolerance(10.0 * UT_TOL))
    {
    int nbNod = 4;
    std::vector<Nodes::Node> node(nbNod);

    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);

    Eigen::Vector3d p0(0, 0, 0), p1(1, 0, 0), p2(0, 1, 0), p3(0, 0, 1), u0(0, 0, 0), v0(0, 0, 0);
    Eigen::Vector3d zero(0,0,0);
    double phi0(0), phi(0), phiv0(0), phiv(0);

    Nodes::Node n0 = {p0,
                      u0,
                      v0,
                      zero,zero,zero,zero,
                      phi0,
                      phi,
                      phiv0,
                      phiv};
    Nodes::Node n1 = {p1,
                      u0,
                      v0,
                      zero,zero,zero,zero,
                      phi0,
                      phi,
                      phiv0,
                      phiv};
    Nodes::Node n2 = {p2,
                      u0,
                      v0,
                      zero,zero,zero,zero,
                      phi0,
                      phi,
                      phiv0,
                      phiv};
    Nodes::Node n3 = {p3,
                      u0,
                      v0,
                      zero,zero,zero,zero,
                      phi0,
                      phi,
                      phiv0,
                      phiv};

    node[0] = n0;
    node[1] = n1;
    node[2] = n2;
    node[3] = n3;

    for (int i = 0; i < nbNod; i++)
        {
        node[i].u0 = rand_vec3d(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
        node[i].v0 = rand_vec3d(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
        }
    // carefull with indices (starting from 1)
    Tetra::Tet t(node, 0, {1, 2, 3, 4});

    double dt = distrib(gen);

    Eigen::Vector3d uk = rand_vec3d(M_PI * distrib(gen), 2 * M_PI * distrib(gen));

    double uk00 = uk.x();
    double uk01 = uk.y();
    double uk02 = uk.z();
    double K = distrib(gen);
    double Js = 0.5 + distrib(gen);  // add 0.5 to center Js around 1

    double Kbis = 2.0 * K / Js;
    std::cout << "uniaxial anisotropy test on a tetrahedron" << std::endl;
    if (!DET_UT) std::cout << "seed =" << sd << std::endl;
    std::cout << "uk =" << uk << std::endl;


    // ref code (with minimal adaptations of integrales method in file MuMag_integrales.cc of
    // src_Tube_scalfmm_thiaville_ec_mu_oersted_thiele_dyn20180903.tgz )
    double u_nod[3][Tetra::N];
    double v_nod[3][Tetra::N];
    double u[3][Tetra::NPI];
    double v[3][Tetra::NPI];

    for (int ie = 0; ie < Tetra::N; ie++)
        {
        int i = t.ind[ie];
        Nodes::Node &nod = node[i];
        for (int d = 0; d < 3; d++)
            {
            u_nod[d][ie] = nod.u0(d);
            v_nod[d][ie] = nod.v0(d);
            }
        }

    tiny::mult<double, 3, Tetra::N, Tetra::NPI>(u_nod, Tetra::a, u);
    tiny::mult<double, 3, Tetra::N, Tetra::NPI>(v_nod, Tetra::a, v);

    Eigen::Matrix<double,Pt::DIM,Tetra::NPI> U;
    Eigen::Matrix<double,Pt::DIM,Tetra::NPI> V;
    for (int npi = 0; npi < Tetra::NPI; npi++)
        {
        U.col(npi) << u[0][npi], u[1][npi], u[2][npi];
        V.col(npi) << v[0][npi], v[1][npi], v[2][npi];
        }  // deep copy instead of re-computing

    //code to test
    Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> H_aniso;
    H_aniso.setZero();
    Eigen::Matrix<double,Tetra::NPI,1> contrib_aniso;
    contrib_aniso.setZero();
    contrib_aniso = t.calc_aniso_uniax(uk, Kbis, THETA * dt, U, V, H_aniso);
    // end of code to test

    for (int npi = 0; npi < Tetra::NPI; npi++)
        {
        double uk0_u = uk00 * u[0][npi] + uk01 * u[1][npi] + uk02 * u[2][npi];
        double uk0_v = uk00 * v[0][npi] + uk01 * v[1][npi] + uk02 * v[2][npi];

        double uHau = 2 * K / Js * uk0_u * uk0_u;

        double Ht[3];
        Ht[0] = 2 * K / Js * uk0_v * uk00;
        Ht[1] = 2 * K / Js * uk0_v * uk01;
        Ht[2] = 2 * K / Js * uk0_v * uk02;

        double H[3];
        H[0] = (2 * K / Js * uk0_u * uk00);
        H[1] = (2 * K / Js * uk0_u * uk01);
        H[2] = (2 * K / Js * uk0_u * uk02);

        Eigen::Vector3d refHval = Eigen::Vector3d( H[0] + THETA * dt * Ht[0], H[1] + THETA * dt * Ht[1],
                                    H[2] + THETA * dt * Ht[2] );
        std::cout << "ref H value = " << refHval << "; H_aniso=" << H_aniso << std::endl;

        BOOST_TEST( (refHval - H_aniso.col(npi)).norm() == 0.0,
                   "mismatch in uniaxial anisotropy field value");

        BOOST_TEST(uHau == contrib_aniso(npi));
        }
    }

BOOST_AUTO_TEST_CASE(anisotropy_cubic, *boost::unit_test::tolerance(10.0 * UT_TOL))
    {
    int nbNod = 4;
    std::vector<Nodes::Node> node(nbNod);

    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);

    Eigen::Vector3d p0(0, 0, 0), p1(1, 0, 0), p2(0, 1, 0), p3(0, 0, 1), u0(0, 0, 0), v0(0, 0, 0);
    Eigen::Vector3d zero(0,0,0);
    double phi0(0), phi(0), phiv0(0), phiv(0);

    Nodes::Node n0 = {p0,
                      u0,
                      v0,
                      zero,zero,zero,zero,
                      phi0,
                      phi,
                      phiv0,
                      phiv};
    Nodes::Node n1 = {p1,
                      u0,
                      v0,
                      zero,zero,zero,zero,
                      phi0,
                      phi,
                      phiv0,
                      phiv};
    Nodes::Node n2 = {p2,
                      u0,
                      v0,
                      zero,zero,zero,zero,
                      phi0,
                      phi,
                      phiv0,
                      phiv};
    Nodes::Node n3 = {p3,
                      u0,
                      v0,
                      zero,zero,zero,zero,
                      phi0,
                      phi,
                      phiv0,
                      phiv};

    node[0] = n0;
    node[1] = n1;
    node[2] = n2;
    node[3] = n3;

    for (int i = 0; i < nbNod; i++)
        {
        node[i].u0 = rand_vec3d(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
        node[i].v0 = rand_vec3d(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
        }
    // carefull with indices (starting from 1)
    Tetra::Tet t(node, 0, {1, 2, 3, 4});

    double dt = distrib(gen);

    Pt::pt3D tmp_rand_v = Pt::pt3D(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
    Eigen::Vector3d rand_vect { tmp_rand_v.x(), tmp_rand_v.y(), tmp_rand_v.z() };
    
    Pt::pt3D tmp(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
    Eigen::Vector3d ex {tmp.x(),tmp.y(),tmp.z()}; 
    Eigen::Vector3d ey = ex.cross(rand_vect);
    Eigen::Vector3d ez = ex.cross(ey);
    ey.normalize();
    ez.normalize();

    double uk00 = ex.x();
    double uk01 = ex.y();
    double uk02 = ex.z();

    double uk10 = ey.x();
    double uk11 = ey.y();
    double uk12 = ey.z();

    double uk20 = ez.x();
    double uk21 = ez.y();
    double uk22 = ez.z();

    double K3 = distrib(gen);
    double Js = 0.5 + distrib(gen);  // add 0.5 to center Js around 1

    double K3bis = 2.0 * K3 / Js;
    std::cout << "cubic anisotropy test on a tetrahedron" << std::endl;
    if (!DET_UT) std::cout << "seed =" << sd << std::endl;
    std::cout << "ex =" << ex << std::endl;
    std::cout << "ey =" << ey << std::endl;
    std::cout << "ez =" << ez << std::endl;

    // ref code (with minimal adaptations of integrales method in file MuMag_integrales.cc of
    // src_Tube_scalfmm_thiaville_ec_mu_oersted_thiele_dyn20180903.tgz )
    double u_nod[3][Tetra::N];
    double v_nod[3][Tetra::N];
    double u[3][Tetra::NPI];
    double v[3][Tetra::NPI];

    for (int ie = 0; ie < Tetra::N; ie++)
        {
        int i = t.ind[ie];
        Nodes::Node &nod = node[i];
        for (int d = 0; d < 3; d++)
            {
            u_nod[d][ie] = nod.u0(d);
            v_nod[d][ie] = nod.v0(d);
            }
        }

    tiny::mult<double, 3, Tetra::N, Tetra::NPI>(u_nod, Tetra::a, u);
    tiny::mult<double, 3, Tetra::N, Tetra::NPI>(v_nod, Tetra::a, v);

    Eigen::Matrix<double,Pt::DIM,Tetra::NPI> U,V;
    for (int npi = 0; npi < Tetra::NPI; npi++)
        {
        U.col(npi) << u[0][npi], u[1][npi], u[2][npi];
        V.col(npi) << v[0][npi], v[1][npi], v[2][npi];
        }  // deep copy instead of re-computing

        //code to test
        Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> H_aniso;
        H_aniso.setZero();
        Eigen::Matrix<double,Tetra::NPI,1> contrib_aniso;
        contrib_aniso.setZero();
        contrib_aniso = t.calc_aniso_cub(ex, ey, ez, K3bis, THETA * dt, U, V, H_aniso);
        // end of code to test

    for (int npi = 0; npi < Tetra::NPI; npi++)
        {
        double uk0_u = uk00 * u[0][npi] + uk01 * u[1][npi] + uk02 * u[2][npi];
        double uk1_u = uk10 * u[0][npi] + uk11 * u[1][npi] + uk12 * u[2][npi];
        double uk2_u = uk20 * u[0][npi] + uk21 * u[1][npi] + uk22 * u[2][npi];

        std::cout << "uk0_u = " << uk0_u << std::endl;
        std::cout << "uk1_u = " << uk1_u << std::endl;
        std::cout << "uk2_u = " << uk2_u << std::endl;

        double uk0_v = uk00 * v[0][npi] + uk01 * v[1][npi] + uk02 * v[2][npi];
        double uk1_v = uk10 * v[0][npi] + uk11 * v[1][npi] + uk12 * v[2][npi];
        double uk2_v = uk20 * v[0][npi] + uk21 * v[1][npi] + uk22 * v[2][npi];

        std::cout << "uk0_v = " << uk0_v << std::endl;
        std::cout << "uk1_v = " << uk1_v << std::endl;
        std::cout << "uk2_v = " << uk2_v << std::endl;

        double uHa3u = -2 * K3 / Js
                       * (uk0_u * (1 - uk0_u * uk0_u) * uk0_u + uk1_u * (1 - uk1_u * uk1_u) * uk1_u
                          + uk2_u * (1 - uk2_u * uk2_u) * uk2_u);

        double Ht[3];
        Ht[0] = -2 * K3 / Js * uk0_v * (1 - 3 * uk0_u * uk0_u) * uk00;
        Ht[1] = -2 * K3 / Js * uk1_v * (1 - 3 * uk1_u * uk1_u) * uk01;
        Ht[2] = -2 * K3 / Js * uk2_v * (1 - 3 * uk2_u * uk2_u) * uk02;

        double H[3];
        H[0] = -2 * K3 / Js * uk0_u * (1 - uk0_u * uk0_u) * uk00
               - 2 * K3 / Js * uk1_u * (1 - uk1_u * uk1_u) * uk10
               - 2 * K3 / Js * uk2_u * (1 - uk2_u * uk2_u) * uk20;
        H[1] = -2 * K3 / Js * uk0_u * (1 - uk0_u * uk0_u) * uk01
               - 2 * K3 / Js * uk1_u * (1 - uk1_u * uk1_u) * uk11
               - 2 * K3 / Js * uk2_u * (1 - uk2_u * uk2_u) * uk21;
        H[2] = -2 * K3 / Js * uk0_u * (1 - uk0_u * uk0_u) * uk02
               - 2 * K3 / Js * uk1_u * (1 - uk1_u * uk1_u) * uk12
               - 2 * K3 / Js * uk2_u * (1 - uk2_u * uk2_u) * uk22;

        Eigen::Vector3d refHval { H[0] + THETA*dt*Ht[0], H[1] + THETA*dt*Ht[1], H[2] + THETA*dt*Ht[2] };
        std::cout << "ref H value = " << refHval << "; H_aniso=" << H_aniso.col(npi) << std::endl;
        double result = (refHval - H_aniso.col(npi)).norm();
        std::cout << "result= " << result << std::endl;
        BOOST_TEST( (refHval - H_aniso.col(npi)).norm() == 0.0, "mismatch in cubic anisotropy field value");
        BOOST_TEST(uHa3u == contrib_aniso(npi), "mismatch in cubic anisotropy contrib_aniso value");
        }
    }

BOOST_AUTO_TEST_CASE(anisotropy_H_aniso, *boost::unit_test::tolerance(10.0 * UT_TOL))
    {
    double Kbis(.5);
    Eigen::Vector3d uk {1,3,5};
    Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> bob,H_a;
    H_a.setZero();
    bob << 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15;
    for(int npi = 0;npi<Tetra::NPI;npi++)
        { H_a.col(npi) += (Kbis * uk.dot( bob.col(npi))) * uk; }
    std::cout << "H_a= " << H_a << std::endl << "Kbis*uk.trans()*bob= " << Kbis*uk.transpose()*bob <<std::endl;

    Eigen::Array<double,Nodes::DIM,Tetra::NPI> tmp;
    tmp.colwise() = uk.array();
    Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> test = tmp.rowwise() * (Kbis*uk.transpose()*bob).array() ;
    std::cout << "test= " << test << std::endl;
    double dist = (H_a - test).norm();
    BOOST_TEST( dist == 0. );
    }

BOOST_AUTO_TEST_SUITE_END()
