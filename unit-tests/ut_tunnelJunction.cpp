#define BOOST_TEST_MODULE tunnelJunctionTest

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <iostream>
#include <random>

#include "tunnelJunction.h"
#include "ut_config.h"

BOOST_AUTO_TEST_SUITE(ut_tunnelJunction)

/** returns a unit vector pointing to (theta,phi) in spherical coordinates */
Eigen::Vector3d unit_vec(double theta, double phi)
    { return Eigen::Vector3d(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)); }

/** random unit vector */
Eigen::Vector3d rand_unit(std::mt19937 &gen)
    {
    std::uniform_real_distribution<> distrib(0.0, 1.0);
    return unit_vec(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
    }

BOOST_AUTO_TEST_CASE(conductance_P_AP)
    {
    /* G_P = 1/RA, G_AP = G_P/(1 + TMR) */
    Tunnel::Settings st;
    st.p = Eigen::Vector3d(0, 0, 1);
    st.RA = 1e-12;
    st.TMR = 1.5;
    const double GP = Tunnel::conductance(st, st.p);
    const double GAP = Tunnel::conductance(st, -st.p);
    BOOST_TEST(GP * st.RA == 1.0, boost::test_tools::tolerance(1e-14));
    BOOST_TEST((GP - GAP) / GAP == st.TMR, boost::test_tools::tolerance(1e-14));
    }

BOOST_AUTO_TEST_CASE(conductance_cosine_law)
    {
    /* G(theta) = (G_P + G_AP)/2 + (G_P - G_AP)/2 cos(theta) */
    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);
    if (!DET_UT) std::cout << "seed =" << sd << std::endl;

    Tunnel::Settings st;
    st.p = rand_unit(gen);
    st.RA = 1e-12 * (1 + distrib(gen));
    st.TMR = 5 * distrib(gen);
    const Eigen::Vector3d m0 = rand_unit(gen);
    const double GP = 1 / st.RA, GAP = GP / (1 + st.TMR);
    const double G = (GP + GAP) / 2 + (GP - GAP) / 2 * st.p.dot(m0);
    BOOST_TEST(Tunnel::conductance(st, m0) / G == 1.0, boost::test_tools::tolerance(1e-13));
    }

BOOST_AUTO_TEST_CASE(polarizations_default)
    {
    /* symmetric junction by default: P1 = P2 = sqrt(TMR/(2 + TMR)), 0.6547 for TMR = 150 % */
    Tunnel::Settings st;
    st.TMR = 1.5;
    const Tunnel::Polarizations P = Tunnel::polarizations(st);
    BOOST_TEST(P.P1 == std::sqrt(1.5 / 3.5), boost::test_tools::tolerance(1e-14));
    BOOST_TEST(P.P2 == P.P1, boost::test_tools::tolerance(1e-14));
    st.P1 = 0.8;
    const Tunnel::Polarizations Q = Tunnel::polarizations(st);
    BOOST_TEST(Q.P1 == 0.8, boost::test_tools::tolerance(1e-14));
    BOOST_TEST(Q.P1 * Q.P2 == 1.5 / 3.5, boost::test_tools::tolerance(1e-14));
    }

BOOST_AUTO_TEST_CASE(polarization_collinear)
    {
    /* Jullière: (P1 + P2)/(1 + P1 P2) along m0 in P, (P2 - P1)/(1 - P1 P2) in AP */
    Tunnel::Settings st;
    st.p = Eigen::Vector3d(0, 0, 1);
    st.TMR = 1.5;
    st.P1 = 0.7;
    const Tunnel::Polarizations P = Tunnel::polarizations(st);
    const Eigen::Vector3d PP = Tunnel::polarization(st, st.p);
    const Eigen::Vector3d PAP = Tunnel::polarization(st, -st.p);
    BOOST_TEST(PP.dot(st.p) == (P.P1 + P.P2) / (1 + P.P1 * P.P2),
               boost::test_tools::tolerance(1e-14));
    BOOST_TEST(PAP.dot(-st.p) == (P.P2 - P.P1) / (1 - P.P1 * P.P2),
               boost::test_tools::tolerance(1e-14));
    BOOST_TEST(PP.cross(st.p).norm() < 1e-15);
    BOOST_TEST(PAP.cross(st.p).norm() < 1e-15);
    }

BOOST_AUTO_TEST_CASE(polarization_angle)
    {
    /* P_T.m0 = (P1 cos + P2)/(1 + P1 P2 cos), transverse part P1 sin/(1 + P1 P2 cos), in the plane
     * (p, m0) */
    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);
    if (!DET_UT) std::cout << "seed =" << sd << std::endl;

    Tunnel::Settings st;
    st.p = rand_unit(gen);
    st.TMR = 0.1 + 5 * distrib(gen);
    const Tunnel::Polarizations P = Tunnel::polarizations(st);
    const Eigen::Vector3d m0 = rand_unit(gen);
    const double c = st.p.dot(m0), s = st.p.cross(m0).norm();
    const Eigen::Vector3d PT = Tunnel::polarization(st, m0);
    const double den = 1 + P.P1 * P.P2 * c;
    BOOST_TEST(PT.dot(m0) == (P.P1 * c + P.P2) / den, boost::test_tools::tolerance(1e-13));
    BOOST_TEST((PT - PT.dot(m0) * m0).norm() == P.P1 * s / den,
               boost::test_tools::tolerance(1e-12));
    BOOST_TEST(std::abs(PT.dot(st.p.cross(m0))) < 1e-14);
    }

BOOST_AUTO_TEST_SUITE_END()
