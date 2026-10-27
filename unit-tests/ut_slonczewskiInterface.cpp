#define BOOST_TEST_MODULE slonczewskiInterfaceTest

#include <boost/test/unit_test.hpp>

#include <cmath>
#include <iostream>
#include <random>

#include "slonczewskiInterface.h"
#include "ut_config.h"

BOOST_AUTO_TEST_SUITE(ut_slonczewskiInterface)

/** returns a unit vector pointing to (theta,phi) in spherical coordinates */
Eigen::Vector3d unit_vec(double theta, double phi)
    { return Eigen::Vector3d(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)); }

/** relative error of a with respect to the reference b (vectors or matrices) */
template <typename T>
double rel_err(const T &a, const T &b)
    { return (a - b).norm() / b.norm(); }

/** random unit vector */
Eigen::Vector3d rand_unit(std::mt19937 &gen)
    {
    std::uniform_real_distribution<> distrib(0.0, 1.0);
    return unit_vec(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
    }

BOOST_AUTO_TEST_CASE(conductance)
    {
    /* G w = Gpar (u.w) u + Gr (w - (u.w) u) - Gi u x w */
    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);
    if (!DET_UT) std::cout << "seed =" << sd << std::endl;

    const Eigen::Vector3d u = rand_unit(gen);
    const Eigen::Vector3d w(distrib(gen), distrib(gen), distrib(gen));
    const double Gpar = distrib(gen), Gr = distrib(gen), Gi = distrib(gen);

    const Eigen::Vector3d result = Slonczewski::conductance(u, Gpar, Gr, Gi) * w;
    const Eigen::Vector3d ref = Gpar * u.dot(w) * u + Gr * (w - u.dot(w) * u) - Gi * u.cross(w);
    BOOST_TEST( rel_err(result, ref) < 1e-14 );
    }

BOOST_AUTO_TEST_CASE(bulkConductance, *boost::unit_test::tolerance(1e-14))
    {
    const double D0 = 2.3e-3, P = 0.4, lsf = 10e-9;
    const double thick = D0 * (1 - P * P) / lsf;
    BOOST_TEST( Slonczewski::bulkConductance(D0, P, lsf, INFINITY) == thick );
    BOOST_TEST( Slonczewski::bulkConductance(D0, P, lsf, 2 * lsf) == thick / std::tanh(2.0) );
    }

BOOST_AUTO_TEST_CASE(renormalizeF1)
    {
    /* the renormalized interface gives the same spin current F1 -> N as the bulk and the interface
     solved for s_F1: beta K - gF s_F1 = P K + Gpar (s_F1 - sN), the flux being P* K - Gpar* sN */
    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);
    if (!DET_UT) std::cout << "seed =" << sd << std::endl;

    const double beta = distrib(gen), gF = 0.1 + distrib(gen), K = distrib(gen), sN = distrib(gen);
    const Slonczewski::F1prm itf = {distrib(gen), 0.1 + distrib(gen)};

    const double sF1 = ((beta - itf.P) * K + itf.Gpar * sN) / (gF + itf.Gpar);
    const double flux = beta * K - gF * sF1;

    const Slonczewski::F1prm r = Slonczewski::renormalizeF1(beta, gF, itf);
    // the flux may be small (two terms of opposite signs): error relative to the size of the terms
    const double scale = std::abs(beta * K) + std::abs(gF * sF1);
    BOOST_TEST( std::abs(r.P * K - r.Gpar * sN - flux) / scale < 1e-13 );

    // a bulk much more conducting than the interface: s_F1 = 0, the interface is unchanged
    const Slonczewski::F1prm r_inf = Slonczewski::renormalizeF1(beta, 1e12, itf);
    BOOST_TEST( r_inf.P == itf.P, boost::test_tools::tolerance(1e-10) );
    BOOST_TEST( r_inf.Gpar == itf.Gpar, boost::test_tools::tolerance(1e-10) );
    }

BOOST_AUTO_TEST_CASE(effectiveBC_spacer_node)
    {
    /* the effective condition gives the spin current of the explicit solution of the spacer node:
     P1 K p + G1 (sF1 - sN) = P2 K m0 + G2 (sN - sF), q = P2 K m0 + G2 (sN - sF) */
    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);
    if (!DET_UT) std::cout << "seed =" << sd << std::endl;

    const Eigen::Vector3d p = rand_unit(gen), m0 = rand_unit(gen);
    const Eigen::Matrix3d G1 = Slonczewski::conductance(p, 0.5 + distrib(gen), 0.5 + distrib(gen),
                                                        0.2 * distrib(gen));
    const Eigen::Matrix3d G2 = Slonczewski::conductance(m0, 0.5 + distrib(gen), 0.5 + distrib(gen),
                                                        0.2 * distrib(gen));
    const double P1 = distrib(gen), P2 = distrib(gen), K = distrib(gen);
    const Eigen::Vector3d sF1(distrib(gen), distrib(gen), distrib(gen));
    const Eigen::Vector3d sF(distrib(gen), distrib(gen), distrib(gen));

    const Eigen::Vector3d sN = (G1 + G2).lu().solve(K * (P1 * p - P2 * m0) + G1 * sF1 + G2 * sF);
    const Eigen::Vector3d q_ref = P2 * K * m0 + G2 * (sN - sF);
    // spin current conservation through the spacer
    const Eigen::Vector3d q1 = P1 * K * p + G1 * (sF1 - sN);
    BOOST_TEST( rel_err(q1, q_ref) < 1e-12 );

    const Slonczewski::BC bc = Slonczewski::effectiveBC(G1, P1, p, G2, P2, m0);
    const Eigen::Vector3d q = K * bc.Peff + bc.Gs * (sF1 - sF);
    BOOST_TEST( rel_err(q, q_ref) < 1e-12 );

    // Gs = G2 A^-1 G1 = G1 A^-1 G2
    const Eigen::Matrix3d Gs_sym = G1 * (G1 + G2).lu().solve(G2);
    BOOST_TEST( rel_err(bc.Gs, Gs_sym) < 1e-12 );
    }

BOOST_AUTO_TEST_CASE(effectiveBC_limits)
    {
    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);
    if (!DET_UT) std::cout << "seed =" << sd << std::endl;

    const Eigen::Vector3d p = rand_unit(gen), m0 = rand_unit(gen);
    const double P1 = distrib(gen), P2 = distrib(gen);
    const double G1par = 0.5 + distrib(gen), G1r = 0.5 + distrib(gen);
    const double G2par = 0.5 + distrib(gen), G2r = 0.5 + distrib(gen);

    /* resistive F1/N interface (G1 << G2): Peff = P1 p, Gs = G1, the flux condition with the
     polarization P1 p */
    const double eps = 1e-8;
    const Eigen::Matrix3d G1 = Slonczewski::conductance(p, eps * G1par, eps * G1r, 0);
    const Eigen::Matrix3d G2 = Slonczewski::conductance(m0, G2par, G2r, 0);
    Slonczewski::BC bc = Slonczewski::effectiveBC(G1, P1, p, G2, P2, m0);
    // absolute error: the polarizations are at most 1, P1 may be small
    BOOST_TEST( (bc.Peff - P1 * p).norm() < 1e-6 );
    BOOST_TEST( rel_err(bc.Gs, G1) < 1e-6 );

    /* resistive N/F2 interface (G2 << G1): Peff = P2 m0 */
    bc = Slonczewski::effectiveBC(Slonczewski::conductance(p, G1par, G1r, 0), P1, p,
                                  Slonczewski::conductance(m0, eps * G2par, eps * G2r, 0), P2, m0);
    BOOST_TEST( (bc.Peff - P2 * m0).norm() < 1e-6 );
    }

BOOST_AUTO_TEST_CASE(effectiveBC_collinear)
    {
    /* m0 = p: two longitudinal conductances in series, and the polarizations weighted by the
     conductances, Peff = (G2 P1 + G1 P2)/(G1 + G2) p; the transverse part does not act on p */
    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);
    if (!DET_UT) std::cout << "seed =" << sd << std::endl;

    const Eigen::Vector3d p = rand_unit(gen);
    const double P1 = distrib(gen), P2 = distrib(gen);
    const double G1par = 0.5 + distrib(gen), G2par = 0.5 + distrib(gen);
    const Eigen::Matrix3d G1 = Slonczewski::conductance(p, G1par, 0.5 + distrib(gen), distrib(gen));
    const Eigen::Matrix3d G2 = Slonczewski::conductance(p, G2par, 0.5 + distrib(gen), distrib(gen));

    const Slonczewski::BC bc = Slonczewski::effectiveBC(G1, P1, p, G2, P2, p);
    const double Peff = (G2par * P1 + G1par * P2) / (G1par + G2par);
    const Eigen::Vector3d Gs_p = bc.Gs * p;
    const double Gs = G1par * G2par / (G1par + G2par);
    BOOST_TEST( rel_err(bc.Peff, Eigen::Vector3d(Peff * p)) < 1e-13 );
    BOOST_TEST( rel_err(Gs_p, Eigen::Vector3d(Gs * p)) < 1e-13 );
    }

BOOST_AUTO_TEST_CASE(stack)
    {
    /* without N/F2, Peff = P1 p and Gs = G1; with a N/F2 interface much more conducting than F1/N,
     the condition tends to the one without N/F2 */
    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);
    if (!DET_UT) std::cout << "seed =" << sd << std::endl;

    Slonczewski::Stack st;
    st.p = rand_unit(gen);
    st.P1 = distrib(gen);
    st.G1par = 0.5 + distrib(gen);
    st.G1r = 0.5 + distrib(gen);
    st.G1i = 0.2 * distrib(gen);
    const Eigen::Vector3d m0 = rand_unit(gen);
    const Eigen::Matrix3d G1 = Slonczewski::conductance(st.p, st.G1par, st.G1r, st.G1i);

    const Slonczewski::BC bcB = st.bc(m0);
    BOOST_TEST( rel_err(bcB.Peff, Eigen::Vector3d(st.P1 * st.p)) < 1e-15 );
    BOOST_TEST( rel_err(bcB.Gs, G1) < 1e-15 );

    st.hasNF2 = true;
    st.P2 = distrib(gen);
    st.G2par = 1e8 * (0.5 + distrib(gen));
    st.G2r = 1e8 * (0.5 + distrib(gen));
    st.G2i = 0;
    const Slonczewski::BC bcA = st.bc(m0);
    BOOST_TEST( (bcA.Peff - bcB.Peff).norm() < 1e-6 );
    BOOST_TEST( rel_err(bcA.Gs, bcB.Gs) < 1e-6 );

    // and with N/F2, the condition is the one of effectiveBC
    const Slonczewski::BC bc = Slonczewski::effectiveBC(G1, st.P1, st.p,
            Slonczewski::conductance(m0, st.G2par, st.G2r, st.G2i), st.P2, m0);
    BOOST_TEST( rel_err(bcA.Peff, bc.Peff) < 1e-15 );
    BOOST_TEST( rel_err(bcA.Gs, bc.Gs) < 1e-15 );
    }

BOOST_AUTO_TEST_SUITE_END()
