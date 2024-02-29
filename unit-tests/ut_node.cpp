#define BOOST_TEST_MODULE nodeTest

#include <boost/test/unit_test.hpp>

#include <cstdlib>
#include <iostream>
#include <random>

#include "node.h"
#include "ut_config.h"  // for tolerance UT_TOL macro

BOOST_AUTO_TEST_SUITE(ut_node)

/*-----------------------------------------------------*/
/* zero lvl tests : direct elementary member functions */
/*-----------------------------------------------------*/

BOOST_AUTO_TEST_CASE(node_get_p_lvl0)
    {
    Nodes::Node n;
    Eigen::Vector3d pPos(1.0, 0.0, 0.0);

    n.p = pPos;

    BOOST_CHECK(n.p.x() == 1.0);
    BOOST_CHECK(n.p.y() == 0.0);
    BOOST_CHECK(n.p.z() == 0.0);
    }

/*---------------------------------------*/
/* first lvl tests : nested calculus,... */
/*---------------------------------------*/

BOOST_AUTO_TEST_CASE(node_get_p_lvl1, *boost::unit_test::tolerance(UT_TOL))
    {
    Nodes::Node n;
    Eigen::Vector3d pPos(1.0, 3.0, 5.0);

    n.p = pPos;
    n.p.normalize();

    BOOST_TEST(n.p.norm() == 1.0);
    }

/*---------------------------------------*/
/* second lvl tests : pure mathematics   */
/*---------------------------------------*/

// Build a unit vector from the cylindrical coordinates (theta, z).
// If theta and z are uniformly distributed in [-pi, pi] and [-1, 1]
// respectively, the resulting vector is isotropically distributed.
static Eigen::Vector3d unit_vector(double theta, double z)
    {
    double r = sqrt(1 - z * z);
    return Eigen::Vector3d(r * cos(theta), r * sin(theta), z);
    }

BOOST_AUTO_TEST_CASE(node_setBasis_eigen_formula, *boost::unit_test::tolerance(10.0 * UT_TOL))
    {
    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(-1.0, 1.0);

    Nodes::Node n;
    n.d[0].u = unit_vector(M_PI * distrib(gen), distrib(gen));
    n.setBasis(M_PI * distrib(gen));

    /* code to test */
    Eigen::Index minIdx;
    n.d[0].u.cwiseAbs().minCoeff(&minIdx);
    n.ep.setUnit(minIdx);
    /* end code to test */

    /* reference code */
    Eigen::Vector3d ep_ref;
    double abs_x = fabs(n.d[0].u.x()), abs_y = fabs(n.d[0].u.y()), abs_z = fabs(n.d[0].u.z());
        if (abs_x < abs_y)
            {
            if (abs_x < abs_z)
                { ep_ref = Eigen::Vector3d::UnitX(); }
            else
                { ep_ref = Eigen::Vector3d::UnitZ(); }
            }
        else
            {
            if (abs_y < abs_z)
                { ep_ref = Eigen::Vector3d::UnitY(); }
            else
                { ep_ref = Eigen::Vector3d::UnitZ(); }
            }
    /* end reference code */

    if (!DET_UT) std::cout << "seed =" << sd << std::endl;
    std::cout << "u= " << n.d[0].u << std::endl;
    std::cout << "ep_ref= " << ep_ref << std::endl;
    std::cout << "ep= " << n.ep << std::endl;
    Eigen::Vector3d result = ep_ref - n.ep;
    BOOST_TEST( result.norm() == 0.0 );
    }


BOOST_AUTO_TEST_CASE(node_e_p, *boost::unit_test::tolerance(10.0 * UT_TOL))
    {
    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(-1.0, 1.0);

    Nodes::Node n;

    // test the orthonormality of the basis (u0, ep, eq)
    n.d[0].u = unit_vector(M_PI * distrib(gen), distrib(gen));
    n.setBasis(M_PI * distrib(gen));

    if (!DET_UT) std::cout << "seed =" << sd << std::endl;
    BOOST_TEST( n.d[0].u.norm() == 1.0 );
    BOOST_TEST( n.ep.norm() == 1.0 );
    BOOST_TEST( n.eq.norm() == 1.0 );
    BOOST_TEST( n.d[0].u.dot(n.ep) == 0.0 );
    BOOST_TEST( n.ep.dot(n.eq) == 0.0 );
    BOOST_TEST( n.eq.dot(n.d[0].u) == 0.0 );
    }

BOOST_AUTO_TEST_CASE(node_evol, *boost::unit_test::tolerance(1e3 * UT_TOL))
    {
    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(-1.0, 1.0);

    Nodes::Node n;

    double vp = distrib(gen);
    double vq = distrib(gen);
    double dt = distrib(gen) + 1.0;

    n.d[0].u = unit_vector(M_PI * distrib(gen), distrib(gen));
    n.setBasis(M_PI * distrib(gen));

    Eigen::Vector3d v = vp * n.ep + vq * n.eq;

    n.make_evol(vp, vq, dt);

    if (!DET_UT) std::cout << "seed =" << sd << std::endl;
    std::cout << "v = " << v << std::endl;
    std::cout << "simplify[v] = " << n.get_v(Nodes::NEXT) << std::endl;

    BOOST_TEST( (n.get_v(Nodes::NEXT) - v).norm() == 0.0);
    }

BOOST_AUTO_TEST_SUITE_END()
