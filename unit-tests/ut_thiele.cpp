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
    BOOST_TEST(true);
    }

BOOST_AUTO_TEST_SUITE_END()

