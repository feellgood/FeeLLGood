#define BOOST_TEST_MODULE tet_chargesTest

#include <boost/test/unit_test.hpp>

#include <random>

#include "ut_tools.h"
#include "ut_config.h"

/*-----------------------------------------------------

class tetra charges is tested here

-------------------------------------------------------*/

BOOST_AUTO_TEST_SUITE(ut_tet_charges)

BOOST_AUTO_TEST_CASE(tet_exchange_lumping, *boost::unit_test::tolerance(UT_TOL))
    {
    int nsrc(0);
    std::vector<double> srcDen;
    srcDen.resize( Tetra::NPI );
    
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
    Tetra::prm p;
    p.J = 1.0;
    auto getter = Nodes::get_u<Nodes::NEXT>;
    t.charges(p, getter, srcDen, nsrc);
    BOOST_CHECK(nsrc == Tetra::NPI);
    BOOST_CHECK(srcDen.size() == nsrc);
    }

BOOST_AUTO_TEST_SUITE_END()
