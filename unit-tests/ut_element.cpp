#define BOOST_TEST_MODULE tetraTest

#include <boost/test/unit_test.hpp>

#include <random>

#include "tetra.h"
#include "ut_config.h"
#include "ut_tools.h"

BOOST_AUTO_TEST_SUITE(ut_element)

/*-----------------------------------------------------*/
/* zero lvl tests : direct elementary member functions */
/*-----------------------------------------------------*/

const int idxPrmToTest = 42;

BOOST_AUTO_TEST_CASE(Tri_full_constructor)
    {
    const int nbNod(3);
    std::vector<Nodes::Node> node(nbNod);
    std::cout << "Tri constructor test\n" << std::endl;
    
    Triangle::Tri f(node, idxPrmToTest, {1, 2, 3});
    /*
    The constructor triangle do modify all indices to fit to zero based index convention
    */
    BOOST_CHECK(f.ind[0] == 0);
    BOOST_CHECK(f.ind[1] == 1);
    BOOST_CHECK(f.ind[2] == 2);
    BOOST_CHECK((f.getN() == Triangle::N) && (f.getNPI() == Triangle::NPI));
    BOOST_CHECK(f.dMs == 0);
    BOOST_CHECK(f.idxPrm == idxPrmToTest);
    BOOST_CHECK(f.Lp.norm() == 0);
    }


BOOST_AUTO_TEST_CASE(Tet_constructor)
    {
    std::cout << "Tet constructor test\n";
    const int nbNod(4);
    std::vector<Nodes::Node> node;
    dummyNodes<nbNod>(node);

    Tetra::Tet tet(node, idxPrmToTest, {1, 2, 3, 4});
    /*
    The constructor tetra do modify all indices to fit to zero based index convention
    */
    BOOST_CHECK(tet.ind[0] == 0);
    BOOST_CHECK(tet.ind[1] == 1);
    BOOST_CHECK(tet.ind[2] == 2);
    BOOST_CHECK(tet.ind[3] == 3);
    BOOST_CHECK((tet.getN() == Tetra::N) && (tet.getNPI() == Tetra::NPI));
    BOOST_CHECK(tet.idxPrm == idxPrmToTest);
    BOOST_CHECK(tet.Lp.norm() == 0);
    }

BOOST_AUTO_TEST_CASE(Tet_existNodes)
    {
    std::cout << "test Tetra::Tet::existNodes()\n";
    const int nbNod = 4;
    std::vector<Nodes::Node> node;
    dummyNodes<nbNod>(node);

    Tetra::Tet tet_good(node, idxPrmToTest, {1,2,3,4});
    Tetra::Tet tet_bad1(node, idxPrmToTest, {0,1,2,3});
    Tetra::Tet tet_bad2(node, idxPrmToTest, {2,3,4,5});
    BOOST_CHECK(tet_good.existNodes());
    BOOST_CHECK(!tet_bad1.existNodes());
    BOOST_CHECK(!tet_bad2.existNodes());
    }

BOOST_AUTO_TEST_SUITE_END()
