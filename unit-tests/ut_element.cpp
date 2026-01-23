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

BOOST_AUTO_TEST_CASE(Fac_full_constructor)
    {
    std::vector<Nodes::Node> node;
    std::cout << "4 param constructor" << std::endl;
    
    Facette::Fac f(node, 0, idxPrmToTest, {0, 0, 0});
    /*
    The constructor facette does not modify all indices to fit to zero based index convention
    it's done later
    */
    const int idx = 0;
    BOOST_CHECK(f.ind[0] == idx);
    BOOST_CHECK(f.ind[1] == idx);
    BOOST_CHECK(f.ind[2] == idx);
    BOOST_CHECK((f.getN() == Facette::N) && (f.getNPI() == Facette::NPI));
    BOOST_CHECK(f.dMs == 0);
    BOOST_CHECK(f.idxPrm == idxPrmToTest);
    BOOST_CHECK(f.Lp.norm() == 0);
    }


BOOST_AUTO_TEST_CASE(Tet_constructor)
    {
    std::cout << "constructor test with empty node vector\n";
    std::vector<Nodes::Node> node(0);

    Tetra::Tet tet(node, idxPrmToTest, {0, 0, 0, 0});
    /*
    The constructor tetra do modify all indices to fit to zero based index convention
    */
    const int idx = -1;
    BOOST_CHECK(tet.ind[0] == idx);
    BOOST_CHECK(tet.ind[1] == idx);
    BOOST_CHECK(tet.ind[2] == idx);
    BOOST_CHECK(tet.ind[3] == idx);
    BOOST_CHECK((tet.getN() == Tetra::N) && (tet.getNPI() == Tetra::NPI));
    BOOST_CHECK(tet.idxPrm == idxPrmToTest);
    BOOST_CHECK(tet.Lp.norm() == 0);
    }

BOOST_AUTO_TEST_CASE(Tet_constructor_with_wrong_init_list)
    {
    std::cout << "constructor test with empty node vector\n";
    std::vector<Nodes::Node> node(0);
    const int extra(0);
    Tetra::Tet tet(node, idxPrmToTest, {0, 0, 0, 0, extra});
    BOOST_CHECK(tet.ind.empty());
    }

BOOST_AUTO_TEST_SUITE_END()
