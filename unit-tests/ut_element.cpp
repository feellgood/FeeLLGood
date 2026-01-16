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

BOOST_AUTO_TEST_CASE(Fac_assemble_vect, *boost::unit_test::tolerance(UT_TOL))
    {
    const int nbNod = 3;
    std::vector<Nodes::Node> node;
    dummyNodes<nbNod>(node);

    Facette::Fac f(node, nbNod, 0, {1, 2, 3});  // carefull with the index shift
    std::cout << "fac.assemble_vect test: test the norm of the output vector built\n";
    
    std::vector<double> L(2*nbNod,1.0);
    f.assemble_vect(L);// init val of f.Lp is zero
    BOOST_TEST(algebra::norm(L) == sqrt(6.0));
    f.Lp[2*Facette::N - 1] -= 1.0;
    f.assemble_vect(L);
    BOOST_TEST(algebra::norm(L) == sqrt(5.0));
    }

BOOST_AUTO_TEST_CASE(Fac_assemble_mat, *boost::unit_test::tolerance(UT_TOL))
    {
    const int nbNod = 3;
    std::vector<Nodes::Node> node;
    dummyNodes<nbNod>(node);

    Facette::Fac f(node, nbNod, 0, {1, 2, 3});  // carefull with the index shift
    std::cout << "fac.assemble_mat test: test the output matrix built\n";
    
    algebra::MatrixShape shape(2*nbNod);
    for(unsigned int i=0;i<2*nbNod;i++)
        for(unsigned int j=0;j<2*nbNod;j++)
            shape[i].insert(j);
    algebra::r_sparseMat Kr(shape);
    f.assemble_mat(Kr);// init val of f.Kp is zero
    
    for(unsigned int i=0;i<2*nbNod;i++)
        for(unsigned int j=0;j<2*nbNod;j++)
            {
            //std::cout <<"K(" << i <<"," << j << ")= " << Kr(i,j) << std::endl; 
            BOOST_CHECK( Kr(i,j) == 0.0); // must be stricly zero
            }
    
    f.Kp(2*Facette::N - 1,2*Facette::N - 1) += 1.0;
    f.assemble_mat(Kr);
    int i = Facette::N - 1;
    BOOST_CHECK(Kr(2*f.ind[i]+1, 2*f.ind[i]+1) == 1.0);
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

BOOST_AUTO_TEST_CASE(Tet_assemble_vect, *boost::unit_test::tolerance(UT_TOL))
    {
    const int nbNod = 4;
    std::vector<Nodes::Node> node;
    dummyNodes<nbNod>(node);

    // carefull with indices (starting from 1)
    Tetra::Tet t(node, 0, {1, 2, 3, 4});
    std::cout << "Tet.assemble_vect test: test the norm of the output vector built\n";
    
    std::vector<double> L(2*nbNod,1.0);
    t.assemble_vect(L);// init val of t.Lp is zero
    BOOST_TEST(algebra::norm(L) == sqrt(8.0));
    t.Lp[2*Tetra::N - 1] -= 1.0;
    t.assemble_vect(L);
    BOOST_TEST(algebra::norm(L) == sqrt(7.0));
    }

BOOST_AUTO_TEST_CASE(Tet_assemble_mat, *boost::unit_test::tolerance(UT_TOL))
    {
    const int nbNod = 4;
    std::vector<Nodes::Node> node;
    dummyNodes<nbNod>(node);

    Tetra::Tet t(node, 0, {1, 2, 3, 4});  // carefull with the index shift
    std::cout << "Tet.assemble_mat test: test the output matrix built\n";
    
    algebra::MatrixShape shape(2*nbNod);
    for(unsigned int i=0;i<2*nbNod;i++)
        for(unsigned int j=0;j<2*nbNod;j++)
            shape[i].insert(j);
    algebra::r_sparseMat Kr(shape);
    t.assemble_mat(Kr);// init val of t.Kp is zero
    
    for(unsigned int i=0;i<2*nbNod;i++)
        for(unsigned int j=0;j<2*nbNod;j++)
            {
            //std::cout <<"K(" << i <<"," << j << ")= " << Kr(i,j) << std::endl; 
            BOOST_CHECK( Kr(i,j) == 0.0); // must be stricly zero
            }
    
    t.Kp(2*Tetra::N - 1,2*Tetra::N - 1) += 1.0;
    t.assemble_mat(Kr);
    int i = Tetra::N - 1;
    BOOST_CHECK(Kr(2*t.ind[i]+1, 2*t.ind[i]+1) == 1.0);
    }

BOOST_AUTO_TEST_SUITE_END()
