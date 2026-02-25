#define BOOST_TEST_MODULE fac_chargesTest

#include <boost/test/unit_test.hpp>

#include <random>
#include "facette.h"
#include "ut_tools.h"
#include "ut_config.h"

/*-----------------------------------------------------

class facette charges is tested here

-------------------------------------------------------*/

BOOST_AUTO_TEST_SUITE(ut_fac_charges)

BOOST_AUTO_TEST_CASE(fac_charges, *boost::unit_test::tolerance(UT_TOL))
    {
    const int nbNod = 3;
    std::vector<Nodes::Node> node;
    dummyNodes<nbNod>(node);

    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);

    for (int i = 0; i < nbNod; i++)
        {
        node[i].d[Nodes::NEXT].u = rand_vec3d(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
        }

    Facette::Fac f(node, nbNod, 0, {1, 2, 3});  // carefull with the index shift
    double dMs = distrib(gen);
    auto getter = Nodes::get_u<Nodes::NEXT>;

    // code to test
    Eigen::Matrix<double,Facette::NPI,1> result = f.charges(dMs, getter);

    //ref code begin
    Eigen::Matrix<double,Facette::NPI,1> resultRef;
    resultRef.setZero();
    
    Eigen::Matrix<double,Nodes::DIM,Facette::N> vec_nod;
    for(int i=0;i<Facette::N;i++)
        { vec_nod.col(i) << node[i].d[Nodes::NEXT].u; }
    Eigen::Matrix<double,Nodes::DIM,Facette::NPI> _u = vec_nod * Facette::eigen_a;
    resultRef = dMs*f.weight.cwiseProduct( _u.transpose()*f.n );

    //ref code end
    for (int j=0;j < Facette::NPI;j++)
        BOOST_TEST(resultRef(j) == result(j));
    }

BOOST_AUTO_TEST_SUITE_END()
