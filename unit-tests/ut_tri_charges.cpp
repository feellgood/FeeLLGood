#define BOOST_TEST_MODULE tri_chargesTest

#include <boost/test/unit_test.hpp>

#include <random>
#include "triangle.h"
#include "ut_tools.h"
#include "ut_config.h"

/*-----------------------------------------------------

class triangle charges is tested here

-------------------------------------------------------*/

BOOST_AUTO_TEST_SUITE(ut_tri_charges)

BOOST_AUTO_TEST_CASE(tri_charges, *boost::unit_test::tolerance(UT_TOL))
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

    Triangle::Tri f(node, nbNod, 0, {1, 2, 3});  // carefull with the index shift
    double dMs = distrib(gen);
    auto getter = Nodes::get_u<Nodes::NEXT>;

    // code to test
    Eigen::Matrix<double,Triangle::NPI,1> result = f.charges(dMs, getter);

    //ref code begin
    Eigen::Matrix<double,Triangle::NPI,1> resultRef;
    resultRef.setZero();
    
    Eigen::Matrix<double,Nodes::DIM,Triangle::N> vec_nod;
    for(int i=0;i<Triangle::N;i++)
        { vec_nod.col(i) << node[i].d[Nodes::NEXT].u; }
    Eigen::Matrix<double,Nodes::DIM,Triangle::NPI> _u = vec_nod * Triangle::eigen_a;
    resultRef = dMs*f.weight.cwiseProduct( _u.transpose()*f.n );

    //ref code end
    for (int j=0;j < Triangle::NPI;j++)
        BOOST_TEST(resultRef(j) == result(j));
    }

BOOST_AUTO_TEST_SUITE_END()
