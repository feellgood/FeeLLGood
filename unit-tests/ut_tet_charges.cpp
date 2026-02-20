#define BOOST_TEST_MODULE tet_chargesTest

#include <boost/test/unit_test.hpp>

#include <random>

#include "ut_tools.h"
#include "ut_config.h"

/*-----------------------------------------------------

class tetra charges is tested here

-------------------------------------------------------*/

BOOST_AUTO_TEST_SUITE(ut_tet_charges)

BOOST_AUTO_TEST_CASE(tet_charges, *boost::unit_test::tolerance(UT_TOL))
    {
    const int nbNod = 4;
    std::vector<Nodes::Node> node;
    dummyNodes<nbNod>(node);

    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);

    for (int i = 0; i < nbNod; i++)
        {
        node[i].d[Nodes::NEXT].u = rand_vec3d(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
        }

    // carefull with indices (starting from 1)
    Tetra::Tet t(node, 0, {1, 2, 3, 4});
    Tetra::prm p;
    p.Ms = distrib(gen);
    auto getter = Nodes::get_u<Nodes::NEXT>;

    // code to test
    Eigen::Matrix<double,Tetra::NPI,1> result = t.charges(p.Ms, getter);

    // ref code begin
    Eigen::Matrix<double,Nodes::DIM,Tetra::N> vec_nod;
    for (int i = 0; i < Tetra::N; i++) vec_nod.col(i) = getter(node[i]);//t.getNode(i) is private

    Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> dudx = vec_nod * (t.da.col(Nodes::IDX_X)).replicate(1,Tetra::NPI);
    Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> dudy = vec_nod * (t.da.col(Nodes::IDX_Y)).replicate(1,Tetra::NPI);
    Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> dudz = vec_nod * (t.da.col(Nodes::IDX_Z)).replicate(1,Tetra::NPI);

    Eigen::Matrix<double,Tetra::NPI,1> resultRef;
    for (int j = 0; j < Tetra::NPI; j++)
        { resultRef(j) = dudx(0,j) + dudy(1,j) + dudz(2,j); }
    resultRef = t.weight.cwiseProduct(resultRef);
    resultRef *= -p.Ms;
    // ref code end
    for (int j = 0; j < Tetra::NPI; j++)
        BOOST_TEST(resultRef(j) == result(j));
    }

BOOST_AUTO_TEST_SUITE_END()
