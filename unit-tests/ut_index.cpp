#define BOOST_TEST_MODULE indexTest

#include <boost/test/unit_test.hpp>
#include <numeric>   // iota
#include "meshUtils.h"

BOOST_AUTO_TEST_SUITE(ut_index)

BOOST_AUTO_TEST_CASE(idx_suppress_copies)
    {
    std::vector<int> my_indices {1,7,666,3,3,3,4,5,-10};
    suppress_copies<int>(my_indices);
    BOOST_CHECK(my_indices.size() == (size_t)7);
    BOOST_CHECK(my_indices[0] == -10);
    BOOST_CHECK(my_indices[1] == 1);
    BOOST_CHECK(my_indices[2] == 3);
    BOOST_CHECK(my_indices[3] == 4);
    BOOST_CHECK(my_indices[4] == 5);
    BOOST_CHECK(my_indices[5] == 7);
    BOOST_CHECK(my_indices[6] == 666);
    }

BOOST_AUTO_TEST_CASE(deprecated_stuff)
    {
    const int NOD = 1000;
    std::vector<int> ld {1,3,5,42,123};
    std::vector<int> lvd;
    lvd.assign(ld.begin(),ld.end());
    lvd.resize(2*ld.size());
    // bind2nd() is deprecated, we want to check if the replacing code is doing the same.
    std::transform(ld.begin(),ld.end(),lvd.begin()+ld.size(), bind2nd(std::plus<int>(),NOD) );

    std::vector<int> lvdNew;
    lvdNew.resize(2*ld.size());
    for(unsigned int i=0;i<ld.size();i++)
        {
        lvdNew[i] = ld[i];
        lvdNew[i + ld.size()] = ld[i] + NOD;
        }
    std::cout << "lvd= ";
    for(unsigned int i=0;i<lvd.size();i++) {std::cout<< lvd[i] << '\t'; }
    std::cout << "\nlvdNew= ";
    for(unsigned int i=0;i<lvdNew.size();i++) {std::cout<< lvdNew[i] << '\t'; }
    std::cout << std::endl;
    BOOST_CHECK(lvd.size() == lvdNew.size());
    for(unsigned int i=0;i<lvd.size();i++)
        { BOOST_CHECK(lvd[i] == lvdNew[i]); }
    }

BOOST_AUTO_TEST_CASE(LLG_index_lvd_preparation)
    {
    /*
     *if diff size is initialized bigger than what is really needed, there are some zero's in the
     resulting diff if set_difference fifth argument is just diff.begin().
     But you do not know in advance diff size, so to avoid setting an arbitrary initial size and
     resizing after calling set_difference it is better not to memory initialize diff and call
     set_difference with the proper inserter
     * */
    const int NOD = 10;
    std::vector<int> v1{0, 0, 1, 2, 5, 5, 5, 9};
    suppress_copies<int>(v1);
    std::cout << "v1= ";
    for(unsigned int i=0;i<v1.size();i++) { std::cout << v1[i] << '\t'; }
    std::cout << std::endl;
    std::vector<int> v2(NOD);
    std::iota(v2.begin(),v2.end(),0);
    std::cout << "v2= ";
    for(unsigned int i=0;i<v2.size();i++) { std::cout << v2[i] << '\t'; }
    std::cout << std::endl;
    std::vector<int> diff;
    std::set_difference(v2.begin(),v2.end(),v1.begin(),v1.end(),std::inserter(diff,diff.begin()));
    std::cout << "v2 - v1 = diff= ";
    for(unsigned int i=0;i<diff.size();i++) { std::cout << diff[i] << '\t'; }
    std::cout << std::endl;
    BOOST_CHECK(diff.size() == (size_t)5);
    }

BOOST_AUTO_TEST_SUITE_END()

