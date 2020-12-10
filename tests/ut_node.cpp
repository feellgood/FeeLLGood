#define BOOST_TEST_MODULE nodeTest

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(ut_node)

/*---------------------------------------*/
/* minus one test: check boost is fine   */
/*---------------------------------------*/

BOOST_AUTO_TEST_CASE(Stupid)
{
float x = 1.0;
BOOST_CHECK(x != 0.0f);
}

BOOST_AUTO_TEST_SUITE_END()
