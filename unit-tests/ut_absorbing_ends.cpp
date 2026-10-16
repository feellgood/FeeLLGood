#define BOOST_TEST_MODULE absorbingEndsTest

#include <boost/test/unit_test.hpp>

#include <cmath>

#include "absorbing_ends.h"
#include "ut_config.h"  // for tolerance UT_TOL macro

BOOST_AUTO_TEST_SUITE(ut_absorbing_ends)

/** zones of length 50 nm at both ends of a wire along z, from -200 nm to 200 nm */
AbsorbingEnds wire()
    {
    AbsorbingEnds ab;
    ab.enable = true;
    ab.axis = Nodes::IDX_Z;
    ab.length = 50e-9;
    ab.alpha_max = 1.0;
    ab.lo = -200e-9;
    ab.hi = 200e-9;
    return ab;
    }

BOOST_AUTO_TEST_CASE(disabled)
    {
    AbsorbingEnds ab = wire();
    ab.enable = false;
    BOOST_TEST(ab.alpha(0.05, 199e-9) == 0.05);
    }

BOOST_AUTO_TEST_CASE(outside_the_zones, *boost::unit_test::tolerance(UT_TOL))
    {
    const AbsorbingEnds ab = wire();
    BOOST_TEST(ab.alpha(0.05, 0.0) == 0.05);
    BOOST_TEST(ab.alpha(0.05, 149e-9) == 0.05);
    BOOST_TEST(ab.alpha(0.05, -149e-9) == 0.05);
    BOOST_TEST(ab.alpha(0.05, 150e-9) == 0.05);  // inner boundary of the zone: s = 0
    }

BOOST_AUTO_TEST_CASE(profile, *boost::unit_test::tolerance(10.0*UT_TOL))
    {
    const AbsorbingEnds ab = wire();
    const double alpha0 = 0.05;
    // end of the regions: alpha_max
    BOOST_TEST(ab.alpha(alpha0, 200e-9) == 1.0);
    BOOST_TEST(ab.alpha(alpha0, -200e-9) == 1.0);
    // beyond the end (a Gauss point cannot be there, but the profile is clamped)
    BOOST_TEST(ab.alpha(alpha0, 210e-9) == 1.0);
    // middle of the zone: sin^2(pi/4) = 1/2
    BOOST_TEST(ab.alpha(alpha0, 175e-9) == alpha0 + 0.5 * (1.0 - alpha0));
    // symmetric zones
    for (double z = 150e-9; z <= 200e-9; z += 5e-9)
        { BOOST_TEST(ab.alpha(alpha0, z) == ab.alpha(alpha0, -z)); }
    }

/* finite differences: their truncation error is far above UT_TOL, hence explicit tolerances */
BOOST_AUTO_TEST_CASE(zero_slope_at_both_boundaries)
    {
    const AbsorbingEnds ab = wire();
    const double alpha0 = 0.05;
    const double dz = 1e-12;
    // slopes relative to (alpha_max - alpha0)/length
    const double scale = (ab.alpha_max - alpha0) / ab.length;
    const double slope_inner = (ab.alpha(alpha0, 150e-9 + dz) - ab.alpha(alpha0, 150e-9)) / dz;
    const double slope_end = (ab.alpha(alpha0, 200e-9) - ab.alpha(alpha0, 200e-9 - dz)) / dz;
    const double slope_middle = (ab.alpha(alpha0, 175e-9 + dz) - ab.alpha(alpha0, 175e-9)) / dz;
    BOOST_TEST(std::abs(slope_inner / scale) < 1e-3);
    BOOST_TEST(std::abs(slope_end / scale) < 1e-3);
    // in the middle, the slope is pi/2 (alpha_max - alpha0)/length
    BOOST_TEST(slope_middle / scale == M_PI / 2, boost::test_tools::tolerance(1e-4));
    }

BOOST_AUTO_TEST_SUITE_END()
