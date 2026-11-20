#define BOOST_TEST_MODULE octreeHeightTest

#include <boost/test/unit_test.hpp>

#include "octree_height.h"

using namespace scal_fmm;

const double boxWidth = 2.01;  // as in fmm_demag.h

/** n^3 points at the centers of the cells of a regular grid filling the cube [a, a + w]^3 */
std::vector<Eigen::Vector3d> grid(const int n, const double a, const double w)
    {
    std::vector<Eigen::Vector3d> pos;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                { pos.push_back(Eigen::Vector3d(a + w * (i + 0.5) / n, a + w * (j + 0.5) / n,
                                                a + w * (k + 0.5) / n)); }
    return pos;
    }

BOOST_AUTO_TEST_SUITE(ut_octree_height)

/* the key of a cell at depth h-1 is the key at depth MaxLevels-1 shifted by 3*(MaxLevels-h) bits */
BOOST_AUTO_TEST_CASE(morton_key_parent)
    {
    const int shift = 3 * (MaxLevels - 2);  // depth 1: 8 cells, key = x + 2 y + 4 z
    BOOST_CHECK_EQUAL(mortonKey(Eigen::Vector3d(-0.5, -0.5, -0.5), boxWidth) >> shift, 0u);
    BOOST_CHECK_EQUAL(mortonKey(Eigen::Vector3d(0.5, -0.5, -0.5), boxWidth) >> shift, 1u);
    BOOST_CHECK_EQUAL(mortonKey(Eigen::Vector3d(-0.5, 0.5, -0.5), boxWidth) >> shift, 2u);
    BOOST_CHECK_EQUAL(mortonKey(Eigen::Vector3d(-0.5, -0.5, 0.5), boxWidth) >> shift, 4u);
    BOOST_CHECK_EQUAL(mortonKey(Eigen::Vector3d(0.5, 0.5, 0.5), boxWidth) >> shift, 7u);
    // points outside the box are put in the boundary cells
    BOOST_CHECK_EQUAL(mortonKey(Eigen::Vector3d(-2, -2, -2), boxWidth), 0u);
    }

/* few particles: the smallest octree */
BOOST_AUTO_TEST_CASE(min_levels)
    {
    BOOST_CHECK_EQUAL(octreeHeight(grid(4, -1, 2), boxWidth), MinLevels);
    }

/* 64^3 particles filling the box: 8^3 leaves at height 4 (512 particles per leaf), 16^3 leaves
at height 5 (64 particles per leaf, less than MinParticlesPerLeaf) */
BOOST_AUTO_TEST_CASE(uniform_cube)
    {
    BOOST_CHECK_EQUAL(octreeHeight(grid(64, -boxWidth / 2, boxWidth), boxWidth), 4);
    }

/* the same particles in a cube 16 times smaller (one cell of depth 4): only the non-empty leaves
count, the octree is 4 levels deeper */
BOOST_AUTO_TEST_CASE(empty_leaves_ignored)
    {
    BOOST_CHECK_EQUAL(octreeHeight(grid(64, 0, boxWidth / 16), boxWidth), 8);
    }

/* all particles at the same point: the deepest octree */
BOOST_AUTO_TEST_CASE(max_levels)
    {
    std::vector<Eigen::Vector3d> pos(1000, Eigen::Vector3d(0.1, 0.2, 0.3));
    BOOST_CHECK_EQUAL(octreeHeight(pos, boxWidth), MaxLevels);
    }

BOOST_AUTO_TEST_SUITE_END()
