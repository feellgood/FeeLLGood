#ifndef OCTREE_HEIGHT_H
#define OCTREE_HEIGHT_H

/** \file octree_height.h
\brief height of the octree of the fast multipole method, adapted to the particles.

The height is the deepest one (at least MinLevels, at most MaxLevels) whose non-empty leaves
contain on average at least MinParticlesPerLeaf particles. A deeper octree only adds cells holding
a few particles: the far field (M2L) then costs much more than the direct computation (P2P) it
replaces. This heuristic comes from st-feeLLGood, where it was calibrated with the rotation kernel
(P = 9) on meshes of 2e4 to 7.5e6 particles: the fastest octree always has ~300 particles per
non-empty leaf, one level deeper (~35 particles) is 4 to 5 times slower, one level less (~2000
particles) 10 to 25 % slower; the threshold 150 selects the fastest one.

The non-empty leaves are counted from the sorted Morton keys of the particles: this is
independent of ScalFMM, and costs a sort of the particles.
*/

#include <algorithm>
#include <cstdint>
#include <vector>

#include <eigen3/Eigen/Dense>

namespace scal_fmm
    {
const int MinLevels = 3;                /**< minimal height of the octree */
const int MaxLevels = 10;               /**< maximal height of the octree */
const double MinParticlesPerLeaf = 150; /**< minimal mean number of particles per non-empty leaf */

/** Morton key of point p of the cubic box of width boxWidth centered on the origin, at depth
MaxLevels-1 (2^(MaxLevels-1) cells per direction). The key of the cell containing p at depth
h-1 (the leaves of an octree of height h) is mortonKey(p) >> 3*(MaxLevels-h). */
inline uint64_t mortonKey(const Eigen::Ref<const Eigen::Vector3d> p /**< [in] */,
                          const double boxWidth /**< [in] */)
    {
    const int n = 1 << (MaxLevels - 1);
    uint64_t key = 0;
    for (int d = 0; d < 3; d++)
        {
        int c = static_cast<int>((p(d) + boxWidth / 2) / boxWidth * n);
        c = std::clamp(c, 0, n - 1);
        for (int b = 0; b < MaxLevels - 1; b++)
            { key |= uint64_t((c >> b) & 1) << (3 * b + d); }
        }
    return key;
    }

/** height of the octree for the particles of positions pos, in the cubic box of width boxWidth
centered on the origin: the deepest octree whose non-empty leaves contain on average at least
MinParticlesPerLeaf particles, clamped to [MinLevels, MaxLevels] */
inline int octreeHeight(const std::vector<Eigen::Vector3d> &pos /**< [in] */,
                        const double boxWidth /**< [in] */)
    {
    std::vector<uint64_t> keys(pos.size());
    std::transform(pos.begin(), pos.end(), keys.begin(),
                   [boxWidth](const Eigen::Vector3d &p) { return mortonKey(p, boxWidth); });
    std::sort(keys.begin(), keys.end());
    for (int h = MaxLevels; h > MinLevels; h--)
        {
        const int shift = 3 * (MaxLevels - h);
        size_t nbLeaves = 0;
        for (size_t i = 0; i < keys.size(); i++)
            {
            if (i == 0 || (keys[i] >> shift) != (keys[i - 1] >> shift))
                { nbLeaves++; }
            }
        if (double(keys.size()) >= MinParticlesPerLeaf * double(nbLeaves))
            { return h; }
        }
    return MinLevels;
    }
    }  // namespace scal_fmm
#endif
