#ifndef FMM_CORE_H
#define FMM_CORE_H

/** \file fmm_core.h
\brief fast multipole method (FMM) of ScalFMM: potentials at target points created by point
charges.

This interface depends neither on ScalFMM nor on Eigen: the ScalFMM headers, which are long to
compile, are only included by fmm_core.cpp.
*/

#include <memory>
#include <vector>

namespace scal_fmm
    {
/** width of the cubic box of the octree, centered on the origin: the particles must be in it */
const double boxWidth = 2.01;

/** \class FmmCore
octree, kernel and algorithm of the FMM for a fixed set of particles: targets, where the potential
is computed, and sources, carrying the charges
*/
class FmmCore
    {
public:
    /** builds the octree of height nbLevels with the particles of positions xyz (x, y, z of the
    nbTargets targets, then of the sources, in the box of width boxWidth centered on the origin),
    and the kernel and the algorithm of the FMM, which uses nbThreads threads */
    FmmCore(const std::vector<double> &xyz /**< [in] */, const int nbTargets /**< [in] */,
            const int nbLevels /**< [in] */, const int nbThreads /**< [in] */);

    /** destructor, defined in fmm_core.cpp where the implementation is complete */
    ~FmmCore();

    /** computes the potentials at the targets, $\sum_j q_j / |x - y_j|$, created by the charges q
    of the sources */
    void compute(const std::vector<double> &q /**< [in] */,
                 std::vector<double> &potentials /**< [out] */);

private:
    /** octree, kernel and algorithm, defined in fmm_core.cpp */
    struct Impl;

    /** implementation */
    std::unique_ptr<Impl> impl;
    };
    }  // namespace scal_fmm
#endif
