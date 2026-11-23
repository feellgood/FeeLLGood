#ifndef FMM_CORE_H
#define FMM_CORE_H

/** \file fmm_core.h
\brief fast multipole method (FMM) of ScalFMM: potentials at target points created by point
charges.

This interface depends neither on ScalFMM nor on Eigen: the ScalFMM headers are only included by
fmm_core.cpp. They are long to compile, and the BLAS and LAPACK declarations of the Chebyshev kernel
conflict with those of the MKL headers, included by Eigen when feeLLGood is built with MKL.
*/

#include <array>
#include <memory>
#include <vector>

#include "config.h"

namespace scal_fmm
    {
/** width of the cubic box of the octree, centered on the origin: the particles must be in it */
const double boxWidth = 2.01;

/** number of right-hand sides of the FMM: the potentials of u and v are computed in a single pass,
the potential of v being needed only if FIRST_ORDER is OFF */
const int NVALS = FIRST_ORDER ? 1 : 2;

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

    /** computes in a single pass, for each right-hand side v, the potentials at the targets,
    $\sum_j q_j / |x - y_j|$, created by the charges q[v] of the sources */
    void compute(const std::array<std::vector<double>, NVALS> &q /**< [in] */,
                 std::array<std::vector<double>, NVALS> &potentials /**< [out] */);

private:
    /** octree, kernel and algorithm, defined in fmm_core.cpp */
    struct Impl;

    /** implementation */
    std::unique_ptr<Impl> impl;
    };
    }  // namespace scal_fmm
#endif
