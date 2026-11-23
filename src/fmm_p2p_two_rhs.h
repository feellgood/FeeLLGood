#ifndef FMM_P2P_TWO_RHS_H
#define FMM_P2P_TWO_RHS_H

/** \file fmm_p2p_two_rhs.h
\brief direct interactions (P2P) of the ScalFMM interpolation kernels with two right-hand sides.

The potentials of u and v are computed in a single FMM pass, with two right-hand sides (NVALS = 2).
ScalFMM does not provide a suitable P2P for this case: depending on the version, it either has none
or a scalar one which also computes the forces. This header specializes the class
DirectInteractionComputer of ScalFMM for a scalar kernel (1/r) and two right-hand sides: each
distance is computed once for both potentials, with the same SIMD type as the P2P of ScalFMM. It
must be included after the ScalFMM header that declares DirectInteractionComputer, and before the
instantiation of the kernel.
*/

#include "Kernels/Interpolation/FInterpP2PKernels.hpp"
#include "Kernels/P2P/FP2P.hpp"
#include "Utils/FMath.hpp"

namespace scal_fmm
    {
/** potentials at the targets of the sources of the neighbors, for two right-hand sides;
ComputeClass is a SIMD type holding NbDouble doubles, the source arrays of ScalFMM being padded
accordingly */
template<class ComputeClass, int NbDouble, class ContainerClass, class MatrixKernelClass>
void fullRemoteTwoRhs(ContainerClass *const FRestrict targets,
                      const ContainerClass *const neighbors[], const int nbNeighbors,
                      const MatrixKernelClass *const matrixKernel)
    {
    const FSize nbTargets = targets->getNbParticles();
    const double *const tx = targets->getPositions()[0];
    const double *const ty = targets->getPositions()[1];
    const double *const tz = targets->getPositions()[2];
    double *const pot0 = targets->getPotentials(0);
    double *const pot1 = targets->getPotentials(1);

    for (int n = 0; n < nbNeighbors; n++)
        {
        if (!neighbors[n])
            { continue; }
        const FSize nbBlocks = (neighbors[n]->getNbParticles() + NbDouble - 1) / NbDouble;
        const ComputeClass *const q0 = (const ComputeClass *)neighbors[n]->getPhysicalValues(0);
        const ComputeClass *const q1 = (const ComputeClass *)neighbors[n]->getPhysicalValues(1);
        const ComputeClass *const sx = (const ComputeClass *)neighbors[n]->getPositions()[0];
        const ComputeClass *const sy = (const ComputeClass *)neighbors[n]->getPositions()[1];
        const ComputeClass *const sz = (const ComputeClass *)neighbors[n]->getPositions()[2];
        for (FSize t = 0; t < nbTargets; t++)
            {
            const ComputeClass x = FMath::ConvertTo<ComputeClass, const double *>(&tx[t]);
            const ComputeClass y = FMath::ConvertTo<ComputeClass, const double *>(&ty[t]);
            const ComputeClass z = FMath::ConvertTo<ComputeClass, const double *>(&tz[t]);
            ComputeClass sum0 = FMath::Zero<ComputeClass>();
            ComputeClass sum1 = FMath::Zero<ComputeClass>();
            for (FSize s = 0; s < nbBlocks; s++)
                {
                const ComputeClass k = matrixKernel->evaluate(x, y, z, sx[s], sy[s], sz[s]);
                sum0 += k * q0[s];
                sum1 += k * q1[s];
                }
            pot0[t] += FMath::ConvertTo<double, ComputeClass>(sum0);
            pot1[t] += FMath::ConvertTo<double, ComputeClass>(sum1);
            }
        }
    }

#if defined(SCALFMM_USE_AVX)
typedef __m256d P2PComputeClass; /**< SIMD type of the P2P, as in ScalFMM */
const int P2PNbDouble = 4;       /**< number of doubles in P2PComputeClass */
#elif defined(SCALFMM_USE_AVX2)
typedef __m512d P2PComputeClass;
const int P2PNbDouble = 8;
#elif defined(SCALFMM_USE_SSE)
typedef __m128d P2PComputeClass;
const int P2PNbDouble = 2;
#else
typedef double P2PComputeClass;
const int P2PNbDouble = 1;
#endif
    }  // namespace scal_fmm

/** P2P of a scalar kernel with two right-hand sides. The target-source FMM algorithm only calls
P2PRemote; P2PInner and P2P are given for completeness. */
template<>
struct DirectInteractionComputer<double, 1, 2>
    {
    /** potentials at the targets of the sources of the neighbors */
    template<typename ContainerClass, typename MatrixKernelClass>
    static void P2PRemote(ContainerClass *const FRestrict targets,
                          const ContainerClass *const neighbors[], const int nbNeighbors,
                          const MatrixKernelClass *const matrixKernel)
        {
        scal_fmm::fullRemoteTwoRhs<scal_fmm::P2PComputeClass, scal_fmm::P2PNbDouble>(
                targets, neighbors, nbNeighbors, matrixKernel);
        }

    /** interactions between the particles of a single container */
    template<typename ContainerClass, typename MatrixKernelClass>
    static void P2PInner(ContainerClass *const FRestrict particles,
                         const MatrixKernelClass *const matrixKernel)
        {
        const FSize nb = particles->getNbParticles();
        const double *const x = particles->getPositions()[0];
        const double *const y = particles->getPositions()[1];
        const double *const z = particles->getPositions()[2];
        for (int v = 0; v < 2; v++)
            {
            const double *const q = particles->getPhysicalValues(v);
            double *const pot = particles->getPotentials(v);
            for (FSize i = 0; i < nb; i++)
                {
                for (FSize j = 0; j < nb; j++)
                    {
                    if (j != i)
                        {
                        pot[i] += matrixKernel->evaluate(x[i], y[i], z[i], x[j], y[j], z[j])
                                  * q[j];
                        }
                    }
                }
            }
        }

    /** mutual interactions between a container and its neighbors */
    template<typename ContainerClass, typename MatrixKernelClass>
    static void P2P(ContainerClass *const FRestrict targets, ContainerClass *const neighbors[],
                    const int nbNeighbors, const MatrixKernelClass *const matrixKernel)
        {
        P2PRemote(targets, neighbors, nbNeighbors, matrixKernel);
        const ContainerClass *const self[1] = {targets};
        for (int n = 0; n < nbNeighbors; n++)
            {
            if (neighbors[n])
                { P2PRemote(neighbors[n], self, 1, matrixKernel); }
            }
        }
    };

#endif
