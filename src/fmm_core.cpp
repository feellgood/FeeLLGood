#include "fmm_core.h"

#include <omp.h>

#include <algorithm>
#include <iostream>
#include <sstream>

#include "Components/FParticleType.hpp"
#include "Components/FTypedLeaf.hpp"
#include "Containers/FOctree.hpp"
#include "Core/FFmmAlgorithmThreadTsm.hpp"
#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"
#include "Kernels/Interpolation/FInterpMatrixKernel.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

/** The BLAS library, used by the Chebyshev kernel, must run on a single thread: the FMM is already
parallelized with OpenMP, and the BLAS threads would compete with it (with MKL, whose OpenMP runtime
differs from the one of GCC, the run does not progress anymore). Weak symbols: null if the BLAS
library is not OpenBLAS, respectively MKL. */
extern "C" void openblas_set_num_threads(int) __attribute__((weak));
extern "C" void MKL_Set_Num_Threads(int) __attribute__((weak));

namespace scal_fmm
    {
namespace
    {
const int ORDER = 6;          /**< number of Chebyshev interpolation points per direction */
const int SizeSubLevels = 3;  /**< maximal size of the sub levels */

typedef double FReal; /**< all computations are made in double precision */
typedef FTypedChebCell<FReal, ORDER> CellClass;
typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
typedef FTypedLeaf<FReal, ContainerClass> LeafClass;
typedef FOctree<FReal, CellClass, ContainerClass, LeafClass> OctreeClass;
typedef FInterpMatrixKernelR<FReal> MatrixKernelClass; /**< interaction kernel 1/r */
typedef FChebSymKernel<FReal, CellClass, ContainerClass, MatrixKernelClass, ORDER> KernelClass;
typedef FFmmAlgorithmThreadTsm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>
        FmmClass;

const FPoint<FReal> boxCenter(0., 0., 0.); /**< center of the box of the octree */
    }  // namespace

/** octree, kernel and algorithm of the FMM. The particles are indexed in the order of xyz: targets
from 0 to nbTargets - 1, then sources. */
struct FmmCore::Impl
    {
    const int nbTargets;                    /**< number of targets */
    OctreeClass tree;                       /**< octree */
    const MatrixKernelClass matrixKernel;   /**< interaction kernel, referenced by the kernel */
    std::unique_ptr<KernelClass> kernels;   /**< kernel */
    std::unique_ptr<FmmClass> algo; /**< algorithm, built once the particles are inserted */

    Impl(const int nbT, const int nbLevels)
        : nbTargets(nbT),
          tree(nbLevels, std::min(SizeSubLevels, nbLevels - 1), boxWidth, boxCenter)
        {
        // The Chebyshev kernel prints the time of its precomputations on std::cout: discard it.
        std::ostringstream discarded;
        struct CoutRedirect
            {
            std::streambuf *saved;
            explicit CoutRedirect(std::streambuf *buf) : saved(std::cout.rdbuf(buf)) {}
            ~CoutRedirect() { std::cout.rdbuf(saved); }
            } redirect(discarded.rdbuf());
        kernels = std::make_unique<KernelClass>(nbLevels, boxWidth, boxCenter, &matrixKernel);
        }
    };

FmmCore::FmmCore(const std::vector<double> &xyz, const int nbTargets, const int nbLevels,
                 const int nbThreads)
    {
    omp_set_num_threads(nbThreads);
    if (openblas_set_num_threads)
        { openblas_set_num_threads(1); }
    if (MKL_Set_Num_Threads)
        { MKL_Set_Num_Threads(1); }
    impl = std::make_unique<Impl>(nbTargets, nbLevels);
    const int nbParticles = xyz.size() / 3;
    for (int i = 0; i < nbParticles; i++)
        {
        const FPoint<FReal> p(xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]);
        if (i < nbTargets)
            { impl->tree.insert(p, FParticleType::FParticleTypeTarget, i); }
        else
            { impl->tree.insert(p, FParticleType::FParticleTypeSource, i, 0.0); }
        }
    // after omp_set_num_threads(): it copies the kernel for each thread
    impl->algo = std::make_unique<FmmClass>(&impl->tree, impl->kernels.get());
    }

FmmCore::~FmmCore() = default;

void FmmCore::compute(const std::vector<double> &q, std::vector<double> &potentials)
    {
    const int nbTargets = impl->nbTargets;
    // physicalValues[idxPart] = q, reset potentials
    impl->tree.forEachLeaf(
            [&q, nbTargets](LeafClass *leaf)
            {
                const int nbParticlesInLeaf = leaf->getSrc()->getNbParticles();
                const auto &indexes = leaf->getSrc()->getIndexes();
                FReal *const physicalValues = leaf->getSrc()->getPhysicalValues();
                for (int idxPart = 0; idxPart < nbParticlesInLeaf; ++idxPart)
                    { physicalValues[idxPart] = q[indexes[idxPart] - nbTargets]; }

                std::fill_n(leaf->getTargets()->getPotentials(),
                            leaf->getTargets()->getNbParticles(), 0);
            });

    impl->tree.forEachCell([](CellClass *cell) { cell->resetToInitialState(); });

    impl->algo->execute();

    potentials.resize(nbTargets);
    impl->tree.forEachLeaf(
            [&potentials](LeafClass *leaf)
            {
                const FReal *const pot = leaf->getTargets()->getPotentials();
                const int nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
                const auto &indexes = leaf->getTargets()->getIndexes();
                for (int idxPart = 0; idxPart < nbParticlesInLeaf; ++idxPart)
                    { potentials[indexes[idxPart]] = pot[idxPart]; }
            });
    }
    }  // namespace scal_fmm
