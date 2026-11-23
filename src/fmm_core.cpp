#include "fmm_core.h"

#include <omp.h>

#include <algorithm>

#include "Components/FParticleType.hpp"
#include "Components/FTypedLeaf.hpp"
#include "Containers/FOctree.hpp"
#include "Core/FFmmAlgorithmThreadTsm.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"
#include "Kernels/Rotation/FRotationCell.hpp"
#include "Kernels/Rotation/FRotationKernel.hpp"

namespace scal_fmm
    {
namespace
    {
const int P = 9;              /**< truncation of the spherical harmonics series */
const int SizeSubLevels = 3;  /**< maximal size of the sub levels */

typedef double FReal; /**< all computations are made in double precision */
typedef FTypedRotationCell<FReal, P> CellClass;
typedef FP2PParticleContainerIndexed<FReal> ContainerClass;
typedef FTypedLeaf<FReal, ContainerClass> LeafClass;
typedef FOctree<FReal, CellClass, ContainerClass, LeafClass> OctreeClass;
typedef FRotationKernel<FReal, CellClass, ContainerClass, P> KernelClass;
typedef FFmmAlgorithmThreadTsm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>
        FmmClass;

const FPoint<FReal> boxCenter(0., 0., 0.); /**< center of the box of the octree */
    }  // namespace

/** octree, kernel and algorithm of the FMM. The particles are indexed in the order of xyz: targets
from 0 to nbTargets - 1, then sources. */
struct FmmCore::Impl
    {
    const int nbTargets; /**< number of targets */
    OctreeClass tree;    /**< octree */
    KernelClass kernels; /**< kernel */
    std::unique_ptr<FmmClass> algo; /**< algorithm, built once the particles are inserted */

    Impl(const int nbT, const int nbLevels)
        : nbTargets(nbT),
          tree(nbLevels, std::min(SizeSubLevels, nbLevels - 1), boxWidth, boxCenter),
          kernels(nbLevels, boxWidth, boxCenter)
        {}
    };

FmmCore::FmmCore(const std::vector<double> &xyz, const int nbTargets, const int nbLevels,
                 const int nbThreads)
    {
    omp_set_num_threads(nbThreads);
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
    impl->algo = std::make_unique<FmmClass>(&impl->tree, &impl->kernels);
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
