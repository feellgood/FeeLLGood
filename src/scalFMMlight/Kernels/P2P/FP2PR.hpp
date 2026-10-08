// See LICENCE file at project root

#ifndef FP2PR_HPP
#define FP2PR_HPP

#include "../../Utils/FGlobal.hpp"
#include "../../Utils/FMath.hpp"

/**
 * @brief FullRemote template function
 */

template <class FReal, class ContainerClass>
static void FullRemote(ContainerClass* const FRestrict inTargets, const ContainerClass* const inNeighbors[], const int limitNeighbors)
    {
    const FSize nbParticlesTargets = inTargets->getNbParticles();
    const FReal*const targetsX = inTargets->getPositions()[0];
    const FReal*const targetsY = inTargets->getPositions()[1];
    const FReal*const targetsZ = inTargets->getPositions()[2];
    FReal*const targetsPotentials = inTargets->getPotentials();

    for(FSize idxNeighbors = 0 ; idxNeighbors < limitNeighbors ; ++idxNeighbors)
        {
        if( inNeighbors[idxNeighbors] )
            {
            const FSize nbParticlesSources = inNeighbors[idxNeighbors]->getNbParticles();
            const FReal*const sourcesPhysicalValues = (const FReal*)inNeighbors[idxNeighbors]->getPhysicalValues();
            const FReal*const sourcesX = (const FReal*)inNeighbors[idxNeighbors]->getPositions()[0];
            const FReal*const sourcesY = (const FReal*)inNeighbors[idxNeighbors]->getPositions()[1];
            const FReal*const sourcesZ = (const FReal*)inNeighbors[idxNeighbors]->getPositions()[2];

            for(FSize idxTarget = 0 ; idxTarget < nbParticlesTargets ; ++idxTarget)
                {
                const FReal tx = targetsX[idxTarget];
                const FReal ty = targetsY[idxTarget];
                const FReal tz = targetsZ[idxTarget];
                FReal tpo(0.0);

                for(FSize idxSource = 0 ; idxSource < nbParticlesSources ; ++idxSource)
                    {
                    FReal dx = tx - sourcesX[idxSource];
                    FReal dy = ty - sourcesY[idxSource];
                    FReal dz = tz - sourcesZ[idxSource];
                    tpo += sourcesPhysicalValues[idxSource] / FMath::Sqrt(dx*dx + dy*dy + dz*dz);
                    }
                targetsPotentials[idxTarget] += tpo;
                }
            }
        }
    }

#endif // FP2PR_HPP
