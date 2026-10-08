// See LICENCE file at project root
#ifndef FTYPEDLEAF_HPP
#define FTYPEDLEAF_HPP

#include "../Utils/FAssert.hpp"
#include "../Utils/FPoint.hpp"
#include "./FParticleType.hpp"

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FTypedLeaf
* @brief
* Please read the license
* This class is used to enable the use of typed particles
* (source XOR target) or simple system (source AND target).
*
* Particles should be typed to enable targets/sources difference.
*/
template<class FReal, class ContainerClass>
class FTypedLeaf
{
     /** The sources containers */
    ContainerClass sources;

     /** The targets containers */
    ContainerClass targets;

public:
    /** do-nothing destructor */
    ~FTypedLeaf(){}

    /**
     * To add a new particle in the leaf
     * @param inParticlePosition the position of the new particle
     * @param type to know if it is a target
     * @param args followed by other param given by the user (variadic template)
     */
    template<typename... Args>
    void push(const FPoint<FReal>& inParticlePosition, const FParticleType type, Args ... args)
        {
        if(type == FParticleType::FParticleTypeTarget) targets.push(inParticlePosition, args...);
        else sources.push(inParticlePosition, args...);
        }

    /**
    * To get all the sources in a leaf
    * @return a pointer to the list of particles that are sources
    */
    ContainerClass* getSrc() { return &this->sources; }

    /**
    * To get all the target in a leaf
    * @return a pointer to the list of particles that are targets
    */
    ContainerClass* getTargets() { return &this->targets; }

};

#endif //FTYPEDLEAF_HPP
