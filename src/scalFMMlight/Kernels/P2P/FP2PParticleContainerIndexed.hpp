// See LICENCE file at project root
#ifndef FP2PPARTICLECONTAINERINDEXED_HPP
#define FP2PPARTICLECONTAINERINDEXED_HPP

#include <vector>

#include "../../Utils/FGlobal.hpp"
#include "../../Utils/FAlignedMemory.hpp"
#include "../../Utils/FMath.hpp"
#include "../../Utils/FPoint.hpp"
#include "../../Components/FParticleType.hpp"

const int NbAttributesPerParticle = 2; // 2 = 1 + 1 : 1 for the scalar charges, 1 for the scalar potential

/** container for indexed particles */
template<class FReal>
class FP2PParticleContainerIndexed
{
protected:
    /** size of a chunck to align memory */
    static const FSize MemoryAlignement = FP2PDefaultAlignement;

    /** the indices of the particles */
    std::vector<FSize> indexes;

    /** The number of particles in the container */
    FSize nbParticles;

    /** 3 pointers to 3 arrays of FReal to store (X,Y,Z) positions */
    FReal* positions[3];

    /** The attributes of the particles */
    FReal* attributes[NbAttributesPerParticle];

    /** The allocated memory */
    FSize allocatedParticles;

    /** Ending call for pushing the attributes */
    template<int index>
    void addParticleValue(const FSize /*insertPosition*/){}

    /** Recursive template : filling call for each attributes values */
    template<int index, typename... Args>
    void addParticleValue(const FSize insertPosition, const FReal value, Args... args)
        {
        // Compile test to ensure indexing
        static_assert(index < NbAttributesPerParticle, "Index to get attributes is out of scope.");
        // insert the value
        attributes[index][insertPosition] = value;
        // Continue for remaining values
        addParticleValue<index+1>( insertPosition, args...);
        }

    /** increase position size if needed. nbParticles increases by step = 1 while allocatedParticles formula has a 1.5 multiplier: it might be unoptimized */
    void increaseSizeIfNeeded(void)
        {
        static const FSize DefaultNbParticles = FSize(MemoryAlignement/sizeof(FReal));
        if( nbParticles >= allocatedParticles )
            {
            // allocate memory
            allocatedParticles = (FMath::Max(DefaultNbParticles,FSize( FReal(nbParticles+1)*1.5 )) + DefaultNbParticles - 1) & ~(DefaultNbParticles-1);
            // init with 0
            const size_t allocatedBytes = sizeof(FReal)*(3 + NbAttributesPerParticle)*allocatedParticles;
            FReal (*newData)[allocatedParticles] = reinterpret_cast<FReal (*)[allocatedParticles]>(FAlignedMemory::AllocateBytes<MemoryAlignement>(allocatedBytes));
            memset( newData, 0, allocatedBytes);
            const char*const toDelete  = reinterpret_cast<const char*>(positions[0]);
            // copy memory, position components and attributes
            if (nbParticles != 0)
                {
                memcpy(newData[0], positions[0], sizeof newData[0]);
                memcpy(newData[1], positions[1], sizeof newData[1]);
                memcpy(newData[2], positions[2], sizeof newData[2]);
                memcpy(newData[3], attributes[0], sizeof newData[3]);
                memcpy(newData[4], attributes[1], sizeof newData[4]);
                }
            positions[0] = newData[0];
            positions[1] = newData[1];
            positions[2] = newData[2];
            attributes[0] = newData[3];
            attributes[1] = newData[4];
            // delete old
            FAlignedMemory::DeallocBytes<MemoryAlignement>(toDelete);
            }
        }

public:
    /** Basic contructor */
    FP2PParticleContainerIndexed() : nbParticles(0), allocatedParticles(0)
        {
        memset(positions, 0, sizeof(positions[0]) * 3);
        memset(attributes, 0, sizeof(attributes[0]) * NbAttributesPerParticle);
        }
    
    /** constructor copy is deleted */
    FP2PParticleContainerIndexed(FP2PParticleContainerIndexed&) = delete;
    
    /** operator copy is deleted */
    FP2PParticleContainerIndexed& operator=(const FP2PParticleContainerIndexed&) = delete;

    /** Destructor: dealloc the memory using first pointer */
    ~FP2PParticleContainerIndexed() { FAlignedMemory::DeallocBytes<MemoryAlignement>(positions[0]); }

    /**
   * @brief getNbParticles
   * @return the number of particles
   */
    FSize getNbParticles() const { return nbParticles; }

    /**
   * @brief getPositions
   * @return a FReal*[3] to get access to the positions
   */
    const FReal*const* getPositions() const { return positions; }

    /**
   * @brief getPositions
   * @return get the position in write mode
   */
    FReal* const* getPositions() { return positions; }

    /** Push to fill both position and corresponding index, called by FSimpleLeaf
     * Should have a particle position, type, index and followed by attributes 
     */
    template<typename... Args>
    void push(const FPoint<FReal>& inParticlePosition,
              const FParticleType particleType,
              const FSize index, Args... args)
        {
        increaseSizeIfNeeded();

        // insert particle data
        positions[0][nbParticles] = inParticlePosition.getX();
        positions[1][nbParticles] = inParticlePosition.getY();
        positions[2][nbParticles] = inParticlePosition.getZ();
        // insert attribute data
        addParticleValue<0>( nbParticles, args...);
        nbParticles += 1;
        indexes.push_back(index);
        }
    
    /** getter for the charges */
    inline FReal* getPhysicalValues(void) { return attributes[0]; }

    /** const getter for the charges */
    inline const FReal* getPhysicalValues(void) const { return attributes[0]; }

    /** getter for the potential */
    inline FReal* getPotentials(void) { return attributes[1]; }

    /** const getter for the potential */
    inline const FReal* getPotentials(void) const { return attributes[1]; }

    /** getter for a reference to the vector of indices */
    const std::vector<FSize>& getIndexes() const { return indexes; }
};

#endif // FP2PPARTICLECONTAINERINDEXED_HPP
