// See LICENCE file at project root
#ifndef FALIGNEDMEMORY_HPP
#define FALIGNEDMEMORY_HPP

#include <cstdint>


/**
 * This should be used to allocate and deallocate aligned memory to AlignementValue.
 */
namespace FAlignedMemory {
/** Specifies the alignment requirement */
template <std::size_t AlignementValue>
struct alignas(AlignementValue) aligned_block {
    /** char array to "alignas" */
    char bytes[AlignementValue];
};

/** alloc memory */
template <std::size_t AlignementValue>
inline void* AllocateBytes(const std::size_t inSize){
    if(inSize == 0){
        return nullptr;
    }

    // Ensure alignment is a power of 2
    static_assert(AlignementValue != 0 && ((AlignementValue-1)&AlignementValue) == 0, "Alignement must be a power of 2");

    std::size_t blockCount = (inSize + AlignementValue - 1) / AlignementValue;
    return new aligned_block<AlignementValue>[blockCount];
}

/** free memory */
template <std::size_t AlignementValue>
inline void DeallocBytes(const void* ptrToFree){
    delete[] reinterpret_cast<const aligned_block<AlignementValue>*>(ptrToFree);
}

}

#endif // FALIGNEDMEMORY_HPP
