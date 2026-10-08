#ifndef FCOORDINATECOMPUTER_HPP
#define FCOORDINATECOMPUTER_HPP

#include "./FTreeCoordinate.hpp"
#include "../Utils/FPoint.hpp"
#include "../Utils/FMath.hpp"
#include "../Utils/FAssert.hpp"

/**
 * @brief FCoordinateComputer namespace is containing four static functions to manipulate coordinate
 * from the simulation box properties.
 */
namespace FCoordinateComputer
{
    /** return coordinate index */
    template <class FReal>
    static inline int GetTreeCoordinate(const FReal inRelativePosition, const FReal boxWidth,
                                        const FReal boxWidthAtLeafLevel, const int treeHeight)
        {
        FAssertLF( (inRelativePosition >= 0 && inRelativePosition <= boxWidth), "inRelativePosition : ",inRelativePosition, " boxWidth ", boxWidth );
        if(inRelativePosition == boxWidth) { return FMath::pow2(treeHeight-1)-1; }
        const FReal indexFReal = inRelativePosition / boxWidthAtLeafLevel;
        return static_cast<int>(indexFReal);
        }

    /** return coordinate index from physical position */
    template <class FReal>
    static inline FTreeCoordinate GetCoordinateFromPosition(const FPoint<FReal>& centerOfBox,
                                                            const FReal boxWidth,
                                                            const int treeHeight,
                                                            const FPoint<FReal>& pos)
        {
        const FPoint<FReal> boxCorner(centerOfBox,-(boxWidth/2));
        const FReal boxWidthAtLeafLevel(boxWidth/FReal(1<<(treeHeight-1)));
        const FReal x = GetTreeCoordinate<FReal>( pos.getX() - boxCorner.getX(), boxWidth, boxWidthAtLeafLevel, treeHeight);
        const FReal y = GetTreeCoordinate<FReal>( pos.getY() - boxCorner.getY(), boxWidth, boxWidthAtLeafLevel, treeHeight);
        const FReal z = GetTreeCoordinate<FReal>( pos.getZ() - boxCorner.getZ(), boxWidth, boxWidthAtLeafLevel, treeHeight);
        return FTreeCoordinate(x,y,z);
        }

    /** return coordinate index from physical position relative to corner */
    template <class FReal>
    static inline FTreeCoordinate GetCoordinateFromPositionAndCorner(const FPoint<FReal>& cornerOfBox,
                                                                     const FReal boxWidth,
                                                                     const int treeHeight,
                                                                     const FPoint<FReal>& pos)
        {
        const FReal boxWidthAtLeafLevel(boxWidth/FReal(1<<(treeHeight-1)));
        // position has to be relative to corner not center
        const FReal x = GetTreeCoordinate<FReal>( pos.getX() - cornerOfBox.getX(), boxWidth, boxWidthAtLeafLevel, treeHeight);
        const FReal y = GetTreeCoordinate<FReal>( pos.getY() - cornerOfBox.getY(), boxWidth, boxWidthAtLeafLevel, treeHeight);
        const FReal z = GetTreeCoordinate<FReal>( pos.getZ() - cornerOfBox.getZ(), boxWidth, boxWidthAtLeafLevel, treeHeight);
        return FTreeCoordinate(x,y,z);
        }

    /** return position from coordinate */
    template <class FReal>
    static inline FPoint<FReal> GetPositionFromCoordinate(const FPoint<FReal>& centerOfBox,
                                                          const FReal boxWidth,
                                                          const int treeHeight,
                                                          const FTreeCoordinate& pos)
        {
        const FPoint<FReal> boxCorner(centerOfBox,-(boxWidth/2));
        const FReal boxWidthAtLeafLevel(boxWidth/FReal(1<<(treeHeight-1)));
        // position has to be relative to corner not center
        const FReal x = pos.getX()*boxWidthAtLeafLevel + boxCorner.getX();
        const FReal y = pos.getY()*boxWidthAtLeafLevel + boxCorner.getY();
        const FReal z = pos.getZ()*boxWidthAtLeafLevel + boxCorner.getZ();
        return FTreeCoordinate(x,y,z);
        }
}

#endif // FCOORDINATECOMPUTER_HPP

