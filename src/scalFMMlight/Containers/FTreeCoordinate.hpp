// See LICENCE file at project root
#ifndef FTREECOORDINATE_HPP
#define FTREECOORDINATE_HPP

#include "../Utils/FPoint.hpp"

/**
 * @author Berenger Bramas (berenger.bramas@inria.fr)
 * @class FTreeCoordinate
 * Please read the license
 *
 * This class represents tree coordinate. It is used to save the position in "box unit" (not system/space unit!).
 * It is directly related to morton index, as interleaves bits from this coordinate make the morton index
 */
class FTreeCoordinate : public FPoint<int> {
public:
    /** Default constructor (position = {0,0,0})*/
    FTreeCoordinate(): FPoint<int>() {}

    /** constructor from Morton index  */
    explicit FTreeCoordinate(const MortonIndex mindex) { setPositionFromMorton(mindex); }

    /** Constructor from args
     * @param inX the x
     * @param inY the y
     * @param inZ the z
     */
    explicit FTreeCoordinate(const int inX,const int inY,const int inZ): FPoint<int>(inX, inY, inZ){}

    /** default copy constructor */
    FTreeCoordinate(const FTreeCoordinate&) = default;

    /**
     * Copy assignment
     * @param other the source class to copy
     * @return this a reference to the current object
     */
    FTreeCoordinate& operator=(const FTreeCoordinate& other) = default;

    /**
     * To get the morton index of the current position
     * @param inLevel the level of the component
     * @return morton index
     */
    MortonIndex getMortonIndex() const
        {
        MortonIndex index = 0x0LL;
        MortonIndex mask = 0x1LL;
        // the order is xyz.xyz...
        MortonIndex mx = FPoint<int>::data()[0] << 2;
        MortonIndex my = FPoint<int>::data()[1] << 1;
        MortonIndex mz = FPoint<int>::data()[2];

        while( (mask <= mz) || ((mask << 1) <= my) || ((mask << 2) <= mx))
            {
            index |= (mz & mask);
            mask <<= 1;
            index |= (my & mask);
            mask <<= 1;
            index |= (mx & mask);
            mask <<= 1;

            mz <<= 2;
            my <<= 2;
            mx <<= 2;
            }

        return index;
        }

    /** This function set the position of the current object using a morton index
     * @param inIndex the morton index to compute position
     */
    void setPositionFromMorton(MortonIndex inIndex) {
        MortonIndex mask = 0x1LL;

        FPoint<int>::data()[0] = 0;
        FPoint<int>::data()[1] = 0;
        FPoint<int>::data()[2] = 0;

        while(inIndex >= mask) {
            FPoint<int>::data()[2] |= int(inIndex & mask);
            inIndex >>= 1;
            FPoint<int>::data()[1] |= int(inIndex & mask);
            inIndex >>= 1;
            FPoint<int>::data()[0] |= int(inIndex & mask);

            mask <<= 1;
        }
    }
};

#endif //FTREECOORDINATE_HPP
