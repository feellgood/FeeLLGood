// See LICENCE file at project root
#ifndef FTREECOORDINATE_HPP
#define FTREECOORDINATE_HPP

/**
 * @author Berenger Bramas (berenger.bramas@inria.fr)
 * @class FTreeCoordinate
 * Please read the license
 *
 * This class represents tree coordinate. It is used to save the position in "box unit" (not system/space unit!).
 * It is directly related to morton index, as interleaves bits from this coordinate make the morton index
 */
class FTreeCoordinate
{
private:
    std::array<int,3> ind;

public:
    /** Default constructor (position = {0,0,0})*/
    FTreeCoordinate() { ind = {0}; }

    /** constructor from Morton index  */
    explicit FTreeCoordinate(const MortonIndex mindex) { setPositionFromMorton(mindex); }

    /** Constructor from args
     * @param inX the x
     * @param inY the y
     * @param inZ the z
     */
    explicit FTreeCoordinate(const int inX,const int inY,const int inZ) { ind = {inX,inY,inZ}; }

    /** default copy constructor */
    FTreeCoordinate(const FTreeCoordinate&) = default;

    /**
     * Copy assignment
     * @param other the source class to copy
     * @return this a reference to the current object
     */
    FTreeCoordinate& operator=(const FTreeCoordinate& other) = default;

    int getX(void) const {return ind[0];}
    int getY(void) const {return ind[1];}
    int getZ(void) const {return ind[2];}

    void setX(const int inX) {ind[0] = inX;}
    void setY(const int inY) {ind[1] = inY;}
    void setZ(const int inZ) {ind[2] = inZ;}


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
        MortonIndex mx = ind[0] << 2;
        MortonIndex my = ind[1] << 1;
        MortonIndex mz = ind[2];

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
        ind[0]=0;
        ind[1]=0;
        ind[2]=0;

        while(inIndex >= mask) {
            ind[2] |= int(inIndex & mask);
            inIndex >>= 1;
            ind[1] |= int(inIndex & mask);
            inIndex >>= 1;
            ind[0] |= int(inIndex & mask);

            mask <<= 1;
        }
    }
};

#endif //FTREECOORDINATE_HPP
