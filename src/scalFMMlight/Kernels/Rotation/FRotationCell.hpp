// See LICENCE file at project root
#ifndef FROTATIONCELL_HPP
#define FROTATIONCELL_HPP

#include <complex>

#include "../../Utils/FGlobal.hpp"
#include "../../Containers/FTreeCoordinate.hpp"
#include "../../Extensions/FExtendCellType.hpp"

/** This class is a cell used for the rotation based kernel
  * The size of the multipole and local vector are based on a template
  * User should choose this parameter P carefully to match with the P of the kernel.
  *
  * Multipole/Local vectors contain value as:
  * {0,0}{1,0}{1,1}...{P,P-1}{P,P}
  * So the size of such vector is (P+2)*(P+1)/2
  */
template <class FReal, int P>
class FRotationCell
{
protected:
    /** Size of multipole vector */
    static const int MultipoleSize = ((P+2)*(P+1))/2;

    /** Size of local vector */
    static const int LocalSize = ((P+2)*(P+1))/2;

    /** Multipole vector (static memory) for multipole extension */
    std::complex<FReal> multipole_exp[MultipoleSize];

    /** Local vector (static memory) for local extension */
    std::complex<FReal> local_exp[LocalSize];

    /** Morton index (need by most elements) */
    MortonIndex mortonIndex;

    /** The position */
    FTreeCoordinate coordinate;

    /** Level in tree */
    std::size_t level;

public:
    /** default constructor */
    FRotationCell() : mortonIndex(0) {}

    /** Copy constructor: Copy the value in the vectors */
    FRotationCell(const FRotationCell& other) { (*this) = other; }

    /** Default destructor */
    virtual ~FRotationCell() {}

    /** Copy operator: copies only the value in the vectors */
    FRotationCell& operator=(const FRotationCell& other)
        {
        for(int idx = MultipoleSize; 0 < idx ; --idx) { multipole_exp[idx] = other.multipole_exp[idx]; }
        for(int idx = LocalSize; 0 < idx ; --idx) { local_exp[idx] = other.local_exp[idx]; }
        return *this;
        }

    /** set the level */
    void setLevel(std::size_t inLevel) { this->level = inLevel; }

    /** To get the morton index */
    MortonIndex getMortonIndex() const { return this->mortonIndex; }

    /** To set the morton index */
    void setMortonIndex(const MortonIndex inMortonIndex) { this->mortonIndex = inMortonIndex; }

    /** To get the position */
    const FTreeCoordinate& getCoordinate() const { return this->coordinate; }

    /** To set the position from 3 indices */
    void setCoordinate(const int inX, const int inY, const int inZ) {
        this->coordinate.setX(inX);
        this->coordinate.setY(inY);
        this->coordinate.setZ(inZ);
    }

    /** Get Multipole array */
    const std::complex<FReal>* getMultipole() const { return multipole_exp; }

    /** Get Local array */
    const std::complex<FReal>* getLocal() const { return local_exp; }

    /** Get Multipole array */
    std::complex<FReal>* getMultipole() { return multipole_exp; }

    /** Get Local array */
    std::complex<FReal>* getLocal() { return local_exp; }

    /** Reset all values to zero */
    void resetToInitialState()
        {
        for(int idx = 0 ; idx < MultipoleSize ; ++idx)
            { multipole_exp[idx] = std::complex<FReal> {0.,0.}; }
        for(int idx = 0 ; idx < LocalSize ; ++idx)
            { local_exp[idx] = std::complex<FReal> {0.,0.}; }
        }
};

template <class FReal, int P>
class FTypedRotationCell : public FRotationCell<FReal, P>, public FExtendCellType {
public:
    /** Reset to initial state the cell */
    void resetToInitialState(){
        FRotationCell<FReal, P>::resetToInitialState();
        FExtendCellType::resetToInitialState();
    }
};

#endif // FROTATIONCELL_HPP
