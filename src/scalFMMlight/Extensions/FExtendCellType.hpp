// See LICENCE file at project root
#ifndef FEXTENDCELLTYPE_HPP
#define FEXTENDCELLTYPE_HPP

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FExtendCellType
* This class is an extension, it proposes a target/source extension for cell.
* Because cells may have child that contains only sources or targets (in Tsm system) then it is important
* to not compute for nothing.
*/
class FExtendCellType {
protected:
    /** Current type = targets */
    bool containsTargets;

    /** Current type = sources */
    bool containsSources;

public:
    /** Default constructor */
    FExtendCellType() : containsTargets(false), containsSources(false) {}

    /** Copy constructor */
    FExtendCellType(const FExtendCellType& other) : containsTargets(other.containsTargets),
            containsSources(other.containsSources){}

    /** Copy operator */
    FExtendCellType& operator=(const FExtendCellType& other) {
        this->containsTargets = other.containsTargets;
        this->containsSources = other.containsSources;
        return *this;
    }

    /** To know if a cell has sources */
    bool hasSrcChild() const { return containsSources; }

    /** To know if a cell has targets */
    bool hasTargetsChild() const { return containsTargets; }

    /** To set cell as sources container */
    void setSrcChildTrue() { containsSources = true; }

    /** To set cell as targets container */
    void setTargetsChildTrue() { containsTargets = true; }

public:
    /** Save current object */
    template <class BufferWriterClass>
    void save(BufferWriterClass& buffer) const {
        buffer << containsTargets;
        buffer << containsSources;
    }
    /** Retrieve current object */
    template <class BufferReaderClass>
    void restore(BufferReaderClass& buffer) {
        buffer >> containsTargets;
        buffer >> containsSources;
    }
    /** reset to unknown type */
    void resetToInitialState(){
        containsTargets = false;
        containsSources = false;
    }

    /** return full size = targets + sources  */
    FSize getSavedSize() const { return FSize(sizeof(containsTargets) + sizeof(containsSources)); }
};


#endif //FEXTENDCELLTYPE_HPP


