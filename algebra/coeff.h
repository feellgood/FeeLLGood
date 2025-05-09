#ifndef COEFF_H
#define COEFF_H

/** \file coeff.h
 * \brief set of class to handle sparse matrix and vector coefficients
 * two dedicated classes vector and matrix coefficients
 it is permissive: indices might be negative. This might be usefull when relative indexation is needed
 */

namespace algebra
{
/**
\class v_coeff
container for pairs of indice nd double value to represent a coefficient of a sparse vector
*/
class v_coeff
{
public:
    /** constructor */
    inline v_coeff(const int i,const double c):_i(i),_c(c) {}

    /** getter for the value of the coefficient */
    inline double getVal(void) const {return _c;}

    /** setter for the value of the coefficient */
    inline void setVal(const double val) {_c = val;}

/** ref to coeff value */
inline double & valRef(void) {return _c;}

    /** increment value with val */
    inline void add(const double val) { _c += val;}
/**
lexicographic order
*/
    inline bool operator< (const v_coeff &c) const
    { return (this->_i < c._i); }

/**
two coeffs are equal when their indices are equals
*/
    inline bool operator== (const v_coeff &c) const
    { return (this->_i == c._i); }

    /** index of the coefficient */
    int _i;

/** printing function */
    inline void print(std::ostream & flux) const { flux << "(" << _i << ";" << _c << ")"; }

private:
/** value of the coefficient */
    double _c;
}; //end class v_coeff

/** operator<< for v_coeff */
inline std::ostream & operator<<(std::ostream & flux, v_coeff const& v) {v.print(flux); return flux;}

} // end namespace algebra

#endif
