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
    v_coeff(const int i,const double c):_i(i),_c(c) {}

    /** getter for the value of the coefficient */
    double getVal() const {return _c;}

    /** index of the coefficient */
    int _i;

/** printing function */
    void print(std::ostream & flux) const { flux << "(" << _i << ";" << _c << ")"; }

private:
/** value of the coefficient */
    double _c;
}; //end class v_coeff

} // end namespace algebra

#endif
