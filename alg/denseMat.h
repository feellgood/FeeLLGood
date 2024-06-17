#ifndef ALG_DENSEMAT_H
#define ALG_DENSEMAT_H

/** \file alg_denseMat.h 
 \brief dense matrix
 */

#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>

namespace alg
{

/**
\class denseMat
usual dense matrix, container for values is a std::vector. Recommended for small dense matrices.
*/

class denseMat : public std::vector<double> {
  protected:
/** number of lines */
	size_t _nrows;

/** number of columns */
	size_t _ncols;

  public:
/** constructor */
inline denseMat(std::vector<double> &v, size_t l, size_t c): std::vector<double>(), _nrows(l), _ncols(c){
      assert(l*c == v.size());
      assign(v.begin(),v.end());
      }

/** constructor */
inline denseMat(size_t l, size_t c): std::vector<double>(c*l), _nrows(l), _ncols(c)  { }

/** constructor with an initial value for all coefficients */
inline denseMat(size_t l, size_t c, double value): std::vector<double>(c*l, value), _nrows(l), _ncols(c)  { }

/** constructor for an empty dense matrix */
inline denseMat(void) { _nrows = _ncols = 0; }

/** getter for the value of the coefficient */
inline const double& operator ()(size_t l, size_t c) const {
      assert(l < _nrows && c < _ncols);
      return *(this->begin() + c*_nrows+l);
    }

/** setter for the value of the coefficient */    
inline double& operator ()(size_t l, size_t c) {
      assert(l < _nrows && c < _ncols);
      return *(this->begin() + c*_nrows+l);
    }

/** convert back a denseMat to a std::vector */
    const std::vector<double> &as_vector(void) const { return *this; }

/** setter for the number of lines */
    inline size_t& nrows(void) { return _nrows; }

/** getter for the number of lines */
    inline const size_t& nrows(void) const { return _nrows; }

/** setter for the number of columns */
    inline size_t& ncols(void) { return _ncols; }

/** getter for the number of columns */
    inline const size_t& ncols(void) const { return _ncols; }

/** set all coefficients to zero */
    inline void clear()
        { std::fill(this->begin(),this->end(), 0.0); }

    /** printing function */
    inline void print(std::ostream & flux) const 
    { 
      for (size_t l=0; l<_nrows; ++l){
          flux<<'{'; 
          for (size_t c=0; c<_ncols; ++c){ 
              flux << *(this->begin() + c*_nrows+l) << " "; 
              }
          flux<<"}\n"; 
          }
    }
}; // end class denseMat

/** operator<< for denseMat */
inline std::ostream & operator<<(std::ostream & flux, denseMat const& m) {m.print(flux); return flux;}
}
#endif //ALG_DENSEMAT_H

