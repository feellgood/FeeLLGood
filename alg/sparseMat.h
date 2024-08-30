#ifndef ALG_SPARSEMAT_H
#define ALG_SPARSEMAT_H

/** \file alg_sparseMat.h 
 \brief read and write sparse matrix
r_sparseMat : read sparse matrix : it is buit calling the constructor with a w_sparseMat as argument
w_sparseMat : write sparse matrix : a std::vector of m_coeff = triplet of two indices and a value (i,j,value)
 */

#include <iostream>
#include <vector>
#include <cassert>

namespace alg
{

/**
\class w_sparseMat
write sparse Matrix, it is a container for objects m_coeff. 
If some m_coeff have the same indices, they will be summed to build the final matrix coefficient
*/
class w_sparseMat
{
	friend class r_sparseMat;

public:
/** constructor */
    inline w_sparseMat(const int _N):N(_N) { m.resize(N); }

/** inserter for a coefficient */
    inline void insert(const m_coeff &co)
        {
        assert(co._i<N && co._j<N );
        m[co._i].insert({co._j, co.getVal()});
        }

/** getter for the number of lines */
	inline int getDim(void) const {return N;}

private:
/** dimension of sparse matrix, N is the number of lines */
    int N;

/** container for the write sparse matrix coefficients */
    std::vector<w_sparseVect> m;
}; // end class w_sparseMat


/** \class r_sparseMat
read sparse matrix	 
The constructor is buiding from a write sparse matrix the data to access efficiently the coefficients values
*/
class r_sparseMat
{
public:
	/** constructor */
	inline r_sparseMat(w_sparseMat &A):N(A.getDim())
		{
		m.resize(N);// N is the number of lines
		if (!A.m.empty())
			{ for(int i=0; i<N; ++i){ m[i].collect(A.m[i]); } }
		}

/** printing function */
	inline void print(void) { std::for_each(m.begin(),m.end(),[](r_sparseVect const& _v) {std::cout << _v;} ); }

/** printing function */
	inline void print(std::ostream & flux) const
	{ std::for_each(m.begin(),m.end(),[&flux](r_sparseVect const& _v) {_v.print(flux);} ); }

	/** getter for the number of lines */
	inline int getDim(void) const {return N;}

/** return true if the coefficient exists */
	inline bool exist(const int &i, const int &j) const { return ( (i<N)&&(m[i].exist(j)) ); }

/** getter for an innner sparse vector */
	inline alg::r_sparseVect & operator() (const int & i) {return m[i];}

/** getter for a coefficient value */
	inline double operator() (const int &i, const int &j) const { return m[i].getVal(j); }

/** setter for a coefficient value */
	inline void setVal (const int &i, const int &j, const double val) { return m[i].setVal(j, val); }

private:
/** dimension of the sparse matrix (nb of lines) */
	const int N;

/** coefficient container */
	std::vector<r_sparseVect> m;
}; // end class r_sparseMat

/** operator<< for r_sparseMat */
inline std::ostream & operator<<(std::ostream & flux, r_sparseMat const& m) {m.print(flux); return flux;}

} // end namespace alg

#endif //ALG_SPARSEMAT_H
