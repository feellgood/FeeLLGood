#ifndef SPARSEMAT_H
#define SPARSEMAT_H

/** \file sparseMat.h 
 \brief read and write sparse matrix
r_sparseMat : read sparse matrix : it is buit calling the constructor with a w_sparseMat as argument
w_sparseMat : write sparse matrix : a std::vector of m_coeff = triplet of two indices and a value (i,j,value)
 */

#include <iostream>
#include <vector>
#include <cassert>

namespace algebra
{

/**
\class w_sparseMat
write sparse Matrix, it is a container for coefficients of a 'line' sparse matrix.
If some m_coeff have the same indices, they will be summed to build the final matrix coefficient
*/
class w_sparseMat
{
	friend class r_sparseMat;

public:
/** constructor */
    inline w_sparseMat(const int _N):N(_N) { m.resize(N); }

/** inserter for a coefficient val at line i col j */
    inline void insert(const int i, const int j, const double val)
        {
        assert(i<N && j<N );
        m[i].insert({j, val});
        }

/** getter for the number of lines */
	inline int getDim(void) const {return N;}

private:
/** dimension of sparse matrix, N is the number of lines */
    const int N;

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
	inline r_sparseVect & operator() (const int & i) {return m[i];}

/** getter for a coefficient value */
	inline double operator() (const int &i, const int &j) const { return m[i].getVal(j); }

private:
/** dimension of the sparse matrix (nb of lines) */
	const int N;

/** coefficient container */
	std::vector<r_sparseVect> m;
}; // end class r_sparseMat



} // end namespace algebra

#endif
