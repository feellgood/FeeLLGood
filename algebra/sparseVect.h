#ifndef SPARSEVECT_H
#define SPARSEVECT_H

/** \file sparseVect.h
 * \brief sparse vector
a write sparse vector is a collection of v_coeff, which is a couple composed of an index and a double value
to populate with coefficients the sparse vector, use insert method
If several v_coeff have the same index, then they are automatically summed up.
a read sparse vector is a std::vector of v_coeff, built by its constructor
 **/

#include <map>
#include <vector>
#include <iostream>
#include <algorithm>

namespace algebra
{

/** overloading of << for printing functions */
template<typename T>
std::ostream& operator << (std::ostream& os, const std::vector<T>& v)
    {
    os << "[";
    for (auto it = v.begin(); it != v.end(); ++it)
        { os << "{" << it->_i << ":" << it->getVal() << "},"; }
    os << "]";
    return os;
    }

/**
\class w_sparseVect
it is a class for sparse vector in writing mode, using std::map<int,double> as coefficient container
*/
class w_sparseVect
{
public:
    /** inserter with a coefficient */
    inline void insert(v_coeff coeff) { coefs[coeff._i] += coeff.getVal(); }

    /** return true if the coefficient exists */
    inline bool exist(const int &idx) const { return coefs.find(idx) != coefs.end(); }

    /** getter for the value of a coefficient of index idx, if several coeffs have the same index then it returns the value of the first occurence */
    inline double getVal(const int idx) const
        {
        auto it = coefs.find(idx);
        if (it == coefs.end()) return 0;
        return it->second;
        }

    /** Insert the coefficients into vector v */
    void inorder_insert(std::vector<v_coeff> &v)
        {
        for (auto it = coefs.begin(); it != coefs.end(); ++it)
            v.push_back(v_coeff{it->first, it->second});
        }

private:
    /** coeffs container */
    std::map<int, double> coefs;
}; // end class w_sparseVect


/**
\class r_sparseVect
read sparse vector : it is a std::vector container for v_coeff, in reading mode
*/
class r_sparseVect: public std::vector<v_coeff>
{
public:
    /** default constructor */
    r_sparseVect(): std::vector<v_coeff>() {}

    /** constructor from a write sparse vector */
    r_sparseVect(w_sparseVect &v): std::vector<v_coeff>() { collect(v); }

    /** return true if the coefficient exists */
    inline bool exist(const int idx) const
        {
        return ( std::find_if(begin(),end(),[this,&idx](v_coeff coeff){return (coeff._i == idx);}) != end() );
        }

    /** collect method is sorting all v_coeffs, eventually with redundant indices, and is summing coeffs with same indices. It removes the coeffs that have been summed. */
    inline void collect(w_sparseVect &v)
        {
        clear();
        v.inorder_insert(*this);
        }

    /** getter for the value of a coefficient of index idx
    if several coeffs have the same index then it returns the value of the first occurence
    return zero if coefficient of index idx does not exist
     */
    inline double getVal(const int idx) const
        {
        double val(0);
        auto it = std::find_if(begin(),end(),[this,&idx](v_coeff coeff){return (coeff._i == idx); } );
        if (it != end()) val = it->getVal();
        return val;
        }

    /** scalar product */
    inline double dot(const std::vector<double> & X) const
        {
        double val(0);
        for(auto it=begin();it!=end();++it)
            { if(it->_i < (int)(X.size()) ) { val += it->getVal()*X[it->_i]; } }
        return val;
        }

    /** printing function */
    inline void print(std::ostream & flux) const
        {
        flux<<'{';
        std::for_each(begin(),end(), [&flux](const v_coeff &c)
            { flux << '{' << c._i << ':' << c.getVal() <<'}';});
        flux<<"}\n";
        }
}; // end class r_sparseVect

} // end namespace algebra

#endif
