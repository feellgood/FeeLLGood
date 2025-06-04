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
        { os << "{" << it->first << ":" << it->second << "},"; }
    os << "]";
    return os;
    }

/**
\class v_coeff
A coefficient of a sparse vector is an (index, value) pair.
*/
using v_coeff = std::pair<int, double>;

/**
\class w_sparseVect
it is a class for sparse vector in writing mode, using std::map<int,double> as coefficient container
*/
class w_sparseVect
{
public:
    /** inserter with a coefficient */
    void insert(v_coeff coeff) { coefs[coeff.first] += coeff.second; }

    /** Insert the coefficients into vector v */
    void inorder_insert(std::vector<v_coeff> &v) const
        {
        for (auto it = coefs.begin(); it != coefs.end(); ++it)
            v.push_back(v_coeff{it->first, it->second});
        }

    /** Return the number of elements */
    size_t size() const { return coefs.size(); }

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
    /** constructor from a write sparse vector */
    r_sparseVect(const w_sparseVect &v): std::vector<v_coeff>()
        {
        reserve(v.size());
        v.inorder_insert(*this);
        }

    /** getter for the value of a coefficient of index idx
    if several coeffs have the same index then it returns the value of the first occurence
    return zero if coefficient of index idx does not exist
     */
    double getVal(const int idx) const
        {
        double val(0);
        auto it = std::find_if(begin(),end(),
                [this,&idx](v_coeff coeff){return (coeff.first == idx); } );
        if (it != end()) val = it->second;
        return val;
        }

    /** scalar product */
    double dot(const std::vector<double> & X) const
        {
        double val(0);
        for(auto it=begin();it!=end();++it)
            { if(it->first < (int)(X.size()) ) { val += it->second * X[it->first]; } }
        return val;
        }

    /** printing function */
    void print(std::ostream & flux) const
        {
        flux<<'{';
        std::for_each(begin(),end(), [&flux](const v_coeff &c)
            { flux << '{' << c.first << ':' << c.second <<'}';});
        flux<<"}\n";
        }
}; // end class r_sparseVect

} // end namespace algebra

#endif
