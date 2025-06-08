#ifndef SPARSEVECT_H
#define SPARSEVECT_H

/** \file sparseVect.h
 * \brief sparse vector
a write sparse vector is a collection of v_coeff, which is a couple composed of an index and a double value
to populate with coefficients the sparse vector, use insert method
If several v_coeff have the same index, then they are automatically summed up.
a read sparse vector is a std::vector of v_coeff, built by its constructor
 **/

#include <set>
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
    /** Constructor. */
    w_sparseVect(int dimension) : N(dimension) {}

    /** inserter with a coefficient */
    void insert(int i, double val)
        {
        if (i < 0 || i >= N)
            { std::cerr << "Error: out-of-range index in sparse vector.\n"; exit(1); }
        coefs[i] += val;
        }

    /** Insert the coefficients into vector v */
    void inorder_insert(std::vector<int> &indices, std::vector<double> &values) const
        {
        for (auto it = coefs.begin(); it != coefs.end(); ++it)
            {
            indices.push_back(it->first);
            values.push_back(it->second);
            }
        }

    /** Return the number of elements */
    size_t size() const { return coefs.size(); }

    /** Dimension of the vector. All indices lie within [0, N). */
    const int N;

private:
    /** coeffs container */
    std::map<int, double> coefs;
}; // end class w_sparseVect


/**
\class r_sparseVect
read sparse vector : it is a std::vector container for v_coeff, in reading mode
*/
class r_sparseVect
{
public:
    /** constructor from a write sparse vector */
    r_sparseVect(const w_sparseVect &v)
        {
        indices.reserve(v.size());
        values.reserve(v.size());
        v.inorder_insert(indices, values);
        }

    /** constructor from size and shape */
    r_sparseVect(const std::set<int>& shape)
        {
        values.resize(shape.size());
        indices.reserve(shape.size());
        for (auto it = shape.begin(); it != shape.end(); ++it)
            indices.push_back(*it);
        }

    /** zero all elements, while preserving the shape */
    void clear()
        {
        for (double& value: values)
            value = 0;
        }

    /** add the value at position i, which must belog to the shape */
    void add(int i, double val)
        {
        for (size_t k = 0; k < indices.size(); ++k)
            if (indices[k] == i)
                {
                values[k] += val;
                return;
                }
        std::cerr << "Error: invalid index in add().\n";
        exit(1);
        }

    /** getter for the value of a coefficient of index idx
    if several coeffs have the same index then it returns the value of the first occurence
    return zero if coefficient of index idx does not exist
     */
    double getVal(const int idx) const
        {
        for (size_t i = 0; i < indices.size(); ++i)
            if (indices[i] == idx)
                { return values[i]; }
        return 0;
        }

    /** Scalar product. The caller must ensure the vector dimensions match. */
    double dot(const std::vector<double> & X) const
        {
        double val(0);
        for (size_t i = 0; i < indices.size(); ++i)
            { val += values[i] * X[indices[i]]; }
        return val;
        }

    /** printing function */
    void print(std::ostream & flux) const
        {
        flux<<'{';
        for (size_t i = 0; i < indices.size(); ++i)
            { flux << '{' << indices[i] << ':' << values[i] << '}'; }
        flux<<"}\n";
        }

private:
    std::vector<int> indices;  /**< array of vector indices. */
    std::vector<double> values;  /**< array of vector values matching the indices */
}; // end class r_sparseVect

} // end namespace algebra

#endif
