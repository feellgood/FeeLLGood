#ifndef SPARSEVECT_H
#define SPARSEVECT_H

/** \file sparseVect.h
 * \brief sparse vector
A write-mode sparse vector is an (index -> value) map.

A read-mode sparse vector holds a list of sorted indices, and a list of matching values. The list
of indices is immutable once the vector has been constructed, but the values can later be changed.
 **/

#include <set>
#include <map>
#include <vector>
#include <iostream>
#include <algorithm>

namespace algebra
{

/**
\class v_coeff
A coefficient of a sparse vector is an (index, value) pair.
*/
using v_coeff = std::pair<int, double>;

/**
\class w_sparseVect
A sparse vector in writing mode is an (index -> value) map.
*/
using w_sparseVect = std::map<int, double>;

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
        for (auto it = v.begin(); it != v.end(); ++it)
            {
            indices.push_back(it->first);
            values.push_back(it->second);
            }
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
        size_t k = 0;
        for (; k < indices.size(); ++k)
            if (indices[k] == i)
                break;
        assert(k < indices.size());
        values[k] += val;
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
