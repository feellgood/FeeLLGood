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
\class w_sparseVect
A sparse vector in writing mode is an (index -> value) map.
*/
using w_sparseVect = std::map<int, double>;

/**
\class r_sparseVect
read sparse vector: holds a list of sorted indices, and a list of matching values.
*/
struct r_sparseVect
{
    std::vector<int> indices;  /**< array of vector indices. */
    std::vector<double> values;  /**< array of vector values matching the indices */
}; // end class r_sparseVect

} // end namespace algebra

#endif
