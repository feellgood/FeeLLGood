#ifndef FG_H5UTILS_H
#define FG_H5UTILS_H

/** \file h5Utils.h
 * HDF5 helpers shared by fg-sol2hdf5 and fg-h5extract. Attributes are written the way h5py
 * writes Python values, so that files written by the Python and the C++ tools are
 * interchangeable: strings as scalar variable length UTF-8 strings, numbers as scalar
 * doubles, tables as two dimensional datasets of doubles.
 */

#include <string>
#include <vector>

#include <H5Cpp.h>

namespace fg
    {
/** Write a string attribute named name on obj. */
void writeStringAttribute(H5::H5Object &obj, const std::string &name, const std::string &value);

/** Write a double attribute named name on obj. */
void writeDoubleAttribute(H5::H5Object &obj, const std::string &name, double value);

/** Read the string attribute named name of obj, fixed or variable length. Throws
 * H5::Exception if it does not exist, std::runtime_error if it is not a string. */
std::string readStringAttribute(const H5::H5Object &obj, const std::string &name);

/** Value of an attribute as text: a string as is, a scalar number formatted as Python prints
 * it, anything else as a short description of its type. */
std::string attributeToString(const H5::Attribute &attr);

/** Create in group a dataset named name holding rows, a table of doubles. All the rows must
 * have the same length; throws std::runtime_error otherwise. */
H5::DataSet writeTable(H5::Group &group, const std::string &name,
                       const std::vector<std::vector<double>> &rows);

/** Read a dataset of numbers as a table of doubles. A one dimensional dataset gives one
 * row per value. */
std::vector<std::vector<double>> readTable(const H5::DataSet &dataset);

/** Names of the attributes of obj, in alphabetical order, as h5py lists them. */
std::vector<std::string> attributeNames(const H5::H5Object &obj);

/** Names of the members of group, in alphabetical order, as h5py lists them. */
std::vector<std::string> memberNames(const H5::Group &group);
    }  // namespace fg

#endif
