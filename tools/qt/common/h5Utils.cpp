#include "h5Utils.h"

#include <algorithm>
#include <stdexcept>

#include "pyRepr.h"

namespace fg
    {
namespace
    {
/** Variable length UTF-8 string type, the one h5py uses for Python str. */
H5::StrType vlenUtf8()
    {
    H5::StrType type(H5::PredType::C_S1, H5T_VARIABLE);
    type.setCset(H5T_CSET_UTF8);
    return type;
    }
    }  // namespace

void writeStringAttribute(H5::H5Object &obj, const std::string &name, const std::string &value)
    {
    H5::Attribute attr = obj.createAttribute(name, vlenUtf8(), H5::DataSpace(H5S_SCALAR));
    attr.write(vlenUtf8(), value);
    }

void writeDoubleAttribute(H5::H5Object &obj, const std::string &name, double value)
    {
    H5::Attribute attr = obj.createAttribute(name, H5::PredType::IEEE_F64LE,
                                             H5::DataSpace(H5S_SCALAR));
    attr.write(H5::PredType::NATIVE_DOUBLE, &value);
    }

std::string readStringAttribute(const H5::H5Object &obj, const std::string &name)
    {
    H5::Attribute attr = obj.openAttribute(name);
    if (attr.getTypeClass() != H5T_STRING)
        {
        throw std::runtime_error("attribute '" + name + "' is not a string");
        }
    std::string value;
    attr.read(attr.getStrType(), value);
    return value;
    }

std::string attributeToString(const H5::Attribute &attr)
    {
    const H5T_class_t cls = attr.getTypeClass();
    const bool scalar = attr.getSpace().getSimpleExtentNpoints() == 1;
    if (cls == H5T_STRING)
        {
        std::string value;
        attr.read(attr.getStrType(), value);
        return value;
        }
    if (cls == H5T_FLOAT && scalar)
        {
        double x;
        attr.read(H5::PredType::NATIVE_DOUBLE, &x);
        return pyRepr(x);
        }
    if (cls == H5T_INTEGER && scalar)
        {
        long long n;
        attr.read(H5::PredType::NATIVE_LLONG, &n);
        return std::to_string(n);
        }
    return "<" + std::to_string(attr.getSpace().getSimpleExtentNpoints()) + " values>";
    }

H5::DataSet writeTable(H5::Group &group, const std::string &name,
                       const std::vector<std::vector<double>> &rows)
    {
    const hsize_t nRows = rows.size();
    const hsize_t nCols = rows.empty() ? 0 : rows[0].size();
    std::vector<double> flat;
    flat.reserve(nRows * nCols);
    for (const auto &row : rows)
        {
        if (row.size() != nCols)
            {
            throw std::runtime_error("dataset '" + name + "': rows of different lengths");
            }
        flat.insert(flat.end(), row.begin(), row.end());
        }
    const hsize_t dims[2] = {nRows, nCols};
    H5::DataSet ds = group.createDataSet(name, H5::PredType::IEEE_F64LE, H5::DataSpace(2, dims));
    if (!flat.empty())
        {
        ds.write(flat.data(), H5::PredType::NATIVE_DOUBLE);
        }
    return ds;
    }

std::vector<std::vector<double>> readTable(const H5::DataSet &dataset)
    {
    H5::DataSpace space = dataset.getSpace();
    const int rank = space.getSimpleExtentNdims();
    std::vector<hsize_t> dims(std::max(rank, 1), 1);
    if (rank > 0)
        {
        space.getSimpleExtentDims(dims.data());
        }
    const hsize_t nRows = dims[0];
    hsize_t nCols = 1;
    for (int k = 1; k < rank; ++k)
        {
        nCols *= dims[k];
        }
    std::vector<double> flat(nRows * nCols);
    if (!flat.empty())
        {
        dataset.read(flat.data(), H5::PredType::NATIVE_DOUBLE);
        }
    std::vector<std::vector<double>> rows(nRows);
    for (hsize_t i = 0; i < nRows; ++i)
        {
        rows[i].assign(flat.begin() + i * nCols, flat.begin() + (i + 1) * nCols);
        }
    return rows;
    }

std::vector<std::string> attributeNames(const H5::H5Object &obj)
    {
    std::vector<std::string> names;
    auto collect = [](hid_t, const char *name, const H5A_info_t *, void *data) -> herr_t
        {
        static_cast<std::vector<std::string> *>(data)->push_back(name);
        return 0;
        };
    hsize_t idx = 0;
    H5Aiterate2(obj.getId(), H5_INDEX_NAME, H5_ITER_INC, &idx, collect, &names);
    return names;
    }

std::vector<std::string> memberNames(const H5::Group &group)
    {
    std::vector<std::string> names;
    const hsize_t n = group.getNumObjs();
    for (hsize_t i = 0; i < n; ++i)
        {
        names.push_back(group.getObjnameByIdx(i));
        }
    std::sort(names.begin(), names.end());
    return names;
    }
    }  // namespace fg
