#include "meshData.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "textTable.h"

namespace
    {
/** Euclidean norm of v. */
double norm(const Vec3 &v) { return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); }
    }  // namespace

std::string loadMesh(const std::string &fileName, MeshData &mesh)
    {
    MeshData m;
    static_cast<fg::Mesh &>(m) = fg::readMesh(fileName);
    if (m.triangles.empty())
        {
        throw std::runtime_error("no triangles found in mesh '" + fileName + "'.");
        }
    for (const auto &tri : m.triangles)
        {
        m.surfaceNodes.insert(tri.begin(), tri.end());
        }
    mesh = std::move(m);
    return "✓ mesh file '" + fileName + "' loaded, " + std::to_string(mesh.triangles.size())
           + " triangles, " + std::to_string(mesh.tetrahedrons.size()) + " tetrahedrons";
    }

std::vector<Vec3> outwardNormals(const MeshData &mesh)
    {
    // sorted node indices of a tetrahedron face -> node opposite to it
    std::map<std::array<std::size_t, 3>, std::size_t> apex;
    for (const auto &tet : mesh.tetrahedrons)
        {
        for (int k = 0; k < 4; ++k)
            {
            std::array<std::size_t, 3> face;
            for (int j = 0, f = 0; j < 4; ++j)
                {
                if (j != k)
                    {
                    face[f++] = tet[j];
                    }
                }
            std::sort(face.begin(), face.end());
            apex.emplace(face, tet[k]);
            }
        }

    std::vector<Vec3> normals;
    normals.reserve(mesh.triangles.size());
    for (const auto &tri : mesh.triangles)
        {
        const Vec3 &p0 = mesh.points[tri[0]], &p1 = mesh.points[tri[1]],
                   &p2 = mesh.points[tri[2]];
        const Vec3 u = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
        const Vec3 v = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
        Vec3 n = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                  u[0] * v[1] - u[1] * v[0]};
        std::array<std::size_t, 3> key = tri;
        std::sort(key.begin(), key.end());
        const auto it = apex.find(key);
        if (it != apex.end())
            {
            const Vec3 &a = mesh.points[it->second];
            const double side = n[0] * (p0[0] - a[0]) + n[1] * (p0[1] - a[1])
                                + n[2] * (p0[2] - a[2]);
            if (side < 0)
                {
                n = {-n[0], -n[1], -n[2]};
                }
            }
        normals.push_back(n);
        }
    return normals;
    }

Magnetization loadSol(const std::string &fileName)
    {
    const fg::TextTable table = fg::readTextTable(fileName);
    Magnetization mag;
    for (const auto &row : table.rows)
        {
        if (row.size() < 4 || row[0] < 0)
            {
            throw std::runtime_error("'" + fileName + "' is not a magnetization .sol file");
            }
        mag[static_cast<std::size_t>(row[0])] = {row[1], row[2], row[3]};
        }
    return mag;
    }

std::vector<double> computeSurfaceCharges(const MeshData &mesh, const Magnetization &m,
                                          const std::vector<Vec3> &normals)
    {
    std::vector<double> charges;
    if (mesh.triangles.size() != normals.size())
        {
        return charges;
        }
    charges.reserve(normals.size());
    for (std::size_t t = 0; t < normals.size(); ++t)
        {
        Vec3 avg = {0, 0, 0};
        for (std::size_t node : mesh.triangles[t])
            {
            auto it = m.find(node);
            if (it == m.end())
                {
                return {};
                }
            for (int k = 0; k < 3; ++k)
                {
                avg[k] += it->second[k];
                }
            }
        const double na = norm(avg);
        const double nn = norm(normals[t]);
        double dot = 0;
        for (int k = 0; k < 3; ++k)
            {
            dot += avg[k] * normals[t][k];
            }
        charges.push_back((na > 0 && nn > 0) ? dot / (na * nn) : 0.0);
        }
    return charges;
    }
