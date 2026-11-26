#ifndef FG3D_MESHDATA_H
#define FG3D_MESHDATA_H

/** \file meshData.h
 * Mesh and magnetization data rendered by fg-3D, and the surface charge computation.
 */

#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "gmshMesh.h"

using Vec3 = fg::Vec3;

/** A mesh with the set of the nodes of its surface. */
struct MeshData : fg::Mesh
    {
    std::set<std::size_t> surfaceNodes;  //!< nodes belonging to at least one triangle
    };

/** Magnetization read from a .sol file: node index -> (mx, my, mz). */
using Magnetization = std::map<std::size_t, Vec3>;

/** Load the triangles and tetrahedrons of a gmsh mesh file, any version, and find its surface
 * nodes. Throws std::runtime_error if the file cannot be read or contains no triangles. A
 * fg::GmshSession must exist. Returns a log message. */
std::string loadMesh(const std::string &fileName, MeshData &mesh);

/** Load the magnetization of a feeLLGood .sol file: node index in the first column, mx, my,
 * mz in the next three. Throws std::runtime_error if the file cannot be read or has less than
 * four columns. */
Magnetization loadSol(const std::string &fileName);

/** Normals of the triangles, of length twice their area, pointing out of the tetrahedron each
 * triangle bounds, as feeLLGood orients its surface elements. A triangle bounding two
 * tetrahedrons, at an interface between two volume regions, points out of the first one
 * found. A triangle bounding no tetrahedron keeps the orientation of its nodes. */
std::vector<Vec3> outwardNormals(const MeshData &mesh);

/** Surface charge of each triangle: dot product of the normalized average magnetization of
 * its three nodes with its unit outward normal. normals[i] is the outward normal of
 * triangles[i]. Returns an empty vector if a triangle node has no magnetization. */
std::vector<double> computeSurfaceCharges(const MeshData &mesh, const Magnetization &m,
                                          const std::vector<Vec3> &normals);

#endif
