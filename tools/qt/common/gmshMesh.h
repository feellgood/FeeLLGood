#ifndef FG_GMSHMESH_H
#define FG_GMSHMESH_H

/** \file gmshMesh.h
 * Reading of gmsh mesh files through the gmsh API, shared by the feeLLGood tools.
 */

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace fg
    {
/** A point or a vector in 3D. */
using Vec3 = std::array<double, 3>;

/** \class GmshSession
 * Initializes gmsh, silenced, for the lifetime of the object, and finalizes it on
 * destruction. Exactly one must exist while gmsh functions are called.
 */
class GmshSession
    {
public:
    /** Initialize gmsh and silence its terminal output. */
    GmshSession();

    /** Finalize gmsh. */
    ~GmshSession();

    GmshSession(const GmshSession &) = delete;
    GmshSession &operator=(const GmshSession &) = delete;
    };

/** Nodes, triangles and tetrahedrons of a mesh. Node indices are gmsh node tags minus one,
 * the convention feeLLGood uses for the node indices written in .sol files. */
struct Mesh
    {
    std::vector<Vec3> points;                              //!< node coordinates
    std::vector<std::array<std::size_t, 3>> triangles;     //!< node indices of triangles
    std::vector<std::array<std::size_t, 4>> tetrahedrons;  //!< node indices of tetrahedrons
    };

/** Read the nodes, the first order triangles and the first order tetrahedrons of a gmsh mesh
 * file, any version. Other element types are ignored. If require3D is true, the model must be
 * 3D and its volumes must contain tetrahedrons only, as feeLLGood requires. The model stays
 * loaded in gmsh until the next call, so that the caller can query it further. A GmshSession
 * must exist. Throws std::runtime_error if the file cannot be read or does not meet these
 * requirements. */
Mesh readMesh(const std::string &fileName, bool require3D = false);

/** Convert a gmsh node tag to a node index (tag minus one). Throws std::runtime_error if the
 * tag is not in [1, nbNodes]. */
std::size_t nodeIndex(std::size_t tag, std::size_t nbNodes);
    }  // namespace fg

#endif
