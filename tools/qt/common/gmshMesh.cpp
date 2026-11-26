#include "gmshMesh.h"

#include <fstream>
#include <stdexcept>

#include <gmsh.h>

namespace fg
    {
namespace
    {
constexpr int GMSH_TRIANGLE = 2;
constexpr int GMSH_TETRAHEDRON = 4;
    }  // namespace

GmshSession::GmshSession()
    {
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 0);
    }

GmshSession::~GmshSession() { gmsh::finalize(); }

std::size_t nodeIndex(std::size_t tag, std::size_t nbNodes)
    {
    if (tag == 0 || tag > nbNodes)
        {
        throw std::runtime_error("node tag " + std::to_string(tag) + " out of range");
        }
    return tag - 1;
    }

Mesh readMesh(const std::string &fileName, bool require3D)
    {
    if (!std::ifstream(fileName))
        {
        throw std::runtime_error("the file '" + fileName + "' does not exist.");
        }

    Mesh m;
    try
        {
        gmsh::clear();
        gmsh::open(fileName);

        if (require3D)
            {
            if (gmsh::model::getDimension() != 3)
                {
                throw std::runtime_error("the model is not 3D");
                }
            std::vector<int> types;
            gmsh::model::mesh::getElementTypes(types, 3, -1);
            for (int t : types)
                {
                if (t != GMSH_TETRAHEDRON)
                    {
                    throw std::runtime_error("only first order tetrahedrons are supported in "
                                             "volumes");
                    }
                }
            }

        std::vector<std::size_t> nodeTags;
        std::vector<double> coords, parametricCoords;
        gmsh::model::mesh::getNodes(nodeTags, coords, parametricCoords);
        const std::size_t nbNodes = nodeTags.size();
        m.points.resize(nbNodes);
        for (std::size_t i = 0; i < nbNodes; ++i)
            {
            m.points[nodeIndex(nodeTags[i], nbNodes)] = {coords[3 * i], coords[3 * i + 1],
                                                         coords[3 * i + 2]};
            }

        std::vector<std::size_t> elemTags, elemNodeTags;
        gmsh::model::mesh::getElementsByType(GMSH_TRIANGLE, elemTags, elemNodeTags);
        for (std::size_t i = 0; i + 2 < elemNodeTags.size(); i += 3)
            {
            m.triangles.push_back({nodeIndex(elemNodeTags[i], nbNodes),
                                   nodeIndex(elemNodeTags[i + 1], nbNodes),
                                   nodeIndex(elemNodeTags[i + 2], nbNodes)});
            }

        gmsh::model::mesh::getElementsByType(GMSH_TETRAHEDRON, elemTags, elemNodeTags);
        for (std::size_t i = 0; i + 3 < elemNodeTags.size(); i += 4)
            {
            m.tetrahedrons.push_back({nodeIndex(elemNodeTags[i], nbNodes),
                                      nodeIndex(elemNodeTags[i + 1], nbNodes),
                                      nodeIndex(elemNodeTags[i + 2], nbNodes),
                                      nodeIndex(elemNodeTags[i + 3], nbNodes)});
            }
        }
    catch (const std::exception &e)
        {
        throw std::runtime_error("error while loading mesh '" + fileName + "': " + e.what());
        }
    catch (...)
        {  // gmsh may throw non std exceptions
        throw std::runtime_error("error while loading mesh '" + fileName + "'");
        }
    return m;
    }
    }  // namespace fg
