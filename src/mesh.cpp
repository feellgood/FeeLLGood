#include "mesh.h"
#include "meshUtils.h"
#include <utility>

using namespace Mesh;

void mesh::infos(void) const
    {
    std::cout << "mesh:\n";
    std::cout << "  nodes:              " << getNbNodes() << '\n';
    std::cout << "  triangles:          " << tri.size() << '\n';
    std::cout << "  tetraedrons:        " << tet.size() << '\n';
    std::cout << "  total volume:       " << vol << '\n';
    }

double mesh::thiele(const int region /** region index, or -1 for all magnetic regions */) const
    {
    double sum = std::transform_reduce(EXEC_POL,magTet.begin(),magTet.end(),0.0,std::plus<>(),
            [this, region](const int idxElem)
                {//the lambda computes sq(sum(grad(u))) on the tetrahedron indexed by idxElem
                double val(0);
                const Tetra::Tet &te = tet[idxElem];
                if((te.idxPrm == region) || (region == -1))
                    {
                    Eigen::Matrix<double,Nodes::DIM,Tetra::N> u_nod;
                    for (int i = 0; i< Tetra::N; i++)
                        u_nod.col(i) = getNode_u(te.ind[i]);
                
                    Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> du_dz = u_nod * (te.da.col(Nodes::IDX_Z)).replicate(1,Tetra::NPI);
                
                    for(int npi=0;npi<Tetra::NPI;npi++)
                        val += du_dz.col(npi).dot(du_dz.col(npi)) * te.weight[npi];
                    }
                return val;
                });
    double volume(totalMagVol);
    if (region != -1)
        { volume = paramTetra[region].volume; }
    double cross_section = volume/l.z();//average cross-section
    return 2.0*cross_section/sum;
    }


double mesh::avg(const std::function<double(Nodes::Node, Nodes::index)>& getter /**< [in] */,
                 Nodes::index d /**< [in] */,
                 int region /**< region index, or -1 for all magnetic regions */) const
    {
    double sum = std::transform_reduce(EXEC_POL, magTet.begin(), magTet.end(), 0.0,
                                       std::plus<>(),
                                       [this, &getter, &d, region](const int &idxElem)
                                       {
                                       const Tetra::Tet &te = tet[idxElem];
                                       if (te.idxPrm != region && region != -1)
                                           return 0.0;
                                       Eigen::Matrix<double,Tetra::NPI,1> val;
                                       te.interpolation(getter, d, val);
                                       return te.weight.dot(val);
                                       });
    double volume = (region == -1) ? totalMagVol : paramTetra[region].volume;
    return sum / volume;
    }

double mesh::doOnNodes(const double init_val, const Nodes::index coord,
                     const std::function<bool(double, double)>& whatToDo) const
    {
    double result(init_val);
    std::for_each(node.begin(), node.end(),
                  [&result, coord, whatToDo](const Nodes::Node &n)
                  {
                  double val = n.p(coord);
                  if (whatToDo(val, result)) result = val;
                  });
    return result;
    }

bool mesh::controlTriangles()
    {
    if (tet.empty())
        {
        std::cerr << "Error: not a single tetrahedron is present in the mesh\n";
        return false;
        }

    // Build a sorted list of all the triangles in the mesh: both the faces of the tetrahedrons
    // and the (triangular) surface elements.
    std::vector<BasicTri> allTriCtnr;

    // First, add the faces of the tetrahedrons.
    for (const Tetra::Tet &tetrahedron : tet)
        {  // Add the faces of this tetrahedron, oriented outwards.
        const std::array<int,Tetra::N>& ind = tetrahedron.ind;
        const int region = tetrahedron.idxPrm;
        allTriCtnr.push_back(BasicTri({ind[0], ind[2], ind[1]}, region));
        allTriCtnr.push_back(BasicTri({ind[1], ind[2], ind[3]}, region));
        allTriCtnr.push_back(BasicTri({ind[2], ind[0], ind[3]}, region));
        allTriCtnr.push_back(BasicTri({ind[3], ind[0], ind[1]}, region));
        }

    // Then add the surface elements.
    for (Triangle::Tri &curTri : tri)
        { allTriCtnr.push_back(BasicTri(curTri)); }

    // Then sort: triangles that are topologically the same (same set of vertices) will end up
    // bunched together.
    std::sort(allTriCtnr.begin(), allTriCtnr.end());

    // Update dMs on each triangular element.
    for (Triangle::Tri &curTri : tri)
        {
        BasicTri triangle(curTri);
        // Find the tetrahedron faces that match this surface element.
        auto face = std::lower_bound(allTriCtnr.begin(), allTriCtnr.end(), triangle);
        for (; face != allTriCtnr.end() && face->nodesInd == triangle.nodesInd; face++)
            {
            if (!face->isSurfaceElement)
                {
                double Ms = paramTetra[face->idRegion].Ms;
                bool isFlipped = triangle.isFlipped ^ face->isFlipped;
                curTri.dMs += isFlipped ? -Ms : Ms;
                }
            }
        }

    // Now walk the list of all the triangles in search for meshing errors.
    std::vector<size_t> indNodTriToAdd;
    std::vector<int> idxPrmTriToAdd;

    size_t inf = 0;
    // For each triangle plus one time after
    for (size_t sup = 1; sup <= allTriCtnr.size(); sup++)
        {
        // If the triangle is different
        if (sup == allTriCtnr.size() || allTriCtnr[inf].nodesInd != allTriCtnr[sup].nodesInd)
            {
            if (!diffTriHandler(inf, sup, allTriCtnr, indNodTriToAdd, idxPrmTriToAdd))
                { return false; }

            inf = sup;
            }
        }

    // Creates the missings elements in tri
    for (size_t i = 0; i + 2 < indNodTriToAdd.size(); i = i + 3)
        {
        int curIdxPrm = idxPrmTriToAdd[i/3];
        int i0 = indNodTriToAdd[i];
        int i1 = indNodTriToAdd[i+1];
        int i2 = indNodTriToAdd[i+2];
        tri.push_back(Triangle::Tri(node, curIdxPrm, {i0, i1, i2}));
        }
    return true;
    }

bool Mesh::mesh::diffTriHandler(const size_t inf, const size_t sup, const std::vector<BasicTri> &allTriCtnr,
                                std::vector<size_t> &indNodTriToAdd, std::vector<int> &idxPrmTriToAdd)
    {
    int idCurSurfReg = -1;
    int nbTetraFaces = 0;
    int nbSurfTri = 0;
    std::pair<int,int> pairIdCurVolRegs(-1, -1);

    // Gets the triangle's stats
    for (size_t i = inf; i < sup; i++)
        {
        if (allTriCtnr[i].isSurfaceElement)
            {
            idCurSurfReg = allTriCtnr[inf].idRegion;
            nbSurfTri++;
            }
        else
            {
            nbTetraFaces++;
            if (pairIdCurVolRegs.first == -1)
                { pairIdCurVolRegs.first = allTriCtnr[inf].idRegion; }
            else
                { pairIdCurVolRegs.second = allTriCtnr[inf].idRegion; }
            }
        }

    if (nbSurfTri == 0 && (nbTetraFaces == 1 ||
            (nbTetraFaces == 2 && pairIdCurVolRegs.first != pairIdCurVolRegs.second)))
        {   // If it is an interface or a surface triangle not in tri
        indNodTriToAdd.push_back(allTriCtnr[inf].nodesInd[0]);
        indNodTriToAdd.push_back(allTriCtnr[inf].nodesInd[1]);
        indNodTriToAdd.push_back(allTriCtnr[inf].nodesInd[2]);
        idxPrmTriToAdd.push_back(allTriCtnr[inf].idRegion);
        }

    if (findNoErrInTriangle(nbSurfTri, nbTetraFaces, pairIdCurVolRegs, idCurSurfReg))
        { return true; }
    else
        { return false; }
    }

bool mesh::findNoErrInTriangle(const int nbSurfTri, const int nbTetraFaces,
                           const std::pair<int,int> pairIdVolRegs, const int idSurfReg)
    {
    if (nbSurfTri > 1)
        {
        std::cerr << "Error: bad mesh. " << nbSurfTri << " instances of"
                     " the same surface triangle have been found\n";
        return false;
        }
    else if (nbTetraFaces == 0)
        {
        std::cerr << "Error: bad mesh. A triangle which belongs to 1 "
                     "surface region and no tetrahedron has been found\n";
        return false;
        }
    else if (nbTetraFaces > 2)
        {
        if (nbSurfTri == 1)
            {
            std::cerr << "Error: bad mesh. A triangular face shared by " << nbTetraFaces
                      << " tetrahedrons has been found in the surface region " << idSurfReg << "\n";
            }
        else
            {
            std::cerr << "Error: bad mesh. A triangular face shared by "
                      << nbTetraFaces << " tetrahedrons has been found\n";
            }
        return false;
        }
    else if (nbTetraFaces == 2 && nbSurfTri == 1 && pairIdVolRegs.first == pairIdVolRegs.second)
        {
        std::cerr << "Error: bad mesh. An internal triangle has been "
                     "found in the surface region " << idSurfReg << "\n";
        return false;
        }
    else
        {
        return true;
        }
    }

void mesh::sortNodes(Nodes::index long_axis)
    {
    // Sort the nodes along this axis, indirectly through an array of indices.
    std::vector<int> permutation(node.size());
    std::iota(permutation.begin(), permutation.end(), 0);
    std::sort(permutation.begin(), permutation.end(),
              [this, long_axis](const int a, const int b)
              { return node[a].p(long_axis) < node[b].p(long_axis); });
    node_index.resize(node.size());
    for (size_t i = 0; i < node.size(); i++)
        node_index[permutation[i]] = i;

    // Actually sort the array of nodes.
    std::vector<Nodes::Node> node_copy(node);
    for (size_t i = 0; i < node.size(); i++)
        node[i] = node_copy[permutation[i]];

    // Update the indices stored in the elements.
    std::for_each(tet.begin(), tet.end(),
                  [this](Tetra::Tet &tetrahedron)
                  {
                      tetrahedron.ind[0] = node_index[tetrahedron.ind[0]];
                      tetrahedron.ind[1] = node_index[tetrahedron.ind[1]];
                      tetrahedron.ind[2] = node_index[tetrahedron.ind[2]];
                      tetrahedron.ind[3] = node_index[tetrahedron.ind[3]];
                  });
    std::for_each(tri.begin(), tri.end(),
                  [this](Triangle::Tri &triangle)
                  {
                      triangle.ind[0] = node_index[triangle.ind[0]];
                      triangle.ind[1] = node_index[triangle.ind[1]];
                      triangle.ind[2] = node_index[triangle.ind[2]];
                  });
    }

double mesh::surface(std::vector<int> &triIndices) const
    {
    double S(0);
    std::for_each(triIndices.begin(),triIndices.end(),[this,&S](const int idx)
                  { S += tri[idx].calc_surf(); });
    return S;
    }

void mesh::setExtSpaceField(Settings &s /**< [in] */)
    {
    extSpaceField.resize(magTet.size());
    int k(0);
    std::for_each(magTet.begin(), magTet.end(), [this,&s,&k](const int idx)
        {
        const Tetra::Tet &t = tet[idx];
        Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> pg = t.getPtGauss();
        for(int i=0;i<Tetra::NPI;i++)
            { extSpaceField[k].col(i) = s.getFieldSpace(pg.col(i)); }
        k++;
        });
    }
