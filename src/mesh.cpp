#include "mesh.h"
#include "meshUtils.h"
#include <utility>

using namespace Mesh;

/** \class BasicTri
BasicTri is a simplified version of the class Tri which contains only
the nodes, the surface region and the region.
It is used when checking the validity of the mesh is needed but
not all the triangle elements are needed.
*/
class BasicTri
    {
public:
    /** Constructor using the indices vector and the region. sRegion is false by default.
     * The nodes are sorted because it makes it easier to compare two BasicTri.
    */
    BasicTri(const std::vector<int> &inds, const int idReg, const bool surface = false)
        : nodesInd({inds[0], inds[1], inds[2]}), idRegion(idReg), isSurfaceElement(surface)
        { std::sort(nodesInd.begin(), nodesInd.end()); }

    /** Constructor used to copy a surface region triangle */
    explicit BasicTri(const Triangle::Tri & tri)
        : BasicTri(tri.ind, tri.idxPrm, true)
        {}

    /** Compare the nodes indice lexicographically */
    bool operator<(const BasicTri & o) const
        { return this->nodesInd < o.nodesInd; }

    /** The 3 node indices composing the triangle */
    std::array<int,3> nodesInd;

    /** index of the region of the triangle */
    int idRegion;

    /** If the triangle is in a surface region */
    bool isSurfaceElement;

    }; // class BasicTri

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
                Tetra::Tet const &te = tet[idxElem];
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
                                       [this, getter, &d, region](const int &idxElem)
                                       {
                                       Tetra::Tet const &te = tet[idxElem];
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
                  [&result, coord, whatToDo](Nodes::Node const &n)
                  {
                  double val = n.p(coord);
                  if (whatToDo(val, result)) result = val;
                  });
    return result;
    }

void mesh::indexReorder()
    {
    std::set<Triangle::Tri> sf;  // implicit use of operator< redefined in class Tri

    std::for_each(tet.begin(), tet.end(), [this,&sf](Tetra::Tet const &te)
                  {
                  const int ia = te.ind[0];
                  const int ib = te.ind[1];
                  const int ic = te.ind[2];
                  const int id = te.ind[3];

                  sf.insert(Triangle::Tri(node, 0, te.idxPrm, {ia, ic, ib} ));
                  sf.insert(Triangle::Tri(node, 0, te.idxPrm, {ib, ic, id} ));
                  sf.insert(Triangle::Tri(node, 0, te.idxPrm, {ia, id, ic} ));
                  sf.insert(Triangle::Tri(node, 0, te.idxPrm, {ia, ib, id} ));
                  });  // end for_each

    std::for_each(tri.begin(), tri.end(), [this, &sf](Triangle::Tri &fa)
                  {
                  int i0 = fa.ind[0], i1 = fa.ind[1], i2 = fa.ind[2];
                  auto it = sf.end();
                  for (int perm = 0; perm < 2; perm++)
                      {
                      for (int nrot = 0; nrot < 3; nrot++)
                          {
                          Triangle::Tri fc(node, 0, 0, {0, 0, 0} );
                          fc.ind[(0 + nrot) % 3] = i0;
                          fc.ind[(1 + nrot) % 3] = i1;
                          fc.ind[(2 + nrot) % 3] = i2;
                          it = sf.find(fc);
                          if (it != sf.end()) break;
                          }

                      if (it != sf.end())
                          {  // found
                          Eigen::Vector3d p0p1 = node[it->ind[1]].p - node[it->ind[0]].p;
                          Eigen::Vector3d p0p2 = node[it->ind[2]].p - node[it->ind[0]].p;
                                
                          // fa.Ms will have the magnitude of first arg of copysign, with second arg
                          // sign
                          fa.dMs = std::copysign( paramTetra[it->idxPrm].Ms,
                               p0p1.dot(p0p2.cross(fa.calc_norm())) );  // carefull, calc_norm
                                                                        // computes the normal to
                                                                        // the triangle before idx swap
                          }
                      std::swap(i1, i2);  // it seems from ref archive we do not want to swap
                                        // inner tri indices but local i1 and i2
                      fa.n = fa.calc_norm(); // update normal vector
                      }                   // end perm
                  });  // end for_each
    }

void mesh::checkTriangles()
    {
    std::vector<BasicTri> allTriCtnr;
    // Put all triangles into allTriCtnr after doing the conversion
    std::for_each(tet.begin(), tet.end(),       // For each tetrahedron
            [this, &allTriCtnr](const Tetra::Tet &tetrahedron)
            {
            for (int i = 0; i < Tetra::N; i++)    // For each 4 triangles
                {
                BasicTri curTri({tetrahedron.ind[i], tetrahedron.ind[(i+1) % Tetra::N],
                                        tetrahedron.ind[(i+2) % Tetra::N]},
                                 tetrahedron.idxPrm);
                allTriCtnr.push_back(curTri);
                }
            });
    std::for_each(tri.begin(), tri.end(), [this, &allTriCtnr](const Triangle::Tri &curTri)
            {   // For each surface element
            allTriCtnr.push_back(BasicTri(curTri));
            });

    if (allTriCtnr.empty())
        {
        std::cerr << "Error: not a single triangle is present in the mesh\n";
        exit(1);
        }

    std::sort(allTriCtnr.begin(), allTriCtnr.end());
    BasicTri *prevTri = &allTriCtnr[0];
    std::pair<int,int> pairIdCurVolRegs(-1, -1);
    int idCurSurfReg = -1;
    int nbTotTri = 0;
    int nbSurfReg = 0;

    for (size_t i = 0; i < allTriCtnr.size(); i++)
        {   // For each triangle
        if (allTriCtnr[i].nodesInd != prevTri->nodesInd)
            {   // If the triangle changes
            analyzeTriangle(nbTotTri, nbSurfReg, pairIdCurVolRegs, idCurSurfReg);
            nbTotTri = 0;
            nbSurfReg = 0;
            idCurSurfReg = -1;
            pairIdCurVolRegs.first = -1;
            pairIdCurVolRegs.second = -1;
            }

        if (allTriCtnr[i].isSurfaceElement)
            {
            idCurSurfReg = allTriCtnr[i].idRegion;
            nbSurfReg++;
            }
        else
            {
            if (pairIdCurVolRegs.first == -1)
                { pairIdCurVolRegs.first = allTriCtnr[i].idRegion; }
            else
                { pairIdCurVolRegs.second = allTriCtnr[i].idRegion; }
            }
        nbTotTri++;
        prevTri = &allTriCtnr[i];
        }
    analyzeTriangle(nbTotTri, nbSurfReg, pairIdCurVolRegs, idCurSurfReg);
    }

void mesh::analyzeTriangle(const int nbTotTri, const int nbSurfReg,
                           const std::pair<int,int> &pairIdVolRegs, const int idSurfReg)
    {
    if (nbSurfReg > 1)
        {
        std::cerr << "Error: wrong mesh generation. A triangle which belongs to "
                     "multiple surface regions has been found\n";
        exit(1);
        }
    else if (nbSurfReg == nbTotTri)
        {
        std::cout << "Error: wrong mesh generation. A triangle which only belongs "
                     "to 1 surface region has been found\n";
        exit(1);
        }
    else if (nbTotTri > 3 || (nbTotTri == 3 && nbSurfReg == 0))
        {
        std::cout << "Error: wrong mesh generation. Three "
                     "identical triangles have been found\n";
        if (nbSurfReg == 1)
            {
            std::cout << "In surface region " << idSurfReg << "\n";
            }
        exit(1);
        }
    else if (nbTotTri == 3 && nbSurfReg == 1 && pairIdVolRegs.first == pairIdVolRegs.second)
        {
        std::cout << "Error: wrong mesh generation. An internal triangle belonging "
                     "to the surface region " << idSurfReg << " has been found\n";
        exit(1);
        }
    }

void mesh::sortNodes(Nodes::index long_axis)
    {
    // Sort the nodes along this axis, indirectly through an array of indices.
    std::vector<int> permutation(node.size());
    std::iota(permutation.begin(), permutation.end(), 0);
    std::sort(permutation.begin(), permutation.end(),
              [this, long_axis](int a, int b)
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
    std::for_each(triIndices.begin(),triIndices.end(),[this,&S](int idx)
                  { S += tri[idx].calc_surf(); });
    return S;
    }

void mesh::setExtSpaceField(Settings &s /**< [in] */)
    { // see here for ref code
// /data/jc/st_feellgood_2024/src_Tube_scalfmm_zhang_ec_mu_oersted_spinHall_thiele_dyn20240320_dev
    extSpaceField.resize(magTet.size());
    int k(0);
    std::for_each(magTet.begin(), magTet.end(), [this,&s,&k](const int idx)
        {
        Tetra::Tet const &t = tet[idx];
        Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> pg;// gauss points
        t.getPtGauss(pg);
        for(int i=0;i<Tetra::NPI;i++)
            { extSpaceField[k].col(i) = s.getFieldSpace(pg.col(i)); }
        k++;
        });
    }
