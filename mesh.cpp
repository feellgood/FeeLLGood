#include "mesh.h"
#include "meshUtils.h"
#include <utility>

using namespace Mesh;

void mesh::infos(void) const
    {
    std::cout << "mesh:\n";
    std::cout << "  nodes:              " << getNbNodes() << '\n';
    std::cout << "  faces:              " << fac.size() << '\n';
    std::cout << "  tetraedrons:        " << tet.size() << '\n';
    std::cout << "  total volume:       " << vol << '\n';
    }

double mesh::avg(std::function<double(Nodes::Node, Nodes::index)> getter /**< [in] */,
                 Nodes::index d /**< [in] */,
                 int region /**< region index, or -1 for all regions */) const
    {
    double sum = std::transform_reduce(EXEC_POL, tet.begin(), tet.end(), 0.0,
                                       std::plus<>(),
                                       [getter, &d, region](Tetra::Tet const &te)
                                       {
                                       if (te.idxPrm != region && region != -1)
                                           return 0.0;
                                       Eigen::Matrix<double,Tetra::NPI,1> val;
                                       te.interpolation(getter, d, val);
                                       return te.weight.dot(val);
                                       });
    double volume = (region == -1) ? vol : paramTetra[region].volume;
    return sum / volume;
    }

double mesh::doOnNodes(const double init_val, const Nodes::index coord,
                     std::function<bool(double, double)> whatToDo) const
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
    std::set<Facette::Fac> sf;  // implicit use of operator< redefined in class Fac

    std::for_each(tet.begin(), tet.end(), [this,&sf](Tetra::Tet const &te)
                  {
                  const int ia = te.ind[0];
                  const int ib = te.ind[1];
                  const int ic = te.ind[2];
                  const int id = te.ind[3];

                  sf.insert(Facette::Fac(node, 0, te.idxPrm, {ia, ic, ib} ));
                  sf.insert(Facette::Fac(node, 0, te.idxPrm, {ib, ic, id} ));
                  sf.insert(Facette::Fac(node, 0, te.idxPrm, {ia, id, ic} ));
                  sf.insert(Facette::Fac(node, 0, te.idxPrm, {ia, ib, id} ));
                  });  // end for_each

    std::for_each(fac.begin(), fac.end(), [this, &sf](Facette::Fac &fa)
                  {
                  int i0 = fa.ind[0], i1 = fa.ind[1], i2 = fa.ind[2];
                  std::set<Facette::Fac>::iterator it = sf.end();
                  for (int perm = 0; perm < 2; perm++)
                      {
                      for (int nrot = 0; nrot < 3; nrot++)
                          {
                          Facette::Fac fc(node, 0, 0, {0, 0, 0} );
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
                                
                          // fa.Ms will have the magnitude of first arg of copysign, with second arg sign
                          fa.dMs = std::copysign( paramTetra[it->idxPrm].J/mu0,
                               p0p1.dot(p0p2.cross(fa.calc_norm())) );  // carefull, calc_norm computes the normal to the face before idx swap
                          }
                      std::swap(i1, i2);  // it seems from ref archive we do not want to swap
                                        // inner fac indices but local i1 and i2
                      fa.n = fa.calc_norm(); // update normal vector
                      }                   // end perm
                  });  // end for_each
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
    std::for_each(fac.begin(), fac.end(),
                  [this](Facette::Fac &facette)
                  {
                      facette.ind[0] = node_index[facette.ind[0]];
                      facette.ind[1] = node_index[facette.ind[1]];
                      facette.ind[2] = node_index[facette.ind[2]];
                  });
    }

void mesh::build_lvd(std::vector<int> &lvd)
    {
    const int NOD = node.size();
    std::vector<int> ld(NOD);// list of Dirichlet nodes where to solve LLG
    std::vector<int> all(NOD);
    std::vector<int> dofs; // list of the nodes where LLG does apply
    std::iota(all.begin(),all.end(),0); // all = { 0,1,...,NOD-1}
    std::for_each(tet.begin(),tet.end(),[this,&dofs](const Tetra::Tet &t)
                  {
                  if (isMagnetic(t))
                      {
                      dofs.push_back(t.ind[0]);
                      dofs.push_back(t.ind[1]);
                      dofs.push_back(t.ind[2]);
                      dofs.push_back(t.ind[3]);
                      }
                  });

    suppress_copies<int>(dofs);
    auto it = std::set_difference(all.begin(),all.end(),dofs.begin(),dofs.end(),ld.begin());
    auto nb_coeffs = it - ld.begin();
    lvd.resize(2*nb_coeffs);
    for(int i=0;i<nb_coeffs;i++)
        {
        lvd[2*i] = ld[i];
        lvd[2*i+1] = ld[i] + NOD;
        }
    }

double mesh::surface(std::vector<int> &facIndices)
    {
    double S(0);
    std::for_each(facIndices.begin(),facIndices.end(),[this,&S](int idx)
                  { S += fac[idx].calc_surf(); });
    return S;
    }

void mesh::buildBoundaryConditions(std::vector<Facette::prm> &paramFacette,
                                   Mesh::boundaryCondition<Eigen::Vector3d> &BC_spin)
    {
    for(unsigned int i=0;i<fac.size();i++)
        {
        Eigen::Vector3d Pu = paramFacette[fac[i].idxPrm].Pu;
        bool Pu_finite = std::isfinite(Pu.norm());
        double V = paramFacette[fac[i].idxPrm].V;
        bool V_finite = std::isfinite(V);
        double J = paramFacette[fac[i].idxPrm].J;
        bool J_finite = std::isfinite(J);
        if (!Pu_finite && !V_finite && !J_finite)
            { /* magnetic surface */ }
        else if (!Pu_finite && V_finite && !J_finite)
            {
            // Pu is undefined, f is not part of a surface boundary condition for spin acc
            //this is an isopotential V for boundary condition for electrostat
            
            }
        else if (!Pu_finite && !V_finite && J_finite)
            {
            // Pu is undefined, f is not part of a surface boundary condition for spin acc
            //this is an isopotential J for boundary condition for electrostat

            }
        else if (Pu_finite && !V_finite && J_finite)
            {
            // we have both J and Pu on the same facette, we can compute Q for boundary
            // condition of spin accumulation
            // warning: we mix the the two(or more) surfaces where to define Q for spin acc BC's
            BC_spin.facIdx.push_back(i);
            BC_spin.value = J*Pu;
            }
        else if (Pu_finite && V_finite && !J_finite)
            {
            // Pu and V are finite, we have to compute J = current density normal to f
            // to define Q for boundary conditions

            }
        else if (V_finite && J_finite)
            {
            std::cout << "Fatal Error: Boundary conditions overdefined: J and V set on the same surface.\n";
            exit(1);
            }
        }
    }
