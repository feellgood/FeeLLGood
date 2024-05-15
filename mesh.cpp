#include "mesh.h"

using namespace Mesh;

void mesh::infos(void) const
    {
    std::cout << "mesh:\n";
    std::cout << "  nodes:              " << getNbNodes() << '\n';
    std::cout << "  faces:              " << fac.size() << '\n';
    std::cout << "  tetraedrons:        " << tet.size() << '\n';
    std::cout << "  total volume:       " << vol << '\n';
    }

double mesh::updateNodes(Eigen::Ref<Eigen::VectorXd> X, const double dt)
    {
    double v2max = 0.0;
    const unsigned int NOD = node.size();

    for (unsigned int i = 0; i < NOD; i++)
        {
        const double & vp = X(i) * gamma0;
        const double & vq = X(NOD + i) * gamma0;
        double v2 = vp*vp + vq*vq;
        if (v2 > v2max)
            {
            v2max = v2;
            }
        node[i].make_evol(vp, vq, dt);
        }

    return sqrt(v2max);
    }

void mesh::buildInitGuess(Eigen::Ref<Eigen::VectorXd> G) const
    {
    const unsigned int NOD = node.size();

    for (unsigned int i = 0; i < NOD; i++)
        {
        Eigen::Vector3d const& n_v = getNode_v(i);//node[i].v;
            
        G(i) = n_v.dot(node[i].ep) / gamma0;
        G(NOD + i) = n_v.dot(node[i].eq) / gamma0;
        }
    }

double mesh::avg(std::function<double(Nodes::Node, Nodes::index)> getter /**< [in] */,
                 Nodes::index d /**< [in] */) const
    {
    double sum = std::transform_reduce(EXEC_POL, tet.begin(), tet.end(), 0.0, std::plus{},
                                       [getter, &d](Tetra::Tet const &te)
                                       {
                                       Eigen::Matrix<double,Tetra::NPI,1> val;
                                       te.interpolation(getter, d, val);
                                       return te.weight.dot(val);
                                       });
    return sum / vol;
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

void mesh::indexReorder(std::vector<Tetra::prm> const &prmTetra)
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

    std::for_each(fac.begin(), fac.end(), [this, &prmTetra, &sf](Facette::Fac &fa)
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
                          fa.dMs = std::copysign( prmTetra[it->idxPrm].J/mu0,
                               p0p1.dot(p0p2.cross(fa.calc_norm())) );  // carefull, calc_norm computes the normal to the face before idx swap
                          }
                      std::swap(i1, i2);  // it seems from ref archive we do not want to swap
                                        // inner fac indices but local i1 and i2
                      fa.n = fa.calc_norm(); // update normal vector
                      }                   // end perm
                  });  // end for_each
    }

void mesh::sortNodes()
    {
    // Find the longest axis of the sample.
    Nodes::index long_axis;
    if (l.x() > l.y())
        {
        if (l.x() > l.z())
            long_axis = Nodes::IDX_X;
        else
            long_axis = Nodes::IDX_Z;
        }
    else
        {
        if (l.y() > l.z())
            long_axis = Nodes::IDX_Y;
        else
            long_axis = Nodes::IDX_Z;
        }

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
    std::for_each(s.begin(), s.end(),
                  [this](Mesh::Surf &surface)
                  {
                      std::vector<Mesh::Triangle> &elements = surface.elem;
                      std::for_each(elements.begin(), elements.end(),
                                    [this](Mesh::Triangle &tri)
                                    {
                                        tri.ind[0] = node_index[tri.ind[0]];
                                        tri.ind[1] = node_index[tri.ind[1]];
                                        tri.ind[2] = node_index[tri.ind[2]];
                                    });
                  });
    }
