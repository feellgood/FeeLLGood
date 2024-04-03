#ifndef mesh_h
#define mesh_h

/** \file mesh.h
\brief class mesh, readMesh is expecting a mesh file in gmsh text format 2.2, with first order tetraedrons and triangular facettes.
*/

#include <algorithm>
#include <execution>
#include <map>
#include <numeric>
#include <set>

#include "facette.h"
#include "node.h"
#include "surface.h"
#include "tetra.h"

#include "feellgoodSettings.h"

namespace Mesh
    {

/** \class mesh
class for storing the mesh, including mesh geometry values, containers for the nodes, triangular
faces and tetrahedrons. nodes data are private. They are accessible only through getter and setter.
*/
class mesh
    {
public:
    /** constructor : read mesh file, reorder indices and computes some values related to the mesh :
     center and length along coordinates and diameter = max(l(x|y|z)), volume and surface */
    inline mesh(Settings const &mySets /**< [in] */)
        {
        readMesh(mySets);
        indexReorder(mySets.paramTetra);  // reordering of index nodes for facette orientation and defines fac::dMs
        if (mySets.verbose)
            { std::cout << "  reindexed\n"; }

        double xmin = minNodes(Nodes::IDX_X);
        double xmax = maxNodes(Nodes::IDX_X);

        double ymin = minNodes(Nodes::IDX_Y);
        double ymax = maxNodes(Nodes::IDX_Y);

        double zmin = minNodes(Nodes::IDX_Z);
        double zmax = maxNodes(Nodes::IDX_Z);

        l = Eigen::Vector3d(xmax - xmin, ymax - ymin, zmax - zmin);
        diam = l.maxCoeff();
        c = Eigen::Vector3d(0.5 * (xmax + xmin), 0.5 * (ymax + ymin), 0.5 * (zmax + zmin));

        sortNodes();

        vol = std::transform_reduce(EXEC_POL, tet.begin(), tet.end(), 0.0, std::plus{},
                                    [](Tetra::Tet const &te) { return te.calc_vol(); });

        surf = std::transform_reduce(EXEC_POL, fac.begin(), fac.end(), 0.0, std::plus{},
                                     [](Facette::Fac const &fa) { return fa.surf; });
        }

    /** return number of nodes  */
    inline int getNbNodes(void) const { return node.size(); }

    /** return number of triangular fac */
    inline int getNbFacs(void) const { return fac.size(); }

    /** return number of tetrahedrons */
    inline int getNbTets(void) const { return tet.size(); }

    /** getter : return node.p */
    inline const Eigen::Vector3d getNode_p(const int i) const { return node[i].p; }
    
    /** getter : return node.u */
    inline const Eigen::Vector3d getNode_u(const int i) const { return node[i].get_u(Nodes::NEXT); }

    /** getter : return node.v */
    inline const Eigen::Vector3d getNode_v(const int i) const { return node[i].get_v(Nodes::NEXT); }

    /** setter for u0 */
    inline void set_node_u0(const int i, Eigen::Vector3d const &val)
        {
        node[i].d[Nodes::CURRENT].u = val;
        }

    /** fix to zero node[i].v */
    inline void set_node_zero_v(const int i) { node[i].d[Nodes::NEXT].v.setZero(); }

    /** basic informations on the mesh */
    void infos(void) const
        {
        std::cout << "mesh:\n";
        std::cout << "  bounding box diam:  " << diam << '\n';
        std::cout << "  nodes:              " << getNbNodes() << '\n';
        std::cout << "  faces:              " << fac.size() << '\n';
        std::cout << "  tetraedrons:        " << tet.size() << '\n';
        std::cout << "  total surface:      " << surf << '\n';
        std::cout << "  total volume:       " << vol << '\n';
        }

    /** call setBasis for all nodes, and update P matrix for all elements */
    void setBasis(const double r)
        {
        std::for_each(EXEC_POL, node.begin(), node.end(),
                      [&r](Nodes::Node &nod) { nod.setBasis(r); });

        std::for_each(EXEC_POL, tet.begin(), tet.end(),
                      [](Tetra::Tet &t) { t.buildMatP();} );

        std::for_each(EXEC_POL, fac.begin(), fac.end(),
                      [](Facette::Fac &f) { f.buildMatP();} );
        }

    /** make_evol on all nodes, and returns v_max */
    double updateNodes(Eigen::Ref<Eigen::VectorXd> X, const double dt)
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

    /** build init guess for bicgstab */
    void buildInitGuess(Eigen::Ref<Eigen::VectorXd> G) const
        {
        const unsigned int NOD = node.size();

        for (unsigned int i = 0; i < NOD; i++)
            {
            Eigen::Vector3d const& n_v = getNode_v(i);//node[i].v;
            
            G(i) = n_v.dot(node[i].ep) / gamma0;
            G(NOD + i) = n_v.dot(node[i].eq) / gamma0;
            }
        }

    /** call evolution for all the nodes */
    inline void evolution(void)
        {
        std::for_each(EXEC_POL, node.begin(), node.end(),
                      [](Nodes::Node &nod) { nod.evolution(); });
        }

    /** isobarycenter */
    Eigen::Vector3d c;

    /** lengths along x,y,z axis */
    Eigen::Vector3d l;

    /** max of l coordinates, to define a bounding box */
    double diam;

    /** total surface */
    double surf;

    /** total volume of the mesh */
    double vol;

    /** face container */
    std::vector<Facette::Fac> fac;

    /** tetrahedron container */
    std::vector<Tetra::Tet> tet;

    /** surface container */
    std::vector<Mesh::Surf> s;

    /** read a solution from a file (tsv formated) and initialize fem struct to restart computation
     * from that distribution, return time
     */
    double readSol(bool VERBOSE /**< [in] */,
                   const std::string fileName /**< [in] input .sol text file */);

    /** computes an analytical initial magnetization distribution as a starting point for the
     * simulation */
    inline void init_distrib(Settings const &mySets /**< [in] */)
        {
        std::for_each(node.begin(),node.end(),[&mySets](Nodes::Node &n)
                                              {
                                              n.d[Nodes::CURRENT].u = mySets.getMagnetization(n.p);
                                              n.d[Nodes::NEXT].u = n.d[Nodes::CURRENT].u;
                                              n.d[Nodes::NEXT].phi = 0.;
                                              n.d[Nodes::NEXT].phiv = 0.;
                                              } );
        }

    /**
    average component of either u or v through getter on the whole set of tetetrahedron
    */
    double avg(std::function<double(Nodes::Node, Nodes::index)> getter /**< [in] */,
               Nodes::index d /**< [in] */) const
        {
        double sum =
                std::transform_reduce(EXEC_POL, tet.begin(), tet.end(), 0.0, std::plus{},
                                      [getter, &d](Tetra::Tet const &te)
                                      {
                                          Eigen::Matrix<double,Tetra::NPI,1> val;
                                          te.interpolation(getter, d, val);
                                          return te.weight.dot(val);
                                      });

        return sum / vol;
        }

    /** text file (tsv) writing function for a solution */
    void savesol(const int precision /**< [in] numeric precision in .sol output text file */,
                 const std::string fileName /**< [in] */,
                 std::string const &metadata /**< [in] */) const;

    /** text file (tsv) writing function for a solution of a side problem, used by electrostatSolver
     */
    bool savesol(const int precision /**< [in] */, const std::string fileName /**< [in] */,
                 std::string const &metadata /**< [in] */,
                 std::vector<double> const &val /**< [in] */) const;

    /** setter for node[i]; what_to_set will fix what is the part of the node struct to set (usefull
     * for fmm_demag.h) */
    inline void set(const int i, std::function<void(Nodes::Node &, const double)> what_to_set,
                    const double val)
        {
        what_to_set(node[i], val);
        }

private:
    /** node container: not initialized by constructor, but later while reading the mesh by member
     * function init_node */
    std::vector<Nodes::Node> node;

    /** Index of a node in the `node` vector. This vector is itself indexed by the node position in
     * the *.msh and *.sol files. In other words, the node found at `file_idx` in a file is stored
     * as `node[node_index[file_idx]]`.
     *
     * This is the inverse of the permutation we applied when sorting the nodes. */
    std::vector<int> node_index;

    /** map of the surface region physical names from mesh file */
    std::map<int, std::string> surfRegNames;

    /** map of the volume region physical names from mesh file */
    std::map<int, std::string> volRegNames;

    /** memory allocation for the nodes */
    inline void init_node(const int Nb) { node.resize(Nb); }

    /** reading mesh format 2.2 text file function */
    void readMesh(Settings const &mySets);

    /** loop on nodes to apply predicate 'whatTodo'  */
    double doOnNodes(const double init_val, const Nodes::index coord,
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

    /** return the minimum of all nodes coordinate along coord axis */
    inline double minNodes(const Nodes::index coord) const
        {
        return doOnNodes(__DBL_MAX__, coord, [](double a, double b) { return a < b; });
        }

    /** return the maximum of all nodes coordinate along coord axis */
    inline double maxNodes(const Nodes::index coord) const
        {
        return doOnNodes(__DBL_MIN__, coord, [](double a, double b) { return a > b; });
        }

    /** redefine orientation of triangular faces in accordance with the tetrahedron
    * reorientation of the tetrahedrons if needed; definition of Ms on facette elements
    Indices and orientation convention :

                        v
                      .
                    ,/
                   /
                2(ic)                                 2
              ,/|`\                                 ,/|`\
            ,/  |  `\                             ,/  |  `\
          ,/    '.   `\                         ,6    '.   `5
        ,/       |     `\                     ,/       8     `\
      ,/         |       `\                 ,/         |       `\
     0(ia)-------'.--------1(ib) --> u     0--------4--'.--------1
      `\.         |      ,/                 `\.         |      ,/
         `\.      |    ,/                      `\.      |    ,9
            `\.   '. ,/                           `7.   '. ,/
               `\. |/                                `\. |/
                  `3(id)                                `3
                     `\.
                        ` w

*/
    void indexReorder(std::vector<Tetra::prm> const &prmTetra)
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

    /** Sort the nodes along the longest axis of the sample. This should reduce the bandwidth of
     * the matrix we will have to solve for. */
    void sortNodes()
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
    };

    }  // end namespace Mesh

#endif
