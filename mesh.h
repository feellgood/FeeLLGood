#ifndef mesh_h
#define mesh_h

/** \file mesh.h
\brief class mesh, readMesh is expecting a mesh file in gmsh text format 2.2, with first order tetraedrons and triangular facettes.
*/

#include <algorithm>
#include <execution>

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
     center and length along coordinates,full volume */
    inline mesh(Settings const &mySets /**< [in] */)
        {
        readMesh(mySets);
        indexReorder(mySets.paramTetra);

        if (mySets.verbose)
            { std::cout << "  reindexed\n"; }

        double xmin = minNodes(Nodes::IDX_X);
        double xmax = maxNodes(Nodes::IDX_X);

        double ymin = minNodes(Nodes::IDX_Y);
        double ymax = maxNodes(Nodes::IDX_Y);

        double zmin = minNodes(Nodes::IDX_Z);
        double zmax = maxNodes(Nodes::IDX_Z);

        l = Eigen::Vector3d(xmax - xmin, ymax - ymin, zmax - zmin);
        c = Eigen::Vector3d(0.5 * (xmax + xmin), 0.5 * (ymax + ymin), 0.5 * (zmax + zmin));

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
        sortNodes(long_axis);

        vol = std::transform_reduce(EXEC_POL, tet.begin(), tet.end(), 0.0, std::plus{},
                                    [](Tetra::Tet const &te) { return te.calc_vol(); });
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

    /** return projection of speed at node i along ep */
    inline double getProj_ep(const int i) const {return node[i].proj_ep();}

    /** return projection of speed at node i along eq */
    inline double getProj_eq(const int i) const {return node[i].proj_eq();}

    /** setter for u0 */
    inline void set_node_u0(const int i, Eigen::Vector3d const &val)
        { node[i].d[Nodes::CURRENT].u = val; }

    /** fix to zero node[i].v */
    inline void set_node_zero_v(const int i) { node[i].d[Nodes::NEXT].v.setZero(); }

    /** basic informations on the mesh */
    void infos(void) const;

    /** call setBasis for all nodes, and update P matrix for all elements */
    inline void setBasis(const double r)
        {
        std::for_each(EXEC_POL, node.begin(), node.end(),
                      [&r](Nodes::Node &nod) { nod.setBasis(r); });
        }

    /** make_evol on all nodes, and returns v_max */
    double updateNodes(Eigen::Ref<Eigen::VectorXd> X, const double dt);

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
               Nodes::index d /**< [in] */) const;

    /** text file (tsv) writing function for a solution */
    void savesol(const int precision /**< [in] numeric precision in .sol output text file */,
                 const std::string fileName /**< [in] */,
                 std::string const &metadata /**< [in] */) const;

    /** text file (tsv) writing function for a solution of a side problem, used by electrostatSolver
     */
    bool savesol(const int precision /**< [in] */,
                 const std::string fileName /**< [in] */,
                 std::string const &metadata /**< [in] */,
                 std::vector<double> const &val /**< [in] */) const;

    /** setter for node[i]; what_to_set will fix what is the part of the node struct to set (usefull
     * for fmm_demag.h) */
    inline void set(const int i /**< [in] */,
                    std::function<void(Nodes::Node &, const double)> what_to_set /**< [in] */,
                    const double val /**< [in] */)
        { what_to_set(node[i], val); }

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

    /** test if mesh file contains surfaces and regions mentionned in yaml settings and their dimensions */
    void checkMeshFile(Settings const &mySets /**< [in] */);
    
    /** read Nodes from mesh file */
    void readNodes(Settings const &mySets /**< [in] */);

    /** read tetraedrons of the settings volume regions */
    void readTetraedrons(Settings const &mySets /**< [in] */);

    /** read facettes of the settings surface regions */
    void readTriangles(Settings const &mySets /**< [in] */);

    /** reading mesh format 2.2 text file function */
    void readMesh(Settings const &mySets /**< [in] */);

    /** loop on nodes to apply predicate 'whatTodo'  */
    double doOnNodes(const double init_val /**< [in] */,
                     const Nodes::index coord /**< [in] */,
                     std::function<bool(double, double)> whatToDo /**< [in] */) const;

    /** return the minimum of all nodes coordinate along coord axis */
    inline double minNodes(const Nodes::index coord /**< [in] */) const
        {
        return doOnNodes(__DBL_MAX__, coord, [](double a, double b) { return a < b; });
        }

    /** return the maximum of all nodes coordinate along coord axis */
    inline double maxNodes(const Nodes::index coord /**< [in] */) const
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
    void indexReorder(std::vector<Tetra::prm> const &prmTetra);

    /** Sort the nodes along the longest axis of the sample. This should reduce the bandwidth of
     * the matrix we will have to solve for. */
    void sortNodes(Nodes::index long_axis /**< [in] */);

    }; // end class mesh

    }  // end namespace Mesh
#endif
