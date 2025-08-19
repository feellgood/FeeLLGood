#ifndef mesh_h
#define mesh_h

/** \file mesh.h
\brief class mesh, readMesh is expecting a mesh file in gmsh format either text or binary, from version 2.2 to the latest 4.1. The mesh has to use only first order tetraedrons and triangular facettes, mixed meshes are not allowed.
*/

#include <cmath>
#include <algorithm>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include <execution>
#pragma GCC diagnostic pop

#include "facette.h"
#include "node.h"
#include "tetra.h"
#include "feellgoodSettings.h"

namespace Mesh
    {
    /** template class to describe boundary conditions
     * T should be either a 3D vector or a double
     * indices stored in facIdx are the indices defining a surface in mesh facette list where the
     * boundary condition does apply
     * */
template <class T>
class boundaryCondition
    {
    public:
        boundaryCondition(const std::string &n, const int idx, T &v): name(n), value(v)
            { facIdx.push_back(idx); }
        std::string name;        /**< name of the corresponding region */
        std::vector<int> facIdx; /**< list of the facette indices defining a surface stored in mesh */
        T value;                 /**< value on the surface defined by the list facIdx */
    };

/** template class to regroup altogether the boundary conditions of a problem to solve
 * the vector of boundaryCondition<T> stores the various surfaces where to apply the boundary
 * condition value (type T should be double or EigenVector3d)
 * Nota Bene: values stored in each boundaryCondition<T> have all the same type, for example are
 * double, but might have different physical nature, for example a potential V and a normal current
 * density Jn. The name in each boundary condition is there to discriminate such a situation
 * */
template <class T>
class allBoundCond
    {
    public:
        std::vector<boundaryCondition<T> > BC;
        void push_back(const std::string name, const int idx, T &v)
            {
            auto it = std::find_if(BC.begin(),BC.end(),[&name](boundaryCondition<T> & bc)
                    { return (bc.name == name); });
            if (it !=BC.end())
                it->facIdx.pusk_back(idx);
            else
                { // create a new bound cond to add to BC
                boundaryCondition<T> bc(name,idx,v);
                BC.push_back(bc); // emplace_back?
                }
            }
    };

/** \class mesh
class for storing the mesh, including mesh geometry values, containers for the nodes, triangular
faces and tetrahedrons. nodes data are private. They are accessible only through getter and setter.
*/
class mesh
    {
public:
    /** constructor : read mesh file, reorder indices and computes some values related to the mesh :
     center and length along coordinates,full volume */
    mesh(Settings &mySets /**< [in] */)
        : paramTetra(mySets.paramTetra)
        {
        readMesh(mySets);
        indexReorder();

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

        // Compute the per-region volumes and the total volume.
        std::for_each(tet.begin(), tet.end(),
                [this](Tetra::Tet const &te)
                    {
                    paramTetra[te.idxPrm].volume += te.calc_vol();
                    });
        vol = std::transform_reduce(paramTetra.begin(), paramTetra.end(), 0.0,
                std::plus<>(), [](const Tetra::prm &region){ return region.volume; });

        // Build the list of all the mesh edges.
        edges.reserve(tet.size() * 6);  // 6 (non unique) edges per tetrahedron
        std::for_each(tet.begin(), tet.end(),
                [this](Tetra::Tet const &te)
                    {
                    for (int i = 0; i < 3; ++i)
                        {
                        for (int j = i + 1; j < 4; ++j)
                            { edges.push_back(std::minmax(te.ind[i], te.ind[j])); }
                        }
                    });
        std::sort(EXEC_POL, edges.begin(), edges.end());
        auto last = std::unique(EXEC_POL, edges.begin(), edges.end());
        edges.erase(last, edges.end());
        edges.shrink_to_fit();  // save memory, as this could be quite big
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

    /** make_evol on i^th node */
    inline void updateNode(int i, double vp, double vq, const double dt)
        { node[i].make_evol(vp*gamma0, vq*gamma0, dt); }

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

    /** Reference to the volume regions in Settings. */
    std::vector<Tetra::prm> &paramTetra;

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
               Nodes::index d /**< [in] */,
               int region = -1 /**< region index, or -1 for all regions */) const;

    /** Compute the maximum angle of the magnetization between two adjacent nodes. */
    double max_angle() const
        {
        double min_dot_product = std::transform_reduce(EXEC_POL, edges.begin(), edges.end(), 1.0,
                [](double a, double b){ return std::min(a, b); },
                [this](const std::pair<int, int> edge)
                    {
                    Eigen::Vector3d m1 = getNode_u(edge.first);
                    Eigen::Vector3d m2 = getNode_u(edge.second);
                    return m1.dot(m2);
                    });
        return std::acos(min_dot_product);
        }

    /** text file (tsv) writing function for a solution, node indices are zero based */
    void savesol(const int precision /**< [in] numeric precision in .sol output text file */,
                 const std::string fileName /**< [in] */,
                 std::string const &metadata /**< [in] */) const;

    /** setter for node[i]; what_to_set will fix what is the part of the node struct to set (usefull
     * for fmm_demag.h) */
    inline void set(const int i /**< [in] */,
                    std::function<void(Nodes::Node &, const double)> what_to_set /**< [in] */,
                    const double val /**< [in] */)
        { what_to_set(node[i], val); }

    /** getter for the node_index value at position i */
    inline int getNodeIndex(const int i) const { return node_index[i]; }

    /** build the list of the indices where to solve LLG from input mesh, output in lvd */
    void build_lvd(std::vector<int> &lvd);

    /** return true if this tetraedron is magnetic */
    inline bool isMagnetic(const Tetra::Tet &t) { return (paramTetra[t.idxPrm].J > 0); }

    /** Edges of all the tetrahedrons, i.e. all the unique pairs of adjacent node indices.
     * Each pair is sorted: first < second.
     * The list is sorted lexicographically, as per std::pair::operator<(). */
    std::vector<std::pair<int, int>> edges;

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
    void indexReorder();

    /** Sort the nodes along the longest axis of the sample. This should reduce the bandwidth of
     * the matrix we will have to solve for. */
    void sortNodes(Nodes::index long_axis /**< [in] */);

    /** returns the surface defined by the set of facets of indices in facIndices
     * each elementary surface triangle defined by points p0,p1,p2 is computed using norm(cross(p0p1,p0p2))/2, it is always
     * positive  */
    double surface(std::vector<int> &facIndices);

    /** prepare boundary conditions surfaces for spin accumulation problem (electrostatic + spin
     * diffusion ) */
    void buildBoundaryConditions(std::vector<Facette::prm> &paramFacette,
            Mesh::boundaryCondition<Eigen::Vector3d> &BC_spin);
    }; // end class mesh

    }  // end namespace Mesh
#endif
