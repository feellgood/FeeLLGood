#ifndef mesh_h
#define mesh_h

/** \file mesh.h
\brief class mesh, readMesh is expecting a mesh file in gmsh format either text or binary, from
version 2.2 to the latest 4.1. The mesh has to use only first order tetraedrons and triangles,
mixed meshes are not allowed.
*/

#include <cmath>
#include <algorithm>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include <execution>
#pragma GCC diagnostic pop

#include "triangle.h"
#include "node.h"
#include "tetra.h"
#include "settings.h"

class BasicTri;

namespace Mesh
    {
    /** A mesh edge is a sorted pair of adjacent node indices: first < second. */
    using Edge = std::pair<int, int>;

/** \class mesh
class for storing the mesh, including mesh geometry values, containers for the nodes, triangles
and tetrahedrons. nodes data are private. They are accessible only through getter and setter.
*/
class mesh
    {
public:
    /** constructor : read mesh file, reorder indices and computes some values related to the mesh :
     center and length along coordinates,full volume */
    mesh(Settings &mySets /**< [in] */,
         bool simplified = false /** pass true in the unit tests */)
        : paramTriangle(mySets.paramTriangle), paramTetra(mySets.paramTetra), volumeRegions(mySets.paramTetra.size())
        {
        readMesh(mySets);
        if (simplified)
            { return; }
        if (!controlTriangles())
            { exit(1); }

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

        totalMagVol = 0;
        // Compute the per-region volumes and the total volume.
        std::for_each(tet.begin(), tet.end(),
                [this](const Tetra::Tet &te)
                    {
                    double vol_tet = te.calc_vol();
                    paramTetra[te.idxPrm].volume += vol_tet;
                    if(isMagnetic(te))
                        totalMagVol += vol_tet;
                    });
        vol = std::transform_reduce(paramTetra.begin(), paramTetra.end(), 0.0,
                std::plus<>(), [](const Tetra::prm &region){ return region.volume; });

        // Build the list of tetrahedrons for each region.
        std::for_each(tet.begin(), tet.end(), [this](const Tetra::Tet &te)
            {
            volumeRegions[te.idxPrm].push_back(te.idx);
            });

        // Build the list of all the mesh edges.
        edges.reserve(tet.size() * 6);  // 6 (non-unique) edges per tetrahedron
        std::for_each(tet.begin(), tet.end(),
                [this](const Tetra::Tet &te)
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
        magNode.resize(node.size());
        std::fill(magNode.begin(),magNode.end(),false);
        std::for_each(tet.begin(),tet.end(),[this](Tetra::Tet &te)
                {
                if(isMagnetic(te))
                    {
                    magTet.push_back(te.idx);
                    for (int i=0;i<4;i++)
                        { magNode[te.ind[i]] = true; }
                    }
                });

        for(unsigned int i=0;i<tri.size();i++)
            {
            if (isMagnetic(tri[i]) && !mySets.paramTriangle[tri[i].idxPrm].suppress_charges)
                { magTri.push_back(i); }
            }

        if(mySets.getFieldType() == R4toR3)
            { setExtSpaceField(mySets); }
        }

    /** operator copy deleted */
    mesh(const mesh &) = delete;

    /** operator= deleted */
    mesh &operator=(const mesh &) = delete;

    /** return number of nodes  */
    inline int getNbNodes(void) const { return node.size(); }

    /** return number of triangles */
    inline int getNbTris(void) const { return tri.size(); }

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
    inline void set_node_u0(const int i, const Eigen::Vector3d &val)
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
    inline void updateNode(const int i, const double vp, const double vq, const double dt)
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

    /** total magnetic volume of the mesh */
    double totalMagVol;

    /** triangle container. Contains only those on a surface region. */
    std::vector<Triangle::Tri> tri;

    /** Reference to the surface regions in Settings. */
    std::vector<Triangle::prm> &paramTriangle;

    /** tetrahedron container */
    std::vector<Tetra::Tet> tet;

    /** Reference to the volume regions in Settings. */
    std::vector<Tetra::prm> &paramTetra;

    /** read a solution from a file (tsv formated) and initialize fem struct to restart computation
     * from that distribution, return time
     */
    double readSol(bool VERBOSE /**< [in] */,
                   const std::string& fileName /**< [in] input .sol text file */);

    /** computes an analytical initial magnetization distribution as a starting point for the
     * simulation. If the node is not magnetic then it is set to NAN. */
    inline void init_distrib(const Settings &mySets /**< [in] */)
        {
        for (int nodeIdx = 0; nodeIdx < int(node.size()); ++nodeIdx)
            {
            Nodes::Node &n = node[nodeIdx];
            n.d[Nodes::NEXT].phi = 0.;
            n.d[Nodes::NEXT].phiv = 0.;

            // A non-magnetic node's magnetization is a vector of NAN.
            if (!magNode[nodeIdx])
                {
                n.d[Nodes::CURRENT].u = Eigen::Vector3d(NAN, NAN, NAN);
                n.d[Nodes::NEXT].u = n.d[Nodes::CURRENT].u;
                continue;
                }

            // If the initial magnetization depends only on the node position, we do not need to
            // build the list of region names.
            if (mySets.getMagType() == POSITION_ONLY)
                {
                n.d[Nodes::CURRENT].u = mySets.getMagnetization(n.p);
                n.d[Nodes::NEXT].u = n.d[Nodes::CURRENT].u;
                continue;
                }

            // Get the list of region indices this node belongs to.
            std::set<int> nodeRegions;
            // skip region 0 (__default__) which should be empty
            for (size_t regIdx = 1; regIdx < volumeRegions.size(); ++regIdx)
                {
                const std::vector<int> &regionTetras = volumeRegions[regIdx];
                bool node_in_region = std::any_of(EXEC_POL,
                    regionTetras.begin(), regionTetras.end(),
                    [this, nodeIdx](const int tetIdx)
                    {
                    const Tetra::Tet &tetrahedron = tet[tetIdx];
                    for (int i = 0; i < Tetra::N; ++i) // node of tetrahedron tetIdx
                        {
                        if (tetrahedron.ind[i] == nodeIdx)
                            { return true; }
                        }
                    return false;
                    });
                if (node_in_region)
                    { nodeRegions.insert(regIdx); }
                }

            // Get the list of region names this node belongs to.
            std::vector<std::string> region_names;
            region_names.resize(nodeRegions.size());
            std::transform(nodeRegions.begin(), nodeRegions.end(), region_names.begin(),
                [this](const int regIdx){ return paramTetra[regIdx].regName; });

            n.d[Nodes::CURRENT].u = mySets.getMagnetization(n.p, region_names);
            n.d[Nodes::NEXT].u = n.d[Nodes::CURRENT].u;
            }
        }

    /**
    returns Thiele length in the magnetic volume region of index region.
    this is an integral on a volume region of sq(grad(u).grad(u)) with usual weighted dot on each
    tetra. It is normalized by a surface: the average cross-section
    
    */
    double thiele(const int region /** region index, or -1 for all magnetic regions */) const;

    /**
    average component of either u or v through getter on the whole set of tetetrahedron
    */
    double avg(const std::function<double(Nodes::Node, Nodes::index)>& getter /**< [in] */,
               Nodes::index d /**< [in] */,
               int region = -1 /**< region index, or -1 for all magnetic regions */) const;

    /** Compute the maximum angle of the magnetization between two adjacent nodes. */
    double max_angle() const
        {
        double min_dot_product = std::transform_reduce(EXEC_POL, edges.begin(), edges.end(), 1.0,
                [](const double a, const double b){ return std::min(a, b); },
                [this](const Edge &edge)
                    {
                    Eigen::Vector3d m1 = getNode_u(edge.first);
                    Eigen::Vector3d m2 = getNode_u(edge.second);
                    return m1.dot(m2);
                    });
        return std::acos(min_dot_product);
        }

    /** text file (tsv) writing function for a solution, node indices are zero based
     * first column is index, then magnetization components, then scalar magnetic potential.
     * If spin accumulation was involved, there are three extra columns for the componants of the
     * spin diffusion vector.
     * Since some volume region might be non magnetic, magnetization is undefined on those nodes and
     * nan is used */
    void savesol(const int precision /**< [in] numeric precision in .sol output text file */,
            const std::string& fileName /**< [in] */,
             const std::string &metadata /**< [in] */,
             bool withSpinAcc /**< [in] */,
             std::vector<Eigen::Vector3d> &s /**< [in] spin accumulation (might be empty) */) const;

    /** setter for node[i]; what_to_set will fix what is the part of the node struct to set (usefull
     * for fmm_demag.h) */
    inline void set(const int i /**< [in] */,
                    const std::function<void(Nodes::Node &, const double)>& what_to_set /**< [in] */,
                    const double val /**< [in] */)
        { what_to_set(node[i], val); }

    /** getter for the node_index value at position i */
    inline int getNodeIndex(const int i) const { return node_index[i]; }

    /** return true if this tetraedron is magnetic */
    inline bool isMagnetic(const Tetra::Tet &t) const { return (paramTetra[t.idxPrm].Ms > 0); }

    /** return true if this triangle is magnetic */
    inline bool isMagnetic(const Triangle::Tri &f)
        { return (magNode[f.ind[0]] && magNode[f.ind[1]] && magNode[f.ind[2]]); }

    /** Update dMs on each triangular element, with a sign that is consistent with the triangle's
     * orientation.
     * Check whether each triangle of tet and tri satisfies one of the three following conditions :
     *      - Is part of no surface region, but part of 2 tetraedrons of the same volume region
     *      - Is part of at most one surface region, and 2 tetraedrons of differents volume region
     *      - Is part of 1 tetraedron and at most one surface region
     * Returns true if all triangles in one of those three cases, false otherwise.
     */
    bool controlTriangles();

    /** Edges of all the tetrahedrons. Each edge is a sorted pair of indices.
     * The list is sorted lexicographically, as per std::pair::operator<(). */
    std::vector<Edge> edges;

    /** list of the magnetic nodes, using inner indices (non gmsh indices). If true it is magnetic.
     * a node is magnetic if it belongs to a magnetic tetrahedron. Consequently any node on a
     * magnetic/non magnetic interface is set to true in magNode. */
    std::vector<bool> magNode;

     /** list of the indices of all magnetic tetrahedrons from all volume regions */
    std::vector<int> magTet;

    /** List of the indices of all surface triangles which can hold magnetic charges,
     * i.e. the triangles where:
     *   * all the vertices are magnetic nodes
     *   * the `suppress_charges` flag is unset. */
    std::vector<int> magTri;

    /** external applied space field, values on gauss points, size is number of tetraedrons */
    std::vector< Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> > extSpaceField;

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

    /** List of tetrahedrons making each volume region. */
    std::vector<std::vector<int>> volumeRegions;

    /** test if mesh file contains surfaces and regions mentionned in yaml settings and their
     * dimensions */
    void checkMeshFile(const Settings &mySets /**< [in] */);
    
    /** read Nodes from mesh file */
    void readNodes(const Settings &mySets /**< [in] */);

    /** read tetraedrons of the settings volume regions */
    void readTetraedrons(const Settings &mySets /**< [in] */);

    /** read triangles of the settings surface regions */
    void readTriangles(const Settings &mySets /**< [in] */);

    /** reading mesh format 2.2 text file function */
    void readMesh(const Settings &mySets /**< [in] */);

    /** loop on nodes to apply predicate 'whatTodo'  */
    double doOnNodes(const double init_val /**< [in] */,
                     const Nodes::index coord /**< [in] */,
                     const std::function<bool(double, double)>& whatToDo /**< [in] */) const;

    /** return the minimum of all nodes coordinate along coord axis */
    inline double minNodes(const Nodes::index coord /**< [in] */) const
        {
        return doOnNodes(__DBL_MAX__, coord, [](const double a, const double b)
                        { return a < b; });
        }

    /** return the maximum of all nodes coordinate along coord axis */
    inline double maxNodes(const Nodes::index coord /**< [in] */) const
        {
        return doOnNodes(__DBL_MIN__, coord, [](const double a, const double b)
                        { return a > b; });
        }

    /** Sort the nodes along the longest axis of the sample. This should reduce the bandwidth of
     * the matrix we will have to solve for. */
    void sortNodes(Nodes::index long_axis /**< [in] */);

    /** Called in controlTriangles. Returns true if a generation error is detected, false otherwise. */
    bool findErrInTriangle(const int nbSurfTri, const int nbTetraFaces,
                         const std::pair<int,int> pairIdVolRegs, const int idSurfReg);

    /** Called in controlTriangles. If a minor generation error is detected, updates indNodTriToAdd and
     * regsTriToAdd. Returns true if no major generation error is detected, false otherwise. */
    bool diffTriHandler(const size_t inf, const size_t sup, const std::vector<BasicTri> &allTriCtnr,
                        std::vector<size_t> &indNodTriToAdd, std::vector<std::pair<int,int>> &regsTriToAdd);

    /** returns the surface defined by the set of triangle of indices in triIndices
     * each elementary surface triangle defined by points p0,p1,p2 is computed using
     * norm(cross(p0p1,p0p2))/2, it is always positive  */
    double surface(std::vector<int> &triIndices) const;

    /** when external applied field is of field_type R4toR3 values of field_space are stored in
     * spaceField */
    void setExtSpaceField(Settings &s /**< [in] */);
    }; // end class mesh
    }  // end namespace Mesh
#endif
