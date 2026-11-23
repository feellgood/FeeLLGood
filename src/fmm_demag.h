#ifndef FMM_DEMAG_H
#define FMM_DEMAG_H

/** \file fmm_demag.h
\brief computation of the scalar magnetic potential by the fast multipole method. This header
prepares the particles (magnetic nodes and Gauss points) and their charges; the octree, the kernel
and the algorithm of ScalFMM are in the class FmmCore (fmm_core.h).
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <vector>

#include "fmm_core.h"
#include "mesh.h"
#include "octree_height.h"

/** \namespace scal_fmm
to grab altogether the templates and functions using scalfmm for the computation of the demag field
*/

namespace scal_fmm
    {
/** \class fmm
to initialize the FMM for the computation of the demagnetizing field, and launch the computation
easily with calc_demag public member
*/
class fmm
    {
public:
    /** constructor, initialize memory for sources and corrections, build the FMM with all the
     * particles */
    inline fmm(Mesh::mesh &msh /**< [in] */,
               std::vector<Tetra::prm> & prmTet /**< [in] */,
               std::vector<Triangle::prm> & prmTri /**< [in] */,
               const int ScalfmmNbThreads /**< [in] */)
        : prmTetra(prmTet), prmTriangle(prmTri), norm(2. / msh.l.maxCoeff())
        {
        const std::vector<Eigen::Vector3d> pos = particlePositions(msh);
        nbLevels = octreeHeight(pos, boxWidth);
        for (size_t i = 0; i < msh.magNode.size(); i++)
            {
            if (msh.magNode[i])
                { targetNodes.push_back(i); }
            }
        std::vector<double> xyz;
        xyz.reserve(3 * pos.size());
        for (const Eigen::Vector3d &p : pos)
            { xyz.insert(xyz.end(), {p.x(), p.y(), p.z()}); }
        core = std::make_unique<FmmCore>(xyz, targetNodes.size(), nbLevels, ScalfmmNbThreads);

        for (int v = 0; v < NVALS; v++)
            {
            srcDen[v].resize(msh.magTri.size() * Triangle::NPI + msh.magTet.size() * Tetra::NPI);
            corr[v].resize(msh.magNode.size());
            }
        }

    /**
    Compute the demagnetizing field. Include the second order corrections if FIRST_ORDER=OFF
    (which is the default).
    */
    void calc_demag(Mesh::mesh &msh /**< [in] */)
        {
        calc_charges(Nodes::get_u<Nodes::NEXT>, msh, srcDen[0], corr[0]);
        if constexpr (NVALS == 2)
            { calc_charges(Nodes::get_v<Nodes::NEXT>, msh, srcDen[1], corr[1]); }
        core->compute(srcDen, potentials);

        // potentials of u (phi) and, if NVALS = 2, of v (phiv) at the magnetic nodes
        const std::array<void (*)(Nodes::Node &, const double), 2> setters = {Nodes::set_phi,
                                                                             Nodes::set_phiv};
        for (int v = 0; v < NVALS; v++)
            {
            for (size_t t = 0; t < targetNodes.size(); t++)
                {
                const int node = targetNodes[t];
                msh.set(node, setters[v], (potentials[v][t] * norm + corr[v][node]) / (4 * M_PI));
                }
            }
        }

    /** height of the octree */
    int getNbLevels(void) const { return nbLevels; }

    /** corrections associated to the nodes, contributions only due to the triangles, for u and v */
    std::array<std::vector<double>, NVALS> corr;

    /** all volume region parameters for the tetraedrons */
    const std::vector<Tetra::prm> &prmTetra;

    /** all surface region parameters for the triangles */
    const std::vector<Triangle::prm> &prmTriangle;

private:
    /** sources: both surface and volume charges, for u and v */
    std::array<std::vector<double>, NVALS> srcDen;

    /** potentials at the targets, computed by the FMM, for u and v */
    std::array<std::vector<double>, NVALS> potentials;

    /** normalization coefficient */
    const double norm;

    /** height of the octree, adapted to the particles (see octree_height.h) */
    int nbLevels;

    /** indices of the magnetic nodes, which are the targets of the FMM, in the order of the
    targets */
    std::vector<int> targetNodes;

    /** octree, kernel and algorithm of the FMM, built once by the constructor */
    std::unique_ptr<FmmCore> core;

    /** normalized positions of all the particles to insert in the tree: magnetic nodes (targets),
    Gauss points of the magnetic tetrahedrons and triangles (sources) */
    std::vector<Eigen::Vector3d> particlePositions(const Mesh::mesh &msh) const
        {
        std::vector<Eigen::Vector3d> pos;
        pos.reserve(msh.magNode.size() + msh.magTet.size() * Tetra::NPI
                    + msh.magTri.size() * Triangle::NPI);
        for (size_t i = 0; i < msh.magNode.size(); i++)
            {
            if (msh.magNode[i])
                { pos.push_back(norm * (msh.getNode_p(i) - msh.c)); }
            }
        for (const int idx : msh.magTet)
            {
            Eigen::Matrix<double, Nodes::DIM, Tetra::NPI> gauss = msh.tet[idx].getPtGauss();
            for (int j = 0; j < Tetra::NPI; j++)
                { pos.push_back(norm * (gauss.col(j) - msh.c)); }
            }
        for (const int idx : msh.magTri)
            {
            Eigen::Matrix<double, Nodes::DIM, Triangle::NPI> gauss = msh.tri[idx].getPtGauss();
            for (int j = 0; j < Triangle::NPI; j++)
                { pos.push_back(norm * (gauss.col(j) - msh.c)); }
            }
        return pos;
        }

    /** computes all charges from tetraedrons and triangles for the demag field to feed a tree in
     * the fast multipole algo (scalfmm), and the corrections of the triangles at the nodes
     */
    void calc_charges(const std::function<const Eigen::Vector3d(const Nodes::Node&)>& getter,
            Mesh::mesh &msh, std::vector<double> &charges, std::vector<double> &corrections)
        {
        int nsrc(0);
        std::fill(charges.begin(),charges.end(),0);
        std::for_each(msh.magTet.begin(), msh.magTet.end(),
                [this, &msh, &getter, &nsrc, &charges](const int idx)
                {
                Tetra::Tet &t = msh.tet[idx];
                Eigen::Matrix<double,Tetra::NPI,1> result =
                        t.charges(prmTetra[t.idxPrm].Ms, getter);
                for(int i=0;i<Tetra::NPI;i++)
                    { charges[nsrc+i] = result(i); }
                nsrc += Tetra::NPI;
                });
        std::fill(corrections.begin(),corrections.end(),0);
        std::for_each(msh.magTri.begin(), msh.magTri.end(),
                [&msh, &getter, &nsrc, &charges, &corrections](const int idx)
                {
                Triangle::Tri &f = msh.tri[idx];
                Eigen::Matrix<double,Triangle::NPI,1> result = f.charges(f.dMs, getter);
                for(int i=0;i<Triangle::NPI;i++)
                    { charges[nsrc+i] = result(i); }
                nsrc += Triangle::NPI;
                f.correctionCharges(getter,result,corrections);
                });
        }

    };  // end class fmm

    }  // namespace scal_fmm
#endif
