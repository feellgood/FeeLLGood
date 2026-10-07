/** \file fmm_demag.cpp
\brief implementation of the interface to scalfmm 3. This is the only translation unit that includes
scalfmm headers.
*/

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include <omp.h>

#include "fmm_demag.h"

// THETA (config.h) is a macro, and a parameter name in the lapack interface used by scalfmm
#pragma push_macro("THETA")
#undef THETA
#include "scalfmm/algorithms/fmm.hpp"
#include "scalfmm/container/particle.hpp"
#include "scalfmm/container/point.hpp"
#include "scalfmm/interpolation/interpolation.hpp"
#include "scalfmm/matrix_kernels/laplace.hpp"
#include "scalfmm/operators/fmm_operators.hpp"
#include "scalfmm/tree/box.hpp"
#include "scalfmm/tree/cell.hpp"
#include "scalfmm/tree/for_each.hpp"
#include "scalfmm/tree/group_tree_view.hpp"
#include "scalfmm/tree/leaf_view.hpp"
#pragma pop_macro("THETA")

namespace scal_fmm
    {
const int DIM = 3;  /**< space dimension */

/** automatic tree height: minimum and maximum values (root included) */
const std::size_t minTreeHeight = 3;
const std::size_t maxTreeHeight = 12;

typedef double FReal; /**< all computations are made in double precision */

/** matrix kernel 1/r, both for near and far fields */
using MatrixKernelClass = scalfmm::matrix_kernels::laplace::one_over_r;

/** near field: direct particle to particle interactions */
using NearFieldClass = scalfmm::operators::near_field_operator<MatrixKernelClass>;

/** far field approximation: uniform interpolation, M2L accelerated by FFT */
using InterpolatorClass =
        scalfmm::interpolation::interpolator<FReal, DIM, MatrixKernelClass,
                                             scalfmm::options::uniform_<scalfmm::options::fft_>>;

/** far field operator */
using FarFieldClass = scalfmm::operators::far_field_operator<InterpolatorClass>;

/** near and far fields altogether */
using FmmOperatorClass = scalfmm::operators::fmm_operators<NearFieldClass, FarFieldClass>;

/** convenient typedef for the positions */
using PointClass = scalfmm::container::point<FReal, DIM>;

/** convenient typedef for the simulation box */
using BoxClass = scalfmm::component::box<PointClass>;

/** a source has one input (its charge) and a dummy output, its variable is the index in srcDen */
using SourceClass = scalfmm::container::particle<FReal, DIM, FReal, 1, FReal, 1, std::size_t>;

/** a target has a dummy input and one output (the potential), its variable is the node index */
using TargetClass = scalfmm::container::particle<FReal, DIM, FReal, 1, FReal, 1, std::size_t>;

/** convenient typedef for the cells */
using CellClass = scalfmm::component::cell<typename InterpolatorClass::storage_type>;

/** convenient typedef for the leaves of the source tree */
using SourceLeafClass = scalfmm::component::leaf_view<SourceClass>;

/** convenient typedef for the leaves of the target tree */
using TargetLeafClass = scalfmm::component::leaf_view<TargetClass>;

/** convenient typedef for the tree of the sources */
using SourceTreeClass = scalfmm::component::group_tree_view<CellClass, SourceLeafClass, BoxClass>;

/** convenient typedef for the tree of the targets */
using TargetTreeClass = scalfmm::component::group_tree_view<CellClass, TargetLeafClass, BoxClass>;

const double boxWidth = 2.01;            /**< bounding box max dimension */
const PointClass boxCenter{0., 0., 0.};  /**< center of the bounding box */

/** implementation of class fmm: trees, operators, charges and corrections */
struct fmm::impl
    {
    /** constructor, build the trees of sources and targets, and the fmm operators */
    impl(Mesh::mesh &msh, std::vector<Tetra::prm> &prmTet, std::vector<Triangle::prm> &prmTri,
         const int order, const int height, const int groupSize)
        : prmTetra(prmTet), prmTriangle(prmTri), norm(2. / msh.l.maxCoeff()), order(order),
          treeHeight(height > 0 ? height : autoTreeHeight(msh)), groupSize(groupSize),
          sourceTree(buildSourceTree(msh)), targetTree(buildTargetTree(msh)),
          interpolator(order, treeHeight, boxWidth), nearField(false), farField(interpolator),
          fmmOperator(nearField, farField)
        {
        srcDen.resize( msh.magTri.size()*Triangle::NPI + msh.magTet.size()*Tetra::NPI );
        corr.resize(msh.magNode.size());
        }

    /** corrections associated to the nodes, contributions only due to the triangles */
    std::vector<double> corr;

    /** all volume region parameters for the tetraedrons */
    const std::vector<Tetra::prm> &prmTetra;

    /** all surface region parameters for the triangles */
    const std::vector<Triangle::prm> &prmTriangle;

    /** sources: both surface and volume charges */
    std::vector<double> srcDen;

    /** normalization coefficient */
    double norm;

    /** order of the interpolation polynomials of the far field */
    const std::size_t order;

    /** height of the trees, root included */
    const std::size_t treeHeight;

    /** number of leaves and cells per group in the trees */
    const std::size_t groupSize;

    /** tree of the sources (Gauss points of the magnetic tetrahedrons and triangles) */
    SourceTreeClass sourceTree;

    /** tree of the targets (magnetic nodes) */
    TargetTreeClass targetTree;

    /** interpolator of the far field. The scalfmm operators below hold references: interpolator,
     * nearField and farField must outlive fmmOperator, and be declared before it. */
    InterpolatorClass interpolator;

    /** near field operator, without mutual interactions since sources and targets differ */
    NearFieldClass nearField;

    /** far field operator */
    FarFieldClass farField;

    /** near and far field operators */
    FmmOperatorClass fmmOperator;

    /** smallest tree height such that the non empty leaves hold on average at most order^3
     * particles, the number of interpolation points in a cell. This is a heuristic, to be
     * calibrated with ci-tests/benchmark-fmm.py. The number of non empty leaves is estimated assuming the
     * particles fill the bounding box of the mesh: along each direction, the number of non empty
     * leaves is max(1, l_i/max(l) * 2^(h-1)).
     */
    std::size_t autoTreeHeight(const Mesh::mesh &msh) const
        {
        const double nbParticles =
                msh.magTet.size()*Tetra::NPI + msh.magTri.size()*Triangle::NPI
                + std::count(msh.magNode.begin(), msh.magNode.end(), true);
        const double particlesPerLeaf = std::pow(order, DIM);
        const Eigen::Vector3d relativeSize = msh.l / msh.l.maxCoeff();
        std::size_t h = minTreeHeight;
        for (; h < maxTreeHeight; h++)
            {
            const double nbLeavesPerDim = std::pow(2.0, h - 1);
            double nbLeaves = 1.0;
            for (int i = 0; i < DIM; i++)
                { nbLeaves *= std::max(1.0, relativeSize(i)*nbLeavesPerDim); }
            if (nbParticles / nbLeaves <= particlesPerLeaf)
                { break; }
            }
        return h;
        }

    /** returns the normalized position of a point */
    PointClass normalized(const Eigen::Ref<const Eigen::Vector3d> p,
                          const Eigen::Ref<const Eigen::Vector3d> c) const
        {
        Eigen::Vector3d q = norm*(p - c);
        return PointClass{q.x(), q.y(), q.z()};
        }

    /** build the tree of the magnetic nodes, the particle variable is the node index */
    TargetTreeClass buildTargetTree(const Mesh::mesh &msh) const
        {
        std::vector<TargetClass> container;
        for (std::size_t idx = 0; idx < msh.magNode.size(); ++idx)
            {
            if (msh.magNode[idx])
                {
                TargetClass p;
                p.position() = normalized(msh.getNode_p(idx), msh.c);
                p.inputs(0) = 0.0;
                p.outputs(0) = 0.0;
                p.variables(idx);
                container.push_back(p);
                }
            }
        return TargetTreeClass(treeHeight, order, BoxClass(boxWidth, boxCenter), groupSize, groupSize,
                               container, false);
        }

    /** build the tree of the Gauss points of the magnetic tetrahedrons and triangles, the particle
     * variable is the index in srcDen. Tetrahedrons first, then triangles, as in calc_charges.
     */
    SourceTreeClass buildSourceTree(const Mesh::mesh &msh) const
        {
        std::vector<SourceClass> container;
        container.reserve(msh.magTet.size()*Tetra::NPI + msh.magTri.size()*Triangle::NPI);
        insertCharges<Tetra::Tet, Tetra::NPI>(msh.tet, msh.magTet, msh.c, container);
        insertCharges<Triangle::Tri, Triangle::NPI>(msh.tri, msh.magTri, msh.c, container);
        return SourceTreeClass(treeHeight, order, BoxClass(boxWidth, boxCenter), groupSize, groupSize,
                               container, false);
        }

    /**
    function template to insert volume or surface charges in a container of sources. class T is
    Tet or Tri, it must have getPtGauss() method to get the Gauss points, second template parameter
    is NPI of the namespace containing class T.
    idxContainer is the list of indices of the magnetic T elements stored in container.
    */
    template<class T, const int NPI>
    void insertCharges(const std::vector<T> &container, const std::vector<int> &idxContainer,
                       const Eigen::Ref<const Eigen::Vector3d> c,
                       std::vector<SourceClass> &sources) const
        {
        std::for_each(idxContainer.begin(), idxContainer.end(),
                      [this, &container, c, &sources](const int idxElem)
                      {
                      const T &elem = container[idxElem];
                      Eigen::Matrix<double,Nodes::DIM,NPI> gauss = elem.getPtGauss();

                      for (int j = 0; j < NPI; j++)
                          {
                          SourceClass p;
                          p.position() = normalized(gauss.col(j), c);
                          p.inputs(0) = 0.0;
                          p.outputs(0) = 0.0;
                          p.variables(sources.size());
                          sources.push_back(p);
                          }
                      });
        }

    /** computes all charges from tetraedrons and triangles for the demag field to feed a tree in the
     * fast multipole algo (scalfmm)
     */
    void calc_charges(const std::function<const Eigen::Vector3d(const Nodes::Node&)>& getter,
            Mesh::mesh &msh)
        {
        int nsrc(0);
        std::fill(srcDen.begin(),srcDen.end(),0);
        std::for_each(msh.magTet.begin(),msh.magTet.end(),[this, &msh, &getter, &nsrc](const int idx)
                {
                Tetra::Tet &t = msh.tet[idx];
                Eigen::Matrix<double,Tetra::NPI,1> result =
                        t.charges(prmTetra[t.idxPrm].Ms, getter);
                for(int i=0;i<Tetra::NPI;i++)
                    { srcDen[nsrc+i] = result(i); }
                nsrc += Tetra::NPI;
                });
        std::fill(corr.begin(),corr.end(),0);
        std::for_each(msh.magTri.begin(),msh.magTri.end(),[this, &msh, &getter, &nsrc](const int idx)
                {
                Triangle::Tri &f = msh.tri[idx];
                Eigen::Matrix<double,Triangle::NPI,1> result = f.charges(f.dMs, getter);
                for(int i=0;i<Triangle::NPI;i++)
                    { srcDen[nsrc+i] = result(i); }
                nsrc += Triangle::NPI;
                f.correctionCharges(getter,result,corr);
                });
        }

    /**
    computes the demag field, with (getter  = u,setter = phi) or (getter = v,setter = phi_v)
    */
    void demag(const std::function<const Eigen::Vector3d(const Nodes::Node&)>& getter,
               const std::function<void(Nodes::Node &, const double)>& setter, Mesh::mesh &msh)
        {
        calc_charges(getter, msh);

        // physical values of the sources are the charges
        scalfmm::component::for_each_leaf(sourceTree.begin(), sourceTree.end(),
                [this](SourceLeafClass &leaf)
                {
                for (auto p_ref : leaf)
                    {
                    auto p = typename SourceLeafClass::proxy_type(p_ref);
                    p.inputs(0) = srcDen[std::get<0>(p.variables())];
                    }
                });

        // reset potentials, multipoles and local expansions
        targetTree.reset_outputs();
        sourceTree.reset_multipoles();
        targetTree.reset_locals();

        scalfmm::algorithms::omp::task_dep(sourceTree, targetTree, fmmOperator);

        scalfmm::component::for_each_leaf(targetTree.begin(), targetTree.end(),
                [this, &msh, &setter](TargetLeafClass &leaf)
                {
                for (auto const p_ref : leaf)
                    {
                    const auto p = typename TargetLeafClass::const_proxy_type(p_ref);
                    const std::size_t idx = std::get<0>(p.variables());
                    msh.set(idx, setter, (p.outputs(0) * norm + corr[idx]) / (4 * M_PI));
                    }
                });
        }
    };  // end struct fmm::impl

fmm::fmm(Mesh::mesh &msh, std::vector<Tetra::prm> &prmTet, std::vector<Triangle::prm> &prmTri,
         const int ScalfmmNbThreads, const int order, const int treeHeight, const int groupSize)
    {
    omp_set_num_threads(ScalfmmNbThreads);
    pImpl = std::make_unique<impl>(msh, prmTet, prmTri, order, treeHeight, groupSize);
    }

fmm::~fmm() = default;

int fmm::treeHeight() const { return pImpl->treeHeight; }

void fmm::calc_demag(Mesh::mesh &msh)
    {
    pImpl->demag(Nodes::get_u<Nodes::NEXT>, Nodes::set_phi, msh);
    if (!FIRST_ORDER)
        { pImpl->demag(Nodes::get_v<Nodes::NEXT>, Nodes::set_phiv, msh); }
    }

    }  // namespace scal_fmm
