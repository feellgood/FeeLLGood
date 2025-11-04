#ifndef FMM_DEMAG_H
#define FMM_DEMAG_H

/** \file fmm_demag.h
\brief this header is the interface to scalfmm. Its purpose is to prepare an octree for the application of
the fast multipole algorithm, and to compute the scalar magnetic potential and the demagnetizing field.
*/

#include "Components/FParticleType.hpp"
#include "Components/FTypedLeaf.hpp"
#include "Containers/FOctree.hpp"
#include "Core/FFmmAlgorithmThreadTsm.hpp"
#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"
#include "Kernels/Rotation/FRotationCell.hpp"
#include "Kernels/Rotation/FRotationKernel.hpp"

#include "mesh.h"

/** \namespace scal_fmm
to grab altogether the templates and functions using scalfmm for the computation of the demag field
*/

namespace scal_fmm
    {
const int P = 9;              /**< truncation of the spherical harmonics series */
const int NbLevels = 6;       /**< number of levels in the tree */
const int SizeSubLevels = 3;  /**< size of the sub levels  */

typedef double FReal; /**< parameter of scalfmm templates, all computations are made in double precision */

typedef FTypedRotationCell<FReal, P>
        CellClass; /**< convenient typedef for the definition of cell type in scalfmm  */

typedef FP2PParticleContainerIndexed<FReal>
        ContainerClass; /**< convenient typedef for the definition of container for scalfmm */

typedef FTypedLeaf<FReal, ContainerClass>
        LeafClass; /**< convenient typedef for the definition of leaf for scalfmm  */

typedef FOctree<FReal, CellClass, ContainerClass, LeafClass>
        OctreeClass; /**< convenient typedef for the definition of the octree for scalfmm */

typedef FRotationKernel<FReal, CellClass, ContainerClass, P>
        KernelClass; /**< convenient typedef for the kernel for scalfmm */

typedef FFmmAlgorithmThreadTsm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass>
        FmmClass; /**< convenient typedef for handling altogether the differents scalfmm object
                     templates used in feellgood  */

const double boxWidth = 2.01;              /**< bounding box max dimension */
const FPoint<FReal> boxCenter(0., 0., 0.); /**< center of the bounding box */

/** \class fmm
to initialize a tree and a kernel for the computation of the demagnetizing field, and launch the
computation easily with calc_demag public member
*/
class fmm
    {
public:
    /** constructor, initialize memory for tree, kernel, sources corrections, initialize all sources
     */
    inline fmm(Mesh::mesh &msh /**< [in] */,
               std::vector<Tetra::prm> & prmTet /**< [in] */,
               std::vector<Facette::prm> & prmFac /**< [in] */,
               const int ScalfmmNbThreads /**< [in] */)
        : prmTetra(prmTet), prmFacette(prmFac), NOD(msh.getNbNodes()),
          tree(NbLevels, SizeSubLevels, boxWidth, boxCenter), kernels(NbLevels, boxWidth, boxCenter)
        {
        omp_set_num_threads(ScalfmmNbThreads);
        norm = 2. / msh.l.maxCoeff();

        FSize idxPart = 0;
        for (idxPart = 0; idxPart < NOD; ++idxPart)
            {
            Eigen::Vector3d pTarget = norm*(msh.getNode_p(idxPart) - msh.c);
            tree.insert(FPoint<FReal>(pTarget.x(), pTarget.y(), pTarget.z()),
                        FParticleType::FParticleTypeTarget, idxPart);
            }

        insertCharges<Tetra::Tet, Tetra::NPI>(msh.tet, idxPart, msh.c);
        insertCharges<Facette::Fac, Facette::NPI>(msh.fac, idxPart, msh.c);

        srcDen.resize( msh.getNbFacs()*Facette::NPI + msh.getNbTets()*Tetra::NPI );
        corr.resize(NOD);
        }

    /**
    launch the calculation of the demag field with second order corrections
    */
    void calc_demag(Mesh::mesh &msh /**< [in] */)
        {
        demag(Nodes::get_u<Nodes::NEXT>, Nodes::set_phi, msh);
        demag(Nodes::get_v<Nodes::NEXT>, Nodes::set_phiv, msh);
        }

    /** sources */
    std::vector<double> srcDen;
    
    /** corrections associated to the nodes, contributions only due to the facettes */
    std::vector<double> corr;

    /** all volume region parameters for the tetraedrons */
    std::vector<Tetra::prm> prmTetra;

    /** all surface region parameters for the facettes */
    std::vector<Facette::prm> prmFacette;

private:
    const int NOD; /**< number of nodes */

    OctreeClass tree;    /**< tree initialized by constructor */

    KernelClass kernels; /**< kernel initialized by constructor */

    double norm; /**< normalization coefficient */

    /**
    function template to insert volume or surface charges in tree for demag computation. class T is
    Tet or Fac, it must have getPtGauss() method, second template parameter is NPI of the namespace
    containing class T
    */
    template<class T, const int NPI>
    void insertCharges(std::vector<T> const &container, FSize &idx, Eigen::Ref<Eigen::Vector3d> const c)
        {
        std::for_each(container.begin(), container.end(),
                      [this, c, &idx](T const &elem)
                      {
                          Eigen::Matrix<double,Nodes::DIM,NPI> gauss;
                          elem.getPtGauss(gauss);

                          for (int j = 0; j < NPI; j++, idx++)
                              {
                              double x = norm*(gauss(0,j) - c.x());
                              double y = norm*(gauss(1,j) - c.y());
                              double z = norm*(gauss(2,j) - c.z());
                              tree.insert(FPoint<FReal>(x,y,z), FParticleType::FParticleTypeSource, idx, 0.0);
                              }
                      });  // end for_each
        }

    /** computes all charges from tetraedrons and facettes for the demag field to feed a tree in the fast multipole algo (scalfmm)
     */
    void calc_charges(std::function<const Eigen::Vector3d(Nodes::Node)> getter, Mesh::mesh &msh)
        {
        int nsrc(0);
        std::fill(srcDen.begin(),srcDen.end(),0);

        std::for_each(msh.tet.begin(), msh.tet.end(),
                      [this, getter, &nsrc](Tetra::Tet const &tet)
                          { tet.charges(prmTetra[tet.idxPrm], getter, srcDen, nsrc); });
        std::fill(corr.begin(),corr.end(),0);
        std::for_each(msh.fac.begin(), msh.fac.end(),
                      [this, getter, &nsrc](Facette::Fac const &fac)
                          { fac.charges(prmFacette[fac.idxPrm], getter, srcDen, nsrc, corr); });
        }

    /**
    computes the demag field, with (getter  = u,setter = phi) or (getter = v,setter = phi_v)
    */
    void demag(std::function<const Eigen::Vector3d(Nodes::Node)> getter,
               std::function<void(Nodes::Node &, const double)> setter, Mesh::mesh &msh)
        {
        FmmClass algo(&tree, &kernels);
        calc_charges(getter, msh);

        // reset potentials and forces - physicalValues[idxPart] = Q
        tree.forEachLeaf(
                [this](LeafClass *leaf)
                {
                    const int nbParticlesInLeaf = leaf->getSrc()->getNbParticles();
                    const auto &indexes = leaf->getSrc()->getIndexes();
                    FReal *const physicalValues = leaf->getSrc()->getPhysicalValues();
                    for (int idxPart = 0; idxPart < nbParticlesInLeaf; ++idxPart)
                        {
                        physicalValues[idxPart] = srcDen[indexes[idxPart] - NOD];
                        }

                    std::fill_n(leaf->getTargets()->getPotentials(),
                                leaf->getTargets()->getNbParticles(), 0);
                });

        tree.forEachCell([](CellClass *cell) { cell->resetToInitialState(); });

        algo.execute();

        tree.forEachLeaf(
                [this, &msh, setter](LeafClass *leaf)
                {
                    const FReal *const potentials = leaf->getTargets()->getPotentials();
                    const int nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
                    const auto &indexes = leaf->getTargets()->getIndexes();
                    for (int idxPart = 0; idxPart < nbParticlesInLeaf; ++idxPart)
                        {
                        const int indexPartOrig = indexes[idxPart];
                        msh.set(indexPartOrig, setter,
                                (potentials[idxPart] * norm + corr[indexPartOrig]) / (4 * M_PI));
                        }
                });
        }
    };  // end class fmm

    }  // namespace scal_fmm
#endif
