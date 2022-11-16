#ifndef FMM_DEMAG_H
#define FMM_DEMAG_H

/** \file fmm_demag.h
this header is the interface to scalfmm. Its purpose is to prepare an octree for the application of the fast multipole algorithm, and to compute the demag field.
*/

// scalFMM includes

#include "Utils/FTic.hpp"
#include "Utils/FParameters.hpp"

#include "Components/FTypedLeaf.hpp"
#include "Components/FParticleType.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"

#include "Core/FFmmAlgorithmThreadTsm.hpp"

#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"
#include "Kernels/Rotation/FRotationKernel.hpp"
#include "Kernels/Rotation/FRotationCell.hpp"

#include "mesh.h"

/**
\namespace scal_fmm to grab altogether the templates and functions using scalfmm for the computation of the demag field 
*/

namespace scal_fmm{
    
    const int P = 9;/**< constant parameter for some scalfmm templates */
    
    typedef double FReal; /**< parameter of scalfmm templates */
    
    typedef FTypedRotationCell<FReal, P>            CellClass; /**< convenient typedef for the definition of cell type in scalfmm  */

    typedef FP2PParticleContainerIndexed<FReal>         ContainerClass; /**< convenient typedef for the definition of container for scalfmm */

    typedef FTypedLeaf<FReal, ContainerClass >                      LeafClass;/**< convenient typedef for the definition of leaf for scalfmm  */    

    typedef FOctree< FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;/**< convenient typedef for the definition of the octree for scalfmm */

    typedef FRotationKernel< FReal, CellClass, ContainerClass, P >          KernelClass;/**< convenient typedef for the kernel for scalfmm */

    typedef FFmmAlgorithmThreadTsm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;/**< convenient typedef for handling altogether the differents scalfmm object templates used in feellgood  */
    
    const int NbLevels = 8; /**< number of levels in the tree */
    const int SizeSubLevels = 6; /**< size of the sub levels  */
    const double boxWidth=2.01;/**< bounding box max dimension */
    const FPoint<FReal> boxCenter(0., 0., 0.);/**< center of the bounding box */

/**
\class fmm 
to initialize a tree and a kernel for the computation of the demagnetizing field, and launch the computation easily with calc_demag public member
*/ 
class fmm 
    {
    public:
        /** constructor, initialize memory for tree, kernel, sources corrections, initialize all sources */
        inline fmm(mesh &msh,bool VERBOSE,const int ScalfmmNbThreads): NOD(msh.getNbNodes()) ,FAC( msh.getNbFacs()) , TET( msh.getNbTets()) , SRC( FAC * Facette::NPI + TET * Tetra::NPI) , tree(NbLevels, SizeSubLevels, boxWidth, boxCenter), kernels( NbLevels, boxWidth, boxCenter)
            {
            omp_set_num_threads(ScalfmmNbThreads);
            norm = 1./(2.*msh.diam);
            
            FTic counter;
            counter.tic();

            FSize idxPart=0;
            for(idxPart=0; idxPart< NOD; ++idxPart)
                {
                Pt::pt3D pTarget = norm*(msh.getNode(idxPart).p - msh.c);
                tree.insert( FPoint<FReal>(pTarget.x(), pTarget.y(), pTarget.z()) , FParticleType::FParticleTypeTarget, idxPart);//, 0.0);    
		//tree.insert( FPoint<FReal>(pTarget.x(), pTarget.y(), pTarget.z()) , FParticleType::target, idxPart);                
		}
            
            insertCharges<Tetra::Tet,Tetra::NPI>(msh.tet,idxPart,msh.c);
            insertCharges<Facette::Fac,Facette::NPI>(msh.fac,idxPart,msh.c);
            
            counter.tac();
            if(VERBOSE) 
                { std::cout << "Nodes, volume & surface charges inserted.\nDone (Creating and Inserting Particles = " 
                    << counter.elapsed() << "s), using " << ScalfmmNbThreads << " threads.\n" << std::endl; }  
            }
        
        /**
        launch the calculation of the demag field with second order corrections
        */
        void calc_demag(mesh &msh,Settings &mySettings)
        {
        FTic counter;
        counter.tic();
        
        demag(Nodes::get_u, Nodes::set_phi, msh, mySettings);
        demag(Nodes::get_v, Nodes::set_phiv, msh, mySettings);

        counter.tac();
        if(mySettings.verbose) { std::cout << "Magnetostatics done in " << counter.elapsed() << " s." << std::endl; }
        }
        
    private:
        const int NOD;/**< number of nodes */
        const int FAC;/**< number of facettes */
        const int TET;/**< number of tetrahedrons */
        const int SRC;/**< number of charge sources */
        
        OctreeClass tree;/**< tree initialized by constructor */
        KernelClass kernels;/**< kernel initialized by constructor */
        
        double norm;/**< normalization coefficient */
        
        /**
        function template to insert volume or surface charges in tree for demag computation. class T is Tet or Fac, it must have interpolation method, second template parameter is NPI of the namespace containing class T 
        */
    template <class T,const int NPI> void insertCharges(std::vector<T> const& container,FSize &idx,Pt::pt3D const& c)
        {
        std::for_each(container.begin(),container.end(),[this,c,&idx](T const& elem)              
            {  
            Pt::pt3D gauss[NPI];
            elem.interpolation(Nodes::get_p,gauss);// for facette it is interpolation<Pt::pt3D> called here
        
            for (int j=0; j<NPI; j++, idx++)
                {
                Pt::pt3D pSource = norm*(gauss[j] - c);
                tree.insert( FPoint<FReal>(pSource.x(), pSource.y(), pSource.z()) , FParticleType::FParticleTypeSource, idx, 0.0);
                //tree.insert( FPoint<FReal>(pSource.x(), pSource.y(), pSource.z()) , FParticleType::source, idx, 0.0);
		}
            });//end for_each    
        }
        
    /**
    computes the demag field, with (getter  = u,setter = phi) or (getter = v,setter = phi_v)
    */

    void demag(std::function<const Pt::pt3D (Nodes::Node)> getter,std::function<void (Nodes::Node &,const double)> setter,mesh &msh,Settings &settings)
        {
        FmmClass algo(&tree, &kernels);
        
        std::vector<double> srcDen(SRC,0);
        std::vector<double> corr(NOD,0);
        
        msh.calc_charges(getter,srcDen,corr,settings);
        
        // reset potentials and forces - physicalValues[idxPart] = Q
        tree.forEachLeaf([this,&srcDen](LeafClass* leaf)
            {
            const int nbParticlesInLeaf = leaf->getSrc()->getNbParticles();
            const FVector<long long>& indexes = leaf->getSrc()->getIndexes();
            FReal* const physicalValues = leaf->getSrc()->getPhysicalValues();
            for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart)
                { physicalValues[idxPart]=srcDen[ indexes[idxPart] - NOD ]; }
            
            std::fill_n(leaf->getTargets()->getPotentials(),leaf->getTargets()->getNbParticles(),0);    
            });

        tree.forEachCell([](CellClass* cell){ cell->resetToInitialState(); });
        
        algo.execute();

        tree.forEachLeaf([this,&corr,&msh,setter](LeafClass* leaf){
            const FReal* const potentials = leaf->getTargets()->getPotentials();
            const int nbParticlesInLeaf  = leaf->getTargets()->getNbParticles();
            const FVector<long long>& indexes  = leaf->getTargets()->getIndexes();

            for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart)
                {
                const int indexPartOrig = indexes[idxPart];
                //setter(msh.setNode(indexPartOrig), (potentials[idxPart]*norm + corr[indexPartOrig])/(4*M_PI));
                msh.set(indexPartOrig,setter,(potentials[idxPart]*norm + corr[indexPartOrig])/(4*M_PI) );
                //setter(msh.node[indexPartOrig], (potentials[idxPart]*norm + corr[indexPartOrig])/(4*M_PI));
                }
            });
        }
    };//end class fmm
    
}//end namespace
#endif
