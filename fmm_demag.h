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

#include "fem.h"

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
    const FPoint<FReal> centerOfBox(0., 0., 0.);/**< center of the bounding box */

/**
\class fmm 
to initialize a tree and a kernel for the computation of the demagnetizing field, and launch the computation easily with calc_demag public member
*/ 
class fmm 
    {
    public:
        /** constructor, initialize memory for tree, kernel, sources corrections, initialize all sources */
        inline fmm(Fem &fem,bool VERBOSE,const int ScalfmmNbThreads): NOD(fem.node.size()) ,FAC( fem.fac.size()) , TET( fem.tet.size()) , SRC( FAC * Facette::NPI + TET * Tetra::NPI)
            {
            omp_set_num_threads(ScalfmmNbThreads);
            
            tree = new OctreeClass(NbLevels, SizeSubLevels, boxWidth, centerOfBox);
            if (!tree) SYSTEM_ERROR;
            
            kernels=new KernelClass( NbLevels, boxWidth, centerOfBox);// kernelPreArgs... , NbLevels, boxWidth, centerOfBox);
            if (!kernels) SYSTEM_ERROR;
            
            srcDen = (FReal*) new FReal[SRC]; if (!srcDen) SYSTEM_ERROR;
            corr=(FReal*) new FReal[NOD]; if (!corr) SYSTEM_ERROR;
            
            norm = 1./(2.*fem.diam);
            
            FTic counter;

            if(VERBOSE)
                { std::cout << "Creating & Inserting particles ...\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl; }
            counter.tic();

            Pt::pt3D c = fem.c;
            FSize idxPart=0;

            std::for_each(fem.node.begin(),fem.node.end(),[this,&c,&idxPart](Nodes::Node const& n)
                {
                Pt::pt3D pTarget = norm*(n.p - c);
                const FPoint<FReal> particlePosition(pTarget.x(), pTarget.y(), pTarget.z());
                tree->insert(particlePosition, FParticleType::FParticleTypeTarget, idxPart, 0.0);
                idxPart++;
                });//end for_each sur node
            if(VERBOSE) { std::cout << "Physical nodes inserted." << std::endl; }
            
            insertCharges<Tetra::Tet,Tetra::NPI>(fem.tet,idxPart,c);
            if(VERBOSE) { std::cout << "Volume charges inserted." << std::endl; }
            
            insertCharges<Facette::Fac,Facette::NPI>(fem.fac,idxPart,c);
            if(VERBOSE) { std::cout << "Surface charges inserted." << std::endl; }
            
            counter.tac();
            if(VERBOSE) 
                { std::cout << "Done (Creating and Inserting Particles = " 
                    << counter.elapsed() << "s).\n>> ScalFMM initialized, using " 
                    << ScalfmmNbThreads << " threads.\n" << std::endl; }  
            }
        /** destructor */
        ~fmm ()
            {
            delete tree;
            delete kernels; 
            delete [] srcDen;
            delete [] corr;    
            }
        
        /**
        launch the calculation of the demag field with second order corrections
        */
        void calc_demag(Fem &fem,Settings &mySettings)
        {
        if(mySettings.verbose) { std::cout << "\t magnetostatics ..................... "; }
        FTic counter;
        counter.tic();
        
        demag<0> (fem,mySettings); // Hd(u)
        demag<1> (fem,mySettings); // Hd(v)
        counter.tac();
        if(mySettings.verbose) { std::cout << "Done in " << counter.elapsed() << " s." << std::endl; }
        }
        
    private:
        const int NOD;/**< number of nodes */
        const int FAC;/**< number of facettes */
        const int TET;/**< number of tetrahedrons */
        const int SRC;/**< number of charge sources */
        
        OctreeClass *tree    = nullptr;/**< tree initialized by constructor */
        KernelClass *kernels = nullptr;/**< kernel initialized by constructor */
        
        FReal *srcDen = nullptr;/**< source buffer */
        FReal *corr = nullptr;/**< correction coefficients buffer */
        
        double norm;/**< normalization coefficient */
        
        /**
        function template to insert volume or surface charges in tree for demag computation. class T is Tet or Fac, it must have interpolation method, second template parameter is NPI of the namespace containing class T 
        */
    template <class T,const int NPI> void insertCharges(std::vector<T> const& container,FSize &idx,Pt::pt3D const& c)
        {
        std::for_each(container.begin(),container.end(),[this,c,&idx](T const& elem)              
            {  
            double gauss[Pt::DIM][NPI];
            elem.interpolation(Nodes::get_p,gauss);
        
            for (int j=0; j<NPI; j++, idx++)
                {
                Pt::pt3D pSource = norm*(Pt::pt3D(gauss[0][j],gauss[1][j],gauss[2][j]) - c);
                const FPoint<FReal> particlePosition(pSource.x(), pSource.y(), pSource.z());
                tree->insert(particlePosition, FParticleType::FParticleTypeSource, idx, 0.0);
                }
            });//end for_each    
        }
        
    /**
    template to computes the demag field, template parameter is either 0 or 1
    */
    template <int Hv> void demag(Fem &fem,Settings &settings)
        {
        FmmClass algo(tree, kernels);
        
        std::fill_n(srcDen,SRC,0);
        std::fill_n(corr,NOD,0);
        
        int nsrc = 0;
        std::function<const Pt::pt3D (Nodes::Node)> getter;
        std::function<void (Nodes::Node &,const double)> setter;
        
        if(Hv)
            { getter = Nodes::get_v; setter = Nodes::set_phiv;}
        else
            { getter = Nodes::get_u; setter = Nodes::set_phi;}
        
        std::for_each(fem.tet.begin(),fem.tet.end(),[this,getter,&nsrc,&settings](Tetra::Tet const& tet)              
            { tet.charges(getter,srcDen,nsrc, nu0 * settings.paramTetra[tet.idxPrm].J ); });//end for_each on tet

        std::for_each(fem.fac.begin(),fem.fac.end(),[this,getter,&nsrc,&settings](Facette::Fac const& fac)
            { fac.charges(getter,srcDen,corr,nsrc); });// end for_each on fac

        fflush(NULL);

        { // reset potentials and forces - physicalValues[idxPart] = Q
        tree->forEachLeaf([this](LeafClass* leaf)
            {
            const int nbParticlesInLeaf = leaf->getSrc()->getNbParticles();
            const FVector<long long>& indexes = leaf->getSrc()->getIndexes();
            FReal* const physicalValues = leaf->getSrc()->getPhysicalValues();
            for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart)
                { physicalValues[idxPart]=srcDen[ indexes[idxPart] - NOD ]; }
            
            std::fill_n(leaf->getTargets()->getPotentials(),leaf->getTargets()->getNbParticles(),0);    
            });

        tree->forEachCell([](CellClass* cell){ cell->resetToInitialState(); });
        }// end reset

        algo.execute();

        tree->forEachLeaf([this,&fem,setter](LeafClass* leaf){
            const FReal* const potentials = leaf->getTargets()->getPotentials();
            const int nbParticlesInLeaf  = leaf->getTargets()->getNbParticles();
            const FVector<long long>& indexes  = leaf->getTargets()->getIndexes();

            for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart)
                {
                const int indexPartOrig = indexes[idxPart];
                setter(fem.node[indexPartOrig], (potentials[idxPart]*norm + corr[indexPartOrig])/(4*M_PI));
                }
            });
        }
    };//end class fmm
    
}//end namespace
#endif
