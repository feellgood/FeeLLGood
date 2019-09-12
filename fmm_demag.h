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



/** double redefinition for the parametrization of some scalfmm templates */
#define FReal double



/**
\namespace scal_fmm to grab altogether the templates and functions using scalfmm for the computation of the demag field 
*/

namespace scal_fmm{
    
    /** constant parameter for some scalfmm templates */
    const int P = 9;
    
    // pb avec ces templates : ils prennent un argument de plus en 1.5.0 et 1.4.148 //ct
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
        /** constructor */
        inline fmm(Fem &fem,bool VERBOSE,const int ScalfmmNbThreads): NOD(fem.node.size()) ,FAC( fem.fac.size()) , TET( fem.tet.size()) , SRC( FAC * Facette::NPI + TET * Tetra::NPI)
            {
            omp_set_num_threads(ScalfmmNbThreads);
            
            tree = new OctreeClass(NbLevels, SizeSubLevels, boxWidth, centerOfBox);
            if (!tree) SYSTEM_ERROR;
            
            kernels=new KernelClass( NbLevels, boxWidth, centerOfBox);// kernelPreArgs... , NbLevels, boxWidth, centerOfBox);
            if (!kernels) SYSTEM_ERROR;
            
            srcDen = (FReal*) new FReal[SRC]; if (!srcDen) SYSTEM_ERROR;
            corr=(FReal*) new FReal[NOD]; if (!corr) SYSTEM_ERROR;
            
            norm = fem.fmm_normalizer;
            init(fem, VERBOSE);
            if(VERBOSE) { std::cout << "\n>> ScalFMM initialized, using " << ScalfmmNbThreads << " threads.\n" << std::endl; }  
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
        demag<1> (fem,mySettings); // Hd(v), second order contribution
        counter.tac();
        if(mySettings.verbose) { std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl; }
        }
        
    private:
        const int NOD;/**< number of nodes */
        const int FAC;/**< number of facettes */
        const int TET;/**< number of tetrahedrons */
        const int SRC;/**< number of charge sources */
        
        OctreeClass *tree    = nullptr;/**< tree initialized by init private member */
        KernelClass *kernels = nullptr;/**< kernel initialized by init private member */
        
        FReal *srcDen;/**< source buffer */
        FReal *corr;/**< correction coefficients buffer */
        
        double norm;/**< normalization coefficient */
        
        /**
        initialization function for the building of an octree to compute the demag field
        */
        void init(Fem &fem,bool VERBOSE)
            {
            FTic counter;

            if(VERBOSE)
                { std::cout << "Creating & Inserting particles ...\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl; }
            counter.tic();

            Pt::pt3D c = fem.c;
            FSize idxPart=0;

            std::for_each(fem.node.begin(),fem.node.end(),[this,c,&idxPart](Nodes::Node const& n)
                {
                Pt::pt3D pTarget = n.p - c; pTarget *= norm;
                const FPoint<FReal> particlePosition(pTarget.x(), pTarget.y(), pTarget.z());
                tree->insert(particlePosition, FParticleTypeTarget, idxPart, 0.0);//ct 1 pour target    
                idxPart++;
                });//end for_each sur node
  
            if(VERBOSE) { std::cout << "Physical nodes inserted." << std::endl; }

            std::for_each(fem.tet.begin(),fem.tet.end(),[this,c,&idxPart](Tetra::Tet const& tet)              
                {       // sources de volume
                double gauss[DIM][Tetra::NPI];
                tet.interpolation(Nodes::get_p,gauss);
        
                for (int j=0; j<Tetra::NPI; j++, idxPart++)
                    {
                    Pt::pt3D pSource = Pt::pt3D(gauss[0][j],gauss[1][j],gauss[2][j]) - c; pSource *= norm;
                    const FPoint<FReal> particlePosition(pSource.x(), pSource.y(), pSource.z());
                    tree->insert(particlePosition, FParticleTypeSource, idxPart, 0.0);
                    }
                });//end for_each on tet

            if(VERBOSE) { std::cout << "Volume charges inserted." << std::endl; }
    
            std::for_each(fem.fac.begin(),fem.fac.end(),[this,c,&idxPart](Facette::Fac const& fac)
                {        // sources de surface
                double gauss[DIM][Facette::NPI];
                fac.interpolation(Nodes::get_p,gauss);
        
                for (int j=0; j<Facette::NPI; j++, idxPart++)
                    {
                    Pt::pt3D pSource = Pt::pt3D(gauss[0][j],gauss[1][j],gauss[2][j]) - c; pSource *= norm;
                    const FPoint<FReal> particlePosition(pSource.x(), pSource.y(), pSource.z());
                    tree->insert(particlePosition, FParticleTypeSource, idxPart, 0.0);
                    }
                });//end for_each on fac
            counter.tac();
            if(VERBOSE)
                { std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s)." << std::endl; }
            }
        
    /**
    template to computes the demag field, template parameter is either 0 or 1
    */
        template <int Hv>
        void demag(Fem &fem,Settings &settings)
        {
        FmmClass algo(tree, kernels);
        
        //FReal *srcDen=(FReal*) new FReal[SRC]; if (!srcDen) SYSTEM_ERROR;
        
        memset(srcDen, 0, SRC*sizeof(FReal));

        //FReal *corr=(FReal*) new FReal[NOD]; if (!corr) SYSTEM_ERROR;
        memset(corr, 0, NOD*sizeof(FReal));

        int nsrc = 0;
        std::function<Pt::pt3D (Nodes::Node)> getter;

        if(Hv)
            { getter = Nodes::get_v; }
        else
            { getter = Nodes::get_u;}

// *********************** TETRAS ********************
        std::for_each(fem.tet.begin(),fem.tet.end(),[this,getter,&nsrc,&settings](Tetra::Tet const& tet)              
            {
            double Ms = nu0 * settings.paramTetra[tet.idxPrm].J;
            /*---------------- INTERPOLATION ---------------*/
            double dudx[DIM][Tetra::NPI], dudy[DIM][Tetra::NPI], dudz[DIM][Tetra::NPI];
            tet.interpolation(getter,dudx,dudy,dudz);
            /*-----------------------------------------------*/

            for (int j=0; j<Tetra::NPI; j++, nsrc++)
                { srcDen[nsrc] = -Ms * ( dudx[0][j] + dudy[1][j] + dudz[2][j] ) * tet.weight[j]; }
            });//end for_each on tet


//      ************************ FACES **************************
        const bool pot_corr = settings.analytic_corr; 

        std::for_each(fem.fac.begin(),fem.fac.end(),[this,&fem,pot_corr,getter,&nsrc](Facette::Fac const& fac)
            {
            double Ms = fac.Ms;
            Pt::pt3D n = fac.n;
            double u[DIM][Facette::NPI];
    
            fac.interpolation(getter,u);
    
            for (int j=0; j<Facette::NPI; j++, nsrc++)
                { srcDen[nsrc] = Ms * ( u[0][j]*n.x() + u[1][j]*n.y() + u[2][j]*n.z() ) * fac.weight[j]; }

            if (pot_corr)
                {// calc coord gauss
                double gauss[DIM][Facette::NPI];
        
                fac.interpolation(Nodes::get_p,gauss);
                // calc corr node by node
                for (int i=0; i<Facette::N; i++)
                    {
                    int i_ = fac.ind[i];
                    Pt::pt3D p_i_ = fem.node[i_].p;	      
                    for (int j=0; j<Facette::NPI; j++)
                        {
                        Pt::pt3D pg = Pt::pt3D(gauss[Pt::IDX_X][j], gauss[Pt::IDX_Y][j], gauss[Pt::IDX_Z][j]);
                        double rij = Pt::dist(p_i_,pg);
                        double sj = Ms* ( u[0][j]*n.x() + u[1][j]*n.y() + u[2][j]*n.z() );
                        corr[i_]-= sj/rij * fac.weight[j];
                        }
                    corr[i_]+= fac.potential(getter,i);//potential<Hv>(fem.node, fac, i);
                    }
                }
            });// end for_each on fac

        fflush(NULL);

        { // reset potentials and forces - physicalValues[idxPart] = Q
        tree->forEachLeaf([this](LeafClass* leaf)
            {
            const int nbParticlesInLeaf = leaf->getSrc()->getNbParticles();
            const FVector<long long>& indexes = leaf->getSrc()->getIndexes(); // pas int mais long long  *ct*
            FReal* const physicalValues = leaf->getSrc()->getPhysicalValues();
            memset(physicalValues, 0, nbParticlesInLeaf*sizeof(FReal));

            for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
                const int indexPartOrig = indexes[idxPart];
                const int nsrc = indexPartOrig-NOD;
                assert((nsrc>=0) && (nsrc<SRC));
                physicalValues[idxPart]=srcDen[nsrc];
                }
            });

        tree->forEachLeaf([](LeafClass* leaf){
        FReal*const potentials = leaf->getTargets()->getPotentials();
        const int nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
        memset(potentials, 0, nbParticlesInLeaf*sizeof(FReal));
        });
    
        tree->forEachCell([](CellClass* cell){ cell->resetToInitialState(); });
        }// end reset

        algo.execute();

        tree->forEachLeaf([this,&fem](LeafClass* leaf){
            const FReal*const potentials = leaf->getTargets()->getPotentials();
            const int nbParticlesInLeaf  = leaf->getTargets()->getNbParticles();
            //const FVector<int>& indexes  = leaf->getTargets()->getIndexes(); // *ct*
            const FVector<long long>& indexes  = leaf->getTargets()->getIndexes(); // int -> long long *ct*

            for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart)
                {
                const int indexPartOrig = indexes[idxPart];

                if (Hv) 
                    fem.node[indexPartOrig].phiv = (potentials[idxPart]*norm + corr[indexPartOrig])/(4*M_PI);
                else 
                    fem.node[indexPartOrig].phi  = (potentials[idxPart]*norm + corr[indexPartOrig])/(4*M_PI);
                }
            });

        //delete [] srcDen;
        //delete [] corr;
        }
    
};//end class fmm
    
}//end namespace
#endif
