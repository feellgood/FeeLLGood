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
    
    const int NbLevels = 8; 
    const int SizeSubLevels = 6;
    const double boxWidth=2.01;//2.01 ct
    const FPoint<FReal> centerOfBox(0., 0., 0.);// manque le typename du template FPoint *ct*

class fmm 
    {
    public:
        inline fmm() {}
        
    private:
        
    };
    
/**
initialization function for the building of an octree and a kernel passed to scalfmm to compute the demag field
*/
template <class CellClass, class ContainerClass, class LeafClass, class OctreeClass, class KernelClass, class FmmClass, typename... Args>
int init(Fem &fem,bool VERBOSE,const int ScalfmmNbThreads, OctreeClass* &tree, KernelClass* &kernels, Args... kernelPreArgs)
{
    FTic counter;
    
    omp_set_num_threads(ScalfmmNbThreads);
    if(VERBOSE) { std::cout << "\n>> ScalFMM using " << ScalfmmNbThreads << " threads.\n" << std::endl; }
    
    // -----------------------------------------------------
    tree=new OctreeClass(NbLevels, SizeSubLevels, boxWidth, centerOfBox);
    if (!tree) SYSTEM_ERROR;
    // -----------------------------------------------------

fem.fmm_normalizer = 1./(2.*fem.diam);
double norm= fem.fmm_normalizer;

if(VERBOSE){
    std::cout << "fmm_normalizer = " << norm << std::endl;
    std::cout << "Creating & Inserting particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    }
    counter.tic();

    Pt::pt3D c = fem.c;
    FSize idxPart=0;

    std::for_each(fem.node.begin(),fem.node.end(),
    [c,norm,&tree,&idxPart](Nodes::Node const& n)
        {
        Pt::pt3D pTarget = n.p - c; pTarget *= norm;
        const FPoint<FReal> particlePosition(pTarget.x(), pTarget.y(), pTarget.z());
        tree->insert(particlePosition, FParticleTypeTarget, idxPart, 0.0);//ct 1 pour target    
        idxPart++;
        });//end for_each sur node
  
    if(VERBOSE) { std::cout << "Physical nodes inserted." << std::endl; }

    std::for_each(fem.tet.begin(),fem.tet.end(),
    [&fem,c,norm,&tree,&idxPart](Tetra::Tet const& tet)              
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
    
    std::for_each(fem.fac.begin(),fem.fac.end(),
    [&fem,c,norm,&tree,&idxPart](Facette::Fac const& fac)
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
    if(VERBOSE){
    std::cout << "Done  " << "(@Creating and Inserting Particles = " << counter.elapsed() << "s).\nCreate kernel ..." << std::endl;
    }
    counter.tic();

    kernels=new KernelClass( kernelPreArgs... , NbLevels, boxWidth, centerOfBox);
    if (!kernels) SYSTEM_ERROR;

    counter.tac();
    if(VERBOSE) { std::cout << "Done  " << " in " << counter.elapsed() << "s)." << std::endl; }

return 0;
}

/**
template to computes the demag field, first template parameter is either 0 or 1
*/
template <int Hv, class CellClass, class ContainerClass, class LeafClass, class OctreeClass,
          class KernelClass, class FmmClass, typename... Args>
void demag(Fem &fem,Settings &settings, OctreeClass *tree, KernelClass *kernels, Args... kernelPreArgs)
{
FmmClass algo(tree, kernels);

const int NOD = fem.node.size();
const int FAC = fem.fac.size();
const int TET = fem.tet.size();
const int SRC = FAC * Facette::NPI + TET * Tetra::NPI;

FReal *srcDen=(FReal*) new FReal[SRC]; if (!srcDen) SYSTEM_ERROR;
memset(srcDen, 0, SRC*sizeof(FReal));

FReal *corr=(FReal*) new FReal[NOD]; if (!corr) SYSTEM_ERROR;
memset(corr, 0, NOD*sizeof(FReal));

int nsrc = 0;
std::function<Pt::pt3D (Nodes::Node)> getter;

if(Hv)
    { getter = Nodes::get_v; }
else
    { getter = Nodes::get_u;}

/*********************** TETRAS *********************/
std::for_each(fem.tet.begin(),fem.tet.end(),[getter,&nsrc,&srcDen,&settings](Tetra::Tet const& tet)              
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

std::for_each(fem.fac.begin(),fem.fac.end(),[pot_corr,getter,&nsrc,&srcDen,&corr,&fem](Facette::Fac const& fac)
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
    tree->forEachLeaf([NOD,SRC,&srcDen](LeafClass* leaf)
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

double norm = fem.fmm_normalizer;

    tree->forEachLeaf([&fem,norm,corr](LeafClass* leaf){
        const FReal*const potentials = leaf->getTargets()->getPotentials();
        const int nbParticlesInLeaf  = leaf->getTargets()->getNbParticles();
        //const FVector<int>& indexes  = leaf->getTargets()->getIndexes(); // *ct*
	const FVector<long long>& indexes  = leaf->getTargets()->getIndexes(); // int -> long long *ct*

        for(int idxPart = 0 ; idxPart < nbParticlesInLeaf ; ++idxPart){
	    const int indexPartOrig = indexes[idxPart];

	    if (Hv) 
       		fem.node[indexPartOrig].phiv = (potentials[idxPart]*norm + corr[indexPartOrig])/(4*M_PI);
	    else 
       		fem.node[indexPartOrig].phi  = (potentials[idxPart]*norm + corr[indexPartOrig])/(4*M_PI);
        }
    });

delete [] srcDen;
delete [] corr;
}

}//end namespace
#endif
