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

#include "config.h" //for NbThreads

/** constant parameter for some scalfmm templates */
static const int P = 9;

/** double redefinition for the parametrization of some scalfmm templates */
#define FReal double

// pb avec ces templates : ils prennent un argument de plus en 1.5.0 et 1.4.148 //ct
//typedef FTypedRotationCell<P>            CellClass; ct
typedef FTypedRotationCell<FReal, P>            CellClass; /**< convenient typedef for the definition of cell type in scalfmm  */

//typedef FP2PParticleContainerIndexed<>         ContainerClass; //ct
typedef FP2PParticleContainerIndexed<FReal>         ContainerClass; /**< convenient typedef for the definition of container for scalfmm */

typedef FTypedLeaf<FReal, ContainerClass >                      LeafClass;/**< convenient typedef for the definition of leaf for scalfmm  */

typedef FOctree< FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;/**< convenient typedef for the definition of the octree for scalfmm */

typedef FRotationKernel< FReal, CellClass, ContainerClass, P >          KernelClass;/**< convenient typedef for the kernel for scalfmm */

typedef FFmmAlgorithmThreadTsm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;/**< convenient typedef for handling altogether the differents scalfmm object templates used in feellgood  */

/**
\namespace fmm to grab altogether the templates and functions using scalfmm for the computation of the demag field 
*/

namespace fmm{

/**
initialization function for the building of an octree and a kernel passed to scalfmm to compute the demag field
*/
template <class CellClass, class ContainerClass, class LeafClass, class OctreeClass, class KernelClass, class FmmClass, typename... Args>
int init(Fem &fem, OctreeClass* &tree, KernelClass* &kernels, Args... kernelPreArgs)
{
    FTic counter;
    const int NbLevels = 8; 
    const int SizeSubLevels = 6;
    //const unsigned int NbThreads  =  omp_get_max_threads();

    omp_set_num_threads(ScalfmmNbThreads);
    if(VERBOSE) { std::cout << "\n>> ScalFMM using " << ScalfmmNbThreads << " threads.\n" << std::endl; }

    const double boxWidth=2.01;//2.01 ct
    const FPoint<FReal> centerOfBox(0., 0., 0.);// manque le typename du template FPoint *ct*

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
    [c,norm,&tree,&idxPart](Node const& n)
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
        double nod[3][Tetra::N], gauss[3][Tetra::NPI];
        for (int i=0; i<Tetra::N; i++)
            {
            int i_= tet.ind[i];
            nod[0][i] = fem.node[i_].p.x();
            nod[1][i] = fem.node[i_].p.y();
            nod[2][i] = fem.node[i_].p.z();
            }
        tiny::mult<double, 3, Tetra::N, Tetra::NPI> (nod, Tetra::a, gauss);

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
        double nod[3][Facette::N], gauss[3][Facette::NPI];
        for (int i=0; i<Facette::N; i++)
            {
            int i_= fac.ind[i];
            nod[0][i] = fem.node[i_].p.x();
            nod[1][i] = fem.node[i_].p.y();
            nod[2][i] = fem.node[i_].p.z();
            }
        tiny::mult<double, 3, Facette::N, Facette::NPI> (nod, Facette::a, gauss);

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
computes correction on potential on facettes using fem struct
*/
template <int Hv> double potential(std::vector<Node> const& myNode, Facette::Fac const& fac, int i);


/**
computes the demag field
*/
template <int Hv, class CellClass, class ContainerClass, class LeafClass, class OctreeClass,
          class KernelClass, class FmmClass, typename... Args>
void demag(Fem &fem,Settings &settings, OctreeClass *tree, KernelClass *kernels, Args... kernelPreArgs)
{
FTic counter;
FmmClass algo(tree, kernels);

if(VERBOSE) { std::cout << "\t magnetostatics ..................... "; }

const int NOD = fem.NOD;
const int FAC = fem.fac.size();
const int TET = fem.tet.size();
const int SRC = FAC * Facette::NPI + TET * Tetra::NPI;

FReal *srcDen=(FReal*) new FReal[SRC]; if (!srcDen) SYSTEM_ERROR;
memset(srcDen, 0, SRC*sizeof(FReal));

FReal *corr=(FReal*) new FReal[NOD]; if (!corr) SYSTEM_ERROR;
memset(corr, 0, NOD*sizeof(FReal));

int nsrc = 0;

/*********************** TETRAS *********************/
std::for_each(fem.tet.begin(),fem.tet.end(),
[&nsrc,&srcDen,&settings,&fem](Tetra::Tet const& tet)              
    {
    double Ms = nu0 * settings.paramTetra[tet.idxPrm].J;
   /*---------------- INTERPOLATION ---------------*/
    double u_nod[3][Tetra::N];
    double dudx[3][Tetra::NPI], dudy[3][Tetra::NPI], dudz[3][Tetra::NPI];

    for (int i=0; i<Tetra::N; i++)
        {
        Node &node = fem.node[ tet.ind[i] ];
        u_nod[Pt::IDX_X][i] = (Hv? node.v.x(): node.u.x());
        u_nod[Pt::IDX_Y][i] = (Hv? node.v.y(): node.u.y());
        u_nod[Pt::IDX_Z][i] = (Hv? node.v.z(): node.u.z());        
        }

    tiny::mult<double, 3, Tetra::N, Tetra::NPI> (u_nod, tet.dadx, dudx);
	tiny::mult<double, 3, Tetra::N, Tetra::NPI> (u_nod, tet.dady, dudy);
	tiny::mult<double, 3, Tetra::N, Tetra::NPI> (u_nod, tet.dadz, dudz);
   /*-----------------------------------------------*/

    for (int j=0; j<Tetra::NPI; j++, nsrc++){
        double div_u = dudx[0][j] + dudy[1][j] + dudz[2][j];
        srcDen[nsrc] = -Ms * div_u * tet.weight[j];
        }
    });//end for_each on tet


/************************ FACES **************************/
const bool pot_corr = settings.analytic_corr; 

std::for_each(fem.fac.begin(),fem.fac.end(),
[pot_corr,&nsrc,&srcDen,&corr,&fem](Facette::Fac const& fac)
    {
    double Ms = fac.Ms;
    Pt::pt3D n = fac.n;
        /** calc u gauss **/  
    double u_nod[3][Facette::N], u[3][Facette::NPI];
    for (int i=0; i<Facette::N; i++)
        {
        Node &node = fem.node[ fac.ind[i] ];
        u_nod[Pt::IDX_X][i] = (Hv? node.v.x(): node.u.x());
        u_nod[Pt::IDX_Y][i] = (Hv? node.v.y(): node.u.y());
        u_nod[Pt::IDX_Z][i] = (Hv? node.v.z(): node.u.z());
        }

    tiny::mult<double, 3, Facette::N, Facette::NPI> (u_nod, Facette::a, u);

        /** calc sigma, fill distrib.alpha **/
    for (int j=0; j<Facette::NPI; j++, nsrc++)
        {
        double un = u[0][j]*n.x() + u[1][j]*n.y() + u[2][j]*n.z();
        srcDen[nsrc] = Ms * un * fac.weight[j];
        }

    if (pot_corr)
        {/** calc coord gauss **/
        double nod[3][Facette::N], gauss[3][Facette::NPI];
        for (int i=0; i<Facette::N; i++)
            {
            int i_= fac.ind[i];
            nod[0][i] = fem.node[i_].p.x();
            nod[1][i] = fem.node[i_].p.y();
            nod[2][i] = fem.node[i_].p.z();
            }
        tiny::mult<double, 3, Facette::N, Facette::NPI> (nod, Facette::a, gauss);

      /** calc corr node by node **/
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
        corr[i_]+= potential<Hv>(fem.node, fac, i);
        }
      }
   });// end for_each on fac

fflush(NULL);

    { // reset potentials and forces - physicalValues[idxPart] = Q

    tree->forEachLeaf([&](LeafClass* leaf)
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

    tree->forEachLeaf([&](LeafClass* leaf){
	FReal*const potentials = leaf->getTargets()->getPotentials();
	const int nbParticlesInLeaf = leaf->getTargets()->getNbParticles();
	memset(potentials, 0, nbParticlesInLeaf*sizeof(FReal));
        });


    tree->forEachCell([&](CellClass* cell){ cell->resetToInitialState(); });

    }// end reset

counter.tic();
algo.execute();
 counter.tac();
if(VERBOSE) { std::cout << "Done  " << "(@Algorithm = " << counter.elapsed() << "s)." << std::endl; }

double norm = fem.fmm_normalizer;

    tree->forEachLeaf([&](LeafClass* leaf){
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

template <int Hv>
double potential(std::vector<Node> const& myNode, Facette::Fac const& fac, int i) // template, mais Hv est utilisé comme un booleen 
{
  double Ms = fac.Ms;
  Pt::pt3D n = fac.n;

 int ii  = (i+1)%3;
 int iii = (i+2)%3;

 int i_,ii_,iii_;
 i_=fac.ind[i];  ii_=fac.ind[ii];  iii_=fac.ind[iii];

Node const& node1 = myNode[i_];
Node const& node2 = myNode[ii_];
Node const& node3 = myNode[iii_];

Pt::pt3D p1p2 = node2.p - node1.p;
Pt::pt3D p1p3 = node3.p - node1.p;

//double b = sqrt( sq(x2-x1) + sq(y2-y1) + sq(z2-z1) );
double b = p1p2.norm();

//double t = (x2-x1)*(x3-x1) + (y2-y1)*(y3-y1) + (z2-z1)*(z3-z1);
double t = Pt::pScal(p1p2,p1p3);
double h = 2.*fac.surf;
 t/=b;  h/=b;
 double a = t/h;  double c = (t-b)/h;

double s1, s2, s3;
if (Hv) {
	s1 = Pt::pScal(node1.v,n);//node1.v[0]*nx + node1.v[1]*ny + node1.v[2]*nz;
	s2 = Pt::pScal(node2.v,n);//node2.v[0]*nx + node2.v[1]*ny + node2.v[2]*nz;
	s3 = Pt::pScal(node3.v,n);//node3.v[0]*nx + node3.v[1]*ny + node3.v[2]*nz;
   }
else {
	s1 = Pt::pScal(node1.u,n);//node1.u[0]*nx + node1.u[1]*ny + node1.u[2]*nz;
	s2 = Pt::pScal(node1.u,n);//node2.u[0]*nx + node2.u[1]*ny + node2.u[2]*nz;
	s3 = Pt::pScal(node1.u,n);//node3.u[0]*nx + node3.u[1]*ny + node3.u[2]*nz;
   }

 double l = s1;
 double j = (s2-s1)/b;
 double k = t/b/h*(s1-s2) + (s3-s1)/h;

 double cc1 = c*c+1;
 double r = sqrt(h*h + (c*h+b)*(c*h+b));
 double ll = log( (cc1*h + c*b + sqrt(cc1)*r) / (b*(c+sqrt(cc1))) );

 double pot1, pot2, pot3, pot;
 pot1 = b*b/pow(cc1,1.5)*ll + c*b*r/cc1 + h*r - c*b*b/cc1 - sqrt(a*a+1)*h*h;
 pot1*= j/2.;

 pot2 = -c*b*b/pow(cc1,1.5)*ll + b*r/cc1 - h*h/2. + h*h*log(c*h+b+r) - b*b/cc1;
 pot2*= k/2.;

 pot3 = h*log(c*h+b+r) - h + b/sqrt(cc1)*ll;
 pot3*= l;

 pot = pot1 + pot2 + pot3 + h*(k*h/2.+l)*(1-log(h*(a+sqrt(a*a+1)))) - k*h*h/4.;

 return Ms*pot;
}
}
#endif
