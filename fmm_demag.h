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


#ifndef FMM_DEMAG_H
#define FMM_DEMAG_H

static const int P = 9;/**< constant parameter for some scalfmm templates */

/** double redefinition for the parametrization of some scalfmm templates */
#define FReal double


/**
replace enum FParticleType::FParticleTypeSource
*/
#define typeSource 0

/**
replace enum FParticleType::FParticleTypeTarget
*/
#define typeTarget 1


// pb avec ces templates : ils prennent un argument de plus en 1.5.0 //ct
//typedef FTypedRotationCell<P>            CellClass; ct
typedef FTypedRotationCell<FReal, P>            CellClass; /**< convenient typedef for the definition of cell type in scalfmm  */

//typedef FP2PParticleContainerIndexed<>         ContainerClass; //ct
typedef FP2PParticleContainerIndexed<FReal>         ContainerClass; /**< convenient typedef for the definition of container for scalfmm */

typedef FTypedLeaf<FReal, ContainerClass >                      LeafClass;/**< convenient typedef for the definition of leaf for scalfmm  */

typedef FOctree< FReal, CellClass, ContainerClass , LeafClass >  OctreeClass;/**< convenient typedef for the definition of the octree for scalfmm */

typedef FRotationKernel< FReal, CellClass, ContainerClass, P >          KernelClass;/**< convenient typedef for the kernel for scalfmm */

typedef FFmmAlgorithmThreadTsm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;/**< convenient typedef for handling altogether the differents scalfmm object templates used in feellgood  */


#ifndef VERBOSE  // same goes for SYSTEM_ERROR
#error "fem.h" must be #included before "fmm_demag.h"
#endif

/**
\namespace fmm to grab altogether the templates and functions using scalfmm for the computation of the demag field 
*/
namespace fmm{
template <class CellClass, class ContainerClass, class LeafClass, class OctreeClass,
          class KernelClass, class FmmClass, typename... Args>

/**
initialization function for the building of an octree and a kernel passed to scalfmm to compute the demag field
*/
int init(Fem &fem, OctreeClass* &tree, KernelClass* &kernels, Args... kernelPreArgs)
{
    FTic counter;
    const int NbLevels = 8; 
    const int SizeSubLevels = 6;
    const unsigned int NbThreads  =  omp_get_max_threads(); // open mp function

    omp_set_num_threads(NbThreads);
    if(VERBOSE) { std::cout << "\n>> Using " << NbThreads << " threads.\n" << std::endl; }

    const double boxWidth=2.01;//2.01 ct

//const FPoint centerOfBox(0., 0., 0.);
    const FPoint<FReal> centerOfBox(0., 0., 0.);// manque le typename du template FPoint *ct*

    // -----------------------------------------------------
    tree=new OctreeClass(NbLevels, SizeSubLevels, boxWidth, centerOfBox);
    if (!tree) SYSTEM_ERROR;
    // -----------------------------------------------------

fem.fmm_normalizer = 2./fem.diam*0.999999; // pourquoi 2 ?
double norm= fem.fmm_normalizer;

const int NOD = fem.NOD;
const int FAC = fem.FAC;
const int TET = fem.TET;

    if(VERBOSE){
    std::cout << "Creating & Inserting particles ..." << std::endl;
    std::cout << "\tHeight : " << NbLevels << " \t sub-height : " << SizeSubLevels << std::endl;
    }
    counter.tic();

int idxPart=0;
for (int i=0; i<NOD; i++, idxPart++){       // cibles (noeuds)
    double xTarget = (fem.node[i].p.x()-fem.c.x()) * norm;
    double yTarget = (fem.node[i].p.y()-fem.c.y()) * norm;
    double zTarget = (fem.node[i].p.z()-fem.c.z()) * norm;
    const FPoint<FReal> particlePosition(xTarget, yTarget, zTarget);// manque le typename du template FPoint ct
    //tree->insert(particlePosition, FParticleType::FParticleTypeTarget, idxPart, 0.0);//ct
	tree->insert(particlePosition, typeTarget, idxPart, 0.0);//ct 1 pour target
    }

    std::cout << "Physical nodes inserted." << std::endl;
    
for (int t=0; t<TET; t++){       // sources de volume
    Tetra::Tet &tet = fem.tet[t];
    double nod[3][Tetra::N], gauss[3][Tetra::NPI];
    for (int i=0; i<Tetra::N; i++){
        int i_= tet.ind[i];
        nod[0][i] = fem.node[i_].p.x();
        nod[1][i] = fem.node[i_].p.y();
		nod[2][i] = fem.node[i_].p.z();
	}
    tiny::mult<double, 3, Tetra::N, Tetra::NPI> (nod, tet.a, gauss);

    for (int j=0; j<Tetra::NPI; j++, idxPart++){
        double xSource = (gauss[0][j]-fem.c.x()) * norm;
        double ySource = (gauss[1][j]-fem.c.y()) * norm;
		double zSource = (gauss[2][j]-fem.c.z()) * norm;
        const FPoint<FReal> particlePosition(xSource, ySource, zSource);// manque le typename du template FPoint ct
        //tree->insert(particlePosition, FParticleType::FParticleTypeSource, idxPart, 0.0);//ct
	tree->insert(particlePosition, typeSource, idxPart, 0.0);//ct 0 pour source
	}
    }

for (int f=0; f<FAC; f++){        // sources de surface
    Facette::Fac &fac = fem.fac[f];
    double nod[3][Facette::N], gauss[3][Facette::NPI];
    for (int i=0; i<Facette::N; i++){
        int i_= fac.ind[i];
        nod[0][i] = fem.node[i_].p.x();
        nod[1][i] = fem.node[i_].p.y();
		nod[2][i] = fem.node[i_].p.z();
	}
    tiny::mult<double, 3, Facette::N, Facette::NPI> (nod, fac.a, gauss);

    for (int j=0; j<Facette::NPI; j++, idxPart++){
        double xSource = (gauss[0][j]-fem.c.x()) * norm;
        double ySource = (gauss[1][j]-fem.c.y()) * norm;
		double zSource = (gauss[2][j]-fem.c.z()) * norm;
        const FPoint<FReal> particlePosition(xSource, ySource, zSource);// manque le typename du template FPoint ct
        //tree->insert(particlePosition, FParticleType::FParticleTypeSource, idxPart, 0.0);//ct
	tree->insert(particlePosition, typeSource, idxPart, 0.0);//ct 0 pour source
	}
    }
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
template <int Hv> double potential(Fem &fem, Facette::Fac &fac, int i);


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
const int FAC = fem.FAC;
const int TET = fem.TET;
const int SRC = fem.SRC;

FReal *srcDen=(FReal*) new FReal[SRC]; if (!srcDen) SYSTEM_ERROR;
memset(srcDen, 0, SRC*sizeof(FReal));

FReal *corr=(FReal*) new FReal[NOD]; if (!corr) SYSTEM_ERROR;
memset(corr, 0, NOD*sizeof(FReal));

int nsrc = 0;

/*********************** TETRAS *********************/
for (int t=0; t<TET; t++){
    Tetra::Tet &tet = fem.tet[t];
    //double Ms = nu0 * settings.param[ std::make_pair("Js",tet.reg) ];
    double Ms = nu0 * settings.paramTetra[tet.idxPrm].J;
   /*---------------- INTERPOLATION ---------------*/
    double u_nod[3][Tetra::N];
    double dudx[3][Tetra::NPI], dudy[3][Tetra::NPI], dudz[3][Tetra::NPI];

    for (int i=0; i<Tetra::N; i++) {
        int i_= tet.ind[i];
        Node &node = fem.node[i_];
        //for (int d=0; d<3; d++) u_nod[d][i] = (Hv? node.v[d]: node.u[d]);
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
        double q = -Ms * div_u * tet.weight[j];

        srcDen[nsrc] = q;
        }
 }


/************************ FACES **************************/
for (int f=0; f<FAC; f++){
    Facette::Fac &fac = fem.fac[f];
    //const int N    = Fac::N;
    //const int NPI  = Fac::NPI;
    double Ms = fac.Ms;
    Pt::pt3D n = fac.n;//nx=fac.nx; ny=fac.ny; nz=fac.nz;

    /** calc u gauss **/  
    double u_nod[3][Facette::N], u[3][Facette::NPI];
        for (int i=0; i<Facette::N; i++){
        int i_= fac.ind[i];
        Node &node = fem.node[i_];
        //for (int d=0; d<3; d++) u_nod[d][i] = (Hv? node.v[d]: node.u[d]);
	u_nod[Pt::IDX_X][i] = (Hv? node.v.x(): node.u.x());
	u_nod[Pt::IDX_Y][i] = (Hv? node.v.y(): node.u.y());
	u_nod[Pt::IDX_Z][i] = (Hv? node.v.z(): node.u.z());
        }

    tiny::mult<double, 3, Facette::N, Facette::NPI> (u_nod, fac.a, u);

    /** calc sigma, fill distrib.alpha **/
    for (int j=0; j<Facette::NPI; j++, nsrc++){
        double un = u[0][j]*n.x() + u[1][j]*n.y() + u[2][j]*n.z();
        double s = Ms * un * fac.weight[j];
        srcDen[nsrc] =  s; 
        }

    if (settings.analytic_corr) {
      /** calc coord gauss **/
      double nod[3][Facette::N], gauss[3][Facette::NPI];
      for (int i=0; i<Facette::N; i++) {
	      int i_= fac.ind[i];
	      //Node &node = fem.node[i_]; // inutilisé *ct*
	      nod[0][i] = fem.node[i_].p.x();
	      nod[1][i] = fem.node[i_].p.y();
	      nod[2][i] = fem.node[i_].p.z();
          }
      tiny::mult<double, 3, Facette::N, Facette::NPI> (nod, fac.a, gauss);

      /** calc corr node by node **/
      for (int i=0; i<Facette::N; i++) {
	      int i_= fac.ind[i];
	      
		//Node &node = fem.node[i_];
		Pt::pt3D p_i_ = fem.node[i_].p;	      
		double x,y,z;
	      //x=node.p.x();  y=node.p.y();  z=node.p.z();
		x=p_i_.x(); y=p_i_.y(); z=p_i_.z();	      
		for (int j=0; j<Facette::NPI; j++) {
	          double xg,yg,zg;
	          xg=gauss[0][j];  yg=gauss[1][j];  zg=gauss[2][j];
	          double rij = sqrt( sq(x-xg) + sq(y-yg) + sq(z-zg) );
	          double sj = Ms* ( u[0][j]*n.x() + u[1][j]*n.y() + u[2][j]*n.z() );
	          corr[i_]-= sj/rij * fac.weight[j];
	          }
	      corr[i_]+= potential<Hv>(fem, fac, i);
          }
      }
   }

fflush(NULL);

    { // reset potentials and forces - physicalValues[idxPart] = Q

    tree->forEachLeaf([&](LeafClass* leaf){
	const int nbParticlesInLeaf = leaf->getSrc()->getNbParticles();
	//const FVector<int>& indexes = leaf->getSrc()->getIndexes(); // *ct*
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


    tree->forEachCell([&](CellClass* cell){
	cell->resetToInitialState();
        });

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
double potential(Fem &fem, Facette::Fac &fac, int i) // template, mais Hv est utilisé comme un booleen 
{
  double Ms = fac.Ms;
  Pt::pt3D n = fac.n;

 int ii  = (i+1)%3;
 int iii = (i+2)%3;

 int i_,ii_,iii_;
 i_=fac.ind[i];  ii_=fac.ind[ii];  iii_=fac.ind[iii];

 //double x1,x2,x3, y1,y2,y3, z1,z2,z3;
Node &node1 = fem.node[i_];
Node &node2 = fem.node[ii_];
Node &node3 = fem.node[iii_];
 //x1=node1.x;  x2=node2.x;  x3=node3.x;
 //y1=node1.y;  y2=node2.y;  y3=node3.y;
 //z1=node1.z;  z2=node2.z;  z3=node3.z;

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
