/** \file fem.h
\brief principal header, contains the structs node fac tet , and fem <br>
This file is called by almost all source files in feellgood project <br>
It does also contains the definition of many constants for the solver, and for scalfmm 
*/

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <list>

#include <time.h>
#include <sys/times.h>
#include <unistd.h>

/* only pull GSL if needed */
#ifdef STAT
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#endif

#include "ANN.h" // ct "ANN/ANN.h"		// ANN declarations
#include "gmm/gmm_kernel.h"  // ct gmm_kernel.h plutôt que gmm.h , qui appelle des fichiers de getfem
#include "gmm/gmm_precond_diagonal.h" //ct

// scalFMM includes

#include "Utils/FTic.hpp"
#include "Utils/FParameters.hpp"

#include "Components/FTypedLeaf.hpp"
#include "Components/FParticleType.hpp"

#include "Containers/FOctree.hpp"
#include "Containers/FVector.hpp"

#include "Core/FFmmAlgorithmThreadTsm.hpp"

/*
#include "Kernels/Spherical/FSphericalKernel.hpp"
#include "Kernels/Spherical/FSphericalCell.hpp"
*/

#include "Kernels/P2P/FP2PParticleContainerIndexed.hpp"

#include "Kernels/Rotation/FRotationKernel.hpp"
#include "Kernels/Rotation/FRotationCell.hpp"

/*
#include "Kernels/Chebyshev/FChebCell.hpp"
#include "Kernels/Chebyshev/FChebMatrixKernel.hpp"
#include "Kernels/Chebyshev/FChebKernel.hpp"
#include "Kernels/Chebyshev/FChebSymKernel.hpp"
*/

using namespace std;

/* macros for messages and errors */
#ifdef LIBRARY
#include <stdexcept>  // for runtime_error
#include <cerrno>
#include <system_error>
const pair<string,int> VERBOSE_KEY {"verbose", -1};
#define IF_VERBOSE(fem) if ((fem).param[VERBOSE_KEY] != 0) /**< macro for the definition of a verbose mode if feellgood compiled as a library */
#define SYSTEM_ERROR {throw system_error(errno,generic_category());}/**< macro for error handling if feellgood compiled as a library */
#else
#define IF_VERBOSE(fem) /**< macro for the definition of a verbose mode if feellgood compiled as an executable  */
#define SYSTEM_ERROR {exit(1);} /**< macro to exit executable on some errors */
#endif

const int D = 3;         /**< dimension, required by ANN for kdTree */
const double invD = 1./(double)D;/**< convenient const value */
const size_t NCLASSES=100;/**< number of BIN for statistics */
const double HMAX=1e6;/**< maximum field value ?? (to check) */

/* for statics 
const double EPSILON = 1e-40;
const double DUMAX   = 0.1;
const double DTMIN   = 1e-14;
const double DTMAX   = 1.e-5;
*/

/* for dynamics */

/** \f$ \epsilon \f$ is smallest possible value (for what value ? to check) */
const double EPSILON = 1e-40;

/** min for du step */
const double DUMIN   = 1e-9;	//1e-6

/** max for du step */
const double DUMAX   = 0.02;	// 0.02

/** minimum step time for time integrator */
const double DTMIN   = 1e-14;

/** maximum step time for time integrator */
const double DTMAX   = 1e-7;	// 1e-7

/** no idea */
const double TAUR    = 100.*DTMAX;


const double mu0 = 4.*M_PI*1e-7;/**< \f$ \mu_0 = 4 \pi 10^{-7} \f$ */
const double nu0 = 1./mu0;/**< \f$ \nu_0 = 1/\mu_0 \f$ */

static const int P = 9;/**< constant parameter for some scalfmm templates */

// pb avec ces templates : ils prennent un argument de plus en 1.5.0 //ct
//typedef FTypedRotationCell<P>            CellClass; ct
typedef FTypedRotationCell<double,P>            CellClass; /**< convenient typedef for the definition of cell type in scalfmm  */

//typedef FP2PParticleContainerIndexed<>         ContainerClass; //ct
typedef FP2PParticleContainerIndexed<double>         ContainerClass; /**< convenient typedef for the definition of container for scalfmm */

typedef FTypedLeaf<double, ContainerClass >                      LeafClass;/**< convenient typedef for the definition of leaf for scalfmm  */

typedef FOctree<double, CellClass, ContainerClass , LeafClass >  OctreeClass;/**< convenient typedef for the definition of the octree for scalfmm */

typedef FRotationKernel<double, CellClass, ContainerClass, P >          KernelClass;/**< convenient typedef for the kernel for scalfmm */

typedef FFmmAlgorithmThreadTsm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;/**< convenient typedef for handling altogether the differents scalfmm object templates used in feellgood  */

typedef gmm::wsvector <double>   write_vector;/**< convenient macro */
typedef gmm::rsvector <double>   read_vector;/**< convenient macro */

typedef gmm::row_matrix	<write_vector>   write_matrix;/**< convenient macro */
typedef gmm::row_matrix	<read_vector>    read_matrix;/**< convenient macro */

typedef double triple[3];/**< a 3D point */

/** \struct Node
Node is containing physical point of coordinates \f$ (x,y,z) \f$, magnetization value at p. 
Many other values for the computation of the scalar potential \f$ \phi \f$
*/
struct Node {
    double x;/**< Physical position x of the node */
double y;/**< Physical position y  of the node */
double z;/**< Physical position z  of the node */
    triple u0;/**< magnetization initial value ? */
triple v0;/**< no idea */
triple u;/**< magnetization value */
triple v;/**< no idea */
    triple ep;/**< no idea */
triple eq;/**< no idea */
    double phi0;/**< scalar potential initial value ? */
double phi;/**< scalar potential value */
    double phiv0;/**< no idea */
double phiv;/**< no idea */
    };

/** \struct Fac
Face is a struct containing the index references to nodes, it has a triangular shape and should not be degenerated 
*/
struct Fac{
	static const int N = 3; /**< number of sommits */
	static const int NPI = 4; /**< number of weights  */
	int reg;/**< .msh region number */
	double surf; /**< surface of the face */
	double Ms; /**< magnetization at saturation of the face */    
	double nx;/**< x component of the normal vector */
	double ny;/**< y component of the normal vector */
	double nz;/**< z component of the normal vector */
	int ind[N];/**< indices table */
	double weight[NPI];/**< weights table */
	double a[N][NPI];          /**< hat functions table */
    };

/** \struct Tet
Tet is a tetrahedron, containing the index references to nodes, must not be flat 
   */ 
struct Tet{
static const int N = 4;/**< number of sommits */
static const int NPI = 5;/**< number of weights  */
int reg;/**< .msh region number */
double vol;/**< volume of the tetrahedron */
int ind[N];/**< indices to the nodes */
double weight[NPI];/**< weights */
double a[N][NPI];/**< hat functions */
double dadx[N][NPI];/**< variations of hat function along x directions */
double dady[N][NPI];/**< variations of hat function along y directions */
double dadz[N][NPI];/**< variations of hat function along z directions */
};

/** \struct Seq
Seq describe a sequence of field from \f$ B_{ini} \f$ to \f$ B_{fin} \f$ by steps \f$ dB \f$
*/    
struct Seq{
    double Bini;/**< starting value */
double Bfin;/**< ending value */
double dB;/**< step field */
    triple a;/**< direction (should be normalized) */
    };

/** \struct Stat
used to build some statistics, with GSL library
*/
struct Stat{
#ifdef STAT
    gsl_histogram* h;/**< pointer to GSL histogram */
#else
    void* h;  /**< placeholder, we don't include GSL */
#endif
    double M; /**< M defined by \f$ M = \alpha/ \tau_R |log(dt/ \tau_R)| \f$ */
    double R; /**< R defined by \f$ R = dt/\tau_R*|log(dt/\tau_R)| \f$ */
    double r; /**< r defined by \f$ r = 0.5*dt/\tau_R*|log(dt/\tau_R)| \f$ */
    };

/** \struct Fem
massive container to grab altogether all parameters of a simulation, including mesh geometry, containers for the mesh
*/
struct Fem{
	string  simname;/**< simulation name */
	string pbname;     /**< file pro */
	int REG; /**< number of regions in the msh (to check) */
	int NOD;/**< number of nodes in the corresponding container */
	int FAC;/**< number of faces in the corresponding container */
	int TET;/**< number of tetrahedron in the corresponding container */
	int SRC;/**< number of sources for scalfmm */
	int SEQ;/**< number of sequences ? (to check) */
	double cx;/**< x axis of what? */
	double cy;/**< y axis of what? */
	double cz;/**< z axis of what? */
	double lx;/**< length along x axis (to check) */
	double ly;/**< length along y axis (to check) */
	double lz;/**< length along z axis (to check) */
	double diam;/**< diameter of the mesh (if a wire) */
	double surf;/**< total surface */
	double vol;/**< total volume of the mesh */
    double scale;/**< scaling factor from gmsh files to feellgood */
double fmm_normalizer;/**< no idea; probably normalizing constant for the computation of the demag field */
    double as[3];/**< no idea */
double locmax;/**< no idea */
    double t;/**< physical current time of the simulation */
double tf;/**< end time of the simulation */
double dt;/**< step time of the simulation */
double vmax;/**< maximum speed of what ? */

double E0[4];/**< table to store initial energy values <br> 
index convention : 0-exchange 1-anisotropy 2-demagnetizing 3-applied */

double E[4]; /**< table to store energy values at time t <br>
index convention : 0-exchange 1-anisotropy 2-demagnetizing 3-applied */

double DW_vz0;/**< initial speed of the domain wall (to check) */
double DW_vz;/**< speed of the domain wall along z */
double DW_dir;/**< direction of the domain wall (to check) */
double DW_z; /**< domain wall displacement along Oz */

double Etot0;/**< initial total energy (to check) */
double Etot;/**< total energy */
double evol;/**< no idea */
double phy;/**< no idea */
    vector <Node> node; /**< node container */
    vector <Fac>  fac; /**< face container */
    vector <Tet>  tet; /**< tetrahedron container */
    double Bext;/**< amplitude of the applied field (to check) */
    triple Hext;/**< external applied field direction (should be normalized) */
    gmm::diagonal_precond <read_matrix>  *prc;/**< diagonal preconditionner */
    //gmm::ilutp_precond < read_matrix > *prc;
    map <pair<string,int>,double> param;/**< a map to store different parameters of a simulation, such as suppression of the charges on some surfaces or not, the material parameter values as \f$ Ms \f$ or \f$ \alpha \f$  */
    
Stat stat;/**< to build some histograms  */

    ANNkd_tree* kdtree;/**< a kdtree to find efficiently the closest set of nodes to a physical point in the mesh  */
    ANNpointArray pts;/**< container for the building of the kdtree */
    };

/** \struct Regions
contains two maps of the volume regions and surfaces numbers
*/
struct Regions{
    map <int,int> surfaces;/**< map for the surfaces region generated by gmsh */
    map <int,int> volumes;/**< map for the volumes region generated by gmsh */
    };

double cputime();/**< convenient function to compute the duration of a simulation */
void extract_comment(istream &flux);/**< parser */

/** 
\return square of a number \f$ x^2 \f$
*/
inline double sq(double x) {return x*x;}

void dialog(Fem& fem, vector<Seq> &seq);/**< some text send to terminal to see what step in the sequence the sover is computing */

void lecture(Fem &fem, double scale, Regions *regions);/**< reading file function */

void alloc_nodal(Fem &fem);/**< allocating function of what ? */

void femutil(Fem &fem);/**< utilitary function to affect fem struct */

void femutil_node(Fem &fem);/**< read the nodes from a .msh file to the node container in fem struct */
void femutil_tet(Fem &fem);/**< read the tetrahedron from a .msh file to the tetrahedron container in fem struct */

void femutil_fac(Fem &fem);/**< read the faces from a .msh file to the face container in fem struct */
void femutil_facMs(Fem &fem);/**< no idea */

void chapeaux(Fem &fem);/**< computes the hat functions for all containers */ 
void affichage(Fem &fem);/**< printing function? (to check) */
void direction(Fem &fem);/**< no idea */

void restoresol(Fem& fem, string *filename);/**< read a solution from a file (tsv formated) and initialize fem struct  */

void init_distrib(Fem &fem);/**< computes an analytical initial magnetization distribution as a starting point for the simulation */

// void dirichlet(Fem &fem);
// void periodic(Fem &fem);

void normalize(triple &a);/**< in place normalizing function of vector a */
double u_moy(Fem &fem, int d);/**< computes the average magnetization */
double v_moy(Fem &fem, int d);/**< computes what ? */

void energy(Fem &fem);/**< computes energies */
bool recentrage(Fem &fem, double thres);/**< recentering algorithm for the study of the motion of an object, for example a domain wall. Mesh must be adequate. */

void saver(Fem &fem, ofstream &fout, int nt);/**< saving function for a solution */
void savecfg_vtk(Fem &fem, int nt, string *filename);/**< text file (vtk) writing function for a solution */
void savesol(Fem &fem, int nt, string *filename);/**< text file (tsv) writing function for a solution */
void savestat(Fem &fem, int nt);/**< file writing function for the statistics (if any) */
void saveH(Fem &fem, int nt);/**< save the field values */

int  vsolve(Fem &fem, long nt);/**< solver */
void base_projection(Fem &fem);/**< computes the projection of the llg operators on the elements */

/**
computes the contribution of the tetrahedron to the integrals
*/
void integrales(Fem &fem, Tet &tet, gmm::dense_matrix <double> &AE, vector <double> &BE);

/**
computes the contribution of the surface to the integrals
*/
void integrales(Fem &fem, Fac &fac, gmm::dense_matrix <double> &AE, vector <double> &BE);

/**
projections on the tetrahedrons
*/
void projection(Fem &fem, Tet &elt,
           gmm::dense_matrix <double> &A,  vector <double> &B,
           gmm::dense_matrix <double> &Ap, vector <double> &Bp);

/**
projection on the faces
*/
void projection(Fem &fem, Fac &elt,
           gmm::dense_matrix <double> &A,  vector <double> &B,
           gmm::dense_matrix <double> &Ap, vector <double> &Bp);

/**
matrix assembly with all the contributions of the tetrahedrons
*/
void assemblage(Fem &fem, Tet &elt,
           gmm::dense_matrix <double> &Ke, vector <double> &Le,
           write_matrix &K, write_vector &L);

/**
matrix assembly with all the contributions of the faces
*/
void assemblage(Fem &fem, Fac &elt,
           gmm::dense_matrix <double> &Ke, vector <double> &Le,
           write_matrix &K, write_vector &L);

/**
time evolution solver
*/
void evolution(Fem& fem);

/**
reset the fem struct to restart another simulation (to check)
*/
void reset(Fem& fem);
