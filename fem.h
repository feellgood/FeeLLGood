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
#define IF_VERBOSE(fem) if ((fem).param[VERBOSE_KEY] != 0)
#define SYSTEM_ERROR {throw system_error(errno,generic_category());}
#else
#define IF_VERBOSE(fem)
#define SYSTEM_ERROR {exit(1);}
#endif

const int D = 3;         // dimension
const double invD = 1./(double)D;
const size_t NCLASSES=100;
const double HMAX=1e6;

/* for statics 
const double EPSILON = 1e-40;
const double DUMAX   = 0.1;
const double DTMIN   = 1e-14;
const double DTMAX   = 1.e-5;
*/

/* for dynamics */
const double EPSILON = 1e-40;
const double DUMIN   = 1e-9;	//1e-6
const double DUMAX   = 0.02;	// 0.02
const double DTMIN   = 1e-14;
const double DTMAX   = 1e-7;	// 1e-7
const double TAUR    = 100.*DTMAX;
/*
*/

const double mu0 = 4.*M_PI*1e-7;
const double nu0 = 1./mu0;

static const int P = 9;
// pb avec ce template : il doit prendre deux arguments //ct
//typedef FTypedRotationCell<P>            CellClass; ct
typedef FTypedRotationCell<double,P>            CellClass; //ct

//typedef FP2PParticleContainerIndexed<>         ContainerClass; //ct
typedef FP2PParticleContainerIndexed<double>         ContainerClass; //ct

typedef FTypedLeaf<double, ContainerClass >                      LeafClass;
typedef FOctree<double, CellClass, ContainerClass , LeafClass >  OctreeClass;
typedef FRotationKernel<double, CellClass, ContainerClass, P >          KernelClass;

typedef FFmmAlgorithmThreadTsm<OctreeClass, CellClass, ContainerClass, KernelClass, LeafClass > FmmClass;

typedef gmm::wsvector <double>   write_vector;
typedef gmm::rsvector <double>   read_vector;

typedef gmm::row_matrix	<write_vector>   write_matrix;
typedef gmm::row_matrix	<read_vector>    read_matrix;

typedef double triple[3];

struct Node {
    double x, y, z;
    triple u0,v0,u,v;
    triple ep, eq;
    double phi0, phi;
    double phiv0, phiv;
    };


struct Fac{
    static const int N = 3, NPI = 4;
    int reg;
    double surf, Ms;
    double nx, ny, nz;                // normale ext
    int ind[N];
    double weight[NPI];
    double a[N][NPI];          // fct chapeaux
    };
    
struct Tet{
    static const int N = 4, NPI = 5;
    int reg;
    double vol;
    int ind[N];
    double weight[NPI];
    double a[N][NPI];
    double dadx[N][NPI], dady[N][NPI], dadz[N][NPI];
    };
    
struct Seq{
    double Bini, Bfin, dB;
    triple a;
    };

struct Stat{
#ifdef STAT
    gsl_histogram* h;
#else
    void* h;  // placeholder, we don't include GSL
#endif
    double M; // M := alpha/tauR*|log(dt/tauR)|
    double R; // R :=    dt/tauR*|log(dt/tauR)|
    double r; // r := 0.5*dt/tauR*|log(dt/tauR)|
    };

struct Fem{
    string  simname, pbname;     // nom simul, fichier pro
    int REG, NOD, FAC, TET, SRC, SEQ;
    double cx,cy,cz, lx,ly,lz, diam, surf, vol;
    double scale, fmm_normalizer;
    double as[3], locmax;
    double t, tf, dt, vmax;
    double E0[4], E[4]; //0-exchange 1-anisotropy 2-demagnetizing 3-applied
    double DW_vz0, DW_vz, DW_dir, DW_z; /* vitesse de chgt de referentiel et deplacement selon Oz */
    double Etot0, Etot, evol, phy;
    vector <Node> node;
    vector <Fac>  fac;
    vector <Tet>  tet; 
    double Bext;
    triple Hext;
    gmm::diagonal_precond <read_matrix>  *prc;
    //gmm::ilutp_precond < read_matrix > *prc;
    map <pair<string,int>,double> param;
    Stat stat;

    ANNkd_tree* kdtree;
    ANNpointArray pts;
    };

struct Regions{
    map <int,int> surfaces;
    map <int,int> volumes;
    };

double cputime();
void extract_comment(istream &flux);
inline double sq(double x) {return x*x;}

void dialog(Fem& fem, vector<Seq> &seq);
void lecture(Fem &fem, double scale, Regions *regions);
void alloc_nodal(Fem &fem);

void femutil(Fem &fem);
void femutil_node(Fem &fem);
void femutil_tet(Fem &fem);
void femutil_fac(Fem &fem);
void femutil_facMs(Fem &fem);

void chapeaux(Fem &fem);
void affichage(Fem &fem);
void direction(Fem &fem);

void restoresol(Fem& fem, string *filename);
void init_distrib(Fem &fem);

// void dirichlet(Fem &fem);
// void periodic(Fem &fem);

void normalize(triple &a);
double u_moy(Fem &fem, int d);
double v_moy(Fem &fem, int d);

void energy(Fem &fem);
bool recentrage(Fem &fem, double thres);

void saver(Fem &fem, ofstream &fout, int nt);
void savecfg_vtk(Fem &fem, int nt, string *filename);
void savesol(Fem &fem, int nt, string *filename);
void savestat(Fem &fem, int nt);
void saveH(Fem &fem, int nt);

int  vsolve(Fem &fem, long nt);
void base_projection(Fem &fem);

void integrales(Fem &fem, Tet &tet, gmm::dense_matrix <double> &AE, vector <double> &BE);
void integrales(Fem &fem, Fac &fac, gmm::dense_matrix <double> &AE, vector <double> &BE);

void projection(Fem &fem, Tet &elt,
           gmm::dense_matrix <double> &A,  vector <double> &B,
           gmm::dense_matrix <double> &Ap, vector <double> &Bp);

void projection(Fem &fem, Fac &elt,
           gmm::dense_matrix <double> &A,  vector <double> &B,
           gmm::dense_matrix <double> &Ap, vector <double> &Bp);

void assemblage(Fem &fem, Tet &elt,
           gmm::dense_matrix <double> &Ke, vector <double> &Le,
           write_matrix &K, write_vector &L);

void assemblage(Fem &fem, Fac &elt,
           gmm::dense_matrix <double> &Ke, vector <double> &Le,
           write_matrix &K, write_vector &L);

void evolution(Fem& fem);
void reset(Fem& fem);
