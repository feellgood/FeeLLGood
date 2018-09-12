#ifndef tetra_h
#define tetra_h

/** \file tetra.h
  \brief namespace Tetra
  header containing Tet class, some constants, and integrales
 */

#include "gmm/gmm_kernel.h" // pour dense_matrix dans namespace Tetra

//#include "feellgoodSettings.h"

#include "config.h"

#include "node.h"

typedef double triple[DIM];/**< a 3D point */

/** 
\return square of a number \f$ x^2 \f$
*/
inline double sq(double x /**< [in] */ ) {return x*x;}

/**
in place normalizing function of triple
a */
inline void normalize(triple &a /**< [in,out] */)
{
double norme=sqrt(sq(a[0])+sq(a[1])+sq(a[2]));
a[0]/= norme;a[1]/= norme;a[2]/= norme;
}


/** \namespace Tetra
 to grab altogether some constants for struct Tet
 */
namespace Tetra
{
const int N = 4;/**< number of sommits */
const int NPI = 5;/**< number of weights  */

const double A=1./4.;/**< some constant to build hat functions */
const double B=1./6.;/**< some constant to build hat functions */
const double C=1./2.;/**< some constant to build hat functions */
const double D=-2./15.;/**< some constant to build hat functions */
const double E=3./40.;/**< some constant to build hat functions */
const double u[NPI]   = {A,B,B,B,C};/**< some constants to build hat functions */
const double v[NPI]   = {A,B,B,C,B};/**< some constants to build hat functions */
const double w[NPI]   = {A,B,C,B,B};/**< some constants to build hat functions */
const double pds[NPI] = {D,E,E,E,E};/**< some constant weights to build hat functions */


/** \class prm
region number and material constants
*/
struct prm
	{
	int reg;/**< region number */	
	double alpha;/**< \f$ \alpha \f$ damping parameter */
	double A;/**< constant */
	double J;/**< exchange */
	double K;/**< uniaxial anisotropy constant */	
	double K3;/**< third order uniaxial anisotropy constant */	
	double uk[DIM][DIM]; /**< anisotropy tensor 3*3 */	
	double Uz;/**< for spin polarized current */
	double beta;/**< non adiabatic constant \f$ \beta \f$ for spin polarization current */	
	};

/**
initialisation of a constant matrix
*/
void init_dadu(gmm::dense_matrix <double> &X); // ça pourrait faire partie d'un constructeur si Tetra devient une classe

/** \class Tet
Tet is a tetrahedron, containing the index references to nodes, must not be flat <br>
indices convention is<br>
```
                        v
                      .
                    ,/
                   /
                2(ic)                                 2
              ,/|`\                                 ,/|`\
            ,/  |  `\                             ,/  |  `\
          ,/    '.   `\                         ,6    '.   `5
        ,/       |     `\                     ,/       8     `\
      ,/         |       `\                 ,/         |       `\
     0(ia)-------'.--------1(ib) --> u     0--------4--'.--------1
      `\.         |      ,/                 `\.         |      ,/
         `\.      |    ,/                      `\.      |    ,9
            `\.   '. ,/                           `7.   '. ,/
               `\. |/                                `\. |/
                  `3(id)                                `3
                     `\.
                        ` w
```
*/
class Tet{
	public:
		inline Tet() {reg = 0;} /**< default constructor */
		int reg;/**< .msh region number */
		int idxPrm;/**< index of the material parameters of the tetrahedron */		
		double vol;/**< volume of the tetrahedron */
		int ind[N];/**< indices to the nodes */
		double weight[NPI];/**< weights */
		double a[N][NPI];/**< hat functions */
		double dadx[N][NPI];/**< variations of hat function along x directions */
		double dady[N][NPI];/**< variations of hat function along y directions */
		double dadz[N][NPI];/**< variations of hat function along z directions */
	

		/**
		computes matrix u = unod*a; with unod = {node[ind].u0}_N and a from tetrahedron directly from vector<node>, without any copy		
		*/		
		void calc_u(double u[DIM][NPI],std::vector <Node> const& myNode);
		
		/**
		computes dudx = unod * dadx with unod = {node[ind].u0}_N		
		*/
		void calc_dudx(double dudx[DIM][NPI],std::vector <Node> const& myNode);

		/**
		computes dudy = unod * dady with unod = {node[ind].u0}_N		
		*/
		void calc_dudy(double dudy[DIM][NPI],std::vector <Node> const& myNode);

		/**
		computes dudz = unod * dadz with unod = {node[ind].u0}_N		
		*/
		void calc_dudz(double dudz[DIM][NPI],std::vector <Node> const& myNode);

		/**
		computes matrix v = vnod*a; with vnod = {node[ind].v0}_N and a from tetrahedron directly from vector<node>, without any copy		
		*/
		void calc_v(double v[DIM][NPI],std::vector <Node> const& myNode);


		/**
		computes dvdx = vnod * dadx with vnod = {node[ind].v0}_N		
		*/
		void calc_dvdx(double dvdx[DIM][NPI],std::vector <Node> const& myNode);

		/**
		computes dvdy = vnod * dady with vnod = {node[ind].v0}_N		
		*/
		void calc_dvdy(double dvdy[DIM][NPI],std::vector <Node> const& myNode);

		/**
		computes dvdz = vnod * dadz with vnod = {node[ind].v0}_N		
		*/
		void calc_dvdz(double dvdz[DIM][NPI],std::vector <Node> const& myNode);

		/**
		computes the integral contribution of the tetrahedron to the evolution of the magnetization
		*/		
		void integrales(std::vector<Tetra::prm> const& params,std::vector<Node> const& myNode,double Hext[DIM],double Vz,double theta,double dt,gmm::dense_matrix <double> &AE, std::vector <double> &BE);

		/**
		convenient getter for N, usefull for templates projection and assemblage
		*/
		inline int getN(void) {return N;}
		
		/**
		computes volume		
		*/
		void calc_vol(std::vector<Node> const& myNode);	
	};
}

#endif /* tetra_h */
