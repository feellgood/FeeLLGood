#ifndef tetra_h
#define tetra_h

/** \file tetra.h
  \brief namespace Tetra
  header containing Tet class, some constants, and integrales
 */

#include "gmm_kernel.h" // pour dense_matrix dans namespace Tetra

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

/** constant matrix */
const double dadu[N][DIM] = {{-1.,-1.,-1.},{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};


/** constant matrix \f$ j:0..4 \f$
a[0][j]   = 1.-u[j]-v[j]-w[j];
a[1][j]   = u[j];
a[2][j]   = v[j];
a[3][j]   = w[j];
*/
const double a[N][NPI] = {{1.-u[0]-v[0]-w[0],1.-u[1]-v[1]-w[1],1.-u[2]-v[2]-w[2],1.-u[3]-v[3]-w[3],1.-u[4]-v[4]-w[4]},
{u[0],u[1],u[2],u[3],u[4]}, {v[0],v[1],v[2],v[3],v[4]}, {w[0],w[1],w[2],w[3],w[4]}};


/** \class prm
region number and material constants
*/
struct prm
	{
	int reg;/**< region number */	
	double alpha;/**< \f$ \alpha \f$ damping parameter */
	double A;/**< exchange constant stiffness */
	double J;/**< \f$ M_s = \nu_0 J \f$ */
	double K;/**< uniaxial anisotropy constant */	
	double K3;/**< third order uniaxial anisotropy constant */	
	double uk[DIM][DIM]; /**< anisotropy tensor 3*3 */	
	double Uz;/**< for spin polarized current */
	double beta;/**< non adiabatic constant \f$ \beta \f$ for spin polarization current */	
	
	/**
	 print the struct parameters
	 */
	inline void infos()
		{
		std::cout<< "volume region number = " << reg <<std::endl;
		std::cout<< "alpha = " << alpha <<std::endl;
		std::cout<< "A = " << A <<std::endl;
		std::cout<< "J = " << J <<std::endl;
		
		if(K!=0)
			{std::cout<< "K = " << K <<std::endl;
			std::cout<< "a = [ [ " << uk[0][0] << "," << uk[0][1] <<"," << uk[0][2] << "]" << std::endl;
			std::cout<< "      [ " << uk[1][0] << "," << uk[1][1] <<"," << uk[1][2] << "]" << std::endl;
			std::cout<< "      [ " << uk[2][0] << "," << uk[2][1] <<"," << uk[2][2] << "] ]" << std::endl;
			}
		else std::cout << "no anisotropy" << std::endl;
		};	
	};

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
		inline Tet() {reg = 0;idxPrm=-1;} /**< default constructor */
		int reg;/**< .msh region number */
		int idxPrm;/**< index of the material parameters of the tetrahedron */		
		double vol;/**< volume of the tetrahedron */
		int ind[N];/**< indices to the nodes */
		double weight[NPI];/**< weights \f$ w_i = |J| p_i  \f$ with  \f$ p_i = pds[i] = (D,E,E,E,E) \f$ */
		double dadx[N][NPI];/**< variations of hat function along x directions */
		double dady[N][NPI];/**< variations of hat function along y directions */
		double dadz[N][NPI];/**< variations of hat function along z directions */
	
		/** initializes weight and dad(x|y|z) */
		void init(std::vector<Node> const& myNode,double epsilon);
		
		/** basic region infos */		
		inline void infos(){std::cout<< reg << ":" << idxPrm <<std::endl;};
		
		/**
		computes the integral contribution of the tetrahedron to the evolution of the magnetization
		*/		
		void integrales(std::vector<Tetra::prm> const& params,std::vector<Node> const& myNode,double Hext[DIM],double Vz,double theta,double dt,double tau_r,gmm::dense_matrix <double> &AE, std::vector <double> &BE);

		/**
		convenient getter for N, usefull for templates projection and assemblage
		*/
		inline int getN(void) {return N;}
		
		/**
        initializes nod matrix from vector myNode
        */
		void getNod(gmm::dense_matrix <double> &nod,std::vector <Node> const& myNode);
		
        /**
        \return \f$ |J| \f$ build Jacobian \f$ J \f$
        */
        double Jacobian(double J[DIM][DIM],std::vector <Node> const& myNode);
        
		/**
		computes volume		
		*/
		void calc_vol(std::vector<Node> const& myNode);
    };//end class Tetra
}

#endif /* tetra_h */
