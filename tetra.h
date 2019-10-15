#ifndef tetra_h
#define tetra_h

/** \file tetra.h
  \brief namespace Tetra
  header containing Tet class, some constants, and integrales
 */
#include <functional>

#include "gmm/gmm_kernel.h"

#include "config.h"
#include "node.h"

#include "tiny.h"

/** a 3D point */
typedef double triple[DIM];

/** 
\return square of a number \f$ x^2 \f$
*/
inline double sq(double x /**< [in] */ ) {return x*x;}

/** in place normalizing function of triple a */
inline void normalize(triple &a /**< [in,out] */)
{
double norme=sqrt(sq(a[0])+sq(a[1])+sq(a[2]));
if(norme>0)
    { a[0]/= norme;a[1]/= norme;a[2]/= norme; }
else {std::cerr<< "division by zero : normalize." <<std::endl;}
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
constexpr double a[N][NPI] = {{1.-u[0]-v[0]-w[0],1.-u[1]-v[1]-w[1],1.-u[2]-v[2]-w[2],1.-u[3]-v[3]-w[3],1.-u[4]-v[4]-w[4]},
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
	
	/** print the struct parameters */
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

    /** \class Obj
     container for buffering the projection(Tetrahedron) matrix results when mutex is locked
     */
    class Obj
    {
    public:
        /** constructor */
        inline Obj(const int _idx)
        {
        idx = _idx;  
        }
    int idx;/**< index of the corresponding tetrahedron */
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
		inline Tet(int _NOD): NOD(_NOD), Kp(2*N,2*N), Lp(2*N) 
        { reg = 0; idxPrm=-1; treated = false;} /**< default constructor */
		
		/** constructor for readMesh */
		inline Tet(const std::vector<Nodes::Node>  *_p_node /**< [in] pointer to the nodes */,
                   const int _reg /**< [in] region number */,
                   const int _idx /**< [in] region index in region vector */,
                   const int i0 /**< [in] node index */,
                   const int i1 /**< [in] node index */,
                   const int i2 /**< [in] node index */,
                   const int i3 /**< [in] node index */,
                    const double epsilon/**< for degeneracy test */) : idxPrm(_idx),reg(_reg),refNode(_p_node), Kp(2*N,2*N), Lp(2*N) 
            {
            NOD = _p_node->size();
            ind[0] = i0; ind[1] = i1; ind[2] = i2; ind[3] = i3;
            for (int i=0; i<N; i++) ind[i]--;           // convention Matlab/msh -> C++
            treated = false;
            calc_vol();
            init(epsilon);
            } 
		
		
		int idxPrm;/**< index of the material parameters of the tetrahedron */		
		double vol;/**< volume of the tetrahedron */
		int ind[N];/**< indices to the nodes */
		double weight[NPI];/**< weights \f$ w_i = |J| p_i  \f$ with  \f$ p_i = pds[i] = (D,E,E,E,E) \f$ */
		double dadx[N][NPI];/**< variations of hat function along x directions */
		double dady[N][NPI];/**< variations of hat function along y directions */
		double dadz[N][NPI];/**< variations of hat function along z directions */
        
        bool treated;/**< flag */
        
		/** initializes weight hat function and dad(x|y|z) if \f$ | detJ | < \epsilon \f$ jacobian is considered degenerated */
		void init(double epsilon);
        
        
        /** weighted scalar product */
        inline double weightedScalarProd(const double X[NPI]) const
            {return (X[0]*weight[0] + X[1]*weight[1] + X[2]*weight[2] + X[3]*weight[3] +X[4]*weight[4]);}
		
		
		
		/** interpolation for 3D vector field: the getter function is given as a parameter in order to know what part of the node you want to interpolate */
		inline void interpolation(std::function<Pt::pt3D (Nodes::Node)> getter,double result[DIM][NPI]) const
        {
        double vec_nod[DIM][N];
        getVecDataFromNode(getter,vec_nod);
        
        tiny::mult<double, DIM, N, NPI> (vec_nod, a, result);
        }
		
		/** interpolation for a tensor : the getter function is given as a parameter in order to know what part of the node you want to interpolate */
        inline void interpolation(std::function<Pt::pt3D (Nodes::Node)> getter,double Tx[DIM][NPI],double Ty[DIM][NPI],double Tz[DIM][NPI]) const
        {
        double vec_nod[DIM][N];
        getVecDataFromNode(getter,vec_nod);
        
        tiny::mult<double, DIM, N, NPI> (vec_nod, dadx, Tx);
        tiny::mult<double, DIM, N, NPI> (vec_nod, dady, Ty);
        tiny::mult<double, DIM, N, NPI> (vec_nod, dadz, Tz);
        }
		
		
		/** interpolation for 3D vector field and a tensor : the getter function is given as a parameter in order to know what part of the node you want to interpolate */
        inline void interpolation(std::function<Pt::pt3D (Nodes::Node)> getter,double result[DIM][NPI],
                                  double Tx[DIM][NPI],double Ty[DIM][NPI],double Tz[DIM][NPI]) const
        {
        double vec_nod[DIM][N];
        getVecDataFromNode(getter,vec_nod);
        
        tiny::mult<double, DIM, N, NPI> (vec_nod, a, result);
        tiny::mult<double, DIM, N, NPI> (vec_nod, dadx, Tx);
        tiny::mult<double, DIM, N, NPI> (vec_nod, dady, Ty);
        tiny::mult<double, DIM, N, NPI> (vec_nod, dadz, Tz);
        }
		
		/** interpolation for components of a field : the getter function is given as a parameter in order to know what part of the node you want to interpolate */
        inline void interpolation(std::function<double (Nodes::Node)> getter,double Xx[NPI],double Xy[NPI],double Xz[NPI]) const
        {
        double scalar_nod[N];    
        getScalDataFromNode(getter,scalar_nod);
        tiny::neg_transposed_mult<double, N, NPI> (scalar_nod, dadx, Xx);
        tiny::neg_transposed_mult<double, N, NPI> (scalar_nod, dady, Xy);
        tiny::neg_transposed_mult<double, N, NPI> (scalar_nod, dadz, Xz);
        }
		
		/** interpolation for scalar field : the getter function is given as a parameter in order to know what part of the node you want to interpolate */
        inline void interpolation(std::function<double (Nodes::Node)> getter,double result[NPI]) const
        {
        double scalar_nod[N];    
        getScalDataFromNode(getter,scalar_nod);
        tiny::transposed_mult<double, N, NPI> (scalar_nod, a, result);
        }
		
		/** interpolation for scalar field : the getter function is given as a parameter in order to know what part of the node you want to interpolate */
        inline void interpolation(std::function<double (Nodes::Node,Pt::index)> getter,Pt::index idx,double result[NPI]) const
        {
        double scalar_nod[N];    
        getScalDataFromNode(getter,idx,scalar_nod);
        tiny::transposed_mult<double, N, NPI> (scalar_nod, a, result);
        }
		
		/** basic infos */		
		inline void infos() const {std::cout<< "reg="<< reg << ":" << idxPrm << "ind:"<< ind[0]<< "\t"<< ind[1]<< "\t"<< ind[2]<< "\t"<< ind[3] <<std::endl;};
		
        /** more infos */		
		inline void infos(Nodes::index idx) const 
            {
            std::cout<< "reg="<< reg << ":" << idxPrm << "ind:"<< ind[0]<< "\t"<< ind[1]<< "\t"<< ind[2]<< "\t"<< ind[3] <<std::endl;
            switch(idx)
                {
                case Nodes::IDX_p : for(int i=0;i<N;i++) {std::cout<<"p_" << i<<  ((*refNode)[ ind[i] ]).p  <<std::endl;} break;
                case Nodes::IDX_u : for(int i=0;i<N;i++) {std::cout<<"m_" << i<< ((*refNode)[ ind[i] ]).u  <<std::endl;} break;
                case Nodes::IDX_phi : for(int i=0;i<N;i++) {std::cout<<"phi" << i<< ((*refNode)[ ind[i] ]).phi  <<std::endl;} break;
                default:break;
                }
            };
        
		/** computes the integral contribution of the tetrahedron to the evolution of the magnetization */		
		void integrales(std::vector<Tetra::prm> const& params,double Hext[DIM],double Vz,double theta,double dt,double tau_r,gmm::dense_matrix <double> &AE, std::vector <double> &BE)  const;

        /** computes all the contributions to the energy of the tetrahedron */
        //void energy(Tetra::prm const& param,double E[4],const double Hext[DIM],double uz_drift) const;
        
        /** exchange energy of the tetrahedron */
        double exchangeEnergy(Tetra::prm const& param,const double dudx[DIM][NPI],const double dudy[DIM][NPI],const double dudz[DIM][NPI]) const;
        
        /** anisotropy energy of the tetrahedron */
        double anisotropyEnergy(Tetra::prm const& param,const double u[DIM][NPI]) const;
        
        /** volume charges  */
        void charges(std::function<Pt::pt3D (Nodes::Node)> getter,double *srcDen,int &nsrc,double Ms) const;
        
        /** demagnetizing energy of the tetrahedron */
        double demagEnergy(Tetra::prm const& param,const double dudx[DIM][NPI],const double dudy[DIM][NPI],const double dudz[DIM][NPI],const double phi[NPI]) const;
        
        /** zeeman energy of the tetrahedron */
        double zeemanEnergy(Tetra::prm const& param,double uz_drift,const double Hext[DIM],const double u[DIM][NPI]) const;
        
        /** computes projection of a tetrahedron using inner matrix in tetra object */
        void projection(gmm::dense_matrix <double> const& A,  std::vector <double> const& B);
        
        /** matrix assembly using inner matrix in tetra */
        void assemblage_mat(write_matrix &K) const;
        
        /** vector assembly using inner matrix in tetra */
        void assemblage_vect(write_vector &L) const;
        
        /** getter for N */
		inline int getN(void) const {return N;}
		
		/** getter for NPI */
		inline int getNPI(void) const {return NPI;}
		
		/** getter for region */
		inline int getRegion(void) const {return reg;}
		
        /** \return \f$ |J| \f$ build Jacobian \f$ J \f$ */
        double Jacobian(double J[DIM][DIM]);
        
    private:
        int NOD;/**< total number of nodes, also an offset for filling sparseMatrix */
        int reg;/**< .msh region number */
        
        const std::vector<Nodes::Node>  *refNode;/**< direct access to the Nodes */
        
        /** computes volume	of the tetrahedron */
		void calc_vol(void);
        
        /** getter to access and copy some vector parts of the node vector */
		inline void getVecDataFromNode(std::function<Pt::pt3D (Nodes::Node)> getter,double vecData[DIM][N]) const
		{
        for (int i=0; i<N; i++)
            {
            Nodes::Node const& node = (*refNode)[ ind[i] ];
    
            vecData[Pt::IDX_X][i]   = getter(node).x();
            vecData[Pt::IDX_Y][i]   = getter(node).y();
            vecData[Pt::IDX_Z][i]   = getter(node).z();
            }
        }
		
		/** getter to access and copy some scalar parts of the node vector */
		inline void getScalDataFromNode(std::function<double (Nodes::Node)> getter,double scalData[N]) const
		{
        for (int i=0; i<N; i++)
            { scalData[i] =  getter( (*refNode)[ ind[i] ] ); }    
        }
		
		/** getter to access and copy some scalar parts (vector components) of the node vector */
		inline void getScalDataFromNode(std::function<double (Nodes::Node,Pt::index)> getter,Pt::index idx,double scalData[N]) const
		{
        for (int i=0; i<N; i++)
            { scalData[i] =  getter( (*refNode)[ ind[i] ] ,idx); }    
        }
        
        gmm::dense_matrix <double> Kp;/**< Kp(2*N,2*N) initialized by constructor */
        std::vector <double>  Lp;/**< Lp(2*N) initialized by constructor  */
};//end class Tetra
}//end namespace Tetra

#endif /* tetra_h */
