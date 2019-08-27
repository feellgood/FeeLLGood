#ifndef facette_h
#define facette_h

/** \file facette.h
  \brief contains namespace Facette
  header containing Fac class, and some constants and a less_than operator to redo orientation of triangular faces
 */

#include "gmm/gmm_kernel.h"

typedef gmm::wsvector <double>   write_vector;/**< gmm write vector */
typedef gmm::rsvector <double>   read_vector; /**< gmm read vector */

typedef gmm::row_matrix	<write_vector>   write_matrix; /**< gmm write sparse matrix */
typedef gmm::row_matrix	<read_vector>    read_matrix; /**< gmm read sparse matrix */

#include "config.h"
#include "node.h"

/** convenient typedef for mtl4 */
//typedef typename mtl::Collection< mtl::compressed2D<double> >::value_type v_type;

/** convenient typedef for mtl4 inserter */
//typedef mtl::mat::inserter< mtl::compressed2D<double>,mtl::update_plus<v_type> > sparseInserter;

/** \namespace Facette
 to grab altogether some constants and calculation functions for class Fac
 */
namespace Facette
{
const int N = 3; /**< number of sommits */
const int NPI = 4; /**< number of weights  */

const double u[NPI]   = {   1/3.,   1/5.,   3/5.,   1/5.};/**< some constants to build hat functions */
const double v[NPI]   = {   1/3.,   1/5.,   1/5.,   3/5.};/**< some constants  to build hat functions */
const double pds[NPI] = {-27/96., 25/96., 25/96., 25/96.};/**< some constant weights  to build hat functions */

/** hat function constants */
const double a[N][NPI] ={
    {1.-u[0]-v[0],1.-u[1]-v[1],1.-u[2]-v[2],1.-u[3]-v[3]},{u[0],u[1],u[2],u[3]},{v[0],v[1],v[2],v[3]}};

// for (int j=0; j<NPI; j++) { a[0][j] = 1.-u[j]-v[j]; a[1][j] = u[j]; a[2][j] = v[j];}


/** \class prm
region number and material constants
*/
struct prm
	{
	int reg;/**< region number */	
	double Js;/**< surface exchange */
	double Ks;/**< uniaxial surface anisotropy constant */	
	double uk[DIM]; /**< anisotropy axis */	
	
	/**
	 print the struct parameters
	 */
	inline void infos()
		{
		std::cout<< "surface region number = " << reg <<std::endl;
		std::cout<< "Js = " << Js <<std::endl;
		
		if(Ks!=0)
			{std::cout<< "Ks = " << Ks <<std::endl;
			std::cout<< "a = [ " << uk[0] << "," << uk[1] <<"," << uk[2] << "]" << std::endl;
			}
		else std::cout << "no surface anisotropy" << std::endl;
		};	
		
	};

    /** \class Obj
     container for buffering the projection(facettes) matrix results when mutex is locked
     */
    
class Obj
    {
    public:
        /** constructor */
        inline Obj(const int _idx)//,gmm::dense_matrix <double> const& K,std::vector <double> const& L)
        {
        idx = _idx;   
        //Ke=K;Le=L;
        }
        
        int idx;/**< index of the corresponding facette */
        //gmm::dense_matrix <double> Ke;/**< small matrix resulting from projection */
        //std::vector <double> Le;/**< small vector resulting from projection */
    };
    
    
/** \class Fac
Face is a class containing the index references to nodes, it has a triangular shape and should not be degenerated 
*/
class Fac{
	public:
		inline Fac(int _NOD):NOD(_NOD),Ksp(2*N,2*N), Lsp(2*N)  /**< constructor */
            {reg = 0; treated = false;}
        
		/** constructor from a region number and three indices */		
		inline Fac(int r,int i0,int i1,int i2) {reg = r; ind[0]=i0;ind[1]=i1;ind[2]=i2;}
		
		/** constructor from a region number, idxPrm and three indices */		
		inline Fac(int r,int idx,int i0,int i1,int i2) {reg = r; idxPrm=idx; ind[0]=i0;ind[1]=i1;ind[2]=i2;}
		
		int reg;/**< .msh region number */
		int idxPrm;/**< index of the material parameters of the facette */		
		double surf; /**< surface of the face */
		double Ms; /**< magnetization at saturation of the face */    
		Pt::pt3D n;/**< normal vector to the face */	
		int ind[N];/**< indices table of the nodes */
		double weight[NPI];/**< weights table */
		
		bool treated;/**< flag */
		
        /** initialize weight  */
        inline void init(void)
            {for (int j=0; j<NPI; j++) {weight[j] = 2.*surf*pds[j]; }}// detJ = 2*surf;
    
        /** weighted scalar product */
        inline double weightedScalarProd(const double X[NPI]) const
            {return (X[0]*weight[0] + X[1]*weight[1] + X[2]*weight[2] + X[3]*weight[3] );}
        
        /** interpolation */
        void interpolation(double u[DIM][NPI]) const;
        
		/** computes the integral contribution of the triangular face */
		void integrales(std::vector<Facette::prm> const& params, std::vector <double> &BE) const;
		
        /** total energy of the facette = anisotropy + demag */
        void energy(Facette::prm const& param,double E[5]);
        
        /** anisotropy energy of the facette */
        double anisotropyEnergy(Facette::prm const& param) const;
        
        /** demagnetizing energy of the facette */
        double demagEnergy(void) const;
        
        /** compute projection of a face */
        void projection(gmm::dense_matrix <double> const& A, std::vector <double> const& B,gmm::dense_matrix <double> &Ap, std::vector <double> &Bp) const;
        
        /** compute projection of a face from inner object matrix */
        void projection(gmm::dense_matrix <double> const& A, std::vector <double> const& B);
        
        /** assemblage of the matrix and vector elements */
        void assemblage(const int NOD,gmm::dense_matrix <double> const& Ke, std::vector <double> const& Le,write_matrix &K,write_vector &L) const;
		
        /** assemblage of the matrix and vector elements from inner matrix in facette object */
        void assemblage(write_matrix &K,write_vector &L) const;
        
        /**
		convenient getter for N
		*/		
		inline int getN(void) {return N;}	

		/**
		computes normal vector n and surface surf		
		*/		
		void calc_surf(void);
        
        /** pointer to the nodes */
        inline void setRefNode(std::vector<Nodes::Node>  *_p_node) {refNode = _p_node;}
        
    private:
        int NOD;/**< number of nodes */
        std::vector<Nodes::Node>  *refNode;/**< direct access to the Nodes */
        gmm::dense_matrix <double> Ksp;/**< matrix initialized by constructor */
        std::vector <double> Lsp;/**< vector initialized by constructor */
	};

/**
operator less_than for the orientation of the facette, lexicographic order
*/
struct less_than
{
/**
operator() for the comparison of two faces with lexicographical order
*/
bool operator()(Fac f1, Fac f2) const
  {
  if (f1.ind[0]<f2.ind[0]) return true;
  else
     if ((f1.ind[0]==f2.ind[0]) && (f1.ind[1]<f2.ind[1])) return true;
     else
        if ((f1.ind[0]==f2.ind[0]) && (f1.ind[1]==f2.ind[1]) && (f1.ind[2]<f2.ind[2])) return true;

  return false;
  }
};

}

#endif /* facette_h */
