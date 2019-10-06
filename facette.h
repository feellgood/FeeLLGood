#ifndef facette_h
#define facette_h

/** \file facette.h
  \brief contains namespace Facette
  header containing Fac class, and some constants and a less_than operator to redo orientation of triangular faces
 */

#include <functional>

#include "config.h"
#include "node.h"

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
constexpr double a[N][NPI] = {{1.-u[0]-v[0],1.-u[1]-v[1],1.-u[2]-v[2],1.-u[3]-v[3]},{u[0],u[1],u[2],u[3]},{v[0],v[1],v[2],v[3]}};

/** \class prm
region number and material constants
*/
struct prm
	{
	int reg;/**< region number */	
	double Js;/**< surface exchange */
	double Ks;/**< uniaxial surface anisotropy constant */	
	double uk[DIM]; /**< anisotropy axis */	
	
	/** print the struct parameters */
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
        inline Obj(const int _idx)
        {
        idx = _idx;   
        }
        
        int idx;/**< index of the corresponding facette */
     };
    
    
/** \class Fac
Face is a class containing the index references to nodes, it has a triangular shape and should not be degenerated 
*/
class Fac{
	public:
		inline Fac(int _NOD /**< [in] */):NOD(_NOD),Ksp(2*N,2*N), Lsp(2*N)  /**< constructor */
            {reg = 0; treated = false;}
        
        /** constructor used by readMesh */
        inline Fac(const std::vector<Nodes::Node>  *_p_node /**< [in] pointer to the nodes */,
                   const int _reg /**< [in] region number */,
                   const int _idx /**< [in] region index in region vector */,
                   const int i0 /**< [in] node index */,
                   const int i1 /**< [in] node index */,
                   const int i2 /**< [in] node index */) : idxPrm(_idx),reg(_reg),refNode(_p_node),Ksp(2*N,2*N), Lsp(2*N)
            {
                NOD = refNode->size();
                if((0<i0)&&(i0<=NOD)&&(0<i1)&&(i1<=NOD)&&(0<i2)&&(i2<=NOD))
                    {
                    ind[0] = i0; ind[1] = i1; ind[2] = i2;
                    for (int i=0; i<N; i++) ind[i]--; // to force index to start from 0 (C++) instead of Matlab/msh convention
                    treated = false;
                    calc_surf();
                    init();
                    }
                else {std::cout<< "wrong indices in triangular facette." <<std::endl;SYSTEM_ERROR;}
            }
        
		/** constructor from a region number and three indices */		
		inline Fac(int r,int i0,int i1,int i2) {reg = r; ind[0]=i0;ind[1]=i1;ind[2]=i2;}
		
		/** constructor from a region number, idxPrm and three indices */		
		inline Fac(int r,int idx,int i0,int i1,int i2) {reg = r; idxPrm=idx; ind[0]=i0;ind[1]=i1;ind[2]=i2;}
		
		int idxPrm;/**< index of the material parameters of the facette */	
			
		double surf; /**< surface of the face */
		double Ms; /**< magnetization at saturation of the face */    
		Pt::pt3D n;/**< normal vector to the face */	
		int ind[N];/**< indices table of the nodes */
		double weight[NPI];/**< weights table */
		
		bool treated;/**< flag */
		
        
    
        /** weighted scalar product : factorized formulation */
        inline double weightedScalarProd(const double X[NPI] /**< [in] */) const
            {return ( X[0]*weight[0] + (X[1] +X[2] + X[3])*weight[1] );}
            //{return (X[0]*weight[0] + X[1]*weight[1] + X[2]*weight[2] + X[3]*weight[3] );}
        
        /** interpolation for 3D vector field : the getter function is given as a parameter in order to know what part of the node you want to interpolate */
        inline void interpolation(std::function<Pt::pt3D (Nodes::Node)> getter /**< [in] */,double result[DIM][NPI] /**< [out] */) const
        {
        double vec_nod[DIM][N];
        for (int i=0; i<N; i++)
            {
            Nodes::Node const& node = (*refNode)[ ind[i] ];
    
            vec_nod[0][i]   = getter(node).x();
            vec_nod[1][i]   = getter(node).y();
            vec_nod[2][i]   = getter(node).z();
            }
        double sx = (vec_nod[Pt::IDX_X][0] + vec_nod[Pt::IDX_X][1] + vec_nod[Pt::IDX_X][2] )/5.0;
        double sy = (vec_nod[Pt::IDX_Y][0] + vec_nod[Pt::IDX_Y][1] + vec_nod[Pt::IDX_Y][2] )/5.0;
        double sz = (vec_nod[Pt::IDX_Z][0] + vec_nod[Pt::IDX_Z][1] + vec_nod[Pt::IDX_Z][2] )/5.0;
        
        result[0][0] = result[0][1] = result[0][2] = result[0][3] = sx;
        result[1][0] = result[1][1] = result[1][2] = result[1][3] = sy;
        result[2][0] = result[2][1] = result[2][2] = result[2][3] = sz;
        
        result[0][0] *= (5.0/3.0); result[1][0] *= (5.0/3.0); result[2][0] *= (5.0/3.0);
        result[0][1] += vec_nod[0][0]*2.0/5.0; result[0][2] += vec_nod[0][1]*2.0/5.0; result[0][3] += vec_nod[0][2]*2.0/5.0;
        result[1][1] += vec_nod[1][0]*2.0/5.0; result[1][2] += vec_nod[1][1]*2.0/5.0; result[1][3] += vec_nod[1][2]*2.0/5.0;
        result[2][1] += vec_nod[2][0]*2.0/5.0; result[2][2] += vec_nod[2][1]*2.0/5.0; result[2][3] += vec_nod[2][2]*2.0/5.0;
        
        //tiny::mult<double, DIM, N, NPI> (vec_nod, a, result);
        }
        
        /** interpolation for scalar field : the getter function is given as a parameter in order to know what part of the node you want to interpolate */
        inline void interpolation(std::function<double (Nodes::Node)> getter /**< [in] */,double result[NPI] /**< [out] */) const
        {
        double scalar_nod[N];    
        for (int i=0; i<N; i++)
            {
            Nodes::Node const& n = (*refNode)[ ind[i] ];
            scalar_nod[i] =  getter(n);
            }
        
        double sn = scalar_nod[0] + scalar_nod[1] + scalar_nod[2];
        
        result[0] = sn/3.0; result[1] = (sn + 2.0*scalar_nod[0])/5.0; result[2] = (sn + 2.0*scalar_nod[1])/5.0; result[2] = (sn + 2.0*scalar_nod[2])/5.0;
        //tiny::transposed_mult<double, N, NPI> (scalar_nod, a, result);
        }
        
		/** computes the integral contribution of the triangular face */
		void integrales(std::vector<Facette::prm> const& params /**< [in] */, std::vector <double> &BE /**< [out] */) const;
		
        /** anisotropy energy of the facette */
        double anisotropyEnergy(Facette::prm const& param /**< [in] */,const double u[DIM][NPI] /**< [in] */) const;
        
        /** demagnetizing energy of the facette */
        double demagEnergy(const double u[DIM][NPI] /**< [in] */,const double phi[NPI] /**< [in] */) const;
        
        /** compute projection of a face from inner object matrix */
        void projection(gmm::dense_matrix <double> const& A, std::vector <double> const& B);
        
        /** assemblage of the matrix and vector elements from inner matrix in facette object */
        void assemblage(write_matrix &K,write_vector &L) const;
        
        /** getter for N */		
		inline int getN(void) const {return N;}	

		/** getter for NPI */		
		inline int getNPI(void) const {return NPI;}	
		
        /**
        computes correction on potential
        */
        double potential(std::function<Pt::pt3D (Nodes::Node)> getter, int i) const;
        
        void calcCorr(std::function<const Pt::pt3D (Nodes::Node)> getter,double *corr,double u[DIM][NPI]) const;
        
    private:
        int NOD;/**< number of nodes */
        int reg;/**< .msh region number */
        
        
        const std::vector<Nodes::Node>  *refNode;/**< direct access to the Nodes */
        gmm::dense_matrix <double> Ksp;/**< matrix initialized by constructor */
        std::vector <double> Lsp;/**< vector initialized by constructor */
        
        /** computes normal vector n and surface surf		*/		
		void calc_surf(void);
        
        /** initialize weight hat function */
        inline void init(void)
            {for (int j=0; j<NPI; j++) {weight[j] = 2.*surf*pds[j]; }}// detJ = 2*surf;
};//end class Fac

/*    
inline bool operator< (const Fac &f1, const Fac &f2)
    {
    if (f1.ind[0]<f2.ind[0]) return true;
    else
        if ((f1.ind[0]==f2.ind[0]) && (f1.ind[1]<f2.ind[1])) return true;
        else
            if ((f1.ind[0]==f2.ind[0]) && (f1.ind[1]==f2.ind[1]) && (f1.ind[2]<f2.ind[2])) return true;

    return false;
    }
*/

/** operator less_than for the orientation of the facette, lexicographic order */
struct less_than
{
/** operator() for the comparison of two faces with lexicographical order */

bool operator()(const Fac &f1, const Fac &f2) const
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
