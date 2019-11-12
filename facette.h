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
	Pt::pt3D uk; /**< anisotropy axis */	
	
	/** print the struct parameters */
	inline void infos()
		{
		std::cout<< "surface region number = " << reg <<std::endl;
		std::cout<< "Js = " << Js <<std::endl;
		
		if(Ks!=0)
			{ std::cout<< "Ks*a = "<<Ks << "*[ " << uk << "]" << std::endl; }
		else std::cout << "no surface anisotropy" << std::endl;
		};	
		
	};

/** \class Fac
Face is a class containing the index references to nodes, it has a triangular shape and should not be degenerated 
*/
class Fac{
	public:
		inline Fac(int _NOD /**< [in] */):NOD(_NOD),reg(0)  /**< constructor */
            { treated = false; }
        
        /** constructor used by readMesh */
        inline Fac(const std::vector<Nodes::Node>  *_p_node /**< [in] pointer to the nodes */,
                   const int _reg /**< [in] region number */,
                   const int _idx /**< [in] region index in region vector */,
                   const int i0 /**< [in] node index */,
                   const int i1 /**< [in] node index */,
                   const int i2 /**< [in] node index */) : idxPrm(_idx),NOD(_p_node->size()),reg(_reg),refNode(_p_node)
            {
		ind[0] = i0; ind[1] = i1; ind[2] = i2;
                for (int i=0; i<N; i++) ind[i]--; // to force index to start from 0 (C++) instead of Matlab/msh convention
                treated = false;
		double surf = calc_surf();
                for (int j=0; j<NPI; j++) {weight[j] = 2.*surf*pds[j]; }
            }
        
		/** constructor from a region number and three indices */		
		inline Fac(int r,int i0,int i1,int i2): NOD(0),reg(r)
		{ind[0]=i0;ind[1]=i1;ind[2]=i2;}
		
		/** constructor from a region number, idxPrm and three indices */		
		inline Fac(int r,int idx,int i0,int i1,int i2): NOD(0),reg(r) 
		{idxPrm=idx; ind[0]=i0;ind[1]=i1;ind[2]=i2;}
		
		int idxPrm;/**< index of the material parameters of the facette */	
			
		double Ms; /**< magnetization at saturation of the face */    
		int ind[N];/**< indices table of the nodes */
		double weight[NPI];/**< weights table */
		
		bool treated;/**< flag */
		
        /** weighted scalar product : factorized formulation */
        inline double weightedScalarProd(const double X[NPI] /**< [in] */) const
            {return ( X[0]*weight[0] + (X[1] +X[2] + X[3])*weight[1] );}
            //{return (X[0]*weight[0] + X[1]*weight[1] + X[2]*weight[2] + X[3]*weight[3] );}
        
        /** interpolation for 3D vector field : the getter function is given as a parameter in order to know what part of the node you want to interpolate */
        inline void interpolation(std::function<Pt::pt3D (Nodes::Node)> getter /**< [in] */,double result[Pt::DIM][NPI] /**< [out] */) const
        {
        double vec_nod[Pt::DIM][N];
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
        
        /** basic infos */		
	inline void infos() const {std::cout<< "reg="<< reg << ":" << idxPrm << "ind:"<< ind[0]<< "\t"<< ind[1]<< "\t"<< ind[2] <<std::endl;};
        
	/** computes the integral contribution of the triangular face */
	void integrales(std::vector<Facette::prm> const& params /**< [in] */, Pt::pt3D BE[N] /**< [out] */) const;
		
        /** anisotropy energy of the facette */
        double anisotropyEnergy(Facette::prm const& param /**< [in] */,const double u[Pt::DIM][NPI] /**< [in] */) const;
        
        /** surface charges  */
        void charges(std::function<Pt::pt3D (Nodes::Node)> getter,double *srcDen,double *corr,int &nsrc) const;
        
        /** demagnetizing energy of the facette */
        double demagEnergy(const double u[Pt::DIM][NPI] /**< [in] */,const double phi[NPI] /**< [in] */) const;
        
        /** assemblage of the matrix elements from inner matrix in facette object */
        void assemblage_mat(write_matrix &K) const;
        
        /** assemblage of the vector elements from inner matrix in facette object */
        inline void assemblage_vect(std::vector<double> &L) const
            { for (int i=0; i < N; i++) { L[NOD+ind[i]] += Lp[i]; L[ind[i]] += Lp[N+i]; } }
        
        /** getter for N */		
	inline int getN(void) const {return N;}	

	/** getter for NPI */		
	inline int getNPI(void) const {return NPI;}	
		
	/** computes correction on potential*/
        double potential(std::function<Pt::pt3D (Nodes::Node)> getter, int i) const;
        
        /** computes correction values */
        void calcCorr(std::function<const Pt::pt3D (Nodes::Node)> getter,double *corr,double u[Pt::DIM][NPI]) const;
        
        /** lexicographic order on indices */
        inline bool operator< (const Fac &f) const
            {
            if (this->ind[0] < f.ind[0]) return true;
            else if ((this->ind[0] == f.ind[0]) && (this->ind[1] < f.ind[1])) return true;
                else if ((this->ind[0] == f.ind[0]) && (this->ind[1] == f.ind[1]) && (this->ind[2] < f.ind[2])) return true;

            return false;
            }
    
	/** small matrix for integrales */
	double Kp[2*N][2*N];        
	
    	/** small vector for integrales */
    	double Lp[2*N];

	/** computes the norm to the face */
	Pt::pt3D calc_norm(void) const;

	/** computes surface of the face		*/		
	double calc_surf(void) const;
    
    private:
        const int NOD;/**< number of nodes */
        const int reg;/**< .msh region number */
        
        const std::vector<Nodes::Node>  *refNode;/**< direct access to the Nodes */
        
};//end class Fac


/** operator less_than for the orientation of the facette, lexicographic order */
struct less_than 
{
    /** operator () for ordering facettes with operator< overloaded */
    bool operator()(const Fac &f1, const Fac &f2) const {return (f1<f2);} 
};

}//end namespace

#endif /* facette_h */
