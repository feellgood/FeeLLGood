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
	bool suppress_charges; /**< suppress charges if true */
	double Ks;/**< uniaxial surface anisotropy constant */	
	Pt::pt3D uk; /**< anisotropy axis */	
	
	/** print the struct parameters */
	inline void infos()
		{
		std::cout<< "surface region number = " << reg << "suppress charges = " << suppress_charges <<std::endl;
		
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
        /** constructor */
	inline Fac(int _NOD /**< [in] */):treated(false),NOD(_NOD),reg(0) { }

        /** constructor used by readMesh */
        inline Fac(const std::shared_ptr<Nodes::Node[]> _p_node /**< [in] pointer to the node */,
                   const int _NOD /**< [in] nb nodes */,
                   const int _reg /**< [in] region number */,
                   const int _idx /**< [in] region index in region vector */,
                   const int i0 /**< [in] node index */,
                   const int i1 /**< [in] node index */,
                   const int i2 /**< [in] node index */) : idxPrm(_idx),treated(false),NOD(_NOD),reg(_reg),refNode(_p_node)
        {
        ind[0] = i0; ind[1] = i1; ind[2] = i2;
        for (int i=0; i<N; i++) ind[i]--; // to force index to start from 0 (C++) instead of Matlab/msh convention
        //treated = false;
        if(NOD>0)
            { surf = calc_surf(); }
        else
            { surf = 0.0; Ms = 0.0; }
        }
		
		/** constructor from a region number, idxPrm and three indices */		
		inline Fac(const int r,const int idx,const int i0,const int i1,const int i2): idxPrm(idx),surf(0),Ms(0),treated(false),NOD(0),reg(r)
            {ind[0]=i0;ind[1]=i1;ind[2]=i2;}
		
		int idxPrm;/**< index of the material parameters of the facette */	
		
		double surf;/**< surface of the element */
		double Ms; /**< magnetization at saturation of the face */    
		int ind[N];/**< indices table of the nodes */
		bool treated;/**< flag */

        /** weighted scalar product : factorized formulation: weight(1)=weight(2)=weight(3) */
        inline double weightedScalarProd(const double X[NPI] /**< [in] */) const
            { return ( X[0]*weight(0) + (X[1] +X[2] + X[3])*weight(1) ); }
        
        /** interpolation template; T == 3D vector field or T == double .The getter function is given as a parameter in order to know what part of the node you want to interpolate */
        //To check with reference code : is there a missing transposition ?
        //tiny::mult<double, DIM, N, NPI> (vec_nod, a, result); //if T == double
        //tiny::transposed_mult<double, N, NPI> (scalar_nod, a, result); //if T == PT::pt3D
        template <class T>
        void interpolation(std::function< T (Nodes::Node)> getter /**< [in] */,T result[NPI] /**< [out] */) const
        {
        T X_nod[N];
        for (int i=0; i<N; i++)
            { X_nod[i] = getter( refNode[ ind[i] ] ); }
        
        T s = X_nod[0] + X_nod[1] + X_nod[2];
        result[0] = s/3.0;
        result[1] = (s + 2.0*X_nod[0])/5.0;
        result[2] = (s + 2.0*X_nod[1])/5.0;
        result[3] = (s + 2.0*X_nod[2])/5.0;
        }
        
        /** interpolation for 3D vector field : the getter function is given as a parameter in order to know what part of the node you want to interpolate */
        inline void interpolation(std::function<Pt::pt3D (Nodes::Node)> getter /**< [in] */,Pt::pt3D result[NPI] /**< [out] */) const
            { interpolation<Pt::pt3D>(getter,result); }
     
        /** interpolation for scalar field : the getter function is given as a parameter in order to know what part of the node you want to interpolate */
        inline void interpolation(std::function<double (Nodes::Node)> getter /**< [in] */,double result[NPI] /**< [out] */) const
            { interpolation<double>(getter,result); }

        /** basic infos */		
        inline void infos() const {std::cout<< "reg="<< reg << ":" << idxPrm << "ind:"<< ind[0]<< "\t"<< ind[1]<< "\t"<< ind[2] <<std::endl;};
        
        /** computes the integral contribution of the triangular face */
        void integrales(std::vector<Facette::prm> const& params /**< [in] */, Pt::pt3D BE[N] /**< [out] */) const;

        /** anisotropy energy of the facette */
        double anisotropyEnergy(Facette::prm const& param /**< [in] */,const Pt::pt3D u[NPI] /**< [in] */) const;
        
        /** surface charges  */
        void charges(std::function<Pt::pt3D (Nodes::Node)> getter,std::vector<double> &srcDen,std::vector<double> &corr,int &nsrc) const;
        
        /** demagnetizing energy of the facette */
        double demagEnergy(const Pt::pt3D u[NPI] /**< [in] */,const double phi[NPI] /**< [in] */) const;
        
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
        void calcCorr(std::function<const Pt::pt3D (Nodes::Node)> getter,std::vector<double> &corr,Pt::pt3D u[NPI]) const;
        
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

    /** computes surface of the face */
    double calc_surf(void) const;
    
    private:
        const int NOD;/**< number of nodes */
        const int reg;/**< .msh region number */
        
        /** direct access to the Nodes */
        //const std::vector<Nodes::Node>  *refNode;
        const std::shared_ptr<Nodes::Node[]> refNode;
        
        /** computes weight coefficients */
        inline double weight(const int i) const { return 2.0*surf*Facette::pds[i]; }
};//end class Fac

}//end namespace

#endif /* facette_h */
