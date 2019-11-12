/** \file linear_algebra.h
\brief secondary header, it grabs altogether the linear algebra by the solver to apply fem method  <br>
It encapsulates the calls to GMM, the assemblage and projection of the matrix for all elements <br> 
projection and matrix assembly is multithreaded for tetrahedron, monothread for facette
*/
#include <deque>
#include <queue>
#include <thread>
#include <mutex>
#include <future>
#include <random>

#include "gmm/gmm_precond_diagonal.h"
#include "gmm/gmm_iter.h"
#include "gmm/gmm_solver_bicgstab.h"
#include "gmm/gmm_solver_gmres.h"

#include "config.h"

#include "feellgoodSettings.h"
#include "pt3D.h"
#include "tetra.h"
#include "facette.h"
#include "node.h"


#ifndef linear_algebra_h
#define linear_algebra_h


/** \class LinAlgebra
convenient class to grab altogether some part of the calculations involved using gmm solver at each timestep
*/
class LinAlgebra
{
public:
	/** constructor */	
    inline LinAlgebra(Settings & s /**< [in] */,
                      std::vector<Nodes::Node> & myNode /**< [in] */,
                      std::vector <Tetra::Tet> &myTet /**< [in] */,
                      std::vector <Facette::Fac> & myFace /**< [in] */) :  NbTH(s.solverNbTh),NOD(myNode.size()),refNode(&myNode),refFac(&myFace),refTet(&myTet),settings(s)
    {
    tab_TH.resize(NbTH+1);
    refTetIt.resize(NbTH);
    prepareItTet(myTet);
    base_projection(!RAND_DETERMINIST);
    }
    
    /** destructor */
    ~LinAlgebra ()
        { if (prc != nullptr) {delete prc;} }
    
	/** pointer to diagonal preconditionner  */
	gmm::diagonal_precond <read_matrix> *prc = nullptr;

    /** mono thread solver , debug only */ 
    int monoThreadSolver(timing const& t_prm,long nt);
    
    /**  solver, uses bicgstab and gmres, sparse matrix and vector are filled with multiThreading */
	int  solver(timing const& t_prm /**< [in] */,long nt /**< [in] */);

    /** setter for DW_dz */
    inline void set_DW_vz(double vz /**< [in] */){DW_vz = vz;}    

    /** getter for v_max */
    inline double get_v_max(void) {return v_max;}
    
private:
    /** number of threads, initialized by constructor */ 
    const int NbTH;
    
    const int NOD;/**< total number of nodes, also an offset for filling sparseMatrix, initialized by constructor */
    std::vector<Nodes::Node>  *refNode;/**< direct access to the Nodes */
	std::vector <Facette::Fac> *refFac; /**< direct access to the faces */
	std::vector <Tetra::Tet> *refTet; /**< direct access to the tet */
	
    /** vector of pair of iterators for the tetrahedrons for multithreading */
	std::vector < std::pair<std::vector<Tetra::Tet>::iterator,std::vector<Tetra::Tet>::iterator> > refTetIt; 
	
	double DW_vz;/**< speed of the domain wall */
	const Settings &settings;/**< settings */
    double v_max;/**< maximum speed */
	
	/** mutex to avoid improper access to sparse matrix */
    std::mutex my_mutex;
    
    /** thread vector */
    std::vector<std::thread> tab_TH;
    
    /** will be used to obtain a seed for the random number generator engine */
    std::random_device rd;
    
    /** update nodes.u and v values and find the maximum speed value  */
    void updateNodes(std::vector<double> const& X,const double dt);
    
    /** prepare refTetIt pairs, used by constructor */
    void prepareItTet(std::vector <Tetra::Tet> &myTet);
    
    /** computes the local vector basis {ep,eq} in the tangeant plane for projection on the elements */
    void base_projection(bool determinist);
    
    
    /** template function to provide P matrix coefficients, with respect to its block diagonal structure */
    template<class T,int N> double Pcoeff(T const& x,int i,int j) const
    {
    double val = 0;
    int node_i = i%N;
    
    /*
	double P[2*N][3*N] = { {0} }; // P must be filled with zero
	
    for (int i=0; i<N; i++){
  	  Nodes::Node const& n = (*refNode)[x.ind[i]];
	P[i][i]  = n.ep.x();  P[i][N+i]  = n.ep.y();  P[i][2*N+i]  = n.ep.z();
	P[N+i][i]= n.eq.x();  P[N+i][N+i]= n.eq.y();  P[N+i][2*N+i]= n.eq.z();
    	}
    */
    
    if(node_i == (j%N))
        {
        Nodes::Node const& n = (*refNode)[x.ind[node_i]];
            
        if(i<N)
            { val = n.ep(j/N); }
        else
            { val = n.calc_eq()(j/N); }
        }
    //std::cout <<"("<< i<<"; "<< j<<"):"<<val <<"\t";
    return val;
    }
    
    template<class T,int N> void projection_vect(T &x, Pt::pt3D B[N])
    {
        //tiny::mult<double,2*N,3*N>(P,B,x.Lp);

    for (int i=0; i<(2*N); i++)
        {
        x.Lp[i] = 0;
        for (int k=0; k<N; k++) { x.Lp[i] += Pcoeff<T,N>(x,i,k)*B[k].x(); }
        for (int k=0; k<N; k++) { x.Lp[i] += Pcoeff<T,N>(x,i,N+k)*B[k].y(); }
        for (int k=0; k<N; k++) { x.Lp[i] += Pcoeff<T,N>(x,i,2*N+k)*B[k].z(); }
        }
    x.treated = false;
    }
    
    
/** template to make projection for T, tetra or facette. It computes Ap = (P*A)*trans(P) and Bp = P*B and stores resuls in inner matrix Kp and vector Lp of class T*/
	template <class T,int N> void projection_mat(T &x,double A[3*N][3*N])//,  double B[3*N])
	{
    double PA[2*N][3*N]; // no need to initialize with zeros
//tiny::mult<double,2*N,3*N,3*N>(P,A,PA);
	for (int i=0; i<(2*N); i++) 
        {
        for (int k=0; k<(3*N); k++)
            {
            PA[i][k]=0;
            for (int j=0; j<(3*N); j++) { PA[i][k] += Pcoeff<T,N>(x,i,j)*A[j][k]; }
            }
        }
        
//tiny::mult<double,2*N,3*N>(P,B,x.Lp);
/*
    for (int i=0; i<(2*N); i++)
        {
        x.Lp[i] = 0;
        for (int k=0; k<(3*N); k++) { x.Lp[i] += Pcoeff<T,N>(x,i,k)*B[k]; }
        }*/

//tiny::direct_transposed_mult<double,2*N,3*N,2*N>(PA,P,x.Kp);
   for (int i=0; i<(2*N); i++) 
        for (int k=0; k<(2*N); k++)
        {
       x.Kp[i][k] = 0;
       for (int j=0; j<(3*N); j++) { x.Kp[i][k] += PA[i][j]* Pcoeff<T,N>(x,k,j); }
       }
    
    x.treated = false;
    }

    /** template to insert coeff in sparse matrix K_TH and vector L_TH, T is Tetra or Facette */
    template <class T> void insertCoeff(std::vector<T> & container, write_matrix &K_TH, std::vector<double> &L_TH)
    {
    std::for_each( container.begin(), container.end(), [&K_TH,&L_TH](T & my_elem)
        { if(!my_elem.treated) {my_elem.assemblage_mat(K_TH);my_elem.assemblage_vect(L_TH);my_elem.treated = true;} } ); 
    }
    
}; // fin class linAlgebra

#endif
