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
    inline LinAlgebra(Settings & s /**< [in] */,const int _NOD /**< [in] */,
                      std::vector<Nodes::Node> & myNode /**< [in] */,
                      std::vector <Tetra::Tet> const& myTet /**< [in] */,
                      std::vector <Facette::Fac> & myFace /**< [in] */) :  NbTH(s.solverNbTh),NOD(_NOD),refNode(&myNode),refFac(&myFace)
    {
    settings = s;
    my_lock = new std::mutex;
    tab_TH.resize(NbTH+1);
    refTet.resize(NbTH);
    deepCopyTet(myTet);
    if(VERBOSE) { std::cout << NbTH+1 << " threads for assembling matrix." << std::endl; }
    }
    
	/** pointer to diagonal preconditionner  */
	gmm::diagonal_precond <read_matrix> *prc;

    /** solver, uses bicgstab and gmres */
	int  solver(long nt /**< [in] */);

    /** setter for dt */
    inline void set_dt(double _dt /**< [in] */){dt = _dt;}
    
    /** setter for DW_dz */
    inline void set_DW_vz(double vz /**< [in] */){DW_vz = vz;}    

    /** setter for Hext */
    inline void set_Hext(double Hx /**< [in] */,double Hy /**< [in] */,double Hz /**< [in] */){Hext[0]=Hx;Hext[1]=Hy;Hext[2]=Hz;}

    /** getter for v_max */
    inline double get_v_max(void) {return v_max;}
    
private:
    /** number of threads, initialized by constructor */ 
    const int NbTH;
    
    const int NOD;/**< total number of nodes, also an offset for filling sparseMatrix, initialized by constructor */
    std::vector<Nodes::Node>  *refNode;/**< direct access to the Nodes */
	std::vector <Facette::Fac> *refFac; /**< direct access to the faces */
	
	std::vector < std::vector <Tetra::Tet> > refTet; /**< splitted copy of the tetrahedrons for multithreading */
	
	double Hext[DIM];/**< applied field */
    double dt;/**< timestep */
    double DW_vz;/**< speed of the domain wall */
	Settings settings;/**< settings */
    double v_max;/**< maximum speed */
	
	/** mutex to avoid improper access to sparse matrix */
    std::mutex *my_lock;
    
    /** buffer fr facette::Obj when try_lock fail */
    std::queue<Facette::Obj> buff_fac;
    
    /** thread vector */
    std::vector<std::thread> tab_TH;
    
    /** will be used to obtain a seed for the random number generator engine */
    std::random_device rd;
    
    /** deep copy of all the tetrahedrons to refTet container, used by constructor */
    void deepCopyTet(std::vector <Tetra::Tet> const& myTet);
    
    /** computes the local vector basis {ep,eq} in the tangeant plane for projection on the elements */
    void base_projection(void);
	
}; // fin class linAlgebra

#endif
