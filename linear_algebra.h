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
    
	/** pointer to diagonal preconditionner  */
	gmm::diagonal_precond <read_matrix> *prc;

    /** solver, uses bicgstab and gmres */
	int  solver(long nt /**< [in] */);

    /** setter for dt */
    inline void set_dt(double _dt /**< [in] */){dt = _dt;}
    
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
	
	double dt;/**< timestep */
    double DW_vz;/**< speed of the domain wall */
	const Settings &settings;/**< settings */
    double v_max;/**< maximum speed */
	
	/** mutex to avoid improper access to sparse matrix */
    std::mutex my_mutex;
    
    /** buffer fr facette::Obj when try_lock fail */
    std::queue<Facette::Obj> buff_fac;
    
    /** thread vector */
    std::vector<std::thread> tab_TH;
    
    /** will be used to obtain a seed for the random number generator engine */
    std::random_device rd;
    
    /** prepare refTetIt pairs, used by constructor */
    void prepareItTet(std::vector <Tetra::Tet> &myTet);
    
    /** computes the local vector basis {ep,eq} in the tangeant plane for projection on the elements */
    void base_projection(bool determinist);
    
    template <class T> void insertCoeff(std::vector<T> container, write_matrix &K_TH, write_vector &L_TH)
    {
    std::for_each( container.begin(), container.end(), [&K_TH,&L_TH](T & my_elem)
        {
        if(!my_elem.treated) {my_elem.assemblage_mat(K_TH);my_elem.assemblage_vect(L_TH);my_elem.treated = true;} 
        } ); 
    }
    
}; // fin class linAlgebra

#endif
