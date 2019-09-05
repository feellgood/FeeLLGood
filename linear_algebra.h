/** \file linear_algebra.h
\brief secondary header, it grabs altogether the linear algebra by the solver to apply fem method  <br>
It encapsulates the calls to MTL4 , the assemblage and projection of the matrix for all elements <br> 
projection and matrix assembly is multithreaded for tetrahedron, monothread for facette
*/
#include <deque>
#include <queue>
#include <thread>
#include <mutex>
#include <future>

#include "gmm/gmm_precond_diagonal.h"
#include "gmm/gmm_iter.h"
#include "gmm/gmm_solver_bicgstab.h"
#include "gmm/gmm_solver_gmres.h"

#include "fem.h"

#ifndef linear_algebra_h
#define linear_algebra_h


typedef gmm::wsvector <double>   write_vector;/**< gmm write vector build on std::map, log(n) for read and write access */
typedef gmm::rsvector <double>   read_vector; /**< gmm read vector */

typedef gmm::row_matrix	<write_vector>   write_matrix; /**< gmm write sparse matrix */
typedef gmm::row_matrix	<read_vector>    read_matrix; /**< gmm read sparse matrix */

/** \class LinAlgebra
convenient class to grab altogether some part of the calculations involved using gmm solver at each timestep
*/
class LinAlgebra
{
public:
	/** constructor */	
    inline LinAlgebra(Settings & s /**< [in] */,const int _NOD /**< [in] */,
                      std::vector<Nodes::Node> & myNode /**< [in] */,
                      std::vector <Tetra::Tet> & myTet /**< [in] */,
                      std::vector <Facette::Fac> & myFace /**< [in] */) :  NOD(_NOD),refNode(&myNode),refFac(&myFace)
    {
    settings = s;
    NbTH = s.solverNbTh;
    my_lock = new std::mutex;
    tab_TH.resize(NbTH+1);
    refTet.resize(NbTH);
    const unsigned long block_size = std::distance(myTet.begin(),myTet.end())/NbTH;

    std::vector<Tetra::Tet>::iterator it_begin = myTet.begin();

    for(int i=0;i<(NbTH-1);i++) 
        {
        std::vector<Tetra::Tet>::iterator it_end = it_begin;
        std::advance(it_end,block_size);
        refTet[i].resize(block_size,Tetra::Tet(_NOD));
        std::copy( it_begin, it_end, refTet[i].begin() );
        it_begin = it_end;
        }
    const unsigned long last_block_size = std::distance(it_begin,myTet.end());
    refTet[NbTH-1].resize(last_block_size,Tetra::Tet(_NOD));
    std::copy( it_begin, myTet.end(), refTet[NbTH-1].begin() );
    
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
    inline double get_v_max() {return v_max;}
    
    /** getter node */
    inline Nodes::Node getNode(int i /**< [in] */) {return (*refNode)[i];}
    
    /** getter node physical position */
    inline Pt::pt3D getNodePhysPos(int i /**< [in] */) {return (*refNode)[i].p;} 
    
    
private:
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
    
    /** number of threads, initialized by constructor */ 
    int NbTH;
    
    /** buffer fr facette::Obj when try_lock fail */
    std::queue<Facette::Obj> buff_fac;
    
    /** thread vector */
    std::vector<std::thread> tab_TH;
    
/** computes the local vector basis {ep,eq} in the tangeant plane for projection on the elements */
inline void base_projection(void)
	{ std::for_each(refNode->begin(),refNode->end(),[](Nodes::Node &n) { n.buildBase_epeq();}); }
	
}; // fin class linAlgebra

#endif
