/** \file electrostatSolver.h
  \brief solver for electrostatic problem when STT is required
  header containing electrostatSolver class. It uses biconjugate stabilized gradient with diagonal preconditioner. The solver is only called once to compute voltages V for each nodes of the mesh, when STT computation is involved.
 */

#include <map>

#include "gmm/gmm_iter.h"
#include "gmm/gmm_solver_bicgstab.h"
#include "gmm/gmm_solver_cg.h"


#include "fem.h"
#include "config.h"
#include "mesh.h"

#include "tetra.h"
#include "facette.h"
#include "tiny.h"

#include "spinTransferTorque.h"

/** \class electrostatSolver
this class is containing both data and a solver to compute potential from dirichlet boundary conditions problem for the current density flowing in the sample.
*/
class electrostatSolver {
public:
    /** constructor */
    inline electrostatSolver(Mesh::mesh const& _msh /**< [in] reference to the mesh */, STT const& _p_stt /**< all spin transfer torque parameters */,
                             const double _tol /**< [in] tolerance for solvers */,
                             const bool v /**< [in] verbose bool */,
                             const int max_iter /**< [in] maximum number of iteration */ ): msh(_msh), p_stt(_p_stt), verbose(v), MAXITER(max_iter) 
                             {
                             if(verbose) { std::cout << "Dirichlet boundary conditions..." << std::endl; infos(); }
                             solve(_tol);
                             }

/** computes the gradient(V) for tetra tet */
    void calc_gradV(Tetra::Tet const& tet,Pt::pt3D (&gradV)[Tetra::NPI]) const
    	{
    	if(V.size()>0)
		{
		for (int npi=0; npi<Tetra::NPI; npi++) 
    			{ 
    			double vx(0),vy(0),vz(0);
    
    			for (int i=0; i<Tetra::N; i++)
        			{
        			const double V_tet = V[ tet.ind[i] ];
        			vx += V_tet*tet.dadx[i][npi];
        			vy += V_tet*tet.dady[i][npi];
        			vz += V_tet*tet.dadz[i][npi];
        			}
    
    			gradV[npi] = Pt::pt3D(vx,vy,vz);
    			}
		}
//tiny::transposed_mult<double, N, NPI> (V_nod, dadx, dVdx);
//tiny::transposed_mult<double, N, NPI> (V_nod, dady, dVdy);
//tiny::transposed_mult<double, N, NPI> (V_nod, dadz, dVdz);
    	}

void calc_Hm_STT(Tetra::Tet const& tet, Pt::pt3D (&gradV)[Tetra::NPI], Pt::pt3D (&Hm)[Tetra::NPI]) const
{
Pt::pt3D p_g[Tetra::NPI];
tet.interpolation(Nodes::get_p,p_g);

for (int npi=0; npi<Tetra::NPI; npi++)
    { Hm[npi] = -p_stt.sigma*gradV[npi]*p_g[npi]; }
}

/** electrostatic potential values for boundary conditions, V.size() is the size of the vector of nodes */ 
    std::vector<double> V;

private:
    /** mesh object to store nodes, fac, tet, and others geometrical values related to the mesh ( const ref ) */
	Mesh::mesh msh;
	
	/** spin transfer torque parameters */
	STT p_stt;
    
    
    /** if verbose set to true, some printing are sent to terminal */
    const bool verbose;

    /** maximum number of iteration for biconjugate stabilized gradient */
    const int MAXITER; //fixed to 5000 in ref code
    
    /** basic informations on boundary conditions */
    inline void infos(void) 
        {
        std::cout << "sigma: " << p_stt.sigma << std::endl;
        
        std::for_each(p_stt.boundaryCond.begin(),p_stt.boundaryCond.end(),
                      [](std::pair<std::string,double> const& p)
                      { std::cout << "regName: " << p.first << "\tV :" << p.second << std::endl; } );
        }
        
    /** assemble the matrix K from tet and Ke inputs */
    inline void assembling_mat(Tetra::Tet const& tet, gmm::dense_matrix <double> &Ke, write_matrix &K)
    {
    for (int ie=0; ie<Tetra::N; ie++)
        { for (int je=0; je<Tetra::N; je++) { K(tet.ind[ie], tet.ind[je]) += Ke(ie, je); } }
    }
    
    /** assemble the vector L from fac and Le inputs */
    inline void assembling_vect(Facette::Fac const& fac, std::vector <double> &Le, write_vector &L)
    { for (int ie=0; ie<Facette::N; ie++) { L[ fac.ind[ie] ] += Le[ie]; } }
    
    
    /** compute integrales for matrix coefficients,input from tet ; sigma is the region conductivity */
    void integrales(Tetra::Tet const& tet, double sigma, gmm::dense_matrix <double> &AE)
    {
    
    
        for (int npi=0; npi<Tetra::NPI; npi++)
            {
            double w, dai_dx, dai_dy, dai_dz, daj_dx, daj_dy, daj_dz;
            w = tet.weight[npi];

            for (int ie=0; ie<Tetra::N; ie++)
                {
                dai_dx= tet.dadx[ie][npi];  dai_dy= tet.dady[ie][npi];  dai_dz= tet.dadz[ie][npi];

                for (int je=0; je<Tetra::N; je++)
                    {
                    daj_dx= tet.dadx[je][npi];  daj_dy= tet.dady[je][npi];  daj_dz= tet.dadz[je][npi];
                    AE(ie, je) += sigma*(dai_dx*daj_dx + dai_dy*daj_dy + dai_dz*daj_dz)*w;		
                    }
                }
            }
    }
    
    /** compute integrales for vector coefficients, input from facette */
void integrales(Facette::Fac const& fac,double pot_val, std::vector <double> &BE)
    {
    for (int npi=0; npi<Facette::NPI; npi++)
            { for (int ie=0; ie<Facette::N; ie++) { BE[ie] -= Facette::a[ie][npi]*pot_val*fac.weight(npi); } }
    }

    /** fill matrix and vector to solve potential values on each node */
void prepareData(write_matrix &Kw, write_vector & Lw)
{
for (unsigned int ne=0; ne < msh.tet.size(); ne++){
    Tetra::Tet const& tet = msh.tet[ne];
    gmm::dense_matrix <double> K(Tetra::N, Tetra::N);
    double sigma = p_stt.sigma;
    integrales(tet,sigma, K);
    assembling_mat(tet, K, Kw);
    }

for (unsigned int ne=0; ne < msh.fac.size(); ne++){
    Facette::Fac const& fac = msh.fac[ne];
    std::vector <double> L(Facette::N);
    double pot_val =0; // to find in V_values
    integrales(fac, pot_val, L);     
    assembling_vect(fac, L, Lw);
    }    
}

/** solver, using biconjugate stabilized gradient, with diagonal preconditionner */
int solve(const double iter_tol)
{
const int NOD = msh.getNbNodes();

write_matrix Kw(NOD, NOD);
write_vector Lw(NOD);
write_vector Xw(NOD);

prepareData(Kw,Lw);

read_matrix  Kr(NOD, NOD);    gmm::copy(Kw, Kr);

//Dirichlet boundary conditions
std::for_each( p_stt.boundaryCond.begin(),p_stt.boundaryCond.end(), [this,&Kw,&Lw,&Xw](auto const& it) 
	{ 
	double V = it.second;
	std::string name = it.first;
	
	auto surf_it = std::find_if(msh.s.begin(),msh.s.end(), [name] (Mesh::Surf const& _s) { return (_s.getName() == name); }  );
	
	if(surf_it != msh.s.end())
		{
		std::for_each( surf_it->elem.begin(),surf_it->elem.end(), [V,&Kw,&Lw,&Xw](Mesh::Triangle const& tri)
			{
			for (int ie=0; ie<Facette::N ; ie++) // should be Triangle::N
				{
				const int i= tri.ind[ie];
				std::for_each(mat_row(Kw, i).begin(),mat_row(Kw, i).end(), [](std::pair<const long unsigned int, double> & _it) { _it.second = 0.0; } );
            			Kw(i, i) =  1e9;
            			Lw[i]    =  V*1e9;  
            			Xw[i]    =  V;
				}
			});
		}
	} );
/*
for (unsigned int ne=0; ne<msh.fac.size(); ne++)
    {
    Facette::Fac const& fac = msh.fac[ne];
    std::map<std::string,double>::iterator it = V_values.find(fac.idxPrm); // might not be idxPrm 
        
    if (it != V_values.end() ) 
        {
        double V = it->second;
        for (int ie=0; ie<Facette::N; ie++)
            {
            const int i= fac.ind[ie];
            std::for_each(mat_row(Kw, i).begin(),mat_row(Kw, i).end(), [](std::pair<const long unsigned int, double> & it) { it.second = 0.0; } );
            Kw(i, i) =  1e9;
            Lw[i]    =  V*1e9;  
            Xw[i]    =  V;
            }
        }
    }
*/

if(verbose) { std::cout << "line weighting..." << std::endl; }
/* equilibrage des lignes */
for (int i=0; i<NOD; i++)
    {
    double norme=gmm::vect_norminf(mat_row(Kw, i));
    Lw[i]/=norme;
    std::for_each(mat_row(Kw, i).begin(),mat_row(Kw, i).end(), [norme](std::pair<const long unsigned int, double> & it) { it.second/=norme; } );
    }

gmm::copy(Kw, Kr);
read_vector  Lr(NOD);        gmm::copy(Lw, Lr);

if(verbose) { std::cout << "solving ..." << std::endl; }

gmm::iteration iter(iter_tol);
iter.set_maxiter(MAXITER);
iter.set_noisy(verbose);

gmm::bicgstab(Kr, Xw, Lr, gmm::diagonal_precond <read_matrix>(Kr), iter);

V.resize(NOD);
gmm::copy(Xw, V);
return 0;
}

}; //end class electrostatSolver
