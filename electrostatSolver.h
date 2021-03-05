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

/** \class electrostatSolver
this class is containing both data and a solver to compute potential from dirichlet boundary conditions problem for the current density flowing in the sample.
*/
class electrostatSolver {
public:
    /** constructor */
    inline electrostatSolver(mesh const& _msh /**< [in] reference to the mesh */, Tetra::STT const& p_stt /**< all spin transfer torque parameters */,
                             const double _tol /**< [in] tolerance for solvers */,
                             const bool v /**< [in] verbose bool */,
                             const bool fastSolve /**< [in] if true fast_solve is used */,
                             const int max_iter /**< [in] maximum number of iteration */ ): msh(_msh), NOD(msh.getNbNodes()), TET(msh.getNbTets()), FAC(msh.getNbFacs()), verbose(v), MAXITER(max_iter) 
                             {
                             if (p_stt.bc != Tetra::boundary_conditions::Undef)
                                {
                                sigma_values.insert( std::pair<int,double>(p_stt.reg,p_stt.sigma) );
                                
                                std::for_each(p_stt.full_bc_data.begin(),p_stt.full_bc_data.end(), 
                                              [this](Tetra::bc_data const& bc_d)
                                                {
                                                switch(bc_d.typ)
                                                    {
                                                    case Tetra::type_val_reg::potV: 
                                                        V_values.insert(std::pair<int,double>(bc_d.reg,bc_d.val) );
                                                    break;
                                                    case Tetra::type_val_reg::densJ:
                                                        J_values.insert(std::pair<int,double>(bc_d.reg,bc_d.val) );
                                                    break;
                                                    default: std::cout << "type error for physical constant in boundary conditions" << std::endl; exit(1); break;
                                                    }
                                                    
                                                } 
                                             ); //end for_each
                                }
                            else
                                { std::cout << "warning : undefined boundary conditions for STT" << std::endl; exit(1);}
                            if(verbose) { infos(); }
                            
                            if (fastSolve)
                                {
                                Vd.resize(NOD);
                                
                                for (int ne=0; ne<FAC; ne++)
                                    {
                                    Facette::Fac const& fac = msh.fac[ne]; 
                                    std::map<int,double>::iterator it = V_values.find( fac.getRegion() ); 
                                    if (it != V_values.end() ) 
                                        { // found
                                        for (int ie=0; ie<Facette::N; ie++)
                                            {
                                            int i= fac.ind[ie];  
                                            Vd[i]    =  it->second;
                                            ld.push_back(i);
                                            }
                                        }
                                    }
                                
                                std::sort( ld.begin(), ld.end() );
                                ld.erase( std::unique(ld.begin(), ld.end() ) , ld.end() );
                                
                                dofs.resize(NOD);
                                std::vector<size_t> all(NOD);
                                std::iota (all.begin(), all.end(), 0); // Fill with 0, 1, ..., NOD-1.
                                std::vector<size_t>::iterator it = std::set_difference (all.begin(), all.end(), ld.begin(), ld.end(), dofs.begin());
                                dofs.resize(it-dofs.begin());
                                    
                                fast_solve(_tol);
                                }
                            else solve(_tol);
                            }

private:
    /** mesh object to store nodes, fac, tet, and others geometrical values related to the mesh */
	mesh msh;
    
    /** number of nodes */
    const int NOD;
    /** number of tet */
    const int TET;
    /** number of fac */
    const int FAC;
    /** if verbose set to true, some printing are sent to terminal */
    const bool verbose;

    /** maximum number of iteration for biconjugate stabilized gradient */
    const int MAXITER; //fixed to 5000 in ref code
    
    std::map<int,double> sigma_values;/**< conductivity region volume table */
    
    /** boundary conditions : table of current densities (int is a surface region) */
    std::map<int,double> J_values;
    
    /** boundary conditions : table of voltage (int is a surface region) */
    std::map<int,double> V_values; 
    
/* ------ data structure for fast solver --- */
    /** freedom degree list for the potential */
    std::vector <size_t> dofs;
    
    /** dirichlet index node list */
    std::vector <size_t> ld;
    
    /** Dirichlet potential values at the nodes, 0 elsewhere */
    std::vector <double> Vd;
/* ------ end data structure for fast solver --- */
    
    /** basic informations on boundary conditions */
    inline void infos(void) 
        {
        std::for_each(sigma_values.begin(),sigma_values.end(),
                      [](std::pair<int,double> const& p)
                      { std::cout << "reg: " << p.first << "\tsigma :" << p.second << std::endl; } );
        
        std::for_each(J_values.begin(),J_values.end(),
                      [](std::pair<int,double> const& p)
                      { std::cout << "reg: " << p.first << "\tJ :" << p.second << std::endl; } );
        
        std::for_each(V_values.begin(),V_values.end(),
                      [](std::pair<int,double> const& p)
                      { std::cout << "reg: " << p.first << "\tV :" << p.second << std::endl; } );
        }
    
    /** sigma value getter */
    inline double getSigma(const int reg)
        {
        double val(0);
        std::map<int,double>::iterator it = sigma_values.find(reg); 
        if (it != sigma_values.end() ) val = it->second;
        
        return val;
        }
        
    /** Jn value getter */
    inline double getJ(const int reg)
        {
        double val(0);
        std::map<int,double>::iterator it = J_values.find(reg); 
        if (it != J_values.end() ) val = it->second;
        
        return val;
        }
    
    /** potential value getter */
    inline double getV(const int reg)
        {
        double val(0);
        std::map<int,double>::iterator it = V_values.find(reg); 
        if (it != V_values.end() ) val = it->second;
        
        return val;
        }
    
    /** assemble the matrix K from tet and Ke inputs */
    void assembling_mat(Tetra::Tet const& tet, gmm::dense_matrix <double> &Ke, write_matrix &K)
    {
    for (int ie=0; ie<Tetra::N; ie++)
        {
        const int i= tet.ind[ie];             
        for (int je=0; je<Tetra::N; je++) { K(i, tet.ind[je]) += Ke(ie, je); }
        }
    }
    
    /** assemble the vector L from fac and Le inputs */
    void assembling_vect(Facette::Fac const& fac, std::vector <double> &Le, write_vector &L)
    {
    for (int ie=0; ie<Facette::N; ie++) { L[ fac.ind[ie] ] += Le[ie]; }
    }
    
    
    /** compute integrales for matrix coefficients,input from tet */
    void integrales(Tetra::Tet const& tet, gmm::dense_matrix <double> &AE)
    { //sigma is the region conductivity
    double sigma = getSigma(tet.getRegion());

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
void integrales(Facette::Fac const& fac, std::vector <double> &BE)
{
double J = getJ(fac.getRegion());

for (int npi=0; npi<Facette::NPI; npi++)
    { for (int ie=0; ie<Facette::N; ie++) { BE[ie] -= Facette::a[ie][npi]*J*fac.weight(npi); } }
}

/** solver, using biconjugate stabilized gradient, with diagonal preconditionner */
int solve(const double iter_tol)
{
write_matrix Kw(NOD, NOD);
write_vector Lw(NOD);
write_vector Xw(NOD);

if(verbose) { std::cout << "assembling..." << std::endl; }

for (int ne=0; ne<TET; ne++){
    Tetra::Tet const& tet = msh.tet[ne];
    gmm::dense_matrix <double> K(Tetra::N, Tetra::N);
    integrales(tet, K);
    assembling_mat(tet, K, Kw);
    }

for (int ne=0; ne<FAC; ne++){
    Facette::Fac const& fac = msh.fac[ne];
    std::vector <double> L(Facette::N);
    integrales(fac, L);     
    assembling_vect(fac, L, Lw);
    }

read_matrix  Kr(NOD, NOD);    gmm::copy(Kw, Kr);

if(verbose) { std::cout << "boundary conditions..." << std::endl; }
// conditions de Dirichlet
for (int ne=0; ne<FAC; ne++)
    {
    Facette::Fac const& fac = msh.fac[ne];

    std::map<int,double>::iterator it = V_values.find(fac.getRegion()); 
        
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

read_vector Xr(NOD); gmm::copy(Xw, Xr);

msh.setNodesPotential(Xr);
return 0;
}
    
/* ----------------------------------------------------- */
/* experimental Ohms'law solver */
/* ----------------------------------------------------- */
    
    
    
    
    
/** other solver (using reduced sub space), using conjugate gradient, with diagonal preconditionner with mixt boundary conditions */
int fast_solve(const double iter_tol)
{ 
write_matrix Kw(NOD, NOD);
write_vector Lw(NOD);

FTic counter;

counter.tic();
for (int ne=0; ne<TET; ne++){
    Tetra::Tet const& tet = msh.tet[ne];
    gmm::dense_matrix <double> K(Tetra::N, Tetra::N);
    integrales(tet, K);   
    assembling_mat(tet, K, Kw);
    }

for (int ne=0; ne<FAC; ne++){
    Facette::Fac const& fac = msh.fac[ne];
    std::vector <double> L(Facette::N);
    integrales(fac, L);     
    assembling_vect(fac, L, Lw);
    }
counter.tac();
std::cout << "\t integrales & assembling total computing time: " << counter.elapsed() << " s\n" << std::endl;


counter.tic();
// Modification du second membre pour tenir compte des valeurs du potentiel aux noeuds de dirichlet
    {
    read_vector  Lr(NOD);
    gmm::copy(Lw,  Lr);

    read_matrix  Kr(NOD, NOD);  
    gmm::copy(Kw,  Kr);

    gmm::mult(Kr, gmm::scaled(Vd, -1.0), Lr, Lw);  // Kr * (-Vd) + Lr --> Lw
    }

write_matrix subKw(dofs.size(), dofs.size());      
gmm::copy(gmm::sub_matrix(Kw, gmm::sub_index(dofs)    , gmm::sub_index(dofs))    , subKw);
read_matrix  Kr(dofs.size(), dofs.size());
gmm::copy(subKw,  Kr);

write_vector subLw(dofs.size());
gmm::copy(gmm::sub_vector(Lw, gmm::sub_index(dofs))    , subLw);
read_vector  Lr(dofs.size());
gmm::copy(subLw,  Lr);

write_vector Xw(dofs.size());
counter.tac();
std::cout << "\t many copies with subindices vector total computing time: " << counter.elapsed() << " s\n" << std::endl;



counter.tic();
std::cout << "\t solving .......................... " << std::endl;

gmm::iteration iter(iter_tol);
iter.set_maxiter(MAXITER);
iter.set_noisy(verbose);

gmm::identity_matrix PS;   // Optional scalar product for cg

gmm::cg(Kr, Xw, Lr, PS, gmm::diagonal_precond <read_matrix>(Kr), iter); // Conjugate gradient

read_vector Xr(dofs.size());
gmm::copy(Xw, Xr);

for (size_t ip=0; ip<dofs.size(); ip++)
    { msh.set_elec_pot( dofs[ip] , Xr[ip] ); }


for (size_t ip=0; ip<ld.size(); ip++)
    { msh.set_elec_pot( ld[ip] , Vd[ ld[ip] ] ); }

counter.tac();
std::cout << "\t total cg solving time: " << counter.elapsed() << " s\n" << std::endl;
    
return 0;
}

}; //end class electrostatSolver
