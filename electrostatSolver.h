/** \file electrostatSolver.h
  \brief solver for electrostatic problem when STT is required
  header containing electrostatSolver class. It uses biconjugate stabilized gradient with diagonal preconditioner. The solver is only called once to compute voltages V for each nodes of the mesh, when STT computation is involved.
 */

#include <map>

//#include "gmm/gmm.h"
//#include "gmm/gmm_MUMPS_interface.h"
//#include "gmm/gmm_precond_ilu.h"
#include "gmm/gmm_iter.h"
#include "gmm/gmm_solver_bicgstab.h"

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
                             const bool v /**< [in] verbose bool */,
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
                            
                            std::cout << "new implementation of spin transfer torque not yet supported, coming soon..." << std::endl; 
                            infos();
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
    
    /** template to assemble the matrix and vector, T is Tet or Fac */
    template <class T>
    void assembling(T const& obj, gmm::dense_matrix <double> &Ke, std::vector <double> &Le, write_matrix &K, write_vector &L)
    {
    const int N = obj.getN();
    for (int ie=0; ie<N; ie++)
        {
        int i= obj.ind[ie];             
        for (int je=0; je<N; je++)
            {
            int j= obj.ind[je];
            K(i, j)+= Ke(ie, je);
            }
        L[i]+= Le[ie];
        }
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
void integrales(Facette::Fac &fac, std::vector <double> &BE)
{
double J = getJ(fac.getRegion());

for (int npi=0; npi<Facette::NPI; npi++)
    {
    double w =fac.weight(npi);
    for (int ie=0; ie<Facette::N; ie++) { BE[ie] -= Facette::a[ie][npi]*J *w; }
    }
}

/** solver, using biconjugate stabilized gradient, with diagonal preconditionner */
int solve(void)
{
write_matrix Kw(NOD, NOD);
write_vector Lw(NOD);
write_vector Xw(NOD);

if(verbose) { std::cout << "assembling..." << std::endl; }

for (int ne=0; ne<TET; ne++){
    Tetra::Tet &tet = msh.tet[ne];
    gmm::dense_matrix <double> K(Tetra::N, Tetra::N);
    std::vector <double> L(Tetra::N);
    integrales(tet, K);
    assembling<Tetra::Tet>(tet, K, L, Kw, Lw);
    }

for (int ne=0; ne<FAC; ne++){
    Facette::Fac &fac = msh.fac[ne];
    gmm::dense_matrix <double> K(Facette::N, Facette::N);
    std::vector <double> L(Facette::N);
    integrales(fac, L);     
    assembling<Facette::Fac>(fac, K, L, Kw, Lw);
    }

read_matrix  Kr(NOD, NOD);    gmm::copy(Kw, Kr);

// conditions de Dirichlet
for (int ne=0; ne<FAC; ne++)
    {
    Facette::Fac &fac = msh.fac[ne];
    const int reg = fac.getRegion();
    
    double V = getV(reg);
       for (int ie=0; ie<Facette::N; ie++)
            {
            const int i= fac.ind[ie];
            auto row = mat_const_row(Kw, i);
            //gmm::linalg_traits<gmm::row_matrix<write_vector > >::const_sub_row_type row = mat_const_row(Kw, i);
            std::for_each(row.begin(),row.end(),[&Kw,i]( std::pair<const long unsigned int, double>& it){ Kw(i,it.first) = 0.0; } );
            //for (write_vector::const_iterator it=vect_const_begin(row); it!=vect_const_end(row); ++it) Kw(i, it->first)=0;
            Kw(i, i) =  1;
            Lw[i]    =  V;  
            Xw[i]    =  V;
            }
    }

/* equilibrage des lignes */
for (int i=0; i<NOD; i++)
    {
    auto row = mat_const_row(Kw, i);
    //gmm::linalg_traits<gmm::row_matrix<write_vector > >::const_sub_row_type row = mat_const_row(Kw, i);
    double norme=gmm::vect_norminf(row);
    Lw[i]/=norme;
    
    std::for_each(row.begin(),row.end(),[&Kw,norme,i]( std::pair<const long unsigned int, double>& it){ Kw(i,it.first)/=norme; } );
    //for (write_vector::const_iterator it=vect_const_begin(row); it!=vect_const_end(row); ++it) Kw(i, it->first)/=norme;
    }

gmm::copy(Kw, Kr);
read_vector  Lr(NOD);        gmm::copy(Lw, Lr);

if(verbose) { std::cout << "preconditionning ..." << std::endl; }

gmm::diagonal_precond <read_matrix> prc(Kr);
//gmm::ilu_precond <read_matrix> prc(Kr);

if(verbose) { std::cout << "solving ..." << std::endl; }

gmm::iteration iter(1e-8);
iter.set_maxiter(MAXITER);
iter.set_noisy(verbose);

gmm::bicgstab(Kr, Xw, Lr, prc, iter);

read_vector Xr(NOD); gmm::copy(Xw, Xr);

for (int i=0; i<NOD; i++) { msh.setNode(i).V = Xr[i]; }

return 0;
}
    
    
}; //end class electrostatSolver
