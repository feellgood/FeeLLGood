#include "Utils/FTic.hpp"
#include "linear_algebra.h"

int LinAlgebra::solver(timing const& t_prm,long nt)
{
FTic counter;

counter.tic();

write_matrix K_TH(2*NOD, 2*NOD);
std::vector<double> L_TH(2*NOD,0);

insertCoeff<Tetra::Tet>(refMsh->tet,K_TH,L_TH);
insertCoeff<Facette::Fac>(refMsh->fac,K_TH,L_TH);

counter.tac();

if(settings.verbose) { std::cout << "matrix assembly done in " << counter.elapsed() << " s.\n"; }

read_matrix Kr(2*NOD,2*NOD);
gmm::copy(K_TH, Kr);

std::vector<double> X(2*NOD);// filled with zero by default

gmm::iteration bicg_iter(1e-6);
gmm::iteration gmr_iter(1e-6);
bicg_iter.set_maxiter(settings.MAXITER);
gmr_iter.set_maxiter(settings.MAXITER);
bicg_iter.set_noisy(false);//VERBOSE
gmr_iter.set_noisy(false);//VERBOSE

if (!nt) 
    { prc = new gmm::diagonal_precond <read_matrix> (Kr); }
else if (!(nt % (settings.REFRESH_PRC)))
    {
    delete prc;
    prc = new gmm::diagonal_precond <read_matrix> (Kr);
    } 

counter.tic();
//gmm::clear(X); // no need to clear X,it is zero initialized 
gmm::bicgstab(Kr, X, L_TH, *prc, bicg_iter);

if (!(bicg_iter.converged() )) 
    {
    counter.tac();
	if(settings.verbose)
        {
        std::cout << "solver: bicg (prc "
                  << (nt % (settings.REFRESH_PRC)) << ") FAILED after "
                  << bicg_iter.get_iteration() << " iterations, in "
                  << counter.elapsed() << " s.\n";
        }
    gmm::diagonal_precond <read_matrix>  gmr_prc (Kr);
    counter.tic();
    
    gmm::clear(X);// X is cleared for safety, it should be useless
    gmm::gmres(Kr, X, L_TH, gmr_prc, 50, gmr_iter);

if (!(gmr_iter.converged() )) {
	counter.tac();
    if(settings.verbose)
        {
        std::cout << "solver: gmres FAILED after "
                  << gmr_iter.get_iteration() << " iterations, in "
                  << counter.elapsed() << " s.\n";
        }
    return 1;
    }
    else {counter.tac();
        if(settings.verbose)
            {
            std::cout << "solver: gmres converged in "
                      << gmr_iter.get_iteration() << " iterations, "
                      << counter.elapsed() << " s.\n";
            }
        }
    }
else 
    {counter.tac();
    if(settings.verbose)
        {
        std::cout << "solver: bicg (prc "
                  << (nt % (settings.REFRESH_PRC)) << ") converged in "
                  << bicg_iter.get_iteration() << " iterations, "
                  << counter.elapsed() << " s.\n";
        }
    }

updateNodes(X,t_prm.get_dt());
return 0;
}

