#include <algorithm>
#include <functional>
#include <mutex>

#include "Utils/FTic.hpp"
#include "linear_algebra.h"


int LinAlgebra::solver(timing const& t_prm,long nt)
{
FTic counter;

counter.tic();

write_matrix K_TH(2*NOD, 2*NOD);
std::vector<double> L_TH(2*NOD,0);

std::mutex my_mutex;

for(int i=0;i<NbTH;i++) 
    {
    tab_TH[i] = std::thread( [this,&my_mutex,&K_TH,&L_TH,i]() //,&Hext,&t_prm
        {
        std::for_each(refTetIt[i].first,refTetIt[i].second, [&my_mutex,&K_TH,&L_TH](Tetra::Tet & tet)
            {
            if(my_mutex.try_lock())
                {
                tet.assemblage_mat(K_TH);tet.assemblage_vect(L_TH);tet.treated = true;
                my_mutex.unlock();    
                }
            });//end for_each
        }); //end thread
    }

tab_TH[NbTH] = std::thread( [this,&my_mutex,&K_TH,&L_TH]()
    {
    std::for_each(refMsh->fac.begin(),refMsh->fac.end(), [&my_mutex,&K_TH,&L_TH](Facette::Fac & fac)
        {
        
        if(my_mutex.try_lock())
            {
            fac.assemblage_mat(K_TH);fac.assemblage_vect(L_TH);fac.treated =true;
            my_mutex.unlock();    
            }
        });
    }
);//end thread

for(int i=0;i<(NbTH+1);i++) {tab_TH[i].join();}

insertCoeff<Tetra::Tet>(refMsh->tet,K_TH,L_TH);
insertCoeff<Facette::Fac>(refMsh->fac,K_TH,L_TH);

counter.tac();

if(settings.verbose) { std::cout << "..matrix assembly done in " << counter.elapsed() << "s" << std::endl; }

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
	if(settings.verbose) { std::cout << "bicg FAILED in " << bicg_iter.get_iteration() << "iterations, duration: " << counter.elapsed() << " s" << std::endl; }
    gmm::diagonal_precond <read_matrix>  gmr_prc (Kr);
    counter.tic();
    
    gmm::clear(X);// X is cleared for safety, it should be useless
    gmm::gmres(Kr, X, L_TH, gmr_prc, 50, gmr_iter);

if (!(gmr_iter.converged() )) {
	counter.tac();
    if(settings.verbose) { std::cout << "gmres FAILED in " << gmr_iter.get_iteration() << "iterations, duration: " << counter.elapsed() << " s" << std::endl; }
    return 1;
    }
    else {counter.tac();
        if(settings.verbose) { std::cout << "v-solve in " << gmr_iter.get_iteration() << " iterations (gmres) duration: " << counter.elapsed() << " s" << std::endl; }
        }
    }
else 
    {counter.tac();
    if(settings.verbose) { std::cout << "v-solve converged in " << bicg_iter.get_iteration() << " (bicg,prc-" << (nt % (settings.REFRESH_PRC)) << ") :duration: " 
<< counter.elapsed() << " s" << std::endl; }
    }

updateNodes(X,t_prm.get_dt());
return 0;
}

