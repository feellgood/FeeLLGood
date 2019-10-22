#include <algorithm>
#include <functional>

#include "Utils/FTic.hpp"
#include "linear_algebra.h"
    
int LinAlgebra::solver(timing const& t_prm,long nt)
{
FTic counter;

counter.tic();

base_projection(!RAND_DETERMINIST);
write_matrix K_TH(2*NOD, 2*NOD);
double L_TH[2*NOD] = {0};

std::cout << ", starting matrix assembly...";

for(int i=0;i<NbTH;i++) 
    {
    tab_TH[i] = std::thread( [this,&t_prm,&K_TH,&L_TH,i]() 
        {
        std::for_each(refTetIt[i].first,refTetIt[i].second, [this,&t_prm,&K_TH,&L_TH](Tetra::Tet & tet)
            {
            //gmm::dense_matrix <double> K(3*Tetra::N,3*Tetra::N); 
            double K[3*Tetra::N][3*Tetra::N] = { {0} }; 
                //std::vector <double> L(3*Tetra::N);
                double L[3*Tetra::N] = {0};
                
            tet.integrales(settings.paramTetra,t_prm.dt,settings.Hext,DW_vz,t_prm.TAUR,K, L);     
            tet.projection( K, L);
            tet.treated = false;
            if(my_mutex.try_lock())
                {
                tet.assemblage_mat(K_TH);
                tet.assemblage_vect(L_TH);
                tet.treated = true;
                my_mutex.unlock();    
                }
            else { tet.treated = false;}
            });//end for_each
        }); //end thread
    }
    
tab_TH[NbTH] = std::thread( [this,&K_TH,&L_TH]()
    {
    std::for_each(refFac->begin(),refFac->end(), [this,&K_TH,&L_TH](Facette::Fac & fac)
        {
        //gmm::dense_matrix <double> Ks(3*Facette::N,3*Facette::N);
        //std::vector <double> Ls(3*Facette::N);
        double Ks[3*Facette::N][3*Facette::N] = { {0} };
        double Ls[3*Facette::N] = {0};
        
            fac.integrales(settings.paramFacette,Ls);     
        fac.projection( Ks, Ls);
        fac.treated = false;
        if(my_mutex.try_lock())
            {
            fac.assemblage_mat(K_TH);
            fac.assemblage_vect(L_TH);
            fac.treated =true;
            my_mutex.unlock();    
            }
        else { fac.treated = false; }
        });
    }
);//end thread

for(int i=0;i<(NbTH+1);i++) {tab_TH[i].join();}


insertCoeff<Tetra::Tet>(*refTet,K_TH,L_TH);
insertCoeff<Facette::Fac>(*refFac,K_TH,L_TH);

counter.tac();

if(settings.verbose) { std::cout << "..matrix assembly done in " << counter.elapsed() << "s" << std::endl; }

read_matrix Kr(2*NOD,2*NOD);
gmm::copy(K_TH, Kr);

std::vector<double> Lr(2*NOD);
for(int i=0;i<2*NOD;i++) {Lr[i] = L_TH[i];} // should be replaced by a std::copy 

std::vector<double> Xw(2*NOD);

gmm::iteration bicg_iter(1e-6);
gmm::iteration gmr_iter(1e-6);
bicg_iter.set_maxiter(settings.MAXITER);
gmr_iter.set_maxiter(settings.MAXITER);
bicg_iter.set_noisy(false);//VERBOSE
gmr_iter.set_noisy(false);//VERBOSE

counter.tic();

if (!nt) 
    {
    if(settings.verbose) { std::cout << "computing prc.";std::fflush(NULL); }
    prc = new gmm::diagonal_precond <read_matrix> (Kr);
	counter.tac();    
	if(settings.verbose) { std::cout << "elapsed time = " << counter.elapsed() << "s" << std::endl; }
    }
else if (!(nt % (settings.REFRESH_PRC)))
    {
    delete prc;
    if(settings.verbose) { std::cout << "computing prc.";std::fflush(NULL); }
    prc = new gmm::diagonal_precond <read_matrix> (Kr);
    counter.tac();    
	if(settings.verbose) { std::cout << "elapsed time = " << counter.elapsed() << "s" << std::endl; }
    } 

counter.tic();
gmm::bicgstab(Kr, Xw, Lr, *prc, bicg_iter);

if (!(bicg_iter.converged() )) 
    {
    counter.tac();
	if(settings.verbose) { std::cout << "bicg FAILED in " << bicg_iter.get_iteration() << "iterations, duration: " << counter.elapsed() << " s" << std::endl; }
    gmm::diagonal_precond <read_matrix>  gmr_prc (Kr);
    counter.tic();
    
    gmm::clear(Xw);
    gmm::gmres(Kr, Xw, Lr, gmr_prc, 50, gmr_iter);

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

//read_vector Xr(2*NOD);    gmm::copy(Xw, Xr);
double v2max = 0.0;

int i=0;
std::for_each(refNode->begin(),refNode->end(),
    [this,&t_prm,&i,&v2max,&Xw](Nodes::Node &n)
        {
        double vp = Xw[i];
        double vq = Xw[NOD+i];
        double v2 = vp*vp + vq*vq;
        if (v2>v2max) { v2max = v2; }
        n.make_evol(vp,vq,t_prm.dt);    
        i++;    
        }
);//end for_each

v_max = sqrt(v2max);
return 0;
}
