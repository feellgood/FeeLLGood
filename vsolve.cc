#include <algorithm>


#include <functional>

#include "Utils/FTic.hpp"

#include "linear_algebra.h"

    
int LinAlgebra::vsolve(long nt)
{
FTic counter;

counter.tic();

base_projection();   // definit plan tangent

write_matrix K_TH(2*NOD, 2*NOD);
write_vector L_TH(2*NOD);

std::cout << "solver iteration nt#" << nt <<  ", starting matrix assembly..." << std::endl;

for(int i=0;i<NbTH;i++) 
    {
    tab_TH[i] = std::thread( [this,&K_TH,&L_TH,i]() 
        {
            thread_local gmm::dense_matrix <double> K(3*Tetra::N,3*Tetra::N);
            thread_local std::vector <double> L(3*Tetra::N);
            
            std::for_each(refTet[i].begin(),refTet[i].end(), [this,&K_TH,&L_TH](Tetra::Tet & tet)
                {
                tet.integrales(settings.paramTetra,Hext,DW_vz,settings.theta,dt,settings.TAUR,K, L);     
                tet.projection( K, L);
                if(my_lock->try_lock())
                    {
                    tet.assemblage(K_TH, L_TH);
                    tet.treated = true;
                    my_lock->unlock();    
                    return;
                    }
                else { tet.treated = false;}
                });//end for_each
        }); //end thread
    }
    
tab_TH[NbTH] = std::thread( [this,&K_TH,&L_TH]()
    {
    thread_local gmm::dense_matrix <double> Ks(3*Facette::N,3*Facette::N);
    //thread_local gmm::dense_matrix <double> Ksp(2*Facette::N,2*Facette::N);
    thread_local std::vector <double> Ls(3*Facette::N);
    //thread_local std::vector <double> Lsp(2*Facette::N);
    
    std::for_each(refFac->begin(),refFac->end(), [this,&K_TH,&L_TH](Facette::Fac & fac)
        {
        fac.integrales(settings.paramFacette, Ls);     
        fac.projection( Ks, Ls);
        
        if(my_lock->try_lock())
            {
            fac.assemblage(K_TH,L_TH);
            fac.treated =true;
            my_lock->unlock();    
            return;
            }
        else { fac.treated = false; }
        }
        );
    }
);//end thread

for(int i=0;i<(NbTH+1);i++) {tab_TH[i].join();}

for(int i=0;i<(NbTH);i++)
    { std::for_each(refTet[i].begin(),refTet[i].end(),[&K_TH,&L_TH](Tetra::Tet const& tet){if(!tet.treated) tet.assemblage(K_TH,L_TH);}); }

std::for_each( (*refFac).begin(), (*refFac).end(), [&K_TH,&L_TH](Facette::Fac const& fac){if(!fac.treated) fac.assemblage(K_TH,L_TH); } );    

/*
while (!buff_fac.empty())
    {
        Facette::Obj const& x = buff_fac.front();
        (*refFac)[x.idx].assemblage(NOD,x.Ke,x.Le,K_TH,L_TH);
        buff_fac.pop();
    }
*/
counter.tac();

if(VERBOSE) { std::cout << "Matrix assembly done, elapsed time = " << counter.elapsed() << "s" << std::endl; }

read_matrix Kr(2*NOD,2*NOD);    gmm::copy(K_TH, Kr);
read_vector Lr(2*NOD);          gmm::copy(L_TH, Lr);
write_vector Xw(2*NOD);

gmm::iteration bicg_iter(1e-6);
gmm::iteration gmr_iter(1e-6);
bicg_iter.set_maxiter(settings.MAXITER);
gmr_iter.set_maxiter(settings.MAXITER);
bicg_iter.set_noisy(false);//VERBOSE
gmr_iter.set_noisy(false);//VERBOSE

counter.tic();

if (!nt) 
    {
    if(VERBOSE) { std::cout << "computing prc.";std::fflush(NULL); }
    prc = new gmm::diagonal_precond <read_matrix> (Kr);
	counter.tac();    
	if(VERBOSE) { std::cout << "elapsed time = " << counter.elapsed() << "s" << std::endl; }
    }
else if (!(nt % (settings.REFRESH_PRC)))
    {
    delete prc;
    if(VERBOSE) { std::cout << "computing prc.";std::fflush(NULL); }
    prc = new gmm::diagonal_precond <read_matrix> (Kr);
    counter.tac();    
	if(VERBOSE) { std::cout << "elapsed time = " << counter.elapsed() << "s" << std::endl; }
    } 


counter.tic();
gmm::bicgstab(Kr, Xw, L_TH, *prc, bicg_iter);

if (!(bicg_iter.converged() )) 
    {
    counter.tac();
	if(VERBOSE) { std::cout << "bicg FAILED in " << bicg_iter.get_iteration() << "iterations, duration: " << counter.elapsed() << " s" << std::endl; }
    gmm::diagonal_precond <read_matrix>  gmr_prc (Kr);
    counter.tic();
    
    gmm::clear(Xw);
    gmm::gmres(Kr, Xw, Lr, gmr_prc, 50, gmr_iter);

    if (!(gmr_iter.converged() )) {
	counter.tac();
    if(VERBOSE) { std::cout << "gmres FAILED in " << gmr_iter.get_iteration() << "iterations, duration: " << counter.elapsed() << " s" << std::endl; }
    return 1;
    }
    else {counter.tac();
        if(VERBOSE) { std::cout << "v-solve in " << gmr_iter.get_iteration() << " iterations (gmres) duration: " << counter.elapsed() << " s" << std::endl; }
        }
    }
else 
    {counter.tac();
    if(VERBOSE) { std::cout << "v-solve converged in " << bicg_iter.get_iteration() << " (bicg,prc-" << (nt % (settings.REFRESH_PRC)) << ") :duration: " 
<< counter.elapsed() << " s" << std::endl; }
    }

read_vector Xr(2*NOD);    gmm::copy(Xw, Xr);
double v2max = 0.0;

int i=0;
std::for_each(refNode->begin(),refNode->end(),
    [this,&i,&v2max,&Xr](Nodes::Node &n)
        {
        double vp = Xr[i];
        double vq = Xr[NOD+i];
        double v2 = vp*vp + vq*vq;
        if (v2>v2max) { v2max = v2; }
        n.make_evol(vp,vq,dt);    
        i++;    
        }
);//end for_each

v_max = sqrt(v2max);
return 0;
}
