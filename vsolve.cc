#include <algorithm>


#include <functional>

#include "Utils/FTic.hpp"

#include "linear_algebra.h"


void LinAlgebra::assemblage(const int NOD,const int N,const int ind[],
           mtl::dense2D <double> const& Ke, mtl::dense_vector <double> const& Le,
           mtl::dense_vector<double> &L)//mtl::compressed2D<double> &K, avant dernier
    {
    for (int i=0; i < N; i++)
        {
        int i_= ind[i];             
        my_lock->lock();
        for (int j=0; j < N; j++)
            {
            int j_= ind[j];
            (*ins)(NOD+i_,j_) << Ke(i,j);      (*ins)(NOD+i_, NOD+j_) << Ke(  i,N+j);
            (*ins)(    i_,j_) << Ke(N+i,j);    (*ins)(    i_, NOD+j_) << Ke(N+i,N+j);
            }
        L(NOD+i_) += Le(i);//L[NOD+i_]+= Le[  i];
        L(i_) += Le(N+i);//L[    i_]+= Le[N+i];
        
        my_lock->unlock();   
        }
    }


int LinAlgebra::vsolve(double dt,long nt)
{
FTic counter;

const int NOD = refNode->size();

counter.tic();

base_projection();   // definit plan tangent

mtl::compressed2D<double> K_TH(2*NOD, 2*NOD);
mtl::mat::set_to_zero( K_TH );

ins = new sparseInserter(K_TH,64);

mtl::dense_vector<double> L_TH(2*NOD);
mtl::vec::set_to_zero( L_TH );

const unsigned long block_size = std::distance(refTet->begin(),refTet->end())/NbTH;

std::vector<Tetra::Tet>::iterator it_begin = refTet->begin();

for(int i=0;i<(NbTH-1);i++) 
    {
    std::vector<Tetra::Tet>::iterator it_end = it_begin;
    std::advance(it_end,block_size);
    
    tab_TH[i] = std::thread( [this,NOD,dt,&L_TH,it_begin,it_end]() 
        {
            mtl::dense2D <double> K(3*Tetra::N,3*Tetra::N), Kp(2*Tetra::N,2*Tetra::N);
            mtl::dense_vector <double> L(3*Tetra::N), Lp(2*Tetra::N);
            
            std::for_each(it_begin,it_end, [this,dt,&L_TH,&K,&L,&Kp,&Lp,NOD](Tetra::Tet & tet) //,&P
                {
                mtl::mat::set_to_zero(K); mtl::mat::set_to_zero(Kp);
                mtl::vec::set_to_zero(L); mtl::vec::set_to_zero(Lp);
                tet.integrales(settings->paramTetra,Hext,DW_vz,settings->theta,dt,settings->TAUR,K, L);     
                tet.projection( K, L, Kp, Lp);//P
                assemblage(NOD, Tetra::N, tet.ind,Kp,Lp, L_TH );
                }
                );//end for_each
        } 
    ); //end thread
    
    it_begin = it_end;
    }

    //last thread must be treated differently 
tab_TH[NbTH-1] = std::thread( [this,NOD,dt,&L_TH,it_begin]() 
    {
        mtl::dense2D <double> K(3*Tetra::N,3*Tetra::N), Kp(2*Tetra::N,2*Tetra::N);
        mtl::dense_vector <double> L(3*Tetra::N), Lp(2*Tetra::N);
            
        std::for_each(it_begin,refTet->end(), [this,dt,&L_TH,&K,&L,&Kp,&Lp,NOD](Tetra::Tet & tet) //,&P
            {
            mtl::mat::set_to_zero(K); mtl::mat::set_to_zero(Kp);
            mtl::vec::set_to_zero(L); mtl::vec::set_to_zero(Lp);
            tet.integrales(settings->paramTetra,Hext,DW_vz,settings->theta,dt,settings->TAUR,K, L);     
            tet.projection( K, L, Kp, Lp);//P
            assemblage(NOD, Tetra::N, tet.ind,Kp,Lp, L_TH );
            }
            );//end for_each
    } 
);//end last thread
    
for(int i=0;i<NbTH;i++) {tab_TH[i].join();}

mtl::dense2D <double> Ps(2*Facette::N,3*Facette::N);

mtl::dense2D <double> Ks(3*Facette::N,3*Facette::N), Ksp(2*Facette::N,2*Facette::N);
mtl::dense_vector <double> Ls(3*Facette::N), Lsp(2*Facette::N);

for_each(refFac->begin(),refFac->end(),
    [this,&L_TH,&Ks,&Ls,&Ksp,&Lsp,NOD,&Ps](Facette::Fac & fac)
        {
        mtl::mat::set_to_zero(Ks); mtl::mat::set_to_zero(Ksp);
        mtl::vec::set_to_zero(Ls); mtl::vec::set_to_zero(Lsp);
        fac.integrales(settings->paramFacette, Ls);     
        fac.projection(Ps, Ks, Ls, Ksp, Lsp);
        fac.assemblage(ins,NOD, Ksp, Lsp, L_TH);    
        }
);

delete ins;//Kw should be ready

counter.tac();

if(VERBOSE) { std::cout << "elapsed time = " << counter.elapsed() << "s" << std::endl; }

mtl::dense_vector<double> Xw(2*NOD);
mtl::vec::set_to_zero(Xw);//should be useless

itl::noisy_iteration<double> bicg_iter(L_TH,settings->MAXITER,1e-6);

bicg_iter.set_quite(true);

counter.tic();

if (!nt) 
    {
    if(VERBOSE) { std::cout << "computing prc.";std::fflush(NULL); }
    prc = new itl::pc::diagonal < mtl::compressed2D<double> >(K_TH); //mind the constructor call syntax
	counter.tac();    
	if(VERBOSE) { std::cout << "elapsed time = " << counter.elapsed() << "s" << std::endl; }
    }
else if (!(nt % (settings->REFRESH_PRC)))
    {
    delete prc;
    if(VERBOSE) { std::cout << "computing prc.";std::fflush(NULL); }
    prc = new itl::pc::diagonal < mtl::compressed2D<double> >(K_TH); //mind the constructor call syntax
    counter.tac();    
	if(VERBOSE) { std::cout << "elapsed time = " << counter.elapsed() << "s" << std::endl; }
    } 


counter.tic();

bicgstab(K_TH, Xw, L_TH, *prc, bicg_iter);

itl::noisy_iteration<double> gmr_iter(L_TH,settings->MAXITER,1e-6);
gmr_iter.set_quite(true);

if (!(bicg_iter.is_converged() )) 
    {
    counter.tac();
	if(VERBOSE) { std::cout << "bicg FAILED in " << bicg_iter.iterations() << "iterations, duration: " << counter.elapsed() << " s" << std::endl; }
    itl::pc::diagonal < mtl::compressed2D<double> >  gmr_prc (K_TH);
    counter.tic();
    
    itl::gmres(K_TH, Xw, L_TH, gmr_prc, gmr_prc,gmr_iter, gmr_iter);
    if (!(gmr_iter.is_converged() )) {
	counter.tac();
    if(VERBOSE) { std::cout << "gmres FAILED in " << gmr_iter.iterations() << "iterations, duration: " << counter.elapsed() << " s" << std::endl; }
    return 1;
    }
    else {counter.tac();
        if(VERBOSE) { std::cout << "v-solve in " << gmr_iter.iterations() << " iterations (gmres) duration: " << counter.elapsed() << " s" << std::endl; }
        }
    }
else 
    {counter.tac();
    if(VERBOSE) { std::cout << "v-solve converged in " << bicg_iter.iterations() << " (bicg,prc-" << (nt % (settings->REFRESH_PRC)) << ") :duration: " 
<< counter.elapsed() << " s" << std::endl; }
    }

double v2max = 0.0;

int i=0;
for_each(refNode->begin(),refNode->end(),
    [&i,&v2max,&Xw,NOD,dt](Node &n)
        {
        double vp = Xw[i];
        double vq = Xw[NOD+i];
        double v2 = vp*vp + vq*vq;
        if (v2>v2max) { v2max = v2; }
        n.make_evol(vp,vq,dt);    
        i++;    
        }
);//end for_each

v_max = sqrt(v2max);
return 0;
}
