#include <algorithm>


#include <functional>

#include "Utils/FTic.hpp"

#include "linear_algebra.h"

void LinAlgebra::assemblage_monoThread(const int N, const int ind[],
           mtl::dense2D <double> const& Ke, mtl::dense_vector <double> const& Le,
           mtl::dense_vector<double> &L)
{
    for (int i=0; i < N; i++)
        {
        int i_= ind[i];             
            for (int j=0; j < N; j++)
                {
                int j_= ind[j];
                (*ins)(NOD+i_,j_) << Ke(i,j);      (*ins)(NOD+i_, NOD+j_) << Ke(  i,N+j);
                (*ins)(    i_,j_) << Ke(N+i,j);    (*ins)(    i_, NOD+j_) << Ke(N+i,N+j);
                }
            L(NOD+i_) += Le(i);//L[NOD+i_]+= Le[  i];
            L(i_) += Le(N+i);//L[    i_]+= Le[N+i];
        }
}

void LinAlgebra::assemblageTet(const int i, const int ind[],
           mtl::dense2D <double> const& Ke, mtl::dense_vector <double> const& Le,
           mtl::dense_vector<double> &L)
{
    if(my_lock->try_lock())
        {
        for (int i=0; i < Tetra::N; i++)
            {
            int i_= ind[i];             
            
                for (int j=0; j < Tetra::N; j++)
                    {
                    int j_= ind[j];
                    (*ins)(NOD+i_,j_) << Ke(i,j);      (*ins)(NOD+i_, NOD+j_) << Ke(  i,Tetra::N+j);
                    (*ins)(    i_,j_) << Ke(Tetra::N+i,j);    (*ins)(    i_, NOD+j_) << Ke(Tetra::N+i,Tetra::N+j);
                    }
                L(NOD+i_) += Le(i);//L[NOD+i_]+= Le[  i];
                L(i_) += Le(Tetra::N+i);//L[    i_]+= Le[N+i];
            }
        my_lock->unlock();    
        return;
        }
    else
        {
        buff_tet[i].push(Tetra::Obj(ind,Ke,Le));    
        }
}

void LinAlgebra::assemblageFac(const int ind[],
           mtl::dense2D <double> const& Ke, mtl::dense_vector <double> const& Le,
           mtl::dense_vector<double> &L)
{
    if(my_lock->try_lock())
        {
        for (int i=0; i < Facette::N; i++)
            {
            int i_= ind[i];             
            
                for (int j=0; j < Facette::N; j++)
                    {
                    int j_= ind[j];
                    (*ins)(NOD+i_,j_) << Ke(i,j);      (*ins)(NOD+i_, NOD+j_) << Ke(  i,Facette::N+j);
                    (*ins)(    i_,j_) << Ke(Facette::N+i,j);    (*ins)(    i_, NOD+j_) << Ke(Facette::N+i,Facette::N+j);
                    }
                L(NOD+i_) += Le(i);//L[NOD+i_]+= Le[  i];
                L(i_) += Le(Facette::N+i);//L[    i_]+= Le[N+i];
            }
        my_lock->unlock();    
        return;
        }
    else
        {
        buff_fac.push(Facette::Obj(ind,Ke,Le));    
        }
}

    
int LinAlgebra::vsolve(long nt)
{
FTic counter;

counter.tic();

base_projection();   // definit plan tangent

mtl::compressed2D<double> K_TH(2*NOD, 2*NOD);
mtl::mat::set_to_zero( K_TH );

ins = new sparseInserter(K_TH,64);

mtl::dense_vector<double> L_TH(2*NOD);
mtl::vec::set_to_zero( L_TH );

for(int i=0;i<NbTH;i++) 
    {
    tab_TH[i] = std::thread( [this,&L_TH,i]() 
        {
            thread_local mtl::dense2D <double> K(3*Tetra::N,3*Tetra::N);
            thread_local mtl::dense2D <double> Kp(2*Tetra::N,2*Tetra::N);
            thread_local mtl::dense_vector <double> L(3*Tetra::N);
            thread_local mtl::dense_vector <double>  Lp(2*Tetra::N);
            
            std::for_each(refTet[i].begin(),refTet[i].end(), [this,&L_TH,i](Tetra::Tet const& tet)
                {
                mtl::mat::set_to_zero(K); mtl::mat::set_to_zero(Kp);
                mtl::vec::set_to_zero(L); mtl::vec::set_to_zero(Lp);
                tet.integrales(settings->paramTetra,Hext,DW_vz,settings->theta,dt,settings->TAUR,K, L);     
                tet.projection( K, L, Kp, Lp);
                assemblageTet(i,tet.ind,Kp ,Lp, L_TH );
                });//end for_each
        }); //end thread
    }
    
tab_TH[NbTH] = std::thread( [this,&L_TH]()
    {
    //thread_local mtl::dense2D <double> Ps(2*Facette::N,3*Facette::N);
    thread_local mtl::dense2D <double> Ks(3*Facette::N,3*Facette::N);
    thread_local mtl::dense2D <double> Ksp(2*Facette::N,2*Facette::N);
    thread_local mtl::dense_vector <double> Ls(3*Facette::N);
    thread_local mtl::dense_vector <double> Lsp(2*Facette::N);

    std::for_each(refFac->begin(),refFac->end(),
    [this,&L_TH](Facette::Fac const& fac)
        {
        mtl::mat::set_to_zero(Ks); mtl::mat::set_to_zero(Ksp);
        mtl::vec::set_to_zero(Ls); mtl::vec::set_to_zero(Lsp);
        fac.integrales(settings->paramFacette, Ls);     
        fac.projection( Ks, Ls, Ksp, Lsp);//(Ps, Ks, Ls, Ksp, Lsp);
        assemblageFac(fac.ind, Ksp, Lsp, L_TH);    
        }
        );
    }
);//end thread

for(int i=0;i<(NbTH+1);i++) {tab_TH[i].join();}

for(int i=0;i<(NbTH);i++)
{
    while (!buff_tet[i].empty())
    {
        Tetra::Obj const& x = buff_tet[i].front();
        assemblage_monoThread(Tetra::N,x.ind,x.Ke,x.Le,L_TH);
        //delete [] x.ind; 
        buff_tet[i].pop();
    }
/*
    std::for_each(buff_TH[i].begin(),buff_TH[i].end(), [this,&L_TH](Obj const& x) 
    { 
    assemblage_monoThread(x.N,x.ind,x.Ke,x.Le,L_TH);    
    delete [] x.ind; 
    } );
*/
//buff_TH[i].clear();
    
}

while (!buff_fac.empty())
    {
        Facette::Obj const& x = buff_fac.front();
        assemblage_monoThread(Facette::N,x.ind,x.Ke,x.Le,L_TH);
        //delete [] x.ind; 
        buff_fac.pop();
    }

delete ins;//K_TH should be ready

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
std::for_each(refNode->begin(),refNode->end(),
    [this,&i,&v2max,&Xw](Node &n)
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
