#include <algorithm>


#include <functional>

#include "Utils/FTic.hpp"

#include "linear_algebra.h"


void LinAlgebra::assemblageTet(const int i_buf, const int i_tet,
           gmm::dense_matrix <double> const& Ke, std::vector <double> const& Le,
           write_matrix &K,write_vector &L)
{
    if(my_lock->try_lock())
        {
        refTet[i_buf][i_tet].assemblage(NOD,Ke,Le,K,L);
        my_lock->unlock();    
        return;
        }
    else { buff_tet[i_buf].push(Tetra::Obj(i_tet,Ke,Le)); }
}

void LinAlgebra::assemblageFac(const int i_fac,
           gmm::dense_matrix <double> const& Ke, std::vector <double> const& Le,
            write_matrix &K,write_vector &L)
{
    if(my_lock->try_lock())
        {
        (*refFac)[i_fac].assemblage(NOD,Ke,Le,K,L);
        my_lock->unlock();    
        return;
        }
    else
        { buff_fac.push(Facette::Obj(i_fac,Ke,Le)); }
}

    
int LinAlgebra::vsolve(long nt)
{
FTic counter;

counter.tic();

base_projection();   // definit plan tangent

//mtl::compressed2D<double> K_TH(2*NOD, 2*NOD);
//mtl::mat::set_to_zero( K_TH );

//ins = new sparseInserter(K_TH,64);

//mtl::dense_vector<double> L_TH(2*NOD);
//mtl::vec::set_to_zero( L_TH );

write_matrix K_TH(2*NOD, 2*NOD);
write_vector L_TH(2*NOD);

std::cout << "solver #" << nt <<  "\nstarting matrix assembly..." << std::endl;

for(int i=0;i<NbTH;i++) 
    {
    tab_TH[i] = std::thread( [this,&K_TH,&L_TH,i]() 
        {
            thread_local gmm::dense_matrix <double> K(3*Tetra::N,3*Tetra::N);
            thread_local gmm::dense_matrix <double> Kp(2*Tetra::N,2*Tetra::N);
            thread_local std::vector <double> L(3*Tetra::N);
            thread_local std::vector <double>  Lp(2*Tetra::N);
            int i_tet =0;
            std::for_each(refTet[i].begin(),refTet[i].end(), [this,&K_TH,&L_TH,i,&i_tet](Tetra::Tet const& tet)
                {
                //mtl::mat::set_to_zero(K); mtl::mat::set_to_zero(Kp);
                //mtl::vec::set_to_zero(L); mtl::vec::set_to_zero(Lp);
                tet.integrales(settings->second_order,settings->paramTetra,Hext,DW_vz,settings->theta,dt,settings->TAUR,K, L);     
                tet.projection( K, L, Kp, Lp);
                assemblageTet(i, i_tet, Kp, Lp, K_TH, L_TH );
                /*
                if(my_lock->try_lock())
                    {
                    refTet[i][i_tet].assemblage(NOD, Kp, Lp, K_TH, L_TH);
                    my_lock->unlock();    
                    return;
                    }
                else { buff_tet[i].push(Tetra::Obj(i_tet,Kp,Lp)); }
                */
                i_tet++;
                });//end for_each
        }); //end thread
    }
    
tab_TH[NbTH] = std::thread( [this,&K_TH,&L_TH]()
    {
    thread_local gmm::dense_matrix <double> Ks(3*Facette::N,3*Facette::N);
    thread_local gmm::dense_matrix <double> Ksp(2*Facette::N,2*Facette::N);
    thread_local std::vector <double> Ls(3*Facette::N);
    thread_local std::vector <double> Lsp(2*Facette::N);
    int i_fac=0;
    std::for_each(refFac->begin(),refFac->end(),
    [this,&K_TH,&L_TH,&i_fac](Facette::Fac const& fac)
        {
        //mtl::mat::set_to_zero(Ks); mtl::mat::set_to_zero(Ksp);
        //mtl::vec::set_to_zero(Ls); mtl::vec::set_to_zero(Lsp);
        fac.integrales(settings->paramFacette, Ls);     
        fac.projection( Ks, Ls, Ksp, Lsp);
        assemblageFac(i_fac, Ksp, Lsp,K_TH, L_TH);    
        i_fac++;
            
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
        refTet[i][x.idx].assemblage(NOD,x.Ke,x.Le,K_TH,L_TH);
        buff_tet[i].pop();
    }
}

while (!buff_fac.empty())
    {
        Facette::Obj const& x = buff_fac.front();
        (*refFac)[x.idx].assemblage(NOD,x.Ke,x.Le,K_TH,L_TH);
        buff_fac.pop();
    }

//delete ins;//K_TH should be ready

counter.tac();

std::cout << "ending matrix assembly." << std::endl;
if(VERBOSE) { std::cout << "elapsed time = " << counter.elapsed() << "s" << std::endl; }

/*
mtl::dense_vector<double> Xw(2*NOD);
mtl::vec::set_to_zero(Xw);//should be useless
itl::noisy_iteration<double> bicg_iter(L_TH,settings->MAXITER,1e-6);
bicg_iter.set_quite(true);
*/

read_matrix Kr(2*NOD,2*NOD);    gmm::copy(K_TH, Kr);
read_vector Lr(2*NOD);          gmm::copy(L_TH, Lr);
write_vector Xw(2*NOD);

gmm::iteration bicg_iter(1e-6);
gmm::iteration gmr_iter(1e-6);
bicg_iter.set_maxiter(settings->MAXITER);
gmr_iter.set_maxiter(settings->MAXITER);
bicg_iter.set_noisy(VERBOSE);
gmr_iter.set_noisy(VERBOSE);

counter.tic();

if (!nt) 
    {
    if(VERBOSE) { std::cout << "computing prc.";std::fflush(NULL); }
    //prc = new itl::pc::diagonal < mtl::compressed2D<double> >(K_TH); //mind the constructor call syntax
	prc = new gmm::diagonal_precond <read_matrix> (Kr);
	counter.tac();    
	if(VERBOSE) { std::cout << "elapsed time = " << counter.elapsed() << "s" << std::endl; }
    }
else if (!(nt % (settings->REFRESH_PRC)))
    {
    delete prc;
    if(VERBOSE) { std::cout << "computing prc.";std::fflush(NULL); }
    //prc = new itl::pc::diagonal < mtl::compressed2D<double> >(K_TH); //mind the constructor call syntax
    prc = new gmm::diagonal_precond <read_matrix> (Kr);
    counter.tac();    
	if(VERBOSE) { std::cout << "elapsed time = " << counter.elapsed() << "s" << std::endl; }
    } 


counter.tic();
gmm::bicgstab(Kr, Xw, L_TH, *prc, bicg_iter);

//itl::noisy_iteration<double> gmr_iter(L_TH,settings->MAXITER,1e-6);
//gmr_iter.set_quite(true);

if (!(bicg_iter.converged() )) 
    {
    counter.tac();
	if(VERBOSE) { std::cout << "bicg FAILED in " << bicg_iter.get_iteration() << "iterations, duration: " << counter.elapsed() << " s" << std::endl; }
    //itl::pc::diagonal < mtl::compressed2D<double> >  gmr_prc (K_TH);
    gmm::diagonal_precond <read_matrix>  gmr_prc (Kr);
    counter.tic();
    
    //itl::gmres(K_TH, Xw, L_TH, gmr_prc, gmr_prc,gmr_iter, gmr_iter);
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
    if(VERBOSE) { std::cout << "v-solve converged in " << bicg_iter.get_iteration() << " (bicg,prc-" << (nt % (settings->REFRESH_PRC)) << ") :duration: " 
<< counter.elapsed() << " s" << std::endl; }
    }

read_vector Xr(2*NOD);    gmm::copy(Xw, Xr);
double v2max = 0.0;

int i=0;
std::for_each(refNode->begin(),refNode->end(),
    [this,&i,&v2max,&Xr](Node &n)
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
