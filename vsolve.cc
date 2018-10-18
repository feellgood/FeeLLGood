#include <algorithm>

#include <thread>
#include <functional>

#include "Utils/FTic.hpp"

#include "linear_algebra.h"

void LinAlgebra::feedMat(const int NOD,double dt, mtl::compressed2D<double> *K_T, mtl::dense_vector<double> *L_T, std::vector<Tetra::Tet>::iterator it_b, std::vector<Tetra::Tet>::iterator it_e)
{
mtl::dense2D <double> K(3*Tetra::N,3*Tetra::N), Kp(2*Tetra::N,2*Tetra::N);
mtl::dense_vector <double> L(3*Tetra::N), Lp(2*Tetra::N);

//mtl::dense2D <double> P(2*Tetra::N,3*Tetra::N);

sparseInserter *ins;
ins = new sparseInserter(*(K_T),256); //64 maybe too small

for_each(it_b,it_e,
    [this,dt,L_T,&K,&L,&Kp,&Lp,&ins,NOD](Tetra::Tet & tet) //,&P
        {
        mtl::mat::set_to_zero(K); mtl::mat::set_to_zero(Kp);
        mtl::vec::set_to_zero(L); mtl::vec::set_to_zero(Lp);
        tet.integrales(settings->paramTetra,Hext,DW_vz,settings->theta,dt,settings->TAUR,K, L);     
        tet.projection( K, L, Kp, Lp);//P
        tet.assemblage( ins,NOD, Kp, Lp,*(L_T) );// on passe l'inserter plutot que Kw
        }
);
delete ins;    
}

int LinAlgebra::vsolve(double dt,long nt)
{
FTic counter;

const int NOD = refNode->size();

counter.tic();

base_projection();   // definit plan tangent

const int NbTH=8;

std::vector < mtl::compressed2D<double> * > K_TH;
std::vector < mtl::dense_vector<double> * > L_TH;
std::thread tab_TH[NbTH];

const unsigned long block_size = std::distance(refTet->begin(),refTet->end())/NbTH;


for(int i=0;i<NbTH;i++)
    {
    K_TH.push_back( new mtl::compressed2D<double>(2*NOD, 2*NOD) );
    mtl::mat::set_to_zero( *(K_TH[i]) );
    
    L_TH.push_back( new mtl::dense_vector<double>(2*NOD)  );
    mtl::vec::set_to_zero( *(L_TH[i]) );
    } 

std::vector<Tetra::Tet>::iterator it_begin = refTet->begin();

for(int i=0;i<(NbTH-1);i++) 
    {
    std::vector<Tetra::Tet>::iterator it_end = it_begin;
    std::advance(it_end,block_size);
    mtl::compressed2D<double> *K = K_TH[i];
    mtl::dense_vector<double> *L = L_TH[i];

    tab_TH[i] = std::thread( [this,NOD,dt,K,L,it_begin,it_end]() {feedMat(NOD,dt,K,L,it_begin,it_end);} ); 
    it_begin = it_end;
    }

mtl::compressed2D<double> *K = K_TH[NbTH-1];
mtl::dense_vector<double> *L = L_TH[NbTH-1];
tab_TH[NbTH-1] = std::thread( [this,NOD,dt,K,L,it_begin]() {feedMat(NOD,dt,K,L,it_begin,refTet->end());} );
    
for(int i=0;i<NbTH;i++) {tab_TH[i].join();}

mtl::compressed2D<double>  Kw(2*NOD, 2*NOD);
mtl::mat::set_to_zero(Kw);
mtl::dense_vector<double> Lw(2*NOD);
mtl::vec::set_to_zero(Lw);

sparseInserter *ins;
ins = new sparseInserter(Kw,256); //64 maybe too small

mtl::dense2D <double> Ps(2*Facette::N,3*Facette::N);

mtl::dense2D <double> Ks(3*Facette::N,3*Facette::N), Ksp(2*Facette::N,2*Facette::N);
mtl::dense_vector <double> Ls(3*Facette::N), Lsp(2*Facette::N);

for_each(refFac->begin(),refFac->end(),
    [this,&Lw,&Ks,&Ls,&Ksp,&Lsp,&ins,NOD,&Ps](Facette::Fac & fac)
        {
        mtl::mat::set_to_zero(Ks); mtl::mat::set_to_zero(Ksp);
        mtl::vec::set_to_zero(Ls); mtl::vec::set_to_zero(Lsp);
        fac.integrales(settings->paramFacette, Ls);     
        fac.projection(Ps, Ks, Ls, Ksp, Lsp);
        fac.assemblage(ins, NOD, Ksp, Lsp, Lw);    
        }
);

delete ins;//Kw should be ready

for(int n=0;n<NbTH;n++)
    {
    Kw += *(K_TH[n]);
    delete K_TH[n];
    Lw += *(L_TH[n]);
    delete L_TH[n];
    }
    
counter.tac();

if(VERBOSE) { std::cout << "elapsed time = " << counter.elapsed() << "s" << std::endl; }

mtl::dense_vector<double> Xw(2*NOD);
mtl::vec::set_to_zero(Xw);//should be useless

itl::noisy_iteration<double> bicg_iter(Lw,settings->MAXITER,1e-6);

bicg_iter.set_quite(true);

counter.tic();

if (!nt) 
    {
    if(VERBOSE) { std::cout << "computing prc.";std::fflush(NULL); }
    prc = new itl::pc::diagonal < mtl::compressed2D<double> >(Kw); //mind the constructor call syntax
	counter.tac();    
	if(VERBOSE) { std::cout << "elapsed time = " << counter.elapsed() << "s" << std::endl; }
    }
else if (!(nt % (settings->REFRESH_PRC)))
    {
    delete prc;
    if(VERBOSE) { std::cout << "computing prc.";std::fflush(NULL); }
    prc = new itl::pc::diagonal < mtl::compressed2D<double> >(Kw); //mind the constructor call syntax
    counter.tac();    
	if(VERBOSE) { std::cout << "elapsed time = " << counter.elapsed() << "s" << std::endl; }
    } 


counter.tic();

bicgstab(Kw, Xw, Lw, *prc, bicg_iter);

itl::noisy_iteration<double> gmr_iter(Lw,settings->MAXITER,1e-6);
gmr_iter.set_quite(true);

if (!(bicg_iter.is_converged() )) 
    {
    counter.tac();
	if(VERBOSE) { std::cout << "bicg FAILED in " << bicg_iter.iterations() << "iterations, duration: " << counter.elapsed() << " s" << std::endl; }
    itl::pc::diagonal < mtl::compressed2D<double> >  gmr_prc (Kw);
    counter.tic();
    
    itl::gmres(Kw, Xw, Lw, gmr_prc, gmr_prc,gmr_iter, gmr_iter);
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
