#include <algorithm>

#include "Utils/FTic.hpp"

#include "linear_algebra.h"


int LinAlgebra::vsolve(double dt,long nt)
{
FTic counter;

const int NOD = refNode->size();

counter.tic();

base_projection();   // definit plan tangent

mtl::compressed2D<double>  Kw(2*NOD, 2*NOD);
mtl::mat::set_to_zero(Kw);
mtl::dense_vector<double> Lw(2*NOD);
mtl::vec::set_to_zero(Lw);

if(VERBOSE) { std::cout <<"DW speed = " << DW_vz << "\nmatrix size " << 2*NOD << "; assembling ..." << std::endl; }

sparseInserter *ins;
ins = new sparseInserter(Kw,256); //64 maybe too small

mtl::dense2D <double> K(3*Tetra::N,3*Tetra::N), Kp(2*Tetra::N,2*Tetra::N);
mtl::dense_vector <double> L(3*Tetra::N), Lp(2*Tetra::N);

mtl::dense2D <double> P(2*Tetra::N,3*Tetra::N);

Settings *mySettings = settings;
double H[DIM];
H[0] = Hext[0];H[1] = Hext[1];H[2] = Hext[2];
double vz = DW_vz;


for_each(refTet->begin(),refTet->end(),
    [mySettings,&H,vz,dt,&Lw,&K,&L,&Kp,&Lp,&ins,NOD,&P](Tetra::Tet & tet)
        {
        mtl::mat::set_to_zero(K); mtl::mat::set_to_zero(Kp);
        mtl::vec::set_to_zero(L); mtl::vec::set_to_zero(Lp);
        tet.integrales(mySettings->paramTetra,H,vz,mySettings->theta,dt,mySettings->TAUR,K, L);     
        tet.projection( P, K, L, Kp, Lp);
        tet.assemblage( ins,NOD, Kp, Lp, Lw);
        //projection<Tetra::Tet>(tet, P, K, L, Kp, Lp);
        //assemblage<Tetra::Tet>(tet, ins,NOD, Kp, Lp, Lw);//Kw avant dernier
        }
);

counter.tac();

if(VERBOSE) { std::cout << "vector<tetra> done. elapsed time = " << counter.elapsed() << "s" << std::endl; }

counter.tic();

mtl::dense2D <double> Ps(2*Facette::N,3*Facette::N);

mtl::dense2D <double> Ks(3*Facette::N,3*Facette::N), Ksp(2*Facette::N,2*Facette::N);
mtl::dense_vector <double> Ls(3*Facette::N), Lsp(2*Facette::N);

for_each(refFac->begin(),refFac->end(),
    [this,&Lw,&Ks,&Ls,&Ksp,&Lsp,&ins,NOD,&Ps](Facette::Fac &fac)
        {
        mtl::mat::set_to_zero(Ks); mtl::mat::set_to_zero(Ksp);
        mtl::vec::set_to_zero(Ls); mtl::vec::set_to_zero(Lsp);
        fac.integrales(settings->paramFacette,*refNode, Ls);     
        projection<Facette::Fac>(fac, Ps, Ks, Ls, Ksp, Lsp);
        assemblage<Facette::Fac>(fac, ins, NOD, Ksp, Lsp, Lw);//Kw avant dernier    
        }
);

delete ins;//Kw should be ready

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
