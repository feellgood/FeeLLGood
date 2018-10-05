#include <algorithm>

#include "linear_algebra.h"

#include "gmm_iter.h"
#include "gmm_solver_bicgstab.h"
#include "gmm_solver_gmres.h"

int LinAlgebra::vsolve(double dt,long nt)
{
time_t timeStart;

const int NOD = refNode->size();

base_projection();   // definit plan tangent

write_matrix Kw(2*NOD, 2*NOD);
write_vector Lw(2*NOD);

if(VERBOSE) { std::cout <<"DW speed = " << DW_vz << "\nmatrix size " << 2*NOD << "; assembling ..." << std::endl; }

time(&timeStart);

for_each(refTet->begin(),refTet->end(),
    [this,dt,&Kw,&Lw](Tetra::Tet & tet)
        { 
        gmm::dense_matrix <double> K(3*Tetra::N,3*Tetra::N), Kp(2*Tetra::N,2*Tetra::N);
        std::vector <double> L(3*Tetra::N), Lp(2*Tetra::N);
        tet.integrales(settings->paramTetra,*refNode,Hext,DW_vz,settings->theta,dt,settings->TAUR,K, L);     
        projection<Tetra::Tet>(tet, K, L, Kp, Lp);
        assemblage<Tetra::Tet>(tet, Kp, Lp, Kw, Lw);
        }
);

for_each(refFac->begin(),refFac->end(),
    [this,&Kw,&Lw](Facette::Fac &fac)
        {
        gmm::dense_matrix <double> K(3*Facette::N,3*Facette::N), Kp(2*Facette::N,2*Facette::N);
        std::vector <double> L(3*Facette::N), Lp(2*Facette::N);
        fac.integrales(settings->paramFacette,*refNode, L);     
        projection<Facette::Fac>(fac, K, L, Kp, Lp);
        assemblage<Facette::Fac>(fac, Kp, Lp, Kw, Lw);    
        }
);

time_t timeEnd;
time(&timeEnd);

if(VERBOSE) { std::cout << "elapsed time = " << difftime(timeEnd,timeStart) << "s" << std::endl; }

read_matrix Kr(2*NOD,2*NOD);    gmm::copy(Kw, Kr);
read_vector Lr(2*NOD);          gmm::copy(Lw, Lr);
write_vector Xw(2*NOD);

gmm::iteration bicg_iter(1e-6);
gmm::iteration gmr_iter(1e-6);
bicg_iter.set_maxiter(MAXITER);
gmr_iter.set_maxiter(MAXITER);
bicg_iter.set_noisy(VERBOSE);
gmr_iter.set_noisy(VERBOSE);

time(&timeStart);

if (!nt) 
    {
    if(VERBOSE) { std::cout << "computing prc.";std::fflush(NULL); }

    prc = new gmm::diagonal_precond <read_matrix> (Kr); // prc est public dans linAlgebra
	time(&timeEnd);    
	if(VERBOSE) { std::cout << "elapsed time = " << difftime(timeEnd,timeStart) << "s" << std::endl; }
    }
else if (!(nt % REFRESH_PRC))
    {
    delete prc;
    if(VERBOSE) { std::cout << "computing prc.";std::fflush(NULL); }
    prc = new gmm::diagonal_precond <read_matrix> (Kr);
    time(&timeEnd);    
	if(VERBOSE) { std::cout << "elapsed time = " << difftime(timeEnd,timeStart) << "s" << std::endl; }
    } 

time(&timeStart);
gmm::bicgstab(Kr, Xw, Lr, *prc, bicg_iter);

if(VERBOSE) {std::cout<< "bi-conjugate gradient stabilized done." <<std::endl;}

if (!(bicg_iter.converged())) {
    time(&timeEnd);
	if(VERBOSE) { std::cout << "%5t  bicg FAILED in " << bicg_iter.get_iteration() << " %30T. " << difftime(timeEnd,timeStart) << " s" << std::endl; }
    gmm::diagonal_precond <read_matrix>  gmr_prc (Kr);
    time(&timeStart);
    gmm::clear(Xw);
    gmm::gmres(Kr, Xw, Lr, gmr_prc, 50, gmr_iter);
    if (!(gmr_iter.converged())) {
	time(&timeEnd);
    if(VERBOSE) { std::cout << "gmres FAILED in " << gmr_iter.get_iteration() << " difftime: " << difftime(timeEnd,timeStart) << " s" << std::endl; }
    return 1;
    }
    else {time(&timeEnd);
        if(VERBOSE) { std::cout << "v-solve in " << gmr_iter.get_iteration() << " (gmres) difftime= " << difftime(timeEnd,timeStart) << " s" << std::endl; }
        }
    }
else {time(&timeEnd);
    if(VERBOSE) { std::cout << "%5t v-solve in " << bicg_iter.get_iteration() << " (bicg,prc-" << (nt % REFRESH_PRC) << ") .......... " 
<< difftime(timeEnd,timeStart) << " s" << std::endl; }
    }

read_vector Xr(2*NOD);    gmm::copy(Xw, Xr);

double v2max = 0.0;

int i=0;
for_each(refNode->begin(),refNode->end(),
    [&i,&v2max,&Xr,NOD,dt](Node &n)
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
