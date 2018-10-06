#include <algorithm>

#include "linear_algebra.h"


int LinAlgebra::vsolve(double dt,long nt)
{
time_t timeStart;

const int NOD = refNode->size();

base_projection();   // definit plan tangent

mtl::compressed2D<double>  Kw(2*NOD, 2*NOD);
mtl::mat::set_to_zero(Kw);
mtl::dense_vector<double> Lw(2*NOD);
mtl::vec::set_to_zero(Lw);


if(VERBOSE) { std::cout <<"DW speed = " << DW_vz << "\nmatrix size " << 2*NOD << "; assembling ..." << std::endl; }

time(&timeStart);

{
mtl::mat::inserter< mtl::compressed2D<double>,mtl::update_plus<v_type> > ins(Kw,64); //64 maybe too small

for_each(refTet->begin(),refTet->end(),
    [this,dt,&Kw,&Lw,&ins](Tetra::Tet & tet)
        { 
        mtl::dense2D <double> K(3*Tetra::N,3*Tetra::N), Kp(2*Tetra::N,2*Tetra::N);
        mtl::dense_vector <double> L(3*Tetra::N), Lp(2*Tetra::N);
        tet.integrales(settings->paramTetra,*refNode,Hext,DW_vz,settings->theta,dt,settings->TAUR,K, L);     
        projection<Tetra::Tet>(tet, K, L, Kp, Lp);
        assemblage<Tetra::Tet>(tet,ins, Kp, Lp, Kw, Lw);
        }
);

for_each(refFac->begin(),refFac->end(),
    [this,&Kw,&Lw,&ins](Facette::Fac &fac)
        {
        mtl::dense2D <double> K(3*Facette::N,3*Facette::N), Kp(2*Facette::N,2*Facette::N);
        mtl::dense_vector <double> L(3*Facette::N), Lp(2*Facette::N);
        fac.integrales(settings->paramFacette,*refNode, L);     
        projection<Facette::Fac>(fac, K, L, Kp, Lp);
        assemblage<Facette::Fac>(fac,ins, Kp, Lp, Kw, Lw);    
        }
);
}//inserter destructor called 
time_t timeEnd;
time(&timeEnd);

if(VERBOSE) { std::cout << "elapsed time = " << difftime(timeEnd,timeStart) << "s" << std::endl; }

mtl::dense_vector<double> Xw(2*NOD);

itl::noisy_iteration<double> bicg_iter(Lw,MAXITER,1e-6);

itl::noisy_iteration<double> gmr_iter(Lw,MAXITER,1e-6);

//bicg_iter.set_maxiter(MAXITER);
//gmr_iter.set_maxiter(MAXITER);

bicg_iter.set_quite(VERBOSE);
gmr_iter.set_quite(VERBOSE);

time(&timeStart);

itl::pc::diagonal < mtl::compressed2D<double> > prc(Kw);

/*
if (!nt) 
    {
    if(VERBOSE) { std::cout << "computing prc.";std::fflush(NULL); }

    itl::pc::diagonal < mtl::compressed2D<double> > prc(Kw); 
	time(&timeEnd);    
	if(VERBOSE) { std::cout << "elapsed time = " << difftime(timeEnd,timeStart) << "s" << std::endl; }
    }
else if (!(nt % REFRESH_PRC))
    {
    delete prc;
    if(VERBOSE) { std::cout << "computing prc.";std::fflush(NULL); }
    itl::pc::diagonal < mtl::compressed2D<double> > prc(Kw); 
    time(&timeEnd);    
	if(VERBOSE) { std::cout << "elapsed time = " << difftime(timeEnd,timeStart) << "s" << std::endl; }
    } 
*/


time(&timeStart);
bicgstab(Kw, Xw, Lw, prc, bicg_iter);

if(VERBOSE) {std::cout<< "bi-conjugate gradient stabilized done." <<std::endl;}

if (!(bicg_iter.is_converged() )) {
    time(&timeEnd);
	if(VERBOSE) { std::cout << "%5t  bicg FAILED in " << bicg_iter.iterations() << " %30T. " << difftime(timeEnd,timeStart) << " s" << std::endl; }
    itl::pc::diagonal < mtl::compressed2D<double> >  gmr_prc (Kw);
    time(&timeStart);
    //gmm::clear(Xw);
    itl::gmres(Kw, Xw, Lw, gmr_prc, gmr_prc,gmr_iter, gmr_iter);
    if (!(gmr_iter.is_converged() )) {
	time(&timeEnd);
    if(VERBOSE) { std::cout << "gmres FAILED in " << gmr_iter.iterations() << " difftime: " << difftime(timeEnd,timeStart) << " s" << std::endl; }
    return 1;
    }
    else {time(&timeEnd);
        if(VERBOSE) { std::cout << "v-solve in " << gmr_iter.iterations() << " (gmres) difftime= " << difftime(timeEnd,timeStart) << " s" << std::endl; }
        }
    }
else {time(&timeEnd);
    if(VERBOSE) { std::cout << "%5t v-solve in " << bicg_iter.iterations() << " (bicg,prc-" << (nt % REFRESH_PRC) << ") .......... " 
<< difftime(timeEnd,timeStart) << " s" << std::endl; }
    }

//read_vector Xr(2*NOD);    gmm::copy(Xw, Xr);

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
