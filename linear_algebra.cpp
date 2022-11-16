#include "linear_algebra.h"

void LinAlgebra::updateNodes(std::vector<double> const& X,const double dt)
{
double v2max = 0.0;

for(int i=0; i < NOD ; i++)
    {
    double vp = X[i];
    double vq = X[NOD+i];
    double v2 = vp*vp + vq*vq;
    if (v2>v2max) { v2max = v2; }
    refMsh->setNode(i).make_evol(vp,vq,dt);    
    }

v_max = sqrt(v2max);
}

void LinAlgebra::base_projection(bool determinist)
{
double r;
    
if(!determinist)
    {    
    std::mt19937 gen(rd());// random number generator: standard Mersenne twister initialized with seed rd()
    std::uniform_real_distribution<> distrib(0.0,1.0);
    r = distrib(gen);
    }
else
    {
    r = rand() / (RAND_MAX+1.);
    }

for(int i=0; i < NOD ; i++) refMsh->setNode(i).setBasis(M_2_PI*r);
}

void LinAlgebra::prepareElements(Pt::pt3D const& Hext /**< [in] applied field */, timing const& t_prm /**< [in] */)
{
base_projection(!RAND_DETERMINIST);

std::for_each(std::execution::par,refMsh->tet.begin(),refMsh->tet.end(),[this,&Hext,&t_prm](Tetra::Tet & tet)
	{
	double K[3*Tetra::N][3*Tetra::N] = { {0} }; 
        Pt::pt3D L[Tetra::N];
            
        tet.integrales(settings.paramTetra,t_prm,Hext,settings.recentering_direction,DW_vz,K, L);
        tet.projection(K,L);
	});

std::for_each(std::execution::par,refMsh->fac.begin(),refMsh->fac.end(),[this](Facette::Fac & fac)
	{
	double Ks[3*Facette::N][3*Facette::N] = { {0} };
        Pt::pt3D Ls[Facette::N];
        
        fac.integrales(settings.paramFacette,Ls);
        fac.projection(Ks,Ls);
	});
}

