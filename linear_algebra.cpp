#include "linear_algebra.h"

void LinAlgebra::updateNodes(std::vector<double> const& X,const double dt)
{
double v2max = 0.0;

int i=0;
std::for_each(refNode->begin(),refNode->end(),
    [this,&v2max,&i,&dt,&X](Nodes::Node &n)
        {
        double vp = X[i];
        double vq = X[NOD+i];
        double v2 = vp*vp + vq*vq;
        if (v2>v2max) { v2max = v2; }
        n.make_evol(vp,vq,dt);    
        i++;    
        }
);//end for_each

v_max = sqrt(v2max);
}

void LinAlgebra::prepareItTet(std::vector <Tetra::Tet> &myTet)
{
const size_t block_size = myTet.size()/NbTH;
const int extra_blocks = myTet.size()%NbTH;
std::vector<Tetra::Tet>::iterator it = myTet.begin();

for (int i = 0; i < NbTH; i++)
    {
    refTetIt[i].first = it;
    std::advance(it, block_size + (i<extra_blocks ? 1 : 0));
    refTetIt[i].second = it;
    }
assert(refTetIt[NbTH-1].second == myTet.end());
}

void LinAlgebra::base_projection(bool determinist)
{
if(!determinist)
    {    
    std::mt19937 gen(rd());// random number generator: standard Mersenne twister initialized with seed rd()
    std::uniform_real_distribution<> distrib(0.0,1.0);
    
    std::for_each(refNode->begin(),refNode->end(),[&gen,&distrib](Nodes::Node &n) 
    { 
        n.theta_sph = M_PI * distrib(gen);
        n.phi_sph = M_2_PI * distrib(gen);
    }); 
    }
else
    {
    std::for_each(refNode->begin(),refNode->end(),[](Nodes::Node &n) 
    { 
        n.theta_sph = M_PI * rand() / (RAND_MAX+1.);
        n.phi_sph = M_2_PI * rand() / (RAND_MAX+1.);
    });
    }
}
