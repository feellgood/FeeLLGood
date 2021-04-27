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

/*
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
*/
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
double r1,r2;
    
if(!determinist)
    {    
    std::mt19937 gen(rd());// random number generator: standard Mersenne twister initialized with seed rd()
    std::uniform_real_distribution<> distrib(0.0,1.0);
    
    r1 = distrib(gen);
    r2 = distrib(gen);
    }
else
    {
    r1 = rand() / (RAND_MAX+1.);
    r2 = rand() / (RAND_MAX+1.);
    }

for(int i=0; i < NOD ; i++)
    {
    refMsh->setNode(i).theta_sph = M_PI*r1;
    refMsh->setNode(i).phi_sph = M_2_PI*r2;
    refMsh->setNode(i).calc_ep();
    refMsh->setNode(i).calc_eq();
    }
}
