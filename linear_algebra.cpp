#include "linear_algebra.h"

void LinAlgebra::prepareItTet(std::vector <Tetra::Tet> &myTet)
{
const unsigned long block_size = std::distance(myTet.begin(),myTet.end())/NbTH;
std::vector<Tetra::Tet>::iterator it = myTet.begin();

std::for_each(refTetIt.begin(),refTetIt.end(),[&it,myTet,block_size] (std::pair<std::vector<Tetra::Tet>::iterator,std::vector<Tetra::Tet>::iterator> & p_it)
    {
    p_it.first = it;
    std::advance(it,block_size);
    p_it.second = it;            
    }
    );
refTetIt[NbTH-1].second = std::prev(myTet.end());
}

void LinAlgebra::base_projection(void)
{
    std::mt19937 gen(rd());// random number generator: standard Mersenne twister initialized with seed rd()
    std::uniform_real_distribution<> distrib(0.0,1.0);
    std::for_each(refNode->begin(),refNode->end(),[&gen,&distrib](Nodes::Node &n) 
    { 
        double theta = M_PI * distrib(gen);
        double phi = M_2_PI * distrib(gen);
        n.ep = Pt::pt3D(theta,phi)*n.u0;
        n.ep.normalize(); // required because u0 is not necessarily unit vector ?? to check wilth initial expression with exprTK       
        n.eq = n.u0*n.ep; 
        n.eq.normalize();
    }); 
}
