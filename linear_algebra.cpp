#include "linear_algebra.h"


void LinAlgebra::deepCopyTet(std::vector <Tetra::Tet> const& myTet)
{
const unsigned long block_size = std::distance(myTet.begin(),myTet.end())/NbTH;

std::vector<Tetra::Tet>::const_iterator it_begin = myTet.begin();

for(int i=0;i<(NbTH-1);i++) 
    {
    std::vector<Tetra::Tet>::const_iterator it_end = it_begin;
    std::advance(it_end,block_size);
    refTet[i].resize(block_size,Tetra::Tet(NOD));
    std::copy( it_begin, it_end, refTet[i].begin() );
    it_begin = it_end;
    }
const unsigned long last_block_size = std::distance(it_begin,myTet.end());
refTet[NbTH-1].resize(last_block_size,Tetra::Tet(NOD));
std::copy( it_begin, myTet.end(), refTet[NbTH-1].begin() );    
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
