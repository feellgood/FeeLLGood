#include "fem.h"

void Fem::direction(void)
{
if(VERBOSE) { std::cerr << "direction " << std::endl; }
const int NPS=1;

ANNidxArray nnIdx = new ANNidx[NPS];    if(!nnIdx) SYSTEM_ERROR;
ANNdistArray dists = new ANNdist[NPS];  if(!dists) SYSTEM_ERROR;
ANNpoint queryPt= annAllocPt(3);

/* bord gauche */
queryPt[0]=c.x();
queryPt[1]=c.y();
queryPt[2]=c.z()-l.z()/2.;

kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
int ns=nnIdx[0];
if(VERBOSE) { std::cerr << "ns : " << ns << std::endl; }

double u2L= node[ns].u.z();
if(VERBOSE) { std::cout << "z : " << queryPt[2] << " u2 : "<< u2L << std::endl; }

/* bord droit */
queryPt[0]=c.x();
queryPt[1]=c.y();
queryPt[2]=c.z()+l.z()/2.;

kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
ns=nnIdx[0];
double u2R= node[ns].u.z();
if(VERBOSE) { std::cout << "z : " << queryPt[2] << " u2 : "<< u2R << std::endl; }

if ((u2L*u2R>0.)&&VERBOSE){ std::cout << "Warning apparently no DW!" << std::endl; }

DW_dir=(u2L>0? 1.: -1.); /* sens de deplacement de la paroi +Oz ou -Oz */
if(VERBOSE) {std::cout << "DW dir   : " << DW_dir << std::endl; }
}

