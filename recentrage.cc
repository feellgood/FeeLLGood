#include "fem.h"

/*-----------------------------------------*/
/*   recentrage de la paroi selon Oz       */
/*-----------------------------------------*/

bool Fem::recentrage(double thres,double mz) // abs(thres) < 1
{
time_t timeStart;
time(&timeStart);

thres = std::min(abs(thres), 1.);
if (fabs(mz)<thres) return false;

if(VERBOSE) {std::cout << "centering ..." << std::flush;}

const int NPS=1;
int ns;

ANNidxArray nnIdx = new ANNidx[NPS];    if(!nnIdx) SYSTEM_ERROR;
ANNdistArray dists = new ANNdist[NPS];  if(!dists) SYSTEM_ERROR;
ANNpoint queryPt= annAllocPt(3);

queryPt[0]=c.x();
queryPt[1]=c.y();

// bord gauche
queryPt[2]=c.z()-l.z()/2.;
kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
ns=nnIdx[0]; 
//double u0L=fem.node[ns].u[0];
//double u1L=fem.node[ns].u[1];
double u2L=node[ns].u.z();

// centre
queryPt[2]=c.z();
kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
ns=nnIdx[0]; 
//double u0C=fem.node[ns].u[0];
//double u1C=fem.node[ns].u[1];
//double u2C=fem.node[ns].u[2];

// bord droit
queryPt[2]=c.z()+l.z()/2.;
kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
ns=nnIdx[0]; 
//double u0R=fem.node[ns].u[0];
//double u1R=fem.node[ns].u[1];
double u2R=node[ns].u.z();

if (u2L*u2R>0) {
#ifdef LIBRARY
    throw runtime_error("Error No Domain Wall");
#else
    std::cout << "Error No Domain Wall" << std::endl;
    exit(1);
#endif
   }

assert(u2L*u2R<0);
double Dz= mz*l.z()/2.*u2L;  // decalage avec signe OK

/* cas ou Dz>0				cas ou Dz<0

<----------------|------->		------->|<----------------	mz = <uz> < 0

ou					ou

---------------->|<-------		<-------|---------------->	mz = <uz> > 0

*/


for (int i=0; i<NOD; i++){    
    Node &tgt_node = node[i];
    double x=tgt_node.p.x();
    double y=tgt_node.p.y();
    double z=tgt_node.p.z()+Dz;

    if (z-c.z() > +l.z()/2.) z=c.z()+l.z()/2.;
    if (z-c.z() < -l.z()/2.) z=c.z()-l.z()/2.;

    queryPt[0]=x;
    queryPt[1]=y;
    queryPt[2]=z;
    kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);

    int ns=nnIdx[0]; 
    tgt_node.u0 = node[ns].u;     
    tgt_node.v  = Pt::pt3D(0.,0.,0.);     
    }

delete queryPt;
delete [] nnIdx;
delete [] dists;

time_t timeEnd;
time(&timeEnd);
if(VERBOSE) { std::cout << "elapsed time = "<< difftime(timeEnd,timeStart) << "s" << std::endl; }
return true;
}
