#include "fem.h"
#include "tiny.h"

/*-----------------------------------------*/
/*   recentrage de la paroi selon Oz       */
/*-----------------------------------------*/

bool recentrage(Fem &fem, double thres) // abs(thres) < 1
{
time_t timeStart;
time(&timeStart);

double mz =u_moy(fem, 2);
thres=min(abs(thres), 1.);
if (fabs(mz)<thres) return false;

//IF_VERBOSE(fem) cout << boost::format("%5t centering %30T.") <<flush;// *ct*
IF_VERBOSE(fem) cout << "%5t centering %30T." << flush;

const int NOD = fem.NOD;
const int NPS=1;
int ns;

ANNidxArray nnIdx = new ANNidx[NPS];    if(!nnIdx) SYSTEM_ERROR;
ANNdistArray dists = new ANNdist[NPS];  if(!dists) SYSTEM_ERROR;
ANNpoint queryPt= annAllocPt(3);

queryPt[0]=fem.cx;
queryPt[1]=fem.cy;

// bord gauche
queryPt[2]=fem.cz-fem.lz/2.;
fem.kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
ns=nnIdx[0]; 
//double u0L=fem.node[ns].u[0];
//double u1L=fem.node[ns].u[1];
double u2L=fem.node[ns].u[2];

// centre
queryPt[2]=fem.cz;
fem.kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
ns=nnIdx[0]; 
//double u0C=fem.node[ns].u[0];
//double u1C=fem.node[ns].u[1];
//double u2C=fem.node[ns].u[2];

// bord droit
queryPt[2]=fem.cz+fem.lz/2.;
fem.kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
ns=nnIdx[0]; 
//double u0R=fem.node[ns].u[0];
//double u1R=fem.node[ns].u[1];
double u2R=fem.node[ns].u[2];

if (u2L*u2R>0) {
#ifdef LIBRARY
    throw runtime_error("Error No Domain Wall");
#else
    cout << "Error No Domain Wall" << endl;
    exit(1);
#endif
   }

assert(u2L*u2R<0);
double Dz= mz*fem.lz/2.*u2L;  // decalage avec signe OK

/* cas ou Dz>0				cas ou Dz<0

<----------------|------->		------->|<----------------	mz = <uz> < 0

ou					ou

---------------->|<-------		<-------|---------------->	mz = <uz> > 0

*/


for (int i=0; i<NOD; i++){    
    Node &tgt_node = fem.node[i];
    double x=tgt_node.x;
    double y=tgt_node.y;
    double z=tgt_node.z+Dz;

    if (z-fem.cz>+fem.lz/2.) z=fem.cz+fem.lz/2.;
    if (z-fem.cz<-fem.lz/2.) z=fem.cz-fem.lz/2.;

    queryPt[0]=x;
    queryPt[1]=y;
    queryPt[2]=z;
    fem.kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);

    int ns=nnIdx[0]; 
    Node &src_node = fem.node[ns];
    tgt_node.u0[0] = src_node.u[0];     
    tgt_node.u0[1] = src_node.u[1];
    tgt_node.u0[2] = src_node.u[2];
    tgt_node.v[0]  = 0.;     
    tgt_node.v[1]  = 0.;
    tgt_node.v[2]  = 0.;
    }

delete queryPt;
delete [] nnIdx;
delete [] dists;

time_t timeEnd;
time(&timeEnd);
IF_VERBOSE(fem) cout << "elapsed time = "<< difftime(timeEnd,timeStart) << "s" <<endl;
return true;
}
