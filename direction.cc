#include "fem.h"

using namespace std;

void Fem::direction(void)
{
//IF_VERBOSE() cerr << "direction " << endl;
const int NPS=1;

ANNidxArray nnIdx = new ANNidx[NPS];    if(!nnIdx) SYSTEM_ERROR;
ANNdistArray dists = new ANNdist[NPS];  if(!dists) SYSTEM_ERROR;
ANNpoint queryPt= annAllocPt(3);

/* bord gauche */
queryPt[0]=cx;
queryPt[1]=cy;
queryPt[2]=cz-lz/2.;

kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
int ns=nnIdx[0];
//IF_VERBOSE() cerr << "ns : " << ns << endl;

double u2L= node[ns].u[2];
IF_VERBOSE() cout << "z : " << queryPt[2] << " u2 : "<< u2L << endl;

/* bord droit */
queryPt[0]=cx;
queryPt[1]=cy;
queryPt[2]=cz+lz/2.;

kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
ns=nnIdx[0];
double u2R= node[ns].u[2];
IF_VERBOSE() cout << "z : " << queryPt[2] << " u2 : "<< u2R << endl;

if (u2L*u2R>0.){
   IF_VERBOSE() cout << "Warning apparently no DW!" << endl;
   }

DW_dir=(u2L>0? 1.: -1.); /* sens de deplacement de la paroi +Oz ou -Oz */
IF_VERBOSE(){cout << "DW dir   : " << DW_dir << endl;//cerr << "fin direction " << endl;
	}
}

