#include "fem.h"

/*-----------------------------------------*/
/*   recentering of the DW along idx_dir       */
/*-----------------------------------------*/

using namespace Pt;

void Fem::direction(enum index idx_dir)
{
if(VERBOSE) { std::cerr << "direction " << std::endl; }
const int NPS=1;

ANNidxArray nnIdx = new ANNidx[NPS];    if(!nnIdx) SYSTEM_ERROR;
ANNdistArray dists = new ANNdist[NPS];  if(!dists) SYSTEM_ERROR;
ANNpoint queryPt= annAllocPt(Pt::DIM);

pt3D p_dir = pt3D(idx_dir);//unit vector

/* bord gauche */
queryPt[0]=c.x()-0.5*pScal(p_dir,l);
queryPt[1]=c.y()-0.5*pScal(p_dir,l);
queryPt[2]=c.z()-0.5*pScal(p_dir,l);

kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
int ns=nnIdx[0];
if(VERBOSE) { std::cerr << "ns : " << ns << std::endl; }

double u2L= node[ns].u(idx_dir);
if(VERBOSE) { std::cout << "z : " << queryPt[2] << " u2 : "<< u2L << std::endl; }

/* bord droit */
queryPt[0]=c.x()+0.5*pScal(p_dir,l);
queryPt[1]=c.y()+0.5*pScal(p_dir,l);
queryPt[2]=c.z()+0.5*pScal(p_dir,l);

kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
ns=nnIdx[0];
double u2R= node[ns].u(idx_dir);
if(VERBOSE) { std::cout << "z : " << queryPt[idx_dir] << " u2 : "<< u2R << std::endl; }

if ((u2L*u2R>0.)&&VERBOSE){ std::cout << "Warning apparently no DW!" << std::endl; }

DW_dir=(u2L>0? 1.: -1.); /* sens de deplacement de la paroi +Oz ou -Oz */
if(VERBOSE) {std::cout << "DW dir   : " << DW_dir << std::endl; }
}



bool Fem::recentrage(double thres,enum index idx_dir,double m_i)
{
time_t timeStart;
time(&timeStart);

thres = std::min(abs(thres), 1.);
if (fabs(m_i)<thres) return false;

if(VERBOSE) {std::cout << "centering ..." << std::flush;}

const int NPS=1;
int ns;

ANNidxArray nnIdx = new ANNidx[NPS];    if(!nnIdx) SYSTEM_ERROR;
ANNdistArray dists = new ANNdist[NPS];  if(!dists) SYSTEM_ERROR;
ANNpoint queryPt= annAllocPt(Pt::DIM);

pt3D p_dir = pt3D(idx_dir);//unit vector

/* bord gauche */
queryPt[0]=c.x()-0.5*pScal(p_dir,l);
queryPt[1]=c.y()-0.5*pScal(p_dir,l);
queryPt[2]=c.z()-0.5*pScal(p_dir,l);

kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
ns=nnIdx[0]; 
double u2L=node[ns].u(idx_dir);

// centre
queryPt[0]=c.x();
queryPt[1]=c.y();
queryPt[2]=c.z();
kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
ns=nnIdx[0];

// bord droit
queryPt[0]=c.x()+0.5*pScal(p_dir,l);
queryPt[1]=c.y()+0.5*pScal(p_dir,l);
queryPt[2]=c.z()+0.5*pScal(p_dir,l);

kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
ns=nnIdx[0]; 
double u2R=node[ns].u(idx_dir);

if (u2L*u2R>0) {
#ifdef LIBRARY
    throw runtime_error("Error No Domain Wall");
#else
    std::cout << "Error No Domain Wall" << std::endl;
    exit(1);
#endif
   }

assert(u2L*u2R<0);
double D_i= m_i*0.5*l(idx_dir)*u2L;  // decalage avec signe OK



for (int i=0; i<NOD; i++){    
    Node &tgt_node = node[i];
    pt3D p = tgt_node.p + D_i*p_dir;
    
    if (p(idx_dir)-c(idx_dir) > +0.5*l(idx_dir)) { p = c + 0.5*l(idx_dir)*p_dir;}
    if (p(idx_dir)-c(idx_dir) < -0.5*l(idx_dir)) {p = c - 0.5*l(idx_dir)*p_dir;}

    queryPt[0]=p.x();
    queryPt[1]=p.y();
    queryPt[2]=p.z();
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
