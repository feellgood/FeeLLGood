#include <functional>

#include "fem.h"

using namespace Pt;

void Fem::direction(enum index idx_dir)
{
const int NPS=1;

ANNidxArray nnIdx = new ANNidx[NPS];    if(!nnIdx) SYSTEM_ERROR;
ANNdistArray dists = new ANNdist[NPS];  if(!dists) SYSTEM_ERROR;
ANNpoint queryPt= annAllocPt(Pt::DIM);

pt3D qPt = c - 0.5*pDirect( pt3D(idx_dir),l);

/* left */
queryPt[0]=c.x(); queryPt[1]=c.y(); queryPt[2]=c.z();

kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
int ns=nnIdx[0];

double u2L= node[ns].u(idx_dir);

/* right */
qPt = c + 0.5*pDirect( pt3D(idx_dir),l);
queryPt[0]=c.x(); queryPt[1]=c.y(); queryPt[2]=c.z();

kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
ns=nnIdx[0];
double u2R= node[ns].u(idx_dir);

if (u2L*u2R>0.) DW_dir = 0.0; else DW_dir=(u2L>0? 1.: -1.); /* sens de deplacement de la paroi +Oz ou -Oz */
delete[] dists;
delete[] nnIdx;
annDeallocPt(queryPt);
}



bool Fem::recentrage(double thres,enum index idx_dir)
{
thres = abs(thres);
double m_i = Fem::avg(Nodes::get_u_comp,idx_dir);
if (fabs(m_i)<thres) return false;

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

if (u2L*u2R>0)
    {
    std::cout << "Error No Domain Wall" << std::endl;
    SYSTEM_ERROR;
   }

assert(u2L*u2R<0);
double D_i= m_i*0.5*l(idx_dir)*u2L;  // decalage avec signe OK

const int NOD = node.size();
for (int i=0; i<NOD; i++){    
    Nodes::Node &tgt_node = node[i];
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
    
annDeallocPt(queryPt);
delete [] nnIdx;
delete [] dists;

return true;
}
