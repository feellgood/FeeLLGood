#include <functional>

#include "fem.h"

using namespace Pt;

void Fem::direction(enum index idx_dir)
{
const int NPS=1;

ANNidxArray nnIdx = new ANNidx[NPS];    if(!nnIdx) SYSTEM_ERROR;
ANNdistArray dists = new ANNdist[NPS];  if(!dists) SYSTEM_ERROR;
ANNpoint queryPt= annAllocPt(Pt::DIM);

pt3D qPt = msh.c - 0.5*pDirect( pt3D(idx_dir),msh.l);

/* left */
queryPt[0]=qPt.x(); queryPt[1]=qPt.y(); queryPt[2]=qPt.z();

kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
int ns=nnIdx[0];
double u2L= msh.getNode(ns).u(idx_dir);

/* right */
qPt = msh.c + 0.5*pDirect( pt3D(idx_dir),msh.l);
queryPt[0]=qPt.x(); queryPt[1]=qPt.y(); queryPt[2]=qPt.z();

kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
ns=nnIdx[0];
double u2R= msh.getNode(ns).u(idx_dir);

if (u2L*u2R>0.) DW_dir = 0.0; else DW_dir=(u2L>0? 1.: -1.); /* sens de deplacement de la paroi +Oz ou -Oz */
delete[] dists;
delete[] nnIdx;
annDeallocPt(queryPt);
}


bool Fem::recenter(double thres,char recentering_direction)
{
enum index idx_dir;
    
switch(recentering_direction)
    {
    case 'X':idx_dir = Pt::IDX_X;break;
    case 'Y':idx_dir = Pt::IDX_Y;break;
    case 'Z':idx_dir = Pt::IDX_Z;break;
    default:idx_dir = Pt::IDX_Z;break;
    }
            
thres = abs(thres);
double m_i = msh.avg(Nodes::get_u_comp,idx_dir);
if (fabs(m_i)<thres) return false;

const int NPS=1;
int ns;

ANNidxArray nnIdx = new ANNidx[NPS];    if(!nnIdx) SYSTEM_ERROR;
ANNdistArray dists = new ANNdist[NPS];  if(!dists) SYSTEM_ERROR;
ANNpoint queryPt= annAllocPt(Pt::DIM);

pt3D p_dir = pt3D(idx_dir);//unit vector
pt3D qPt = msh.c - 0.5*pDirect(p_dir,msh.l);

/* bord gauche */
queryPt[0]= qPt.x(); queryPt[1]= qPt.y(); queryPt[2]= qPt.z();

kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
ns=nnIdx[0]; 
double u2L=msh.getNode(ns).u(idx_dir);

// bord droit
qPt = msh.c + 0.5*pDirect(p_dir,msh.l);
queryPt[0]= qPt.x(); queryPt[1]= qPt.y(); queryPt[2]= qPt.z();

kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);
ns=nnIdx[0]; 
double u2R=msh.getNode(ns).u(idx_dir);

if (u2L*u2R>0)
    {
    std::cout << "Error No Domain Wall" << std::endl;
    SYSTEM_ERROR;
   }

assert(u2L*u2R<0);
double D_i= m_i*0.5*msh.l(idx_dir)*u2L;  // decalage avec signe OK

for (int i=0; i<msh.getNbNodes(); i++)
    {
    pt3D p = msh.getNode(i).p + D_i*p_dir;
    
    if (p(idx_dir)-msh.c(idx_dir) > +0.5*msh.l(idx_dir)) { p = msh.c + 0.5*msh.l(idx_dir)*p_dir;}
    if (p(idx_dir)-msh.c(idx_dir) < -0.5*msh.l(idx_dir)) {p = msh.c - 0.5*msh.l(idx_dir)*p_dir;}

    queryPt[0]=p.x(); queryPt[1]=p.y(); queryPt[2]=p.z();
    kdtree->annkSearch(queryPt, NPS, nnIdx, dists, 0.);

    int ns=nnIdx[0]; 
    //msh.setNode(i).u0 = msh.getNode(ns).u;     
    msh.set_node_u0(i, msh.getNode(ns).u );
    //msh.setNode(i).v  = Pt::pt3D(0.,0.,0.);     
    msh.set_node_zero_v(i);
    }
    
annDeallocPt(queryPt);
delete [] nnIdx;
delete [] dists;

return true;
}
