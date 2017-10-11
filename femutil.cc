#include "fem.h"
#include <set>

struct less_than
{
  bool operator()(Fac f1, Fac f2) const
  {
  if (f1.ind[0]<f2.ind[0]) return true;
  else
     if ((f1.ind[0]==f2.ind[0]) && (f1.ind[1]<f2.ind[1])) return true;
     else
        if ((f1.ind[0]==f2.ind[0]) && (f1.ind[1]==f2.ind[1]) && (f1.ind[2]<f2.ind[2])) return true;

  return false;
  }
};

void femutil_node(Fem &fem)
{
const int NOD = fem.NOD;

fem.pts= annAllocPts(NOD, 3);

// calcul du diametre et du centrage
double xmin, xmax, ymin, ymax, zmin, zmax, lx,ly,lz;
xmin = ymin = zmin = +HUGE;
xmax = ymax = zmax = -HUGE;

for (int i=0; i<NOD; i++){
    double xi,yi,zi;
    xi = fem.node[i].x;      yi = fem.node[i].y;    zi = fem.node[i].z;
    fem.pts[i][0]=xi;	     fem.pts[i][1]=yi;	    fem.pts[i][2]=zi;
    if (xi<xmin) xmin=xi;    if (xi>xmax) xmax=xi;
    if (yi<ymin) ymin=yi;    if (yi>ymax) ymax=yi;
    if (zi<zmin) zmin=zi;    if (zi>zmax) zmax=zi;
    }

// allocation de l'arbre de recherche
fem.kdtree = new ANNkd_tree(fem.pts, NOD, 3);
if (!fem.kdtree) SYSTEM_ERROR;

lx=xmax-xmin; fem.lx=lx;
ly=ymax-ymin; fem.ly=ly;
lz=zmax-zmin; fem.lz=lz;

fem.diam = lx;
if (fem.diam<ly) fem.diam=ly;
if (fem.diam<lz) fem.diam=lz;

fem.cx = 0.5*(xmax+xmin);
fem.cy = 0.5*(ymax+ymin);
fem.cz = 0.5*(zmax+zmin);

fem.as[0] = lx/fem.diam;
fem.as[1] = ly/fem.diam;
fem.as[2] = lz/fem.diam;
}

void femutil_tet(Fem &fem)
{
const int TET = fem.TET;

/*
                        v
                      .
                    ,/
                   /
                2(ic)                                 2
              ,/|`\                                 ,/|`\
            ,/  |  `\                             ,/  |  `\
          ,/    '.   `\                         ,6    '.   `5
        ,/       |     `\                     ,/       8     `\
      ,/         |       `\                 ,/         |       `\
     0(ia)-------'.--------1(ib) --> u     0--------4--'.--------1
      `\.         |      ,/                 `\.         |      ,/
         `\.      |    ,/                      `\.      |    ,9
            `\.   '. ,/                           `7.   '. ,/
               `\. |/                                `\. |/
                  `3(id)                                `3
                     `\.
                        ` w

*/

// calcul des volumes et reorientation des tetraedres si necessaire
double voltot = 0.;

for (int t=0; t<TET; t++){
   Tet &tet = fem.tet[t];
   int i0,i1,i2,i3;
   i0=tet.ind[0];   i1=tet.ind[1];   i2=tet.ind[2];   i3=tet.ind[3];
   
   double x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3;
   x0 = fem.node[i0].x;   y0 = fem.node[i0].y;   z0 = fem.node[i0].z;
   x1 = fem.node[i1].x;   y1 = fem.node[i1].y;   z1 = fem.node[i1].z;
   x2 = fem.node[i2].x;   y2 = fem.node[i2].y;   z2 = fem.node[i2].z;
   x3 = fem.node[i3].x;   y3 = fem.node[i3].y;   z3 = fem.node[i3].z;

   double vecx,vecy,vecz,vol;
   vecx = (y1-y0)*(z2-z0)-(y2-y0)*(z1-z0);
   vecy = (z1-z0)*(x2-x0)-(z2-z0)*(x1-x0);
   vecz = (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);
//   vol  = 1./6.* fabs(vecx*(x3-x0) + vecy*(y3-y0) + vecz*(z3-z0));
   vol  = 1./6.* (vecx*(x3-x0) + vecy*(y3-y0) + vecz*(z3-z0));
   if (vol<0.) {
      tet.ind[3]=i2; tet.ind[2]=i3;
      vol=-vol;
      IF_VERBOSE(fem) cout << "ill-oriented tetrahedron: " << t << " now corrected!"<< endl;
      }
   tet.vol = vol;
   voltot+= vol;}
fem.vol = voltot;

}

void femutil_facMs(Fem &fem)
{
const int FAC = fem.FAC;
const int TET = fem.TET;
pair <string,int> p;
map <pair<string,int>,double> &param = fem.param;

// decomposition des tetraedres en elements de surface
set<Fac, less_than> sf;
IF_VERBOSE(fem) cout << "Nb de Tet " << TET << endl;
for (int t=0; t<TET; t++){
    Tet &tet = fem.tet[t];
    int ia,ib,ic,id;
    ia=tet.ind[0];  ib=tet.ind[1];  ic=tet.ind[2];  id=tet.ind[3];
//    cout << "tet " << t <<"/"<<TET<< endl;
    {
    Fac fac; fac.reg=tet.reg; 
    fac.ind[0]=ia; fac.ind[1]=ic; fac.ind[2]=ib;
    sf.insert(fac);
//    cout << fac.ind[0] << " " << fac.ind[1] << " " << fac.ind[2] << endl;
    }
 
    {
    Fac fac; fac.reg=tet.reg;
    fac.ind[0]=ib; fac.ind[1]=ic; fac.ind[2]=id;
    sf.insert(fac);
//    cout << fac.ind[0] << " " << fac.ind[1] << " " << fac.ind[2] << endl;
    }

    {
    Fac fac; fac.reg=tet.reg; 
    fac.ind[0]=ia; fac.ind[1]=id; fac.ind[2]=ic;
    sf.insert(fac);
//    cout << fac.ind[0] << " " << fac.ind[1] << " " << fac.ind[2] << endl;
    }

    {
    Fac fac; fac.reg=tet.reg; 
    fac.ind[0]=ia; fac.ind[1]=ib; fac.ind[2]=id;
    sf.insert(fac); 
//    cout << fac.ind[0] << " " << fac.ind[1] << " " << fac.ind[2] << endl;
    }
}

/*
cout << "===============================================================" << endl;   

  for (set<Fac, less_than>::const_iterator it=sf.begin(); it!=sf.end(); ++it) {
      Fac f1=*it;
      for (int i=0; i<Fac::N; i++) 
          cout << f1.ind[i] << " ";     
      cout << endl;
      }

cout << "===============================================================" << endl;   
*/

// calcul des normales aux faces
int done = 0;
for (int f=0; f<FAC; f++){
    int progress = 100*double(f)/FAC;
    if (progress>done && !(progress%5)) {
        IF_VERBOSE(fem) cout << progress << "--"; fflush(NULL);
        done = progress;
    }

    Fac &fac = fem.fac[f];
    fac.Ms = 0.;
    p = make_pair("Js", fac.reg);
    double Js = param[p];
    if (Js<0.) continue;  // elimination des facettes a Js<0

    int i0 = fac.ind[0],  i1 = fac.ind[1],  i2 = fac.ind[2];

    set< Fac, less_than >::iterator it=sf.end();
    for (int perm=0; perm<2; perm++) {
        for (int nrot=0; nrot<3; nrot++) {
            Fac fc;

            fc.ind[(0+nrot)%3]=i0; fc.ind[(1+nrot)%3]=i1; fc.ind[(2+nrot)%3]=i2;
            it=sf.find(fc);
            if (it!=sf.end()) break;
        }
      
        if (it!=sf.end()) { // found
           Fac fc = *it;
           int i0=fac.ind[0], i1=fac.ind[1], i2=fac.ind[2];
//           cout << "fac " << i0 << " " << i1 << " " << i2 <<endl;
           i0=fc.ind[0], i1=fc.ind[1], i2=fc.ind[2];
//           cout << "fc  " << i0 << " " << i1 << " " << i2 <<endl;

           double x0 = fem.node[i0].x,  y0 = fem.node[i0].y,  z0 = fem.node[i0].z;
           double x1 = fem.node[i1].x,  y1 = fem.node[i1].y,  z1 = fem.node[i1].z;
           double x2 = fem.node[i2].x,  y2 = fem.node[i2].y,  z2 = fem.node[i2].z;

           double nx   = (y1-y0)*(z2-z0)-(y2-y0)*(z1-z0);
           double ny   = (z1-z0)*(x2-x0)-(z2-z0)*(x1-x0);
           double nz   = (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);

           p = make_pair("Js", fc.reg);
           double Ms = nu0*param[p];
//           cout << "Ms : " << Ms << " " << fac.Ms << endl;
           if (nx*fac.nx+ny*fac.ny+nz*fac.nz > 0) {
               fac.Ms = fac.Ms + Ms;   // la face trouvee a la meme orientation que la face traitee
           }
           else {
               fac.Ms = fac.Ms - Ms;   // la face trouvee a une orientation opposee
           }

//           double xbar=(x0+x1+x2)/3, ybar=(y0+y1+y2)/3., zbar=(z0+z1+z2)/3.;
//           cout << "f "<< f << " loc : " << xbar << " " << ybar << " " << zbar << endl;
//           cout << "Js "<< fac.Ms*mu0 << " n   : " << nx << " " << ny << " " << nz << endl;
        }
    int tmp=i1; i1=i2; i2=tmp;
    }//fin perm
}
IF_VERBOSE(fem) cout << "100" << endl;
}

void femutil_fac(Fem &fem)
{
const int FAC = fem.FAC;

// calcul des surfaces
double surftot = 0.;
for (int f=0; f<FAC; f++){
    Fac &fac = fem.fac[f];
    int i0,i1,i2;
    i0 = fac.ind[0];    i1 = fac.ind[1];    i2 = fac.ind[2];
    
    double x0,y0,z0, x1,y1,z1, x2,y2,z2;
    x0 = fem.node[i0].x;   y0 = fem.node[i0].y;   z0 = fem.node[i0].z;
    x1 = fem.node[i1].x;   y1 = fem.node[i1].y;   z1 = fem.node[i1].z;
    x2 = fem.node[i2].x;   y2 = fem.node[i2].y;   z2 = fem.node[i2].z;
     
    double vecx,vecy,vecz,norm,surf;
    vecx = (y1-y0)*(z2-z0)-(y2-y0)*(z1-z0);
    vecy = (z1-z0)*(x2-x0)-(z2-z0)*(x1-x0);
    vecz = (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);
    norm = sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
    fac.nx = vecx/norm,  fac.ny = vecy/norm,  fac.nz = vecz/norm;
    surf = 0.5* norm;
    fac.surf = surf;
    surftot+= surf;
    }
fem.surf = surftot;
}

void femutil(Fem &fem)
{
femutil_node(fem);
femutil_tet(fem);
femutil_fac(fem);
femutil_facMs(fem);
cout << "surface  : " << fem.surf << endl;
cout << "volume   : " << fem.vol << endl;
}

