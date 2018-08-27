#include "fem.h"
#include <set>

using namespace std;


void Fem::femutil_node(void)
{
pts= annAllocPts(NOD, 3);

// calcul du diametre et du centrage
double xmin, xmax, ymin, ymax, zmin, zmax;
xmin = ymin = zmin = +HUGE;
xmax = ymax = zmax = -HUGE;

for (int i=0; i<NOD; i++){
    double xi,yi,zi;
    xi = node[i].p.x();      yi = node[i].p.y();    zi = node[i].p.z();
    pts[i][0]=xi;	     pts[i][1]=yi;	    pts[i][2]=zi;
    if (xi<xmin) xmin=xi;    if (xi>xmax) xmax=xi;
    if (yi<ymin) ymin=yi;    if (yi>ymax) ymax=yi;
    if (zi<zmin) zmin=zi;    if (zi>zmax) zmax=zi;
    }

// allocation de l'arbre de recherche
kdtree = new ANNkd_tree(pts, NOD, 3);
if (!kdtree) SYSTEM_ERROR;

l = Pt::pt3D(xmax-xmin,ymax-ymin,zmax-zmin);//lx=xmax-xmin;ly=ymax-ymin;lz=zmax-zmin;

diam = l.x();
if (diam<l.y()) diam=l.y();
if (diam<l.z()) diam=l.z();

c = Pt::pt3D(0.5*(xmax+xmin),0.5*(ymax+ymin),0.5*(ymax+ymin));
//cx = 0.5*(xmax+xmin);cy = 0.5*(ymax+ymin);cz = 0.5*(zmax+zmin);

as[0] = l.x()/diam;
as[1] = l.y()/diam;
as[2] = l.z()/diam;
}


void Fem::femutil_tet(void)
{
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

for (int i_t=0; i_t<TET; i_t++){
   Tetra::Tet &te = tet[i_t];
   int i0,i1,i2,i3;
   i0=te.ind[0];   i1=te.ind[1];   i2=te.ind[2];   i3=te.ind[3];
   
//   double x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3;
//   x0 = node[i0].p.x();   y0 = node[i0].p.y();   z0 = node[i0].p.z();
//   x1 = node[i1].p.x();   y1 = node[i1].p.y();   z1 = node[i1].p.z();
//   x2 = node[i2].p.x();   y2 = node[i2].p.y();   z2 = node[i2].p.z();
//   x3 = node[i3].p.x();   y3 = node[i3].p.y();   z3 = node[i3].p.z();

//   double vecx,vecy,vecz;
//   vecx = (y1-y0)*(z2-z0)-(y2-y0)*(z1-z0);
//   vecy = (z1-z0)*(x2-x0)-(z2-z0)*(x1-x0);
//   vecz = (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);

	Pt::pt3D p0 = node[i0].p;
	Pt::pt3D p1 = node[i1].p;
	Pt::pt3D p2 = node[i2].p;
	Pt::pt3D p3 = node[i3].p;
	Pt::pt3D vec = (p1-p0)*(p2-p0);

   double i_vol  = 1./6.* pScal(vec,p3-p0);//(vecx*(x3-x0) + vecy*(y3-y0) + vecz*(z3-z0));
   if (i_vol<0.) {
      te.ind[3]=i2; te.ind[2]=i3;
      i_vol=-i_vol;
      IF_VERBOSE() cout << "ill-oriented tetrahedron: " << i_t << " now corrected!"<< endl;
      }
   te.vol = i_vol;
   voltot+= i_vol;}
vol = voltot;
}

void Fem::femutil_facMs(Settings &settings /**< [in] */)
{
pair <string,int> p;
map <pair<string,int>,double> &param = settings.param;

// decomposition des tetraedres en elements de surface
set<Facette::Fac, Facette::less_than> sf;
IF_VERBOSE() cout << "Nb de Tet " << TET << endl;
for (int i_t=0; i_t<TET; i_t++){
    Tetra::Tet &te = tet[i_t];
    int ia,ib,ic,id;
    ia=te.ind[0];  ib=te.ind[1];  ic=te.ind[2];  id=te.ind[3];
//    cout << "tet " << t <<"/"<<TET<< endl;
    {
    Facette::Fac f; f.reg=te.reg; 
    f.ind[0]=ia; f.ind[1]=ic; f.ind[2]=ib;
    sf.insert(f);
//    cout << fac.ind[0] << " " << fac.ind[1] << " " << fac.ind[2] << endl;
    }
 
    {
    Facette::Fac f; f.reg=te.reg;
    f.ind[0]=ib; f.ind[1]=ic; f.ind[2]=id;
    sf.insert(f);
//    cout << fac.ind[0] << " " << fac.ind[1] << " " << fac.ind[2] << endl;
    }

    {
    Facette::Fac f; f.reg=te.reg; 
    f.ind[0]=ia; f.ind[1]=id; f.ind[2]=ic;
    sf.insert(f);
//    cout << fac.ind[0] << " " << fac.ind[1] << " " << fac.ind[2] << endl;
    }

    {
    Facette::Fac f; f.reg=te.reg; 
    f.ind[0]=ia; f.ind[1]=ib; f.ind[2]=id;
    sf.insert(f); 
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
for (int i_f=0; i_f<FAC; i_f++){
    int progress = 100*double(i_f)/FAC;
    if (progress>done && !(progress%5)) {
        IF_VERBOSE() cout << progress << "--"; fflush(NULL);
        done = progress;
    }

    Facette::Fac &fa = fac[i_f];
    fa.Ms = 0.;
    p = make_pair("Js", fa.reg);
    double Js = param[p];
    if (Js<0.) continue;  // elimination des facettes a Js<0

    int i0 = fa.ind[0],  i1 = fa.ind[1],  i2 = fa.ind[2];

    set< Facette::Fac, Facette::less_than >::iterator it=sf.end();
    for (int perm=0; perm<2; perm++) {
        for (int nrot=0; nrot<3; nrot++) {
            Facette::Fac fc;

            fc.ind[(0+nrot)%3]=i0; fc.ind[(1+nrot)%3]=i1; fc.ind[(2+nrot)%3]=i2;
            it=sf.find(fc);
            if (it!=sf.end()) break;
        }
      
        if (it!=sf.end()) { // found
           Facette::Fac fc = *it;
           int i0=fa.ind[0], i1=fa.ind[1], i2=fa.ind[2];
//           cout << "fac " << i0 << " " << i1 << " " << i2 <<endl;
           i0=fc.ind[0], i1=fc.ind[1], i2=fc.ind[2];
//           cout << "fc  " << i0 << " " << i1 << " " << i2 <<endl;

	Pt::pt3D p0 = node[i0].p; 
	Pt::pt3D p1 = node[i1].p;
	Pt::pt3D p2 = node[i2].p;
        Pt::pt3D n = (p1-p0)*(p2-p0);
	
	//double x0 = node[i0].x,  y0 = node[i0].y,  z0 = node[i0].z;
        //double x1 = node[i1].x,  y1 = node[i1].y,  z1 = node[i1].z;
        //double x2 = node[i2].x,  y2 = node[i2].y,  z2 = node[i2].z;

        //double nx   = (y1-y0)*(z2-z0)-(y2-y0)*(z1-z0);
        //double ny   = (z1-z0)*(x2-x0)-(z2-z0)*(x1-x0);
        //double nz   = (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);

           p = make_pair("Js", fc.reg);
           double Ms = nu0*param[p];
//           cout << "Ms : " << Ms << " " << fac.Ms << endl;
           //if (n.x()*fa.n.x()+n.y()*fa.n.y()+n.z()*fa.n.z() > 0) 
	if (Pt::pScal(n,fa.n) > 0) {
         	fa.Ms = fa.Ms + Ms;   // la face trouvee a la meme orientation que la face traitee
           	}
           else {
               fa.Ms = fa.Ms - Ms;   // la face trouvee a une orientation opposee
           }

//           double xbar=(x0+x1+x2)/3, ybar=(y0+y1+y2)/3., zbar=(z0+z1+z2)/3.;
//           cout << "f "<< f << " loc : " << xbar << " " << ybar << " " << zbar << endl;
//           cout << "Js "<< fac.Ms*mu0 << " n   : " << nx << " " << ny << " " << nz << endl;
        }
    int tmp=i1; i1=i2; i2=tmp;
    }//fin perm
}
IF_VERBOSE() cout << "100" << endl;
}

void Fem::femutil_fac(void)
{
// calcul des surfaces
double surftot = 0.;
for (int f=0; f<FAC; f++){
    Facette::Fac &fa = fac[f];
    int i0,i1,i2;
    i0 = fa.ind[0];    i1 = fa.ind[1];    i2 = fa.ind[2];
    
	    
	//double x0,y0,z0, x1,y1,z1, x2,y2,z2;
    //x0 = node[i0].x;   y0 = node[i0].y;   z0 = node[i0].z;
    //x1 = node[i1].x;   y1 = node[i1].y;   z1 = node[i1].z;
    //x2 = node[i2].x;   y2 = node[i2].y;   z2 = node[i2].z;

	Pt::pt3D p0 = node[i0].p;
     	Pt::pt3D p1 = node[i1].p;
	Pt::pt3D p2 = node[i2].p;

	Pt::pt3D vec = (p1-p0)*(p2-p0);

    //double vecx = (y1-y0)*(z2-z0)-(y2-y0)*(z1-z0);
    //double vecy = (z1-z0)*(x2-x0)-(z2-z0)*(x1-x0);
    //double vecz = (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);
    double norm = vec.norm();//sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
    fa.n = vec/norm;//fa.nx = vec.x()/norm,  fa.ny = vec.y()/norm,  fa.nz = vec.z()/norm;
    fa.surf = 0.5*norm;
    surftot+= fa.surf;
    }
surf = surftot;
}

void Fem::femutil(Settings &settings)
{
femutil_node();
femutil_tet();
femutil_fac();
femutil_facMs(settings);
cout << "surface  : " << surf << std::endl;
cout << "volume   : " << vol << std::endl;
}

