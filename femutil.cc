#include "fem.h"
#include <set>
#include <algorithm>
#include <numeric>

using namespace std;


void Fem::femutil_node(void)
{
pts= annAllocPts(NOD, 3);

int i=0;
std::for_each(node.begin(),node.end(),
[this,&i](Node const& n) 
	{ this->pts[i][0] = n.p.x();this->pts[i][1] = n.p.y();this->pts[i][2] = n.p.z();i++; } 
);// end for_each

// allocation de l'arbre de recherche
kdtree = new ANNkd_tree(pts, NOD, 3);
if (!kdtree) SYSTEM_ERROR;

auto extrema = std::minmax_element(node.begin(),node.end(),[](Node const& n1,Node const& n2) {return !(n2.p.x()<n1.p.x())?n1.p.x():n2.p.x(); } );
double xmin = extrema.first->p.x();
double xmax = extrema.second->p.x();

extrema = std::minmax_element(node.begin(),node.end(),[](Node const& n1,Node const& n2) {return !(n2.p.y()<n1.p.y())?n1.p.y():n2.p.y(); } );
double ymin = extrema.first->p.y();
double ymax = extrema.second->p.y();

extrema = std::minmax_element(node.begin(),node.end(),[](Node const& n1,Node const& n2) {return !(n2.p.z()<n1.p.z())?n1.p.z():n2.p.z(); } );
double zmin = extrema.first->p.z();
double zmax = extrema.second->p.z();

// calcul du diametre et du centrage
l = Pt::pt3D(xmax-xmin,ymax-ymin,zmax-zmin);

diam = l.x();
if (diam<l.y()) diam=l.y();
if (diam<l.z()) diam=l.z();

c = Pt::pt3D(0.5*(xmax+xmin),0.5*(ymax+ymin),0.5*(ymax+ymin));
}

void Fem::femutil_facMs(Settings &settings /**< [in] */)
{
pair <string,int> p;
map <pair<string,int>,double> &param = settings.param;

// decomposition des tetraedres en elements de surface
set<Facette::Fac, Facette::less_than> sf;
IF_VERBOSE() cout << "Nb de Tet " << TET << endl;

std::for_each(tet.begin(),tet.end(),
[&sf](Tetra::Tet const& te) //subtil, le capture oblige à passer une ref du set sf
	{
	int ia,ib,ic,id;
    	ia=te.ind[0];  ib=te.ind[1];  ic=te.ind[2];  id=te.ind[3];
	sf.insert( Facette::Fac(te.reg,ia,ic,ib) );
    	sf.insert( Facette::Fac(te.reg,ib,ic,id) );
    	sf.insert( Facette::Fac(te.reg,ia,id,ic) );
    	sf.insert( Facette::Fac(te.reg,ia,ib,id) ); 
	});//fin for_each


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
           int i0 = it->ind[0];
	   int i1 = it->ind[1];
	   int i2 = it->ind[2];

	   Pt::pt3D p0 = node[i0].p; 
	   Pt::pt3D p1 = node[i1].p;
	   Pt::pt3D p2 = node[i2].p;
           Pt::pt3D n = (p1-p0)*(p2-p0);
	
	   p = make_pair("Js", it->reg);
           double Ms = nu0*param[p];
//           cout << "Ms : " << Ms << " " << fac.Ms << endl;
        if (Pt::pScal(n,fa.n) > 0) { fa.Ms += Ms; }  // la face trouvee a la meme orientation que la face traitee
        else { fa.Ms -= Ms; }  // la face trouvee a une orientation opposee
        }
    int tmp=i1; i1=i2; i2=tmp;
    }//fin perm
}
IF_VERBOSE() cout << "100" << endl;
}

void Fem::femutil(Settings &settings)
{
femutil_node();

//ici on calcule les volumes elementaires et on reoriente si besoin par une permutation d'indice
std::for_each(tet.begin(),tet.end(),[this](Tetra::Tet &te) {te.calc_vol(this->node);} );

//volume total
vol = std::accumulate(tet.begin(),tet.end(),0.0,[](double x,Tetra::Tet &te){return x+te.vol;} );

//ici on calcule les normales et surfaces elementaires
std::for_each(fac.begin(),fac.end(),[this](Facette::Fac &fa) {fa.calc_surf(this->node);} );

//ici la somme des surfaces elementaires = surface totale
surf = std::accumulate(fac.begin(),fac.end(),0.0,[](double x,Facette::Fac &fa){return x+fa.surf;} );

femutil_facMs(settings);
cout << "surface  : " << surf << std::endl;
cout << "volume   : " << vol << std::endl;
}

