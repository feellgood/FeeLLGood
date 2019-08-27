#include "fem.h"
#include <set>
#include <algorithm>
#include <numeric>

void Fem::femutil(Settings &settings)
{
const long nb_nod = node.size();
    
pts= annAllocPts(nb_nod, 3);
kdtree = new ANNkd_tree(pts, nb_nod, 3);
if (!kdtree) SYSTEM_ERROR;

int i=0;
std::for_each(node.begin(),node.end(),
[this,&i](Nodes::Node const& n) 
	{ this->pts[i][0] = n.p.x();this->pts[i][1] = n.p.y();this->pts[i][2] = n.p.z();i++; } 
);// end for_each

const auto minX = std::min_element(node.begin(),node.end(),[](Nodes::Node const& n1,Nodes::Node const& n2) {return (n1.p.x()<n2.p.x()); } );
double xmin = minX->p.x(); 
const auto maxX = std::max_element(node.begin(),node.end(),[](Nodes::Node const& n1,Nodes::Node const& n2) {return (n1.p.x()<n2.p.x()); } );
double xmax = maxX->p.x();

const auto minY = std::min_element(node.begin(),node.end(),[](Nodes::Node const& n1,Nodes::Node const& n2) {return (n1.p.y()<n2.p.y()); } );
double ymin = minY->p.y();
const auto maxY = std::max_element(node.begin(),node.end(),[](Nodes::Node const& n1,Nodes::Node const& n2) {return (n1.p.y()<n2.p.y()); } );
double ymax = maxY->p.y();

const auto minZ = std::min_element(node.begin(),node.end(),[](Nodes::Node const& n1,Nodes::Node const& n2) {return (n1.p.z()<n2.p.z()); } );
double zmin = minZ->p.z();
const auto maxZ = std::max_element(node.begin(),node.end(),[](Nodes::Node const& n1,Nodes::Node const& n2) {return (n1.p.z()<n2.p.z()); } );
double zmax = maxZ->p.z();

// calcul du diametre et du centrage
l = Pt::pt3D(xmax - xmin,ymax - ymin,zmax - zmin);
diam = l.maxLength();
c = Pt::pt3D(0.5*(xmax + xmin),0.5*(ymax + ymin),0.5*(zmax + zmin));

//ici on affecte refNode des tetraedres
std::for_each(tet.begin(),tet.end(),[this](Tetra::Tet &te) {te.setRefNode( &(this->node) );} );

//ici on calcule les volumes elementaires et on reoriente si besoin par une permutation d'indice
std::for_each(tet.begin(),tet.end(),[](Tetra::Tet &te) {te.calc_vol();} );

//volume total
vol = std::accumulate(tet.begin(),tet.end(),0.0,[](double x,Tetra::Tet &te){return x+te.vol;} );

//ici on affecte refNode des facettes
std::for_each(fac.begin(),fac.end(),[this](Facette::Fac &fa) {fa.setRefNode( &(this->node) );} );

//ici on calcule les normales et surfaces elementaires
std::for_each(fac.begin(),fac.end(),[](Facette::Fac &fa) {fa.calc_surf();} );

//ici la somme des surfaces elementaires = surface totale
surf = std::accumulate(fac.begin(),fac.end(),0.0,[](double x,Facette::Fac &fa){return x+fa.surf;} );

//on construit un set de facettes à partir des tetraèdres
std::set<Facette::Fac, Facette::less_than> sf;

std::for_each(tet.begin(),tet.end(),
[&sf](Tetra::Tet const& te)
	{
	int ia=te.ind[0];int ib=te.ind[1];int ic=te.ind[2];int id=te.ind[3];
	sf.insert( Facette::Fac(te.reg,te.idxPrm,ia,ic,ib) );
    sf.insert( Facette::Fac(te.reg,te.idxPrm,ib,ic,id) );
    sf.insert( Facette::Fac(te.reg,te.idxPrm,ia,id,ic) );
    sf.insert( Facette::Fac(te.reg,te.idxPrm,ia,ib,id) ); 
	});//fin for_each

std::for_each(fac.begin(),fac.end(),
[this,&settings,&sf,nb_nod](Facette::Fac &fa)
    {
    fa.Ms = 0.;
    double Js = settings.paramFacette[fa.idxPrm].Js;
    if (Js >= 0.)  // we ignore faces with Js<0
        {
        int i0 = fa.ind[0],  i1 = fa.ind[1],  i2 = fa.ind[2];

        std::set< Facette::Fac, Facette::less_than >::iterator it=sf.end();
        for (int perm=0; perm<2; perm++) 
            {
            for (int nrot=0; nrot<3; nrot++)
                {
                Facette::Fac fc(nb_nod);
                fc.ind[(0+nrot)%3]=i0; fc.ind[(1+nrot)%3]=i1; fc.ind[(2+nrot)%3]=i2;
                it=sf.find(fc);
                if (it!=sf.end()) break;
                }
      
            if (it!=sf.end()) 
                { // found
                Pt::pt3D p0 = node[ it->ind[0] ].p; 
                Pt::pt3D p1 = node[ it->ind[1] ].p;
                Pt::pt3D p2 = node[ it->ind[2] ].p;
                Pt::pt3D n = (p1-p0)*(p2-p0);
	
                //double Ms = nu0*param[ make_pair("Js", it->reg) ];
                double Ms = nu0*settings.paramTetra[it->idxPrm].J;
                if (Pt::pScal(n,fa.n) > 0) { fa.Ms += Ms; }  // la face trouvee a la meme orientation que la face traitee
                else { fa.Ms -= Ms; }  // la face trouvee a une orientation opposee
                }
            int tmp=i1; i1=i2; i2=tmp;
            }//end perm
        }
    }); //end for_each
}

