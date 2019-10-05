#include "fem.h"
#include <set>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <string>

void Fem::infos(void) const
{
std::cout << "This is feeLLGood SHA1= " + std::string(SHAnumber) << std::endl;
std::cout << "\t diam bounding box ="<< diam << std::endl;
std::cout << "\t nodes\t\t\t" << node.size() << std::endl;
std::cout << "\t faces\t\t\t" << fac.size() << std::endl;
std::cout << "\t tetraedrons\t\t" << tet.size() << std::endl;
std::cout << "\t Total surface\t\t"  << surf << std::endl;
std::cout << "\t Total volume\t\t\t" << vol << std::endl;
std::cout << "\t fmm_normalizer = " << fmm_normalizer << std::endl;
}

void Fem::femutil(Settings const& settings)
{
double xmin = minNodes(Pt::IDX_X);
double xmax = maxNodes(Pt::IDX_X);

double ymin = minNodes(Pt::IDX_Y);
double ymax = maxNodes(Pt::IDX_Y);

double zmin = minNodes(Pt::IDX_Z);
double zmax = maxNodes(Pt::IDX_Z);

// calcul du diametre et du centrage
l = Pt::pt3D(xmax - xmin,ymax - ymin,zmax - zmin);
diam = l.maxLength();
fmm_normalizer = 1./(2.*diam);
c = Pt::pt3D(0.5*(xmax + xmin),0.5*(ymax + ymin),0.5*(zmax + zmin));

//volume total
vol = std::accumulate(tet.begin(),tet.end(),0.0,[](double x,Tetra::Tet &te){return x+te.vol;} );

// somme des surfaces elementaires = surface totale
surf = std::accumulate(fac.begin(),fac.end(),0.0,[](double x,Facette::Fac &fa){return x+fa.surf;} );

//on construit un set de facettes avec l'ordre less_than
std::set<Facette::Fac, Facette::less_than> sf;

std::for_each(tet.begin(),tet.end(),[&sf](Tetra::Tet const& te)
	{
	int ia=te.ind[0];int ib=te.ind[1];int ic=te.ind[2];int id=te.ind[3];
    const int reg = te.getRegion();
	sf.insert( Facette::Fac(reg,te.idxPrm,ia,ic,ib) );
    sf.insert( Facette::Fac(reg,te.idxPrm,ib,ic,id) );
    sf.insert( Facette::Fac(reg,te.idxPrm,ia,id,ic) );
    sf.insert( Facette::Fac(reg,te.idxPrm,ia,ib,id) ); 
	});//fin for_each

std::for_each(fac.begin(),fac.end(),[this,&settings,&sf](Facette::Fac &fa)
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
                Facette::Fac fc(node.size());
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


