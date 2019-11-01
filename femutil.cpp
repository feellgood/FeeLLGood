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
}

void Fem::geometry(void)
{
double xmin = minNodes(Pt::IDX_X);
double xmax = maxNodes(Pt::IDX_X);

double ymin = minNodes(Pt::IDX_Y);
double ymax = maxNodes(Pt::IDX_Y);

double zmin = minNodes(Pt::IDX_Z);
double zmax = maxNodes(Pt::IDX_Z);

l = Pt::pt3D(xmax - xmin,ymax - ymin,zmax - zmin);
diam = l.maxLength();
c = Pt::pt3D(0.5*(xmax + xmin),0.5*(ymax + ymin),0.5*(zmax + zmin));
vol = std::accumulate(tet.begin(),tet.end(),0.0,[](double x,Tetra::Tet const& te){return x + te.calc_vol();} );
surf = std::accumulate(fac.begin(),fac.end(),0.0,[](double x,Facette::Fac const& fa){return x+fa.surf;} );
}

void Fem::femutil(Settings const& settings)
{
//facettes set with order less_than
std::set<Facette::Fac, Facette::less_than> sf;
    
std::for_each(tet.begin(),tet.end(),[&sf](Tetra::Tet const& te)
	{
	int ia=te.ind[0];int ib=te.ind[1];int ic=te.ind[2];int id=te.ind[3];
    const int reg = te.getRegion();
	sf.insert( Facette::Fac(reg,te.idxPrm,ia,ic,ib) );
    sf.insert( Facette::Fac(reg,te.idxPrm,ib,ic,id) );
    sf.insert( Facette::Fac(reg,te.idxPrm,ia,id,ic) );
    sf.insert( Facette::Fac(reg,te.idxPrm,ia,ib,id) ); 
	});//end for_each

std::for_each(fac.begin(),fac.end(),[this,&settings,&sf](Facette::Fac &fa)
    {
    fa.Ms = 0.;
    if ( settings.paramFacette[fa.idxPrm].Js >= 0.)  // we ignore faces with Js<0
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
                
                //fa.Ms will have the magnitude of first arg of copysign, with the sign of second arg
                fa.Ms = std::copysign(nu0*settings.paramTetra[it->idxPrm].J , Pt::pTriple( p1-p0 , p2-p0 ,fa.n) );
                }
            std::swap(i1,i2);
            }//end perm
        }
    }); //end for_each
}


