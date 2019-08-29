#include "fem.h"

void Fem::energy(Settings &settings)
{
Etot = 0.0;
double _E[5] = {0.0,0.0,0.0,0.0,0.0};
double uz_drift=2.*DW_z/l.z()*DW_dir;

std::for_each(tet.begin(),tet.end(),[this,settings,&_E,uz_drift](Tetra::Tet const& te) 
    {
    Tetra::prm const& param = settings.paramTetra[te.idxPrm];
    te.energy(param,_E,Hext,uz_drift);    
    }
);

std::for_each(fac.begin(),fac.end(),[settings,&_E](Facette::Fac const& fa)
    {
    Facette::prm const& param = settings.paramFacette[fa.idxPrm];    
    fa.energy(param,_E);
    }
);

for (int e=0; e<4; e++)
    { E[e] = _E[e]; Etot += _E[e]; }

evol = Etot-Etot0;
}
