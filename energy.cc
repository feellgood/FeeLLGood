#include "fem.h"

/*index convention : 0-exchange 1-anisotropy 2-demagnetizing 3-applied*/

void Fem::energy(Settings &settings)
{
Etot = 0.0;
double _E[4] = {0.0,0.0,0.0,0.0};
double uz_drift=2.*DW_z/l.z()*DW_dir;

std::for_each(tet.begin(),tet.end(),[this,settings,&_E,uz_drift](Tetra::Tet const& te) 
    {
    Tetra::prm const& param = settings.paramTetra[te.idxPrm];
    te.energy(param,_E,Hext,uz_drift);    
    }
);

double contribDemag = 0.0;
double contribAnisotropy = 0.0;

std::for_each(fac.begin(),fac.end(),[settings,&contribDemag,&contribAnisotropy](Facette::Fac const& fa)
    {
    Facette::prm const& param = settings.paramFacette[fa.idxPrm];    
    double phi[Facette::NPI];
    double u[DIM][Facette::NPI];

    fa.interpolation(Nodes::get_u,u);
    fa.interpolation(Nodes::get_phi,phi);
    
    if(param.Ks != 0.0)
        { contribAnisotropy += fa.anisotropyEnergy(param,u); }
    contribDemag += fa.demagEnergy(u,phi);
    }
);

for (int e=0; e<4; e++) { E[e] = _E[e]; }
E[1] += contribAnisotropy;
E[2] += contribDemag;
Etot = E[0] + E[1] + E[2] + E[3];
evol = Etot-Etot0;
}
