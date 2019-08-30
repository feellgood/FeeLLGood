#include "fem.h"

/*index convention : 0-exchange 1-anisotropy 2-demagnetizing 3-applied*/

void Fem::energy(Settings &settings)
{
Etot = 0.0;
double uz_drift=2.*DW_z/l.z()*DW_dir;

double exchange = 0.0;
double demag = 0.0;
double anisotropy = 0.0;
double zeeman = 0.0;

std::for_each(tet.begin(),tet.end(),[this,settings,&exchange,&demag,&anisotropy,&zeeman,uz_drift](Tetra::Tet const& te) 
    {
    Tetra::prm const& param = settings.paramTetra[te.idxPrm];
    double u[DIM][Tetra::NPI],dudx[DIM][Tetra::NPI], dudy[DIM][Tetra::NPI], dudz[DIM][Tetra::NPI];
    double phi[Tetra::NPI];

    te.interpolation(Nodes::get_u,u,dudx,dudy,dudz);
    te.interpolation(Nodes::get_phi,phi);
    
    exchange += te.exchangeEnergy(param,dudx,dudy,dudz);
    demag += te.demagEnergy(param,dudx,dudy,dudz,phi);
    if((param.K != 0.0)||(param.K3 != 0.0))
        { anisotropy += te.anisotropyEnergy(param,u); }
    zeeman += te.zeemanEnergy(param,uz_drift,Hext,u);
    }
);

std::for_each(fac.begin(),fac.end(),[settings,&demag,&anisotropy](Facette::Fac const& fa)
    {
    Facette::prm const& param = settings.paramFacette[fa.idxPrm];    
    double phi[Facette::NPI];
    double u[DIM][Facette::NPI];

    fa.interpolation(Nodes::get_u,u);
    fa.interpolation(Nodes::get_phi,phi);
    
    if(param.Ks != 0.0)
        { anisotropy += fa.anisotropyEnergy(param,u); }
    demag += fa.demagEnergy(u,phi);
    }
);

E[0] = exchange;
E[1] = anisotropy;
E[2] = demag;
E[3] = zeeman;
Etot = E[0] + E[1] + E[2] + E[3];
evol = Etot-Etot0;
}
