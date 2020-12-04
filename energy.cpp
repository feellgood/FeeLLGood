#include "fem.h"


void Fem::energy(double const t,Settings & settings)
{
double uz_drift=2.*DW_z/l.z()*DW_dir;

zeroEnergy();

const Pt::pt3D Hext = settings.getValue(t);

std::for_each(tet.begin(),tet.end(),[this,&Hext,&settings,uz_drift](Tetra::Tet const& te) 
    {
    Tetra::prm const& param = settings.paramTetra[te.idxPrm];
    double u[Pt::DIM][Tetra::NPI],dudx[Pt::DIM][Tetra::NPI], dudy[Pt::DIM][Tetra::NPI], dudz[Pt::DIM][Tetra::NPI];
    double phi[Tetra::NPI];

    te.interpolation(Nodes::get_u,u,dudx,dudy,dudz);
    te.interpolation(Nodes::get_phi,phi);
    
    E_exch += te.exchangeEnergy(param,dudx,dudy,dudz);
    
    E_demag += te.demagEnergy(param,dudx,dudy,dudz,phi);
    
    if((param.K != 0.0)||(param.K3 != 0.0))
        { E_aniso += te.anisotropyEnergy(param,u); }
    
    E_zeeman += te.zeemanEnergy(param,uz_drift,Hext,u);
    }
);

std::for_each(fac.begin(),fac.end(),[this,&settings](Facette::Fac const& fa)
    {
    Facette::prm const& param = settings.paramFacette[fa.idxPrm];    
    double phi[Facette::NPI];
    Pt::pt3D u[Facette::NPI];

    fa.interpolation(Nodes::get_u,u);
    fa.interpolation(Nodes::get_phi,phi);
    
    if(param.Ks != 0.0)
        { E_aniso += fa.anisotropyEnergy(param,u); }
    E_demag += fa.demagEnergy(u,phi);
    }
);

calc_Etot();

if (settings.verbose && (Etot > Etot0))
    { std::cout << "Warning energy increases! : Etot= " << Etot << " ; Etot0 = "<< Etot0 << std::endl; }
}
