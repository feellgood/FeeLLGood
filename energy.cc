#include "fem.h"

void Fem::energy(Settings &settings)
{
    const int FAC = fac.size();
    const int TET = tet.size();
Etot = 0.0;
double _E[5] = {0.0,0.0,0.0,0.0,0.0};

double uz_drift=2.*DW_z/l.z()*DW_dir;

/* Contribution des tetraedres */
for (int i_t=0; i_t<TET; i_t++)
    {
	Tetra::Tet &te = tet[i_t];
    Tetra::prm & param = settings.paramTetra[te.idxPrm];
    te.energy(param,_E,Hext,uz_drift);
    }    

/* Contribution des facettes triangulaires a l'energie d'anisotropie */

for (int i_t=0; i_t<FAC; i_t++) {
	Facette::Fac &fa = fac[i_t];
    Facette::prm & param = settings.paramFacette[fa.idxPrm];
    fa.energy(param,_E);
    }
	
for (int e=0; e<4; e++){
	E[e] = _E[e];
	Etot+= _E[e];
    }

evol = Etot-Etot0;
phy = 0;
}
