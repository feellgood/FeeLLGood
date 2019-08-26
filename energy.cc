#include "fem.h"

void Fem::energy(Settings &settings)
{
    const int FAC = fac.size();
    const int TET = tet.size();
Etot = 0.0;
double _E[5] = {0.0,0.0,0.0,0.0,0.0};

double uz_drift=2.*DW_z/l.z()*DW_dir;

/* Contribution des tetraedres */
for (int i_t=0; i_t<TET; i_t++) {
	
    Tetra::Tet &te = tet[i_t];
    Tetra::prm & param = settings.paramTetra[te.idxPrm];
    te.energy(param,_E,Hext,uz_drift);
    /*
    double Ae = param.A;
    double Js = param.J;
    double Ms = nu0 * Js;

    double K = param.K;
    double K3 = param.K3;
    
    double uk00 = param.uk[0][0];
    double uk01 = param.uk[0][1];
    double uk02 = param.uk[0][2];
    double uk10 = param.uk[1][0];
    double uk11 = param.uk[1][1];
    double uk12 = param.uk[1][2];
    double uk20 = param.uk[2][0];
    double uk21 = param.uk[2][1];
    double uk22 = param.uk[2][2];   
   */
   /*-------------------- INTERPOLATION --------------------*/
   /* 
   double u_nod[3][Tetra::N], u[3][Tetra::NPI];
    double dudx[3][Tetra::NPI], dudy[3][Tetra::NPI], dudz[3][Tetra::NPI];
    double q[Tetra::NPI],  phi[Tetra::NPI];
    double phi_nod[Tetra::N], negphi_nod[Tetra::N], Hdx[Tetra::NPI], Hdy[Tetra::NPI], Hdz[Tetra::NPI];

    for (int i=0; i<Tetra::N; i++)
        {
        int i_= te.ind[i];
        Node &n = node[i_];
        u_nod[Pt::IDX_X][i] = n.u.x(); u_nod[Pt::IDX_Y][i] = n.u.y(); u_nod[Pt::IDX_Z][i] = n.u.z();
        phi_nod[i] =  n.phi;
        negphi_nod[i] = -n.phi;
        }

	tiny::transposed_mult<double, Tetra::N, Tetra::NPI> (phi_nod, Tetra::a, phi);

	tiny::mult<double, 3, Tetra::N, Tetra::NPI> (u_nod, Tetra::a, u);
	tiny::mult<double, 3, Tetra::N, Tetra::NPI> (u_nod, te.dadx, dudx);
	tiny::mult<double, 3, Tetra::N, Tetra::NPI> (u_nod, te.dady, dudy);
	tiny::mult<double, 3, Tetra::N, Tetra::NPI> (u_nod, te.dadz, dudz);

    for (int npi=0; npi<Tetra::NPI; npi++){
        //double div_u = dudx[0][npi] + dudy[1][npi] + dudz[2][npi];
        q[npi] = -Ms*(dudx[0][npi] + dudy[1][npi] + dudz[2][npi]);
        }

	tiny::transposed_mult<double, Tetra::N, Tetra::NPI> (negphi_nod, te.dadx, Hdx);
	tiny::transposed_mult<double, Tetra::N, Tetra::NPI> (negphi_nod, te.dady, Hdy);
	tiny::transposed_mult<double, Tetra::N, Tetra::NPI> (negphi_nod, te.dadz, Hdz);

    double dens[5][Tetra::NPI];
    double Eelem[5];

    for (int npi=0; npi<Tetra::NPI; npi++) {

        // cosinus directeurs
        double al0=uk00*u[0][npi] + uk01*u[1][npi] + uk02*u[2][npi];
        double al1=uk10*u[0][npi] + uk11*u[1][npi] + uk12*u[2][npi];
        double al2=uk20*u[0][npi] + uk21*u[1][npi] + uk22*u[2][npi];

//        cout << dudx[0][npi] << "\t" << dudx[1][npi] << "\t" << dudx[2][npi] << endl;
        dens[0][npi] = Ae*(sq( dudx[0][npi] ) + sq( dudy[0][npi] ) + sq( dudz[0][npi] )+
                           sq( dudx[1][npi] ) + sq( dudy[1][npi] ) + sq( dudz[1][npi] )+
                           sq( dudx[2][npi] ) + sq( dudy[2][npi] ) + sq( dudz[2][npi] ) );

        dens[1][npi] = -K  * sq(al0);      // uniaxe
        dens[1][npi]+= +K3 * (sq(al0) * sq(al1) + sq(al1) * sq(al2) + sq(al2) * sq(al0)); //cubique
 
 //       dens[2][npi] = -0.5*Js* (u[0][npi]*Hdx[npi] + u[1][npi]*Hdy[npi] + u[2][npi]*Hdz[npi]);

        dens[2][npi] = 0.5*mu0*q[npi]*phi[npi];

        dens[3][npi] = -Js*(u[0][npi]*Hext[0] + u[1][npi]*Hext[1] + u[2][npi]*Hext[2] + uz_drift*Hext[2]);

        dens[4][npi] = 0.;
        }
    tiny::mult<double, 5, Tetra::NPI> (dens, te.weight, Eelem);
    _E[0] += Eelem[0];
    _E[1] += Eelem[1];
    _E[2] += Eelem[2];
    _E[3] += Eelem[3];
    _E[4] += Eelem[4];
    
    */
}

/* Contribution des facettes triangulaires a l'energie d'anisotropie */

	
for (int i_t=0; i_t<FAC; i_t++) {
	Facette::Fac &fa = fac[i_t];
    double Ms = fa.Ms;
    Pt::pt3D n = fa.n;

    Facette::prm & param = settings.paramFacette[fa.idxPrm];
    double K = param.Ks;
	double uk00 = param.uk[0];
	double uk01 = param.uk[1];
	double uk02 = param.uk[2];    
		
	/*-------------------- INTERPOLATION --------------------*/
	double u_nod[3][Facette::N], u[3][Facette::NPI];
    double phi_nod[Facette::N];
    double q[Facette::NPI],  phi[Facette::NPI];
		
	for (int i=0; i<Facette::N; i++)
		{
		int i_= fa.ind[i];
		Node &n = node[i_];
		u_nod[Pt::IDX_X][i] = n.u.x(); u_nod[Pt::IDX_Y][i] = n.u.y(); u_nod[Pt::IDX_Z][i] = n.u.z();	        
		phi_nod[i] =  n.phi;
        }

	tiny::transposed_mult<double, Facette::N, Facette::NPI> (phi_nod, Facette::a, phi);
	tiny::mult<double, 3, Facette::N, Facette::NPI> (u_nod, Facette::a, u);
	
    for (int npi=0; npi<Facette::NPI; npi++)
        { q[npi] = Ms * (u[0][npi]*n.x() + u[1][npi]*n.y() + u[2][npi]*n.z()); }
	
	/*-------------------------------------------------------*/	    
	double dens[Facette::NPI];
	double Eelem;
		
	for (int npi=0; npi<Facette::NPI; npi++)
        {// cosinus directeurs
		double al0=uk00*u[0][npi] + uk01*u[1][npi] + uk02*u[2][npi];
		dens[npi] = -K  * sq(al0);      // uniaxe
		}
	Eelem=tiny::sp<double, Facette::NPI> (dens, fa.weight);
	_E[1]+=Eelem;

	for (int npi=0; npi<Facette::NPI; npi++) { dens[npi] = 0.5*mu0*q[npi]*phi[npi]; }
	Eelem=tiny::sp<double, Facette::NPI> (dens, fa.weight);
	_E[2]+=Eelem;
    }
	
for (int e=0; e<4; e++){
	E[e] = _E[e];
	Etot+= _E[e];
    }

evol = Etot-Etot0;
phy = 0;
}
