#include "fem.h"

using namespace std;

void Fem::energy(Settings &settings)
{
pair <string,int> p;
//map <pair<string,int>,double> &param = settings.param;

Etot = 0.0;
double _E[5] = {0.0,0.0,0.0,0.0,0.0};

double uz_drift=2.*DW_z/lz*DW_dir;

/* Contribution des tetraedres */
for (int i_t=0; i_t<TET; i_t++) {
	
    Tetra::Tet &te = tet[i_t];
    const int reg = te.reg;
    p = make_pair("Ae", reg);   double Ae = settings.param[p];
    p = make_pair("Js", reg);   double Js = settings.param[p];  double Ms = nu0 * settings.param[p];

    p=make_pair("Ka",reg);      double K = settings.param[p];		//cout << ", Ka=" << K;
    p=make_pair("Ka3",reg);     double K3 = settings.param[p];   	//cout << ", Ka3=" << K3;

    p=make_pair("a1",reg);      double uk00 = settings.param[p];  	//cout << ", a1=" << k0;
    p=make_pair("a2",reg);      double uk01 = settings.param[p];  	//cout << ", a2=" << k1;
    p=make_pair("a3",reg);      double uk02 = settings.param[p];  	//cout << ", a3=" << k2;
    p=make_pair("b1",reg);      double uk10 = settings.param[p];  	//cout << ", a1=" << k0;
    p=make_pair("b2",reg);      double uk11 = settings.param[p];  	//cout << ", a2=" << k1;
    p=make_pair("b3",reg);      double uk12 = settings.param[p];  	//cout << ", a3=" << k2;
    p=make_pair("c1",reg);      double uk20 = settings.param[p];  	//cout << ", a1=" << k0;
    p=make_pair("c2",reg);      double uk21 = settings.param[p];  	//cout << ", a2=" << k1;
    p=make_pair("c3",reg);      double uk22 = settings.param[p];  	//cout << ", a3=" << k2;	
    //p = make_pair("alpha", reg);double alpha = fem.param[p];    
   
   /*-------------------- INTERPOLATION --------------------*/
    double u_nod[3][Tetra::N], u[3][Tetra::NPI];
    double dudx[3][Tetra::NPI], dudy[3][Tetra::NPI], dudz[3][Tetra::NPI];
    double q[Tetra::NPI],  phi[Tetra::NPI];
    double phi_nod[Tetra::N], negphi_nod[Tetra::N], Hdx[Tetra::NPI], Hdy[Tetra::NPI], Hdz[Tetra::NPI];

    for (int i=0; i<Tetra::N; i++){
        int i_= te.ind[i];
        Node &n = node[i_];
        for (int d=0; d<3; d++) {
            u_nod[d][i] = n.u[d];
            }
           phi_nod[i] =  n.phi;
        negphi_nod[i] = -n.phi;
        }

	tiny::transposed_mult<double, Tetra::N, Tetra::NPI> (phi_nod, te.a, phi);

	tiny::mult<double, 3, Tetra::N, Tetra::NPI> (u_nod, te.a, u);
	tiny::mult<double, 3, Tetra::N, Tetra::NPI> (u_nod, te.dadx, dudx);
	tiny::mult<double, 3, Tetra::N, Tetra::NPI> (u_nod, te.dady, dudy);
	tiny::mult<double, 3, Tetra::N, Tetra::NPI> (u_nod, te.dadz, dudz);

    for (int npi=0; npi<Tetra::NPI; npi++){
        double div_u = dudx[0][npi] + dudy[1][npi] + dudz[2][npi];
        q[npi] = -Ms * div_u;
        }

	tiny::transposed_mult<double, Tetra::N, Tetra::NPI> (negphi_nod, te.dadx, Hdx);
	tiny::transposed_mult<double, Tetra::N, Tetra::NPI> (negphi_nod, te.dady, Hdy);
	tiny::transposed_mult<double, Tetra::N, Tetra::NPI> (negphi_nod, te.dadz, Hdz);

   /*-------------------------------------------------------*/
	    
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
    tiny::add<double, 5> (Eelem, _E);
    }

/* Contribution des triangles a l'energie d'anisotropie */

	
for (int i_t=0; i_t<FAC; i_t++) {
	Facette::Fac &fa = fac[i_t];
const int reg = fa.reg;
    double Ms = fa.Ms;
    double nx,ny,nz;
    nx=fa.nx; ny=fa.ny; nz=fa.nz;

	p=make_pair("Ks",reg);      double K = settings.param[p];	//cout << ", Ks=" << K;
	p=make_pair("a1",reg);      double uk00 = settings.param[p];  	//cout << ", a1=" << k0;
	p=make_pair("a2",reg);      double uk01 = settings.param[p];  	//cout << ", a2=" << k1;
	p=make_pair("a3",reg);      double uk02 = settings.param[p];  	//cout << ", a3=" << k2;    
		
	/*-------------------- INTERPOLATION --------------------*/
	double u_nod[3][Facette::N], u[3][Facette::NPI];
    double phi_nod[Facette::N];
    double q[Facette::NPI],  phi[Facette::NPI];
		
	for (int i=0; i<Facette::N; i++){
		int i_= fa.ind[i];
			Node &n = node[i_];
			for (int d=0; d<3; d++) {
				u_nod[d][i] = n.u[d];
            }
        phi_nod[i] =  n.phi;
        }

	tiny::transposed_mult<double, Facette::N, Facette::NPI> (phi_nod, fa.a, phi);
		
	tiny::mult<double, 3, Facette::N, Facette::NPI> (u_nod, fa.a, u);
	
    for (int npi=0; npi<Facette::NPI; npi++) { q[npi] = Ms * (u[0][npi]*nx + u[1][npi]*ny + u[2][npi]*nz); }
	
	/*-------------------------------------------------------*/	    
	double dens[Facette::NPI];
	double Eelem;
		
	for (int npi=0; npi<Facette::NPI; npi++) {
		// cosinus directeurs
		double al0=uk00*u[0][npi] + uk01*u[1][npi] + uk02*u[2][npi];
		dens[npi] = -K  * sq(al0);      // uniaxe
		}
	Eelem=tiny::sp<double, Facette::NPI> (dens, fa.weight);
	_E[1]+=Eelem;

	for (int npi=0; npi<Facette::NPI; npi++) {
		dens[npi] = 0.5*mu0*q[npi]*phi[npi];
		}
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
