#include "fem.h"
#include "tiny.h"

using namespace std;

void energy(Fem &fem,Settings &settings)
{
const int TET = fem.TET;
const int FAC = fem.FAC;

pair <string,int> p;
map <pair<string,int>,double> &param = settings.param;
triple &Hext=fem.Hext;

fem.Etot = 0.0;
double E[5] = {0.0,0.0,0.0,0.0,0.0};

double uz_drift=2.*fem.DW_z/fem.lz*fem.DW_dir;

/* Contribution des tetraedres */
for (int t=0; t<TET; t++) {
	const int N   = Tet::N;
	const int NPI = Tet::NPI;
	
    Tet &tet = fem.tet[t];
    const int reg = tet.reg;
    p = make_pair("Ae", reg);   double Ae = settings.param[p];
    p = make_pair("Js", reg);   double Js = settings.param[p];  double Ms = nu0 * settings.param[p];

    p=make_pair("Ka",reg);      double K = param[p];		//cout << ", Ka=" << K;
    p=make_pair("Ka3",reg);     double K3 = param[p];   	//cout << ", Ka3=" << K3;

    p=make_pair("a1",reg);      double uk00 = param[p];  	//cout << ", a1=" << k0;
    p=make_pair("a2",reg);      double uk01 = param[p];  	//cout << ", a2=" << k1;
    p=make_pair("a3",reg);      double uk02 = param[p];  	//cout << ", a3=" << k2;
    p=make_pair("b1",reg);      double uk10 = param[p];  	//cout << ", a1=" << k0;
    p=make_pair("b2",reg);      double uk11 = param[p];  	//cout << ", a2=" << k1;
    p=make_pair("b3",reg);      double uk12 = param[p];  	//cout << ", a3=" << k2;
    p=make_pair("c1",reg);      double uk20 = param[p];  	//cout << ", a1=" << k0;
    p=make_pair("c2",reg);      double uk21 = param[p];  	//cout << ", a2=" << k1;
    p=make_pair("c3",reg);      double uk22 = param[p];  	//cout << ", a3=" << k2;
    //p = make_pair("alpha", reg);double alpha = fem.param[p];    
   
   /*-------------------- INTERPOLATION --------------------*/
    double u_nod[3][N], u[3][NPI];
    double dudx[3][NPI], dudy[3][NPI], dudz[3][NPI];
    double q[NPI],  phi[NPI];
    double phi_nod[N], negphi_nod[N], Hdx[NPI], Hdy[NPI], Hdz[NPI];

    for (int i=0; i<N; i++){
        int i_= tet.ind[i];
        Node &node = fem.node[i_];
        for (int d=0; d<3; d++) {
            u_nod[d][i] = node.u[d];
            }
           phi_nod[i] =  node.phi;
        negphi_nod[i] = -node.phi;
        }

	tiny::transposed_mult<double, N, NPI> (phi_nod, tet.a, phi);

	tiny::mult<double, 3, N, NPI> (u_nod, tet.a, u);
	tiny::mult<double, 3, N, NPI> (u_nod, tet.dadx, dudx);
	tiny::mult<double, 3, N, NPI> (u_nod, tet.dady, dudy);
	tiny::mult<double, 3, N, NPI> (u_nod, tet.dadz, dudz);

    for (int npi=0; npi<NPI; npi++){
        double div_u = dudx[0][npi] + dudy[1][npi] + dudz[2][npi];
        q[npi] = -Ms * div_u;
        }

	tiny::transposed_mult<double, N, NPI> (negphi_nod, tet.dadx, Hdx);
	tiny::transposed_mult<double, N, NPI> (negphi_nod, tet.dady, Hdy);
	tiny::transposed_mult<double, N, NPI> (negphi_nod, tet.dadz, Hdz);

   /*-------------------------------------------------------*/
	    
    double dens[5][NPI];
    double Eelem[5];

    for (int npi=0; npi<NPI; npi++) {

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
    tiny::mult<double, 5, NPI> (dens, tet.weight, Eelem);
    tiny::add<double, 5> (Eelem, E);
    }

/* Contribution des triangles a l'energie d'anisotropie */

	
for (int t=0; t<FAC; t++) {
	//const int N   = Fac::N;
	//const int NPI = Fac::NPI;
	
	Facette::Fac &fac = fem.fac[t];
	const int reg = fac.reg;
    double Ms = fac.Ms;
    double nx,ny,nz;
    nx=fac.nx; ny=fac.ny; nz=fac.nz;

	p=make_pair("Ks",reg);      double K = param[p];	//cout << ", Ks=" << K;
	p=make_pair("a1",reg);      double uk00 = param[p];  	//cout << ", a1=" << k0;
	p=make_pair("a2",reg);      double uk01 = param[p];  	//cout << ", a2=" << k1;
	p=make_pair("a3",reg);      double uk02 = param[p];  	//cout << ", a3=" << k2;    
		
	/*-------------------- INTERPOLATION --------------------*/
	double u_nod[3][Facette::N], u[3][Facette::NPI];
    double phi_nod[Facette::N];
    double q[Facette::NPI],  phi[Facette::NPI];
		
	for (int i=0; i<Facette::N; i++){
		int i_= fac.ind[i];
			Node &node = fem.node[i_];
			for (int d=0; d<3; d++) {
				u_nod[d][i] = node.u[d];
            }
        phi_nod[i] =  node.phi;
        }

	tiny::transposed_mult<double, Facette::N, Facette::NPI> (phi_nod, fac.a, phi);
		
	tiny::mult<double, 3, Facette::N, Facette::NPI> (u_nod, fac.a, u);
	
    for (int npi=0; npi<Facette::NPI; npi++){
        double un = u[0][npi]*nx + u[1][npi]*ny + u[2][npi]*nz;
        q[npi] = Ms * un;
        }
	
	/*-------------------------------------------------------*/	    
	double dens[Facette::NPI];
	double Eelem;
		
	for (int npi=0; npi<Facette::NPI; npi++) {
		// cosinus directeurs
		double al0=uk00*u[0][npi] + uk01*u[1][npi] + uk02*u[2][npi];
		dens[npi] = -K  * sq(al0);      // uniaxe
		}
	Eelem=tiny::sp<double, Facette::NPI> (dens, fac.weight);
	E[1]+=Eelem;

	for (int npi=0; npi<Facette::NPI; npi++) {
		dens[npi] = 0.5*mu0*q[npi]*phi[npi];
		}
	Eelem=tiny::sp<double, Facette::NPI> (dens, fac.weight);
	E[2]+=Eelem;
    }
	
for (int e=0; e<4; e++){
	fem.E[e] = E[e];
	fem.Etot+= E[e];
    }

fem.evol = fem.Etot-fem.Etot0;
fem.phy = 0;
}
