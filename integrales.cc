#include "fem.h"
#include "tiny.h"

/**
  Calcul des matrices elementaires sur un element
  Schema en theta pr ech
 */

void integrales(Fem &fem,Settings &settings, Tet &tet, gmm::dense_matrix <double> &AE, vector <double> &BE)
{
const int N   = Tet::N;
const int NPI = Tet::NPI;
const int reg = tet.reg;

pair <string,int> p;
map <pair<string,int>,double> &param = settings.param;
triple &Hext=fem.Hext;

//cout << "Hext " << Hext[0] << "\t" << Hext[1] << "\t" << Hext[2] << endl;
//cout << "fem.Hext " << fem.Hext[0] << "\t" << fem.Hext[1] << "\t" << fem.Hext[2] << endl;

p=make_pair("theta",-1);    double s = param[p];   	//cout << ", theta= " << s;

p=make_pair("alpha",reg);   double alpha = param[p];  	//cout << ", alpha= " << alpha;
p=make_pair("Ae",reg);      double A = param[p];    	//cout << ", Ae= " << A ;
p=make_pair("Js",reg);      double J = param[p]+EPSILON;//cout <<", Js=" << J;


p=make_pair("Ka",reg);      double K = param[p];   	//cout << ", Ka=" << K;
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

/* deplacement de paroi par courant polarise en spin*/
p=make_pair("UzDW",reg);      double Uz   = param[p];
p=make_pair("betaDW",reg);    double beta = param[p];

/* ces constantes permettent de factoriser beaucoup d'expressions  */
double Abis = 2.0*A/J;
double Kbis = 2.0*K/J;
double K3bis = 2.0*K3/J;

double dt = settings.dt;

/*-------------------- INTERPOLATION --------------------*/
double u_nod[3][N], u[3][NPI];
double dudx[3][NPI], dudy[3][NPI], dudz[3][NPI];
double negphi0_nod[N], Hdx[NPI], Hdy[NPI], Hdz[NPI];

double v_nod[3][N], v[3][NPI];
double dvdx[3][NPI], dvdy[3][NPI], dvdz[3][NPI];
double negphiv0_nod[N], Hvx[NPI], Hvy[NPI], Hvz[NPI];

for (int i=0; i<N; i++){
    int i_= tet.ind[i];
    Node &node = fem.node[i_];
    for (int d=0; d<3; d++){
        u_nod[d][i]  = node.u0[d];
        v_nod[d][i]  = node.v0[d];
        }
    negphi0_nod[i]  = -node.phi0;
    negphiv0_nod[i] = -node.phiv0;
    }

tiny::mult<double, 3, N, NPI> (u_nod, tet.a, u);
tiny::mult<double, 3, N, NPI> (u_nod, tet.dadx, dudx);
tiny::mult<double, 3, N, NPI> (u_nod, tet.dady, dudy);
tiny::mult<double, 3, N, NPI> (u_nod, tet.dadz, dudz);

tiny::mult<double, 3, N, NPI> (v_nod, tet.a, v);
tiny::mult<double, 3, N, NPI> (v_nod, tet.dadx, dvdx);
tiny::mult<double, 3, N, NPI> (v_nod, tet.dady, dvdy);
tiny::mult<double, 3, N, NPI> (v_nod, tet.dadz, dvdz);

tiny::transposed_mult<double, N, NPI> (negphi0_nod, tet.dadx, Hdx);
tiny::transposed_mult<double, N, NPI> (negphi0_nod, tet.dady, Hdy);
tiny::transposed_mult<double, N, NPI> (negphi0_nod, tet.dadz, Hdz);

tiny::transposed_mult<double, N, NPI> (negphiv0_nod, tet.dadx, Hvx);
tiny::transposed_mult<double, N, NPI> (negphiv0_nod, tet.dady, Hvy);
tiny::transposed_mult<double, N, NPI> (negphiv0_nod, tet.dadz, Hvz);

/* on se place dans le referentiel mobile */
double Vz=fem.DW_vz;

/*-------------------------------------------------------*/
for (int npi=0; npi<NPI; npi++){
    double ai, ai_w, dai_dx, dai_dy, dai_dz, daj_dx, daj_dy, daj_dz;
    double Dai_Daj, Dai_Du0, Dai_Du1, Dai_Du2;
    
	double w = tet.weight[npi];
    double uk0_u = uk00*u[0][npi] + uk01*u[1][npi] + uk02*u[2][npi]; 
    double uk1_u = uk10*u[0][npi] + uk11*u[1][npi] + uk12*u[2][npi]; 
    double uk2_u = uk20*u[0][npi] + uk21*u[1][npi] + uk22*u[2][npi]; 

    double uk0_v = uk00*v[0][npi] + uk01*v[1][npi] + uk02*v[2][npi]; 
    double uk1_v = uk10*v[0][npi] + uk11*v[1][npi] + uk12*v[2][npi]; 
    double uk2_v = uk20*v[0][npi] + uk21*v[1][npi] + uk22*v[2][npi]; 

    double Du2 = sq(dudx[0][npi]) + sq(dudy[0][npi]) + sq(dudz[0][npi]) +
                     sq(dudx[1][npi]) + sq(dudy[1][npi]) + sq(dudz[1][npi]) +
                     sq(dudx[2][npi]) + sq(dudy[2][npi]) + sq(dudz[2][npi]) ;

    double uHext= u[0][npi]*Hext[0]  + u[1][npi]*Hext[1]  + u[2][npi]*Hext[2];
    double uHdu = u[0][npi]*Hdx[npi] + u[1][npi]*Hdy[npi] + u[2][npi]*Hdz[npi];

    double uHau = Kbis* uk0_u*uk0_u;
    double uHa3u = -K3bis*(uk0_u*(1-uk0_u*uk0_u)*uk0_u + uk1_u*(1-uk1_u*uk1_u)*uk1_u + uk2_u*(1-uk2_u*uk2_u)*uk2_u);

    double uHeff = -Abis*Du2 +uHext +uHdu +uHau +uHa3u;

    double alfa=alpha; // seulement pour l'ordre 1 en temps
    double R=0.;

#ifdef STAT
    gsl_histogram_increment (fem.stat.h, uHeff);
#endif

#ifdef ORD2
    double r = 0.1;	     			

    double M = 2.*alpha*r/dt;  			
R = dt/TAUR*abs(log(dt/TAUR));    	

#ifdef STAT
fem.stat.r = r;
fem.stat.M = M;
fem.stat.R = R;
#endif

    if (uHeff>0.){ 
       if (uHeff>M) alfa=alpha+dt/2.*M;
       else alfa=alpha+dt/2.*uHeff;
       }
    else{
       if (uHeff<-M) alfa=alpha/(1.+dt/(2.*alpha)*M);
       else alfa=alpha/(1.-dt/(2.*alpha)*uHeff);
       }
#endif

    triple Ht; // derivee de Hr
    Ht[0]= Hvx[npi] + (Kbis* uk0_v - K3bis* uk0_v*(1-3*uk0_u*uk0_u) )*uk00;   
    Ht[1]= Hvy[npi] + (Kbis* uk0_v - K3bis* uk1_v*(1-3*uk1_u*uk1_u) )*uk01;   
    Ht[2]= Hvz[npi] + (Kbis* uk0_v - K3bis* uk2_v*(1-3*uk2_u*uk2_u) )*uk02;  

    for (int i=0; i<N; i++){
        ai = tet.a[i][npi];
        dai_dx= tet.dadx[i][npi];  dai_dy= tet.dady[i][npi];  dai_dz= tet.dadz[i][npi];
        Dai_Du0 = dai_dx * dudx[0][npi] + dai_dy * dudy[0][npi] + dai_dz * dudz[0][npi];
        Dai_Du1 = dai_dx * dudx[1][npi] + dai_dy * dudy[1][npi] + dai_dz * dudz[1][npi];
        Dai_Du2 = dai_dx * dudx[2][npi] + dai_dy * dudy[2][npi] + dai_dz * dudz[2][npi];

        BE[i]    += (-Abis* Dai_Du0 + ( Kbis* uk0_u*uk00 -K3bis*( uk0_u*(1-uk0_u*uk0_u)*uk00 + uk1_u*(1-uk1_u*uk1_u)*uk10 + uk2_u*(1-uk2_u*uk2_u)*uk20 ) + Hdx[npi] + Hext[0] )*ai) *w;
        BE[N+i]  += (-Abis* Dai_Du1 + ( Kbis* uk0_u*uk01 -K3bis*( uk0_u*(1-uk0_u*uk0_u)*uk01 + uk1_u*(1-uk1_u*uk1_u)*uk11 + uk2_u*(1-uk2_u*uk2_u)*uk21 ) + Hdy[npi] + Hext[1] )*ai) *w;
        BE[2*N+i]+= (-Abis* Dai_Du2 + ( Kbis* uk0_u*uk02 -K3bis*( uk0_u*(1-uk0_u*uk0_u)*uk02 + uk1_u*(1-uk1_u*uk1_u)*uk12 + uk2_u*(1-uk2_u*uk2_u)*uk22 ) + Hdz[npi] + Hext[2] )*ai ) *w;

	ai_w = ai*w;
/* changement de referentiel */
        BE[i]    += +Vz*(u[1][npi]*dudz[2][npi]-u[2][npi]*dudz[1][npi]+alpha*dudz[0][npi]) *ai_w;
        BE[N+i]  += +Vz*(u[2][npi]*dudz[0][npi]-u[0][npi]*dudz[2][npi]+alpha*dudz[1][npi]) *ai_w;
        BE[2*N+i]+= +Vz*(u[0][npi]*dudz[1][npi]-u[1][npi]*dudz[0][npi]+alpha*dudz[2][npi]) *ai_w;

/* second membre pour les termes de courant polarise en spin pour une paroi */
	BE[i]    += -Uz*(u[1][npi]*dudz[2][npi]-u[2][npi]*dudz[1][npi]+beta*dudz[0][npi]) *ai_w;
	BE[N+i]  += -Uz*(u[2][npi]*dudz[0][npi]-u[0][npi]*dudz[2][npi]+beta*dudz[1][npi]) *ai_w;
	BE[2*N+i]+= -Uz*(u[0][npi]*dudz[1][npi]-u[1][npi]*dudz[0][npi]+beta*dudz[2][npi]) *ai_w;


#ifdef ORD2
        BE[    i]+= Ht[0] *ai_w *s*dt; // ordre 2 en temps
        BE[  N+i]+= Ht[1] *ai_w *s*dt;
        BE[2*N+i]+= Ht[2] *ai_w *s*dt;

/* changement de referentiel */
        BE[i]    += +Vz*(u[1][npi]*dvdz[2][npi]-u[2][npi]*dvdz[1][npi]+v[1][npi]*dudz[2][npi]-v[2][npi]*dudz[1][npi]+alpha*dvdz[0][npi]) *ai_w *s*dt;
        BE[N+i]  += +Vz*(u[2][npi]*dvdz[0][npi]-u[0][npi]*dvdz[2][npi]+v[2][npi]*dudz[0][npi]-v[0][npi]*dudz[2][npi]+alpha*dvdz[1][npi]) *ai_w *s*dt;
        BE[2*N+i]+= +Vz*(u[0][npi]*dvdz[1][npi]-u[1][npi]*dvdz[0][npi]+v[0][npi]*dudz[1][npi]-v[1][npi]*dudz[0][npi]+alpha*dvdz[2][npi]) *ai_w *s*dt;

/* second membre pour les termes de courant polarise en spin pour une paroi  pour ordre 2 en temps*/
	BE[i]    += -Uz*(u[1][npi]*dvdz[2][npi]-u[2][npi]*dvdz[1][npi]+v[1][npi]*dudz[2][npi]-v[2][npi]*dudz[1][npi]+beta*dvdz[0][npi]) *ai_w *s*dt;
	BE[N+i]  += -Uz*(u[2][npi]*dvdz[0][npi]-u[0][npi]*dvdz[2][npi]+v[2][npi]*dudz[0][npi]-v[0][npi]*dudz[2][npi]+beta*dvdz[1][npi]) *ai_w *s*dt;
	BE[2*N+i]+= -Uz*(u[0][npi]*dvdz[1][npi]-u[1][npi]*dvdz[0][npi]+v[0][npi]*dudz[1][npi]-v[1][npi]*dudz[0][npi]+beta*dvdz[2][npi]) *ai_w *s*dt;
#endif

        AE(    i,    i)+=  alfa* ai_w;  //lumping
        AE(  N+i,  N+i)+=  alfa* ai_w;
        AE(2*N+i,2*N+i)+=  alfa* ai_w;

        AE(0*N+i,2*N+i)+= +u_nod[1][i]* ai_w; //lumping
        AE(0*N+i,1*N+i)+= -u_nod[2][i]* ai_w;
        AE(1*N+i,0*N+i)+= +u_nod[2][i]* ai_w;
        AE(1*N+i,2*N+i)+= -u_nod[0][i]* ai_w;
        AE(2*N+i,1*N+i)+= +u_nod[0][i]* ai_w;
        AE(2*N+i,0*N+i)+= -u_nod[1][i]* ai_w;

        for (int j=0; j<N; j++){
            //aj = tet.a[j][npi];
            daj_dx= tet.dadx[j][npi];  daj_dy= tet.dady[j][npi];  daj_dz= tet.dadz[j][npi];
            Dai_Daj = dai_dx*daj_dx + dai_dy*daj_dy + dai_dz*daj_dz;

            AE(i,j)        +=  s*dt*(1.+R)* Abis* Dai_Daj *w;
            AE(N+i,N+j)    +=  s*dt*(1.+R)* Abis* Dai_Daj *w;
            AE(2*N+i,2*N+j)+=  s*dt*(1.+R)* Abis* Dai_Daj *w;
	    }
	}
    }

//cout << AE << endl;
//cout << "temps integrales: "<< diff_t << endl;
}

void integrales(Fem &fem,Settings &settings, Fac &fac, gmm::dense_matrix <double> &AE, vector <double> &BE)
{
const int N   = Fac::N;
const int NPI = Fac::NPI;
const int reg = fac.reg;

pair <string,int> p;
map <pair<string,int>,double> &param = settings.param;
//triple &Hext=fem.Hext;
//cout << "Hext " << Hext[0] << "\t" << Hext[1] << "\t" << Hext[2] << endl;
//cout << "fem.Hext " << fem.Hext[0] << "\t" << fem.Hext[1] << "\t" << fem.Hext[2] << endl;

p=make_pair("theta",-1);    //double s = param[p];   	cout << ", theta= " << s;
p=make_pair("alpha",reg);   //double alpha = param[p];  	cout << ", alpha= " << alpha;

p=make_pair("Js",reg);      double J = abs(param[p])+EPSILON;//cout << ", Js=" << J;
p=make_pair("Ks",reg);      double K = param[p];   	//cout << ", Ks=" << K;

p=make_pair("a1",reg);      double uk00 = param[p];  	//cout << ", a1=" << k0;
p=make_pair("a2",reg);      double uk01 = param[p];  	//cout << ", a2=" << k1;
p=make_pair("a3",reg);      double uk02 = param[p];  	//cout << ", a3=" << k2;

double Kbis = 2.0*K/J;
/*-------------------- INTERPOLATION --------------------*/
double u_nod[3][N], u[3][NPI];

for (int i=0; i<N; i++){
    int i_= fac.ind[i];
    Node &node = fem.node[i_];
    for (int d=0; d<3; d++)
        u_nod[d][i]   = node.u0[d];
    }

tiny::mult<double, 3, N, NPI> (u_nod, fac.a, u);

/*-------------------------------------------------------*/

for (int npi=0; npi<NPI; npi++){
    double Kbis_ai;
    double w_uk0_u = fac.weight[npi]*(uk00*u[0][npi] + uk01*u[1][npi] + uk02*u[2][npi]); 

    for (int i=0; i<N; i++){
        Kbis_ai = Kbis*fac.a[i][npi];

        BE[i]    += (Kbis_ai* w_uk0_u*uk00);
        BE[N+i]  += (Kbis_ai* w_uk0_u*uk01);
        BE[2*N+i]+= (Kbis_ai* w_uk0_u*uk02);
	}
    }

//cout << AE << endl;
//cout << "temps integrales: "<< diff_t << endl;
}
