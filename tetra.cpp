/**
  Elementary matrix Calculation for a tetrahedron element 
  Theta Scheme for ech
 */ 

#include "config.h" //pour macro if_verbose

#include "tetra.h"
#include "pt3D.h"
#include "tiny.h"



using namespace Tetra;

void Tet::calc_u(double u[DIM][NPI],std::vector <Node> const& myNode)
{
for(int k=0;k<NPI;k++)
	{ 
	u[0][k]= myNode[ind[0]].u0.x()*a[0][k] + myNode[ind[1]].u0.x()*a[1][k] + myNode[ind[2]].u0.x()*a[2][k] + myNode[ind[3]].u0.x()*a[3][k]; 
	u[1][k]= myNode[ind[0]].u0.y()*a[0][k] + myNode[ind[1]].u0.y()*a[1][k] + myNode[ind[2]].u0.y()*a[2][k] + myNode[ind[3]].u0.y()*a[3][k];
	u[2][k]= myNode[ind[0]].u0.z()*a[0][k] + myNode[ind[1]].u0.z()*a[1][k] + myNode[ind[2]].u0.z()*a[2][k] + myNode[ind[3]].u0.z()*a[3][k];
	}
}

void Tet::calc_dudx(double dudx[DIM][NPI],std::vector <Node> const& myNode)
{
for(int k=0;k<NPI;k++)
	{ 
	dudx[0][k]= myNode[ind[0]].u0.x()*dadx[0][k] + myNode[ind[1]].u0.x()*dadx[1][k] + myNode[ind[2]].u0.x()*dadx[2][k] + myNode[ind[3]].u0.x()*dadx[3][k]; 
	dudx[1][k]= myNode[ind[0]].u0.y()*dadx[0][k] + myNode[ind[1]].u0.y()*dadx[1][k] + myNode[ind[2]].u0.y()*dadx[2][k] + myNode[ind[3]].u0.y()*dadx[3][k];
	dudx[2][k]= myNode[ind[0]].u0.z()*dadx[0][k] + myNode[ind[1]].u0.z()*dadx[1][k] + myNode[ind[2]].u0.z()*dadx[2][k] + myNode[ind[3]].u0.z()*dadx[3][k];
	}
}

void Tet::calc_dudy(double dudy[DIM][NPI],std::vector <Node> const& myNode)
{
for(int k=0;k<NPI;k++)
	{ 
	dudy[0][k]= myNode[ind[0]].u0.x()*dady[0][k] + myNode[ind[1]].u0.x()*dady[1][k] + myNode[ind[2]].u0.x()*dady[2][k] + myNode[ind[3]].u0.x()*dady[3][k]; 
	dudy[1][k]= myNode[ind[0]].u0.y()*dady[0][k] + myNode[ind[1]].u0.y()*dady[1][k] + myNode[ind[2]].u0.y()*dady[2][k] + myNode[ind[3]].u0.y()*dady[3][k];
	dudy[2][k]= myNode[ind[0]].u0.z()*dady[0][k] + myNode[ind[1]].u0.z()*dady[1][k] + myNode[ind[2]].u0.z()*dady[2][k] + myNode[ind[3]].u0.z()*dady[3][k];
	}
}

void Tet::calc_dudz(double dudz[DIM][NPI],std::vector <Node> const& myNode)
{
for(int k=0;k<NPI;k++)
	{ 
	dudz[0][k]= myNode[ind[0]].u0.x()*dadz[0][k] + myNode[ind[1]].u0.x()*dadz[1][k] + myNode[ind[2]].u0.x()*dadz[2][k] + myNode[ind[3]].u0.x()*dadz[3][k]; 
	dudz[1][k]= myNode[ind[0]].u0.y()*dadz[0][k] + myNode[ind[1]].u0.y()*dadz[1][k] + myNode[ind[2]].u0.y()*dadz[2][k] + myNode[ind[3]].u0.y()*dadz[3][k];
	dudz[2][k]= myNode[ind[0]].u0.z()*dadz[0][k] + myNode[ind[1]].u0.z()*dadz[1][k] + myNode[ind[2]].u0.z()*dadz[2][k] + myNode[ind[3]].u0.z()*dadz[3][k];
	}
}
//////////////////////////////////////

void Tet::calc_v(double v[DIM][NPI],std::vector <Node> const& myNode)
{
for(int k=0;k<NPI;k++)
	{ 
	v[0][k]= myNode[ind[0]].v0.x()*a[0][k] + myNode[ind[1]].v0.x()*a[1][k] + myNode[ind[2]].v0.x()*a[2][k] + myNode[ind[3]].v0.x()*a[3][k]; 
	v[1][k]= myNode[ind[0]].v0.y()*a[0][k] + myNode[ind[1]].v0.y()*a[1][k] + myNode[ind[2]].v0.y()*a[2][k] + myNode[ind[3]].v0.y()*a[3][k];
	v[2][k]= myNode[ind[0]].v0.z()*a[0][k] + myNode[ind[1]].v0.z()*a[1][k] + myNode[ind[2]].v0.z()*a[2][k] + myNode[ind[3]].v0.z()*a[3][k];
	}
}

void Tet::calc_dvdx(double dvdx[DIM][NPI],std::vector <Node> const& myNode)
{
for(int k=0;k<NPI;k++)
	{ 
	dvdx[0][k]= myNode[ind[0]].v0.x()*dadx[0][k] + myNode[ind[1]].v0.x()*dadx[1][k] + myNode[ind[2]].v0.x()*dadx[2][k] + myNode[ind[3]].v0.x()*dadx[3][k]; 
	dvdx[1][k]= myNode[ind[0]].v0.y()*dadx[0][k] + myNode[ind[1]].v0.y()*dadx[1][k] + myNode[ind[2]].v0.y()*dadx[2][k] + myNode[ind[3]].v0.y()*dadx[3][k];
	dvdx[2][k]= myNode[ind[0]].v0.z()*dadx[0][k] + myNode[ind[1]].v0.z()*dadx[1][k] + myNode[ind[2]].v0.z()*dadx[2][k] + myNode[ind[3]].v0.z()*dadx[3][k];
	}
}

void Tet::calc_dvdy(double dvdy[DIM][NPI],std::vector <Node> const& myNode)
{
for(int k=0;k<NPI;k++)
	{ 
	dvdy[0][k]= myNode[ind[0]].v0.x()*dady[0][k] + myNode[ind[1]].v0.x()*dady[1][k] + myNode[ind[2]].v0.x()*dady[2][k] + myNode[ind[3]].v0.x()*dady[3][k]; 
	dvdy[1][k]= myNode[ind[0]].v0.y()*dady[0][k] + myNode[ind[1]].v0.y()*dady[1][k] + myNode[ind[2]].v0.y()*dady[2][k] + myNode[ind[3]].v0.y()*dady[3][k];
	dvdy[2][k]= myNode[ind[0]].v0.z()*dady[0][k] + myNode[ind[1]].v0.z()*dady[1][k] + myNode[ind[2]].v0.z()*dady[2][k] + myNode[ind[3]].v0.z()*dady[3][k];
	}
}

void Tet::calc_dvdz(double dvdz[DIM][NPI],std::vector <Node> const& myNode)
{
for(int k=0;k<NPI;k++)
	{ 
	dvdz[0][k]= myNode[ind[0]].v0.x()*dadz[0][k] + myNode[ind[1]].v0.x()*dadz[1][k] + myNode[ind[2]].v0.x()*dadz[2][k] + myNode[ind[3]].v0.x()*dadz[3][k]; 
	dvdz[1][k]= myNode[ind[0]].v0.y()*dadz[0][k] + myNode[ind[1]].v0.y()*dadz[1][k] + myNode[ind[2]].v0.y()*dadz[2][k] + myNode[ind[3]].v0.y()*dadz[3][k];
	dvdz[2][k]= myNode[ind[0]].v0.z()*dadz[0][k] + myNode[ind[1]].v0.z()*dadz[1][k] + myNode[ind[2]].v0.z()*dadz[2][k] + myNode[ind[3]].v0.z()*dadz[3][k];
	}
}

void Tet::integrales(std::vector<Tetra::prm> const& params,std::vector <Node> const& myNode,double Hext[DIM],double Vz,double theta,double dt,gmm::dense_matrix <double> &AE, std::vector <double> &BE)
{
/*
std::map < std::pair<std::string,int>,double> &param = mySets.param;

double s = param[ std::make_pair("theta",-1) ];   	//cout << ", theta= " << s;
double alpha = param[ std::make_pair("alpha",reg) ];  	//cout << ", alpha= " << alpha;
double A = param[ std::make_pair("Ae",reg) ];    	//cout << ", Ae= " << A ;
double J = param[ std::make_pair("Js",reg) ]+EPSILON;//cout <<", Js=" << J;

double K = param[std::make_pair("Ka",reg)];   	//cout << ", Ka=" << K;
double K3 = param[std::make_pair("Ka3",reg)];   	//cout << ", Ka3=" << K3;

double uk00 = param[std::make_pair("a1",reg)];  	//cout << ", a1=" << k0;
double uk01 = param[std::make_pair("a2",reg)];  	//cout << ", a2=" << k1;
double uk02 = param[std::make_pair("a3",reg)];  	//cout << ", a3=" << k2;
double uk10 = param[std::make_pair("b1",reg)];  	//cout << ", a1=" << k0;
double uk11 = param[std::make_pair("b2",reg)];  	//cout << ", a2=" << k1;
double uk12 = param[std::make_pair("b3",reg)];  	//cout << ", a3=" << k2;
double uk20 = param[std::make_pair("c1",reg)];  	//cout << ", a1=" << k0;
double uk21 = param[std::make_pair("c2",reg)];  	//cout << ", a2=" << k1;
double uk22 = param[std::make_pair("c3",reg)];  	//cout << ", a3=" << k2;

double Uz   = param[std::make_pair("UzDW",reg)];//deplacement de paroi par courant polarise en spin
double beta = param[std::make_pair("betaDW",reg)];

*/

double alpha = params[idxPrm].alpha;
double A = params[idxPrm].A;
double J = params[idxPrm].J;

double K = params[idxPrm].K;
double K3 = params[idxPrm].K3;

double uk00 = params[idxPrm].uk[0][0];
double uk01 = params[idxPrm].uk[0][1];
double uk02 = params[idxPrm].uk[0][2];
double uk10 = params[idxPrm].uk[1][0];
double uk11 = params[idxPrm].uk[1][1];
double uk12 = params[idxPrm].uk[1][2];
double uk20 = params[idxPrm].uk[2][0];
double uk21 = params[idxPrm].uk[2][1];
double uk22 = params[idxPrm].uk[2][2];

double Uz = params[idxPrm].Uz;
double beta = params[idxPrm].beta;

/* ces constantes permettent de factoriser beaucoup d'expressions  */
double Abis = 2.0*A/J;
double Kbis = 2.0*K/J;
double K3bis = 2.0*K3/J;
double s_dt = theta*dt;//theta du theta schema, defini dans config.h

/*-------------------- INTERPOLATION --------------------*/
double u_nod[3][N]; 
double u[3][NPI];
double dudx[3][NPI], dudy[3][NPI], dudz[3][NPI];
double negphi0_nod[N], Hdx[NPI], Hdy[NPI], Hdz[NPI];

double v_nod[3][N];
double v[3][NPI];
double dvdx[3][NPI], dvdy[3][NPI], dvdz[3][NPI];
double negphiv0_nod[N], Hvx[NPI], Hvy[NPI], Hvz[NPI];

for (int i=0; i<N; i++){
    Node const& node = myNode[ ind[i] ];
    u_nod[Pt::IDX_X][i]  = node.u0.x(); u_nod[Pt::IDX_Y][i] = node.u0.y(); u_nod[Pt::IDX_Z][i]  = node.u0.z();
    v_nod[Pt::IDX_X][i]  = node.v0.x(); v_nod[Pt::IDX_Y][i] = node.v0.y(); v_nod[Pt::IDX_Z][i]  = node.v0.z();			
//for (int d=0; d<3; d++) { u_nod[d][i]  = node.u0(d); v_nod[d][i]  = node.v0(d); }
    negphi0_nod[i]  = -node.phi0;
    negphiv0_nod[i] = -node.phiv0;
    }

//calc_u(u,myNode);
tiny::mult<double, 3, N, NPI> (u_nod, a, u);
//calc_dudx(dudx,myNode);
tiny::mult<double, 3, N, NPI> (u_nod, dadx, dudx);
//calc_dudy(dudy,myNode);
tiny::mult<double, 3, N, NPI> (u_nod, dady, dudy);
//calc_dudz(dudz,myNode);
tiny::mult<double, 3, N, NPI> (u_nod, dadz, dudz);

//calc_v(v,myNode);
tiny::mult<double, 3, N, NPI> (v_nod, a, v);
//calc_dudx(dvdx,myNode);
tiny::mult<double, 3, N, NPI> (v_nod, dadx, dvdx);
//calc_dudx(dvdy,myNode);
tiny::mult<double, 3, N, NPI> (v_nod, dady, dvdy);
//calc_dudx(dvdz,myNode);
tiny::mult<double, 3, N, NPI> (v_nod, dadz, dvdz);

tiny::transposed_mult<double, N, NPI> (negphi0_nod, dadx, Hdx);
tiny::transposed_mult<double, N, NPI> (negphi0_nod, dady, Hdy);
tiny::transposed_mult<double, N, NPI> (negphi0_nod, dadz, Hdz);

tiny::transposed_mult<double, N, NPI> (negphiv0_nod, dadx, Hvx);
tiny::transposed_mult<double, N, NPI> (negphiv0_nod, dady, Hvy);
tiny::transposed_mult<double, N, NPI> (negphiv0_nod, dadz, Hvz);

/* on se place dans le referentiel mobile */
//double Vz=fem.DW_vz;

/*-------------------------------------------------------*/
for (int npi=0; npi<NPI; npi++){
    double ai, ai_w, dai_dx, dai_dy, dai_dz, daj_dx, daj_dy, daj_dz;
    double Dai_Daj, Dai_Du0, Dai_Du1, Dai_Du2;
    
	double w = weight[npi];
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

    triple Ht; // derivee de Hr : y a t'il une erreur ? on dirait que ce devrait etre Kbis*uk{0|1|2}_v et pas Kbis*ok0_v
    Ht[0]= Hvx[npi] + (Kbis* uk0_v - K3bis* uk0_v*(1-3*uk0_u*uk0_u) )*uk00;   
    Ht[1]= Hvy[npi] + (Kbis* uk0_v - K3bis* uk1_v*(1-3*uk1_u*uk1_u) )*uk01;   
    Ht[2]= Hvz[npi] + (Kbis* uk0_v - K3bis* uk2_v*(1-3*uk2_u*uk2_u) )*uk02;  

    for (int i=0; i<N; i++){
        ai = a[i][npi];
        dai_dx= dadx[i][npi];  dai_dy= dady[i][npi];  dai_dz= dadz[i][npi];
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
        BE[    i]+= Ht[0] *ai_w*s_dt; // ordre 2 en temps
        BE[  N+i]+= Ht[1] *ai_w*s_dt;
        BE[2*N+i]+= Ht[2] *ai_w*s_dt;

/* changement de referentiel */
        BE[i]    += +Vz*(u[1][npi]*dvdz[2][npi]-u[2][npi]*dvdz[1][npi]+v[1][npi]*dudz[2][npi]-v[2][npi]*dudz[1][npi]+alpha*dvdz[0][npi]) *ai_w*s_dt;
        BE[N+i]  += +Vz*(u[2][npi]*dvdz[0][npi]-u[0][npi]*dvdz[2][npi]+v[2][npi]*dudz[0][npi]-v[0][npi]*dudz[2][npi]+alpha*dvdz[1][npi]) *ai_w*s_dt;
        BE[2*N+i]+= +Vz*(u[0][npi]*dvdz[1][npi]-u[1][npi]*dvdz[0][npi]+v[0][npi]*dudz[1][npi]-v[1][npi]*dudz[0][npi]+alpha*dvdz[2][npi]) *ai_w*s_dt;

/* second membre pour les termes de courant polarise en spin pour une paroi  pour ordre 2 en temps*/
	BE[i]    += -Uz*(u[1][npi]*dvdz[2][npi]-u[2][npi]*dvdz[1][npi]+v[1][npi]*dudz[2][npi]-v[2][npi]*dudz[1][npi]+beta*dvdz[0][npi]) *ai_w*s_dt;
	BE[N+i]  += -Uz*(u[2][npi]*dvdz[0][npi]-u[0][npi]*dvdz[2][npi]+v[2][npi]*dudz[0][npi]-v[0][npi]*dudz[2][npi]+beta*dvdz[1][npi]) *ai_w*s_dt;
	BE[2*N+i]+= -Uz*(u[0][npi]*dvdz[1][npi]-u[1][npi]*dvdz[0][npi]+v[0][npi]*dudz[1][npi]-v[1][npi]*dudz[0][npi]+beta*dvdz[2][npi]) *ai_w*s_dt;
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
            daj_dx= dadx[j][npi];  daj_dy= dady[j][npi];  daj_dz= dadz[j][npi];
            Dai_Daj = dai_dx*daj_dx + dai_dy*daj_dy + dai_dz*daj_dz;

            AE(i,j)        +=  s_dt*(1.+R)* Abis* Dai_Daj *w;
            AE(N+i,N+j)    +=  s_dt*(1.+R)* Abis* Dai_Daj *w;
            AE(2*N+i,2*N+j)+=  s_dt*(1.+R)* Abis* Dai_Daj *w;
	    }
	}
    }
}

void Tet::calc_vol(std::vector<Node> const& myNode)
{
int i0,i1,i2,i3;
i0=ind[0];   i1=ind[1];   i2=ind[2];   i3=ind[3];
   
Pt::pt3D p0 = myNode[i0].p;
Pt::pt3D p1 = myNode[i1].p;
Pt::pt3D p2 = myNode[i2].p;
Pt::pt3D p3 = myNode[i3].p;
Pt::pt3D vec = (p1-p0)*(p2-p0);

double i_vol  = 1./6.* pScal(vec,p3-p0);
   if (i_vol<0.) {
      ind[3]=i2; ind[2]=i3;
      i_vol=-i_vol;
      IF_VERBOSE() std::cout << "ill-oriented tetrahedron, now corrected!"<< std::endl;
      }
vol = i_vol;
}



void Tetra::init_dadu(gmm::dense_matrix <double> &X)
	{
	X(0,0)= -1.0;   X(0,1)= -1.0;   X(0,2)= -1.0;
        X(1,0)= +1.0;   X(1,1)=  0.0;   X(1,2)=  0.0;
        X(2,0)=  0.0;   X(2,1)= +1.0;   X(2,2)=  0.0;
        X(3,0)=  0.0;   X(3,1)=  0.0;   X(3,2)= +1.0;
	}




