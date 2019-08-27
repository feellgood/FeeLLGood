/**
  Elementary matrix Calculation for a tetrahedron element 
 */ 

#include "config.h" //pour macro if_verbose

#include "tetra.h"
#include "pt3D.h"
#include "tiny.h"

#include "matBlocDiag.h"

using namespace Tetra;

void Tet::init(double epsilon)
{
double J[Pt::DIM][Pt::DIM];
double detJ = Jacobian(J);
double da[N][Pt::DIM];
    
if (fabs(detJ) < epsilon){
        #ifdef LIBRARY
            ostringstream what;
            what << "Singular jacobian in tetrahedron ";
            throw runtime_error(what.str());
        #else
            std::cerr << "jacobienne singuliere ds le tetraedre " << std::endl;
			infos();
            SYSTEM_ERROR;
        #endif
            }
Pt::inverse(J,detJ);
tiny::mult<double, N, Pt::DIM, Pt::DIM> (Tetra::dadu, J, da);
    
for (int j=0; j<NPI; j++)
    {
    for (int i=0; i<N; i++)
        { dadx[i][j]=da[i][0]; dady[i][j]=da[i][1]; dadz[i][j]=da[i][2]; }
    weight[j]    = detJ * Tetra::pds[j];
    }    
}

void Tet::integrales(std::vector<Tetra::prm> const& params,double Hext[DIM],double Vz,
                     double theta,double dt,double tau_r,gmm::dense_matrix <double> &AE, std::vector <double> &BE) const
{
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
double s_dt = theta*dt;//theta du theta schema

/*-------------------- INTERPOLATION --------------------*/
double u_nod[3][N]; 
double u[3][NPI], dudx[3][NPI], dudy[3][NPI], dudz[3][NPI];//should be Pt::pt3D X[NPI]
double negphi0_nod[N], Hdx[NPI], Hdy[NPI], Hdz[NPI];

double v_nod[3][N];
double v[3][NPI];
double dvdx[3][NPI], dvdy[3][NPI], dvdz[3][NPI];
double negphiv0_nod[N], Hvx[NPI], Hvy[NPI], Hvz[NPI];

for (int i=0; i<N; i++){
    Nodes::Node const& node = (*refNode)[ ind[i] ];
    u_nod[Pt::IDX_X][i]  = node.u0.x(); u_nod[Pt::IDX_Y][i] = node.u0.y(); u_nod[Pt::IDX_Z][i]  = node.u0.z();
    v_nod[Pt::IDX_X][i]  = node.v0.x(); v_nod[Pt::IDX_Y][i] = node.v0.y(); v_nod[Pt::IDX_Z][i]  = node.v0.z();			
    negphi0_nod[i]  = -node.phi0;
    negphiv0_nod[i] = -node.phiv0;
    }

    
tiny::mult<double, 3, N, NPI> (u_nod, a, u);// a is const matrix from namespace tetra 

tiny::mult<double, 3, N, NPI> (u_nod, dadx, dudx);
tiny::mult<double, 3, N, NPI> (u_nod, dady, dudy);
tiny::mult<double, 3, N, NPI> (u_nod, dadz, dudz);

tiny::mult<double, 3, N, NPI> (v_nod, a, v);

tiny::mult<double, 3, N, NPI> (v_nod, dadx, dvdx);
tiny::mult<double, 3, N, NPI> (v_nod, dady, dvdy);
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
    double ai, ai_w, dai_dx, dai_dy, dai_dz;
    double Dai_Daj, Dai_Du0, Dai_Du1, Dai_Du2, contrib;
    
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
    double uk0_uuu = (1-uk0_u*uk0_u)*uk0_u;
    double uk1_uuu = (1-uk1_u*uk1_u)*uk1_u;
    double uk2_uuu = (1-uk2_u*uk2_u)*uk2_u;
    
    //double uHa3u = -K3bis*(uk0_u*(1-uk0_u*uk0_u)*uk0_u + uk1_u*(1-uk1_u*uk1_u)*uk1_u + uk2_u*(1-uk2_u*uk2_u)*uk2_u);
    
    double uHa3u = -K3bis*(uk0_u*uk0_uuu + uk1_u*uk1_uuu + uk2_u*uk2_uuu);
    
    double uHeff = -Abis*Du2 +uHext +uHdu +uHau +uHa3u;

    double alfa=alpha; // seulement pour l'ordre 1 en temps
    double R=0.;

//second order corrections
double r = 0.1;	     			
double M = 2.*alpha*r/dt;  			
R = dt/tau_r*abs(log(dt/tau_r));    	

    if (uHeff>0.){ 
       if (uHeff>M) alfa=alpha+dt/2.*M;
       else alfa=alpha+dt/2.*uHeff;
       }
    else{
       if (uHeff<-M) alfa=alpha/(1.+dt/(2.*alpha)*M);
       else alfa=alpha/(1.-dt/(2.*alpha)*uHeff);
       }
// end second order corrections

    for (int i=0; i<N; i++){
        ai = a[i][npi];
        dai_dx= dadx[i][npi];  dai_dy= dady[i][npi];  dai_dz= dadz[i][npi];
        Dai_Du0 = dai_dx * dudx[0][npi] + dai_dy * dudy[0][npi] + dai_dz * dudz[0][npi];
        Dai_Du1 = dai_dx * dudx[1][npi] + dai_dy * dudy[1][npi] + dai_dz * dudz[1][npi];
        Dai_Du2 = dai_dx * dudx[2][npi] + dai_dy * dudy[2][npi] + dai_dz * dudz[2][npi];

/*        BE[i]    += (-Abis* Dai_Du0 + ( Kbis* uk0_u*uk00 -K3bis*( uk0_u*(1-uk0_u*uk0_u)*uk00 + uk1_u*(1-uk1_u*uk1_u)*uk10 + uk2_u*(1-uk2_u*uk2_u)*uk20 ) + Hdx[npi] + Hext[0] )*ai) *w;
        
        BE[N+i]  += (-Abis* Dai_Du1 + ( Kbis* uk0_u*uk01 -K3bis*( uk0_u*(1-uk0_u*uk0_u)*uk01 + uk1_u*(1-uk1_u*uk1_u)*uk11 + uk2_u*(1-uk2_u*uk2_u)*uk21 ) + Hdy[npi] + Hext[1] )*ai) *w;
        
        BE[2*N+i]+= (-Abis* Dai_Du2 + ( Kbis* uk0_u*uk02 -K3bis*( uk0_u*(1-uk0_u*uk0_u)*uk02 + uk1_u*(1-uk1_u*uk1_u)*uk12 + uk2_u*(1-uk2_u*uk2_u)*uk22 ) + Hdz[npi] + Hext[2] )*ai ) *w;
*/
        BE[i]    += (-Abis* Dai_Du0 + ( Kbis* uk0_u*uk00 -K3bis*( uk0_uuu*uk00 + uk1_uuu*uk10 + uk2_uuu*uk20 ) + Hdx[npi] + Hext[0] )*ai) *w;
        
        BE[N+i]  += (-Abis* Dai_Du1 + ( Kbis* uk0_u*uk01 -K3bis*( uk0_uuu*uk01 + uk1_uuu*uk11 + uk2_uuu*uk21 ) + Hdy[npi] + Hext[1] )*ai) *w;
        
        BE[2*N+i]+= (-Abis* Dai_Du2 + ( Kbis* uk0_u*uk02 -K3bis*( uk0_uuu*uk02 + uk1_uuu*uk12 + uk2_uuu*uk22 ) + Hdz[npi] + Hext[2] )*ai) *w;
        
	ai_w = ai*w;
/* changement de referentiel */
        BE[i]    += +Vz*(u[1][npi]*dudz[2][npi]-u[2][npi]*dudz[1][npi]+alpha*dudz[0][npi]) *ai_w;
        BE[N+i]  += +Vz*(u[2][npi]*dudz[0][npi]-u[0][npi]*dudz[2][npi]+alpha*dudz[1][npi]) *ai_w;
        BE[2*N+i]+= +Vz*(u[0][npi]*dudz[1][npi]-u[1][npi]*dudz[0][npi]+alpha*dudz[2][npi]) *ai_w;

/* second membre pour les termes de courant polarise en spin pour une paroi */
	BE[i]    += -Uz*(u[1][npi]*dudz[2][npi]-u[2][npi]*dudz[1][npi]+beta*dudz[0][npi]) *ai_w;
	BE[N+i]  += -Uz*(u[2][npi]*dudz[0][npi]-u[0][npi]*dudz[2][npi]+beta*dudz[1][npi]) *ai_w;
	BE[2*N+i]+= -Uz*(u[0][npi]*dudz[1][npi]-u[1][npi]*dudz[0][npi]+beta*dudz[2][npi]) *ai_w;

//second order corrections
    triple Ht; //derivee de Hr : y a t'il une erreur ? on dirait que ce devrait etre Kbis*uk{0|1|2}_v et pas Kbis*ok0_v
    Ht[0]= Hvx[npi] + (Kbis* uk0_v - K3bis* uk0_v*(1-3*uk0_u*uk0_u) )*uk00;   
    Ht[1]= Hvy[npi] + (Kbis* uk0_v - K3bis* uk1_v*(1-3*uk1_u*uk1_u) )*uk01;   
    Ht[2]= Hvz[npi] + (Kbis* uk0_v - K3bis* uk2_v*(1-3*uk2_u*uk2_u) )*uk02; 
        
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
// end second order corrections

        
        AE(    i,    i)+=  alfa* ai_w;  //lumping
        AE(  N+i,  N+i)+=  alfa* ai_w;
        AE(2*N+i,2*N+i)+=  alfa* ai_w;

        AE(0*N+i,2*N+i)+= +u_nod[1][i]* ai_w; //lumping
        AE(0*N+i,1*N+i)+= -u_nod[2][i]* ai_w;
        AE(1*N+i,0*N+i)+= +u_nod[2][i]* ai_w;
        AE(1*N+i,2*N+i)+= -u_nod[0][i]* ai_w;
        AE(2*N+i,1*N+i)+= +u_nod[0][i]* ai_w;
        AE(2*N+i,0*N+i)+= -u_nod[1][i]* ai_w;

        for (int j=0; j<N; j++)
            {
            Dai_Daj = dai_dx*dadx[j][npi] + dai_dy*dady[j][npi] + dai_dz*dadz[j][npi];
            contrib = s_dt*(1.+R)* Abis* Dai_Daj *w;
            AE(i,j)        +=  contrib;
            AE(N+i,N+j)    +=  contrib;
            AE(2*N+i,2*N+j)+=  contrib;
            }
        
        
        }
    }
}

void Tet::energy(Tetra::prm const& param,double E[5],const double Hext[DIM],double uz_drift) const
{
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
   
   /*-------------------- INTERPOLATION --------------------*/
double u_nod[3][N], u[3][NPI];
double dudx[3][NPI], dudy[3][NPI], dudz[3][NPI];
double q[NPI],  phi[NPI];
double phi_nod[N], negphi_nod[N], Hdx[NPI], Hdy[NPI], Hdz[NPI];

for (int i=0; i<N; i++)
    {
    int i_= ind[i];
    Nodes::Node &n = (*refNode)[i_];
    u_nod[Pt::IDX_X][i] = n.u.x(); u_nod[Pt::IDX_Y][i] = n.u.y(); u_nod[Pt::IDX_Z][i] = n.u.z();
    phi_nod[i] =  n.phi;
    negphi_nod[i] = -n.phi;
    }

tiny::transposed_mult<double, N, NPI> (phi_nod, a, phi);
tiny::mult<double, 3, N, NPI> (u_nod, a, u);
tiny::mult<double, 3, N, NPI> (u_nod, dadx, dudx);
tiny::mult<double, 3, N, NPI> (u_nod, dady, dudy);
tiny::mult<double, 3, N, NPI> (u_nod, dadz, dudz);

for (int npi=0; npi<NPI; npi++){
        //double div_u = dudx[0][npi] + dudy[1][npi] + dudz[2][npi];
    q[npi] = -Ms*(dudx[0][npi] + dudy[1][npi] + dudz[2][npi]);
    }

tiny::transposed_mult<double, N, NPI> (negphi_nod, dadx, Hdx);
tiny::transposed_mult<double, N, NPI> (negphi_nod, dady, Hdy);
tiny::transposed_mult<double, N, NPI> (negphi_nod, dadz, Hdz);

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
    tiny::mult<double, 5, NPI> (dens, weight, Eelem);
    E[0] += Eelem[0];
    E[1] += Eelem[1];
    E[2] += Eelem[2];
    E[3] += Eelem[3];
    E[4] += Eelem[4];
}

void Tet::projection(//mtl::dense2D <double> &P,
           gmm::dense_matrix <double> const& A,  std::vector <double> const& B,
           gmm::dense_matrix <double> &Ap, std::vector <double> &Bp) const
{
thread_local gmm::dense_matrix <double> P(2*N,3*N);
thread_local gmm::dense_matrix <double> PA(2*N,3*N);
//mtl::mat::set_to_zero(P);
//matBlocDiag::matBloc myP;

//myP.D[0][0] = matBlocDiag::BlocElem(1,2,3,4);

for (int i=0; i<N; i++){
    Nodes::Node const& n = (*refNode)[ind[i]];
    P(i,i)  = n.ep.x();  P(i,N+i)  = n.ep.y();  P(i,2*N+i)  = n.ep.z();
    P(N+i,i)= n.eq.x();  P(N+i,N+i)= n.eq.y();  P(N+i,2*N+i)= n.eq.z();
    }

    gmm::mult(P,A,PA);
//Ap = (P*A)*trans(P);
//Bp = P*B;
gmm::mult(PA, gmm::transposed(P), Ap);

gmm::mult(P,B,Bp);
    
}

void Tet::projection(gmm::dense_matrix <double> const& A,  std::vector <double> const& B)
{
thread_local gmm::dense_matrix <double> P(2*N,3*N);
thread_local gmm::dense_matrix <double> PA(2*N,3*N);

for (int i=0; i<N; i++){
    Nodes::Node const& n = (*refNode)[ind[i]];
    P(i,i)  = n.ep.x();  P(i,N+i)  = n.ep.y();  P(i,2*N+i)  = n.ep.z();
    P(N+i,i)= n.eq.x();  P(N+i,N+i)= n.eq.y();  P(N+i,2*N+i)= n.eq.z();
    }

    gmm::mult(P,A,PA);
//Ap = (P*A)*trans(P);
//Bp = P*B;
gmm::mult(PA, gmm::transposed(P), Kp);

gmm::mult(P,B,Lp);
}

void Tet::assemblage(gmm::dense_matrix <double> const& Ke, std::vector <double> const& Le,
                     write_matrix &K,write_vector &L) const
    {
    for (int i=0; i < N; i++)
        {
        int i_= ind[i];             
        
        for (int j=0; j < N; j++)
            {
            int j_= ind[j];
            K(NOD+i_,j_) += Ke(i,j);      K(NOD+i_, NOD+j_) += Ke(  i,N+j);
            K(    i_,j_) += Ke(N+i,j);    K(    i_, NOD+j_) += Ke(N+i,N+j);
            }
        L[NOD+i_] += Le[i];
        L[i_] += Le[N+i];
        
        }
    }

void Tet::assemblage(write_matrix &K,write_vector &L) const
{
for (int i=0; i < N; i++)
    {
    int i_= ind[i];             
        
    for (int j=0; j < N; j++)
        {
        int j_= ind[j];
        K(NOD+i_,j_) += Kp(i,j);      K(NOD+i_, NOD+j_) += Kp(  i,N+j);
        K(    i_,j_) += Kp(N+i,j);    K(    i_, NOD+j_) += Kp(N+i,N+j);
        }
    L[NOD+i_] += Lp[i];
    L[i_] += Lp[N+i];
    }    
}
    
void Tet::getNod(gmm::dense_matrix <double> &nod)
{
for (int i=0; i<N; i++)
    {
    int i_= ind[i];
    nod(0,i) = (*refNode)[i_].p.x();
    nod(1,i) = (*refNode)[i_].p.y();
    nod(2,i) = (*refNode)[i_].p.z();
    }
}

double Tet::Jacobian(double J[DIM][DIM])
{
Pt::pt3D p0 = (*refNode)[ ind[0] ].p;
Pt::pt3D p1 = (*refNode)[ ind[1] ].p;
Pt::pt3D p2 = (*refNode)[ ind[2] ].p;
Pt::pt3D p3 = (*refNode)[ ind[3] ].p;
J[0][0] = p1.x()-p0.x(); J[0][1] = p2.x()-p0.x(); J[0][2] = p3.x()-p0.x();   
J[1][0] = p1.y()-p0.y(); J[1][1] = p2.y()-p0.y(); J[1][2] = p3.y()-p0.y();
J[2][0] = p1.z()-p0.z(); J[2][1] = p2.z()-p0.z(); J[2][2] = p3.z()-p0.z();
    
return Pt::det(J);
}


void Tet::calc_vol(void)
{
int i0,i1,i2,i3;
i0=ind[0];   i1=ind[1];   i2=ind[2];   i3=ind[3];
   
Pt::pt3D p0 = (*refNode)[i0].p;
Pt::pt3D p1 = (*refNode)[i1].p;
Pt::pt3D p2 = (*refNode)[i2].p;
Pt::pt3D p3 = (*refNode)[i3].p;
Pt::pt3D vec = (p1-p0)*(p2-p0);

vol  = 1./6.* pScal(vec,p3-p0);
   if (vol<0.) {
      ind[3]=i2; ind[2]=i3;
      vol *= -1;
      if(VERBOSE) { std::cout << "ill-oriented tetrahedron, now corrected!"<< std::endl; }
      }
}
