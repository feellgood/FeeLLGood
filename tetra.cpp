/**
  Elementary matrix Calculation for a tetrahedron element 
 */ 

#include "config.h" //pour macro if_verbose

#include "tetra.h"
#include "pt3D.h"
#include "tiny.h"
#include "matBlocDiag.h"

using namespace Tetra;
using namespace Pt;

void Tet::integrales(std::vector<Tetra::prm> const& params,double dt,Pt::pt3D const& Hext,double tau_r,double Vz, double AE[3*N][3*N], double *BE) const
{
double alpha = params[idxPrm].alpha;
pt3D uk[DIM] = { params[idxPrm].uk[0], params[idxPrm].uk[1], params[idxPrm].uk[2]}; 
double Uz = params[idxPrm].Uz;
double beta = params[idxPrm].beta;

double Abis = 2.0*(params[idxPrm].A)/(params[idxPrm].J);
double Kbis = 2.0*(params[idxPrm].K)/(params[idxPrm].J);
double K3bis = 2.0*(params[idxPrm].K3)/(params[idxPrm].J);
double s_dt = THETA*dt;//theta from theta scheme in config.h.in

/*-------------------- INTERPOLATION --------------------*/
double u_nod[DIM][N]; 
pt3D Hd[NPI], dUdx[NPI], dUdy[NPI], dUdz[NPI];
pt3D Hv[NPI], dVdx[NPI], dVdy[NPI], dVdz[NPI];

// u_nod(u0) is required for lumping in AE matrix
getVecDataFromNode(Nodes::get_u0,u_nod);

pt3D U[NPI],V[NPI];
interpolation(Nodes::get_u0,U,dUdx,dUdy,dUdz);
interpolation(Nodes::get_v0,V,dVdx,dVdy,dVdz);
interpolation(Nodes::get_phi0,Hd);
interpolation(Nodes::get_phiv0,Hv);

for (int npi=0; npi<NPI; npi++){
    pt3D uk_u = pt3D(pScal(uk[0],U[npi]), pScal(uk[1],U[npi]), pScal(uk[2],U[npi]));
    pt3D uk_v = pt3D(pScal(uk[0],V[npi]), pScal(uk[1],V[npi]), pScal(uk[2],V[npi]));

Pt::pt3D uk_uuu = pDirect(pt3D(1,1,1) - pDirect(uk_u,uk_u), uk_u);
    
    double uHeff = -Abis*(norme2(dUdx[npi]) + norme2(dUdy[npi]) + norme2(dUdz[npi])); 
	uHeff +=  pScal(U[npi], Hext + Hd[npi]) + Kbis*sq(uk_u(0)) - K3bis*pScal(uk_u,uk_uuu);

    double alfa=calc_alpha_eff(alpha,dt,uHeff);
        
    //second order corrections,Ht = derivee de Hr : y a t'il une erreur ? on dirait que ce devrait etre Kbis*uk_v({0|1|2}) et pas Kbis*uk_v(0)
    pt3D truc = Kbis*uk_v(0)*pt3D(1,1,1) -K3bis*pDirect(uk_v , pt3D(1,1,1)-3*pDirect(uk_u,uk_u));
    pt3D Ht = s_dt*(Hv[npi] + pDirect(truc,uk[0]));
    
    pt3D Heff = Kbis*uk_u(0)*uk[0] - K3bis*( uk_uuu(0)*uk[0] + uk_uuu(1)*uk[1] + uk_uuu(2)*uk[2] ) + Hd[npi] + Hext;
    
    double w = weight[npi];
    double prefactor_contrib = w*s_dt*(1.+ dt/tau_r*abs(log(dt/tau_r)))* Abis;
    
    for (int i=0; i<N; i++){
        double ai_w = w*a[i][npi];
        pt3D Dai_Du = dadx[i][npi]*dUdx[npi] + dady[i][npi]*dUdy[npi] + dadz[i][npi]*dUdz[npi];
        
        pt3D X =  -w*Abis*Dai_Du + ai_w*(alpha*Vz - beta*Uz)*(dUdz[npi] + dVdz[npi]*s_dt) ;
        
        X += ai_w*(Heff + Ht + (Vz-Uz)*(U[npi]*(dUdz[npi]+dVdz[npi]*s_dt) +V[npi]*dUdz[npi]*s_dt) );
        
        BE[i]    += X(0);
        BE[N+i]  += X(1);
        BE[2*N+i]+= X(2);
        
        AE[    i][    i] +=  alfa* ai_w;  //lumping
        AE[  N+i][  N+i] +=  alfa* ai_w;
        AE[2*N+i][2*N+i] +=  alfa* ai_w;

        AE[0*N+i][2*N+i] += ai_w*u_nod[1][i]; //lumping
        AE[0*N+i][1*N+i] += -ai_w*u_nod[2][i];
        AE[1*N+i][0*N+i] += ai_w*u_nod[2][i];
        AE[1*N+i][2*N+i] += -ai_w*u_nod[0][i];
        AE[2*N+i][1*N+i] += ai_w*u_nod[0][i];
        AE[2*N+i][0*N+i] += -ai_w*u_nod[1][i];

        for (int j=0; j<N; j++)
            {
            double contrib = prefactor_contrib*(dadx[i][npi]*dadx[j][npi] + dady[i][npi]*dady[j][npi] + dadz[i][npi]*dadz[j][npi]);
            
            AE[i][j]        +=  contrib;
            AE[N+i][N+j]    +=  contrib;
            AE[2*N+i][2*N+j] +=  contrib;
            }
        }
    }
}

double Tet::exchangeEnergy(Tetra::prm const& param,const double dudx[DIM][NPI],const double dudy[DIM][NPI],const double dudz[DIM][NPI]) const
{
double dens[NPI];

for (int npi=0; npi<NPI; npi++)
    {
    dens[npi] = sq( dudx[0][npi] ) + sq( dudy[0][npi] ) + sq( dudz[0][npi] )+
                sq( dudx[1][npi] ) + sq( dudy[1][npi] ) + sq( dudz[1][npi] )+
                sq( dudx[2][npi] ) + sq( dudy[2][npi] ) + sq( dudz[2][npi] ) ;
    }
return ( param.A * weightedScalarProd(dens) );
}

double Tet::anisotropyEnergy(Tetra::prm const& param,const double u[DIM][NPI]) const
{
double dens[NPI];

double K = param.K;
double K3 = param.K3;

double uk00 = param.uk[0](0);
double uk01 = param.uk[0](1);
double uk02 = param.uk[0](2);
double uk10 = param.uk[1](0);
double uk11 = param.uk[1](1);
double uk12 = param.uk[1](2);
double uk20 = param.uk[2](0);
double uk21 = param.uk[2](1);
double uk22 = param.uk[2](2);   

for (int npi=0; npi<NPI; npi++)
    {
        // cosinus directeurs
    double al0=uk00*u[0][npi] + uk01*u[1][npi] + uk02*u[2][npi];
    double al1=uk10*u[0][npi] + uk11*u[1][npi] + uk12*u[2][npi];
    double al2=uk20*u[0][npi] + uk21*u[1][npi] + uk22*u[2][npi];
    
    dens[npi] = -K  * sq(al0) + K3 * (sq(al0) * sq(al1) + sq(al1) * sq(al2) + sq(al2) * sq(al0)); // uniaxial(K) + cubic(K3)
    }

return weightedScalarProd(dens);
}

void Tet::charges(std::function<Pt::pt3D (Nodes::Node)> getter,double *srcDen,int &nsrc,double Ms) const
{
double dudx[DIM][NPI], dudy[DIM][NPI], dudz[DIM][NPI];
interpolation(getter,dudx,dudy,dudz);
            
for (int j=0; j<NPI; j++, nsrc++)
    { srcDen[nsrc] = -Ms * ( dudx[0][j] + dudy[1][j] + dudz[2][j] ) * weight[j]; }
}

double Tet::demagEnergy(Tetra::prm const& param,const double dudx[DIM][NPI],const double dudy[DIM][NPI],const double dudz[DIM][NPI],const double phi[NPI]) const
{
double dens[NPI];
double Ms = nu0 * param.J;

for (int npi=0; npi<NPI; npi++)
    { dens[npi] = (dudx[0][npi] + dudy[1][npi] + dudz[2][npi])*phi[npi]; }
return ( -0.5*mu0*Ms*weightedScalarProd(dens) );
}

double Tet::zeemanEnergy(Tetra::prm const& param,double uz_drift,Pt::pt3D const& Hext,const double u[DIM][NPI]) const
{
double dens[NPI];
double Js = param.J;

for (int npi=0; npi<NPI; npi++)
    {
    dens[npi] = u[0][npi]*Hext(Pt::IDX_X) + u[1][npi]*Hext(Pt::IDX_Y) + u[2][npi]*Hext(Pt::IDX_Z) + uz_drift*Hext(Pt::IDX_Z);
    }
return ( -Js*weightedScalarProd(dens) );
}

/*
void Tet::projection(double A[3*N][3*N],  double B[3*N])
{
double P[2*N][3*N] = { {0} }; // P must be filled with zero
double PA[2*N][3*N]; // no need to initialize with zeros

for (int i=0; i<N; i++){
    Nodes::Node const& n = (*refNode)[ind[i]];
	P[i][i]  = n.ep.x();  P[i][N+i]  = n.ep.y();  P[i][2*N+i]  = n.ep.z();
	P[N+i][i]= n.eq.x();  P[N+i][N+i]= n.eq.y();  P[N+i][2*N+i]= n.eq.z();
    }
//Ap = (P*A)*trans(P);
//Bp = P*B;
tiny::mult<double,2*N,3*N,3*N>(P,A,PA);
tiny::direct_transposed_mult<double,2*N,3*N,2*N>(PA,P,Kp);
tiny::mult<double,2*N,3*N>(P,B,Lp);
}
*/

void Tet::assemblage_mat(write_matrix &K) const
{
for (int i=0; i < N; i++)
    {
    int i_= ind[i];             
        
    for (int j=0; j < N; j++)
        {
        int j_= ind[j];
        K(NOD+i_,j_) += Kp[i][j];      K(NOD+i_, NOD+j_) += Kp[  i][N+j];
        K(    i_,j_) += Kp[N+i][j];    K(    i_, NOD+j_) += Kp[N+i][N+j];
        }
    }    
}

double Tet::Jacobian(double J[DIM][DIM])
{
Pt::pt3D const & p0 = (*refNode)[ ind[0] ].p;
Pt::pt3D const & p1 = (*refNode)[ ind[1] ].p;
Pt::pt3D const & p2 = (*refNode)[ ind[2] ].p;
Pt::pt3D const & p3 = (*refNode)[ ind[3] ].p;
J[0][0] = p1.x()-p0.x(); J[0][1] = p2.x()-p0.x(); J[0][2] = p3.x()-p0.x();   
J[1][0] = p1.y()-p0.y(); J[1][1] = p2.y()-p0.y(); J[1][2] = p3.y()-p0.y();
J[2][0] = p1.z()-p0.z(); J[2][1] = p2.z()-p0.z(); J[2][2] = p3.z()-p0.z();
    
return Pt::det(J);
}

double Tet::calc_vol(void) const
{
Pt::pt3D const & p0 = (*refNode)[ ind[0] ].p;
Pt::pt3D const & p1 = (*refNode)[ ind[1] ].p;
Pt::pt3D const & p2 = (*refNode)[ ind[2] ].p;
Pt::pt3D const & p3 = (*refNode)[ ind[3] ].p;

return Pt::pTriple(p1-p0,p2-p0,p3-p0)/6.0;
}
