/**
  Elementary matrix Calculation for a tetrahedron element 
 */ 

#include <set>

#include "tetra.h"
#include "pt3D.h"
#include "tiny.h"
#include "matBlocDiag.h"
#include "time_integration.h"

#include "facette.h"

using namespace Tetra;
using namespace Pt;


void Tet::lumping(int const& npi,double alpha_eff,double prefactor, double AE[3*N][3*N]) const
{
const double w = weight[npi];

for (int i=0; i<N; i++)
    {
    const double ai_w = w*a[i][npi];
    const pt3D ai_w_u0 = ai_w*Nodes::get_u0(refNode[ ind[i] ]);

    AE[    i][    i] +=  alpha_eff * ai_w;
    AE[  N+i][  N+i] +=  alpha_eff * ai_w;
    AE[2*N+i][2*N+i] +=  alpha_eff * ai_w;

    AE[0*N+i][2*N+i] += ai_w_u0(Pt::IDX_Y);
    AE[0*N+i][1*N+i] -= ai_w_u0(Pt::IDX_Z);
    AE[1*N+i][0*N+i] += ai_w_u0(Pt::IDX_Z);
    AE[1*N+i][2*N+i] -= ai_w_u0(Pt::IDX_X);
    AE[2*N+i][1*N+i] += ai_w_u0(Pt::IDX_X);
    AE[2*N+i][0*N+i] -= ai_w_u0(Pt::IDX_Y);

    for (int j=0; j<N; j++)
        {
        double contrib = w*prefactor*(dadx[i][npi]*dadx[j][npi] + dady[i][npi]*dady[j][npi] + dadz[i][npi]*dadz[j][npi]);
            
        AE[i][j]        +=  contrib;
        AE[N+i][N+j]    +=  contrib;
        AE[2*N+i][2*N+j] +=  contrib;
        }
    }
}

void Tet::build_BE(int const& npi, Pt::pt3D const & Ht, Pt::pt3D const & Heff,double alpha, double beta, double Abis, double s_dt, double Uz, double Vz,
                   Pt::pt3D U[NPI], Pt::pt3D V[NPI],
            Pt::pt3D dUdx[NPI], Pt::pt3D dUdy[NPI], Pt::pt3D dUdz[NPI], Pt::pt3D dVdx[NPI], Pt::pt3D dVdy[NPI], Pt::pt3D dVdz[NPI], Pt::pt3D BE[N]) const
{// the artificial drift from eventual recentering is only along z !! it should be generalized
for (int i=0; i<N; i++)
    {
    const double ai_w = weight[npi]*a[i][npi];
    BE[i] -= weight[npi]*Abis*(dadx[i][npi]*dUdx[npi] + dady[i][npi]*dUdy[npi] + dadz[i][npi]*dUdz[npi]);
    BE[i] += ai_w*(alpha*Vz - beta*Uz)*(dUdz[npi] + dVdz[npi]*s_dt) ;
    BE[i] += ai_w*(Heff + Ht*s_dt + (Vz-Uz)*(U[npi]*(dUdz[npi]+dVdz[npi]*s_dt) +V[npi]*dUdz[npi]*s_dt) );
    }
}

void Tet::integrales(std::vector<Tetra::prm> const& params, timing const& prm_t,Pt::pt3D const& Hext,double Vz, double AE[3*N][3*N], Pt::pt3D BE[N]) const
{
double alpha_LLG = params[idxPrm].alpha_LLG;
pt3D ex = params[idxPrm].ex;
pt3D ey = params[idxPrm].ey;
pt3D ez = params[idxPrm].ez;

double Abis = 2.0*(params[idxPrm].A)/(params[idxPrm].J);
double Kbis = 2.0*(params[idxPrm].K)/(params[idxPrm].J);
double K3bis = 2.0*(params[idxPrm].K3)/(params[idxPrm].J);
const double s_dt = THETA*prm_t.get_dt();//theta from theta scheme in config.h.in

/*-------------------- INTERPOLATION --------------------*/
pt3D Hd[NPI], dUdx[NPI], dUdy[NPI], dUdz[NPI];
pt3D Hv[NPI], dVdx[NPI], dVdy[NPI], dVdz[NPI];
pt3D U[NPI],V[NPI];

interpolation(Nodes::get_u0,U,dUdx,dUdy,dUdz);
interpolation(Nodes::get_v0,V,dVdx,dVdy,dVdz);
interpolation(Nodes::get_phi0,Hd);
interpolation(Nodes::get_phiv0,Hv);

for (int npi=0; npi<NPI; npi++)
    {
    pt3D uk_u = pt3D(pScal(ex,U[npi]), pScal(ey,U[npi]), pScal(ez,U[npi]));
    pt3D uk_v = pt3D(pScal(ex,V[npi]), pScal(ey,V[npi]), pScal(ez,V[npi]));

    Pt::pt3D cube_uk_u = directCube(uk_u);
    Pt::pt3D uk_uuu = uk_u - cube_uk_u;

    double uHeff = -Abis*(norme2(dUdx[npi]) + norme2(dUdy[npi]) + norme2(dUdz[npi])); 
	uHeff +=  pScal(U[npi], Hext + Hd[npi]) + Kbis*sq(pScal(params[idxPrm].uk,U[npi])) - K3bis*pScal(uk_u,uk_uuu);

    lumping(npi, prm_t.calc_alpha_eff(alpha_LLG,uHeff), prm_t.prefactor*s_dt*Abis, AE);
    
    pt3D Ht = Hv[npi];
    Ht += Kbis*pScal(params[idxPrm].uk,V[npi])*params[idxPrm].uk; 
    Ht += -K3bis*pDirect(uk_v,uk_u-3*cube_uk_u);

    pt3D Heff = Kbis*pScal(params[idxPrm].uk,U[npi])*params[idxPrm].uk;
    Heff += -K3bis*( uk_uuu(0)*ex + uk_uuu(1)*ey + uk_uuu(2)*ez ) + Hd[npi] + Hext;

    build_BE(npi, Ht, Heff, alpha_LLG, params[idxPrm].beta_sc, Abis, s_dt, params[idxPrm].Uz, Vz, U, V, dUdx, dUdy, dUdz, dVdx, dVdy, dVdz, BE);
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

for (int npi=0; npi<NPI; npi++)
    {
    Pt::pt3D m = Pt::pt3D(u[0][npi],u[1][npi],u[2][npi]);
    // uniaxial magnetocrystalline anisotropy constant K, anisotropy axis uk 
    dens[npi] = -param.K*sq( Pt::pScal( param.uk, m) );
        
        // cosinus directeurs
    double al0= Pt::pScal( param.ex, m );
    double al1= Pt::pScal( param.ey, m );
    double al2= Pt::pScal( param.ez, m );
    
    dens[npi] += param.K3*(sq(al0*al1) + sq(al1*al2) + sq(al2*al0)); // cubic anisotropy (K3)
    }

return weightedScalarProd(dens);
}

void Tet::charges(std::function<Pt::pt3D (Nodes::Node)> getter,std::vector<double> &srcDen,int &nsrc,double Ms) const
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

double Tet::zeemanEnergy(Tetra::prm const& param,double uz_drift,Pt::pt3D const& Hext,double const u[DIM][NPI]) const
{
double dens[NPI];

for (int npi=0; npi<NPI; npi++)
    {
    Pt::pt3D u_npi= Pt::pt3D(u[0][npi],u[1][npi],u[2][npi]);
    dens[npi] = Pt::pScal(u_npi,Hext) + uz_drift*Hext(Pt::IDX_Z);
    }
return ( -param.J*weightedScalarProd(dens) );
}

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
Pt::pt3D const & p0 = refNode[ ind[0] ].p;
Pt::pt3D const & p1 = refNode[ ind[1] ].p;
Pt::pt3D const & p2 = refNode[ ind[2] ].p;
Pt::pt3D const & p3 = refNode[ ind[3] ].p;
J[0][0] = p1.x()-p0.x(); J[0][1] = p2.x()-p0.x(); J[0][2] = p3.x()-p0.x();   
J[1][0] = p1.y()-p0.y(); J[1][1] = p2.y()-p0.y(); J[1][2] = p3.y()-p0.y();
J[2][0] = p1.z()-p0.z(); J[2][1] = p2.z()-p0.z(); J[2][2] = p3.z()-p0.z();
    
return Pt::det(J);
}

double Tet::calc_vol(void) const
{
Pt::pt3D const & p0 = refNode[ ind[0] ].p;
Pt::pt3D const & p1 = refNode[ ind[1] ].p;
Pt::pt3D const & p2 = refNode[ ind[2] ].p;
Pt::pt3D const & p3 = refNode[ ind[3] ].p;

return Pt::pTriple(p1-p0,p2-p0,p3-p0)/6.0;
}

std::set<Facette::Fac> Tet::ownedFac() const
{
std::set<Facette::Fac> s;

int ia=ind[0];int ib=ind[1];int ic=ind[2];int id=ind[3];

s.insert( Facette::Fac(nullptr,0,reg,idxPrm,ia,ic,ib) );
s.insert( Facette::Fac(nullptr,0,reg,idxPrm,ib,ic,id) );
s.insert( Facette::Fac(nullptr,0,reg,idxPrm,ia,id,ic) );
s.insert( Facette::Fac(nullptr,0,reg,idxPrm,ia,ib,id) );

return s;
}
