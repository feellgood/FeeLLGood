#include "facette.h"
#include "tiny.h"

using namespace Facette;
using namespace Pt;

void Fac::integrales(std::vector<Facette::prm> const& params, Pt::pt3D BE[N]) const
{
Pt::pt3D const & uk = params[idxPrm].uk;
double Kbis = 2.0*(params[idxPrm].Ks)/Ms; // carefull Ms of the facette here (could lead to div by zero ?)

Pt::pt3D u[NPI];
interpolation<Pt::pt3D>(Nodes::get_u0,u);

for (int npi=0; npi<NPI; npi++)
	{
    double w_uk_u = weight(npi)*pScal(uk,u[npi]);

    for (int i=0; i<N; i++) { BE[i] += Kbis*a[i][npi]*w_uk_u*uk; }
    }
}

double Fac::anisotropyEnergy(Facette::prm const& param,const Pt::pt3D u[NPI]) const
{//surface Neel anisotropy (uk is a uniaxial easy axis)
double dens[NPI];
for (int npi=0; npi<NPI; npi++)
    { dens[npi] = -param.Ks*Pt::sq(pScal(param.uk,u[npi])); }
return weightedScalarProd(dens);
}

void Fac::charges(std::function<Pt::pt3D (Nodes::Node)> getter,std::vector<double> &srcDen,std::vector<double> &corr,int &nsrc) const
{
Pt::pt3D u[NPI];
interpolation<Pt::pt3D>(getter,u);
Pt::pt3D n = calc_norm();

for (int j=0; j<NPI; j++, nsrc++)
    { srcDen[nsrc] = Ms * weight(j) * pScal(u[j],n); }

calcCorr(getter,corr,u);
}

double Fac::demagEnergy(const Pt::pt3D u[NPI],const double phi[NPI]) const
{
double q[NPI];
Pt::pt3D n = calc_norm();

for (int npi=0; npi<NPI; npi++)
        { q[npi] = Ms * pScal(u[npi],n); }
double dens[NPI];
for (int npi=0; npi<NPI; npi++)
    { dens[npi] = 0.5*mu0*q[npi]*phi[npi]; }
return weightedScalarProd(dens);
}

void Fac::assemblage_mat(write_matrix &K) const
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

Pt::pt3D Fac::calc_norm(void) const
{
Pt::pt3D p0 = refNode[ ind[0] ].p;
Pt::pt3D p1 = refNode[ ind[1] ].p;
Pt::pt3D p2 = refNode[ ind[2] ].p;

Pt::pt3D n = (p1-p0)*(p2-p0);
n.normalize();
return n;
}

double Fac::calc_surf(void) const
{
Pt::pt3D p0 = refNode[ ind[0] ].p;
Pt::pt3D p1 = refNode[ ind[1] ].p;
Pt::pt3D p2 = refNode[ ind[2] ].p;

Pt::pt3D vec = (p1-p0)*(p2-p0);

return 0.5*vec.norm();
}

double Fac::potential(std::function<Pt::pt3D (Nodes::Node)> getter, int i) const
{
int ii  = (i+1)%3;
int iii = (i+2)%3;

Nodes::Node const& node1 = refNode[ ind[i] ];
Nodes::Node const& node2 = refNode[ ind[ii] ];
Nodes::Node const& node3 = refNode[ ind[iii] ];

Pt::pt3D p1p2 = node2.p - node1.p;
Pt::pt3D p1p3 = node3.p - node1.p;

double b = p1p2.norm();
double t = Pt::pScal(p1p2,p1p3)/b;
double h = 2.*calc_surf()/b;
double a = t/h;  
double c = (t-b)/h;
double cc1 = 1.0 + sq(c);
double r = sqrt( sq(h) + sq(c*h+b));
double log_1 = log( (cc1*h + c*b + sqrt(cc1)*r) / (b*(c+sqrt(cc1))) );
double log_2 = log(c*h+b+r);

Pt::pt3D n = calc_norm();

double s1 = Pt::pScal(getter(node1),n);
double s2 = Pt::pScal(getter(node2),n);
double s3 = Pt::pScal(getter(node3),n);

double j = (s2-s1)/b;
double k = t/b/h*(s1-s2) + (s3-s1)/h;

double pot1 = b*(b/pow(cc1,1.5)*log_1 + c*(r-b)/cc1) + h*(r - sqrt(a*a+1)*h);

double pot2 = b*(-c*b/pow(cc1,1.5)*log_1 + (r-b)/cc1) + sq(h)*(log_2 - 0.5);

double pot3 = h*(log_2 - 1.0) + b/sqrt(cc1)*log_1;

double pot = 0.5*(j*pot1 + k*pot2) + s1*pot3 + h*(k*h/2.+s1)*(1-log(h*(a+sqrt(a*a+1)))) - 0.25*k*h*h;

return Ms*pot;
}

void Fac::calcCorr(std::function<const Pt::pt3D (Nodes::Node)> getter,std::vector<double> &corr,Pt::pt3D u[NPI]) const
{
// calc coord gauss
Pt::pt3D gauss[NPI];
        
interpolation<Pt::pt3D>(Nodes::get_p,gauss);
// calc corr node by node
Pt::pt3D n = calc_norm();

for (int i=0; i<N; i++)
    {
    const int i_ = ind[i];
    const Pt::pt3D p_i_ = refNode[ i_ ].p;	      
    for (int j=0; j<NPI; j++)
        {
        corr[i_]-= Ms*pScal(u[j],n)*weight(j)/Pt::dist(p_i_, gauss[j]);
        }
    corr[i_]+= potential(getter,i);
    }
}


