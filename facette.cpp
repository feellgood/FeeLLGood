#include "facette.h"
#include "tiny.h"

using namespace Facette;
using namespace Pt;

void Fac::integrales(std::vector<Facette::prm> const& params, Pt::pt3D (&BE)[N]) const
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

double Fac::anisotropyEnergy(Facette::prm const& param,const Pt::pt3D (&u)[NPI]) const
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

double Fac::demagEnergy(const Pt::pt3D (&u)[NPI],const double (&phi)[NPI]) const
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

double Fac::potential(std::function<Pt::pt3D (Nodes::Node)> getter, int i) const
{
int ii  = (i+1)%3;
int iii = (i+2)%3;

Nodes::Node const& node1 = refNode[ ind[i] ];
Nodes::Node const& node2 = refNode[ ind[ii] ];
Nodes::Node const& node3 = refNode[ ind[iii] ];

Pt::pt3D p1p2 = node2.p - node1.p;
Pt::pt3D p1p3 = node3.p - node1.p;

std::function<double(double)> f = [](double x){return sqrt(1.0+x*x);} ;

double b = p1p2.norm();
double t = Pt::pScal(p1p2,p1p3)/b; // carefull with t , if cos(p1p2,p1p3) is negative then t is negative
double _2s = 2.*calc_surf();
double h = _2s/b;

if (_2s<0)
    { std::cout << "facette surface is negative : surface is ill-oriented" << std::endl; exit(1);}

double c = (t-b)/h;
double r = h*f(t/h);// if _2s is positive it is the same as double r = sqrt( sq(_2s/b) + sq(t));
double log_1 = log( (c*t + h + f(c)*r) / (b*(c+f(c))) );
double xi = b*log_1/f(c);

Pt::pt3D n = calc_norm();
double s1 = Pt::pScal(getter(node1),n);
double s2 = Pt::pScal(getter(node2),n);
double s3 = Pt::pScal(getter(node3),n);

double pot = xi*s1 +( (xi*(h + c*t) - b*(r-b))*s2 + b*(r-b-c*xi)*s3 )*b/(_2s*(1+c*c));
return 0.5*Ms*pot;
}

void Fac::calcCorr(std::function<const Pt::pt3D (Nodes::Node)> getter,std::vector<double> &corr,Pt::pt3D (&u)[NPI]) const
{
// calc coord gauss
Pt::pt3D gauss[NPI];
        
interpolation<Pt::pt3D>(Nodes::get_p,gauss);
// calc corr node by node
Pt::pt3D n = calc_norm();

for (int i=0; i<N; i++)
    {
    const int i_ = ind[i];
    Pt::pt3D & p_i_ = refNode[ i_ ].p;	      
    for (int j=0; j<NPI; j++)
        {
        corr[i_]-= Ms*pScal(u[j],n)*weight(j)/Pt::dist(p_i_, gauss[j]);
        }
    corr[i_]+= potential(getter,i);
    }
}


