#include "facette.h"

using namespace Facette;

void Fac::integrales(std::vector<Facette::prm> const& params, std::vector <double> &BE) const
{
double Js = params[idxPrm].Js;
double Ks = params[idxPrm].Ks;

if(Ks!=0.)
{
double uk00 = params[idxPrm].uk[0];
double uk01 = params[idxPrm].uk[1];
double uk02 = params[idxPrm].uk[2];
double Kbis = 2.0*Ks/Js;

double u[DIM][NPI];
interpolation(Nodes::get_u0,u);

    for (int npi=0; npi<NPI; npi++)
        {
        double Kbis_ai;
        double w_uk0_u = weight[npi]*(uk00*u[0][npi] + uk01*u[1][npi] + uk02*u[2][npi]); 

        for (int i=0; i<N; i++)
            {
            Kbis_ai = Kbis*a[i][npi];
            BE[i]    += (Kbis_ai* w_uk0_u*uk00);
            BE[N+i]  += (Kbis_ai* w_uk0_u*uk01);
            BE[2*N+i]+= (Kbis_ai* w_uk0_u*uk02);
            }
        }
    }
}

double Fac::anisotropyEnergy(Facette::prm const& param,const double u[DIM][NPI]) const
{
double uk00 = param.uk[0];
double uk01 = param.uk[1];
double uk02 = param.uk[2];    

double dens[NPI];
for (int npi=0; npi<NPI; npi++)
    {// cosinus directeurs
    double al0=uk00*u[0][npi] + uk01*u[1][npi] + uk02*u[2][npi];
    dens[npi] = -param.Ks*al0*al0;      // uniaxe
    }
return weightedScalarProd(dens);
}

double Fac::demagEnergy(const double u[DIM][NPI],const double phi[NPI]) const
{
double q[NPI];

for (int npi=0; npi<NPI; npi++)
        { q[npi] = Ms * (u[0][npi]*n.x() + u[1][npi]*n.y() + u[2][npi]*n.z()); }
double dens[NPI];
for (int npi=0; npi<NPI; npi++)
    { dens[npi] = 0.5*mu0*q[npi]*phi[npi]; }
return weightedScalarProd(dens);
}

void Fac::projection(gmm::dense_matrix <double> const& A, std::vector <double> const& B)
{
thread_local gmm::dense_matrix <double> P(2*N,3*N);
thread_local gmm::dense_matrix <double> PA(2*N,3*N);

for (int i=0; i<N; i++){
    Nodes::Node const& n = (*refNode)[ind[i]];
    P(i,i)  = n.ep.x();  P(i,N+i)  = n.ep.y();  P(i,2*N+i)  = n.ep.z();
    P(N+i,i)= n.eq.x();  P(N+i,N+i)= n.eq.y();  P(N+i,2*N+i)= n.eq.z();
    }

//Ap = (P*A)*trans(P); Bp = P*B;
gmm::mult(P,A,PA);
gmm::mult(PA, gmm::transposed(P), Ksp);

gmm::mult(P,B,Lsp);

}

void Fac::assemblage(write_matrix &K,write_vector &L) const
{
    for (int i=0; i < N; i++)
        {
        int i_= ind[i];             
        
        for (int j=0; j < N; j++)
            {
            int j_= ind[j];
            K(NOD+i_,j_) += Ksp(i,j);      K(NOD+i_, NOD+j_) += Ksp(  i,N+j);
            K(    i_,j_) += Ksp(N+i,j);    K(    i_, NOD+j_) += Ksp(N+i,N+j);
            }
        L[NOD+i_] += Lsp[i];
        L[i_] += Lsp[N+i];
        
        }
    
}

void Fac::calc_surf(void)
{
int i0,i1,i2;
i0 = ind[0];    i1 = ind[1];    i2 = ind[2];

Pt::pt3D p0 = (*refNode)[i0].p;
Pt::pt3D p1 = (*refNode)[i1].p;
Pt::pt3D p2 = (*refNode)[i2].p;

Pt::pt3D vec = (p1-p0)*(p2-p0);

double norm = vec.norm();//sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
n = vec/norm;//fa.nx = vec.x()/norm,  fa.ny = vec.y()/norm,  fa.nz = vec.z()/norm;
surf = 0.5*norm;
}

double Fac::potential(std::function<Pt::pt3D (Nodes::Node)> getter, int i) const
{
int ii  = (i+1)%3;
int iii = (i+2)%3;

int i_,ii_,iii_;
i_= ind[i];  ii_= ind[ii];  iii_= ind[iii];

Nodes::Node const& node1 = (*refNode)[i_];
Nodes::Node const& node2 = (*refNode)[ii_];
Nodes::Node const& node3 = (*refNode)[iii_];

Pt::pt3D p1p2 = node2.p - node1.p;
Pt::pt3D p1p3 = node3.p - node1.p;

double b = p1p2.norm();
double t = Pt::pScal(p1p2,p1p3);
double h = 2.*surf;
 t/=b;  h/=b;
 double a = t/h;  double c = (t-b)/h;

double s1 = Pt::pScal(getter(node1),n);
double s2 = Pt::pScal(getter(node2),n);
double s3 = Pt::pScal(getter(node3),n);

// formule dégueu à reécrire

double l = s1;
double j = (s2-s1)/b;
double k = t/b/h*(s1-s2) + (s3-s1)/h;

double cc1 = c*c+1;
double r = sqrt(h*h + (c*h+b)*(c*h+b));
double ll = log( (cc1*h + c*b + sqrt(cc1)*r) / (b*(c+sqrt(cc1))) );

double pot1, pot2, pot3, pot;
pot1 = b*b/pow(cc1,1.5)*ll + c*b*r/cc1 + h*r - c*b*b/cc1 - sqrt(a*a+1)*h*h;
pot1*= j/2.;

pot2 = -c*b*b/pow(cc1,1.5)*ll + b*r/cc1 - h*h/2. + h*h*log(c*h+b+r) - b*b/cc1;
pot2*= k/2.;

pot3 = h*log(c*h+b+r) - h + b/sqrt(cc1)*ll;
pot3*= l;

pot = pot1 + pot2 + pot3 + h*(k*h/2.+l)*(1-log(h*(a+sqrt(a*a+1)))) - k*h*h/4.;

return Ms*pot;
}

void Fac::calcCorr(std::function<const Pt::pt3D (Nodes::Node)> getter,double *corr,double u[DIM][NPI]) const
{
// calc coord gauss
double gauss[DIM][NPI];
        
interpolation(Nodes::get_p,gauss);
// calc corr node by node
for (int i=0; i<N; i++)
    {
    int i_ = ind[i];
    Pt::pt3D p_i_ = (*refNode)[i_].p;	      
    for (int j=0; j<NPI; j++)
        {
        Pt::pt3D pg = Pt::pt3D(gauss[Pt::IDX_X][j], gauss[Pt::IDX_Y][j], gauss[Pt::IDX_Z][j]);
        double sj = Ms* ( u[0][j]*n.x() + u[1][j]*n.y() + u[2][j]*n.z() );
        corr[i_]-= sj*weight[j]/Pt::dist(p_i_,pg);
        }
    corr[i_]+= potential(getter,i);//potential<Hv>(fem.node, fac, i);
    }
}


