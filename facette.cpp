#include "facette.h"

#include "tiny.h"

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
/*-------------------- INTERPOLATION --------------------*/
double u_nod[3][N], u[3][NPI];

for (int i=0; i<N; i++){
    Node const& node = (*refNode)[ ind[i] ];
    
    u_nod[0][i]   = node.u0.x();
    u_nod[1][i]   = node.u0.y();
    u_nod[2][i]   = node.u0.z();
    }

tiny::mult<double, 3, N, NPI> (u_nod, a, u);

/*-------------------------------------------------------*/

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


void Fac::projection(gmm::dense_matrix <double> const& A,  std::vector <double> const& B,
           gmm::dense_matrix <double> &Ap, std::vector <double> &Bp) const
{
thread_local gmm::dense_matrix <double> P(2*N,3*N);
thread_local gmm::dense_matrix <double> PA(2*N,3*N);
//mtl::mat::set_to_zero(P);

for (int i=0; i<N; i++){
    Node const& n = (*refNode)[ind[i]];
    P(i,i)  = n.ep.x();  P(i,N+i)  = n.ep.y();  P(i,2*N+i)  = n.ep.z();
    P(N+i,i)= n.eq.x();  P(N+i,N+i)= n.eq.y();  P(N+i,2*N+i)= n.eq.z();
    }

//Ap = (P*A)*trans(P);
//Bp = P*B;
gmm::mult(P,A,PA);
gmm::mult(PA, gmm::transposed(P), Ap);

gmm::mult(P,B,Bp);

}

void Fac::assemblage(const int NOD,gmm::dense_matrix <double> const& Ke, std::vector <double> const& Le,write_matrix &K,write_vector &L) const
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
        L[NOD+i_] += Le[i];//L[NOD+i_]+= Le[  i];
        L[i_] += Le[N+i];//L[    i_]+= Le[N+i];
        
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

