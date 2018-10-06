#include "facette.h"

#include "tiny.h"

using namespace Facette;

void Fac::integrales(std::vector<Facette::prm> const& params, std::vector<Node> const& myNode, mtl::dense_vector <double> &BE)
{
double Js = params[idxPrm].Js;
double Ks = params[idxPrm].Ks;

double uk00 = params[idxPrm].uk[0];
double uk01 = params[idxPrm].uk[1];
double uk02 = params[idxPrm].uk[2];

double Kbis = 2.0*Ks/Js;
/*-------------------- INTERPOLATION --------------------*/
double u_nod[3][N], u[3][NPI];

for (int i=0; i<N; i++){
    Node const& node = myNode[ ind[i] ];
    
    u_nod[0][i]   = node.u0.x();
    u_nod[1][i]   = node.u0.y();
    u_nod[2][i]   = node.u0.z();
    }

tiny::mult<double, 3, N, NPI> (u_nod, a, u);

/*-------------------------------------------------------*/

for (int npi=0; npi<NPI; npi++){
    double Kbis_ai;
    double w_uk0_u = weight[npi]*(uk00*u[0][npi] + uk01*u[1][npi] + uk02*u[2][npi]); 

    for (int i=0; i<N; i++){
        Kbis_ai = Kbis*a[i][npi];

        BE[i]    += (Kbis_ai* w_uk0_u*uk00);
        BE[N+i]  += (Kbis_ai* w_uk0_u*uk01);
        BE[2*N+i]+= (Kbis_ai* w_uk0_u*uk02);
	}
    }
}

void Fac::calc_surf(std::vector<Node> const& myNode)
{
int i0,i1,i2;
i0 = ind[0];    i1 = ind[1];    i2 = ind[2];

Pt::pt3D p0 = myNode[i0].p;
Pt::pt3D p1 = myNode[i1].p;
Pt::pt3D p2 = myNode[i2].p;

Pt::pt3D vec = (p1-p0)*(p2-p0);

double norm = vec.norm();//sqrt(vecx*vecx + vecy*vecy + vecz*vecz);
n = vec/norm;//fa.nx = vec.x()/norm,  fa.ny = vec.y()/norm,  fa.nz = vec.z()/norm;
surf = 0.5*norm;
}

