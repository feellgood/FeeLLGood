#include "linear_algebra.h"
#include "tetra.h"

using namespace std;

void LinAlgebra::integrales(Facette::Fac &fac, gmm::dense_matrix <double> &AE, vector <double> &BE)
{
//const int N   = Fac::N;
//const int NPI = Fac::NPI;
const int reg = fac.reg;

map <pair<string,int>,double> &param = settings.param;

double Js = abs(param[make_pair("Js",reg)])+EPSILON;//cout << ", Js=" << Js;
double Ks = param[make_pair("Ks",reg)];   	//cout << ", Ks=" << Ks;

double uk00 = param[make_pair("a1",reg)];  	//cout << ", a1=" << k0;
double uk01 = param[make_pair("a2",reg)];  	//cout << ", a2=" << k1;
double uk02 = param[make_pair("a3",reg)];  	//cout << ", a3=" << k2;

double Kbis = 2.0*Ks/Js;
/*-------------------- INTERPOLATION --------------------*/
double u_nod[3][Facette::N], u[3][Facette::NPI];

for (int i=0; i<Facette::N; i++){
    int i_= fac.ind[i];
    Node &node = fem.node[i_];
    
    u_nod[0][i]   = node.u0[0];
    u_nod[1][i]   = node.u0[1];
    u_nod[2][i]   = node.u0[2];
    }

tiny::mult<double, 3, Facette::N, Facette::NPI> (u_nod, fac.a, u);

/*-------------------------------------------------------*/

for (int npi=0; npi<Facette::NPI; npi++){
    double Kbis_ai;
    double w_uk0_u = fac.weight[npi]*(uk00*u[0][npi] + uk01*u[1][npi] + uk02*u[2][npi]); 

    for (int i=0; i<Facette::N; i++){
        Kbis_ai = Kbis*fac.a[i][npi];

        BE[i]    += (Kbis_ai* w_uk0_u*uk00);
        BE[Facette::N+i]  += (Kbis_ai* w_uk0_u*uk01);
        BE[2*Facette::N+i]+= (Kbis_ai* w_uk0_u*uk02);
	}
    }
}
