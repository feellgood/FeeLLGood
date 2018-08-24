#include "facette.h"

#include "tiny.h"

using namespace Facette;

void Fac::integrales(Settings &mySets, std::vector<Node> &myNode, std::vector <double> &BE)
{
std::map <std::pair<std::string,int>,double> &param = mySets.param;

double Js = abs(param[std::make_pair("Js",reg)])+EPSILON;//cout << ", Js=" << Js;
double Ks = param[std::make_pair("Ks",reg)];   	//cout << ", Ks=" << Ks;

double uk00 = param[std::make_pair("a1",reg)];  	//cout << ", a1=" << k0;
double uk01 = param[std::make_pair("a2",reg)];  	//cout << ", a2=" << k1;
double uk02 = param[std::make_pair("a3",reg)];  	//cout << ", a3=" << k2;

double Kbis = 2.0*Ks/Js;
/*-------------------- INTERPOLATION --------------------*/
double u_nod[3][N], u[3][NPI];

for (int i=0; i<N; i++){
    Node &node = myNode[ ind[i] ];
    
    u_nod[0][i]   = node.u0[0];
    u_nod[1][i]   = node.u0[1];
    u_nod[2][i]   = node.u0[2];
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
