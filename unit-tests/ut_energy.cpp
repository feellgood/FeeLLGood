#define BOOST_TEST_MODULE tetraTest

#include <boost/test/unit_test.hpp>
#include <iomanip>
#include <random>

#include "tetra.h"
#include "facette.h"
#include "tiny.h"
#include "ut_tools.h"
#include "ut_config.h"

BOOST_AUTO_TEST_SUITE(ut_energy)

/*---------------------------------------*/
/* test: check if demag energy computed through demag field Hd gives same energy as tet+fac demag energy computed from scaler potential phi   */
/*---------------------------------------*/

BOOST_AUTO_TEST_CASE(demagEnergy, *boost::unit_test::tolerance(UT_TOL))
    {
    const int nbNod = 4;
    std::vector<Nodes::Node> node;
    dummyNodes<nbNod>(node);

    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);

    for (int i = 0; i < nbNod; i++)
        {
        node[i].d[1].u = rand_vec3d(M_PI * distrib(gen), 2 * M_PI * distrib(gen));
        node[i].d[1].phi = distrib(gen);
        }

    // carefull with indices (starting from 1)
    Tetra::Tet t(node, 0, {1, 2, 3, 4});

    double t_dadx[Tetra::N][Tetra::NPI],t_dady[Tetra::N][Tetra::NPI],t_dadz[Tetra::N][Tetra::NPI];
    
    for(int i=0;i<Tetra::N;i++)
        for(int j=0;j<Tetra::NPI;j++)
            {
            t_dadx[i][j] = t.dadx(i,j);
            t_dady[i][j] = t.dady(i,j);
            t_dadz[i][j] = t.dadz(i,j);
            }
    // ref code (with minimal adaptations of MuMag::energy method in file MuMag_energy.cc of
    // src_Tube_scalfmm_thiaville_ec_mu_oersted_thiele_dyn20180903.tgz )

    double Ms = distrib(gen);
    double Js = mu0 * Ms;
    double u_nod[3][Tetra::N], u[3][Tetra::NPI];
    double negphi_nod[Tetra::N], Hdx[Tetra::NPI], Hdy[Tetra::NPI], Hdz[Tetra::NPI];
    
    for (int ie=0; ie<Tetra::N; ie++)
        {
        int i= t.ind[ie];
        for (int dim=0; dim<3; dim++)
            { u_nod[dim][ie] = node[i].d[1].u[dim]; }
        negphi_nod[ie] = -node[i].d[1].phi;
        }
    
    tiny::mult<double, 3, Tetra::N, Tetra::NPI> (u_nod, Tetra::a, u);
    tiny::transposed_mult<double, Tetra::N, Tetra::NPI> (negphi_nod, t_dadx, Hdx);
	tiny::transposed_mult<double, Tetra::N, Tetra::NPI> (negphi_nod, t_dady, Hdy);
	tiny::transposed_mult<double, Tetra::N, Tetra::NPI> (negphi_nod, t_dadz, Hdz);
    
    double dens[Tetra::NPI];
    for (int npi=0; npi<Tetra::NPI; npi++) {
        dens[npi] = -0.5*Js* (u[0][npi]*Hdx[npi] + u[1][npi]*Hdy[npi] + u[2][npi]*Hdz[npi]);
        }
    double Edemag(0);
    for (int npi=0; npi<Tetra::NPI; npi++) { Edemag += dens[npi]*t.weight[npi]; }
    // end ref code
    
    //code to test
    Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> _u,dudx,dudy,dudz;
    Eigen::Matrix<double,Tetra::NPI,1> _phi;
    t.interpolation(Nodes::get_u<Nodes::NEXT>, _u, dudx, dudy, dudz);
    t.interpolation(Nodes::get_phi<Nodes::NEXT>, _phi);
    Tetra::prm param;
    param.J = Js;
    double result_to_test = t.demagEnergy(param, dudx, dudy, dudz, _phi);
    std::cout << "E_demag(vol)= " << result_to_test << std::endl;
    
    const int ia = t.ind[0] + 1;
    const int ib = t.ind[1] + 1;
    const int ic = t.ind[2] + 1;
    const int id = t.ind[3] + 1;
    
    std::vector<Facette::Fac> fa;
    fa.push_back( Facette::Fac(node, nbNod, 0, {ia, ic, ib} ));
    fa.push_back( Facette::Fac(node, nbNod, 0, {ib, ic, id} ));
    fa.push_back( Facette::Fac(node, nbNod, 0, {ia, id, ic} ));
    fa.push_back( Facette::Fac(node, nbNod, 0, {ia, ib, id} ));
    std::for_each(fa.begin(),fa.end(),[Ms](Facette::Fac &f){ f.Ms = Ms;});
    
    std::for_each(fa.begin(),fa.end(), [&result_to_test](Facette::Fac &f)
        {
        Eigen::Matrix<double,Facette::NPI,1> phi_fa;
        Eigen::Matrix<double,Nodes::DIM,Facette::NPI> u_fa;
        f.interpolation(Nodes::get_u<Nodes::NEXT>, u_fa);
        f.interpolation(Nodes::get_phi<Nodes::NEXT>, phi_fa);
        result_to_test += f.demagEnergy(u_fa, phi_fa);
        });
    // end code to test
    std::cout << "val ref Edemag= " << std::setprecision(15) << Edemag << std::endl;
    std::cout << "val to test Edemag= " << result_to_test << std::endl;
    BOOST_TEST( Nodes::sq(result_to_test - Edemag) == 0.0 );
    }

BOOST_AUTO_TEST_SUITE_END()
