#define BOOST_TEST_MODULE time_integrationTest

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <random>

#include "ut_config.h" // for tolerance UT_TOL macro
#include "pt3D.h"
#include "time_integration.h"

BOOST_AUTO_TEST_SUITE(ut_time_int)

/*---------------------------------------*/
/* minus one test: check boost is fine   */
/*---------------------------------------*/

BOOST_AUTO_TEST_CASE(Stupid)
{
float x = 1.0;
BOOST_CHECK(x != 0.0f);
}

/*-----------------------------------------------------*/
/* zero lvl tests : direct elementary member functions */ 
/*-----------------------------------------------------*/

BOOST_AUTO_TEST_CASE(timing_constructor, * boost::unit_test::tolerance(UT_TOL))
{
unsigned sd = my_seed();
std::mt19937 gen(my_seed());
std::uniform_real_distribution<> distrib(0.0,1.0);


double a= distrib(gen);
double b= distrib(gen);

timing prm_t = timing(0.0,1.0,std::min(a,b),std::max(a,b));

double dt = prm_t.get_dt();

std::cout << "test timing constructor" << std::endl;
if (!DET_UT) std::cout << "seed =" << sd << std::endl;
BOOST_TEST( dt == sqrt(a*b) );
}


/*---------------------------------------*/
/* first lvl tests : nested calculus,... */
/*---------------------------------------*/

/*---------------------------------------*/
/* second lvl tests : pure mathematics   */
/*---------------------------------------*/

BOOST_AUTO_TEST_CASE(calc_alpha_eff, * boost::unit_test::tolerance(UT_TOL))
{
unsigned sd = my_seed();
std::mt19937 gen(sd);
std::uniform_real_distribution<> distrib(0.0,1.0);

double alpha_LLG = distrib(gen);

double X= distrib(gen);

timing prm_t = timing(0.0,1e-8,1e-14,1e-9);

// ref code (with minimal adaptations of MuMag_Integrales.cc of src_Tube_scalfmm_thiaville_ec_mu_oersted_thiele_dyn20180903.tgz )
double dt = prm_t.get_dt();
double alpha = alpha_LLG;
double uHeff = X;
double r = 0.1;                             
double M = 2.*alpha*r/dt;                   

double alfa=0.;
    if (uHeff>0.){
       if (uHeff>M) alfa=alpha+dt/2.*M;
       else alfa=alpha+dt/2.*uHeff;
       }
    else{
       if (uHeff<-M) alfa=alpha/(1.+dt/(2.*alpha)*M);
       else alfa=alpha/(1.-dt/(2.*alpha)*uHeff);
       }

// end ref code


// code to check
double my_alpha = prm_t.calc_alpha_eff(alpha_LLG,X);
// end code to check

std::cout << "test that calc_alpha_eff gives correct effective damping parameter relatively to ref source code" << std::endl;
if (!DET_UT) std::cout << "seed =" << sd << std::endl;
BOOST_TEST( my_alpha == alfa );
}


BOOST_AUTO_TEST_SUITE_END()
