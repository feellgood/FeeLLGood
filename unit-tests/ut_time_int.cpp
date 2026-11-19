#define BOOST_TEST_MODULE time_integrationTest

#include <boost/test/unit_test.hpp>

#include <cstdlib>
#include <iostream>
#include <random>

#include "tetra.h"
#include "time_integration.h"
#include "ut_config.h"  // for tolerance UT_TOL macro

BOOST_AUTO_TEST_SUITE(ut_time_int)

/*-----------------------------------------------------*/
/* zero lvl tests : direct elementary member functions */
/*-----------------------------------------------------*/

BOOST_AUTO_TEST_CASE(timing_constructor, *boost::unit_test::tolerance(UT_TOL))
    {
    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);

    double a = distrib(gen);
    double b = distrib(gen);

    timing prm_t = timing(1.0, std::min(a, b), std::max(a, b));

    double dt = prm_t.get_dt();

    std::cout << "test timing constructor" << std::endl;
    if (!DET_UT) std::cout << "seed =" << sd << std::endl;
    BOOST_TEST(dt == sqrt(a * b));
    }

/*---------------------------------------*/
/* second lvl tests : pure mathematics   */
/*---------------------------------------*/

BOOST_AUTO_TEST_CASE(calc_alpha_eff, *boost::unit_test::tolerance(UT_TOL))
    {
    /*
    this unit test checks the four branches of the piecewise function defining effective alpha
    ref code comes from research version of st-feellgood MuMag_Integrales.cc (with adaptations)
     */
    unsigned sd = my_seed();
    std::mt19937 gen(sd);
    std::uniform_real_distribution<> distrib(0.0, 1.0);
    double alpha_LLG = distrib(gen);
    timing prm_t = timing(1e-8, 1e-14, 1e-9);

    const double dt = prm_t.get_dt();
    const double r = 0.1;
    double M = 2. * alpha_LLG * r / dt;

    Eigen::Matrix<double,Tetra::NPI,1> alpha_g;
    alpha_g.setConstant(alpha_LLG);

    // x will be a number equal to a dot product of magnetization and Heff
    std::function<double(double)> func = [M,alpha_LLG,dt](const double x)
        {
        //ref code
        double y;
        if (x > 0.)
            {
            if (x > M)
                y = alpha_LLG + dt / 2. * M; //case 1
            else
                y = alpha_LLG + dt / 2. * x; //case 2
            }
        else
            {
            if (x < -M)
                y = alpha_LLG / (1. + dt / (2. * alpha_LLG) * M); //case 3
            else
                y = alpha_LLG / (1. - dt / (2. * alpha_LLG) * x); //case 4
            }
        //end ref code
        return y;
        }; // end lambda

    std::cout << "test that calc_alpha_eff gives correct effective damping parameter\n";
    if (!DET_UT) std::cout << "seed =" << sd << std::endl;
    Eigen::Matrix<double,Tetra::NPI,1> my_alpha;
    double uHeff = distrib(gen);
    Eigen::Matrix<double,Tetra::NPI,1> _uHeff;
    _uHeff.setConstant(uHeff);
    double a_eff = func(uHeff); // case 2
    my_alpha = prm_t.get_alpha_eff<Tetra::NPI>(alpha_g, _uHeff);
    for(int npi=0;npi<Tetra::NPI;npi++) { BOOST_TEST(my_alpha(npi) == a_eff); }

    uHeff *= -1.0;
    a_eff = func(uHeff); // case 4
    _uHeff.setConstant(uHeff);
    my_alpha = prm_t.get_alpha_eff<Tetra::NPI>(alpha_g, _uHeff );
    for(int npi=0;npi<Tetra::NPI;npi++) { BOOST_TEST(my_alpha(npi) == a_eff); }

    uHeff = distrib(gen) + M;
    a_eff = func(uHeff); // case 1
    _uHeff.setConstant(uHeff);
    my_alpha = prm_t.get_alpha_eff<Tetra::NPI>(alpha_g, _uHeff );
    for(int npi=0;npi<Tetra::NPI;npi++) { BOOST_TEST(my_alpha(npi) == a_eff); }

    uHeff *= -1.0;
    a_eff = func(uHeff); // case 3
    _uHeff.setConstant(uHeff);
    my_alpha = prm_t.get_alpha_eff<Tetra::NPI>(alpha_g, _uHeff );
    for(int npi=0;npi<Tetra::NPI;npi++) { BOOST_TEST(my_alpha(npi) == a_eff); }
    }

BOOST_AUTO_TEST_SUITE_END()
