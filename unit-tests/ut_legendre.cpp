#define BOOST_TEST_MODULE legendreTest

#include <boost/test/unit_test.hpp>
#include <iomanip>
#include <random>
#include <cmath> // assoc_legendre()

#include "ut_config.h"
#include "../src/fmm_demag.h"
#include "../src/scalFMMlight/Kernels/Rotation/FRotationKernel.hpp"

/*
 *this unit test is checking the computation of associated Legendre polynom redefined in scalFMM. It
 does also check indexation in the one-D table of the P_{l,m} values, needed for Y_{l,m} Kernel
 * */

BOOST_AUTO_TEST_SUITE(ut_legendre)
BOOST_AUTO_TEST_CASE(legendre, *boost::unit_test::tolerance(10.0*UT_TOL))
    {
/**
 * std::assoc_legendre does compute the associated legendre polynom, but without the Condon-Shortley phase (-1)^m
 * see https://en.cppreference.com/w/cpp/numeric/special_functions/assoc_legendre for more infos.
 * */

    using FReal = scal_fmm::FReal;
    const int P = scal_fmm::P;
    const int SizeArray = ((P+2)*(P+1))/2;

    double values[SizeArray];
    double theta = 1.234;
    values[0] = 1.0;
    values[1] = cos(theta);
    values[2] = -sin(theta);
    double minusOnePow_m;

    auto at_lm = scal_fmm::KernelClass::atLm;
    scal_fmm::KernelClass::computeLegendre(values);
    for(int l=0;l<P;l++)
        {
        minusOnePow_m = 1.0;
        for(int m=0;m<=l;m++)
            {
            BOOST_TEST( values[at_lm(l,m)] == minusOnePow_m*std::assoc_legendre(l,m,cos(theta)) );
            minusOnePow_m *= -1.0;
            }
        }
    }
BOOST_AUTO_TEST_SUITE_END()
