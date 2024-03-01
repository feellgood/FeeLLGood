#define BOOST_TEST_MODULE OMPsolverTest

#include <boost/test/unit_test.hpp>

#include <cstdlib>
#include <iostream>
#include <random>

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

#include "ut_config.h"  // for tolerance UT_TOL macro
#include "ut_solver.h"

BOOST_AUTO_TEST_SUITE(ut_OMP_solver)

BOOST_AUTO_TEST_CASE(class_embedded_solver, *boost::unit_test::tolerance(UT_TOL))
    {
    const int MAX_ITER = 1000;
    const double _TOL = 1e-7;
    const int NOD = 20000;

    omp_set_num_threads(4);
    Eigen::setNbThreads(4);

    DummyLinAlgebra<NOD> bob(MAX_ITER,_TOL);
    
    double t(0);
    
    bob.solve(t);
    double t_end(0.01);
    BOOST_CHECK( t == t_end);
    }

BOOST_AUTO_TEST_SUITE_END()
