#define BOOST_TEST_MODULE algebra_bicg_test

#include <boost/test/unit_test.hpp>

#include <cstdlib>
#include <iostream>
#include <random>

#include "algebra/algebra.h"
#include "algebra/bicg.h"
#include "ut_config.h"  // for tolerance UT_TOL macro
#include "sparse_matrix.h"

BOOST_AUTO_TEST_SUITE(ut_algebra_bicg)

/*
silly_problem_solver tests if algebra::bicg solves Ax=b with A = Id and b=(1  1 ... 1) with default diagonal precond in 0 iteration
*/
BOOST_AUTO_TEST_CASE(silly_problem_solver, *boost::unit_test::tolerance(UT_TOL))
    {
    const int MAX_ITER = 100;
    const double _TOL = 1e-8;
    int N = 10000;
    algebra::iteration iter("bicg",_TOL,false,MAX_ITER);

    std::vector<double> x(N,0.0), b(N,1.0);
    std::vector<MatrixCoefficient> coefficients;
    for(int i=0;i<N;i++) { coefficients.push_back({i, i, 1.0}); }
    algebra::SparseMatrix Ar = buildSparseMat(N, coefficients);
    algebra::bicg<double>(iter,Ar,x,b);
    std::cout << "#iterations:     " << iter.get_iteration() << std::endl;
    std::cout << "estimated error: " << iter.get_res()      << std::endl;
    BOOST_CHECK( iter.get_iteration() == 0);
    std::vector<double> y(N);
    algebra::mult(Ar,x,y);// y = Ar * x
    algebra::sub(b,y); // y -= b;
    BOOST_CHECK( algebra::norm<double>(y) == 0.0 );
    }

/*
rand_sp_mat_problem_solver tests if algebra::bicg solves Ax=b with
 A = Id + extras 1 coefficients randomly placed, with their symmetric 
 and b=(1  1 ... 1) and an initial stupid guess x=2
*/

BOOST_AUTO_TEST_CASE(rand_sp_mat_problem_solver, *boost::unit_test::tolerance(UT_TOL))
    {
    const int MAX_ITER = 2000;
    const double _TOL = 1e-8;
    int N = 10000;
    algebra::iteration iter("bicg",_TOL,false,MAX_ITER);
    std::vector<double> x(N,2.0), b(N,1.0);
    std::vector<MatrixCoefficient> coefficients;
    for(int i=0;i<N;i++) { coefficients.push_back({i, i, 1.0}); }

    std::mt19937 gen(my_seed());
    std::uniform_int_distribution<> distrib(0, N-1);

    for (int nb=0; nb<400; nb++)
        {
        int i = distrib(gen);
        int j = distrib(gen);
        coefficients.push_back({i, j, 1.0});
        coefficients.push_back({j, i, 1.0});
        }

    algebra::SparseMatrix Ar = buildSparseMat(N, coefficients);
    algebra::bicg<double>(iter,Ar,x,b);
    std::cout << "#iterations:     " << iter.get_iteration() << std::endl;
    std::cout << "estimated error: " << iter.get_res()      << std::endl;
    BOOST_CHECK( iter.get_iteration() > 1);
    std::vector<double> y(N);
    algebra::mult(Ar,x,y);// y = Ar * x
    algebra::sub(b,y); // y -= b;
    double err_result = algebra::norm<double>(y);
    std::cout << "norm(Ax - b)= " << err_result << std::endl;
    BOOST_CHECK( err_result < _TOL );
    }

/*
rand_asym_sp_mat_problem_solver tests if algebra::bicg solves Ax=b with
 A = Id + extras 1 coefficients randomly placed ( A is not symmetric) 
 and b=(1  1 ... 1) and an initial stupid guess x=2
*/

BOOST_AUTO_TEST_CASE(rand_asym_sp_mat_problem_solver, *boost::unit_test::tolerance(UT_TOL))
    {
    const int MAX_ITER = 2000;
    const double _TOL = 1e-8;
    const int N = 10000;
    algebra::iteration iter("bicg",_TOL,false,MAX_ITER);
    std::vector<double> x(N,2.0), b(N,1.0);
    std::vector<MatrixCoefficient> coefficients;
    for(int i=0;i<N;i++) { coefficients.push_back({i, i, 1.0}); }
    std::mt19937 gen(my_seed());
    std::uniform_int_distribution<> distrib(0, N-1);

    for (int nb=0; nb<800; nb++)
        {
        int i = distrib(gen);
        int j = distrib(gen);
        coefficients.push_back({i, j, 1.0});
        }

    algebra::SparseMatrix Ar = buildSparseMat(N, coefficients);
    algebra::bicg<double>(iter,Ar,x,b);

    std::cout << "#iterations:     " << iter.get_iteration() << std::endl;
    std::cout << "estimated error: " << iter.get_res()      << std::endl;
    BOOST_CHECK( iter.get_iteration() > 1);
    std::vector<double> y(N);
    algebra::mult(Ar,x,y);// y = Ar * x
    algebra::sub(b,y); // y -= b;
    double err_result = algebra::norm<double>(y);
    std::cout << "norm(Ax - b)= " << err_result << std::endl;
    BOOST_CHECK( err_result < _TOL );
    }

BOOST_AUTO_TEST_SUITE_END()

