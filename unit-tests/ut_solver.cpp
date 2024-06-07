#define BOOST_TEST_MODULE solverTest

#include <boost/test/unit_test.hpp>

#include <cstdlib>
#include <iostream>
#include <random>

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

#include "ut_config.h"  // for tolerance UT_TOL macro
#include "ut_solver.h"

BOOST_AUTO_TEST_SUITE(ut_solver)

/*
check a nested eigen expression to compute v_max in linAlgebra::solver() method
*/
BOOST_AUTO_TEST_CASE(vectorManip, *boost::unit_test::tolerance(UT_TOL))
    {
    const int N = 4;
    Eigen::VectorXd x(2*N);
    x << 1,2,3,-4,5,6,7,8;
    double x_max = (x.reshaped<Eigen::RowMajor>(2,N).colwise().norm()).maxCoeff();
    double result = sqrt(4.0*4.0 + 8.0*8.0);
    std::cout << "x_max= " << x_max << " should be: " << result << std::endl;
    BOOST_CHECK( x_max == result );
    }
/*
silly_problem_solver tests if eigen::bicgstab solves Ax=b with A = Id and b=(1  1 ... 1) with default diagonal precond in 1 iteration
*/
BOOST_AUTO_TEST_CASE(silly_problem_solver, *boost::unit_test::tolerance(UT_TOL))
    {
    const int MAX_ITER = 100;
    const double _TOL = 1e-8;
    int N = 10000;

    Eigen::VectorXd x(N), b(N);
    Eigen::SparseMatrix<double,Eigen::RowMajor> A(N,N);
    std::vector<Eigen::Triplet<double>> coeffs;
    
    for(int i=0;i<N;i++)
        {
        coeffs.push_back(Eigen::Triplet<double>(i,i,1.0) );
        b(i) = 1.0;
        }
    A.setFromTriplets(coeffs.begin(),coeffs.end());
    
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
    solver.setMaxIterations(MAX_ITER);
    solver.setTolerance(_TOL);
    solver.compute(A);
    x = solver.solve(b);
    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error()      << std::endl;
    BOOST_CHECK( solver.iterations() == 1);
    BOOST_CHECK( (x-b).norm() == 0.0 );
    }

/*
rand_sp_mat_problem_solver tests if eigen::bicgstab solves Ax=b with
 A = Id + extras 1 coefficients randomly placed, with their symmetric 
 and b=(1  1 ... 1) and an initial stupid guess x=2
*/
BOOST_AUTO_TEST_CASE(rand_sp_mat_problem_solver, *boost::unit_test::tolerance(UT_TOL))
    {
    const int MAX_ITER = 2000;
    const double _TOL = 1e-8;
    int N = 10000;

    Eigen::VectorXd x(N), b(N);
    Eigen::SparseMatrix<double,Eigen::RowMajor> A(N,N);
    std::vector<Eigen::Triplet<double>> coeffs;

    for(int i=0;i<N;i++)
        {
        coeffs.push_back(Eigen::Triplet<double>(i,i,1.0) );
        b(i) = 1.0;
        }

    std::mt19937 gen(my_seed());
    std::uniform_int_distribution<> distrib(0, N-1);

    for (int nb=0; nb<400; nb++)
        {
        int i = distrib(gen);
        int j = distrib(gen);
        coeffs.push_back(Eigen::Triplet<double>(i,j,1.0) );
        coeffs.push_back(Eigen::Triplet<double>(j,i,1.0) );
        }
    A.setFromTriplets(coeffs.begin(),coeffs.end());
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
    solver.setMaxIterations(MAX_ITER);
    solver.setTolerance(_TOL);
    solver.compute(A);
    x = solver.solveWithGuess(b, Eigen::VectorXd::Constant(N,2.0));
    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error()      << std::endl;
    BOOST_CHECK( solver.iterations() > 1);
    BOOST_CHECK( (x-b).norm() != 0.0 );
    double err_result = ((A*x)-b).norm(); 
    std::cout << "norm(Ax - b)= " << err_result << std::endl;
    BOOST_CHECK( err_result < _TOL );
    }

/*
rand_asym_sp_mat_problem_solver tests if eigen::bicgstab solves Ax=b with
 A = Id + extras 1 coefficients randomly placed ( A is not symmetric) 
 and b=(1  1 ... 1) and an initial stupid guess x=2
*/
BOOST_AUTO_TEST_CASE(rand_asym_sp_mat_problem_solver, *boost::unit_test::tolerance(UT_TOL))
    {
    const int MAX_ITER = 2000;
    const double _TOL = 1e-8;
    const int N = 10000;

    Eigen::VectorXd x(N), b(N);
    Eigen::SparseMatrix<double,Eigen::RowMajor> A(N,N);
    std::vector<Eigen::Triplet<double>> coeffs;

    build_coeffs<N,800>(coeffs,b);
    A.setFromTriplets(coeffs.begin(),coeffs.end());
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
    solver.setMaxIterations(MAX_ITER);
    solver.setTolerance(_TOL);
    solver.compute(A);
    x = solver.solveWithGuess(b, Eigen::VectorXd::Constant(N,2.0));
    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error()      << std::endl;
    BOOST_CHECK( solver.iterations() > 1);
    BOOST_CHECK( (x-b).norm() != 0.0 );
    BOOST_CHECK( ((A*x)-b).norm() <= _TOL );
    }

BOOST_AUTO_TEST_SUITE_END()

