#define BOOST_TEST_MODULE algebra_bicg_test

#include <boost/test/unit_test.hpp>

#include <cstdlib>
#include <iostream>
#include <random>

#include "algebra/algebra.h"
#include "algebra/bicg.h"
#include "algebra/block_precond.h"
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

    /* Id + symmetric ones may be singular (an isolated pair i, j gives the block ((1 1)(1 1))) and
     then bicg may diverge: the number of ones of each row is added to its diagonal, the matrix is
     strictly diagonally dominant, hence invertible */
    std::vector<int> nb_ones(N, 0);
    for (int nb=0; nb<400; nb++)
        {
        int i = distrib(gen);
        int j = distrib(gen);
        coefficients.push_back({i, j, 1.0});
        coefficients.push_back({j, i, 1.0});
        nb_ones[i]++;
        nb_ones[j]++;
        }
    for(int i=0;i<N;i++) { coefficients.push_back({i, i, double(nb_ones[i])}); }

    algebra::SparseMatrix Ar = buildSparseMat(N, coefficients);
    algebra::bicg<double>(iter,Ar,x,b);
    std::cout << "#iterations:     " << iter.get_iteration() << std::endl;
    std::cout << "estimated error: " << iter.get_res()      << std::endl;
    BOOST_CHECK( iter.get_iteration() > 1);
    std::vector<double> y(N);
    algebra::mult(Ar,x,y);// y = Ar * x
    algebra::sub(b,y); // y -= b;
    // bicg stops when norm(Ax - b) <= _TOL * norm(b): relative criterion
    double err_result = algebra::norm<double>(y) / algebra::norm<double>(b);
    std::cout << "norm(Ax - b)/norm(b)= " << err_result << std::endl;
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
    // bicg stops when norm(Ax - b) <= _TOL * norm(b): relative criterion
    double err_result = algebra::norm<double>(y) / algebra::norm<double>(b);
    std::cout << "norm(Ax - b)/norm(b)= " << err_result << std::endl;
    BOOST_CHECK( err_result < _TOL );
    }


/* coefficients of a system with 3 unknowns per node (numbered node by node, 3*i + d), mimicking the
 spin diffusion problem: diagonal blocks d Id + c [e]x (relaxation and strong antisymmetric
 precession coupling, random unit vector e per node), and a symmetric coupling between random pairs
 of nodes (diffusion), the matrix being strictly block diagonally dominant */
static void spin_like_coefficients(const int nbNod, std::mt19937 &gen,
                                   std::vector<MatrixCoefficient> &coefficients)
    {
    std::uniform_real_distribution<> distrib(-1.0, 1.0);
    std::uniform_int_distribution<> node(0, nbNod - 1);
    const double c = 20.0;  // precession >> relaxation: the diagonal preconditioner is poor
    std::vector<double> diag(nbNod, 1.0);
    for (int nb = 0; nb < 4*nbNod; nb++)
        {
        int i = node(gen), j = node(gen);
        if (i == j) continue;
        for (int d = 0; d < 3; d++)
            {
            coefficients.push_back({3*i + d, 3*j + d, -1.0});
            coefficients.push_back({3*j + d, 3*i + d, -1.0});
            }
        diag[i] += 1.0;
        diag[j] += 1.0;
        }
    for (int i = 0; i < nbNod; i++)
        {
        double e[3] = {distrib(gen), distrib(gen), distrib(gen)};
        const double ne = sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
        for (int d = 0; d < 3; d++) e[d] /= ne;
        // row x: (s x e)_x = s_y e_z - s_z e_y, etc.
        coefficients.push_back({3*i, 3*i, diag[i]});
        coefficients.push_back({3*i, 3*i + 1, c*e[2]});
        coefficients.push_back({3*i, 3*i + 2, -c*e[1]});
        coefficients.push_back({3*i + 1, 3*i, -c*e[2]});
        coefficients.push_back({3*i + 1, 3*i + 1, diag[i]});
        coefficients.push_back({3*i + 1, 3*i + 2, c*e[0]});
        coefficients.push_back({3*i + 2, 3*i, c*e[1]});
        coefficients.push_back({3*i + 2, 3*i + 1, -c*e[0]});
        coefficients.push_back({3*i + 2, 3*i + 2, diag[i]});
        }
    }

/*
block_diag_precond_inverse: on a block diagonal matrix, the block diagonal preconditioner is the
exact inverse (M A x = x), and the identity on the Dirichlet components
*/
BOOST_AUTO_TEST_CASE(block_diag_precond_inverse)
    {
    const int nbNod = 100;
    std::mt19937 gen(my_seed());
    std::uniform_real_distribution<> distrib(-1.0, 1.0);
    std::vector<MatrixCoefficient> coefficients;
    for (int i = 0; i < nbNod; i++)
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                { coefficients.push_back({3*i + r, 3*i + c, (r == c ? 5.0 : 0.0) + distrib(gen)}); }
    algebra::SparseMatrix A = buildSparseMat(3*nbNod, coefficients);

    std::vector<int> ld;  // no Dirichlet condition
    algebra::BlockDiagPrecond<3> M(A, nbNod, ld);
    std::vector<double> x(3*nbNod), Ax(3*nbNod), MAx(3*nbNod);
    for (double &xi : x) xi = distrib(gen);
    algebra::mult(A, x, Ax);
    M.apply(Ax, MAx);
    algebra::sub(x, MAx);  // MAx -= x
    const double err = algebra::norm(MAx)/algebra::norm(x);
    std::cout << "norm(M A x - x)/norm(x)= " << err << std::endl;
    BOOST_CHECK( err < 1e-13 );

    // Dirichlet component 3*7+1: row and column replaced by the identity before inversion
    ld = {3*7 + 1};
    algebra::BlockDiagPrecond<3> Md(A, nbNod, ld);
    std::vector<double> e(3*nbNod, 0.0), Me(3*nbNod);
    e[3*7 + 1] = 1.0;
    Md.apply(e, Me);
    // (cofactors divided by the determinant: equal to 1 and 0 up to rounding errors)
    BOOST_CHECK( std::abs(Me[3*7]) < 1e-14 );
    BOOST_CHECK( std::abs(Me[3*7 + 1] - 1.0) < 1e-14 );
    BOOST_CHECK( std::abs(Me[3*7 + 2]) < 1e-14 );
    }

/*
spin_like_block_precond_solver: bicg_dir_prec with the block diagonal preconditioner solves a system
with a strong antisymmetric coupling of the 3 components of each node (as the spin diffusion
problem) and Dirichlet conditions, in fewer iterations than with the diagonal preconditioner
*/
BOOST_AUTO_TEST_CASE(spin_like_block_precond_solver)
    {
    const int MAX_ITER = 5000;
    const double _TOL = 1e-8;
    const int nbNod = 2000;
    const int N = 3*nbNod;
    std::mt19937 gen(my_seed());
    std::vector<MatrixCoefficient> coefficients;
    spin_like_coefficients(nbNod, gen, coefficients);
    algebra::SparseMatrix A = buildSparseMat(N, coefficients);

    std::uniform_real_distribution<> distrib(-1.0, 1.0);
    std::vector<double> b(N), xd(N, 0.0);
    for (double &bi : b) bi = distrib(gen);
    std::vector<int> ld;  // Dirichlet: the 3 components of the first 10 nodes
    for (int i = 0; i < 10; i++)
        for (int d = 0; d < 3; d++)
            {
            ld.push_back(3*i + d);
            xd[3*i + d] = distrib(gen);
            }

    // diagonal preconditioner
    algebra::iteration<double> iter_diag("bicg_dir",_TOL,false,MAX_ITER);
    std::vector<double> x_diag(N, 0.0);
    algebra::bicg_dir(iter_diag, A, x_diag, b, xd, ld);

    // block diagonal preconditioner
    algebra::iteration<double> iter_block("bicg_dir_prec",_TOL,false,MAX_ITER);
    std::vector<double> x(N, 0.0);
    algebra::BlockDiagPrecond<3> M(A, nbNod, ld);
    algebra::bicg_dir_prec(iter_block, A, x, b, xd, ld,
                           [&M](const std::vector<double> &p, std::vector<double> &phat)
                           { M.apply(p, phat); });
    std::cout << "#iterations: diagonal " << iter_diag.get_iteration()
              << ", block diagonal " << iter_block.get_iteration() << std::endl;
    BOOST_CHECK( iter_block.status == algebra::CONVERGED );
    BOOST_CHECK( iter_block.get_iteration() < iter_diag.get_iteration() );

    // Dirichlet values imposed, residual of the free components
    for (const int i : ld)
        { BOOST_CHECK( x[i] == xd[i] ); }
    std::vector<double> r(N), bfree(b);
    algebra::mult(A, x, r);
    algebra::sub(b, r);     // r = Ax - b
    algebra::applyMask(ld, r);
    std::vector<double> Axd(N);
    algebra::mult(A, xd, Axd);
    algebra::sub(Axd, bfree);  // bfree = b - A xd
    algebra::applyMask(ld, bfree);
    const double err = algebra::norm(r)/algebra::norm(bfree);
    std::cout << "norm(Ax - b)/norm(b - A xd) on the free components= " << err << std::endl;
    BOOST_CHECK( err < _TOL );
    }

BOOST_AUTO_TEST_SUITE_END()

