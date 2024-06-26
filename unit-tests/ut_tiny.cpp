#define BOOST_TEST_MODULE tinyTest

#include <boost/test/unit_test.hpp>

#include "tiny.h"
#include "ut_config.h"

BOOST_AUTO_TEST_SUITE(ut_tiny)


BOOST_AUTO_TEST_CASE(tiny_add, *boost::unit_test::tolerance(UT_TOL))
    {
    const int M = 3;
    const int N = 2;
    double A[M][N] = {{1, 2}, {3, 4}, {5, 6}};
    double B[M][N] = {{4, 31}, {-5, 0}, {0, 1}};
    double C[M][N] = {{0}};
    const double val_ref[M][N] = {{5, 33}, {-2, 4}, {5, 7}};
    tiny::add<double, M, N>(A, B, C);

    std::cout << "basic test addition" << std::endl;
    
    for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++)
            BOOST_CHECK(C[i][j] == val_ref[i][j]);
    }

BOOST_AUTO_TEST_CASE(tiny_mult_mat_vect)
    {
    const int M = 3;
    const int N = 2;
    double A[M][N] = {{1, 2}, {3, 4}, {5, 6}};
    double B[N] = {-2, 1};
    double C[M] = {0};
    const double val_ref[M] = {0, -2, -4};
    tiny::mult<double, M, N>(A, B, C);

    std::cout << "basic test mat*vect multiplication" << std::endl;
    for (int j = 0; j < M; j++)
        BOOST_CHECK(C[j] == val_ref[j]);
    }

BOOST_AUTO_TEST_CASE(tiny_mult_mat_mat)
    {
    const int M = 3;
    const int N = 2;
    const int P = 4;
    double A[M][N] = {{1, 2}, {3, 4}, {5, 6}};
    double B[N][P] = {{-2, 1, -3, 0}, {5, 3, 8, 4}};
    double C[M][P] = {{0}};
    const double val_ref[M][P] = {{8, 7, 13, 8}, {14, 15, 23, 16}, {20, 23, 33, 24}};
    tiny::mult<double, M, N, P>(A, B, C);

    std::cout << "basic test mat*mat multiplication" << std::endl;
    for (int i = 0; i < M; i++)
        for (int j = 0; j < P; j++)
            BOOST_CHECK(C[i][j] == val_ref[i][j]);
    }

BOOST_AUTO_TEST_CASE(tiny_sub)
    {
    const int M = 3;
    const int N = 2;
    double A[M][N] = {{1, 2}, {3, 4}, {5, 6}};
    double B[M][N] = {{4, 31}, {-5, 0}, {0, 1}};
    double C[M][N] = {{0}};
    const double val_ref[M][N] = {{-3, -29}, {8, 4}, {5, 5}};
    tiny::sub<double, M, N>(A, B, C);

    std::cout << "basic test substraction" << std::endl;
    for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++)
            BOOST_CHECK(C[i][j] == val_ref[i][j]);
    }

BOOST_AUTO_TEST_CASE(tiny_frobenius, *boost::unit_test::tolerance(UT_TOL))
    {
    const int M = 3;
    const int N = 2;
    double A[M][N] = {{1, 2}, {3, 4}, {5, 6}};

    double result = tiny::frob_norm<double, M, N>(A);

    std::cout << "test frobenius norm" << std::endl;
    BOOST_TEST(result == sqrt(91.0));
    }

BOOST_AUTO_TEST_CASE(tiny_dist, *boost::unit_test::tolerance(UT_TOL))
    {
    const int M = 3;
    const int N = 2;
    double A[M][N] = {{1, 2}, {3, 4}, {5, 6}};
    double B[M][N] = {{4, 31}, {-5, 0}, {1e-17, 1}};
    double C[M][N] = {{0}};
    const double val_ref[M][N] = {{-3, -29}, {8, 4}, {5, 5}};
    tiny::sub<double, M, N>(A, B, C);

    double result = tiny::dist<double, M, N>(C, val_ref);

    std::cout << "test dist = frobenius_norm(A-B)" << std::endl;
    BOOST_TEST(result == 0.0);
    }

BOOST_AUTO_TEST_SUITE_END()
