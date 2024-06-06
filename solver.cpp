#include "chronometer.h"
#include "linear_algebra.h"

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

int LinAlgebra::solver(timing const &t_prm)
    {
    chronometer counter(2);
    std::vector<Eigen::Triplet<double>> w_K_TH;

    std::for_each(refMsh->tet.begin(), refMsh->tet.end(),
                      [this,&w_K_TH](Tetra::Tet &my_elem) { my_elem.assemblage_mat(NOD,w_K_TH); } );
    
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double,Eigen::RowMajor>,Eigen::IncompleteLUT<double>> _solver;
    _solver.setTolerance(TOL);
    _solver.setMaxIterations(MAXITER);

    Eigen::SparseMatrix<double,Eigen::RowMajor> K(2*NOD,2*NOD);
    K.setFromTriplets(w_K_TH.begin(),w_K_TH.end());
    if (verbose)
        {
        std::cout << "matrix assembly done in " << counter.millis() << std::endl;
        counter.reset();
        }
    _solver.preconditioner().setDroptol(ILU_tol);
    _solver.preconditioner().setFillfactor(ILU_fill_factor);
    _solver.analyzePattern(K);// numerical values in K are not used
    _solver.factorize(K);

    if(verbose)
        {
        std::cout << "ILU preconditionner (tolerance;filling factor) = ("<< ILU_tol <<";"<< ILU_fill_factor << ")\n";
        }

    if(_solver.info() != Eigen::Success)
        {
        std::cout <<"sparse matrix decomposition failed" << std::endl;
        exit(1);
        }
    else if (verbose)
        { std::cout << "sparse matrix factorization done in " << counter.millis() << std::endl; }

    Eigen::VectorXd L_TH(2*NOD);// RHS vector of the system to solve
    L_TH.setZero(2*NOD);
    std::for_each(refMsh->tet.begin(), refMsh->tet.end(),
                      [this,&L_TH](Tetra::Tet &my_elem) { my_elem.assemblage_vect(NOD,L_TH); } );
    std::for_each(refMsh->fac.begin(), refMsh->fac.end(),
                      [this,&L_TH](Facette::Fac &my_elem) { my_elem.assemblage_vect(NOD,L_TH); } );

    Eigen::VectorXd X_guess(2*NOD);
    refMsh->buildInitGuess(X_guess);// gamma0 division handled by function buildInitGuess

    counter.reset();

    Eigen::VectorXd sol = _solver.solveWithGuess(L_TH,X_guess);
    int nb_iter = _solver.iterations();
    double solver_error= _solver.error();

    if( (nb_iter > MAXITER) || (solver_error > TOL) )
        {
        if (verbose)
            {
            std::cout << "solver: bicgstab FAILED after " << nb_iter
            << " iterations, in " << counter.millis() << std::endl;
            }
        return 1;
        }

        if (verbose)
            {
            std::cout << "solver: bicgstab converged in " << nb_iter
            << " iterations, " << counter.millis() << std::endl;
            }

    v_max = refMsh->updateNodes(sol, t_prm.get_dt());//gamma0 multiplication handled by updateNodes
    return 0;
    }
