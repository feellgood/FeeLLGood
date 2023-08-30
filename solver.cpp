#include "chronometer.h"
#include "linear_algebra.h"

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

int LinAlgebra::solver(timing const &t_prm)
    {
    const double _TOL(1e-6);
    chronometer counter(2);
    
    std::vector<Eigen::Triplet<double>> w_K_TH;
   
    std::for_each(refMsh->tet.begin(), refMsh->tet.end(),
                      [this,&w_K_TH](Tetra::Tet &my_elem) { my_elem.assemblage_mat(NOD,w_K_TH); } );
    std::for_each(refMsh->fac.begin(), refMsh->fac.end(),
                      [this,&w_K_TH](Facette::Fac &my_elem) { my_elem.assemblage_mat(NOD,w_K_TH); } );
    
    Eigen::VectorXd L_TH = Eigen::VectorXd::Zero(2*NOD);
    
    std::for_each(refMsh->tet.begin(), refMsh->tet.end(),
                      [this,&L_TH](Tetra::Tet &my_elem) { my_elem.assemblage_vect(NOD,L_TH); } );
    std::for_each(refMsh->fac.begin(), refMsh->fac.end(),
                      [this,&L_TH](Facette::Fac &my_elem) { my_elem.assemblage_vect(NOD,L_TH); } );
    

    Eigen::SparseMatrix<double,Eigen::RowMajor> K(2*NOD,2*NOD);
    K.reserve(2*NOD); // initial mem alloc to the size of the diagonal
    K.setFromTriplets(w_K_TH.begin(),w_K_TH.end());

    if (settings.verbose)
        { std::cout << "matrix assembly done in " << counter.millis() << std::endl; }

    Eigen::VectorXd X_guess(2*NOD);
    refMsh->buildInitGuess(X_guess);// gamma0 division handled by function buildInitGuess

    counter.reset();

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double,Eigen::RowMajor>,Eigen::IncompleteLUT<double>> _solver; //,Eigen::DiagonalPreconditioner<double>
    _solver.setTolerance(_TOL);
    _solver.setMaxIterations(settings.MAXITER);
    
    _solver.compute(K);
    if(_solver.info() != Eigen::Success)
        { std::cout <<"sparse matrix decomposition failed.\n";exit(1); }

    Eigen::VectorXd sol = _solver.solveWithGuess(L_TH,X_guess);// solve(L_TH); second parameter of solveWithGuess is guess

    int nb_iter = _solver.iterations();
    double solver_error= _solver.error();
    bool converged = ((nb_iter < settings.MAXITER) && (solver_error < _TOL)); // might be the same as (_solver.info() == Eigen::Success)

    if(!converged)
        {
        if (settings.verbose)
            {
            std::cout << "solver: bicgstab FAILED after " << nb_iter
            << " iterations, in " << counter.millis() << std::endl;
            }
        return 1;
        }
    else
        {
        if (settings.verbose)
            {
            std::cout << "solver: bicgstab converged in " << nb_iter
            << " iterations, " << counter.millis() << std::endl;
            }
        }

    //gamma0 multiplication handled by function updateNodes
    v_max = refMsh->updateNodes(sol, t_prm.get_dt());
    return 0;
    }
