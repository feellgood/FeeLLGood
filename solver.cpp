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
    
    K.setFromTriplets(w_K_TH.begin(),w_K_TH.end());
    _solver.analyzePattern(K);// numerical values in K are not used
    _solver.factorize(K);

    if(_solver.info() != Eigen::Success)
        { std::cout <<"sparse matrix factorize(): decomposition failed.\n";exit(1); }
    else if (settings.verbose)
        { std::cout << "K factorized.\n"; }
    L_TH.setZero(2*NOD);

    std::for_each(refMsh->tet.begin(), refMsh->tet.end(),
                      [this](Tetra::Tet &my_elem) { my_elem.assemblage_vect(NOD,L_TH); } );
    std::for_each(refMsh->fac.begin(), refMsh->fac.end(),
                      [this](Facette::Fac &my_elem) { my_elem.assemblage_vect(NOD,L_TH); } );

    if (settings.verbose)
        { std::cout << "matrix assembly done in " << counter.millis() << std::endl; }

    refMsh->buildInitGuess(X_guess);// gamma0 division handled by function buildInitGuess

    counter.reset();

    Eigen::VectorXd sol = _solver.solveWithGuess(L_TH,X_guess);
    int nb_iter = _solver.iterations();
    double solver_error= _solver.error();

    if( (nb_iter > settings.MAXITER) || (solver_error > settings.TOL) )
        {
        if (settings.verbose)
            {
            std::cout << "solver: bicgstab FAILED after " << nb_iter
            << " iterations, in " << counter.millis() << std::endl;
            }
        return 1;
        }

        if (settings.verbose)
            {
            std::cout << "solver: bicgstab converged in " << nb_iter
            << " iterations, " << counter.millis() << std::endl;
            }

    v_max = refMsh->updateNodes(sol, t_prm.get_dt());//gamma0 multiplication handled by updateNodes
    return 0;
    }
