#include "chronometer.h"
#include "linear_algebra.h"

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

int LinAlgebra::solver(timing const &t_prm, long nt)
    {
    chronometer counter(2);

    write_matrix K_TH(2 * NOD, 2 * NOD);
    //std::vector<Eigen::Triplet<double>> w_K;
   
    std::for_each(refMsh->tet.begin(), refMsh->tet.end(),
                      [this,&K_TH](Tetra::Tet &my_elem) { my_elem.assemblage_mat(NOD,K_TH); } );
    std::for_each(refMsh->fac.begin(), refMsh->fac.end(),
                      [this,&K_TH](Facette::Fac &my_elem) { my_elem.assemblage_mat(NOD,K_TH); } );
    
    std::vector<double> w_L_TH(2 * NOD, 0);
    //Eigen::VectorXd L_TH = Eigen::VectorXd::Zero(2*NOD);
    
    std::for_each(refMsh->tet.begin(), refMsh->tet.end(),
                      [this,&w_L_TH](Tetra::Tet &my_elem) { my_elem.assemblage_vect(NOD,w_L_TH); } );
    std::for_each(refMsh->fac.begin(), refMsh->fac.end(),
                      [this,&w_L_TH](Facette::Fac &my_elem) { my_elem.assemblage_vect(NOD,w_L_TH); } );
    
    if (settings.verbose)
        { std::cout << "matrix assembly done in " << counter.millis() << std::endl; }

    read_matrix Kr(2 * NOD, 2 * NOD);
    gmm::copy(K_TH, Kr);
    std::vector<double> X(2 * NOD);  // filled with zero by default
    refMsh->buildInitGuess(X);// gamma0 division handled by function buildInitGuess

    gmm::iteration bicg_iter(1e-6);
    bicg_iter.set_maxiter(settings.MAXITER);

/*
    Eigen::SparseMatrix<double> K(2*NOD,2*NOD);//,Eigen::RowMajor
    K.reserve(2*NOD); // initial mem alloc to the size of the diagonal
    K.setFromTriplets(w_K.begin(),w_K.end());
*/
    // Let GMM be noisy only at verbosity levels higher than one.
    bicg_iter.set_noisy(settings.verbose > 1);


    if (!prc || (nt - prc_time_step >= settings.REFRESH_PRC))
        {
        if (prc)
            {
            delete prc;
            }
        prc = new gmm::diagonal_precond<read_matrix>(Kr);
        prc_time_step = nt;
        }

//    const double _TOL(1e-6);

    counter.reset();

gmm::bicgstab(Kr, X, w_L_TH, *prc, bicg_iter);
    /*
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> _solver;
    _solver.compute(K);
    if(_solver.info() != Eigen::Success)
        { std::cout <<"sparse matrix decomposition failed.\n"; }

    _solver.setMaxIterations(10000);//settings.MAXITER;
    _solver.setTolerance(_TOL);
    double solver_error(0);

    Eigen::VectorXd L_TH(2*NOD);
    for(int i=0;i< (2*NOD);i++) {L_TH(i) = w_L_TH[i];}
    
    X = _solver.solve(L_TH);//.solveWithGuess(L_TH,L_TH);//second parameter of solveWithGuess is guess
    int nb_iter = _solver.iterations();
    solver_error = _solver.error();

    bool converged = (nb_iter < 10000);
    if (!std::isnan(solver_error))
        { converged = converged && (solver_error < _TOL); }
    else
        {std::cout <<"solver returned NaN.";exit(1);}
*/

    if (!(bicg_iter.converged())) //(!converged)
        {
        if (settings.verbose)
            {
            std::cout << "solver: bicg (prc " << (nt % (settings.REFRESH_PRC)) << ") FAILED after "
            //std::cout << "solver: bicgstab FAILED after " << nb_iter
            << bicg_iter.get_iteration() << " iterations, in " << counter.millis() << std::endl;
            
            //std::cout << "error: " << solver_error << std::endl;
            }
        return 1;
        }
    else
        {
        if (settings.verbose)
            {
            std::cout << "solver: bicg (prc " << (nt % (settings.REFRESH_PRC)) << ") converged in "<< bicg_iter.get_iteration() 
            //std::cout << "solver: bicgstab converged in " << nb_iter
            << " iterations, " << counter.millis() << std::endl;
            }
        }

    // gamma0 multiplication handled within function updateNodes
    v_max = refMsh->updateNodes(X, t_prm.get_dt());
    return 0;
    }
