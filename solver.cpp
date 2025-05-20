#include "chronometer.h"
#include "linear_algebra.h"

#include "algebra/algebra.h"
#include "algebra/bicg.h"

#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>

int LinAlgebra::solver(timing const &t_prm)
    {
    chronometer counter(2);

    algebra::w_sparseMat Kw(2*NOD);
    std::for_each(refMsh->tet.begin(), refMsh->tet.end(),
                      [this,&Kw](Tetra::Tet &my_elem) { my_elem.assemblage_mat(NOD,Kw); } );

    algebra::r_sparseMat Kr(Kw);
    if (verbose)
        {
        std::cout << "matrix assembly done in " << counter.millis() << std::endl;
        counter.reset();
        }

    algebra::iteration iter(TOL);
    iter.set_maxiter(MAXITER);
    iter.set_noisy(verbose);
    
    std::fill(L_rhs.begin(),L_rhs.end(),0);
    std::for_each(refMsh->tet.begin(), refMsh->tet.end(),
                      [this](Tetra::Tet &my_elem) { my_elem.assemblage_vect(NOD,L_rhs); } );
    std::for_each(refMsh->fac.begin(), refMsh->fac.end(),
                      [this](Facette::Fac &my_elem) { my_elem.assemblage_vect(NOD,L_rhs); } );

    Eigen::SparseMatrix<double,Eigen::RowMajor> K(2*NOD,2*NOD);
    for(int i=0;i<(2*NOD);i++)
        {
        for(int j=0;j<(2*NOD);j++)
            {
            double val = Kr(i,j);
            if (val!=0) { K.insert(i,j) = val; }
            }
        }

    if (verbose)
        {
        std::cout << "matrix assembly done in " << counter.millis() << std::endl;
        counter.reset();
        }
//ref code
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double,Eigen::RowMajor>,Eigen::IncompleteLUT<double>> _solver;
    _solver.setTolerance(TOL);
    _solver.setMaxIterations(MAXITER);
    _solver.preconditioner().setDroptol(ILU_tol);
    _solver.preconditioner().setFillfactor(ILU_fill_factor);
    _solver.analyzePattern(K);// numerical values in K are not used
    _solver.factorize(K);

    if(verbose)
        { std::cout << "ILU preconditionner (tolerance;filling factor) = ("<< ILU_tol <<";"<< ILU_fill_factor << ")\n"; }

    if(_solver.info() != Eigen::Success)
        {
        std::cout <<"sparse matrix decomposition failed" << std::endl;
        exit(1);
        }
    else if (verbose)
        { std::cout << "sparse matrix factorization done in " << counter.millis() << std::endl; }

    buildInitGuess(Xw);// gamma0 division handled by function buildInitGuess    

    Eigen::VectorXd L_TH(2*NOD);// RHS vector of the system to solve
    Eigen::VectorXd X_guess(2*NOD);
    for(int i=0;i<(2*NOD);i++)
        {
        L_TH(i) = L_rhs[i];
        X_guess(i) = Xw[i];
        }
    counter.reset();

    Eigen::VectorXd sol = _solver.solveWithGuess(L_TH,X_guess);
// end ref code
    for(int i=0;i<(2*NOD);i++) { Xw[i] = sol(i); }

//    double residu = algebra::bicg<double>(iter, Kr, Xw, L_rhs);
//    if(verbose) std::cout << "residu= " << residu << std::endl;

    int nb_iter = _solver.iterations();//iter.get_iteration();
    double _solver_error= _solver.error();//iter.get_res();

    if( (nb_iter > MAXITER) || (_solver_error > TOL) )
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
        << " iterations, " << counter.millis() << ", residu= "<< _solver_error << std::endl;
        }
    
    double v2max(0.0);
    const double dt = t_prm.get_dt();
    for (int i = 0; i < NOD; i++)
        {
        double vp = Xw[i];
        double vq = Xw[NOD + i];
        double v2 = Nodes::sq(vp) + Nodes::sq(vq);
        if (v2 > v2max)
            { v2max = v2; }
        
        refMsh->updateNode(i, vp, vq, dt);//gamma0 multiplication handled in updateNode
        }
    v_max = sqrt(v2max);
    v_max *= gamma0;
    
    return 0;
    }
