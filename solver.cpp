#include "chronometer.h"
#include "linear_algebra.h"

#include "algebra/algebra.h"
#include "algebra/bicg.h"

int LinAlgebra::solver(timing const &t_prm)
    {
    chronometer counter(2);

    K.clear();
    std::for_each(std::execution::par, refMsh->tet.begin(), refMsh->tet.end(),
                      [this](Tetra::Tet &my_elem) { my_elem.assemblage_mat(K); } );

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
                      [this](Tetra::Tet &my_elem) { my_elem.assemblage_vect(L_rhs); } );
    std::for_each(refMsh->fac.begin(), refMsh->fac.end(),
                      [this](Facette::Fac &my_elem) { my_elem.assemblage_vect(L_rhs); } );

    if (verbose)
        {
        std::cout << "matrix assembly done in " << counter.millis() << std::endl;
        counter.reset();
        }

    buildInitGuess(Xw);// gamma0 division handled by function buildInitGuess
    double residu = algebra::bicg<double>(iter, K, Xw, L_rhs);

    int nb_iter = iter.get_iteration();
    double _solver_error= iter.get_res();

    if( (nb_iter > MAXITER) || (_solver_error > TOL) )
        {
        if (verbose)
            {
            std::cout << "solver: bicgstab FAILED after " << nb_iter
            << " iterations, in " << counter.millis()<< "; residu= " << residu << std::endl;
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
        double vp = Xw[2*i];
        double vq = Xw[2*i+1];
        double v2 = Nodes::sq(vp) + Nodes::sq(vq);
        if (v2 > v2max)
            { v2max = v2; }
        
        refMsh->updateNode(i, vp, vq, dt);//gamma0 multiplication handled in updateNode
        }
    v_max = gamma0*sqrt(v2max);
    
    return 0;
    }
