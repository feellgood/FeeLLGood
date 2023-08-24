#include "chronometer.h"
#include "linear_algebra.h"

int LinAlgebra::solver(timing const &t_prm, long nt)
    {
    chronometer counter(2);

    write_matrix K_TH(2 * NOD, 2 * NOD);
    std::for_each(refMsh->tet.begin(), refMsh->tet.end(),
                      [this,&K_TH](Tetra::Tet &my_elem) { my_elem.assemblage_mat(NOD,K_TH); } );
    std::for_each(refMsh->fac.begin(), refMsh->fac.end(),
                      [this,&K_TH](Facette::Fac &my_elem) { my_elem.assemblage_mat(NOD,K_TH); } );
    
    std::vector<double> L_TH(2 * NOD, 0);
    
    std::for_each(refMsh->tet.begin(), refMsh->tet.end(),
                      [this,&L_TH](Tetra::Tet &my_elem) { my_elem.assemblage_vect(NOD,L_TH); } );
    std::for_each(refMsh->fac.begin(), refMsh->fac.end(),
                      [this,&L_TH](Facette::Fac &my_elem) { my_elem.assemblage_vect(NOD,L_TH); } );
    
    if (settings.verbose)
        {
        std::cout << "matrix assembly done in " << counter.millis() << std::endl;
        }

    read_matrix Kr(2 * NOD, 2 * NOD);
    gmm::copy(K_TH, Kr);

    std::vector<double> X(2 * NOD);  // filled with zero by default
    refMsh->buildInitGuess(X); // gamma0 division handled within function buildInitGuess

    gmm::iteration bicg_iter(1e-6);
    bicg_iter.set_maxiter(settings.MAXITER);

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

    counter.reset();
    gmm::bicgstab(Kr, X, L_TH, *prc, bicg_iter);

    if (!(bicg_iter.converged()))
        {
        if (settings.verbose)
            {
            std::cout << "solver: bicg (prc " << (nt % (settings.REFRESH_PRC)) << ") FAILED after "
                      << bicg_iter.get_iteration() << " iterations, in " << counter.millis()
                      << std::endl;
            }
        return 1;
        }
    else
        {
        if (settings.verbose)
            {
            std::cout << "solver: bicg (prc " << (nt % (settings.REFRESH_PRC)) << ") converged in "
                      << bicg_iter.get_iteration() << " iterations, " << counter.millis()
                      << std::endl;
            }
        }

    // gamma0 multiplication handled within function updateNodes
    v_max = refMsh->updateNodes(X, t_prm.get_dt());
    return 0;
    }
