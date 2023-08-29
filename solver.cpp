#include "chronometer.h"
#include "linear_algebra.h"

int LinAlgebra::solver(timing const &t_prm, long nt)
    {
    chronometer counter(2);

    write_matrix K_TH(2 * NOD, 2 * NOD);
    std::vector<double> L_TH(2 * NOD, 0);

    insertCoeff<Tetra::Tet>(refMsh->tet, K_TH, L_TH);
    insertCoeff<Facette::Fac>(refMsh->fac, K_TH, L_TH);

    if (settings.verbose)
        {
        std::cout << "matrix assembly done in " << counter.millis() << std::endl;
        }

    read_matrix Kr(2 * NOD, 2 * NOD);
    gmm::copy(K_TH, Kr);

    std::vector<double> X(2 * NOD);  // filled with zero by default
    refMsh->buildInitGuess(X); // gamma0 division handled within function buildInitGuess

    gmm::iteration bicg_iter(1e-6);
    gmm::iteration gmr_iter(1e-6);
    bicg_iter.set_maxiter(settings.MAXITER);
    gmr_iter.set_maxiter(settings.MAXITER);

    // Let GMM be noisy only at verbosity levels higher than one.
    bicg_iter.set_noisy(settings.verbose > 1);
    gmr_iter.set_noisy(settings.verbose > 1);

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
    // gmm::clear(X); // no need to clear X,it is zero initialized
    gmm::bicgstab(Kr, X, L_TH, *prc, bicg_iter);

    if (!(bicg_iter.converged()))
        {
        if (settings.verbose)
            {
            std::cout << "solver: bicg (prc " << (nt % (settings.REFRESH_PRC)) << ") FAILED after "
                      << bicg_iter.get_iteration() << " iterations, in " << counter.millis()
                      << std::endl;
            }
        gmm::diagonal_precond<read_matrix> gmr_prc(Kr);
        counter.reset();

        gmm::clear(X);  // X is cleared for safety, it should be useless
        gmm::gmres(Kr, X, L_TH, gmr_prc, 50, gmr_iter);

        if (!(gmr_iter.converged()))
            {
            if (settings.verbose)
                {
                std::cout << "solver: gmres FAILED after " << gmr_iter.get_iteration()
                          << " iterations, in " << counter.millis() << std::endl;
                }
            return 1;
            }
        else
            {
            if (settings.verbose)
                {
                std::cout << "solver: gmres converged in " << gmr_iter.get_iteration()
                          << " iterations, " << counter.millis() << std::endl;
                }
            }
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

    updateNodes(X, t_prm.get_dt());// gamma0 multiplication handled within function updateNodes
    return 0;
    }
