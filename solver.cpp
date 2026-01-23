#include "chronometer.h"
#include "linear_algebra.h"

#include "algebra/bicg.h"

bool LinAlgebra::solve(timing const &t_prm)
    {
    chronometer counter(2);
    iter.reset();

    K.clear();
    std::for_each(EXEC_POL, msh->magTet.begin(), msh->magTet.end(),
                  [this](const int idx)
                      {
                      Tetra::Tet &my_elem = msh->tet[idx];
                      buildMat<Tetra::N>(my_elem.ind, my_elem.Kp);
                      });

    if (verbose)
        {
        std::cout << "matrix assembly done in " << counter.millis() << std::endl;
        counter.reset();
        }

    std::fill(L_rhs.begin(),L_rhs.end(),0);
    std::for_each(msh->magTet.begin(), msh->magTet.end(),
                  [this](const int idx)
                      {
                      Tetra::Tet &my_elem = msh->tet[idx];
                      double *data_start = my_elem.Lp.data();
                      double *data_end = data_start + 2*Tetra::N;
                      std::vector<double> Le(data_start, data_end);
                      buildVect<Tetra::N>(my_elem.ind, Le);
                      });

    std::for_each(msh->magFac.begin(), msh->magFac.end(),
                  [this](const int idx)
                      {
                      Facette::Fac &my_elem = msh->fac[idx];
                      double *data_start = my_elem.Lp.data();
                      double *data_end = data_start + 2*Facette::N;
                      std::vector<double> Le(data_start, data_end);
                      buildVect<Facette::N>(my_elem.ind, Le);
                      });

    /* RHS forced to zero outside mag material defined by mask, corresponding K diagonal coefficients set to 1  */
    algebra::applyMask(lvd,L_rhs);
    std::for_each(lvd.begin(),lvd.end(),[this](int i){ K.set(i,i,1.0); });
    if (verbose)
        {
        std::cout << "matrix assembly done in " << counter.millis() << std::endl;
        counter.reset();
        }
    /*
     *bicg with Dirichlet boundary conditions is used to mask non magnetic mesh regions, the
     Dirichlet values are zeros, so a specific bicg_dir is called assuming those zeros, instead of
     the standard one.
     * */
    buildInitGuess(Xw);// gamma0 division handled by function buildInitGuess
    algebra::bicg_dir<double>(iter, K, Xw, L_rhs, lvd);

    if( (iter.status == algebra::ITER_OVERFLOW) || (iter.status == algebra::CANNOT_CONVERGE) || (iter.get_res() > iter.resmax) )
        {
        if (verbose)
            { std::cout << "solver: " << iter.infos() << std::endl; }
        return true;
        }

    if (verbose)
        { std::cout << "solver: " << iter.infos() <<  " in " << counter.millis() << std::endl; }

    double v2max(0.0);
    const double dt = t_prm.get_dt();
    for (int i = 0; i < NOD; i++)
        {
        if (msh->magNode[i])
            {
            double vp = Xw[2*i];
            double vq = Xw[2*i+1];
            double v2 = Nodes::sq(vp) + Nodes::sq(vq);
            if (v2 > v2max)
                { v2max = v2; }
            msh->updateNode(i, vp, vq, dt);//gamma0 multiplication handled in updateNode
            }
        }
    v_max = gamma0*sqrt(v2max);
    return false;
    }
