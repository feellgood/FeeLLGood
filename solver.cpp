#include "chronometer.h"
#include "linear_algebra.h"

#include "algebra/bicg.h"

int LinAlgebra::solver(timing const &t_prm)
    {
    chronometer counter(2);
    iter.reset();

    K.clear();
    std::for_each(EXEC_POL, refMsh->tet.begin(), refMsh->tet.end(),
                      [this](Tetra::Tet &my_elem) { my_elem.assemblage_mat(K); } );

    if (verbose)
        {
        std::cout << "matrix assembly done in " << counter.millis() << std::endl;
        counter.reset();
        }

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
    algebra::bicg<double>(iter, K, Xw, L_rhs);

    if( (iter.status == algebra::ITER_OVERFLOW) || (iter.status == algebra::CANNOT_CONVERGE) || (iter.get_res() > iter.resmax) )
        {
        if (verbose)
            { std::cout << "solver: " << iter.infos() << std::endl; }
        return 1;
        }

    if (verbose)
        { std::cout << "solver: " << iter.infos() <<  " in " << counter.millis() << std::endl; }
    
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
