#include "linear_algebra.h"

void LinAlgebra::base_projection()
    {
    std::mt19937 gen(rand());  // random number generator: standard Mersenne twister initialized
                               // with pseudo-random seed
    std::uniform_real_distribution<> distrib(0.0, 1.0);
    double r = distrib(gen);

    refMsh->setBasis(M_2_PI * r);
    }

void LinAlgebra::prepareElements(Pt::pt3D const &Hext /**< [in] applied field */,
                                 timing const &t_prm /**< [in] */)
    {
    base_projection();
    Pt::index idx_dir;

    if (!settings.recenter)
        {
        idx_dir = Pt::IDX_UNDEF;
        }
    else
        {
        idx_dir = settings.recentering_direction;
        }

    std::for_each(std::execution::par, refMsh->tet.begin(), refMsh->tet.end(),
                  [this, &Hext, &t_prm, idx_dir](Tetra::Tet &tet)
                  { tet.integrales(settings.paramTetra, t_prm, Hext, idx_dir, DW_vz); });

    std::for_each(std::execution::par, refMsh->fac.begin(), refMsh->fac.end(),
                  [this](Facette::Fac &fac)
                  { fac.integrales(settings.paramFacette); });
    }
