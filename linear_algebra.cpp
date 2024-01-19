#include "linear_algebra.h"

void LinAlgebra::base_projection()
    {
    std::mt19937 gen(rand());  // random number generator: standard Mersenne twister initialized
                               // with pseudo-random seed
    std::uniform_real_distribution<> distrib(0.0, 1.0);
    double r = distrib(gen);

    refMsh->setBasis(M_2_PI * r);
    }

void LinAlgebra::prepareElements(Eigen::Vector3d const &Hext /**< [in] applied field */,
                                 timing const &t_prm /**< [in] */)
    {
    base_projection();
    Nodes::index idx_dir;

    if (!settings.recenter)
        { idx_dir = Nodes::IDX_UNDEF; }
    else
        { idx_dir = settings.recentering_direction; }

    std::for_each(std::execution::par, refMsh->tet.begin(), refMsh->tet.end(),
                  [this, &Hext, &t_prm, idx_dir](Tetra::Tet &tet)
                  { tet.integrales(settings.paramTetra[tet.idxPrm], t_prm, Hext, idx_dir, DW_vz); });

    std::for_each(std::execution::par, refMsh->fac.begin(), refMsh->fac.end(),
                  [this](Facette::Fac &fac)
                  {
                  // the contribution to Lp computed in integrales is due to surface anisotropy
                  if (settings.paramFacette[fac.idxPrm].Ks != 0)
                    { fac.integrales(settings.paramFacette[fac.idxPrm]);  }
                  });
    }
