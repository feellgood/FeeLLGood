#include "linear_algebra.h"

void LinAlgebra::base_projection()
    {
    std::mt19937 gen(rand());  // random number generator: standard Mersenne twister initialized
                               // with pseudo-random seed
    std::uniform_real_distribution<> distrib(0.0, 1.0);
    double r = distrib(gen);

    refMsh->setBasis(M_2_PI * r);
    }

void LinAlgebra::buildInitGuess(algebra::Vector<double> &G) const
    {
    for (int i = 0; i < NOD; i++)
        {
        G[i] = refMsh->getProj_ep(i)/gamma0;
        G[NOD + i] = refMsh->getProj_eq(i)/gamma0;
        }
    }

void LinAlgebra::prepareElements(Eigen::Vector3d const &Hext /**< [in] applied field */,
                                 timing const &t_prm /**< [in] */)
    {
    // it might be more efficient to build calc_Hext inside lambda of the for_each
    std::function<void( Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> H )> calc_Hext =
                  [&Hext](Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> H) { H.colwise() = Hext; }; // set all columns of H to Hext

    std::for_each(EXEC_POL, refMsh->tet.begin(), refMsh->tet.end(),
                  [this, &calc_Hext, &t_prm](Tetra::Tet &tet)
                  { tet.integrales(prmTetra[tet.idxPrm], t_prm, calc_Hext, idx_dir, DW_vz); });

    std::for_each(EXEC_POL, refMsh->fac.begin(), refMsh->fac.end(),
                  [this](Facette::Fac &fac)
                  {
                  // the contribution to Lp computed in integrales is due to surface anisotropy
                  if (prmFacette[fac.idxPrm].Ks != 0)
                    { fac.integrales(prmFacette[fac.idxPrm]);  }
                  });
    }

void LinAlgebra::prepareElements(double const A_Hext /**< [in] amplitude applied field */, timing const &t_prm /**< [in] */)
    {
    std::for_each(EXEC_POL, refMsh->tet.begin(), refMsh->tet.end(),
                  [this, &A_Hext, &t_prm](Tetra::Tet &tet)
                  {
                  Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> sp_H = extSpaceField[tet.idx];
                  std::function<void( Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> H )> calc_Hext =
                  [&sp_H,&A_Hext](Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> H)
                        { H = A_Hext*sp_H; };

                  tet.integrales(prmTetra[tet.idxPrm], t_prm, calc_Hext, idx_dir, DW_vz);
                  });

    std::for_each(EXEC_POL, refMsh->fac.begin(), refMsh->fac.end(),
                  [this](Facette::Fac &fac)
                  {
                  // the contribution to Lp computed in integrales is due to surface anisotropy
                  if (prmFacette[fac.idxPrm].Ks != 0)
                    { fac.integrales(prmFacette[fac.idxPrm]);  }
                  });
    }


void LinAlgebra::setExtSpaceField(Settings &s /**< [in] */)
    { // see here for ref code /data/jc/st_feellgood_2024/src_Tube_scalfmm_zhang_ec_mu_oersted_spinHall_thiele_dyn20240320_dev
    extSpaceField.resize(refMsh->tet.size());
    int k(0);
    std::for_each(refMsh->tet.begin(), refMsh->tet.end(), [this,&s,&k](Tetra::Tet &tet)
        {
        Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> pg;// gauss points
        tet.getPtGauss(pg);
        for(int i=0;i<Tetra::NPI;i++)
            { extSpaceField[k].col(i) = s.getFieldSpace(pg.col(i)); }
        k++;
        });
    }
