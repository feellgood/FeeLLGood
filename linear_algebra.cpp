#include "linear_algebra.h"

void LinAlgebra::base_projection()
    {
    std::mt19937 gen(rand());  // random number generator: standard Mersenne twister initialized
                               // with pseudo-random seed
    std::uniform_real_distribution<> distrib(0.0, 1.0);
    double r = distrib(gen);

    msh->setBasis(M_2_PI * r);
    }

void LinAlgebra::buildInitGuess(std::vector<double> &G) const
    {
    std::fill(G.begin(),G.end(),0);
    for (int i = 0; i < NOD; i++)
        {
        if (msh->magNode[i])
            {
            G[2*i] = msh->getProj_ep(i)/gamma0;
            G[2*i+1] = msh->getProj_eq(i)/gamma0;
            }
        }
    }

void LinAlgebra::prepareElements(Eigen::Vector3d const &Hext /**< [in] applied field */,
                                 timing const &t_prm /**< [in] */)
    {
    // it might be more efficient to build calc_Hext inside lambda of the for_each
    std::function< Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>(void)> calc_Hext =
                  [&Hext](void)
                      {
                      Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> H;
                      H.colwise() = Hext;  // set all columns of H to Hext
                      return H;
                      };

    std::for_each(EXEC_POL, msh->tet.begin(), msh->tet.end(),
                  [this, &calc_Hext, &t_prm](Tetra::Tet &tet)
                  {
                  if (msh->isMagnetic(tet))
                      { tet.integrales(paramTet[tet.idxPrm], t_prm, calc_Hext, idx_dir, DW_vz); }
                  });

    std::for_each(EXEC_POL, msh->fac.begin(), msh->fac.end(),
                  [this](Facette::Fac &fac)
                  {
                  // the contribution to Lp computed in integrales is due to surface anisotropy
                  if ((msh->isMagnetic(fac))&&(paramFac[fac.idxPrm].Ks != 0))
                      { fac.integrales(paramFac[fac.idxPrm]);  }
                  });
    }

void LinAlgebra::prepareElements(double const A_Hext /**< [in] amplitude applied field (might be time dependant)*/,
                                 timing const &t_prm /**< [in] */)
    {
    std::for_each(EXEC_POL, msh->tet.begin(), msh->tet.end(),
                  [this, &A_Hext, &t_prm](Tetra::Tet &tet)
                  {
                  if (msh->isMagnetic(tet))
                      {
                      Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> sp_H = extSpaceField[tet.idx];
                      std::function< Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> (void)> calc_Hext =
                      [&sp_H,&A_Hext](void)
                          {
                          Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> H= A_Hext*sp_H;
                          return H;
                          };
                      tet.integrales(paramTet[tet.idxPrm], t_prm, calc_Hext, idx_dir, DW_vz);
                      }
                  });

    std::for_each(EXEC_POL, msh->fac.begin(), msh->fac.end(),
                  [this](Facette::Fac &fac)
                  {
                  // the contribution to Lp computed in integrales is due to surface anisotropy
                  if ((msh->isMagnetic(fac))&&(paramFac[fac.idxPrm].Ks != 0))
                    { fac.integrales(paramFac[fac.idxPrm]);  }
                  });
    }


void LinAlgebra::setExtSpaceField(Settings &s /**< [in] */)
    { // see here for ref code /data/jc/st_feellgood_2024/src_Tube_scalfmm_zhang_ec_mu_oersted_spinHall_thiele_dyn20240320_dev
    extSpaceField.resize(msh->magTet.size());
    int k(0);
    std::for_each(msh->magTet.begin(), msh->magTet.end(), [this,&s,&k](const int idx)
        {
        Tetra::Tet const &t = msh->tet[idx];
        Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> pg;// gauss points
        t.getPtGauss(pg);
        for(int i=0;i<Tetra::NPI;i++)
            { extSpaceField[k].col(i) = s.getFieldSpace(pg.col(i)); }
        k++;
        });
    }
