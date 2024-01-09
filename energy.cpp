#include "fem.h"

void Fem::energy(double const t, Settings &settings)
    {
    zeroEnergy();
    Eigen::Vector3d Hext = settings.getField(t);

    std::for_each(msh.tet.begin(), msh.tet.end(),
                  [this, &Hext, &settings](Tetra::Tet const &te)
                  {
                      Tetra::prm const &param = settings.paramTetra[te.idxPrm];
                      Eigen::Matrix<double,Pt::DIM,Tetra::NPI> u,dudx,dudy,dudz;
                      //double u[Pt::DIM][Tetra::NPI], dudx[Pt::DIM][Tetra::NPI], dudy[Pt::DIM][Tetra::NPI], dudz[Pt::DIM][Tetra::NPI];
                      Eigen::Vector<double,Tetra::NPI> phi;//double phi[Tetra::NPI];

                      te.interpolation(Nodes::get_u, u, dudx, dudy, dudz);
                      te.interpolation(Nodes::get_phi, phi);

                      E_exch += te.exchangeEnergy(param, dudx, dudy, dudz);

                      E_demag += te.demagEnergy(dudx, dudy, dudz, phi);

                      if ((param.K != 0.0) || (param.K3 != 0.0))
                          {
                          E_aniso += te.anisotropyEnergy(param, u);
                          }

                      E_zeeman += te.zeemanEnergy(param, 0, Hext, u);
                  });

    std::for_each(msh.fac.begin(), msh.fac.end(),
                  [this, &settings](Facette::Fac const &fa)
                  {
                      Facette::prm const &param = settings.paramFacette[fa.idxPrm];
                      Eigen::Vector<double,Facette::NPI> phi;
                      Eigen::Matrix<double,Pt::DIM,Facette::NPI> u;

                      fa.interpolation(Nodes::get_u, u);
                      fa.interpolation(Nodes::get_phi, phi);

                      if (param.Ks != 0.0)
                          {
                          E_aniso += fa.anisotropyEnergy(param, u);
                          }
                      
                      Eigen::Vector<double,Facette::NPI> _phi {phi[0],phi[1],phi[2],phi[3]};
                      E_demag += fa.demagEnergy(u, _phi);
                  });

    calc_Etot();

    if (settings.verbose && (Etot > Etot0))
        {
        std::cout << "WARNING: energy increased from " << Etot0 << " to " << Etot << "\n";
        }
    }
