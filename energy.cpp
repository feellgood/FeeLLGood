#include "fem.h"

void Fem::energy(double const t, Settings &settings)
    {
    zeroEnergy();
    const Pt::pt3D Hext = settings.getField(t);

    std::for_each(msh.tet.begin(), msh.tet.end(),
                  [this, &Hext, &settings](Tetra::Tet const &te)
                  {
                      Tetra::prm const &param = settings.paramTetra[te.idxPrm];
                      double u[Pt::DIM][Tetra::NPI], dudx[Pt::DIM][Tetra::NPI],
                              dudy[Pt::DIM][Tetra::NPI], dudz[Pt::DIM][Tetra::NPI];
                      double phi[Tetra::NPI];

                      te.interpolation(Nodes::get_u, u, dudx, dudy, dudz);
                      te.interpolation(Nodes::get_phi, phi);

                      E_exch += te.exchangeEnergy(param, dudx, dudy, dudz);

                      E_demag += te.demagEnergy(dudx, dudy, dudz, phi);

                      if ((param.K != 0.0) || (param.K3 != 0.0))
                          {
                          E_aniso += te.anisotropyEnergy(param, u);
                          }

                      E_zeeman += te.zeemanEnergy(param, 0, Hext,u);//uz_drift, Hext, u);
                  });

    std::for_each(msh.fac.begin(), msh.fac.end(),
                  [this, &settings](Facette::Fac const &fa)
                  {
                      Facette::prm const &param = settings.paramFacette[fa.idxPrm];
                      double phi[Facette::NPI];
                      Pt::pt3D u[Facette::NPI];

                      fa.interpolation<Pt::pt3D>(Nodes::get_u, u);
                      fa.interpolation<double>(Nodes::get_phi, phi);

                      if (param.Ks != 0.0)
                          {
                          E_aniso += fa.anisotropyEnergy(param, u);
                          }
                      
                      //devNote: these copy will be removed when interpolation rewritten with eigen
                      Eigen::Matrix<double,Pt::DIM,Facette::NPI> _u;
                      for (int i=0;i<Pt::DIM;i++)
                        for(int j=0;j<Facette::NPI;j++) _u(i,j) = u[j](i);
                      
                      Eigen::Vector<double,Facette::NPI> _phi {phi[0],phi[1],phi[2],phi[3]};
                      E_demag += fa.demagEnergy(_u, _phi);
                  });

    calc_Etot();

    if (settings.verbose && (Etot > Etot0))
        {
        std::cout << "WARNING: energy increased from " << Etot0 << " to " << Etot << "\n";
        }
    }
