#include "fem.h"

using namespace Nodes;

void Fem::energy(double const t, Settings &settings)
    {
    std::fill(E.begin(),E.end(),0);
    Etot = 0.0;
    Eigen::Vector3d Hext = settings.getField(t);

    std::for_each(msh.tet.begin(), msh.tet.end(),
                  [this, &Hext, &settings](Tetra::Tet const &te)
                  {
                      Tetra::prm const &param = settings.paramTetra[te.idxPrm];
                      Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> u,dudx,dudy,dudz;
                      Eigen::Matrix<double,Tetra::NPI,1> phi;

                      te.interpolation(Nodes::get_u<NEXT>, u, dudx, dudy, dudz);
                      te.interpolation(Nodes::get_phi<NEXT>, phi);

                      E[EXCHANGE] += te.exchangeEnergy(param, dudx, dudy, dudz);
                      E[DEMAG] += te.demagEnergy(param, dudx, dudy, dudz, phi);

                      if ((param.K != 0.0) || (param.K3 != 0.0))
                          { E[ANISOTROPY] += te.anisotropyEnergy(param, u); }

                      E[ZEEMAN] += te.zeemanEnergy(param, Hext, u);
                  });

    std::for_each(msh.fac.begin(), msh.fac.end(),
                  [this, &settings](Facette::Fac const &fa)
                  {
                      Facette::prm const &param = settings.paramFacette[fa.idxPrm];
                      Eigen::Matrix<double,Facette::NPI,1> phi;
                      Eigen::Matrix<double,Nodes::DIM,Facette::NPI> u;

                      fa.interpolation(Nodes::get_u<NEXT>, u);
                      fa.interpolation(Nodes::get_phi<NEXT>, phi);

                      if (param.Ks != 0.0)
                          { E[ANISOTROPY] += fa.anisotropyEnergy(param, u); }
                      
                      E[DEMAG] += fa.demagEnergy(u, phi);
                  });
    Etot = std::reduce(E.begin(),E.end());//sum of all terms
    if (settings.verbose && (Etot > Etot0))
        { std::cout << "WARNING: energy increased from " << Etot0 << " to " << Etot << "\n"; }
    }
