#include "electrostatSolver.h"

void integrales(Tetra::Tet const &tet, double sigma, double AE[Tetra::N][Tetra::N])
    {
    for (int npi = 0; npi < Tetra::NPI; npi++)
        {
        double w = tet.weight[npi];

        for (int ie = 0; ie < Tetra::N; ie++)
            {
            double dai_dx = tet.da(ie,Nodes::IDX_X);
            double dai_dy = tet.da(ie,Nodes::IDX_Y);
            double dai_dz = tet.da(ie,Nodes::IDX_Z);

            for (int je = 0; je < Tetra::N; je++)
                {
                double daj_dx = tet.da(je,Nodes::IDX_X);
                double daj_dy = tet.da(je,Nodes::IDX_Y);
                double daj_dz = tet.da(je,Nodes::IDX_Z);
                AE[ie][je] += sigma * (dai_dx * daj_dx + dai_dy * daj_dy + dai_dz * daj_dz) * w;
                }
            }
        }
    }

void integrales(Facette::Fac const &fac, double pot_val, std::vector<double> &BE)
    {
    for (int npi = 0; npi < Facette::NPI; npi++)
        for (int ie = 0; ie < Facette::N; ie++)
            {
            BE[ie] -= Facette::a[ie][npi] * pot_val * fac.weight[npi];
            }
    }
