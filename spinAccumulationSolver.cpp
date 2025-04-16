#include "spinAccumulationSolver.h"

double spinAcc::getN0(Tetra::Tet &tet) const
    { return paramTetra[tet.idxPrm].N0; }

double spinAcc::getBeta(Tetra::Tet &tet) const
    { return paramTetra[tet.idxPrm].beta; }

double spinAcc::getLsd(Tetra::Tet &tet) const
    { return paramTetra[tet.idxPrm].lsd; }

double spinAcc::getLsf(Tetra::Tet &tet) const
    { return paramTetra[tet.idxPrm].lsf; }

void spinAcc::prepareExtras(std::vector<Tetra::Tet> &v_tet, electrostatSolver &elec)
    {
    using namespace Nodes;
    using namespace Tetra;

    std::for_each( v_tet.begin(), v_tet.end(), [this, &elec](Tet &tet)
        {
        const int _idx = tet.idx;
        const double sigma = elec.getSigma(tet);
        const double N0 = getN0(tet);
        const double D0 = 2.0*sigma / (sq(CHARGE_ELECTRON) * N0);
        const double beta = getBeta(tet);
        const double lsd = getLsd(tet);
        const double lsf = getLsf(tet);
        const double ksi = sq(lsd/lsf); /**< ksi is in Thiaville notations beta_DW */
        const double Js = 42;// for calc_prefactor(Js); what is that Js ?
        double prefactor = D0 / sq(lsd) / (gamma0*Js/mu0);
        prefactor *= sq(lsd) / (D0 * (1. + sq(ksi))) * BOHRS_MUB * beta / CHARGE_ELECTRON; // the composition of the two formulas simplifies, and is independant from D0, is this correct ?? to check

        tet.extraField = [this, &elec, sigma, _idx](Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> H)
                         { for(int npi = 0; npi<Tetra::NPI; npi++) { H.col(npi) += elec.Hm[_idx].col(npi); } };

        tet.extraCoeffs_BE = [this, &elec, sigma, ksi, prefactor, &tet](Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> U,
                                          Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> dUdx,
                                          Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> dUdy,
                                          Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> dUdz,
                                          Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,N>> BE)
                {
                for (int npi = 0; npi < NPI; npi++)
                    {
                    Eigen::Vector3d const &_gV = elec.gradV[tet.idx].col(npi);
                    Eigen::Vector3d j_grad_u =
                    -sigma * Eigen::Vector3d(_gV.dot( Eigen::Vector3d(dUdx(Nodes::IDX_X,npi), dUdy(Nodes::IDX_X,npi), dUdz(Nodes::IDX_X,npi))),
                                             _gV.dot( Eigen::Vector3d(dUdx(Nodes::IDX_Y,npi), dUdy(Nodes::IDX_Y,npi), dUdz(Nodes::IDX_Y,npi))),
                                             _gV.dot( Eigen::Vector3d(dUdx(Nodes::IDX_Z,npi), dUdy(Nodes::IDX_Z,npi), dUdz(Nodes::IDX_Z,npi))));

                    Eigen::Vector3d m = ksi * j_grad_u + U.col(npi).cross(j_grad_u);
                    for (int i = 0; i < N; i++)
                        {
                        const double ai_w = tet.weight[npi] * a[i][npi];
                        BE.col(i) += ai_w*( elec.Hm[tet.idx].col(npi) + prefactor*m);
                        }
                    } // end loop on npi
                }; //end lambda
        });//end for_each
    }
