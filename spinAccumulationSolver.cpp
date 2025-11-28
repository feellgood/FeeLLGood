#include "algebra/bicg.h"
#include "spinAccumulationSolver.h"
#include "chronometer.h" //date()
#include "tags.h"
#include "solverUtils.h"

using algebra::sq;
using namespace Nodes;

bool spinAcc::checkBoundaryConditions(bool verbose) const
    {
    // nbVolP and nbVolN0 initialized to 1 because of __default__
    unsigned int nbVolP(1);
    unsigned int nbVolN0(1);
    std::for_each(paramTetra.begin(),paramTetra.end(),[&nbVolN0,&nbVolP](Tetra::prm const &p)
        {
        if(p.regName != "__default__")
            {
            if (std::isfinite(p.N0) && (p.N0 != 0)) nbVolN0++;
            if(std::isfinite(p.P) && (0 <= p.P) && (p.P < 1.0)) nbVolP++;
            }
        });
    int nbSurfJ(0);
    int nbSurfS(0);
    std::for_each(paramFacette.begin(),paramFacette.end(),[&nbSurfJ,&nbSurfS,verbose](Facette::prm const &p)
        {
        if(p.regName != "__default__")
            {
            if (verbose)
                {
                std::cout << p.regName << " : J= " << p.J << "; s= {" << p.s[0] << ';' << p.s[1]
                          << ';' << p.s[2] << "}\n";
                }
            if (std::isfinite(p.J)) nbSurfJ++;
            if (std::isfinite(p.s.norm())) nbSurfS++;
            }
        });
    return ( (nbSurfJ == 1)&&(nbSurfS == 1)
           && (nbVolN0 == paramTetra.size())
           && (nbVolP == paramTetra.size()) );
    }

void spinAcc::fillDirichletData(const int k, Eigen::Vector3d &s_value)
    {
    valDirichlet[3*k] = s_value[Nodes::IDX_X];
    idxDirichlet.push_back(3*k);
    valDirichlet[3*k + 1] = s_value[Nodes::IDX_Y];
    idxDirichlet.push_back(3*k + 1);
    valDirichlet[3*k + 2] = s_value[Nodes::IDX_Z];
    idxDirichlet.push_back(3*k + 2);
    }

void spinAcc::boundaryConditions(void)
    {
    std::fill(valDirichlet.begin(),valDirichlet.end(),0.0);
    std::for_each(refMsh->fac.begin(),refMsh->fac.end(),[this](Facette::Fac &f)
        {
        if (std::isnan(paramFacette[f.idxPrm].s.norm()) &&  std::isfinite(paramFacette[f.idxPrm].J))
            {// TODO: check that formula, especially the sign with current convention
            Eigen::Vector3d s_value = -paramFacette[f.idxPrm].J*(BOHRS_MUB/CHARGE_ELECTRON)*paramFacette[f.idxPrm].P;
            for(int j=0;j<Facette::N;j++)
                { fillDirichletData(f.ind[j],s_value); }
            }
        else if (std::isfinite(paramFacette[f.idxPrm].s.norm()) &&  std::isnan(paramFacette[f.idxPrm].J))
            {
            Eigen::Vector3d s_value = paramFacette[f.idxPrm].s;
            for(int j=0;j<Facette::N;j++)
                { fillDirichletData(f.ind[j],s_value); }
            }
        });
    suppress_copies<int>(idxDirichlet);
    }

double spinAcc::getJs(Tetra::Tet const &tet) const
    { return paramTetra[tet.idxPrm].J; }

double spinAcc::getSigma(Tetra::Tet const &tet) const
    { return paramTetra[tet.idxPrm].sigma; }

double spinAcc::getN0(Tetra::Tet &tet) const
    { return paramTetra[tet.idxPrm].N0; }

double spinAcc::getPolarization(Tetra::Tet &tet) const
    { return paramTetra[tet.idxPrm].P; }

double spinAcc::getLsd(Tetra::Tet &tet) const
    { return paramTetra[tet.idxPrm].lsd; }

double spinAcc::getLsf(Tetra::Tet &tet) const
    { return paramTetra[tet.idxPrm].lsf; }

double spinAcc::getSpinHall(Tetra::Tet &tet) const
    { return paramTetra[tet.idxPrm].spinHall; }

void spinAcc::setPotential(std::vector<double> &_V)
    { V = _V; }

void spinAcc::prepareExtras(void)
    {
    using namespace Tetra;

    std::for_each( refMsh->tet.begin(), refMsh->tet.end(), [this](Tet &t)
        {
        const int _idx = t.idx;
        const double sigma = getSigma(t);
        const double P = getPolarization(t);
        const double lsd = getLsd(t);
        const double lsf = getLsf(t);
        const double ksi = sq(lsd/lsf);
        const double Js = getJs(t);// Js = Ms/nu_0
        double prefactor = mu0*BOHRS_MUB*P/(gamma0*Js*CHARGE_ELECTRON*(1.0 + sq(ksi)));

        t.extraField = [this, _idx](Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> H)
                         { for(int npi = 0; npi<Tetra::NPI; npi++) { H.col(npi) += Hm[_idx].col(npi); } };

        t.extraCoeffs_BE = [this, sigma, ksi, prefactor, &t](Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> U,
                                          Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> dUdx,
                                          Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> dUdy,
                                          Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> dUdz,
                                          Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,N>> BE)
                {
                for (int npi = 0; npi < NPI; npi++)
                    {
                    Eigen::Vector3d const &_gV = gradV[t.idx].col(npi);
                    Eigen::Vector3d j_grad_u =
                    -sigma * Eigen::Vector3d(_gV.dot( Eigen::Vector3d(dUdx(IDX_X,npi), dUdy(IDX_X,npi), dUdz(IDX_X,npi))),
                                             _gV.dot( Eigen::Vector3d(dUdx(IDX_Y,npi), dUdy(IDX_Y,npi), dUdz(IDX_Y,npi))),
                                             _gV.dot( Eigen::Vector3d(dUdx(IDX_Z,npi), dUdy(IDX_Z,npi), dUdz(IDX_Z,npi))));

                    Eigen::Vector3d m = ksi * j_grad_u + U.col(npi).cross(j_grad_u);
                    for (int i = 0; i < N; i++)
                        {
                        const double ai_w = t.weight[npi] * a[i][npi];
                        BE.col(i) += ai_w*( Hm[t.idx].col(npi) + prefactor*m);
                        }
                    } // end loop on npi
                }; //end lambda
        });//end for_each
    }

void spinAcc::preCompute(void)
    {
    s.resize(NOD);
    std::for_each(refMsh->tet.begin(), refMsh->tet.end(), [this](Tetra::Tet const &tet)
                 {
                 Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> _gradV;
                 calc_gradV(tet, _gradV);
                 gradV.push_back(_gradV);

                 Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> _Hm;
                 calc_Hm(tet, _gradV, _Hm);
                 Hm.push_back(_Hm);
                 });
    prepareExtras();
    }

bool spinAcc::compute(void)
    {
    bool has_converged = solve();
    if (!has_converged)
        {
        std::cout << "spin accumulation solver: " << iter.infos() << std::endl;
        for(int i=0;i<NOD;i++)
            { s[i].setZero(); }
        }
    else
        { if (verbose) { std::cout << "spin accumulation solved.\n"; } }
    return has_converged;
    }

void spinAcc::calc_gradV(Tetra::Tet const &tet, Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> _gradV)
    {
    for (int npi = 0; npi < Tetra::NPI; npi++)
        {
        Eigen::Vector3d v(0,0,0);
        for (int i = 0; i < Tetra::N; i++)
            { v += V[tet.ind[i]] * tet.da.row(i); }
        _gradV.col(npi) = v;
        }
    }

// Hm is a torque involving current density
void spinAcc::calc_Hm(Tetra::Tet const &tet, Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> _gradV,
                 Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> _Hm)
    {
    Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> p_g;
    tet.getPtGauss(p_g);
    const double sigma = getSigma(tet);
    for (int npi = 0; npi < Tetra::NPI; npi++)
        { _Hm.col(npi) = -sigma * _gradV.col(npi).cross(p_g.col(npi)); }
    }

bool spinAcc::solve(void)
    {
    iter.reset();
    algebra::w_sparseMat Kw(DIM_PROBLEM*NOD);
    std::vector<double> Lw(DIM_PROBLEM*NOD,0);

    std::for_each( refMsh->tet.begin(), refMsh->tet.end(),[this,&Kw,&Lw](Tetra::Tet &elem)
        {
        Eigen::Matrix<double,DIM_PROBLEM*Tetra::N,DIM_PROBLEM*Tetra::N> K;
        K.setZero();
        std::vector<double> L(DIM_PROBLEM*Tetra::N,0.0);
        integrales(elem,K);
        if(refMsh->isMagnetic(elem))
            { integraleMag(elem,K,L); }
        if(paramTetra[elem.idxPrm].spinHall != 0)
            { integraleSpinHall(elem,L); }
        buildMat<Tetra::N,DIM_PROBLEM>(elem.ind, K, Kw);
        buildVect<Tetra::N,DIM_PROBLEM>(elem.ind, L, Lw);
        } );

    // here are the boundary conditions
    std::for_each(refMsh->fac.begin(),refMsh->fac.end(),[this,&Lw](Facette::Fac &f)
        {
        std::vector<double> L(DIM_PROBLEM*Facette::N,0.0);
        Eigen::Vector3d s_value;
        if (paramFacette[f.idxPrm].s.norm() == 0)
            s_value.setZero();
        else if (!std::isfinite( paramFacette[f.idxPrm].s.norm() ) &&
                std::isfinite(paramFacette[f.idxPrm].J)) // we should also test polarization P
            s_value = -paramFacette[f.idxPrm].J*(BOHRS_MUB/CHARGE_ELECTRON)*paramFacette[f.idxPrm].P;
        integrales(f,s_value,L);
        buildVect<Facette::N,DIM_PROBLEM>(f.ind, L, Lw);
        });

    algebra::r_sparseMat Kr(Kw);
    std::vector<double> Xw(DIM_PROBLEM*NOD);
    algebra::bicg_dir(iter, Kr, Xw, Lw, valDirichlet, idxDirichlet);

    for (int i=0; i<NOD; i++)
        for (int j=0; j<3; j++)
        {
        s[i][j] = Xw[DIM_PROBLEM*i + j];
        }
    return (iter.status == algebra::CONVERGED);
    }

void spinAcc::integraleMag(Tetra::Tet &tet,
                           Eigen::Matrix<double,DIM_PROBLEM*Tetra::N,DIM_PROBLEM*Tetra::N> &AE,
                           std::vector<double> &BE)
    {
    using namespace Tetra;

    double sigma = getSigma(tet);
    double N0 = getN0(tet);
    double P = getPolarization(tet);
    double lsd = getLsd(tet);
    double D0=2.0*sigma/(sq(CHARGE_ELECTRON)*N0);
    /* these two constants in a magnetic region are the only parameters involved in the diffusion
       equation for magnetic contribution */
    double cst0 = BOHRS_MUB*P*sigma/CHARGE_ELECTRON;
    double cst1 = D0/sq(lsd);
    Eigen::Matrix<double,N,1> V_nod;
    for (size_t ie=0; ie<N; ie++) { V_nod[ie] = V[ tet.ind[ie] ]; }
    Eigen::Matrix<double,Nodes::DIM,NPI> &_gradV = gradV[tet.idx];

    for (size_t npi=0; npi<NPI; npi++)
        {
        double w = tet.weight[npi];

        for (size_t ie=0; ie<N; ie++)
            {
            const Eigen::Vector3d &m = refMsh->getNode_u(tet.ind[ie]);//magnetization
            Eigen::Vector3d grad_ai = tet.da.row(ie);
            double tmp = cst0*w*grad_ai.dot( _gradV.col(npi) );
            BE[    ie] += tmp*m[0];
            BE[  N+ie] += tmp*m[1];
            BE[2*N+ie] += tmp*m[2];

            tmp = cst1*w*Tetra::a[ie][npi];
            AE(    ie,   N+ie) += tmp*m[2];
            AE(    ie, 2*N+ie) -= tmp*m[1];
            AE(  N+ie,     ie) -= tmp*m[2];
            AE(  N+ie, 2*N+ie) += tmp*m[0];
            AE(2*N+ie,     ie) += tmp*m[1];
            AE(2*N+ie,   N+ie) -= tmp*m[0];
            }
        }
    }

void spinAcc::integraleSpinHall(Tetra::Tet &tet, std::vector<double> &BE)
    {
    using namespace Tetra;
    double spinHall = getSpinHall(tet);
    double cst0 = spinHall*CHARGE_ELECTRON/MASS_ELECTRON;
    Eigen::Matrix<double,N,1> V_nod;
    for (size_t ie=0; ie<N; ie++) { V_nod[ie] = V[ tet.ind[ie] ]; }
    Eigen::Matrix<double,Nodes::DIM,NPI> &_gradV = gradV[tet.idx];
    for (size_t npi=0; npi<NPI; npi++)
        {
        double w = tet.weight[npi];
        for (size_t ie=0; ie<N; ie++)
            {
            Eigen::Vector3d grad_ai = tet.da.row(ie);
            Eigen::Vector3d v = grad_ai.cross(_gradV.col(npi));
            BE[    ie] += cst0*w*v[IDX_X];
            BE[  N+ie] += cst0*w*v[IDX_Y];
            BE[2*N+ie] += cst0*w*v[IDX_Z];
            }
        }
    }

void spinAcc::integrales(Tetra::Tet &tet,
                         Eigen::Matrix<double,DIM_PROBLEM*Tetra::N,DIM_PROBLEM*Tetra::N> &AE)
    {
    using namespace Tetra;

    double sigma = getSigma(tet);
    double N0 = getN0(tet);
    double lsf = getLsf(tet);
    double D0=2.0*sigma/(sq(CHARGE_ELECTRON)*N0);
    double cst0 = D0/sq(lsf);
    for (size_t npi=0; npi<NPI; npi++)
        {
        double w = tet.weight[npi];
        for (size_t ie=0; ie<N; ie++)
            {
            Eigen::Vector3d grad_ai = tet.da.row(ie);
            double tmp = cst0*w*Tetra::a[ie][npi];
            AE(    ie,     ie) += tmp;
            AE(  N+ie,   N+ie) += tmp;
            AE(2*N+ie, 2*N+ie) += tmp;
            for (size_t je=0; je<N; je++)
                {
                tmp = D0*w*(grad_ai.dot( tet.da.row(je) ));
                AE(    ie,     je) += tmp;
                AE(  N+ie,   N+je) += tmp;
                AE(2*N+ie, 2*N+je) += tmp;
                }
            }
        }
    }

void spinAcc::integrales(Facette::Fac &fac, Eigen::Vector3d &Q, std::vector<double> &BE)
    {
    using namespace Nodes;
    using namespace Facette;

    for (int npi=0; npi<NPI; npi++)
        {
        double w = fac.weight[npi];
        for (int ie=0; ie<N; ie++)
            {
            double ai_w = w*a[ie][npi];
            BE[    ie] -= Q[IDX_X]* ai_w;
            BE[  N+ie] -= Q[IDX_Y]* ai_w;
            BE[2*N+ie] -= Q[IDX_Z]* ai_w;
            }
        }
    }

