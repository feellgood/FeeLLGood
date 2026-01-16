#include "algebra/bicg.h"
#include "spinAccumulationSolver.h"
#include "chronometer.h" //date()
#include "tags.h"

using algebra::sq;
using namespace Nodes;

void spinAcc::checkBoundaryConditions(void) const
    {
    // nbVolP and nbVolN0 initialized to 1 because of __default__
    unsigned int nbVolP(1);
    unsigned int nbVolN0(1);
    std::for_each(paramTet.begin(),paramTet.end(),[&nbVolN0,&nbVolP](Tetra::prm const &p)
        {
        if(p.regName != "__default__")
            {
            if (std::isfinite(p.N0) && (p.N0 != 0)) nbVolN0++;
            if(std::isfinite(p.P) && (0 <= p.P) && (p.P < 1.0)) nbVolP++;
            }
        });
    int nbSurfJ(0);
    int nbSurfS(0);
    std::for_each(paramFac.begin(),paramFac.end(),[&nbSurfJ,&nbSurfS](Facette::prm const &p)
        {
        if(p.regName != "__default__")
            {
            if (std::isfinite(p.jn)) nbSurfJ++;
            if (std::isfinite(p.s.norm())) nbSurfS++;
            }
        });

    bool result = ( (nbSurfJ == 1)&&(nbSurfS == 1)
           && (nbVolN0 == paramTet.size())
           && (nbVolP == paramTet.size()) );

    if (!result)
        {
        std::cout << "Error: incorrect boundary conditions for spin diffusion solver.\n";
        exit(1);
        }
    else if (verbose)
        { std::cout << "spin diffusion problem boundary conditions Ok.\n"; }
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
    std::for_each(msh->fac.begin(),msh->fac.end(),[this](Facette::Fac &f)
        {
        if (std::isnan(paramFac[f.idxPrm].s.norm()) &&  std::isfinite(paramFac[f.idxPrm].jn))
            {
             /* units:
             * [jn] = A m^-2; [BOHRS_MUB/CHARGE_ELECTRON] = m^2 ; [P] = 1 ⇒ [s_value] = A
             * check s_value formula, especially the sign with current convention
             * */
            Eigen::Vector3d s_value = -paramFac[f.idxPrm].jn*(BOHRS_MUB/CHARGE_ELECTRON)*paramFac[f.idxPrm].uP;
            for(int j=0;j<Facette::N;j++)
                { fillDirichletData(f.ind[j],s_value); }
            }
        else if (std::isfinite(paramFac[f.idxPrm].s.norm()) &&  std::isnan(paramFac[f.idxPrm].jn))
            {
            Eigen::Vector3d s_value = paramFac[f.idxPrm].s;
            for(int j=0;j<Facette::N;j++)
                { fillDirichletData(f.ind[j],s_value); }
            }
        });
    suppress_copies<int>(idxDirichlet);
    }

double spinAcc::getMs(Tetra::Tet const &tet) const
    { return paramTet[tet.idxPrm].Ms; }

double spinAcc::getSigma(Tetra::Tet const &tet) const
    { return paramTet[tet.idxPrm].sigma; }

double spinAcc::getDiffusionCst(Tetra::Tet &tet) const
    {
    const double N0 = paramTet[tet.idxPrm].N0;
    return 2.0*getSigma(tet)/(sq(CHARGE_ELECTRON)*N0);
    }

double spinAcc::getPolarizationRate(Tetra::Tet &tet) const
    { return paramTet[tet.idxPrm].P; }

double spinAcc::getLsd(Tetra::Tet &tet) const
    { return paramTet[tet.idxPrm].lsd; }

double spinAcc::getLsf(Tetra::Tet &tet) const
    { return paramTet[tet.idxPrm].lsf; }

double spinAcc::getSpinHall(Tetra::Tet &tet) const
    { return paramTet[tet.idxPrm].spinHall; }

void spinAcc::setPotential(std::vector<double> &_V)
    { V = _V; }

void spinAcc::prepareExtras(void)
    {
    using namespace Tetra;

    std::for_each( msh->tet.begin(), msh->tet.end(), [this](Tet &t)
        {
        const double sigma = getSigma(t);
        const double P = getPolarizationRate(t);
        const double lsd = getLsd(t);
        const double lsf = getLsf(t);
        const double ksi = sq(lsd/lsf);
        const double Ms = getMs(t);

        // this formula might be mistaken, mixture of different models, to check
        double prefactor = BOHRS_MUB*P/(gamma0*Ms*CHARGE_ELECTRON*(1.0 + sq(ksi)));

        t.extraField = [this, _idx = t.idx](Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> H)
                         { for(int npi = 0; npi<Tetra::NPI; npi++) { H.col(npi) += Hst[_idx].col(npi); } };

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
                        BE.col(i) += ai_w*( Hst[t.idx].col(npi) + prefactor*m);
                        }
                    } // end loop on npi
                }; //end lambda
        });//end for_each
    }

void spinAcc::preCompute(void)
    {
    s.resize(NOD);
    std::for_each(msh->tet.begin(), msh->tet.end(), [this](Tetra::Tet const &tet)
                 {
                 Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> _gradV;
                 calc_gradV(tet, _gradV);
                 gradV.push_back(_gradV);

                 Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> _Hst;
                 calc_Hst(tet, _gradV, _Hst);
                 Hst.push_back(_Hst);
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

// Hst is a torque involving current density
void spinAcc::calc_Hst(Tetra::Tet const &tet, Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> _gradV,
                 Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> _Hst)
    {
    Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> p_g;
    tet.getPtGauss(p_g);
    const double sigma = getSigma(tet);
    for (int npi = 0; npi < Tetra::NPI; npi++)
        { _Hst.col(npi) = -sigma * _gradV.col(npi).cross(p_g.col(npi)); }
    }

bool spinAcc::solve(void)
    {
    iter.reset();
    algebra::w_sparseMat Kw(DIM_PB*NOD);
    std::vector<double> Lw(DIM_PB*NOD,0);

    std::for_each( msh->tet.begin(), msh->tet.end(),[this,&Kw,&Lw](Tetra::Tet &elem)
        {
        Eigen::Matrix<double,DIM_PB*Tetra::N,DIM_PB*Tetra::N> K;
        K.setZero();
        std::vector<double> L(DIM_PB*Tetra::N,0.0);
        integrales(elem,K);
        if(msh->isMagnetic(elem))
            {
            integraleMag(elem,K);
            integraleMag(elem,L);
            }
        if(paramTet[elem.idxPrm].spinHall != 0)
            { integraleSpinHall(elem,L); }
        buildMat<Tetra::N>(elem.ind, K, Kw);
        buildVect<Tetra::N>(elem.ind, L, Lw);
        } );

    // here are the boundary conditions
    std::for_each(msh->fac.begin(),msh->fac.end(),[this,&Lw](Facette::Fac &f)
        {
        std::vector<double> L(DIM_PB*Facette::N,0.0);
        Eigen::Vector3d s_value;
        if (!std::isfinite(paramFac[f.idxPrm].s.norm()) && std::isfinite(paramFac[f.idxPrm].jn))
            { // we should also test polarization P
            s_value = -paramFac[f.idxPrm].jn*(BOHRS_MUB/CHARGE_ELECTRON)*paramFac[f.idxPrm].uP;
            for (int npi=0; npi<Facette::NPI; npi++)
                {
                const double w = f.weight[npi];
                for (int ie=0; ie<Facette::N; ie++)
                    {
                    double ai_w = w*Facette::a[ie][npi];
                    L[               ie] -= s_value[IDX_X]* ai_w;
                    L[  Facette::N + ie] -= s_value[IDX_Y]* ai_w;
                    L[2*Facette::N + ie] -= s_value[IDX_Z]* ai_w;
                    }
                }
            }
        buildVect<Facette::N>(f.ind, L, Lw);
        });

    algebra::r_sparseMat Kr(Kw);
    std::vector<double> Xw(DIM_PB*NOD);
    algebra::bicg_dir(iter, Kr, Xw, Lw, valDirichlet, idxDirichlet);

    for (int i=0; i<NOD; i++)
        for (int j=0; j<3; j++)
        {
        s[i][j] = Xw[DIM_PB*i + j];
        }
    return (iter.status == algebra::CONVERGED);
    }

void spinAcc::integraleMag(Tetra::Tet &tet,
                           Eigen::Matrix<double,DIM_PB*Tetra::N,DIM_PB*Tetra::N> &AE)
    {
    using namespace Tetra;

    const double D0 = getDiffusionCst(tet);
    /* constant cst1 in a magnetic region is the only LHS parameter involved in the diffusion
     *  equation for magnetic contribution
     *  units: [cst1] = s^-1 : it is 1/tau_sd
     */
    const double cst1 = D0/sq(getLsd(tet));

    for (size_t npi=0; npi<NPI; npi++)
        {
        const double cst1_w = cst1*tet.weight[npi];
        for (size_t ie=0; ie<N; ie++)
            {
            const double tmp = cst1_w*Tetra::a[ie][npi];
            const Eigen::Vector3d &m = msh->getNode_u(tet.ind[ie]);
            AE(    ie,   N+ie) += tmp*m[2];
            AE(    ie, 2*N+ie) -= tmp*m[1];
            AE(  N+ie,     ie) -= tmp*m[2];
            AE(  N+ie, 2*N+ie) += tmp*m[0];
            AE(2*N+ie,     ie) += tmp*m[1];
            AE(2*N+ie,   N+ie) -= tmp*m[0];
            }
        }
    }

void spinAcc::integraleMag(Tetra::Tet &tet, std::vector<double> &BE)
    {
    using namespace Tetra;

    /* constant cst0 in a magnetic region is the only RHS parameter involved in the diffusion
     * equation for magnetic contribution
     * units: [cst0] = [sigma] m^2 = A^2 s^3 m^-1 kg^-1
     */
    const double cst0 = BOHRS_MUB*getPolarizationRate(tet)*getSigma(tet)/CHARGE_ELECTRON;
    Eigen::Matrix<double,Nodes::DIM,NPI> &_gradV = gradV[tet.idx];

    for (size_t npi=0; npi<NPI; npi++)
        {
        const double cst0_w = cst0*tet.weight[npi];

        for (size_t ie=0; ie<N; ie++)
            {
            const Eigen::Vector3d &m = msh->getNode_u(tet.ind[ie]);//magnetization
            Eigen::Vector3d grad_ai = tet.da.row(ie);
            double tmp = cst0_w*grad_ai.dot( _gradV.col(npi) );
            BE[    ie] += tmp*m[0];
            BE[  N+ie] += tmp*m[1];
            BE[2*N+ie] += tmp*m[2];
            }
        }
    }

void spinAcc::integraleSpinHall(Tetra::Tet &tet, std::vector<double> &BE)
    {
    using namespace Tetra;

    const double cst0 = getSpinHall(tet)*CHARGE_ELECTRON/MASS_ELECTRON;
    Eigen::Matrix<double,Nodes::DIM,NPI> &_gradV = gradV[tet.idx];
    for (size_t npi=0; npi<NPI; npi++)
        {
        const double cst0_w = cst0*tet.weight[npi];
        for (size_t ie=0; ie<N; ie++)
            {
            Eigen::Vector3d grad_ai = tet.da.row(ie);
            Eigen::Vector3d v = grad_ai.cross(_gradV.col(npi));
            BE[    ie] += cst0_w*v[IDX_X];
            BE[  N+ie] += cst0_w*v[IDX_Y];
            BE[2*N+ie] += cst0_w*v[IDX_Z];
            }
        }
    }

void spinAcc::integrales(Tetra::Tet &tet,
                         Eigen::Matrix<double,DIM_PB*Tetra::N,DIM_PB*Tetra::N> &AE)
    {//AE has a block structure, it might be usefull to remove loops and fill it with N*N blocks, as exchange
    using namespace Tetra;

    const double lsf = getLsf(tet);
    const double D0 = getDiffusionCst(tet);
    /* units: [D0/sq(lsf)] = s^-1 : it is 1/tau_sf */
    Eigen::Matrix<double,N,N> da_daT = tet.da*tet.da.transpose();

    for (size_t npi=0; npi<NPI; npi++)
        {
        const double D0_w = D0*tet.weight[npi];
        for (size_t ie=0; ie<N; ie++)
            {
            double tmp = D0_w*Tetra::a[ie][npi]/sq(lsf);
            AE(    ie,     ie) += tmp;
            AE(  N+ie,   N+ie) += tmp;
            AE(2*N+ie, 2*N+ie) += tmp;
            for (size_t je=0; je<N; je++)
                {
                AE(    ie,     je) += D0_w*da_daT(ie,je);
                AE(  N+ie,   N+je) += D0_w*da_daT(ie,je);
                AE(2*N+ie, 2*N+je) += D0_w*da_daT(ie,je);
                }
            }
        }
    }

