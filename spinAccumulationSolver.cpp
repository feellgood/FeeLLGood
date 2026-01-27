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
    valDirichlet[DIM_PB*k] = s_value[Nodes::IDX_X];
    idxDirichlet.push_back(DIM_PB*k);
    valDirichlet[DIM_PB*k + 1] = s_value[Nodes::IDX_Y];
    idxDirichlet.push_back(DIM_PB*k + 1);
    valDirichlet[DIM_PB*k + 2] = s_value[Nodes::IDX_Z];
    idxDirichlet.push_back(DIM_PB*k + 2);
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

void spinAcc::preCompute(std::vector<double> &V)
    {
    s.resize(NOD);
    std::for_each(msh->tet.begin(), msh->tet.end(), [this,&V](Tetra::Tet const &tet)
                 {
                 Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> _gradV = calc_gradV(tet,V);
                 gradV.push_back(_gradV);

                 Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> _Hst = calc_Hst(tet, getSigma(tet), _gradV);
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

bool spinAcc::solve(void)
    {
    iter.reset();

    std::for_each( msh->tet.begin(), msh->tet.end(),[this](Tetra::Tet &elem)
        {
        Eigen::Matrix<double,DIM_PB*Tetra::N,DIM_PB*Tetra::N> Ke;
        Ke.setZero();
        std::vector<double> Le(DIM_PB*Tetra::N,0.0);
        integrales(elem,Ke);
        integrales(elem,Le);
        buildMat<Tetra::N>(elem.ind, Ke);
        buildVect<Tetra::N>(elem.ind, Le);
        } );

    /* here are the boundary conditions: either a vector s is defined on a surface, or a normal
     current density and a unit polarization vector */
    std::for_each(msh->fac.begin(),msh->fac.end(),[this](Facette::Fac &f)
        {
        std::vector<double> Le(DIM_PB*Facette::N,0.0);
        Eigen::Vector3d s_value = paramFac[f.idxPrm].s;

        if (std::isfinite(paramFac[f.idxPrm].jn) && std::isfinite(paramFac[f.idxPrm].uP.norm()))
            { s_value = -paramFac[f.idxPrm].jn*(BOHRS_MUB/CHARGE_ELECTRON)*paramFac[f.idxPrm].uP; }

        if (std::isfinite(s_value.norm()))
            {
            for (int npi=0; npi<Facette::NPI; npi++)
                {
                const double w = f.weight[npi];
                for (int ie=0; ie<Facette::N; ie++)
                    {
                    double ai_w = w*Facette::a[ie][npi];
                    Le[               ie] -= s_value[IDX_X]* ai_w;
                    Le[  Facette::N + ie] -= s_value[IDX_Y]* ai_w;
                    Le[2*Facette::N + ie] -= s_value[IDX_Z]* ai_w;
                    }
                }
            }
        buildVect<Facette::N>(f.ind, Le);
        });

    std::vector<double> Xw(DIM_PB*NOD);
    algebra::bicg_dir(iter, K, Xw, L_rhs, valDirichlet, idxDirichlet);

    for (int i=0; i<NOD; i++)
        for (int j=0; j<DIM_PB; j++)
            { s[i][j] = Xw[DIM_PB*i + j]; }
    return (iter.status == algebra::CONVERGED);
    }

void spinAcc::integrales(Tetra::Tet &tet,
                         Eigen::Matrix<double,DIM_PB*Tetra::N,DIM_PB*Tetra::N> &AE)
    {
    /* non magnetic metal contribution to AE has a block diagonal structure:
     * AE = (A 0 0)
            (0 A 0)
            (0 0 A)
    A is a N*N matrix
    A = asMatrixDiagonal(a*w)/tau_sf + c*da*transpose(da) with c = D0*sum(weight)
    */
    using namespace Tetra;

    const double lsf = getLsf(tet);
    const double D0 = getDiffusionCst(tet);

    Eigen::Matrix<double,N,1> a_w = eigen_a*tet.weight;
    Eigen::Matrix<double,N,1> diag = (D0/sq(lsf))*a_w; // units: [D0/sq(lsf)] = s^-1 : it is 1/tau_sf
    Eigen::Matrix<double,N,N> diagBlock = tet.calcDiagBlock(D0,diag);
    AE.block<N,N>(0,0) += diagBlock;
    AE.block<N,N>(N,N) += diagBlock;
    AE.block<N,N>(2*N,2*N) += diagBlock;

//here is the magnetic contribution to AE, it is also block diagonal, and antisymmetric
    if(msh->isMagnetic(tet))
        {
        const double invTau_sd = D0/sq(getLsd(tet)); //units: [D0/sq(lsd)] = s^-1 : it is 1/tau_sd
        diag = invTau_sd * a_w.cwiseProduct(tet.calcOffDiagBlock(IDX_X));
        AE.block<N,N>(N,2*N).diagonal() += diag;
        AE.block<N,N>(2*N,N).diagonal() -= diag;
        diag = invTau_sd * a_w.cwiseProduct(tet.calcOffDiagBlock(IDX_Y));
        AE.block<N,N>(0,2*N).diagonal() -= diag;
        AE.block<N,N>(2*N,0).diagonal() += diag;
        diag = invTau_sd * a_w.cwiseProduct(tet.calcOffDiagBlock(IDX_Z));
        AE.block<N,N>(0,N).diagonal() += diag;
        AE.block<N,N>(N,0).diagonal() -= diag;
        }
    }

void spinAcc::integrales(Tetra::Tet &tet, std::vector<double> &BE)
    {
    using namespace Tetra;

    /* constant cst0 in a magnetic region is the only RHS parameter involved in the diffusion
    * equation for magnetic contribution
    * units: [cst0] = [sigma] m^2 = A^2 s^3 m^-1 kg^-1
    */
    const double cst0 = BOHRS_MUB*getPolarizationRate(tet)*getSigma(tet)/CHARGE_ELECTRON;
    Eigen::Matrix<double,Nodes::DIM,NPI> &_gradV = gradV[tet.idx];

    if(msh->isMagnetic(tet))
        {
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
    if(paramTet[tet.idxPrm].spinHall != 0)
        {
        const double cst0 = getSpinHall(tet)*CHARGE_ELECTRON/MASS_ELECTRON;
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
    }

