#include "algebra/bicg.h"
#include "spinAccumulationSolver.h"
#include "chronometer.h" //date()
#include "tags.h"
#include "solverUtils.h"

bool spinAcc::checkBoundaryConditions(void) const
    {
    int nbSurfJ(0);
    int nbSurfZeroS(0);
    std::for_each(paramFacette.begin(),paramFacette.end(),[&nbSurfJ,&nbSurfZeroS](Facette::prm const &p)
        {
        if (std::isfinite(p.J) && (p.J != 0)) nbSurfJ++;
        if (p.s.norm() == 0) nbSurfZeroS++;
        });
    return ((nbSurfJ == 1)&&(nbSurfZeroS == 1));
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
    using namespace Nodes;
    using namespace Tetra;

    std::for_each( msh.tet.begin(), msh.tet.end(), [this](Tet &t)
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
    msh.s.resize(NOD);
    std::for_each(msh.tet.begin(), msh.tet.end(), [this](Tetra::Tet const &tet)
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
            { msh.s[i].setZero(); }
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
    using namespace Nodes;

    iter.reset();
    algebra::w_sparseMat Kw(DIM_PROBLEM*NOD);
    std::vector<double> Lw(DIM_PROBLEM*NOD,0);

    std::for_each( msh.tet.begin(), msh.tet.end(),[this,&Kw,&Lw](Tetra::Tet &elem)
        {
        Eigen::Matrix<double,DIM_PROBLEM*Tetra::N,DIM_PROBLEM*Tetra::N> K;
        K.setZero();
        std::vector<double> L(DIM_PROBLEM*Tetra::N,0.0);
        integrales(elem,K,L);
        buildMat<Tetra::N,DIM_PROBLEM>(elem.ind, K, Kw);
        buildVect<Tetra::N,DIM_PROBLEM>(elem.ind, L, Lw);
        } );

    // here are the boundary conditions
    std::for_each(msh.fac.begin(),msh.fac.end(),[this,&Lw](Facette::Fac &f)
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
    algebra::bicg(iter, Kr, Xw, Lw);

    for (int i=0; i<NOD; i++)
        {
        msh.s[i] = Eigen::Vector3d(Xw[DIM_PROBLEM*i], Xw[DIM_PROBLEM*i + 1], Xw[DIM_PROBLEM*i +2]);
        }
    return (iter.status == algebra::CONVERGED);
    }

void spinAcc::integrales(Tetra::Tet &tet,
                         Eigen::Matrix<double,DIM_PROBLEM*Tetra::N,DIM_PROBLEM*Tetra::N> &AE,
                         std::vector<double> &BE)
    {
    using algebra::sq;
    using namespace Nodes;
    using namespace Tetra;

    double sigma = getSigma(tet);
    double spinHall = getSpinHall(tet);
    double N0 = getN0(tet);
    double P = getPolarization(tet);
    double lsf = getLsf(tet);
    double D0=2.0*sigma/(sq(CHARGE_ELECTRON)*N0);
    Eigen::Matrix<double,N,1> V_nod;
    Eigen::Matrix<double,Nodes::DIM,NPI> gradV;

    for (size_t ie=0; ie<N; ie++)
        {
        size_t i= tet.ind[ie];
        V_nod[ie] = V[i];
        }

    Eigen::Matrix<double,N,NPI> dadx;
    dadx.colwise() = tet.da.col(IDX_X); // colwise() means da.col(IDX_X) is repeated to build dadx 
    Eigen::Matrix<double,N,NPI> dady;
    dady.colwise() = tet.da.col(IDX_Y);
    Eigen::Matrix<double,N,NPI> dadz;
    dadz.colwise() = tet.da.col(IDX_Z);
    // building explicitely dad(x|y|z) migth be avoided rewritting the following multiplications
    gradV.row(IDX_X) = V_nod.transpose() * dadx;// V_nod^T * dadx
    gradV.row(IDX_Y) = V_nod.transpose() * dady;// V_nod^T * dady
    gradV.row(IDX_Z) = V_nod.transpose() * dadz;// V_nod^T * dadz

    if(msh.isMagnetic(tet))
        {
        double lsd = getLsd(tet);
        double u_nod[DIM][N];

        for (size_t ie=0; ie<N; ie++)
            {
            size_t i= tet.ind[ie];
            for (size_t d=0; d<DIM; d++)
                {
                u_nod[d][ie] = msh.getNode_u(i)(d);//msh.node[i].u[d]; // carefull here do we need u comp of Node(NEXT) or Node(CURRENT) ?
                }
            }

        for (size_t npi=0; npi<NPI; npi++)
            {
            double w = tet.weight[npi];

            for (size_t ie=0; ie<N; ie++)
                {
                Eigen::Vector3d grad_ai = tet.da.row(ie);
                double tmp = BOHRS_MUB*P*sigma/CHARGE_ELECTRON*w*grad_ai.dot( gradV.col(npi) );
                BE[    ie] += tmp*u_nod[0][ie];
                BE[  N+ie] += tmp*u_nod[1][ie];
                BE[2*N+ie] += tmp*u_nod[2][ie];

                tmp = Tetra::a[ie][npi]*D0*w/sq(lsd);
                AE(    ie,   N+ie) += tmp*u_nod[2][ie];
                AE(    ie, 2*N+ie) -= tmp*u_nod[1][ie];
                AE(  N+ie,     ie) -= tmp*u_nod[2][ie];
                AE(  N+ie, 2*N+ie) += tmp*u_nod[0][ie];
                AE(2*N+ie,     ie) += tmp*u_nod[1][ie];
                AE(2*N+ie,   N+ie) -= tmp*u_nod[0][ie];
                }
            }
        }

    for (size_t npi=0; npi<NPI; npi++)
        {
        double w = tet.weight[npi];
        double tmp2 = spinHall*CHARGE_ELECTRON/MASS_ELECTRON *w;

        for (size_t ie=0; ie<N; ie++)
            {
            Eigen::Vector3d grad_ai = tet.da.row(ie);
            Eigen::Vector3d v = grad_ai.cross(gradV.col(npi));
            BE[    ie] += tmp2*v[IDX_X];
            BE[  N+ie] += tmp2*v[IDX_Y];
            BE[2*N+ie] += tmp2*v[IDX_Z];

            double tmp = Tetra::a[ie][npi]*D0*w/sq(lsf);
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

