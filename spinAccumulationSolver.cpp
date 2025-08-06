#include "algebra/bicg.h"
#include "spinAccumulationSolver.h"
#include "chronometer.h" //date()
#include "tags.h"
#include "solverUtils.h"

double spinAcc::getJs(Tetra::Tet const &tet) const
    { return paramTetra[tet.idxPrm].J; }

double spinAcc::getSigma(Tetra::Tet const &tet) const
    { return paramTetra[tet.idxPrm].sigma; }

double spinAcc::getN0(Tetra::Tet &tet) const
    { return paramTetra[tet.idxPrm].N0; }

double spinAcc::getBeta(Tetra::Tet &tet) const
    { return paramTetra[tet.idxPrm].beta; }

double spinAcc::getLsd(Tetra::Tet &tet) const
    { return paramTetra[tet.idxPrm].lsd; }

double spinAcc::getLsf(Tetra::Tet &tet) const
    { return paramTetra[tet.idxPrm].lsf; }

void spinAcc::prepareExtras(void)
    {
    using namespace Nodes;
    using namespace Tetra;

    std::for_each( msh.tet.begin(), msh.tet.end(), [this](Tet &t)
        {
        const int _idx = t.idx;
        const double sigma = elec.getSigma(t);
        const double beta = getBeta(t);
        const double lsd = getLsd(t);
        const double lsf = getLsf(t);
        const double ksi = sq(lsd/lsf);
        const double Js = getJs(t);// Js = Ms/nu_0
        double prefactor = mu0*BOHRS_MUB*beta/(gamma0*Js*CHARGE_ELECTRON*(1.0 + sq(ksi)));

        t.extraField = [this, _idx](Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> H)
                         { for(int npi = 0; npi<Tetra::NPI; npi++) { H.col(npi) += elec.Hm[_idx].col(npi); } };

        t.extraCoeffs_BE = [this, sigma, ksi, prefactor, &t](Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> U,
                                          Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> dUdx,
                                          Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> dUdy,
                                          Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> dUdz,
                                          Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,N>> BE)
                {
                for (int npi = 0; npi < NPI; npi++)
                    {
                    Eigen::Vector3d const &_gV = elec.gradV[t.idx].col(npi);
                    Eigen::Vector3d j_grad_u =
                    -sigma * Eigen::Vector3d(_gV.dot( Eigen::Vector3d(dUdx(IDX_X,npi), dUdy(IDX_X,npi), dUdz(IDX_X,npi))),
                                             _gV.dot( Eigen::Vector3d(dUdx(IDX_Y,npi), dUdy(IDX_Y,npi), dUdz(IDX_Y,npi))),
                                             _gV.dot( Eigen::Vector3d(dUdx(IDX_Z,npi), dUdy(IDX_Z,npi), dUdz(IDX_Z,npi))));

                    Eigen::Vector3d m = ksi * j_grad_u + U.col(npi).cross(j_grad_u);
                    for (int i = 0; i < N; i++)
                        {
                        const double ai_w = t.weight[npi] * a[i][npi];
                        BE.col(i) += ai_w*( elec.Hm[t.idx].col(npi) + prefactor*m);
                        }
                    } // end loop on npi
                }; //end lambda
        });//end for_each
    }

int spinAcc::solve(void)
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

    std::for_each( msh.fac.begin(), msh.fac.end(),[this,&Lw](Facette::Fac &elem)
        {
        std::vector<double> L(DIM_PROBLEM*Facette::N,0.0);
        integrales(elem,L);
        buildVect<Facette::N,DIM_PROBLEM>(elem.ind, L, Lw);
        } );

    algebra::r_sparseMat Kr(Kw);
    std::vector<double> Xw(DIM_PROBLEM*NOD);
    algebra::bicg(iter, Kr, Xw, Lw);

    if ( (iter.status  == algebra::ITER_OVERFLOW) || (iter.status == algebra::CANNOT_CONVERGE) || (iter.get_res() > iter.resmax) )
        {
        if (verbose)
            { std::cout << "spin accumulation solver: " << iter.infos() << std::endl; }
        return 1;
        }

    for (int i=0; i<NOD; i++)
        {
        Qs[i][IDX_X] = Xw[DIM_PROBLEM*i];
        Qs[i][IDX_Y] = Xw[DIM_PROBLEM*i + 1];
        Qs[i][IDX_Z] = Xw[DIM_PROBLEM*i + 2];
        }
    return 0;
    }

void spinAcc::integrales(Tetra::Tet &tet,
                         Eigen::Matrix<double,DIM_PROBLEM*Tetra::N,DIM_PROBLEM*Tetra::N> &AE,
                         std::vector<double> &BE)
    {
    using algebra::sq;
    using namespace Nodes;
    using namespace Tetra;

    double Js = getJs(tet);
    double sigma = getSigma(tet);
    //double CSH = msh.find_param("CSH",reg); // to take into account Spin Hall effects
    double N0 = getN0(tet);
    double beta = getBeta(tet);
    double lsd = getLsd(tet);
    double lsf = getLsf(tet);
    double u_nod[DIM][N] {{0}};
    Eigen::Matrix<double,N,1> V_nod;
    Eigen::Matrix<double,Nodes::DIM,NPI> gradV;

    for (size_t ie=0; ie<N; ie++)
        {
        size_t i= tet.ind[ie];
        if (Js>0.0) // if Js <=0 then it is a non magnetic material, u_nod is zero
            for (size_t d=0; d<DIM; d++)
                {
                u_nod[d][ie] = msh.getNode_u(i)(d);//msh.node[i].u[d]; // carefull here do we need u comp of Node(NEXT) or Node(CURRENT) ?
                }
        V_nod[ie] = elec.V[i];
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

    double D0=2.0*sigma/(sq(CHARGE_ELECTRON)*N0);
    for (size_t npi=0; npi<NPI; npi++)
        {
        double w = tet.weight[npi];

        for (size_t ie=0; ie<N; ie++)
            {
            Eigen::Vector3d grad_ai = tet.da.row(ie);
            double tmp = BOHRS_MUB*beta*sigma/CHARGE_ELECTRON*w*grad_ai.dot( gradV.col(npi) );
            BE[    ie] += tmp*u_nod[0][ie];
            BE[  N+ie] += tmp*u_nod[1][ie];
            BE[2*N+ie] += tmp*u_nod[2][ie];
/*          if (CSH != 0)  // SOT : spin orbit torque CSH is constant Spin Hall
                {
                tmp = CSH*CHARGE_ELECTRON/MASS_ELECTRON *w; // obviously SOT contribution to BE is a cross product, should be rewritten
                BE[    ie] -= tmp*(dai_dz*dVdy[npi]-dai_dy*dVdz[npi]);
                BE[  N+ie] -= tmp*(dai_dx*dVdz[npi]-dai_dz*dVdx[npi]);
                BE[2*N+ie] -= tmp*(dai_dy*dVdx[npi]-dai_dx*dVdy[npi]);
                } */
            tmp = Tetra::a[ie][npi]*D0*w/sq(lsf + EPSILON);
            AE(    ie,     ie) += tmp;
            AE(  N+ie,   N+ie) += tmp;
            AE(2*N+ie, 2*N+ie) += tmp;

            tmp = Tetra::a[ie][npi]*D0*w/sq(EPSILON + lsd);
            AE(    ie,   N+ie) += tmp*u_nod[2][ie];
            AE(    ie, 2*N+ie) -= tmp*u_nod[1][ie];
            AE(  N+ie,     ie) -= tmp*u_nod[2][ie];
            AE(  N+ie, 2*N+ie) += tmp*u_nod[0][ie];
            AE(2*N+ie,     ie) += tmp*u_nod[1][ie];
            AE(2*N+ie,   N+ie) -= tmp*u_nod[0][ie];

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

void spinAcc::integrales(Facette::Fac &fac, std::vector<double> &BE)
    {
    using namespace Nodes;
    using namespace Facette;
// J in the following equation is a normal current density to the facette 
    Eigen::Vector3d Qn = -paramFacette[fac.idxPrm].J*(BOHRS_MUB/CHARGE_ELECTRON)*paramFacette[fac.idxPrm].Pu;
//if(!std::isfinite(paramFacette[fac.idxPrm].J)) {std::cout << "ouch J not finite\n";exit(1); }

    for (int npi=0; npi<NPI; npi++)
        {
        double w = fac.weight[npi];
        for (int ie=0; ie<N; ie++)
            {
            double ai_w = w*a[ie][npi];
            BE[    ie] -= Qn[IDX_X]* ai_w;
            BE[  N+ie] -= Qn[IDX_Y]* ai_w;
            BE[2*N+ie] -= Qn[IDX_Z]* ai_w;
            }
        }
    }

bool spinAcc::save(std::string const &metadata) const
    {
    using namespace Nodes;
    std::ofstream fout(Q_fileName, std::ios::out);
    if (fout.fail())
        {
        std::cout << "cannot open file " << Q_fileName << std::endl;
        SYSTEM_ERROR;
        }

    fout << tags::sol::rw_time << ' ' << date() << '\n' << metadata << std::scientific
         << std::setprecision(precision);

    for (int i = 0; i < NOD; i++)
        {
        int _i = msh.getNodeIndex(i);
        fout << i << '\t' << Qs[_i][IDX_X] << '\t' << Qs[_i][IDX_Y] << '\t' << Qs[_i][IDX_Z] << std::endl;
        }
    fout.close();
    return !(fout.good());
    }
