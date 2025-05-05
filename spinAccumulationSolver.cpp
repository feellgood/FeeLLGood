#include "algebra/bicg.h"
#include "spinAccumulationSolver.h"

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

Eigen::Vector3d spinAcc::get_Qn(Facette::Fac const &fac) const
    {
    Eigen::Vector3d Pu = paramFacette[fac.idxPrm].Pu;
    double J  = paramFacette[fac.idxPrm].J;
    return J*(BOHRS_MUB/CHARGE_ELECTRON)*Pu;
    }

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
        const double Js = getJs(tet);// Js = Ms/nu_0
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

int spinAcc::solve(const double _tol /**< [in] tolerance */ )
    {
    const int DIM_3D = 3;
    algebra::w_sparseMat Kw(DIM_3D*NOD);
    std::vector<double> Lw(DIM_3D*NOD);

    std::for_each( msh.tet.begin(), msh.tet.end(),[this,&Kw,&Lw](Tetra::Tet &elem)
        {
        Eigen::Matrix<double,DIM_3D*Tetra::N,DIM_3D*Tetra::N> K;
        K.setZero();
        Eigen::Matrix<double,DIM_3D*Tetra::N,1> L;
        L.setZero();
        integrales(elem,K,L);
        assembling<Tetra::N>(elem.ind,K,L,Kw,Lw);
        } );

    std::for_each( msh.fac.begin(), msh.fac.end(),[this,&Kw,&Lw](Facette::Fac &elem)
        {
        Eigen::Matrix<double,DIM_3D*Facette::N,DIM_3D*Facette::N> K;
        K.setZero();
        Eigen::Matrix<double,DIM_3D*Facette::N,1> L;
        L.setZero();
        integrales(elem,L);
        assembling<Facette::N>(elem.ind,K,L,Kw,Lw);
        } );

    algebra::iteration iter("bicg",_tol,verbose,MAXITER);

    std::cout << "bicg...\n";
    algebra::r_sparseMat Kr(Kw);
    std::vector<double> Xw(DIM_3D*NOD);
    algebra::bicg(iter, Kr, Xw, Lw);

    if (!(iter.converged()))
        {
        std::cout << "\t bicg FAILED in " << iter.get_iteration() << "; residu= " << iter.get_res() << '\t';
        }

    for (int i=0; i<NOD; i++)
        {
        for (int d=0; d<DIM_3D; d++)
            { Qs[i][d] = Xw[d*NOD+i]; }
        }

    return iter.converged();
    }

void spinAcc::integrales(Tetra::Tet &tet,
                         Eigen::Matrix<double,DIM_PROBLEM*Tetra::N,DIM_PROBLEM*Tetra::N> &AE,
                         Eigen::Matrix<double,DIM_PROBLEM*Tetra::N,1> &BE)
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
            double Dai_DV = grad_ai.dot( gradV.col(npi) );
            double tmp = BOHRS_MUB*beta*sigma/CHARGE_ELECTRON* Dai_DV *w;
            BE[    ie] += tmp*u_nod[0][ie]; // A. de Riz 2019 moins -> plus : to check : c'est quoi ce bazar de signes qui changent ???
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
            AE(    ie,     ie) += tmp;  // lumping
            AE(  N+ie,   N+ie) += tmp;
            AE(2*N+ie, 2*N+ie) += tmp;

            tmp = Tetra::a[ie][npi]*D0*w/sq(EPSILON + lsd); // lsd is lJ  //D0*pow(ilJ, 2.)*ai*w; with ilJ = 1/lJ;
            AE(    ie,   N+ie) += tmp*u_nod[2][ie]; // lumping
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

void spinAcc::integrales(Facette::Fac &fac, Eigen::Matrix<double,DIM_PROBLEM*Facette::N,1> &BE)
    {
    using namespace Nodes;
    Eigen::Vector3d Qn = get_Qn(fac);
    for (int npi=0; npi<Facette::NPI; npi++)
        {
        double w = fac.weight[npi];
        for (int ie=0; ie<Facette::N; ie++)
            {
            double ai_w = w*Facette::a[ie][npi];
            BE[      ie]        -= Qn[IDX_X]* ai_w;
            BE[  Facette::N+ie] -= Qn[IDX_Y]* ai_w;
            BE[2*Facette::N+ie] -= Qn[IDX_Z]* ai_w;
            }
        }
    }
