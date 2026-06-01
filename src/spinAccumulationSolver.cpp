#include "algebra/bicg.h"
#include "spinAccumulationSolver.h"
#include "chronometer.h" //date()

using algebra::sq;
using namespace Nodes;

    spinAcc::spinAcc(const Settings &mySettings /**< [in] */,
            Mesh::mesh &_msh /**< [in] ref to the mesh */,
            const double _tol /**< [in] tolerance for bicg_dir solver */,
            const int max_iter /**< [in] maximum number of iterations */):
            solver<DIM_PB_SPIN_ACC>(_msh, mySettings.paramTetra, mySettings.paramTriangle,
                                    "bicg_dir", _tol, mySettings.verbose, max_iter)
        {
        if (mySettings.spin_acc)
            {
            checkBoundaryConditions();
            valDirichlet.resize(DIM_PB*NOD);
            boundaryConditions();
            electrostatSolver pot_solver(_msh, mySettings.paramTetra, mySettings.paramTriangle,
                                         1e-8, mySettings.verbose, 1000);
            pot_solver.checkBoundaryConditions();
            pot_solver.V.resize(_msh.getNbNodes());
            std::string V_fileName("");
            if(mySettings.V_file)
                V_fileName = mySettings.getSimName() + "_V.sol";
            pot_solver.compute(mySettings.verbose, V_fileName);
            V = pot_solver.V;
            s.resize(NOD);
            if(!compute())
                {
                std::cout << "Error: spin diffusion solver(first try) failed.\n";
                exit(1);
                }
            }
        }

void spinAcc::checkBoundaryConditions(void) const
    {
    // nbVolP and nbVolN0 initialized to 1 because of __default__
    unsigned int nbVolP(1);
    unsigned int nbVolN0(1);
    std::for_each(paramTet.begin(),paramTet.end(),[&nbVolN0,&nbVolP](const Tetra::prm &p)
        {
        if(p.regName != "__default__")
            {
            if (std::isfinite(p.N0) && (p.N0 != 0)) nbVolN0++;
            if(std::isfinite(p.P) && (0 <= p.P) && (p.P < 1.0)) nbVolP++;
            }
        });
    int nbSurfJ(0);
    int nbSurfS(0);
    std::for_each(paramTri.begin(),paramTri.end(),[&nbSurfJ,&nbSurfS](const Triangle::prm &p)
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
    std::for_each(msh->tri.begin(),msh->tri.end(),[this](const Triangle::Tri &f)
        {
        if (std::isnan(paramTri[f.idxPrm].s.norm()) &&  std::isfinite(paramTri[f.idxPrm].jn))
            {
             /* units:
             * [jn] = A m^-2; [BOHRS_MUB/CHARGE_ELECTRON] = m^2 ; [P] = 1 ⇒ [s_value] = A
             * check s_value formula, especially the sign with current convention
             * */
            Eigen::Vector3d s_value =
                    -paramTri[f.idxPrm].jn*(BOHRS_MUB/CHARGE_ELECTRON)*paramTri[f.idxPrm].uP;
            for(int j=0;j<Triangle::N;j++)
                { fillDirichletData(f.ind[j],s_value); }
            }
        else if (std::isfinite(paramTri[f.idxPrm].s.norm()) &&  std::isnan(paramTri[f.idxPrm].jn))
            {
            Eigen::Vector3d s_value = paramTri[f.idxPrm].s;
            for(int j=0;j<Triangle::N;j++)
                { fillDirichletData(f.ind[j],s_value); }
            }
        });
    suppress_copies<int>(idxDirichlet);
    }

double spinAcc::getMs(const Tetra::Tet &tet) const
    { return paramTet[tet.idxPrm].Ms; }

double spinAcc::getSigma(const Tetra::Tet &tet) const
    { return paramTet[tet.idxPrm].sigma; }

double spinAcc::getDiffusionCst(const Tetra::Tet &tet) const
    {
    const double N0 = paramTet[tet.idxPrm].N0;
    return 2.0*getSigma(tet)/(sq(CHARGE_ELECTRON)*N0);
    }

double spinAcc::getPolarizationRate(const Tetra::Tet &tet) const
    { return paramTet[tet.idxPrm].P; }

double spinAcc::getLsd(const Tetra::Tet &tet) const
    { return paramTet[tet.idxPrm].lsd; }

double spinAcc::getLsf(const Tetra::Tet &tet) const
    { return paramTet[tet.idxPrm].lsf; }

void spinAcc::prepareExtraField(void)
    {
    using namespace Tetra;

    std::for_each( msh->tet.begin(), msh->tet.end(), [this](Tet &t)
        {
// this is wierd, since only mag tetrahedrons seems concerned, we should loop over magTet, not tet
        if (msh->isMagnetic(t))
            {
            double D0 = getDiffusionCst(t);
            double prefactor = D0/(sq(getLsd(t))*gamma0*getMs(t));
            t.extraField = [this, &t, prefactor](Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> H)
                        { H += calc_Hst(t, prefactor, s); };
            }
        });//end for_each
    }

bool spinAcc::compute(void)
    {
    prepareExtraField(); // do we have to do that once or each time we want another s computation ?
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
    std::for_each(msh->tri.begin(),msh->tri.end(),[this](Triangle::Tri &f)
        {
        std::vector<double> Le(DIM_PB*Triangle::N,0.0);
        Eigen::Vector3d s_value = paramTri[f.idxPrm].s;

        if (std::isfinite(paramTri[f.idxPrm].jn) && std::isfinite(paramTri[f.idxPrm].uP.norm()))
            { s_value = -paramTri[f.idxPrm].jn*(BOHRS_MUB/CHARGE_ELECTRON)*paramTri[f.idxPrm].uP; }

        if (std::isfinite(s_value.norm()))
            {
            for (int npi=0; npi<Triangle::NPI; npi++)
                {
                const double w = f.weight[npi];
                for (int ie=0; ie<Triangle::N; ie++)
                    {
                    double ai_w = w*Triangle::a[ie][npi];
                    Le[               ie] -= s_value[IDX_X]* ai_w;
                    Le[  Triangle::N + ie] -= s_value[IDX_Y]* ai_w;
                    Le[2*Triangle::N + ie] -= s_value[IDX_Z]* ai_w;
                    }
                }
            }
        buildVect<Triangle::N>(f.ind, Le);
        });

    std::vector<double> Xw(DIM_PB*NOD);
    algebra::bicg_dir(iter, K, Xw, L_rhs, valDirichlet, idxDirichlet);

    for (int i=0; i<NOD; i++)
        for (int j=0; j<DIM_PB; j++)
            { s[i][j] = Xw[DIM_PB*i + j]; }
    return (iter.status == algebra::CONVERGED);
    }

void spinAcc::integrales(const Tetra::Tet &tet,
                         Eigen::Matrix<double,DIM_PB*Tetra::N,DIM_PB*Tetra::N> &AE) const
    {
    /* non-magnetic metal contribution to AE has a block diagonal structure:
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
    Eigen::Matrix<double,N,1> diag = (D0/sq(lsf))*a_w; // units: [D0/sq(lsf)] = s^-1 :
                                                       // it is 1/tau_sf
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
    const Eigen::Matrix<double,Nodes::DIM,NPI> &_gradV = calc_gradV(tet,V);

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
    }

