#include "algebra/bicg.h"
#include "algebra/block_precond.h"
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
            valDirichlet.resize(DIM_PB * NOD);
            boundaryConditions();

            electrostatSolver pot_solver(_msh, mySettings.paramTetra, mySettings.paramTriangle,
                                         1e-8, mySettings.verbose, 1000);
            pot_solver.checkBoundaryConditions();
            pot_solver.V.resize(_msh.getNbNodes());
            std::string V_fileName("");
            if(mySettings.V_file)
                { V_fileName = mySettings.getSimName() + "_V.sol"; }
            pot_solver.compute(mySettings.verbose, V_fileName);

            V = pot_solver.V;
            s.resize(NOD);
            prepareExtraField();
            if(!compute())
                {
                std::cerr << "Error: spin diffusion solver(first try) failed.\n";
                exit(1);
                }
            }
        }

void spinAcc::checkBoundaryConditions(void) const
    {
    // nbVolP and nbVolN0 initialized to 1 because of __default__
    unsigned int nbVolP(1);
    unsigned int nbVolN0(1);
    for (const Tetra::prm &p : paramTet)
        {
        if(p.regName != "__default__")
            {
            if (std::isfinite(p.N0) && (p.N0 != 0))
                { nbVolN0++; }
            if (std::isfinite(p.P) && (0 <= p.P) && (p.P < 1.0))
                { nbVolP++; }
            }
        }
    int nbSurfJ(0);
    int nbSurfS(0);
    for (const Triangle::prm &p : paramTri)
        {
        if(p.regName != "__default__")
            {
            if (std::isfinite(p.jn))
                { nbSurfJ++; }
            if (std::isfinite(p.s.norm()))
                { nbSurfS++; }
            }
        }

    // one surface with jn (as the electrostatic problem), at least one surface with a given s
    bool result = ( (nbSurfJ == 1) && (nbSurfS >= 1)
           && (nbVolN0 == paramTet.size())
           && (nbVolP == paramTet.size()) );

    if (!result)
        {
        std::cerr << "Error: incorrect boundary conditions for spin diffusion solver.\n";
        exit(1);
        }
    else if (verbose)
        { std::cout << "spin diffusion problem boundary conditions Ok.\n"; }
    }

void spinAcc::fillDirichletData(const int k, Eigen::Vector3d &s_value)
    {
    valDirichlet[DIM_PB*k  ] = s_value[Nodes::IDX_X];
    idxDirichlet.push_back(DIM_PB * k    );
    valDirichlet[DIM_PB*k+1] = s_value[Nodes::IDX_Y];
    idxDirichlet.push_back(DIM_PB * k + 1);
    valDirichlet[DIM_PB*k+2] = s_value[Nodes::IDX_Z];
    idxDirichlet.push_back(DIM_PB * k + 2);
    }

void spinAcc::boundaryConditions(void)
    {
    /* Dirichlet conditions on all the surfaces where s is given, including the surface where the
     * current density jn is injected (as in st-feeLLGood). Where jn and uP are given and s is not,
     * the injection is a flux (Neumann) condition, see solve() */
    std::fill(valDirichlet.begin(), valDirichlet.end(), 0.0);
    for (const Triangle::Tri &f : msh->tri)
        {
        if (std::isfinite(paramTri[f.idxPrm].s.norm()))
            {
            Eigen::Vector3d s_value = paramTri[f.idxPrm].s;
            for(int j = 0; j < Triangle::N; j++)
                { fillDirichletData(f.ind[j], s_value); }
            }
        }
    suppress_copies<int>(idxDirichlet);
    }

double spinAcc::getMs(const Tetra::Tet &tet) const
    { return paramTet[tet.idxPrm].Ms; }

double spinAcc::getSigma(const Tetra::Tet &tet) const
    { return paramTet[tet.idxPrm].sigma; }

double spinAcc::getDiffusionCst(const Tetra::Tet &tet) const
    {
    const double N0 = paramTet[tet.idxPrm].N0;
    return 2.0 * getSigma(tet) / (sq(CHARGE_ELECTRON) * N0);
    }

double spinAcc::getPolarizationRate(const Tetra::Tet &tet) const
    { return paramTet[tet.idxPrm].P; }

double spinAcc::getLsd(const Tetra::Tet &tet) const
    { return paramTet[tet.idxPrm].lsd; }

double spinAcc::getLsf(const Tetra::Tet &tet) const
    { return paramTet[tet.idxPrm].lsf; }

Eigen::Matrix<double,Nodes::DIM,Tetra::N> spinAcc::calc_u_nod(const Tetra::Tet &tet) const
    {
    Eigen::Matrix<double,Nodes::DIM,Tetra::N> u_nod;
    for (int ie = 0; ie < Tetra::N; ie++)
        { u_nod.col(ie) = msh->getNode_u(tet.ind[ie]); }
    return u_nod;
    }

void spinAcc::prepareExtraField(void) const
    {
    for (const int idxTet : msh->magTet)
        {
        Tetra::Tet &t = msh->tet[idxTet];
        double D0 = getDiffusionCst(t);
        double prefactor = D0 / (sq(getLsd(t)) * gamma0 * getMs(t));
        t.extraField = [this, &t, prefactor](Eigen::Ref<Eigen::Matrix<double, Nodes::DIM, Tetra::NPI>> H)
                        { H += calc_Hst(t, prefactor, s); };
        }
    }

bool spinAcc::compute(void)
    {
    //prepareExtraField(); // do we have to do that once or each time we want another s computation ?
    bool has_converged = solve();
    if (!has_converged)
        {
        std::cout << "spin accumulation solver: " << iter.infos() << "\n";
        for(int i = 0; i < NOD; i++)
            { s[i].setZero(); }
        }
    else if (verbose)
        { std::cout << "spin accumulation solved.\n"; }
    return has_converged;
    }

bool spinAcc::solve(void)
    {
    iter.reset();

    /* the system depends on the magnetization: it is assembled again at each call, from zero
     * (otherwise the contributions of all the previous calls would be summed up) */
    K.clear();
    std::fill(L_rhs.begin(), L_rhs.end(), 0.0);

    for (Tetra::Tet &elem : msh->tet)
        {
        Eigen::Matrix<double,DIM_PB*Tetra::N,DIM_PB*Tetra::N> Ke;
        Ke.setZero();
        std::vector<double> Le(DIM_PB * Tetra::N, 0.0);
        integrales(elem, Ke);
        integrales(elem, Le);
        buildMat<Tetra::N>(elem.ind, Ke);
        buildVect<Tetra::N>(elem.ind, Le);
        }

    /* flux (Neumann) boundary condition on the surface where the current density jn is injected
     * with the polarization uP: Q_n = -(mu_B/e) jn uP, contribution -int_Gamma Q_n a_i to the RHS
     * (the surfaces where s is given are Dirichlet conditions, see boundaryConditions()) */
    for (Triangle::Tri &f : msh->tri)
        {
        if (std::isfinite(paramTri[f.idxPrm].jn) && std::isfinite(paramTri[f.idxPrm].uP.norm()))
            {
            std::vector<double> Le(DIM_PB * Triangle::N, 0.0);
            const Eigen::Vector3d Qn =
                    -paramTri[f.idxPrm].jn * (BOHRS_MUB/CHARGE_ELECTRON) * paramTri[f.idxPrm].uP;
            for (int npi = 0; npi < Triangle::NPI; npi++)
                {
                const double w = f.weight[npi];
                for (int ie = 0; ie < Triangle::N; ie++)
                    {
                    double ai_w = w * Triangle::a[ie][npi];
                    Le[              ie] -= Qn[IDX_X] * ai_w;
                    Le[  Triangle::N+ie] -= Qn[IDX_Y] * ai_w;
                    Le[2*Triangle::N+ie] -= Qn[IDX_Z] * ai_w;
                    }
                }
            buildVect<Triangle::N>(f.ind, Le);
            }
        }

    std::vector<double> Xw(DIM_PB * NOD);
    /* block diagonal preconditioner (inverse of the 3x3 block of each node), as st-feeLLGood: it
     * captures the precession coupling s x u, local to each node, that the diagonal one misses */
    const algebra::BlockDiagPrecond<DIM_PB> precond(K, NOD, idxDirichlet);
    algebra::bicg_dir_prec(iter, K, Xw, L_rhs, valDirichlet, idxDirichlet,
                           [&precond](const std::vector<double> &p, std::vector<double> &phat)
                           { precond.apply(p, phat); });

    for (int i = 0; i < NOD; i++)
        {
        for (int j = 0; j < DIM_PB; j++)
            { s[i][j] = Xw[DIM_PB*i+j]; }
        }
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

    Eigen::Matrix<double,N,1> a_w = eigen_a * tet.weight;
    Eigen::Matrix<double,N,1> diag = (D0 / sq(lsf)) * a_w; // units: [D0/sq(lsf)] = s^-1 :
                                                       // it is 1/tau_sf
    Eigen::Matrix<double,N,N> diagBlock = tet.calcDiagBlock(D0, diag);
    AE.block<N,N>(    0,     0) += diagBlock;
    AE.block<N,N>(    N,     N) += diagBlock;
    AE.block<N,N>(2 * N, 2 * N) += diagBlock;

//here is the magnetic contribution to AE, it is also block diagonal, and antisymmetric
    if(msh->isMagnetic(tet))
        {
        const double invTau_sd = D0 / sq(getLsd(tet)); //units: [D0/sq(lsd)] = s^-1 : it is 1/tau_sd
        /* nodal magnetization: the same (NEXT, the latest one) as in the RHS */
        const Eigen::Matrix<double,Nodes::DIM,N> u_nod = calc_u_nod(tet);

        diag = invTau_sd * a_w.cwiseProduct(u_nod.row(IDX_X).transpose());
        AE.block<N,N>(    N, 2 * N).diagonal() += diag;
        AE.block<N,N>(2 * N,     N).diagonal() -= diag;
        diag = invTau_sd * a_w.cwiseProduct(u_nod.row(IDX_Y).transpose());
        AE.block<N,N>(    0, 2 * N).diagonal() -= diag;
        AE.block<N,N>(2 * N,     0).diagonal() += diag;
        diag = invTau_sd * a_w.cwiseProduct(u_nod.row(IDX_Z).transpose());
        AE.block<N,N>(    0,     N).diagonal() += diag;
        AE.block<N,N>(    N,     0).diagonal() -= diag;
        }
    }

void spinAcc::integrales(Tetra::Tet &tet, std::vector<double> &BE)
    {
    using namespace Tetra;

    /* constant cst0 in a magnetic region is the only RHS parameter involved in the diffusion
    * equation for magnetic contribution
    * units: [cst0] = [sigma] m^2 = A^2 s^3 m^-1 kg^-1
    */
    const double cst0 = BOHRS_MUB * getPolarizationRate(tet) * getSigma(tet) / CHARGE_ELECTRON;
    Eigen::Matrix<double, Nodes::DIM, NPI> gradV = tet.gradV(V);

    if(msh->isMagnetic(tet))
        {
        /* magnetization at the Gauss points: u varies in the element, a nodal u_i would not be a
         * lumping */
        const Eigen::Matrix<double,Nodes::DIM,NPI> U = calc_u_nod(tet) * eigen_a;
        for (size_t npi = 0; npi < NPI; npi++)
            {
            const Eigen::Vector3d cst0_w_gradV = cst0 * tet.weight[npi] * gradV.col(npi);

            for (size_t ie = 0; ie < N; ie++)
                {
                const double tmp = cst0_w_gradV.dot(tet.da.row(ie));
                BE[    ie] += tmp * U(0,npi);
                BE[  N+ie] += tmp * U(1,npi);
                BE[2*N+ie] += tmp * U(2,npi);
                }
            }
        }
    }

