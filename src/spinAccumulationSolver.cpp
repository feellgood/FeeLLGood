#include <algorithm>
#include <array>
#include <map>
#include <set>
#include <utility>

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
            if (mySettings.usesCurrentDensity() && !mySettings.hasCurrentDensity())
                {
                std::cerr << "Error: a surface region has jn: current_density, but the function "
                             "current_density is not defined.\n";
                exit(1);
                }
            if (mySettings.hasCurrentDensity() && !mySettings.usesCurrentDensity())
                {
                std::cerr << "Warning: current_density is defined, but no surface region has "
                             "jn: current_density: it only acts through Bext.\n";
                }
            checkBoundaryConditions();
            checkDensityOfStates();
            makeStacks();
            valDirichlet.resize(DIM_PB * NOD);
            boundaryConditions();

            electrostatSolver pot_solver(_msh, mySettings.paramTetra, mySettings.paramTriangle,
                                         1e-8, mySettings.verbose, 1000);
            pot_solver.checkBoundaryConditions();
            for (const Triangle::prm &p : mySettings.paramTriangle)
                {
                if (std::isfinite(p.V))
                    { V_D = p.V; }
                }
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

    // at least one surface with jn (as the electrostatic problem), at least one with a given s
    bool result = ( (nbSurfJ >= 1) && (nbSurfS >= 1)
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

void spinAcc::checkDensityOfStates(void) const
    {
    if (!verbose)
        { return; }
    // volume region of the first tetrahedron found with each face, by its sorted node indices
    std::map<std::array<int,Triangle::N>, int> faceRegion;
    std::set<std::pair<int,int>> pairs;  // neighbouring regions of different dens_state
    for (const Tetra::Tet &tet : msh->tet)
        {
        for (int omit = 0; omit < Tetra::N; omit++)  // the face without the node omit
            {
            std::array<int,Triangle::N> key;
            for (int ie = 0, j = 0; ie < Tetra::N; ie++)
                {
                if (ie != omit)
                    { key[j++] = tet.ind[ie]; }
                }
            std::sort(key.begin(), key.end());
            const auto [it, inserted] = faceRegion.emplace(key, tet.idxPrm);
            if (!inserted && it->second != tet.idxPrm
                && paramTet[it->second].N0 != paramTet[tet.idxPrm].N0)
                { pairs.insert(std::minmax(it->second, tet.idxPrm)); }
            }
        }
    for (const auto &[a, b] : pairs)
        {
        std::cout << "Warning: the volume regions " << paramTet[a].regName << " and "
                  << paramTet[b].regName << " have different dens_state (" << paramTet[a].N0
                  << ", " << paramTet[b].N0 << "): the spin accumulation s being continuous, the "
                     "spin chemical potential is not (N0 mu_s is).\n";
        }
    }

void spinAcc::makeStacks(void)
    {
    stacks.resize(paramTri.size());
    /* triangles of the surfaces with a polarizer stack, by their sorted node indices, to find the
     * tetrahedron on the side of the sample, hence its density of states */
    std::map<std::array<int,Triangle::N>, int> faces;
    for (size_t k = 0; k < msh->tri.size(); k++)
        {
        const Triangle::Tri &f = msh->tri[k];
        if (!paramTri[f.idxPrm].slonczewski)
            { continue; }
        if (paramTri[f.idxPrm].stack.hasNF2 && !msh->isMagnetic(msh->tri[k]))
            {
            std::cerr << "Error: surface region " << paramTri[f.idxPrm].regName << " has a "
                         "polarizer stack with a N/F2 interface, it should be on a magnetic region.\n";
            exit(1);
            }
        std::array<int,Triangle::N> key = f.ind;
        std::sort(key.begin(), key.end());
        faces[key] = k;
        }
    if (faces.empty())
        { return; }

    std::vector<double> N0(paramTri.size(), NAN);
    for (const Tetra::Tet &tet : msh->tet)
        {
        for (int omit = 0; omit < Tetra::N; omit++)  // the face without the node omit
            {
            std::array<int,Triangle::N> key;
            for (int ie = 0, j = 0; ie < Tetra::N; ie++)
                {
                if (ie != omit)
                    { key[j++] = tet.ind[ie]; }
                }
            std::sort(key.begin(), key.end());
            const auto found = faces.find(key);
            if (found == faces.end())
                { continue; }
            const int idx = msh->tri[found->second].idxPrm;
            const double n0 = paramTet[tet.idxPrm].N0;
            if (std::isnan(N0[idx]))
                { N0[idx] = n0; }
            else if (N0[idx] != n0)
                {
                std::cerr << "Error: surface region " << paramTri[idx].regName << " has a polarizer "
                             "stack and borders volume regions of different dens_state.\n";
                exit(1);
                }
            }
        }

    for (size_t idx = 0; idx < paramTri.size(); idx++)
        {
        if (!paramTri[idx].slonczewski)
            { continue; }
        stacks[idx] = Slonczewski::makeStack(paramTri[idx].stack, N0[idx]);
        if (verbose)
            {
            const Slonczewski::Stack &s = stacks[idx];
            std::cout << "polarizer stack of " << paramTri[idx].regName << " (dens_state " << N0[idx]
                      << "), in m/s: F1N (renormalized by F1) P = " << s.P1 << ", Gpar = "
                      << s.G1par << ", Gr = " << s.G1r << ", Gi = " << s.G1i;
            if (s.hasNF2)
                {
                std::cout << "; NF2 P = " << s.P2 << ", Gpar = " << s.G2par << ", Gr = " << s.G2r
                          << ", Gi = " << s.G2i << ", transverse: "
                          << (s.transverseSink ? "sink" : "series");
                }
            std::cout << "\n";
            }
        }
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

std::vector<double> spinAcc::potential(void) const
    {
    std::vector<double> Vt(V.size());
    std::transform(V.begin(), V.end(), Vt.begin(),
                   [this](const double v_ref) { return V_D + jn_factor*v_ref; });
    return Vt;
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
     * (the surfaces where s is given are Dirichlet conditions, see boundaryConditions());
     * time-dependent current density: jn -> jn_factor jn. On a surface with a polarizer stack
     * (slonczewski), the flux is replaced by a Robin condition, see stackBoundaryCondition() */
    for (Triangle::Tri &f : msh->tri)
        {
        if (paramTri[f.idxPrm].slonczewski)
            { stackBoundaryCondition(f); }
        else if (std::isfinite(paramTri[f.idxPrm].jn) && std::isfinite(paramTri[f.idxPrm].uP.norm()))
            {
            std::vector<double> Le(DIM_PB * Triangle::N, 0.0);
            const Eigen::Vector3d Qn = -jn_factor * paramTri[f.idxPrm].jn
                                       * (BOHRS_MUB/CHARGE_ELECTRON) * paramTri[f.idxPrm].uP;
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

void spinAcc::stackBoundaryCondition(Triangle::Tri &f)
    {
    /* the spin current entering the sample through the polarizer stack is
     * q = (mu_B/e) jn Peff(m0) - Gs s (jn > 0: electrons entering the sample, flowing from F1 to F2),
     * that is Q_n = -q: contribution -int_Gamma Q_n a_i = int_Gamma q a_i, the source to the RHS
     * and the term in s to the matrix; time-dependent current density: jn -> jn_factor jn */
    using namespace Triangle;
    const prm &param = paramTri[f.idxPrm];
    const double K = jn_factor * param.jn * BOHRS_MUB / CHARGE_ELECTRON;

    Eigen::Matrix<double,DIM_PB*N,DIM_PB*N> Ke;
    Ke.setZero();
    std::vector<double> Le(DIM_PB * N, 0.0);
    for (int npi = 0; npi < NPI; npi++)
        {
        /* magnetization of the sample at the Gauss point (NEXT, as in the volume terms),
         * normalized; only used with a N/F2 interface, the surface being then magnetic */
        Eigen::Vector3d m0 = Eigen::Vector3d::Zero();
        for (int ie = 0; ie < N; ie++)
            { m0 += a[ie][npi] * msh->getNode_u(f.ind[ie]); }
        if (m0.norm() > 0)
            { m0.normalize(); }
        const Slonczewski::BC bc = stacks[f.idxPrm].bc(m0);

        for (int ie = 0; ie < N; ie++)
            {
            const double ai_w = f.weight[npi] * a[ie][npi];
            for (int di = 0; di < DIM_PB; di++)
                {
                Le[di*N + ie] += K * bc.Peff[di] * ai_w;
                for (int je = 0; je < N; je++)
                    for (int dj = 0; dj < DIM_PB; dj++)
                        { Ke(di*N + ie, dj*N + je) += ai_w * a[je][npi] * bc.Gs(di,dj); }
                }
            }
        }
    buildMat<N>(f.ind, Ke);
    buildVect<N>(f.ind, Le);
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

    /* 1/tau_sf, units s^-1: l_sf is the spin diffusion length, longitudinal in a ferromagnet, where
     * the longitudinal diffusion is D0 (1 - P^2) (beta = beta' = P, see the anisotropic diffusion
     * below): 1/tau_sf = D0 (1 - P^2)/l_sf^2, D0/l_sf^2 in a non-magnetic metal */
    double invTau_sf = D0 / sq(lsf);
    if (msh->isMagnetic(tet))
        { invTau_sf *= 1.0 - sq(getPolarizationRate(tet)); }

    Eigen::Matrix<double,N,1> a_w = eigen_a * tet.weight;
    Eigen::Matrix<double,N,1> diag = invTau_sf * a_w;
    Eigen::Matrix<double,N,N> diagBlock = tet.calcDiagBlock(D0, diag);
    AE.block<N,N>(    0,     0) += diagBlock;
    AE.block<N,N>(    N,     N) += diagBlock;
    AE.block<N,N>(2 * N, 2 * N) += diagBlock;

/* here is the magnetic contribution to AE: the precession, block diagonal and antisymmetric, and
 * the anisotropic diffusion, symmetric with full blocks */
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

        /* anisotropic diffusion (beta = beta' = P): -D0 P^2 (u.grad s):(u.grad w), symmetric, it
         * couples the three components; magnetization at the Gauss points, as in the RHS */
        const Eigen::Matrix<double,Nodes::DIM,NPI> U = u_nod * eigen_a;
        tet.add_spinDiffusionAniso(D0 * sq(getPolarizationRate(tet)), U, AE);
        }
    }

void spinAcc::integrales(const Tetra::Tet &tet, std::vector<double> &BE) const
    {
    using namespace Tetra;

    if(msh->isMagnetic(tet))
        {
        /* magnetization at the Gauss points: u varies in the element, a nodal u_i would not be a
         * lumping */
        const Eigen::Matrix<double,Nodes::DIM,NPI> U = calc_u_nod(tet) * eigen_a;

        /* constant cst0 in a magnetic region is the only RHS parameter involved in the diffusion
        * equation for magnetic contribution
        * units: [cst0] = [sigma] m^2 = A^2 s^3 m^-1 kg^-1
        * time-dependent current density: grad V = jn_factor grad V0
        */
        const double cst0 = jn_factor * BOHRS_MUB * getPolarizationRate(tet) * getSigma(tet)
                            / CHARGE_ELECTRON;
        const Eigen::Matrix<double, Nodes::DIM, NPI> gradV = tet.gradV(V);

        for (size_t npi = 0; npi < NPI; npi++)
            {
            const Eigen::Vector3d cst0_w_gradV = cst0 * tet.weight[npi] * gradV.col(npi);

            for (size_t ie = 0; ie < N; ie++)
                {
                const double tmp = cst0_w_gradV.dot(tet.da.row(ie));
                BE[    ie] += tmp * U(IDX_X,npi);
                BE[  N+ie] += tmp * U(IDX_Y,npi);
                BE[2*N+ie] += tmp * U(IDX_Z,npi);
                }
            }
        }
    }

