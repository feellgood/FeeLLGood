#include <algorithm>
#include <fstream>
#include <iostream>
#include "chronometer.h" // date()
#include "tags.h"
#include "electrostatSolver.h"
#include "algebra/cg.h"
#include "solver.h"  // build(Mat|Vect)
#include "meshUtils.h"

void electrostatSolver::checkBoundaryConditions(void) const
    {
    int nbSurfJ(0);
    int nbSurfJ_current_density(0);
    int nbSurfV(0);
    for (const Triangle::prm &p : paramTri)
        {
        if (std::isfinite(p.jn))
            { nbSurfJ++; }
        if (p.jn_current_density)
            { nbSurfJ_current_density++; }
        if (std::isfinite(p.V))
            { nbSurfV++; }
        }
    /* at least one surface with jn, all numeric or all jn: current_density (the potential is then
     * proportional to the single factor current_density(t)/j_ref), and a single surface with V */
    if (!(nbSurfJ >= 1 && nbSurfV == 1))
        {
        std::cerr << "Error: incorrect boundary conditions for potential V solver: at least one "
                     "surface with jn and a single surface with V are required.\n";
        exit(1);
        }
    if (nbSurfJ_current_density != 0 && nbSurfJ_current_density != nbSurfJ)
        {
        std::cerr << "Error: incorrect boundary conditions for potential V solver: the surfaces "
                     "with jn cannot mix jn: current_density and numeric jn.\n";
        exit(1);
        }
    else if (verbose)
        { std::cout << " electrostatic problem boundary conditions Ok.\n"; }
    }

void electrostatSolver::infos(void) const
    {
    std::cout << "Boundary conditions:\n";
    for (const Triangle::prm &p : paramTri)
        {
        if (p.jn_current_density)
            { std::cout << "\tjn= current_density (reference 1 A/m^2)\n"; }
        else if (!std::isnan(p.jn))
            { std::cout << "\tjn= " << p.jn << "\n"; }
        if (!std::isnan(p.V))
            { std::cout << "\tV= " << p.V << "\n"; }
        }
    }

void electrostatSolver::integrales(const Tetra::Tet &tet,
        Eigen::Ref<Eigen::Matrix<double,Tetra::N,Tetra::N> > AE) const
    {
    const double sigma = getSigma(tet);
    for (int npi = 0; npi < Tetra::NPI; npi++)
        {
        double sigma_w = sigma * tet.weight[npi];

        for (int ie = 0; ie < Tetra::N; ie++)
            {
            double dai_dx = tet.da(ie, Nodes::IDX_X);
            double dai_dy = tet.da(ie, Nodes::IDX_Y);
            double dai_dz = tet.da(ie, Nodes::IDX_Z);

            for (int je = 0; je < Tetra::N; je++)
                {
                double daj_dx = tet.da(je, Nodes::IDX_X);
                double daj_dy = tet.da(je, Nodes::IDX_Y);
                double daj_dz = tet.da(je, Nodes::IDX_Z);
                AE(ie, je) += sigma_w * (dai_dx * daj_dx + dai_dy * daj_dy + dai_dz * daj_dz);
                }
            }
        }
    }

void electrostatSolver::integrales(const Triangle::Tri &tri, std::vector<double> &BE) const
    {
    double J_bc = getCurrentDensity(tri); // _bc for Boundary Condition
    for (int npi = 0; npi < Triangle::NPI; npi++)
        {
        double J_w = J_bc * tri.weight[npi];
        for (int ie = 0; ie < Triangle::N; ie++)
            { BE[ie] -= Triangle::a[ie][npi] * J_w; }
        }
    }

void electrostatSolver::compute(const bool verbose, const std::string& V_fileName)
    {
    bool has_converged = solve();
    if (verbose)
        { std::cout << "electrostatic solver: " << iter.infos() << "\n"; }
    if (has_converged)
        {
        if (!V_fileName.empty())
            {
            bool iznogood = save(V_fileName, "## columns: index\tV\n");
            if (verbose && iznogood)
                { std::cout << "file " << V_fileName << " written.\n"; }
            }
        }
    else
        { exit(1); }
    }

bool electrostatSolver::solve(void)
    {
    const int NOD = msh->getNbNodes();

    for (Tetra::Tet &elem : msh->tet)
        {
        Eigen::Matrix<double,DIM_PB_ELEC*Tetra::N,DIM_PB_ELEC*Tetra::N> Ke;
        Ke.setZero();
        integrales(elem, Ke);
        buildMat<Tetra::N>(elem.ind, Ke);
        }

    for (Triangle::Tri &elem : msh->tri)
        {
        std::vector<double> Le(DIM_PB * Triangle::N, 0.0);
        integrales(elem, Le);
        buildVect<Triangle::N>(elem.ind, Le);
        }

    std::vector<int> ld; // vector of the Dirichlet Nodes
    std::vector<double> Vd(NOD); // potential values on Dirichlet nodes, zero on the others

    /* the potential is solved with V = 0 on the surface with V, the imposed V_D being added
     * afterwards (spinAcc::potential, save): the potential is defined up to a constant, and
     * solving around V_D != 0 would add an error of the order of the tolerance times V_D (it
     * drowned the reference potential of jn: current_density, about 5e-14 V for 1 A/m^2) */
    reference_problem = std::any_of(paramTri.begin(), paramTri.end(),
                                    [](const Triangle::prm &p) { return p.jn_current_density; });

    for (const Triangle::Tri &tri : msh->tri)
        {
        double _V = paramTri[tri.idxPrm].V;
        if (!std::isnan(_V))
            {
            V_D = _V;
            for(int ie=0; ie < Triangle::N; ie++)
                {
                int i= tri.ind[ie];
                Vd[i]= 0.0;
                ld.push_back(i);
                }
            }
        }
    suppress_copies<int>(ld);
    iter.reset();
    algebra::cg_dir<double>(iter, K, V, L_rhs, Vd, ld);

    return (iter.status == algebra::CONVERGED);
    }

bool electrostatSolver::save(const std::string &V_fileName, const std::string &metadata) const
    {
    std::ofstream fout(V_fileName, std::ios::out);
    if (fout.fail())
        {
        std::cerr << "cannot open file " << V_fileName << "\n";
        SYSTEM_ERROR;
        }

    fout << tags::common::rw_time << ' ' << date() << '\n' << metadata << std::scientific
         << std::setprecision(precision);
    // V0 = V_D + V for a numeric jn; reference potential V_ref = V for jn: current_density
    const double offset = reference_problem ? 0.0 : V_D;
    const int NOD = msh->getNbNodes();
    for (int i = 0; i < NOD; i++)
        {
        int _i = msh->getNodeIndex(i);
        fout << i << '\t' << offset + V[_i] << "\n";
        }
    fout.close();
    return !(fout.good());
    }

double electrostatSolver::getSigma(const Tetra::Tet &tet) const
    { return paramTet[tet.idxPrm].sigma; }

double electrostatSolver::getCurrentDensity(const Triangle::Tri &tri) const
    {
    double val = paramTri[tri.idxPrm].jn;
    if (!std::isfinite(val))
        { val = 0.0; }
    return val;
    }
