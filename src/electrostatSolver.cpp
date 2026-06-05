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
    int nbSurfV(0);
    std::for_each(paramTri.begin(),paramTri.end(),[&nbSurfJ,&nbSurfV](const Triangle::prm &p)
        {
        if (std::isfinite(p.jn)) nbSurfJ++;
        if (std::isfinite(p.V)) nbSurfV++;
        });
    if (!(nbSurfJ == 1 && nbSurfV == 1))
        {
        std::cerr << "Error: incorrect boundary conditions for potential V solver.\n";
        exit(1);
        }
    else if (verbose)
        { std::cout << " electrostatic problem boundary conditions Ok.\n"; }
    }

void electrostatSolver::infos(void) const
    {
    std::cout << "Boundary conditions:\n";
    std::for_each(paramTri.begin(),paramTri.end(),[](const Triangle::prm &p)
        {
        if (!std::isnan(p.jn)) { std::cout << "\tjn= " << p.jn << "\n"; }
        if (!std::isnan(p.V)) { std::cout << "\tV= " << p.V << "\n"; }
        } );
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
            double dai_dx = tet.da(ie,Nodes::IDX_X);
            double dai_dy = tet.da(ie,Nodes::IDX_Y);
            double dai_dz = tet.da(ie,Nodes::IDX_Z);

            for (int je = 0; je < Tetra::N; je++)
                {
                double daj_dx = tet.da(je,Nodes::IDX_X);
                double daj_dy = tet.da(je,Nodes::IDX_Y);
                double daj_dz = tet.da(je,Nodes::IDX_Z);
                AE(ie,je) += sigma_w*(dai_dx*daj_dx + dai_dy*daj_dy + dai_dz*daj_dz);
                }
            }
        }
    }

void electrostatSolver::integrales(const Triangle::Tri &tri, std::vector<double> &BE) const
    {
    double J_bc = getCurrentDensity(tri);// _bc for Boundary Condition
    for (int npi = 0; npi < Triangle::NPI; npi++)
        {
        double J_w = J_bc *tri.weight[npi];
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
            bool iznogood = save(V_fileName,"## columns: index\tV\n");
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

    std::for_each( msh->tet.begin(),msh->tet.end(),[this](Tetra::Tet &elem)
        {
        Eigen::Matrix<double,DIM_PB_ELEC*Tetra::N,DIM_PB_ELEC*Tetra::N> Ke;
        Ke.setZero();
        integrales(elem,Ke);
        buildMat<Tetra::N>(elem.ind,Ke);
        } );

    std::for_each( msh->tri.begin(), msh->tri.end(),[this](Triangle::Tri &elem)
        {
        std::vector<double> Le(DIM_PB*Triangle::N,0.0);
        integrales(elem,Le);
        buildVect<Triangle::N>(elem.ind,Le);
        } );

    std::vector<int> ld; // vector of the Dirichlet Nodes
    std::vector<double> Vd(NOD); // potential values on Dirichlet nodes, zero on the others

    std::for_each(msh->tri.begin(),msh->tri.end(),[this,&Vd,&ld](const Triangle::Tri &tri)
        {
        double _V = paramTri[tri.idxPrm].V;
        if (!std::isnan(_V))
            {
            for(int ie=0; ie<Triangle::N; ie++)
                {
                int i= tri.ind[ie];
                Vd[i]= _V;
                ld.push_back(i);
                }
            }
        });
    suppress_copies<int>(ld);
    iter.reset();
    algebra::cg_dir<double>(iter, K, V, L_rhs, Vd, ld);

    return ( iter.status == algebra::CONVERGED);
    }

bool electrostatSolver::save(const std::string& V_fileName, const std::string &metadata) const
    {
    std::ofstream fout(V_fileName, std::ios::out);
    if (fout.fail())
        {
        std::cerr << "cannot open file " << V_fileName << "\n";
        SYSTEM_ERROR;
        }

    fout << tags::sol::rw_time << ' ' << date() << '\n' << metadata << std::scientific
         << std::setprecision(precision);
    const int NOD = msh->getNbNodes();
    for (int i = 0; i < NOD; i++)
        {
        int _i = msh->getNodeIndex(i);
        fout << i << '\t' << V[_i] << "\n";
        }
    fout.close();
    return !(fout.good());
    }

double electrostatSolver::getSigma(const Tetra::Tet &tet) const
    { return paramTet[tet.idxPrm].sigma; }

double electrostatSolver::getCurrentDensity(const Triangle::Tri &tri) const
    {
    double val = paramTri[tri.idxPrm].jn;
    if (!std::isfinite(val)) val = 0.0;
    return val;
    }
