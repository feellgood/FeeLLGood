#include <fstream>
#include <iostream>
#include "chronometer.h" // date()
#include "tags.h"
#include "electrostatSolver.h"
#include "algebra/cg.h"
#include "solver.h"  // build(Mat|Vect)
#include "meshUtils.h"

bool electrostatSolver::checkBoundaryConditions(void) const
    {
    int nbSurfJ(0);
    int nbSurfV(0);
    std::for_each(paramFac.begin(),paramFac.end(),[&nbSurfJ,&nbSurfV](Facette::prm const &p)
        {
        if (std::isfinite(p.jn)) nbSurfJ++;
        if (std::isfinite(p.V)) nbSurfV++;
        });
    return ((nbSurfJ == 1)&&(nbSurfV == 1));
    }

void electrostatSolver::infos(void)
    {
    std::cout << "Boundary conditions:\n";
    std::for_each(paramFac.begin(),paramFac.end(),[](const Facette::prm &p)
        {
        if (!std::isnan(p.jn)) { std::cout << "\tjn= " << p.jn << std::endl; }
        if (!std::isnan(p.V)) { std::cout << "\tV= " << p.V << std::endl; }
        } );
    }

void electrostatSolver::integrales(Tetra::Tet const &tet, Eigen::Ref<Eigen::Matrix<double,Tetra::N,Tetra::N> > AE)
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

void electrostatSolver::integrales(Facette::Fac const &fac, std::vector<double> &BE)
    {
    double J_bc = getCurrentDensity(fac);// _bc for Boundary Condition
    for (int npi = 0; npi < Facette::NPI; npi++)
        {
        double J_w = J_bc *fac.weight[npi];
        for (int ie = 0; ie < Facette::N; ie++)
            { BE[ie] -= Facette::a[ie][npi] * J_w; }
        }
    }

void electrostatSolver::compute(const bool verbose, const std::string V_fileName)
    {
    bool has_converged = solve();
    if (verbose)
        { std::cout << "electrostatic solver: " << iter.infos() << std::endl; }
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
    algebra::w_sparseMat Kw(NOD);
    std::vector<double> Lw(NOD, 0.0);

    std::for_each( msh->tet.begin(),msh->tet.end(),[this,&Kw](Tetra::Tet &elem)
        {
        Eigen::Matrix<double,DIM_PB_ELEC*Tetra::N,DIM_PB_ELEC*Tetra::N> K;
        K.setZero();
        integrales(elem,K);
        buildMat<Tetra::N>(elem.ind,K,Kw);
        } );

    std::for_each( msh->fac.begin(), msh->fac.end(),[this,&Lw](Facette::Fac &elem)
        {
        std::vector<double> L(DIM_PB*Facette::N,0.0);
        integrales(elem,L);
        buildVect<Facette::N>(elem.ind,L,Lw);
        } );

    algebra::r_sparseMat Kr(Kw);
    std::vector<int> ld; // vector of the Dirichlet Nodes
    std::vector<double> Vd(NOD); // potential values on Dirichlet nodes, zero on the others

    std::for_each(msh->fac.begin(),msh->fac.end(),[this,&Vd,&ld](Facette::Fac &fac)
        {
        double _V = paramFac[fac.idxPrm].V;
        if (!std::isnan(_V))
            {
            for(int ie=0; ie<Facette::N; ie++)
                {
                int i= fac.ind[ie];
                Vd[i]= _V;
                ld.push_back(i);
                }
            }
        });
    suppress_copies<int>(ld);
    iter.reset();
    algebra::cg_dir<double>(iter, Kr, V, Lw, Vd, ld);

    return ( iter.status == algebra::CONVERGED);
    }

bool electrostatSolver::save(const std::string V_fileName, std::string const &metadata) const
    {
    std::ofstream fout(V_fileName, std::ios::out);
    if (fout.fail())
        {
        std::cout << "cannot open file " << V_fileName << std::endl;
        SYSTEM_ERROR;
        }

    fout << tags::sol::rw_time << ' ' << date() << '\n' << metadata << std::scientific
         << std::setprecision(precision);
    const int NOD = msh->getNbNodes();
    for (int i = 0; i < NOD; i++)
        {
        int _i = msh->getNodeIndex(i);
        fout << i << '\t' << V[_i] << std::endl;
        }
    fout.close();
    return !(fout.good());
    }

double electrostatSolver::getSigma(Tetra::Tet const &tet) const
    { return paramTet[tet.idxPrm].sigma; }

double electrostatSolver::getCurrentDensity(Facette::Fac const &fac) const
    {
    double val = paramFac[fac.idxPrm].jn;
    if (!std::isfinite(val)) val = 0.0;
    return val;
    }
