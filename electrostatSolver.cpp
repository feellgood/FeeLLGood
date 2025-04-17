#include <fstream>
#include <iostream>
#include "chronometer.h" // date()
#include "tags.h"
#include "electrostatSolver.h"

void electrostatSolver::calc_gradV(Tetra::Tet const &tet, Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> _gradV)
    {
    for (int npi = 0; npi < Tetra::NPI; npi++)
        {
        Eigen::Vector3d v(0,0,0);
        for (int i = 0; i < Tetra::N; i++)
            { v += V[tet.ind[i]] * tet.da.row(i); }
        _gradV.col(npi) = v;
        }
    }

void electrostatSolver::calc_Hm(Tetra::Tet const &tet, Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> _gradV,
                 Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> _Hm)
    {
    Eigen::Matrix<double,Nodes::DIM,Tetra::NPI> p_g;
    tet.getPtGauss(p_g);
    const double sigma = getSigma(tet);
    for (int npi = 0; npi < Tetra::NPI; npi++)
        { _Hm.col(npi) = -sigma * _gradV.col(npi).cross(p_g.col(npi)); }
    }

void electrostatSolver::infos(void)
    {
    std::cout << "Boundary conditions:\n";
    std::for_each(paramFacette.begin(),paramFacette.end(),[](Facette::prm &p)
        {
        if (!std::isnan(p.J)) { std::cout << "\tJ= " << p.J << std::endl; }
        if (!std::isnan(p.V)) { std::cout << "\tV= " << p.V << std::endl; }
        } );
    }

void electrostatSolver::assembling_mat(Tetra::Tet const &tet, double Ke[Tetra::N][Tetra::N], std::vector<Eigen::Triplet<double>> &K)
    {
    for (int ie = 0; ie < Tetra::N; ie++)
        {
        for (int je = 0; je < Tetra::N; je++)
            {
            double val = Ke[ie][je];
            if (val != 0)
                { K.push_back( Eigen::Triplet(tet.ind[ie], tet.ind[je], val) ); }
            }
        }
    }

void electrostatSolver::assembling_vect(Facette::Fac const &fac, std::vector<double> const &Le, Eigen::Ref<Eigen::VectorXd> L)
    {
    for (int ie = 0; ie < Facette::N; ie++)
        { L(fac.ind[ie]) += Le[ie]; }
    }

// same formula as sp_acc_llg
void electrostatSolver::integrales(Tetra::Tet const &tet, double AE[Tetra::N][Tetra::N])
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
                AE[ie][je] += sigma_w*(dai_dx*daj_dx + dai_dy*daj_dy + dai_dz*daj_dz);
                }
            }
        }
    }

// same formula as sp_acc_llg
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

void electrostatSolver::prepareData(std::vector<Eigen::Triplet<double>> &Kw, Eigen::Ref<Eigen::VectorXd> Lw)
    {
    std::for_each(msh.tet.begin(), msh.tet.end(), [this, &Kw](Tetra::Tet const &tet)
                  {
                  double K[Tetra::N][Tetra::N] = {{0}};
                  integrales(tet, K);
                  assembling_mat(tet, K, Kw);
                  });

    std::for_each(msh.fac.begin(), msh.fac.end(), [this, &Lw](Facette::Fac const &fac)
                  {
                  std::vector<double> L(Facette::N);
                  integrales(fac, L);
                  assembling_vect(fac, L, Lw);
                  });
    }

int electrostatSolver::solve(const double _tol)
    {
    const int NOD = msh.getNbNodes();

    std::vector<Eigen::Triplet<double>> Kw;
    Eigen::VectorXd Lw = Eigen::VectorXd::Zero(NOD);
    prepareData(Kw, Lw);

    Eigen::VectorXd Xw(NOD);
    /* do something here with boundary conditions */

    if (verbose)
        { std::cout << "line weighting..." << std::endl; }

    Eigen::SparseMatrix<double> Kr(NOD,NOD);
    Kr.setFromTriplets(Kw.begin(),Kw.end());
    std::vector<double> maxRow(NOD,0);

    //here we fill the vector maxRow with infinity norm : maxRow[i] = max{|Kr.row(i)|}
    for(int i=0;i<NOD;i++)
        {
        for(int k=0;k<Kr.outerSize();++k)
            for(Eigen::SparseMatrix<double>::InnerIterator it(Kr,k); it; ++it)
                { if((it.row() == i)&&(fabs(it.value()) > maxRow[i])) { maxRow[i] = fabs(it.value()); } }
        }

    for (int i = 0; i < NOD; i++)
        {
        double norme = maxRow[i];
        Lw[i] /= norme;
        Kr.row(i) /= norme;
        }

    if (verbose)
        { std::cout << "solving ..." << std::endl; }

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> _solver;

    _solver.setTolerance(_tol);
    _solver.setMaxIterations(MAXITER);
    _solver.compute(Kr);

    Eigen::VectorXd sol = _solver.solve(Lw);
    for (int i=0;i<NOD;i++)
        { V[i]= sol(i); }
    return (_solver.iterations() < MAXITER);
    }

bool electrostatSolver::save(std::string const &metadata) const
    {
    std::ofstream fout(fileName, std::ios::out);
    if (fout.fail())
        {
        std::cout << "cannot open file " << fileName << std::endl;
        SYSTEM_ERROR;
        }

    fout << tags::sol::rw_time << ' ' << date() << '\n' << metadata << std::scientific
         << std::setprecision(precision);

    const int NOD = msh.getNbNodes();

    for (int i = 0; i < NOD; i++)
        {
        int _i = msh.getNodeIndex(i);
        fout << i << '\t' << V[_i] << std::endl;
        }
    fout.close();
    return !(fout.good());
    }
