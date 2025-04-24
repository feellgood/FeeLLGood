#include <fstream>
#include <iostream>
#include "chronometer.h" // date()
#include "tags.h"
#include "electrostatSolver.h"
#include "algebra/cg.h"

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

// same formula as sp_acc_llg
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

// same formula as sp_acc_llg
void electrostatSolver::integrales(Facette::Fac const &fac, Eigen::Ref<Eigen::Matrix<double,Facette::N,1> > BE)
    {
    double J_bc = getCurrentDensity(fac);// _bc for Boundary Condition
    for (int npi = 0; npi < Facette::NPI; npi++)
        {
        double J_w = J_bc *fac.weight[npi];
        for (int ie = 0; ie < Facette::N; ie++)
            { BE(ie) -= Facette::a[ie][npi] * J_w; }
        }
    }

int electrostatSolver::solve(const double _tol)
    {
    const int DIM_1D = 1;
    algebra::w_sparseMat Kw(NOD);
    std::vector<double> Lw(NOD, 0.0);
    std::vector<double> Xw(NOD, 0.0);//Eigen::VectorXd Xw(NOD);
    
    if (verbose)
        { std::cout << "line weighting..." << std::endl; }

    std::for_each( msh.tet.begin(),msh.tet.end(),[this,&Kw,&Lw](Tetra::Tet &elem)
        {
        Eigen::Matrix<double,DIM_1D*Tetra::N,DIM_1D*Tetra::N> K;
        K.setZero();
        Eigen::Matrix<double,DIM_1D*Tetra::N,1> L;
        L.setZero();
        integrales(elem,K);
        assembling<Tetra::N>(elem.ind,K,L,Kw,Lw);
        } );

    std::for_each( msh.fac.begin(), msh.fac.end(),[this,&Kw,&Lw](Facette::Fac &elem)
        {
        Eigen::Matrix<double,DIM_1D*Facette::N,DIM_1D*Facette::N> K;
        K.setZero();
        Eigen::Matrix<double,DIM_1D*Facette::N,1> L;
        L.setZero();
        integrales(elem,L);
        assembling<Facette::N>(elem.ind,K,L,Kw,Lw);
        } );

    if (verbose)
        { std::cout << "solving ..." << std::endl; }

    algebra::r_sparseMat Kr(Kw);
    algebra::iteration iter("cg_dir",_tol,verbose,MAXITER);
    std::vector<int> ld; // vector of the Dirichlet Nodes
    std::vector<double> Vd(NOD); // potential values on Dirichlet nodes, zero on the others

    std::for_each(msh.fac.begin(),msh.fac.end(),[this,&Vd,&ld](Facette::Fac &fac)
        {
        double _V = paramFacette[fac.idxPrm].V;
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
    algebra::cg_dir<double>(iter, Kr, Xw, Lw, Vd, ld);

    if (verbose)
        std::cout << "solved in " << iter.get_iteration() <<" ; residu= " << iter.get_res();
    for (int i=0;i<NOD;i++)
        { V[i]= Xw[i]; }
    return ( iter.get_iteration() < (int)MAXITER);
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

    for (int i = 0; i < NOD; i++)
        {
        int _i = msh.getNodeIndex(i);
        fout << i << '\t' << V[_i] << std::endl;
        }
    fout.close();
    return !(fout.good());
    }

double electrostatSolver::getSigma(Tetra::Tet const &tet) const
    { return paramTetra[tet.idxPrm].sigma; }

double electrostatSolver::getCurrentDensity(Facette::Fac const &fac) const
    {
    double val = paramFacette[fac.idxPrm].J;
    if (!std::isfinite(val)) val = 0.0;
    return val;
    }
