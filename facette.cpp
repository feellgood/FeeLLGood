#include "facette.h"

using namespace Facette;
using namespace Nodes;

void Fac::integrales(Facette::prm const &params)
    {
    double Kbis = 2.0*params.Ks / Ms;
    
    Eigen::Matrix<double,DIM,NPI> u;
    interpolation(Nodes::get_u0, u);

    Eigen::Matrix<double,DIM,N> BE;
    BE.setZero();

    for (int npi = 0; npi < NPI; npi++)
        {
        double _prefactor = weight[npi] * Kbis * params.uk.dot(u.col(npi));

        for (int i = 0; i < N; i++)
            for(int k = 0;k<DIM;k++)
                {
                BE(k,i) += _prefactor*a[i][npi]*params.uk(k);
                }
        }
    /*-------------------- PROJECTION --------------------*/
    #if EIGEN_VERSION_AT_LEAST(3,4,0)
        Lp = P * BE.reshaped<Eigen::RowMajor>();
    #else
        Eigen::Matrix<double,DIM*N,1> tmp;
        for(int k=0;k<DIM;k++)
            for(int i=0;i<N;i++)
                tmp(k*N+i) = BE(k,i);
        Lp = P*tmp;
    #endif
    //Lp = P*BE;
    }

double Fac::anisotropyEnergy(Facette::prm const &param,
                             Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> const u) const
    {  // surface Neel anisotropy (uk is a uniaxial easy axis)
    Eigen::Matrix<double,NPI,1> dens = (u.transpose()*param.uk).array().square();
    return -param.Ks*weight.dot(dens);
    }

Eigen::Matrix<double,NPI,1> Fac::charges(std::function<Eigen::Vector3d(Nodes::Node)> getter, std::vector<double> &corr) const
    {
    Eigen::Matrix<double,DIM,N> vec_nod;
    for(int i=0;i<N;i++)
        vec_nod.col(i) << getter(getNode(i));

    Eigen::Matrix<double,DIM,NPI> _u = vec_nod * eigen_a;

    Eigen::Matrix<double,NPI,1> result = Ms*weight.cwiseProduct( _u.transpose()*n );
    Eigen::Matrix<double,DIM,NPI> gauss;
    getPtGauss(gauss);
    // calc corr node by node
    for (int i = 0; i < N; i++)
        {
        const int i_ = ind[i];
        const Eigen::Vector3d &p_i_ = getNode(i).p;
        for (int j = 0; j < NPI; j++)
            {
            double d_ij= (p_i_ - gauss.col(j)).norm();
            corr[i_] -= result(j)/d_ij;//Ms * pScal(u[j], n) * weight(j) / d_ij;
            }
        corr[i_] += potential(getter, i);
        }
    
    return result;
    }

double Fac::demagEnergy(Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> u, Eigen::Ref<Eigen::Matrix<double,NPI,1>> phi) const
    {
    Eigen::Matrix<double,NPI,1> dens = (u.transpose()*n).cwiseProduct(phi);
    return 0.5*mu0*Ms*dens.dot(weight);
    }

double Fac::potential(std::function<Eigen::Vector3d(Nodes::Node)> getter, int i) const
    {
    int ii = (i + 1) % 3;
    int iii = (i + 2) % 3;

    Nodes::Node const &node1 = getNode(i);
    Nodes::Node const &node2 = getNode(ii);
    Nodes::Node const &node3 = getNode(iii);

    Eigen::Vector3d p1p2 = node2.p - node1.p;
    Eigen::Vector3d p1p3 = node3.p - node1.p;

    std::function<double(double)> f = [](double x) { return sqrt(1.0 + x * x); };

    double b = p1p2.norm();
    double t = p1p2.dot(p1p3) / b;  // carefull with t , if cos(p1p2,p1p3) < 0 then t < 0
    double _2s = 2. * surf;
    double h = _2s / b;

    if (_2s < 0)
        {
        std::cout << "facette surface is negative : surface is ill-oriented" << std::endl;
        exit(1);
        }

    double c = (t - b) / h;
    double r = h*f(t/h);  // if _2s > 0 it is the same as double r = sqrt( sq(_2s/b) + sq(t));
    double log_1 = log((c * t + h + f(c) * r) / (b * (c + f(c))));
    double xi = b * log_1 / f(c);

    double s1 = getter(node1).dot(n);
    double s2 = getter(node2).dot(n);
    double s3 = getter(node3).dot(n);

    double pot = xi * s1
                 + ((xi * (h + c * t) - b * (r - b)) * s2 + b * (r - b - c * xi) * s3) * b
                           / (_2s * (1 + c * c));
    return 0.5 * Ms * pot;
    }

