/**
  Elementary matrix Calculation for a tetrahedron element
 */

#include <set>

#include "config.h"  // to get gamma0 constant
#include "tetra.h"
#include "time_integration.h"
#include "facette.h"

using namespace Tetra;
using namespace Nodes;

Eigen::Matrix<double,NPI,1> Tetra::calc_alpha_eff(const double dt, const double alpha,
                                                      Eigen::Ref<Eigen::Matrix<double,NPI,1>> uHeff)
    {
    double reduced_dt = gamma0 * dt;
    Eigen::Matrix<double,NPI,1> a_eff;
    a_eff.setConstant(alpha);
    const double r = 0.1;  // what is that constant, where does it come from ?
    const double M = 2. * alpha * r / reduced_dt;

    for(int npi=0;npi<NPI;npi++)
        {
        double h = uHeff(npi);
        if (h > 0.)
            {
            if (h > M)
                a_eff(npi) = alpha + reduced_dt / 2. * M;
            else
                a_eff(npi) = alpha + reduced_dt / 2. * h;
            }
        else
            {
            if (h < -M)
                a_eff(npi) = alpha / (1. + reduced_dt / (2. * alpha) * M);
            else
                a_eff(npi) = alpha / (1. - reduced_dt / (2. * alpha) * h);
            }
        }
    return a_eff;
    }

void Tet::lumping(Eigen::Ref<Eigen::Matrix<double,NPI,1>> alpha_eff, double prefactor,
                  Eigen::Ref<Eigen::Matrix<double,3*N,3*N>> AE ) const
    {
    Eigen::Matrix<double,N,N> exch_block = da*da.transpose();
    exch_block *= prefactor*weight.sum();
    // contrib is alpha contribution to the diagonal of AE; to help stabilizing the scheme, alpha is modified
    Eigen::Matrix<double,N,1> contrib = eigen_a * weight.cwiseProduct(alpha_eff);
    exch_block.diagonal() += contrib;
    
    AE.block<N,N>(0,0) += exch_block;
    AE.block<N,N>(N,N) += exch_block;
    AE.block<N,N>(2*N,2*N) += exch_block;
      
    // off diagonal blocks are filled with magnetization components weighted products
    // it could be rewritten using block matrix and avoid getNode().get_u() in the loop 
    Eigen::Matrix<double,N,1> a_w = eigen_a * weight; 
    for (int i = 0; i < N; i++)
        {
        const Eigen::Vector3d ai_w_u0 = a_w[i] * getNode(i).get_u(Nodes::CURRENT);
        AE(i, 2*N + i) += ai_w_u0(IDX_Y);
        AE(i, N + i) -= ai_w_u0(IDX_Z);
        AE(N + i, i) += ai_w_u0(IDX_Z);
        AE(N + i, 2*N + i) -= ai_w_u0(IDX_X);
        AE(2*N + i, N + i) += ai_w_u0(IDX_X);
        AE(2*N + i, i) -= ai_w_u0(IDX_Y);
        }
    }

void Tet::add_drift_BE(double alpha, double s_dt, double Vdrift,
                       Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> U,
                       Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> V, 
                       Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> dUd_,
                       Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> dVd_,
                       Eigen::Ref<Eigen::Matrix<double,DIM,N>> BE) const
    {  // the artificial drift from eventual recentering is along x,y or z
    Eigen::Matrix<double,DIM,N> interim;

    for (int npi = 0; npi < NPI; npi++)
        {
        for (int i = 0; i < N; i++)
            {
            interim.col(i) = a[i][npi]*( alpha*dUd_.col(npi) + U.col(npi).cross(dUd_.col(npi))
                    + s_dt*(alpha*dVd_.col(npi) + U.col(npi).cross(dVd_.col(npi)) + V.col(npi).cross(dUd_.col(npi))) );
            }
        BE += Vdrift*weight[npi]*interim;
        }
    }

Eigen::Matrix<double,NPI,1> Tet::calc_aniso_uniax(Eigen::Ref<const Eigen::Vector3d> uk, const double Kbis,
                             const double s_dt,
                             Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> U,
                             Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> V,
                             Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> H_aniso) const
    {
    for(int npi = 0;npi<NPI;npi++)
        { H_aniso.col(npi) += (Kbis * uk.dot( U.col(npi) + s_dt * V.col(npi))) * uk; }
    return Kbis * (U.transpose() * uk).array().square();
    }

Eigen::Matrix<double,NPI,1> Tet::calc_aniso_cub(Eigen::Ref<const Eigen::Vector3d> ex,
                           Eigen::Ref<const Eigen::Vector3d> ey,
                           Eigen::Ref<const Eigen::Vector3d> ez,
                           const double K3bis, const double s_dt,
                           Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> U,
                           Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> V,
                           Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> H_aniso) const
    {
    Eigen::Matrix<double,NPI,1> result;
    for(int npi = 0;npi<NPI;npi++)
        {
        Eigen::Vector3d uk_u = Eigen::Vector3d(ex.dot(U.col(npi)), ey.dot(U.col(npi)), ez.dot(U.col(npi)));
        Eigen::Vector3d uk_v = Eigen::Vector3d(ex.dot(V.col(npi)), ey.dot(V.col(npi)), ez.dot(V.col(npi)));
        Eigen::Vector3d uk_uuu = uk_u.unaryExpr( [](double x){ return x*(1.0 - x*x);} );
        Eigen::Vector3d tmp = uk_v.cwiseProduct(ex);

        H_aniso.col(npi) += -K3bis * (uk_uuu(0) * ex + uk_uuu(1) * ey + uk_uuu(2) * ez
                  + s_dt * tmp.cwiseProduct( Eigen::Vector3d(1, 1, 1) - 3*uk_u.cwiseProduct(uk_u) ));

        result[npi] = uk_u.dot(uk_uuu);
        }
    return -K3bis*result;
    }

void Tet::integrales(Tetra::prm const &param, timing const &prm_t,
                     std::function<void( Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> Hext )> calc_Hext,//Eigen::Vector3d const &Hext,
                     Nodes::index idx_dir, double Vdrift)
    {
    const double alpha = param.alpha_LLG;
    const double Js = param.J;
    const double Abis = 2.0 * param.A / Js;
    const double dt = prm_t.get_dt();
    const double s_dt = THETA * dt * gamma0;  // theta from theta scheme in config.h.in

    /*-------------------- INTERPOLATION --------------------*/
    Eigen::Matrix<double,DIM,NPI> U,dUdx,dUdy,dUdz;
    interpolation(Nodes::get_u<CURRENT>, U, dUdx, dUdy, dUdz);
    Eigen::Matrix<double,Nodes::DIM,N> v_nod;
    for (int i = 0; i < N; i++)
        { v_nod.col(i) = Nodes::get_v<CURRENT>(getNode(i)); }
    Eigen::Matrix<double,DIM,NPI> V = v_nod * eigen_a;
    Eigen::Matrix<double,DIM,NPI> Hd;
    interpolation_field(Nodes::get_phi<CURRENT>, Hd);
    Eigen::Matrix<double,DIM,NPI> Hv;
    interpolation_field(Nodes::get_phiv<CURRENT>, Hv);
    /*-------------------- END INTERPOLATION ----------------*/
    Eigen::Matrix<double,NPI,1> uHeff = -Abis *( dUdx.colwise().squaredNorm() + dUdy.colwise().squaredNorm() + dUdz.colwise().squaredNorm());
    Eigen::Matrix<double,DIM,NPI> H_aniso;
    H_aniso.setZero();

    if(param.K != 0)
        {
        double Kbis = 2.0 * param.K / Js;
        uHeff += calc_aniso_uniax(param.uk, Kbis, s_dt, U, V, H_aniso);
        }
    if(param.K3 != 0)
        {
        double K3bis = 2.0 * param.K3 / Js;
        uHeff += calc_aniso_cub(param.ex, param.ey, param.ez, K3bis, s_dt, U, V, H_aniso);
        }

    Eigen::Matrix<double,DIM,NPI> Hext;
    calc_Hext(Hext);// Hext is defined on gauss points, calc_Hext do a = on Hext

    Eigen::Matrix<double,DIM,NPI> Heff = Hd + Hext;
    Eigen::Matrix<double,DIM,NPI> H = Heff; // we need Hd+Hext for future computations
    extraField(Heff);// extraField do a +=like for STT contrib on Heff
    uHeff += (U.cwiseProduct(Heff)).colwise().sum();//dot product on each col of U and Heff

    Eigen::Matrix<double,NPI,1> a_eff = calc_alpha_eff(dt, alpha, uHeff);
    
    Eigen::Matrix<double,3*N,3*N> AE;
    AE.setZero();
    lumping(a_eff, prm_t.prefactor * s_dt * Abis, AE);

    Eigen::Matrix<double,2*N,3*N> P;
    buildMatP(P);

    /*--------------------   PROJECTION: AE->Kp   --------------------*/
    Kp = P*AE*P.transpose();// with MKL installed this operation should call dgemm_direct
    
    /*********************   calculations on BE   *********************/
    Eigen::Matrix<double,DIM,N> BE;
    BE.setZero();
    if (idx_dir != IDX_UNDEF)
        {
        Eigen::Matrix<double,DIM,NPI>  dVd_dir = v_nod * (da.col(idx_dir)).replicate(1,NPI);
        if (idx_dir == IDX_Z)
            add_drift_BE(alpha, s_dt, Vdrift, U, V, dUdz, dVd_dir, BE);
        else if (idx_dir == IDX_Y)
            add_drift_BE(alpha, s_dt, Vdrift, U, V, dUdy, dVd_dir, BE);
        else if (idx_dir == IDX_X)
            add_drift_BE(alpha, s_dt, Vdrift, U, V, dUdx, dVd_dir, BE);
        }
    H += H_aniso + (s_dt / gamma0) * Hv;
    
    for (int npi = 0; npi < NPI; npi++)
        {
        const double w = weight[npi];
        for (int i = 0; i < N; i++)
            {
            BE.col(i) -= w*Abis*(da(i,0)*dUdx.col(npi) + da(i,1)*dUdy.col(npi) + da(i,2)*dUdz.col(npi));
            BE.col(i) += w*a[i][npi]*H.col(npi);
            }
        }
    extraCoeffs_BE(Js, U, dUdx, dUdy, dUdz, BE);  // STT

    /*--------------------   PROJECTION: BE->Lp   --------------------*/
    #if EIGEN_VERSION_AT_LEAST(3,4,0)
        Lp = P * BE.reshaped<Eigen::RowMajor>();
    #else
        Eigen::Matrix<double,DIM*N,1> tmp;
        for(int k=0;k<DIM;k++)
            for(int i=0;i<N;i++)
                tmp(k*N+i) = BE(k,i);
        Lp = P*tmp;
    #endif
    }

double Tet::exchangeEnergy(Tetra::prm const &param,
                           Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> dudx,
                           Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> dudy,
                           Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> dudz) const
    {// partial reduction on columns with colwise()
    Eigen::Matrix<double,NPI,1> dens = dudx.colwise().squaredNorm()
                                     + dudy.colwise().squaredNorm()
                                     + dudz.colwise().squaredNorm();
    return param.A * weight.dot(dens);
    }

double Tet::anisotropyEnergy(Tetra::prm const &param, Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> u) const
    {
    Eigen::Matrix<double,NPI,1> dens;

    for (int npi = 0; npi < NPI; npi++)
        {
        Eigen::Vector3d m = u.col(npi);
        // uniaxial magnetocrystalline anisotropy constant K, anisotropy axis uk
        dens[npi] = -param.K * sq( param.uk.dot( m ));

        // cosinus directeurs
        double al0 = m.dot(param.ex);
        double al1 = m.dot(param.ey);
        double al2 = m.dot(param.ez);

        dens[npi] += param.K3
                     * (sq(al0 * al1) + sq(al1 * al2) + sq(al2 * al0));  // cubic anisotropy (K3)
        }

    return weight.dot(dens);
    }

void Tet::charges(Tetra::prm const &param,
                  std::function<Eigen::Vector3d(Nodes::Node)> getter,
                  std::vector<double> &srcDen,
                  int &nsrc) const
    {
    Eigen::Matrix<double,Nodes::DIM,N> vec_nod;
    for (int i = 0; i < N; i++)
        { vec_nod.col(i) = getter(getNode(i)); }

    double dud_sum = (vec_nod.row(IDX_X)).dot( da.col(IDX_X) )
                   + (vec_nod.row(IDX_Y)).dot( da.col(IDX_Y) )
                   + (vec_nod.row(IDX_Z)).dot( da.col(IDX_Z) );
    
    dud_sum *= -param.J/mu0;
    for(int i=0;i<Tetra::NPI;i++)
        { srcDen[nsrc+i] = weight(i)*dud_sum; }
    nsrc += Tetra::NPI;
    }

double Tet::demagEnergy(Tetra::prm const &param,
                       Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> dudx,
                       Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> dudy,
                       Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> dudz,
                       Eigen::Ref<Eigen::Matrix<double,NPI,1>> phi) const
    {
    Eigen::Matrix<double,NPI,1> dens;

    for (int npi = 0; npi < NPI; npi++)
        { dens[npi] = (dudx(0,npi) + dudy(1,npi) + dudz(2,npi)) * phi[npi]; }
    return -0.5*param.J*weight.dot(dens);
    }

double Tet::zeemanEnergy(Tetra::prm const &param, Eigen::Ref<Eigen::Vector3d> const Hext,
                        Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> const u) const
    {
    Eigen::Matrix<double,NPI,1> dens = u.transpose() * Hext;

    return -param.J*weight.dot(dens);
    }

double Tet::Jacobian(Eigen::Ref<Eigen::Matrix3d> J)
    {
    Eigen::Vector3d p0p1 = getNode(1).p - getNode(0).p;
    Eigen::Vector3d p0p2 = getNode(2).p - getNode(0).p;
    Eigen::Vector3d p0p3 = getNode(3).p - getNode(0).p;
    J(0,0) = p0p1.x();
    J(0,1) = p0p2.x();
    J(0,2) = p0p3.x();
    J(1,0) = p0p1.y();
    J(1,1) = p0p2.y();
    J(1,2) = p0p3.y();
    J(2,0) = p0p1.z();
    J(2,1) = p0p2.z();
    J(2,2) = p0p3.z();
    return J.determinant();
    }

double Tet::calc_vol(void) const
    {
    Eigen::Vector3d p0p1 = getNode(1).p - getNode(0).p;
    Eigen::Vector3d p0p2 = getNode(2).p - getNode(0).p;
    Eigen::Vector3d p0p3 = getNode(3).p - getNode(0).p;

    return p0p1.dot(p0p2.cross(p0p3))/6.0;
    }

