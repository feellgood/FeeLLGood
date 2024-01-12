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

void Tet::lumping(int const &npi, double alpha_eff, double prefactor,
                  Eigen::Ref<Eigen::Matrix<double,3*N,3*N>> AE ) const
    {
    const double w = weight[npi];

    for (int i = 0; i < N; i++)
        {
        const double ai_w = w * a[i][npi];
        const Eigen::Vector3d ai_w_u0 = ai_w * Nodes::get_u0(refNode[ind[i]]);

        AE(i,i) += alpha_eff * ai_w;
        AE(N + i,N + i) += alpha_eff * ai_w;
        AE(2*N + i,2*N + i) += alpha_eff * ai_w;

        AE(i, 2*N + i) += ai_w_u0(IDX_Y);
        AE(i, N + i) -= ai_w_u0(IDX_Z);
        AE(N + i, i) += ai_w_u0(IDX_Z);
        AE(N + i, 2*N + i) -= ai_w_u0(IDX_X);
        AE(2*N + i, N + i) += ai_w_u0(IDX_X);
        AE(2*N + i, i) -= ai_w_u0(IDX_Y);

        for (int j = 0; j < N; j++)
            {
            double contrib = w * prefactor
                             * (dadx(i,npi) * dadx(j,npi) + dady(i,npi) * dady(j,npi)
                                + dadz(i,npi) * dadz(j,npi));

            AE(i,j) += contrib;
            AE(N + i,N + j) += contrib;
            AE(2*N + i,2*N + j) += contrib;
            }
        }
    }

void Tet::add_drift_BE(int const &npi, double alpha, double s_dt, double Vdrift,
                       Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> U,
                       Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> V, 
                       Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> dUd_,
                       Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> dVd_,
                       Eigen::Ref<Eigen::Matrix<double,DIM,N>> BE) const
    {  // the artificial drift from eventual recentering is along x,y or z
    double w = weight[npi];
    Eigen::Matrix<double,DIM,N> interim;

    for (int i = 0; i < N; i++)
        {
        interim.col(i) = a[i][npi]*( alpha*dUd_.col(npi) + U.col(npi).cross(dUd_.col(npi)) 
                    + s_dt*(alpha*dVd_.col(npi) + U.col(npi).cross(dVd_.col(npi)) + V.col(npi).cross(dUd_.col(npi))) );
        }
    BE += w*Vdrift*interim;
    }

double Tet::calc_aniso_uniax(const int npi, Eigen::Ref<const Eigen::Vector3d> uk, const double Kbis,
                             const double s_dt,
                             Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> U,
                             Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> V,
                             Eigen::Ref<Eigen::Vector3d> H_aniso) const
    {
    H_aniso += (Kbis * uk.dot( U.col(npi) + s_dt * V.col(npi))) * uk;
    return Kbis * sq( uk.dot(U.col(npi)));
    }

double Tet::calc_aniso_cub(const int npi,
                           Eigen::Ref<const Eigen::Vector3d> ex,
                           Eigen::Ref<const Eigen::Vector3d> ey,
                           Eigen::Ref<const Eigen::Vector3d> ez,
                           const double K3bis, const double s_dt,
                           Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> U,
                           Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> V,
                           Eigen::Ref<Eigen::Vector3d> H_aniso) const
    {
    Eigen::Vector3d uk_u = Eigen::Vector3d(ex.dot(U.col(npi)), ey.dot(U.col(npi)), ez.dot(U.col(npi)));
    Eigen::Vector3d uk_v = Eigen::Vector3d(ex.dot(V.col(npi)), ey.dot(V.col(npi)), ez.dot(V.col(npi)));
    
    Eigen::Vector3d uk_uuu = uk_u - uk_u.unaryExpr( [](double x){ return x*x*x;} ); // should be replaced by cube()
    Eigen::Vector3d tmp = uk_v.cwiseProduct(ex);

    H_aniso += -K3bis * (uk_uuu(0) * ex + uk_uuu(1) * ey + uk_uuu(2) * ez
                  + s_dt * tmp.cwiseProduct( Eigen::Vector3d(1, 1, 1) - 3*uk_u.cwiseProduct(uk_u) ));
    return -K3bis *uk_u.dot(uk_uuu);
    }

void Tet::integrales(Tetra::prm const &param, timing const &prm_t,
                     Eigen::Vector3d const &Hext, Nodes::index idx_dir, double Vdrift)
    {
    double alpha = param.alpha_LLG;
    double Js = param.J;
    double Abis = 2.0 * param.A / Js;
    double Kbis = 2.0 * param.K / Js;
    double K3bis = 2.0 * param.K3 / Js;
    const double s_dt = THETA * prm_t.get_dt() * gamma0;  // theta from theta scheme in config.h.in

    Eigen::Matrix<double,3*N,3*N> AE;
    AE.setZero();
    //Eigen::Matrix<double,3*N,1> BE;
    Eigen::Matrix<double,DIM,N> BE;
    BE.setZero();
    
    /*-------------------- INTERPOLATION --------------------*/
    Eigen::Matrix<double,DIM,NPI> U,dUdx,dUdy,dUdz;
    interpolation(Nodes::get_u0, U, dUdx, dUdy, dUdz);
    Eigen::Matrix<double,DIM,NPI> V,dVdx,dVdy,dVdz;
    interpolation(Nodes::get_v0, V, dVdx, dVdy, dVdz);
    Eigen::Matrix<double,DIM,NPI> Hd;
    interpolation_field(Nodes::get_phi0, Hd);
    Eigen::Matrix<double,DIM,NPI> Hv;
    interpolation_field(Nodes::get_phiv0, Hv);
    /*-------------------- END INTERPOLATION ----------------*/

    for (int npi = 0; npi < NPI; npi++)
        {
        Eigen::Vector3d H_aniso;
        H_aniso.setZero();
        
        double contrib_aniso = calc_aniso_uniax(npi, param.uk, Kbis, s_dt, U, V, H_aniso);
        contrib_aniso += calc_aniso_cub(npi, param.ex, param.ey, param.ez, K3bis, s_dt, U, V, H_aniso);

        Eigen::Vector3d H = H_aniso + Hd.col(npi) + Hext + (s_dt / gamma0) * Hv.col(npi);
        const double w = weight[npi];
        for (int i = 0; i < N; i++)
            {
            BE.col(i) += w*(-Abis*(da(i,0)*dUdx.col(npi) + da(i,1)*dUdy.col(npi) + da(i,2)*dUdz.col(npi)) + a[i][npi]*H);
            }

        extraCoeffs_BE(npi, Js, U.col(npi), dUdx.col(npi), dUdy.col(npi), dUdz.col(npi), BE);  // STT
        if (idx_dir != IDX_UNDEF)
            {
            if (idx_dir == IDX_Z)
                add_drift_BE(npi, alpha, s_dt, Vdrift, U, V, dUdz, dVdz, BE);
            else if (idx_dir == IDX_Y)
                add_drift_BE(npi, alpha, s_dt, Vdrift, U, V, dUdy, dVdy, BE);
            else if (idx_dir == IDX_X)
                add_drift_BE(npi, alpha, s_dt, Vdrift, U, V, dUdx, dVdx, BE);
            }

        H = Hext + Hd.col(npi);
        Eigen::Vector3d tmp_H = extraField(npi);  // computes STT contrib
        H += tmp_H;
        
        double uHeff =
                contrib_aniso - Abis *( dUdx.col(npi).squaredNorm() + dUdy.col(npi).squaredNorm() + dUdz.col(npi).squaredNorm());
        uHeff += H.dot( U.col(npi) );

        lumping(npi, prm_t.calc_alpha_eff(alpha, uHeff), prm_t.prefactor * s_dt * Abis, AE);
        }
    /*-------------------- PROJECTIONS --------------------*/
    Kp = P*AE*P.transpose();// with MKL installed this operation should call dgemm_direct

    // the following deep copy might be replaced by reshaped(), but it is in eigen since version 3.4, not 3.3
    Eigen::Matrix<double,DIM*N,1> tmp;
    for(int k=0;k<DIM;k++)
        for(int i=0;i<N;i++)
            tmp(k*N+i) = BE(k,i);
    Lp = P*tmp;
    }

double Tet::exchangeEnergy(Tetra::prm const &param,
                           Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> dudx,
                           Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> dudy,
                           Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> dudz) const
    {
    Eigen::Matrix<double,NPI,1> dens;

    for (int npi = 0; npi < NPI; npi++)
        {
        dens[npi] = dudx.col(npi).squaredNorm()
                  + dudy.col(npi).squaredNorm()
                  + dudz.col(npi).squaredNorm();
        }
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

Eigen::Matrix<double,NPI,1> Tet::charges(std::function<Eigen::Vector3d(Nodes::Node)> getter) const
    {
    Eigen::Matrix<double,DIM,NPI> dudx,dudy,dudz;
    interpolation(getter,dudx,dudy,dudz);
    
    Eigen::Matrix<double,NPI,1> result;
    for (int j = 0; j < NPI; j++)
        { result(j) = -Ms * (dudx(0,j) + dudy(1,j) + dudz(2,j)) * weight[j]; }
    return result;
    }

double Tet::demagEnergy(Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> dudx,
                       Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> dudy,
                       Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> dudz,
                       Eigen::Ref<Eigen::Matrix<double,NPI,1>> phi) const
    {
    Eigen::Matrix<double,NPI,1> dens;

    for (int npi = 0; npi < NPI; npi++)
        { dens[npi] = (dudx(0,npi) + dudy(1,npi) + dudz(2,npi)) * phi[npi]; }
    return -0.5*mu0*Ms*weight.dot(dens);
    }

double Tet::zeemanEnergy(Tetra::prm const &param, double uz_drift, Eigen::Ref<Eigen::Vector3d> const Hext,
                        Eigen::Ref<Eigen::Matrix<double,DIM,NPI>> const u) const
    {
    Eigen::Matrix<double,NPI,1> dens;

    for (int npi = 0; npi < NPI; npi++)
        { dens[npi] = u.col(npi).dot(Hext) + uz_drift * Hext.z(); }
    return -param.J*weight.dot(dens);
    }

double Tet::Jacobian(Eigen::Ref<Eigen::Matrix3d> J)
    {
    Eigen::Vector3d p0p1 = refNode[ind[1]].p - refNode[ind[0]].p;
    Eigen::Vector3d p0p2 = refNode[ind[2]].p - refNode[ind[0]].p;
    Eigen::Vector3d p0p3 = refNode[ind[3]].p - refNode[ind[0]].p;
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
    Eigen::Vector3d p0p1 = refNode[ind[1]].p - refNode[ind[0]].p;
    Eigen::Vector3d p0p2 = refNode[ind[2]].p - refNode[ind[0]].p;
    Eigen::Vector3d p0p3 = refNode[ind[3]].p - refNode[ind[0]].p;

    return p0p1.dot(p0p2.cross(p0p3))/6.0;
    }

std::set<Facette::Fac> Tet::ownedFac() const
    {
    std::set<Facette::Fac> s;

    const int ia = ind[0];
    const int ib = ind[1];
    const int ic = ind[2];
    const int id = ind[3];

    s.insert(Facette::Fac(refNode, 0, idxPrm, {ia, ic, ib} ));
    s.insert(Facette::Fac(refNode, 0, idxPrm, {ib, ic, id} ));
    s.insert(Facette::Fac(refNode, 0, idxPrm, {ia, id, ic} ));
    s.insert(Facette::Fac(refNode, 0, idxPrm, {ia, ib, id} ));

    return s;
    }
