/**
  Elementary matrix Calculation for a tetrahedron element
 */

#include <set>

#include "config.h"  // to get gamma0 constant
#include "pt3D.h"
#include "tetra.h"
#include "time_integration.h"
#include "facette.h"

using namespace Tetra;
using namespace Pt;

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

        AE(i, 2*N + i) += ai_w_u0(Pt::IDX_Y);
        AE(i, N + i) -= ai_w_u0(Pt::IDX_Z);
        AE(N + i, i) += ai_w_u0(Pt::IDX_Z);
        AE(N + i, 2*N + i) -= ai_w_u0(Pt::IDX_X);
        AE(2*N + i, N + i) += ai_w_u0(Pt::IDX_X);
        AE(2*N + i, i) -= ai_w_u0(Pt::IDX_Y);

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
                       Eigen::Ref<Eigen::Matrix<double,Pt::DIM,NPI>> U,
                       Eigen::Ref<Eigen::Matrix<double,Pt::DIM,NPI>> V, 
                       Eigen::Ref<Eigen::Matrix<double,Pt::DIM,NPI>> dUd_,
                       Eigen::Ref<Eigen::Matrix<double,Pt::DIM,NPI>> dVd_,
                       Eigen::Ref<Eigen::Vector<double,3*N>> BE) const
    {  // the artificial drift from eventual recentering is along x,y or z
    double w = weight[npi];
    Eigen::Vector3d interim[N];

    for (int i = 0; i < N; i++)
        {
        interim[i] = a[i][npi]*( alpha*dUd_.col(npi) + U.col(npi).cross(dUd_.col(npi)) 
                    + s_dt*(alpha*dVd_.col(npi) + U.col(npi).cross(dVd_.col(npi)) + V.col(npi).cross(dUd_.col(npi))) );
        }
    
    for (int i = 0; i < N; i++)
        for(int k=0;k<Pt::DIM;k++)
            { BE(k*N + i) += w*Vdrift*interim[i](k); }
    }

double Tet::calc_aniso_uniax(int const &npi, Eigen::Ref<Eigen::Vector3d> const uk, const double Kbis,
                             const double s_dt,
                             Eigen::Ref<Eigen::Matrix<double,Pt::DIM,NPI>> U,
                             Eigen::Ref<Eigen::Matrix<double,Pt::DIM,NPI>> V,
                             Eigen::Ref<Eigen::Vector3d> H_aniso) const
    {
    H_aniso += (Kbis * uk.dot( U.col(npi) + s_dt * V.col(npi))) * uk;
    return Kbis * sq( uk.dot(U.col(npi)));
    }

double Tet::calc_aniso_cub(int const &npi,
                           Eigen::Ref<Eigen::Vector3d> ex,
                           Eigen::Ref<Eigen::Vector3d> ey,
                           Eigen::Ref<Eigen::Vector3d> ez, const double K3bis, const double s_dt,
                           Eigen::Ref<Eigen::Matrix<double,Pt::DIM,NPI>> U,
                           Eigen::Ref<Eigen::Matrix<double,Pt::DIM,NPI>> V,
                           Eigen::Ref<Eigen::Vector3d> H_aniso) const
    {
    Eigen::Vector3d uk_u = Eigen::Vector3d(ex.dot(U.col(npi)), ey.dot(U.col(npi)), ez.dot(U.col(npi)));
    Eigen::Vector3d uk_v = Eigen::Vector3d(ex.dot(V.col(npi)), ey.dot(V.col(npi)), ez.dot(V.col(npi)));
    
    Eigen::Vector3d uk_uuu = uk_u - uk_u.unaryExpr( [](double x){ return x*x*x;} );
    Eigen::Vector3d tmp = uk_v.cwiseProduct(ex);

    H_aniso += -K3bis * (uk_uuu(0) * ex + uk_uuu(1) * ey + uk_uuu(2) * ez
                  + s_dt * tmp.cwiseProduct( Eigen::Vector3d(1, 1, 1) - 3*uk_u.cwiseProduct(uk_u) ));
    return -K3bis *uk_u.dot(uk_uuu);
    }

void Tet::integrales(Tetra::prm const &param, timing const &prm_t,
                     Eigen::Vector3d const &Hext, Pt::index idx_dir, double Vdrift)
    {
    double alpha = param.alpha_LLG;
    double Js = param.J;
    double Abis = 2.0 * param.A / Js;
    double Kbis = 2.0 * param.K / Js;
    double K3bis = 2.0 * param.K3 / Js;
    const double s_dt = THETA * prm_t.get_dt() * gamma0;  // theta from theta scheme in config.h.in

    Eigen::Matrix<double,3*N,3*N> AE;
    AE.setZero();
    Eigen::Vector<double,3*N> BE;
    BE.setZero();
    
    /*-------------------- INTERPOLATION --------------------*/
    
    //interpolation(Nodes::get_u0, U, dUdx, dUdy, dUdz);
    Eigen::Matrix<double,Pt::DIM,N> vec_nod;
    for (int i = 0; i < N; i++)
            vec_nod.col(i) = Nodes::get_u0(refNode[ind[i]]);
    //getDataFromNode<Eigen::Vector3d>(Nodes::get_u0, vec_nod);
    Eigen::Matrix<double,Pt::DIM,NPI> U = vec_nod * eigen_a;
    Eigen::Matrix<double,Pt::DIM,NPI> dUdx = vec_nod * dadx;
    Eigen::Matrix<double,Pt::DIM,NPI> dUdy = vec_nod * dady;
    Eigen::Matrix<double,Pt::DIM,NPI> dUdz = vec_nod * dadz;
    
    //interpolation(Nodes::get_v0, V, dVdx, dVdy, dVdz);
    for (int i = 0; i < N; i++)
            vec_nod.col(i) = Nodes::get_v0(refNode[ind[i]]);
    //getDataFromNode<Eigen::Vector3d>(Nodes::get_v0, vec_nod);
    Eigen::Matrix<double,Pt::DIM,NPI> V = vec_nod * eigen_a;
    Eigen::Matrix<double,Pt::DIM,NPI> dVdx = vec_nod * dadx;
    Eigen::Matrix<double,Pt::DIM,NPI> dVdy = vec_nod * dady;
    Eigen::Matrix<double,Pt::DIM,NPI> dVdz = vec_nod * dadz;
    

    Eigen::Matrix<double,Pt::DIM,NPI> Hd;
    //interpolation(Nodes::get_phi0, Hd);    
    
    Eigen::Vector<double,N> scalar_nod;
    for (int i = 0; i < N; i++)
        scalar_nod[i] = Nodes::get_phi0(refNode[ind[i]]);
    // same as tiny::neg_transposed_mult
    for (int j = 0; j < NPI; j++)
        {
        Hd.col(j).setZero();
            for (int i = 0; i < N; i++)
            { Hd.col(j) -= scalar_nod[i] * Eigen::Vector3d(dadx(i,j), dady(i,j), dadz(i,j)); }
        }
    
    Eigen::Matrix<double,Pt::DIM,NPI> Hv;
    //interpolation(Nodes::get_phiv0, Hv);
    for (int i = 0; i < N; i++)
        scalar_nod[i] = Nodes::get_phiv0(refNode[ind[i]]);
    
    for (int j = 0; j < NPI; j++)
        {
        Hv.col(j).setZero();
            for (int i = 0; i < N; i++)
                { Hv.col(j) -= scalar_nod[i] * Eigen::Vector3d(dadx(i,j), dady(i,j), dadz(i,j)); }
        }
/*-------------------- END INTERPOLATION --------------------*/

    for (int npi = 0; npi < NPI; npi++)
        {
        Eigen::Vector3d H_aniso;
        H_aniso.setZero();
        
        //devNote: for some obscure reason we have to deep copy uk and ex,ey,ez for calc_aniso_xxx functions ??
        Eigen::Vector3d uk { param.uk.x(),param.uk.y() ,param.uk.z() };
        double contrib_aniso = calc_aniso_uniax(npi, uk, Kbis, s_dt, U, V, H_aniso);
        
        Eigen::Vector3d ex { param.ex.x(), param.ex.y(), param.ex.z() };
        Eigen::Vector3d ey { param.ey.x(), param.ey.y(), param.ey.z() };
        Eigen::Vector3d ez { param.ez.x(), param.ez.y(), param.ez.z() };
        contrib_aniso += calc_aniso_cub(npi, ex, ey, ez, K3bis, s_dt, U, V, H_aniso);

        Eigen::Vector3d H = H_aniso + Hd.col(npi) + Hext + (s_dt / gamma0) * Hv.col(npi);
        const double w = weight[npi];
        for (int i = 0; i < N; i++)
            {
            const Eigen::Vector3d interim = -Abis*(da(i,0)*dUdx.col(npi) + da(i,1)*dUdy.col(npi) + da(i,2)*dUdz.col(npi)) + a[i][npi]*H;
            for(int k=0;k<Pt::DIM;k++)
                { BE(k*N + i) += w*interim(k); }
            }

        extraCoeffs_BE(npi, Js, U.col(npi), dUdx.col(npi), dUdy.col(npi), dUdz.col(npi), BE);  // STT
        if (idx_dir != Pt::IDX_UNDEF)
            {
            if (idx_dir == Pt::IDX_Z)
                add_drift_BE(npi, alpha, s_dt, Vdrift, U, V, dUdz, dVdz, BE);
            else if (idx_dir == Pt::IDX_Y)
                add_drift_BE(npi, alpha, s_dt, Vdrift, U, V, dUdy, dVdy, BE);
            else if (idx_dir == Pt::IDX_X)
                add_drift_BE(npi, alpha, s_dt, Vdrift, U, V, dUdx, dVdx, BE);
            }

        H = Hext + Hd.col(npi);
        Pt::pt3D tmp_H(0,0,0);
        extraField(npi, tmp_H);  // add STT contribution to the effective field H through tmp_H
        H += Eigen::Vector3d(tmp_H.x(),tmp_H.y(),tmp_H.z());
        
        double uHeff =
                contrib_aniso - Abis *( dUdx.col(npi).squaredNorm() + dUdy.col(npi).squaredNorm() + dUdz.col(npi).squaredNorm());
        uHeff += H.dot( U.col(npi) );

        lumping(npi, prm_t.calc_alpha_eff(alpha, uHeff), prm_t.prefactor * s_dt * Abis, AE);
        }
    /*-------------------- PROJECTIONS --------------------*/
    Kp = P*AE*P.transpose();// with MKL installed this operation should call dgemm_direct
    Lp = P*BE;
    }

double Tet::exchangeEnergy(Tetra::prm const &param,
                           Eigen::Ref<Eigen::Matrix<double,Pt::DIM,NPI>> dudx,
                           Eigen::Ref<Eigen::Matrix<double,Pt::DIM,NPI>> dudy,
                           Eigen::Ref<Eigen::Matrix<double,Pt::DIM,NPI>> dudz) const
    {
    double dens[NPI];

    for (int npi = 0; npi < NPI; npi++)
        {
        dens[npi] = sq(dudx(0,npi)) + sq(dudy(0,npi)) + sq(dudz(0,npi))
                  + sq(dudx(1,npi)) + sq(dudy(1,npi)) + sq(dudz(1,npi))
                  + sq(dudx(2,npi)) + sq(dudy(2,npi)) + sq(dudz(2,npi));
        }
    return (param.A * weightedScalarProd(dens));
    }

double Tet::anisotropyEnergy(Tetra::prm const &param, Eigen::Ref<Eigen::Matrix<double,Pt::DIM,NPI>> u) const
    {
    double dens[NPI];

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

    return weightedScalarProd(dens);
    }

Eigen::Vector<double,NPI> Tet::charges(std::function<Eigen::Vector3d(Nodes::Node)> getter) const
    {
    //double dudx[DIM][NPI], dudy[DIM][NPI], dudz[DIM][NPI];

    //Eigen::Vector3d vec_nod[N];
    //getDataFromNode<Eigen::Vector3d>(getter, vec_nod);
    Eigen::Matrix<double,DIM,N> vec_nod;
    for (int i = 0; i < N; i++)
            vec_nod.col(i) = getter(refNode[ind[i]]);

    Eigen::Matrix<double,DIM,NPI> dudx = vec_nod * dadx;//tiny::mult<double, N, NPI>(vec_nod, dadx, dudx);
    Eigen::Matrix<double,DIM,NPI> dudy = vec_nod * dady;//tiny::mult<double, N, NPI>(vec_nod, dady, dudy);
    Eigen::Matrix<double,DIM,NPI> dudz = vec_nod * dadz;//tiny::mult<double, N, NPI>(vec_nod, dadz, dudz);

    Eigen::Vector<double,NPI> result;
    for (int j = 0; j < NPI; j++)
        {
        result(j) = -Ms * (dudx(0,j) + dudy(1,j) + dudz(2,j)) * weight[j];
        }
    return result;
    }

double Tet::demagEnergy(Eigen::Ref<Eigen::Matrix<double,Pt::DIM,NPI>> dudx,
                       Eigen::Ref<Eigen::Matrix<double,Pt::DIM,NPI>> dudy,
                       Eigen::Ref<Eigen::Matrix<double,Pt::DIM,NPI>> dudz,
                       Eigen::Ref<Eigen::Vector<double,NPI>> phi) const
    {
    double dens[NPI];
    //double Ms = nu0 * param.J;

    for (int npi = 0; npi < NPI; npi++)
        {
        dens[npi] = (dudx(0,npi) + dudy(1,npi) + dudz(2,npi)) * phi[npi];
        }
    return (-0.5 * mu0 * Ms * weightedScalarProd(dens));
    }

double Tet::zeemanEnergy(Tetra::prm const &param, double uz_drift, Eigen::Ref<Eigen::Vector3d> const Hext,
                        Eigen::Ref<Eigen::Matrix<double,Pt::DIM,NPI>> const u) const
    {
    double dens[NPI];

    for (int npi = 0; npi < NPI; npi++)
        {
        //Pt::pt3D u_npi = Pt::pt3D(u[0][npi], u[1][npi], u[2][npi]);
        dens[npi] = u.col(npi).dot(Hext) + uz_drift * Hext.z();
        }
    return (-param.J * weightedScalarProd(dens));
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

    return p0p1.dot(p0p2.cross(p0p3))/6.0; //Pt::pTriple(p1 - p0, p2 - p0, p3 - p0) / 6.0;
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
