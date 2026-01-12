#ifndef tetra_h
#define tetra_h

/** \file tetra.h
  \brief namespace Tetra
  header containing Tet class, some constants, and integrales
 */

#include <set>

#include "facette.h"
#include "node.h"
#include "time_integration.h"
#include "element.h"

/** \namespace Tetra
 to grab altogether some constants for struct Tet
 */
namespace Tetra
    {
const int N = 4;   /**< number of sommits */
const int NPI = 5; /**< number of weights  */

const double epsilon =
        EPSILON; /**< this constant is defined from a macro in config.h.in, it is used to check the
                    validity of the tetrahedreon, a degeneracy test */

constexpr double A = 1. / 4.;                /**< constant to build hat functions */
constexpr double B = 1. / 6.;                /**< constant to build hat functions */
constexpr double C = 1. / 2.;                /**< constant to build hat functions */
constexpr double D = -2. / 15.;              /**< constant to build hat functions */
constexpr double E = 3. / 40.;               /**< constant to build hat functions */
constexpr double u[NPI] = {A, B, B, B, C};   /**< some constants to build hat functions */
constexpr double v[NPI] = {A, B, B, C, B};   /**< some constants to build hat functions */
constexpr double w[NPI] = {A, B, C, B, B};   /**< some constants to build hat functions */
constexpr double pds[NPI] = {D, E, E, E, E}; /**< some constant weights to build hat functions */

/** constant matrix \f$ j:0..4 \f$
a[0][j]   = 1.-u[j]-v[j]-w[j];
a[1][j]   = u[j];
a[2][j]   = v[j];
a[3][j]   = w[j];
*/
constexpr double a[N][NPI] = {{1. - u[0] - v[0] - w[0], 1. - u[1] - v[1] - w[1],
                               1. - u[2] - v[2] - w[2], 1. - u[3] - v[3] - w[3],
                               1. - u[4] - v[4] - w[4]},
                              {u[0], u[1], u[2], u[3], u[4]},
                              {v[0], v[1], v[2], v[3], v[4]},
                              {w[0], w[1], w[2], w[3], w[4]}};

/** eigen constant matrix a */
const Eigen::Matrix<double,N,NPI> eigen_a = (Eigen::MatrixXd(N,NPI) << a[0][0], a[0][1], a[0][2], a[0][3], a[0][4],
                                                                       a[1][0], a[1][1], a[1][2], a[1][3], a[1][4],
                                                                       a[2][0], a[2][1], a[2][2], a[2][3], a[2][4],
                                                                       a[3][0], a[3][1], a[3][2], a[3][3], a[3][4] ).finished();

/** \class prm
region number and material constants
*/
struct prm
    {
    std::string regName; /**< region name */
    double alpha_LLG;    /**< \f$ \alpha \f$ damping parameter, dimensionless */
    double A;            /**< exchange constant stiffness in [Joule per meter] = kg m^1 s^-2 */
    double Ms;           /**< Magnetization at saturation in [Ms] = A m^-1 if Ms<=0 the the region is non magnetic */
    double K;            /**< uniaxial anisotropy constant */
    Eigen::Vector3d uk;  /**< uniaxial anisotropy axis */

    double K3;           /**< cubic anisotropy constant */
    Eigen::Vector3d ex;  /**< unit vector1 (for cubic anisotropy) */
    Eigen::Vector3d ey;  /**< unit vector2 (for cubic anisotropy) */
    Eigen::Vector3d ez;  /**< unit vector3 (for cubic anisotropy) */

    double P;            /**< spin diffusion polarization rate, dimensionless */
    double N0;           /**< density of states at Fermi level in [Energy]^-1 [Volume]^-1 = kg^-1 m^-5 s^2 */
    double sigma;        /**< electrical conductivity in [Siemens per meter] = kg^-1 m^-3 s^3 A^2  */
    double lsd;          /**< diffusion length related to s-d coupling in a magnetic material m^1 */
    double lsf;          /**< spin diffusion length in a metal (magnetic or not) m^1 */
    double spinHall;     /**< Spin Orbit Torque contribution to spin diffusion due to spin Hall effect */

    double volume = 0; /**< total volume of the region */

    void infos(void);  /**< print the struct parameters */
    };

/** \class Tet
Tet is a tetrahedron, containing the index references to nodes, must not be flat <br>
indices convention is<br>
```
                        v
                      .
                    ,/
                   /
                2(ic)                                 2
              ,/|`\                                 ,/|`\
            ,/  |  `\                             ,/  |  `\
          ,/    '.   `\                         ,6    '.   `5
        ,/       |     `\                     ,/       8     `\
      ,/         |       `\                 ,/         |       `\
     0(ia)-------'.--------1(ib) --> u     0--------4--'.--------1
      `\.         |      ,/                 `\.         |      ,/
         `\.      |    ,/                      `\.      |    ,9
            `\.   '. ,/                           `7.   '. ,/
               `\. |/                                `\. |/
                  `3(id)                                `3
                     `\.
                        ` w
```
*/
class Tet : public element<N,NPI>
    {
public:
    /** constructor. It initializes weight hat function with weights \f$ w_i = |J| p_i \f$
    with  \f$ p_i = pds[i] = (D,E,E,E,E) \f$ and dad(x|y|z) if \f$ | detJ | < \epsilon \f$ jacobian
    is considered degenerated unit tests : Tet_constructor; Tet_inner_tables
    */
    inline Tet(const std::vector<Nodes::Node> &_p_node /**< vector of nodes */,
               const int _idx /**< [in] region index in region vector */,
               std::initializer_list<int> _i /**< [in] node index */)
        : element<N,NPI>(_p_node,_idx,_i), idx(0)
        {
        zeroBasing();
        da.setZero();

        if (existNodes())
            {
            orientate();

            Eigen::Matrix3d J;
            // we have to rebuild the jacobian in case of ill oriented tetrahedron
            double detJ = Jacobian(J);

            if (fabs(detJ) < Tetra::epsilon)
                {
                std::cerr << "Singular jacobian in tetrahedron: |det(J)|= " << fabs(detJ) << std::endl;
                element::infos();
                SYSTEM_ERROR;
                }
            Eigen::Matrix<double,N,Nodes::DIM> dadu;
            dadu << -1., -1., -1., 1., 0., 0., 0., 1., 0., 0., 0., 1.;
            da = dadu * J.inverse();

            for (int j = 0; j < NPI; j++)
                { weight[j] = detJ * Tetra::pds[j]; }
            }
        // do nothing lambda's (usefull for spin transfer torque)
        extraField = [] ( Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> ) {};
        extraCoeffs_BE = [](Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>>,
                            Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>>,
                            Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>>,
                            Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>>,
                            Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,N>>) {};
        }

    /** local hat functions matrix, initialized by constructor: da = dadu * inverse(Jacobian) */
    Eigen::Matrix<double,N,Nodes::DIM> da;

    /** interpolation for scalar field : the getter function is given as a parameter in order to
     * know what part of the node you want to interpolate */
    inline void interpolation(std::function<double(Nodes::Node)> getter,
                              Eigen::Ref<Eigen::Matrix<double,NPI,1>> result) const
        {
        Eigen::Matrix<double,N,1> scalar_nod;
        for (int i = 0; i < N; i++) scalar_nod(i) = getter(getNode(i));
        result = scalar_nod.transpose() * eigen_a;
        }

    /** interpolation for 3D vector field and a tensor : getter function is given as a parameter to
     * know what part of the node you want to interpolate */
    inline void interpolation(std::function<Eigen::Vector3d(Nodes::Node)> getter,
                              Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> result,
                              Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> Tx,
                              Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> Ty,
                              Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> Tz) const
        {
        Eigen::Matrix<double,Nodes::DIM,N> vec_nod;
        for (int i = 0; i < N; i++) vec_nod.col(i) = getter(getNode(i));

        result = vec_nod * eigen_a;
        Tx = vec_nod * (da.col(Nodes::IDX_X)).replicate(1,NPI);
        Ty = vec_nod * (da.col(Nodes::IDX_Y)).replicate(1,NPI);
        Tz = vec_nod * (da.col(Nodes::IDX_Z)).replicate(1,NPI);
        }

    /** interpolation for components of a field : the getter function is given as a parameter in
     * order to know what part of the node you want to interpolate */
    inline void interpolation_field(std::function<double(Nodes::Node)> getter,
                                    Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> X) const
        {
        Eigen::Matrix<double,N,1> scalar_nod;
        for (int i = 0; i < N; i++) scalar_nod(i) = getter(getNode(i));
        
        X.setZero();
        for (int j = 0; j < NPI; j++)
            {
            for (int i = 0; i < N; i++)
                {
                X.col(j) -= (scalar_nod[i] * da.row(i));
                }
            }
        }

    /** interpolation for the component idx of a field : the getter function is given as a parameter in order to
     * know what part of the node you want to interpolate */
    inline void interpolation(std::function<double(Nodes::Node, Nodes::index)> getter, Nodes::index idx,
                              Eigen::Ref<Eigen::Matrix<double,Tetra::NPI,1>> result) const
        {
        Eigen::Matrix<double,N,1> scalar_nod;

        for (int i = 0; i < N; i++)
            {
            scalar_nod[i] = getter(getNode(i),idx);
            }
        result = scalar_nod.transpose() * eigen_a;
        }
    
    /** AE matrix filling with exchange contributions, diagonal contributions and magnetization contribution
    AE has a block structure:
    ( E -Z +Y)
    ( +Z E -X)
    ( -Y +X E)
    E X Y Z are N*N matrices
    E is the sum of the exchange contribution and modified alpha (scheme stabilizer) on its diagonal 
    X Y Z are diagonal matrices with magnetization component contributions
     */
    void lumping(Eigen::Ref<Eigen::Matrix<double,NPI,1>> alpha_eff, double prefactor,
                  Eigen::Ref<Eigen::Matrix<double,3*N,3*N>> AE ) const;

    /** add drift contribution due to eventual recentering to vectors BE */
    void add_drift_BE(double alpha, double s_dt, double Vdrift,
                      Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> U,
                      Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> V, 
                      Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> dUd_,
                      Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> dVd_,
                      Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,N>> BE) const;

    /** append(+=) H_aniso for uniaxial anisotropy contribution, returns contribution to uHeff (used to
     * compute the stabilizing effective damping) */
    Eigen::Matrix<double,NPI,1> calc_aniso_uniax(Eigen::Ref<const Eigen::Vector3d> uk,
                                                 const double Kbis, const double s_dt,
                                                 Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> U,
                                                 Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> V,
                                                 Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> H_aniso) const;

    /** append(+=) H_aniso for cubic anisotropy contribution, returns contribution to uHeff (used to
     * compute the stabilizing effective damping) */
    Eigen::Matrix<double,NPI,1> calc_aniso_cub(Eigen::Ref<const Eigen::Vector3d> ex,
                                               Eigen::Ref<const Eigen::Vector3d> ey,
                                               Eigen::Ref<const Eigen::Vector3d> ez,
                                               const double K3bis, const double s_dt,
                                               Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> U,
                                               Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> V,
                                               Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> H_aniso) const;

    /** computes the integral contribution of the tetrahedron to the evolution of the magnetization
     * calc_Hext is a function that returns external H field defined on gauss points
     */
    void integrales( Tetra::prm const &param, timing const &prm_t,
                    std::function< Eigen::Matrix<double,Nodes::DIM,NPI> (void)> calc_Hext,
                    Nodes::index idx_dir, double Vdrift);

    /** exchange energy of the tetrahedron */
    double exchangeEnergy(Tetra::prm const &param,
                          Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> dudx,
                          Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> dudy,
                          Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> dudz) const;

    /** anisotropy energy of the tetrahedron */
    double anisotropyEnergy(Tetra::prm const &param, Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> u) const;

    /** computes volume charges, the divergence of the magnetization. Result is stored in srcDen, at nsrc position */
    void charges(Tetra::prm const &param,
                 std::function<Eigen::Vector3d(Nodes::Node)> getter,
                 std::vector<double> &srcDen,
                 int &nsrc) const;

    /** demagnetizing energy of the tetrahedron */
    double demagEnergy(Tetra::prm const &param,
                       Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> dudx,
                       Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> dudy,
                       Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> dudz,
                       Eigen::Ref<Eigen::Matrix<double,NPI,1>> phi) const;

    /** zeeman energy of the tetrahedron */
    double zeemanEnergy(Tetra::prm const &param, Eigen::Ref<Eigen::Vector3d> const Hext,
                        Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> const u) const;

    /** \return \f$ |J| \f$ build Jacobian \f$ J \f$ */
    double Jacobian(Eigen::Ref<Eigen::Matrix3d> J);

    /** computes volume of the tetrahedron ; unit test Tet_calc_vol */
    double calc_vol(void) const;

    /** idx is the index of the tetrahedron in the vector of tetrahedron */
    int idx;

    /** for extra contribution to the effective field, such as spin transfert torque Hm */
    std::function<void( Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> H)> extraField;

    /** for extra contribution to the matrix BE, such as spin transfer torque contribs */
    std::function<void(Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> U,
                       Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> dUdx,
                       Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> dUdy,
                       Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> dUdz,
                       Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,N>> BE)> extraCoeffs_BE;

    /** returns gauss points in result = vec_nod*Tetra::a  */
    void getPtGauss(Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> result) const
        {
        Eigen::Matrix<double,Nodes::DIM,N> vec_nod;
        for (int i = 0; i < N; i++)
            {
            vec_nod.col(i) << getNode(i).p;
            }
        result = vec_nod * eigen_a;
        }

private:
    /** orientation redefined through index swaping if needed */
    void orientate(void)
        {
        if (calc_vol() < 0.0)
                std::swap(ind[2], ind[3]);
        }

    };  // end class Tetra

    /** to perform some second order corrections, an effective \f$ \alpha \f$ is computed here with
     * a piecewise formula */
    Eigen::Matrix<double,NPI,1> calc_alpha_eff(const double dt, const double alpha,
                                               Eigen::Ref<Eigen::Matrix<double,NPI,1>> uHeff);
    }   // end namespace Tetra

#endif /* tetra_h */
