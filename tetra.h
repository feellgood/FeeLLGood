#ifndef tetra_h
#define tetra_h

/** \file tetra.h
  \brief namespace Tetra
  header containing Tet class, some constants, and integrales
 */

#include <set>

#include "spinTransferTorque.h"
#include "facette.h"
#include "node.h"
#include "time_integration.h"
#include "tiny.h"

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

/** \class prm
region number and material constants
*/
struct prm
    {
    std::string regName; /**< region name */
    double alpha_LLG;    /**< \f$ \alpha \f$ damping parameter */
    double A;            /**< exchange constant stiffness */
    double J;            /**< \f$ M_s = \nu_0 J \f$ */
    double K;            /**< uniaxial anisotropy constant */
    Pt::pt3D uk;         /**< uniaxial anisotropy axis */

    double K3;   /**< cubic anisotropy constant */
    Pt::pt3D ex; /**< unit vector1 (for cubic anisotropy) */
    Pt::pt3D ey; /**< unit vector2 (for cubic anisotropy) */
    Pt::pt3D ez; /**< unit vector3 (for cubic anisotropy) */
    STT p_STT;   /**< spin transfert torque (thiaville STT) parameters */

    /** print the struct parameters */
    inline void infos()
        {
        std::cout << "volume region name = " << regName << std::endl;
        std::cout << "alpha_LLG = " << alpha_LLG << std::endl;
        std::cout << "A = " << A << std::endl;
        std::cout << "J = " << J << std::endl;

        if (K != 0)
            {
            std::cout << "K*uk =" << K << "*[ " << uk << "]" << std::endl;
            }

        if (K3 != 0)
            {
            std::cout << "K3 = " << K3 << "; ex=[ " << ex << "], ey=[" << ey << "], ez=[" << ez
                      << "]\n";
            }
        };
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

        if (refNode.size() > 0)
            {
            orientate();

            Eigen::Matrix3d J;
            // we have to rebuild the jacobian in case of ill oriented tetrahedron
            double detJ = Jacobian(J);

            if (fabs(detJ) < Tetra::epsilon)
                {
                std::cerr << "Singular jacobian in tetrahedron" << std::endl;
                element::infos();
                SYSTEM_ERROR;
                }
            Eigen::Matrix<double,N,Pt::DIM> dadu;
            dadu << -1., -1., -1., 1., 0., 0., 0., 1., 0., 0., 0., 1.;
            da = dadu * J.inverse();

            for (int j = 0; j < NPI; j++)
                {
                for (int i = 0; i < N; i++)
                    {
                    dadx[i][j] = da(i,0);
                    dady[i][j] = da(i,1);
                    dadz[i][j] = da(i,2);
                    }
                weight[j] = detJ * Tetra::pds[j];
                }
            }
        // do nothing lambda's (usefull for spin transfer torque)
        extraField = [](int, Pt::pt3D &) {};
        extraCoeffs_BE = [](int, double, Pt::pt3D &, Pt::pt3D &, Pt::pt3D &, Pt::pt3D &,
                            Eigen::Ref<Eigen::Vector<double,3*N>>) {};
        }

    /** variations of hat function along x directions */
    double dadx[N][NPI];

    /** variations of hat function along y directions */
    double dady[N][NPI];

    /** variations of hat function along z directions */
    double dadz[N][NPI];

    /** weighted scalar product */
    inline double weightedScalarProd(const double (&X)[NPI]) const
        {
        return (X[0] * weight[0] + X[1] * weight[1] + X[2] * weight[2] + X[3] * weight[3]
                + X[4] * weight[4]);
        }

    /** interpolation for scalar field : the getter function is given as a parameter in order to
     * know what part of the node you want to interpolate */
    inline void interpolation(std::function<double(Nodes::Node)> getter,
                              double (&result)[NPI]) const
        {
        double scalar_nod[N];
        getDataFromNode<double>(getter, scalar_nod);
        tiny::transposed_mult<double, N, NPI>(scalar_nod, a, result);
        }

    /** interpolation for a tensor : the getter function is given as a parameter in order to know
     * what part of the node you want to interpolate */
/*
    inline void interpolation(std::function<Pt::pt3D(Nodes::Node)> getter,
                              double (&Tx)[Pt::DIM][NPI], double (&Ty)[Pt::DIM][NPI],
                              double (&Tz)[Pt::DIM][NPI]) const
        {
        Pt::pt3D vec_nod[N];
        getDataFromNode<Pt::pt3D>(getter, vec_nod);

        tiny::mult<double, N, NPI>(vec_nod, dadx, Tx);
        tiny::mult<double, N, NPI>(vec_nod, dady, Ty);
        tiny::mult<double, N, NPI>(vec_nod, dadz, Tz);
        }
*/

    /** interpolation for 3D vector field and a tensor : getter function is given as a parameter to
     * know what part of the node you want to interpolate */
    inline void interpolation(std::function<Pt::pt3D(Nodes::Node)> getter, Pt::pt3D (&result)[NPI],
                              Pt::pt3D (&Tx)[NPI], Pt::pt3D (&Ty)[NPI], Pt::pt3D (&Tz)[NPI]) const
        {
        double u[Pt::DIM][NPI];
        double dudx[Pt::DIM][NPI], dudy[Pt::DIM][NPI], dudz[Pt::DIM][NPI];

        Pt::pt3D vec_nod[N];
        getDataFromNode<Pt::pt3D>(getter, vec_nod);

        tiny::mult<double, N, NPI>(vec_nod, a, u);
        tiny::mult<double, N, NPI>(vec_nod, dadx, dudx);
        tiny::mult<double, N, NPI>(vec_nod, dady, dudy);
        tiny::mult<double, N, NPI>(vec_nod, dadz, dudz);

        for (int npi = 0; npi < NPI; npi++)
            {
            result[npi] = Pt::pt3D(u[0][npi], u[1][npi], u[2][npi]);
            Tx[npi] = Pt::pt3D(dudx[0][npi], dudx[1][npi], dudx[2][npi]);
            Ty[npi] = Pt::pt3D(dudy[0][npi], dudy[1][npi], dudy[2][npi]);
            Tz[npi] = Pt::pt3D(dudz[0][npi], dudz[1][npi], dudz[2][npi]);
            }  // copie qui pourrait être évitée : à améliorer
        }

    /** interpolation for 3D vector field and a tensor : the getter function is given as a parameter
     * in order to know what part of the node you want to interpolate */
    inline void interpolation(std::function<Pt::pt3D(Nodes::Node)> getter,
                              double (&result)[Pt::DIM][NPI], double (&Tx)[Pt::DIM][NPI],
                              double (&Ty)[Pt::DIM][NPI], double (&Tz)[Pt::DIM][NPI]) const
        {
        Pt::pt3D vec_nod[N];
        getDataFromNode<Pt::pt3D>(getter, vec_nod);

        tiny::mult<double, N, NPI>(vec_nod, a, result);
        tiny::mult<double, N, NPI>(vec_nod, dadx, Tx);
        tiny::mult<double, N, NPI>(vec_nod, dady, Ty);
        tiny::mult<double, N, NPI>(vec_nod, dadz, Tz);
        }

    /** interpolation for components of a field : the getter function is given as a parameter in
     * order to know what part of the node you want to interpolate */
    inline void interpolation(std::function<double(Nodes::Node)> getter, Pt::pt3D (&X)[NPI]) const
        {
        double scalar_nod[N];
        getDataFromNode<double>(getter, scalar_nod);

        // same as tiny::neg_transposed_mult
        for (int j = 0; j < NPI; j++)
            {
            X[j] = Pt::pt3D(0, 0, 0);
            for (int i = 0; i < N; i++)
                {
                X[j] -= (scalar_nod[i] * Pt::pt3D(dadx[i][j], dady[i][j], dadz[i][j]));
                }
            }
        }

    /** interpolation for components of a field : the getter function is given as a parameter in
     * order to know what part of the node you want to interpolate */
    inline void interpolation(std::function<double(Nodes::Node)> getter, double (&Xx)[NPI],
                              double (&Xy)[NPI], double (&Xz)[NPI]) const
        {
        double scalar_nod[N];
        getDataFromNode<double>(getter, scalar_nod);
        tiny::neg_transposed_mult<double, N, NPI>(scalar_nod, dadx, Xx);
        tiny::neg_transposed_mult<double, N, NPI>(scalar_nod, dady, Xy);
        tiny::neg_transposed_mult<double, N, NPI>(scalar_nod, dadz, Xz);
        }

    /** interpolation for scalar field : the getter function is given as a parameter in order to
     * know what part of the node you want to interpolate */
    inline void interpolation(std::function<double(Nodes::Node, Pt::index)> getter, Pt::index idx,
                              double (&result)[NPI]) const
        {
        double scalar_nod[N];

        for (int i = 0; i < N; i++)
            {
            scalar_nod[i] = getter(refNode[ind[i]], idx);
            }
        tiny::transposed_mult<double, N, NPI>(scalar_nod, a, result);
        }

    /** AE matrix filling */
    void lumping(int const &npi, double alpha_eff, double prefactor,
                 Eigen::Ref<Eigen::Matrix<double,3*N,3*N>> AE ) const;

    /** add drift contribution due to eventual recentering to vectors BE */
    void add_drift_BE(int const &npi, double alpha, double s_dt, double Vdrift, Pt::pt3D (&U)[NPI],
                      Pt::pt3D (&V)[NPI], Pt::pt3D (&dUd_)[NPI], Pt::pt3D (&dVd_)[NPI],
                      Eigen::Ref<Eigen::Vector<double,3*N>> BE) const;

    /** append H_aniso for uniaxial anisotropy contribution, returns contribution to uHeff (used to
     * compute the stabilizing effective damping) */
    double calc_aniso_uniax(int const &npi, Pt::pt3D const &uk, const double Kbis,
                            const double s_dt, Pt::pt3D (&U)[NPI], Pt::pt3D (&V)[NPI],
                            Pt::pt3D &H_aniso) const;

    /** append H_aniso for cubic anisotropy contribution, returns contribution to uHeff (used to
     * compute the stabilizing effective damping) */
    double calc_aniso_cub(int const &npi, Pt::pt3D const &ex, Pt::pt3D const &ey,
                          Pt::pt3D const &ez, const double K3bis, const double s_dt,
                          Pt::pt3D (&U)[NPI], Pt::pt3D (&V)[NPI], Pt::pt3D &H_aniso) const;

    /** computes the integral contribution of the tetrahedron to the evolution of the magnetization
     */
    void integrales(std::vector<Tetra::prm> const &params, timing const &prm_t,
                    Pt::pt3D const &Hext, Pt::index idx_dir, double Vdrift);

    /** exchange energy of the tetrahedron */
    double exchangeEnergy(Tetra::prm const &param, const double (&dudx)[Pt::DIM][NPI],
                          const double (&dudy)[Pt::DIM][NPI],
                          const double (&dudz)[Pt::DIM][NPI]) const;

    /** anisotropy energy of the tetrahedron */
    double anisotropyEnergy(Tetra::prm const &param, const double (&u)[Pt::DIM][NPI]) const;

    /** return volume charges  */
    Eigen::Vector<double,NPI> charges(std::function<Pt::pt3D(Nodes::Node)> getter) const;

    /** demagnetizing energy of the tetrahedron */
    double demagEnergy(const double (&dudx)[Pt::DIM][NPI],
                       const double (&dudy)[Pt::DIM][NPI], const double (&dudz)[Pt::DIM][NPI],
                       const double (&phi)[NPI]) const;

    /** zeeman energy of the tetrahedron */
    double zeemanEnergy(Tetra::prm const &param, double uz_drift, Pt::pt3D const &Hext,
                        const double (&u)[Pt::DIM][NPI]) const;

    /** \return \f$ |J| \f$ build Jacobian \f$ J \f$ */
    double Jacobian(Eigen::Ref<Eigen::Matrix3d> J);

    /** computes volume of the tetrahedron ; unit test Tet_calc_vol */
    double calc_vol(void) const;

    /** return a set of the four facettes of the tetrahedron */
    std::set<Facette::Fac> ownedFac() const;

    /** idx is the index of the tetrahedron in the vector of tetrahedron */
    int idx;

    /** for extra contribution to the effective field, such as spin transfert torque Hm */
    std::function<void(int npi, Pt::pt3D &Hm)> extraField;

    /** for extra contribution to the matrix BE, such as spin transfer torque contribs */
    std::function<void(int npi, double Js, Pt::pt3D &U, Pt::pt3D &dUdx, Pt::pt3D &dUdy,
                       Pt::pt3D &dUdz, Eigen::Ref<Eigen::Vector<double,3*N>> BE)> extraCoeffs_BE;

    /** returns gauss points in result = vec_nod*Tetra::a  */
    void getPtGauss(Eigen::Ref<Eigen::Matrix<double,Pt::DIM,NPI>> result) const
        {
        static const Eigen::Matrix<double,N,NPI> eigen_a = [] {
        Eigen::Matrix<double,N,NPI> tmp; tmp << a[0][0], a[0][1], a[0][2], a[0][3], a[0][4],
                                            a[1][0], a[1][1], a[1][2], a[1][3], a[1][4],
                                            a[2][0], a[2][1], a[2][2], a[2][3], a[2][4],
                                            a[3][0], a[3][1], a[3][2], a[3][3], a[3][4];
                                            return tmp; }();
        Eigen::Matrix<double,Pt::DIM,N> vec_nod;
        for (int i = 0; i < N; i++)
            {
            //const Pt::pt3D & tmp = refNode[ind[i]].p;
            vec_nod.col(i) << refNode[ind[i]].p;//tmp.x(), tmp.y(), tmp.z();
            }
        result = vec_nod * eigen_a;
        }

private:
    /** local hat functions matrix, initialized by constructor: da = dadu * inverse(Jacobian) */
    Eigen::Matrix<double,N,Pt::DIM> da;

    void orientate(void)
        {
        if (calc_vol() < 0.0)
                std::swap(ind[2], ind[3]);
        }

    /** template getter to access and copy parts of the node vector of type T= double | Pt::pt3D */
    template<class T>
    void getDataFromNode(std::function<T(Nodes::Node)> getter, T (&X_data)[N]) const
        {
        for (int i = 0; i < N; i++)
            X_data[i] = getter(refNode[ind[i]]);
        }
    };  // end class Tetra
    }   // end namespace Tetra

#endif /* tetra_h */
