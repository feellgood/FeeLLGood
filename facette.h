#ifndef facette_h
#define facette_h

/** \file facette.h
  \brief contains namespace Facette
  header containing Fac class, and some constants and a less_than operator to redo orientation of
  triangular faces
 */

#include "element.h"

/** \namespace Facette
 to grab altogether some constants and calculation functions for class Fac
 */
namespace Facette
    {
const int N = 3;   /**< number of sommits */
const int NPI = 4; /**< number of weights  */

constexpr double u[NPI] = {1 / 3., 1 / 5., 3 / 5.,
                           1 / 5.}; /**< some constants to build hat functions */
constexpr double v[NPI] = {1 / 3., 1 / 5., 1 / 5.,
                           3 / 5.}; /**< some constants  to build hat functions */
constexpr double pds[NPI] = {-27 / 96., 25 / 96., 25 / 96.,
                             25 / 96.}; /**< some constant weights  to build hat functions */

/** hat function constants */
constexpr double a[N][NPI] = {
        {1. - u[0] - v[0], 1. - u[1] - v[1], 1. - u[2] - v[2], 1. - u[3] - v[3]},
        {u[0], u[1], u[2], u[3]},
        {v[0], v[1], v[2], v[3]}};

/** \class prm
region number and material constants
*/
struct prm
    {
    std::string regName;   /**< region name */
    bool suppress_charges; /**< suppress charges if true */
    double Ks;             /**< uniaxial surface anisotropy constant */
    Pt::pt3D uk;           /**< anisotropy axis */

    /** print the struct parameters */
    inline void infos()
        {
        std::cout << "surface region name = " << regName
                  << " ; suppress charges = " << suppress_charges << std::endl;

        if (Ks != 0)
            {
            std::cout << "Ks*a = " << Ks << "*[ " << uk << "]" << std::endl;
            }
        else
            std::cout << "no surface anisotropy" << std::endl;
        };
    };

/** \class Fac
Face is a class containing the index references to nodes, it has a triangular shape and should not
be degenerated
*/
class Fac : public element<N,NPI>
    {
public:
    /** constructor for some unit tests */
    inline Fac(const std::vector<Nodes::Node> &_p_node /**< [in] vector of nodes */)
        : element<N,NPI>(_p_node), idxPrm(0)
        {
        indicesToZero();
        }

    /** constructor used by readMesh */
    inline Fac(const std::vector<Nodes::Node> &_p_node /**< [in] vector of nodes */,
               const int _NOD /**< [in] nb nodes */,
               const int _idx /**< [in] region index in region vector */,
               std::initializer_list<int> _i /**< [in] node index */)
        : element<N,NPI>(_p_node), idxPrm(_idx)
        {
        set_ind(_i);

        if (_NOD > 0)
            {
            zeroBasing();
            surf = calc_surf();
            }
        else
            {
            surf = 0.0;
            Ms = 0.0;
            }  // no index shift here if NOD == 0 : usefull while reordering face indices
        }

    int idxPrm; /**< index of the material parameters of the facette */

    double surf; /**< surface of the element */
    double Ms;   /**< magnetization at saturation of the face */

    /** weighted scalar product : factorized formulation: weight(1)=weight(2)=weight(3) */
    inline double weightedScalarProd(const double (&X)[NPI] /**< [in] */) const
        {
        return (X[0] * weight(0) + (X[1] + X[2] + X[3]) * weight(1));
        }

    /** interpolation template; T == 3D vector field or T == double .The getter function is given as
    a parameter in order to know what part of the node you want to interpolate if T is Pt::pt3D then
    result = vec_nod * a : mind the transposition due to the fact that vec_nod is an array of pt3D
    if T is double then result = transpose(scalar_nod) * a
        */

    // tiny::mult<double, DIM, N, NPI> (vec_nod, a, result); //if T == PT::pt3D
    // tiny::transposed_mult<double, N, NPI> (scalar_nod, a, result); //if T == double
    template<class T>
    void interpolation(std::function<T(Nodes::Node)> getter /**< [in] */,
                       T (&result)[NPI] /**< [out] */) const
        {
        const T x0(getter(refNode[ind[0]])), x1(getter(refNode[ind[1]])),
                x2(getter(refNode[ind[2]]));

        result[0] = x0;
        result[0] += x1;
        result[0] += x2;
        result[1] = (result[0] + 2.0 * x0) / 5.0;
        result[2] = (result[0] + 2.0 * x1) / 5.0;
        result[3] = (result[0] + 2.0 * x2) / 5.0;
        result[0] /= 3.0;
        }

    /** interpolation for 3D vector field : the getter function is given as a parameter in order to
     know what part of the node you want to interpolate This function is only usefull for
     fmm_demag.h, in other parts of the code the generic interpolation template is directly called
     */
    inline void interpolation(std::function<Pt::pt3D(Nodes::Node)> getter /**< [in] */,
                              Pt::pt3D (&result)[NPI] /**< [out] */) const
        {
        interpolation<Pt::pt3D>(getter, result);
        }

    /** basic infos */
    inline void infos() const
        {
        std::cout << "idxPrm:" << idxPrm << "ind: ";
        print_indices();
        };

    /** computes the integral contribution of the triangular face */
    void integrales(std::vector<Facette::prm> const &params /**< [in] */,
                    Pt::pt3D (&BE)[N] /**< [out] */) const;

    /** anisotropy energy of the facette */
    double anisotropyEnergy(Facette::prm const &param /**< [in] */,
                            const Pt::pt3D (&u)[NPI] /**< [in] */) const;

    /** surface charges  */
    void charges(std::function<Pt::pt3D(Nodes::Node)> getter, std::vector<double> &srcDen,
                 std::vector<double> &corr, int &nsrc) const;

    /** demagnetizing energy of the facette */
    double demagEnergy(const Pt::pt3D (&u)[NPI] /**< [in] */,
                       const double (&phi)[NPI] /**< [in] */) const;

    /** projections using Nodes::projection templates */
    inline void projection(double (&K)[3 * N][3 * N], Pt::pt3D L[N])
        {
        Nodes::projection_mat<Facette::Fac, N>(*this, K);
        Nodes::projection_vect<Facette::Fac, N>(*this, L);
        }

    /** computes weight coefficients */
    inline double weight(const int i) const { return 2.0 * surf * Facette::pds[i]; }

    /** computes correction on potential*/
    double potential(std::function<Pt::pt3D(Nodes::Node)> getter, int i) const;

    /** computes correction values */
    void calcCorr(std::function<const Pt::pt3D(Nodes::Node)> getter, std::vector<double> &corr,
                  Pt::pt3D (&u)[NPI]) const;

    /** lexicographic order on indices */
    inline bool operator<(const Fac &f) const
        {
        return (this->ind[0] < f.ind[0])
               || ((this->ind[0] == f.ind[0])
                   && ((this->ind[1] < f.ind[1])
                       || ((this->ind[1] == f.ind[1]) && (this->ind[2] < f.ind[2]))));
        }

    /** computes the norm to the face, returns a unit vector */
    inline Pt::pt3D calc_norm(void) const
        {
        Pt::pt3D n = normal_vect();
        n.normalize();
        return n;
        }

private:

    void orientate(void) {}// orientation is done in mesh::indexReorder

    /** computes surface of the face */
    inline double calc_surf(void) const { return 0.5 * normal_vect().norm(); }

    /** return normal to the triangular face, not normalized */
    inline Pt::pt3D normal_vect() const
        {
        Pt::pt3D p0 = refNode[ind[0]].p;
        Pt::pt3D p1 = refNode[ind[1]].p;
        Pt::pt3D p2 = refNode[ind[2]].p;

        return ((p1 - p0) * (p2 - p0));
        }
    };  // end class Fac

    }  // namespace Facette

#endif /* facette_h */
