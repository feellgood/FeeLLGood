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

static const Eigen::Matrix<double,N,NPI> eigen_a = [] {
    Eigen::Matrix<double,N,NPI> tmp; tmp << a[0][0], a[0][1], a[0][2], a[0][3],
                                            a[1][0], a[1][1], a[1][2], a[1][3],
                                            a[2][0], a[2][1], a[2][2], a[2][3];
                                            return tmp; }();

/** \class prm
region number and material constants
*/
struct prm
    {
    std::string regName;   /**< region name */
    bool suppress_charges; /**< suppress charges if true */
    double Ks;             /**< uniaxial surface anisotropy constant */
    Eigen::Vector3d uk;    /**< anisotropy axis */

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
Face is a class containing the index references to nodes, its surface and its normal unit vector, it has a triangular shape and should not
be degenerated, orientation must be defined in adequation with the mesh
*/
class Fac : public element<N,NPI>
    {
public:
    /** constructor used by readMesh */
    inline Fac(const std::vector<Nodes::Node> &_p_node /**< [in] vector of nodes */,
               const int _NOD /**< [in] nb nodes */,
               const int _idx /**< [in] region index in region vector */,
               std::initializer_list<int> _i /**< [in] node index */)
        : element<N,NPI>(_p_node,_idx,_i)
        {
        Ms = 0;
        if (_NOD > 0)
            {
            zeroBasing();
            surf = calc_surf();
            n = calc_norm();
            }
        else
            {
            surf = 0.0;
            n = Eigen::Vector3d(0,0,0);
            }  // no index shift here if NOD == 0 : usefull while reordering face indices

        for(int i=0;i<NPI;i++)
            { weight[i] = 2.0 * surf * Facette::pds[i]; }
        }

    /** surface of the element */
    double surf;

    /** magnetization at saturation of the facette */
    double Ms;

    /** normal vector (unit vector) */
    Eigen::Vector3d n;

    /** interpolation function on the output of the getter.
    result = vec_nod * a 
    */
    void interpolation(std::function<Eigen::Vector3d(Nodes::Node)> getter /**< [in] */,
                       Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> result /**< [out] */) const
        {
        Eigen::Matrix<double,Nodes::DIM,N> vec_nod;
        for (int i = 0; i < N; i++) vec_nod.col(i) = getter(getNode(i));
        
        result = vec_nod * eigen_a;
        }

    /** interpolation function on the output of the getter,
     mind the transposition: result = transpose(scalar_nod) * a
    */
    void interpolation(std::function<double(Nodes::Node)> getter /**< [in] */,
                       Eigen::Ref<Eigen::Matrix<double,NPI,1>> result /**< [out] */) const
        {
        Eigen::Matrix<double,N,1> scalar_nod;
        for (int i = 0; i < N; i++) scalar_nod(i) = getter(getNode(i));
        
        result = scalar_nod.transpose() * eigen_a;
        }

    /** computes the integral contribution of the triangular face, the only contribution comes from Neel surface anisotropy */
    void integrales(Facette::prm const &params /**< [in] */);

    /** anisotropy energy of the facette */
    double anisotropyEnergy(Facette::prm const &param /**< [in] */,
                            Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> const u /**< [in] */) const;

    /** return surface charges and computes some corrections */
    Eigen::Matrix<double,NPI,1> charges(Facette::prm const &param /**< [in] */,
                                        std::function<Eigen::Vector3d(Nodes::Node)> getter /**< [in] */,
                                        std::vector<double> &corr /**< [in|out]*/ ) const;

    /** demagnetizing energy of the facette */
    double demagEnergy(Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> u /**< [in] */,
                       Eigen::Ref<Eigen::Matrix<double,NPI,1>> phi /**< [in] */) const;

    /** computes correction on potential*/
    double potential(std::function<Eigen::Vector3d(Nodes::Node)> getter, int i) const;

    /** lexicographic order on indices */
    inline bool operator<(const Fac &f) const
        {
        return (this->ind[0] < f.ind[0])
               || ((this->ind[0] == f.ind[0])
                   && ((this->ind[1] < f.ind[1])
                       || ((this->ind[1] == f.ind[1]) && (this->ind[2] < f.ind[2]))));
        }

    /** computes the norm to the face, returns a unit vector */
    inline Eigen::Vector3d calc_norm(void) const
        {
        Eigen::Vector3d _n = normal_vect();
        _n.normalize();
        return _n;
        }

    /** returns Gauss points in result = vec_nod*Facette::a */
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
    void orientate(void) {}// orientation is done in mesh::indexReorder

    /** computes surface of the face */
    inline double calc_surf(void) const { return 0.5 * normal_vect().norm(); }

    /** return normal to the triangular face, not normalized */
    inline Eigen::Vector3d normal_vect() const
        {
        Eigen::Vector3d p0p1 = getNode(1).p - getNode(0).p;
        Eigen::Vector3d p0p2 = getNode(2).p - getNode(0).p;

        return p0p1.cross(p0p2);
        }
    };  // end class Fac

    }  // namespace Facette

#endif /* facette_h */
