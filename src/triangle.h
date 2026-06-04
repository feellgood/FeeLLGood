#ifndef triangle_h
#define triangle_h

/** \file triangle.h
  \brief contains namespace Triangle
  header containing Tri class, and some constants and a less_than operator to redo orientation of
  triangles
 */

#include "element.h"

/** \namespace Triangle
 to grab altogether some constants and calculation functions for class Tri
 */
namespace Triangle
    {
/** number of sommits */
const int N = 3;

#if ONE_GAUSS_POINT  // Single Gauss point at barycenter of triangle
    /** number of Gauss points */
    const int NPI = 1;

    /** u coordinate of single Gauss point at barycenter */
    constexpr double u[NPI] = {1. / 3.};

     /** v coordinate of single Gauss point at barycenter */
    constexpr double v[NPI] = {1. / 3.};

     /** weight for single point integration (area of reference triangle = 1/2) */
    constexpr double pds[NPI] = {1. / 2.};

    /** hat function constants evaluated at barycenter
    All shape functions are equal to 1/3 at barycenter
    */
    constexpr double a[N][NPI] = { {1. - u[0] - v[0]}, {u[0]}, {v[0]}};

    /** eigen matrix a initialization from C array - shape functions evaluated at Gauss point */
    const Eigen::Matrix<double,N,NPI> eigen_a =
        (Eigen::MatrixXd(N,NPI) << a[0][0], a[1][0], a[2][0] ).finished();
#else
    /** number of Gauss points  */
    const int NPI = 4;

    /** some constants to build hat functions */
    constexpr double u[NPI] = {1 / 3., 1 / 5., 3 / 5., 1 / 5.};

    /** some constants  to build hat functions */
    constexpr double v[NPI] = {1 / 3., 1 / 5., 1 / 5., 3 / 5.};

    /** some constant weights  to build hat functions */
    constexpr double pds[NPI] = {-27 / 96., 25 / 96., 25 / 96., 25 / 96.};

    /** hat function constants */
    constexpr double a[N][NPI] = {
        {1. - u[0] - v[0], 1. - u[1] - v[1], 1. - u[2] - v[2], 1. - u[3] - v[3]},
        {u[0], u[1], u[2], u[3]},
        {v[0], v[1], v[2], v[3]}};

    /** eigen matrix a initialization from C array */
    const Eigen::Matrix<double,N,NPI> eigen_a =
            (Eigen::MatrixXd(N,NPI) << a[0][0], a[0][1], a[0][2], a[0][3],
                                       a[1][0], a[1][1], a[1][2], a[1][3],
                                       a[2][0], a[2][1], a[2][2], a[2][3] ).finished();
#endif
/** \class prm
region number and material constants
*/
struct prm
    {
    std::string regName;   /**< region name */
    bool suppress_charges; /**< suppress charges if true */
    double Ks;             /**< uniaxial surface anisotropy constant */
    Eigen::Vector3d uk;    /**< anisotropy axis */
    double V;              /**< fixed potential (boundary condition for electrostatic
                             sub-problem) */
    double jn;             /**< \f$ j \dot \hat{n} \f$ normal current density (boundary condition
                             for electrostatic and spin diffusion sub-problems) */
    Eigen::Vector3d uP;    /**< Polarization spin accumulation orientation (unit vector
                             for spin diffusion problem) */
    Eigen::Vector3d s;     /**< spin diffusion vector (for spin diffusion boundary conditions) */

    /** print the region surface parameters, optional boolean spinAcc for spin accumulation datas */
    inline void infos(const bool spinAcc = false) const
        {
        std::cout << "surface region name = " << regName
                  << " ; suppress charges = " << suppress_charges << std::endl;

        if (Ks != 0)
            {
            std::cout << "Ks*a = " << Ks << "*[ " << uk << "]" << std::endl;
            }
        else
            std::cout << "no surface anisotropy" << std::endl;

        if (spinAcc)
            {
            std::cout<< "V= " << V << "; jn= " << jn << "; uP= " << uP << "; s= " << s << std::endl;
            }
        };
    };

/** \class Tri
Tri is a class containing the index references to nodes, its surface and its normal unit vector,
it has a triangular shape and should not be degenerated, orientation must be defined in adequation
with the mesh
*/
class Tri : public element<N,NPI>
    {
public:
    /** constructor used by readMesh
     * */
    inline Tri(const std::vector<Nodes::Node> &_p_node /**< [in] vector of nodes */,
               const int _idx /**< [in] region index in region vector */,
               std::initializer_list<int> _i /**< [in] node index */)
        : element<N,NPI>(_p_node,_idx,_i),surf(0),dMs(0)
        {
        if (existNodes())
            {
            surf = calc_surf();
            n = calc_norm();
            for(int i=0;i<NPI;i++)
                { weight[i] = 2.0 * surf * Triangle::pds[i]; }
            }
        else
            { std::cerr<<"Error: Tri constructor has out of bound index\n"; }
        }

    /** surface of the element */
    double surf;

    /** difference (d for delta) of magnetization at the triangle between corresponding
     * tetrahedrons */
    double dMs;

    /** normal vector (unit vector) */
    Eigen::Vector3d n;

    /** interpolation function on the output of the getter.
    result = vec_nod * a 
    */
    void interpolation(const std::function<Eigen::Vector3d(Nodes::Node)>& getter /**< [in] */,
                       Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> result /**< [out] */) const
        {
        Eigen::Matrix<double,Nodes::DIM,N> vec_nod;
        for (int i = 0; i < N; i++) vec_nod.col(i) = getter(getNode(i));
        
        result = vec_nod * eigen_a;
        }

    /** interpolation function on the output of the getter,
     mind the transposition: result = transpose(scalar_nod) * a
    */
    void interpolation(const std::function<double(Nodes::Node)>& getter /**< [in] */,
                       Eigen::Ref<Eigen::Matrix<double,NPI,1>> result /**< [out] */) const
        {
        Eigen::Matrix<double,N,1> scalar_nod;
        for (int i = 0; i < N; i++) scalar_nod(i) = getter(getNode(i));
        
        result = scalar_nod.transpose() * eigen_a;
        }

    /** computes the integral contribution of the triangular triangle, the only contribution comes
     * from Neel surface anisotropy */
    void integrales(Triangle::prm const &params /**< [in] */);

    /** anisotropy energy of the triangle */
    double anisotropyEnergy(Triangle::prm const &param /**< [in] */,
            const Eigen::Ref<const Eigen::Matrix<double,Nodes::DIM,NPI>> u /**< [in] */) const;

    /** computes surface charges, return result on NPI */
    Eigen::Matrix<double,NPI,1>  charges(const double &dMs /**< [in] */,
            const std::function<Eigen::Vector3d(const Nodes::Node&)> &getter /**< [in] */) const override;

    /** computes correction on surface charges from localCharges input, result directly stored in
     * corr */
    void correctionCharges(const std::function<Eigen::Vector3d(Nodes::Node)>& getter /**< [in] */,
                           Eigen::Matrix<double,Triangle::NPI,1> &localCharges /**< [in] */,
                           std::vector<double> &corr /**< [in|out]*/) const;

    /** demagnetizing energy of the triangle */
    double demagEnergy(Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,NPI>> u /**< [in] */,
                       const Eigen::Ref<const Eigen::Matrix<double,NPI,1>> phi /**< [in] */) const;

    /** computes correction on potential*/
    double potential(const std::function<Eigen::Vector3d(Nodes::Node)>& getter, int i) const;

    /** Two triangles are equal if they have the same set of vertices, i.e. if their lists of
     * node indices are the same to within a permutation.
     * Warning: two equal triangles may have different orientations. */
    bool operator==(const Tri &f) const
        {
        /* As we rule out degenerate triangles, it is sufficient to test whether all the vertices
         * of f are also vertices of this. */
        return std::find(ind.begin(), ind.end(), f.ind[0]) != ind.end()
            && std::find(ind.begin(), ind.end(), f.ind[1]) != ind.end()
            && std::find(ind.begin(), ind.end(), f.ind[2]) != ind.end();
        }

    /** computes the norm to the triangle, returns a unit vector */
    inline Eigen::Vector3d calc_norm(void) const
        {
        Eigen::Vector3d _n = normal_vect();
        _n.normalize();
        return _n;
        }

    /** returns Gauss points */
    Eigen::Matrix<double,Nodes::DIM,NPI> getPtGauss(void) const override
        {
        Eigen::Matrix<double,Nodes::DIM,N> vec_nod;
        for (int i = 0; i < N; i++)
            {
            vec_nod.col(i) << getNode(i).p;
            }
        return vec_nod*eigen_a;
        }

    /** computes surface of the triangle */
    inline double calc_surf(void) const { return 0.5 * normal_vect().norm(); }

private:
    /** do nothing function: the orientation of the triangles does not change */
    void orientate(void) override {}

    /** return normal to the triangle, not normalized */
    inline Eigen::Vector3d normal_vect() const
        {
        Eigen::Vector3d p0p1 = getNode(1).p - getNode(0).p;
        Eigen::Vector3d p0p2 = getNode(2).p - getNode(0).p;

        return p0p1.cross(p0p2);
        }
    };  // end class Tri

    }  // namespace Triangle

#endif /* triangle_h */
