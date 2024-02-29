#include <eigen3/Eigen/Dense>
#include "tetra.h"

template <int nbNod>
void dummyNodes(std::vector<Nodes::Node> &node)
    {
    node.resize(nbNod);
    Eigen::Vector3d zero(0,0,0),u0(0, 0, 0), v0(0, 0, 0), u(0, 0, 0), v(0, 0, 0);
    double phi0(0), phi(0), phiv0(0), phiv(0);

    if(nbNod ==4)
        {
        Eigen::Vector3d p0(0, 0, 0), p1(1, 0, 0), p2(0, 1, 0), p3(0, 0, 1);
        Nodes::Node n1 = {p0, zero, zero, {{u0, v0, phi0, phiv0}, {u,v,phi,phiv}} };
        Nodes::Node n2 = {p1, zero, zero, {{u0, v0, phi0, phiv0}, {u,v,phi,phiv}} };
        Nodes::Node n3 = {p2, zero, zero, {{u0, v0, phi0, phiv0}, {u,v,phi,phiv}} };
        Nodes::Node n4 = {p3, zero, zero, {{u0, v0, phi0, phiv0}, {u,v,phi,phiv}} };
        node = {n1,n2,n3,n4};
        }
    else if (nbNod == 3)
        {
        Eigen::Vector3d p1(1, 0, 0), p2(0, 1, 0), p3(1, 1, 0);
        Nodes::Node n1 = {p1, zero, zero, {{u0, v0, phi0, phiv0}, {u,v,phi,phiv}} };
        Nodes::Node n2 = {p2, zero, zero, {{u0, v0, phi0, phiv0}, {u,v,phi,phiv}} };
        Nodes::Node n3 = {p3, zero, zero, {{u0, v0, phi0, phiv0}, {u,v,phi,phiv}} };
        node = {n1,n2,n3};
        }
    }

/** returns a unit vector pointing to (theta,phi) in spherical coordinates */
Eigen::Vector3d rand_vec3d(double theta, double phi)
    {
    const double si_t = sin(theta);
    return Eigen::Vector3d( si_t*cos(phi), si_t*sin(phi), cos(theta) );
    }

template<int nbVect>
double sq_frobenius_norm(Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,nbVect>> const X)
    {
    double val(0.0);
    for (int i = 0; i < nbVect; i++)
        {
        val += X.col(i).squaredNorm();
        }
    return val;
    }

double sq_dist(double _x[Nodes::DIM][Tetra::NPI], Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> X)
    {
    double val(0.0);

    for (int i = 0; i < Tetra::NPI; i++)
        for (int j = 0; j < Nodes::DIM; j++)
            { val += Nodes::sq(_x[j][i] - X(j,i)); }
    return val;
    }

double sq_dist(double _x[Tetra::NPI], double _y[Tetra::NPI], double _z[Tetra::NPI],
               Eigen::Ref<Eigen::Matrix<double,Nodes::DIM,Tetra::NPI>> X)
    {
    double val(0.0);

    for (int i = 0; i < Tetra::NPI; i++)
        {
        val += Nodes::sq(_x[i] - X.col(i).x());
        val += Nodes::sq(_y[i] - X.col(i).y());
        val += Nodes::sq(_z[i] - X.col(i).z());
        }

    return val;
    }

double det(const double M[Nodes::DIM][Nodes::DIM])
    {
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
            - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]));
    }

void inverse(double M[Nodes::DIM][Nodes::DIM], double detM)
    {
    double m00 = M[0][0];
    double m01 = M[0][1];
    double m02 = M[0][2];
    double m10 = M[1][0];
    double m11 = M[1][1];
    double m12 = M[1][2];
    double m20 = M[2][0];
    double m21 = M[2][1];
    double m22 = M[2][2];

    M[0][0] = (m11 * m22 - m12 * m21) / detM;
    M[0][1] = (m02 * m21 - m01 * m22) / detM;
    M[0][2] = (m01 * m12 - m02 * m11) / detM;
    M[1][0] = (m12 * m20 - m10 * m22) / detM;
    M[1][1] = (m00 * m22 - m02 * m20) / detM;
    M[1][2] = (m02 * m10 - m00 * m12) / detM;
    M[2][0] = (m10 * m21 - m11 * m20) / detM;
    M[2][1] = (m01 * m20 - m00 * m21) / detM;
    M[2][2] = (m00 * m11 - m01 * m10) / detM;
    }

