#include <eigen3/Eigen/Dense>
#include "pt3D.h"

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

/**
 frobenius norm of a table of pt3D
 */
template<int nbVect>
double sq_frobenius_norm(const Pt::pt3D X[nbVect])
    {
    double val(0.0);
    for (int i = 0; i < nbVect; i++)
        {
        val += X[i].norm2();
        }
    return val;
    }

template<int nbVect>
double sq_frobenius_norm(Eigen::Ref<Eigen::Matrix<double,Pt::DIM,nbVect>> const X)
    {
    double val(0.0);
    for (int i = 0; i < nbVect; i++)
        {
        val += X.col(i).squaredNorm();
        }
    return val;
    }

