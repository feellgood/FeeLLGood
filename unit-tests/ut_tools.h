#include <eigen3/Eigen/Dense>
#include "pt3D.h"

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

