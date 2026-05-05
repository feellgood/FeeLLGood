#ifndef ALGEBRA_CORE_H
#define ALGEBRA_CORE_H

#include <numeric> // inner_product

namespace algebra
{
/** returns scalar product X.Y */
template <typename T>
T dot(const std::vector<T> & X,const std::vector<T> & Y)
	{
	return std::inner_product(X.begin(),X.end(),Y.begin(),(T)0);
	}

/** euclidian norm of vector X */
template <typename T>
T norm(std::vector<T> & X) { return sqrt(fabs( dot<T>(X,X) )); }
}

#endif
