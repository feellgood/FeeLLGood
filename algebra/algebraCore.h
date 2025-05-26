#ifndef ALGEBRA_CORE_H
#define ALGEBRA_CORE_H

#include <numeric> // inner_product

namespace algebra
{
template<typename T> using Vector = std::vector<T,std::allocator<T>>;

/** returns scalar product X.Y */
template <typename T>
T dot(const Vector<T> & X,const Vector<T> & Y)
	{
	return std::inner_product(X.begin(),X.end(),Y.begin(),(T)0);
	}

/** euclidian norm of vector X */
template <typename T>
T norm(Vector<T> & X) { return sqrt(fabs( dot<T>(X,X) )); }
}

#endif
