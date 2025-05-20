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
	T init_val=0;
	/*
	if(std::isnan(val))
	    { std::cout <<"dot:Warning: NaN value.\n"; exit(1); }
	else if(!std::isfinite(val))
	    { std::cout <<"dot:Warning: inf value.\n"; exit(1); }
	*/
	return std::inner_product(X.begin(),X.end(),Y.begin(),init_val);
	}

/** euclidian norm of vector X */
template <typename T>
T norm(Vector<T> & X) { return sqrt(fabs( dot<T>(X,X) )); }
}

#endif
