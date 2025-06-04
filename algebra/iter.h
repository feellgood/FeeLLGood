#ifndef ITER_H
#define ITER_H

#include <iomanip>

/** \file iter.h
\brief iteration class from GMM, with some adaptations and simplifications.

The Iteration object calculates if the solution has reached the desired accuracy,
 or if the maximum number of iterations has been reached.

The method finished() checks the convergence.
The first() method is used to determine the first iteration of the loop.
*/

namespace algebra
{
/**
 \class iteration
 monitor over the successive iterations if a value is converging or not
 */
template <typename T>
class iteration
    {
    protected:

    /** Right hand side norm. */
    T rhsn;

/** Max. number of iterations. */
    size_t maxiter;

/** if noise > 0 iterations are printed. */
    int noise;

/** maximum residu. */
    T resmax;

/** iteration number. */
    size_t nit;

/** last computed residu. */
    T res;

/** true : info was written */
    bool written;

    public :
    /** constructor */
    iteration(T r = 1.0E-8, int noi = 0, size_t mit = (size_t)(-1)): rhsn(1.0), maxiter(mit),
        noise(noi), resmax(r), nit(0), res(std::numeric_limits<T>::max()), written(false) {}

/** increment of the number of iterations */
void operator ++(int) { nit++; written = false; }

/** operator increment */
void operator ++() { (*this)++; }

/** true if iterations are starting */
bool first() { return nit == 0; }

/** set the "noisyness" (verbosity) of the solvers */
void set_noisy(int n) { noise = n; }

/** getter for resmax */
T get_resmax() const { return resmax; }

/** setter for resmax */
void set_resmax(T r) { resmax = r; }

/** getter for residu res */
T get_res() const { return res; }

/** getter for number of iterations */
size_t get_iteration() const { return nit; }

/** getter for the maximum number of iterations */
size_t get_maxiter() const { return maxiter; }

/** setter for the maximum number of iterations */
void set_maxiter(size_t i) { maxiter = i; }

/** getter for the right hand side norm value */
T get_rhsnorm() const { return rhsn; }

/** setter for the right hand side norm */
void set_rhsnorm(T r) { rhsn = r; }

/** monitor the convergence through a number */
bool converged(T nr)
    {
    res = std::fabs(nr);
    if (std::isnan(res))
        { std::cout << "residu is NaN, algo cannot converge.\n"; exit(1); }
    return (res <= rhsn * resmax);
    }

/** returns true if the algo has converged according the convergence criterias through a norm value nr */
bool finished(T nr)
    {
    if (noise > 0 && !written)
        {
        T a = (rhsn == 0) ? 1.0 : rhsn;
        converged(nr);
        std::cout << " iter " << std::setw(3) << nit << " residual " << std::setw(12) << std::fabs(nr) / a << std::endl;
        written = true;
        }
    return converged(nr);
    }
  };
}

#endif
