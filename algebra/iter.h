#ifndef ITER_H
#define ITER_H

#include <iomanip>
#include <vector>

/** \file iter.h
\brief iteration class from GMM, with some adaptations and simplifications.

The Iteration object calculates if the solution has reached the desired accuracy,
 or if the maximum number of iterations has been reached.

The method finished() checks the convergence.
The first() method is used to determine the first iteration of the loop.
*/

namespace algebra
{
double norm(const std::vector<double> & X);

/**
 \class iteration
 monitor over the successive iterations if a value is converging or not
 */
class iteration
    {
    protected :

    /** Right hand side norm. */
    double rhsn;

/** Max. number of iterations. */
    size_t maxiter;

/** if noise > 0 iterations are printed. */
    int noise;

/** maximum residu. */
    double resmax;

/** iteration number. */
    size_t nit;

/** last computed residu. */
    double res;

/** true : info was written */
    bool written;

    public :
    /** constructor */
inline iteration(double r = 1.0E-8, int noi = 0, size_t mit = (size_t)(-1))
      : rhsn(1.0), maxiter(mit), noise(noi), resmax(r), nit(0), res(std::numeric_limits<double>::max()), written(false) {}

/** increment of the number of iterations */
inline void operator ++(int) { nit++; written = false; }

/** operator increment */
inline void operator ++() { (*this)++; }

/** true if iterations are starting */
inline bool first(void) { return nit == 0; }

/** set the "noisyness" (verbosity) of the solvers */
inline void set_noisy(int n) { noise = n; }

/** getter for resmax */
inline double get_resmax(void) const { return resmax; }

/** setter for resmax */
inline void set_resmax(double r) { resmax = r; }

/** getter for residu res */
inline double get_res() const { return res; }

/** getter for number of iterations */
inline size_t get_iteration(void) const { return nit; }

/** setter for the number of iterations */
inline void set_iteration(size_t i) { nit = i; }

/** getter for the maximum number of iterations */
inline size_t get_maxiter(void) const { return maxiter; }

/** setter for the maximum number of iterations */
inline  void set_maxiter(size_t i) { maxiter = i; }

/** getter for the right hand side norm value */
inline double get_rhsnorm(void) const { return rhsn; }

/** setter for the right hand side norm */
void set_rhsnorm(double r) { rhsn = r; }

/** return the monitored algo has converged or not according criteria fixed by right hand side norm rhsn */
inline bool converged(void)
    {
    if (std::isnan(res))
        { std::cout << "residu is NaN, algo cannot converge.\n"; exit(1); }
    return res <= rhsn * resmax;
    }

/** monitor the convergence through a number */
inline bool converged(double nr)
    {
    res = std::fabs(nr);
    return converged();
    }

/** monitor the convergence through a vector */
inline bool converged(const std::vector<double> &v)
    { return converged( norm(v) ); }

/** returns true if the algo has converged according the convergence criterias through a norm value nr */
inline bool finished(double nr)
    {
    if (noise > 0 && !written)
        {
        double a = (rhsn == 0) ? 1.0 : rhsn;
        converged(nr);
        std::cout << " iter " << std::setw(3) << nit << " residual " << std::setw(12) << std::fabs(nr) / a << std::endl;
        written = true;
        }
    return converged(nr);
    }
  };
}

#endif
