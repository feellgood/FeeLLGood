#ifndef ITER_H
#define ITER_H

#include <iomanip>

/** \file iter.h
\brief iteration class from GMM, with some adaptations and simplifications.

The Iteration object calculates if the solution has reached the desired accuracy,
 or if the maximum number of iterations has been reached.

The method finished() checks the convergence.
*/

namespace algebra
{
/** status of iterative algorithm
UNDEFINED means no iteration done
CONVERGED means algo succeeded after some iterations < MAXITER to achieve residu<TOL
ITER_OVERFLOW means number of iterations has exceeded MAXITER
CANNOT_CONVERGE means an algebric operation leads to nan (might happen while computing beta in bicg inner loop)
*/
enum algoStatus
    {
    UNDEFINED = -1,
    CONVERGED = 0,
    ITER_OVERFLOW = 1,
    CANNOT_CONVERGE = 2
    };

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
    const int maxiter;

    /** if true iterations are printed. */
    const bool noise;

    /** iteration number. */
    int nit;

    /** last computed residu. */
    T res;

    /** true : info was written */
    bool written;

    private:
    /** monitored solver name */
    const std::string solver_name;

    /** string version of the status */
    std::string str_status(void) const
        {
        std::string s;
        switch(status)
            {
            case UNDEFINED:
                s="UNDEFINED";
            break;
            case CONVERGED:
                s="CONVERGED";
            break;
            case ITER_OVERFLOW:
                s="ITER_OVERFLOW";
            break;
            case CANNOT_CONVERGE:
                s="CANNOT_CONVERGE";
            break;
            }
        return s;
        }

    public :
    /** constructor */
    iteration(const std::string name, T r, bool _noise, int _maxiter): rhsn(1.0), maxiter(_maxiter),
        noise(_noise), nit(0), res(std::numeric_limits<T>::max()), written(false), solver_name(name), resmax(r)
        { status = UNDEFINED; }

    /** status of the monitored algorithm */
    algebra::algoStatus status;

    /** maximum tolerated residu. */
    const T resmax;

    /** reset all values to monitor a new algorithm exectution */
    void reset(void)
        {
        rhsn=1.0;
        nit = 0;
        res = std::numeric_limits<T>::max();
        written = false;
        status = UNDEFINED;
        }

    /** return a string aggregating the status of the monitored algorithm, the number of iterations and the residu */
    std::string infos(void) const
        {
        std::stringstream sstr;
        sstr << solver_name << " status " << str_status() << " after " << nit << " iterations, residu= " << res;
        return sstr.str();
        }

    /** increment of the number of iterations */
    void operator ++(int)
        {
        nit++;
        written = false;
        if (nit >= maxiter)
            { status = ITER_OVERFLOW; }
        }

    /** operator increment */
    void operator ++() { (*this)++; }

    /** getter for residu res */
    T get_res() const { return res; }

    /** getter for number of iterations */
    int get_iteration() const { return nit; }

    /** getter for the right hand side norm value */
    T get_rhsnorm() const { return rhsn; }

    /** setter for the right hand side norm */
    void set_rhsnorm(T r) { rhsn = r; }

    /** return the monitored algo has converged or not
        according to criteria fixed by right hand side norm rhsn */
    bool converged()
        {
        bool cv(false);
        if (std::isnan(res))
            { status = CANNOT_CONVERGE; }
        else
            {
            cv = (res <= rhsn * resmax);
            if(cv)
                { status = CONVERGED; }
            }
        return cv;
        }

    /** monitor the convergence through a number (the norm of a vector) */
    bool converged(T nr)
        {
        res = std::fabs(nr);
        return converged();
        }

    /** returns true if the algo has converged according the convergence criterias through a norm value nr */
    bool finished(T nr)
        {
        if (noise && !written)
            {
            T a = (rhsn == 0) ? 1.0 : rhsn;
            converged(nr);
            std::cout << " iter " << std::setw(3) << nit << " residual " << std::setw(12) << std::fabs(nr) / a << std::endl;
            written = true;
            }
        return converged(nr);
        }
    };// end class
}//end namespace
#endif
