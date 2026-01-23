#include "algebra.h"

#ifndef BICG_H
#define BICG_H

namespace algebra
{
/** solve A x = rhs. Algo is stabilized biconjugate gradient with diagonal preconditioner,
vectors x and rhs must have the same size.
The status of the convergence is returned in iter.status as well as total number of iterations and error
 */
template <typename T>
void bicg(iteration<T> &iter, SparseMatrix& A, std::vector<T> & x, const std::vector<T> & rhs)
    {
    const size_t DIM = x.size();
    T rho_1(0.0), rho_2(0.0), alpha(0.0), beta(0.0), omega(0.0);
    std::vector<T> p(DIM), phat(DIM,0), shat(DIM), r(DIM), rt(DIM), s(DIM), t(DIM), v(DIM), diag_precond(DIM), b(DIM);
    b.assign(rhs.begin(),rhs.end());// b = rhs;

    A.build_diag_precond<T>(diag_precond);
    iter.set_rhsnorm(norm(b));
    r.assign(b.begin(),b.end());// r = b;
    mult(A, x, v);         // v = A x;
    sub(v, r);             // r -= v; so r = b - A x;

    rt.assign(r.begin(),r.end()); // rt = r;
    p.assign(r.begin(),r.end()); // p = r;
    while (!iter.finished(norm(r)))
        {
        rho_1 = dot(rt,r);

        if(iter.get_iteration() > 0)
            {
            if((rho_2 == 0)||(omega == 0))
                {
                iter.status = CANNOT_CONVERGE;
                break;
                }
            beta = rho_1 / rho_2;
            beta *=  alpha / omega;
            scaled(omega, v);           // v *= omega
            sub(v, p);                  // p -= v; so  p = p - omega v

            scaled(beta, p);            // p *= beta
            add(r, p);                  // p += r; so  p = r + beta p
            }
        p_direct(diag_precond, p, phat);// phat = M p;
        mult(A, phat, v);               //  v = A phat;
        alpha=rho_1/dot(v, rt);         // alpha = rho_1 /(v'*rtilde);
        s.assign(r.begin(), r.end());   // s = r
        scaled_add(v, -alpha, s);       // s += -alpha v; so s = r -alpha v

        if (iter.finished(norm(s)))
            {
            scaled_add(phat, alpha, x); // x += alpha * phat
            break;
            }
        else if ((iter.status == ITER_OVERFLOW)||(iter.status == CANNOT_CONVERGE))
            {break;}

        p_direct(diag_precond, s, shat);// shat = M s;
        mult(A, shat, t);               //  t = A shat;
        omega = dot(t, s)/dot(t,t);    // omega = (t'* s) / (t'*t);
        scaled_add(phat, alpha, x);    // x += alpha * phat;
        scaled_add(shat, omega, x);    // x += omega * shat;
        scaled(omega, t);              // t *= omega
        r.assign(s.begin(), s.end());  // r = s
        sub(t, r);                     // r -= t; so r = s - omega t
        rho_2 = rho_1;
        ++iter;
        }
    }

/** directional stabilized biconjugate gradient (mask) with diagonal preconditioner with Dirichlet conditions, returns residu
 * iter is an iteration object
 * the linear system to solve is A x = rhs
 * ld is a mask: a list of indices that will be zeroed when using applyMask(ld, vector)
 * xd is a vector containing some values on the nodes where Dirichlet applies
 */
template<typename T>
void bicg_dir(iteration<T> &iter, SparseMatrix& A, std::vector<T> & x, const std::vector<T> & rhs,
           const std::vector<T>& xd, const std::vector<int>& ld)
    {
    T rho_1(0.0), rho_2(0.0), alpha(0.0), beta(0.0), omega(0.0);
    const size_t DIM = x.size();
    std::vector<T> p(DIM), phat(DIM), shat(DIM), r(DIM), rt(DIM), s(DIM), t(DIM), v(DIM), diag_precond(DIM), b(rhs);

    A.build_diag_precond<T>(diag_precond);
    mult(A, xd, v);
    sub(v, b);      // b -= A xd; so b = rhs - A xd
    applyMask(ld,b);
    applyMask(ld,diag_precond);

    iter.set_rhsnorm(norm(b));

    r.assign(b.begin(),b.end());// r = b;
    mult(A, x, v);         // v = A x;
    sub(v, r);             // r -= v; so r = b - A x;
    applyMask(ld,r);
    rt.assign(r.begin(),r.end()); // rt = r;
    p.assign(r.begin(),r.end()); // p = r;
    while (!iter.finished(norm(r)))
        {
        rho_1 = dot(rt,r);

        if(iter.get_iteration() > 0)
            {
            if((rho_2 == 0)||(omega == 0))
                {
                iter.status = CANNOT_CONVERGE;
                break;
                }
            beta = (rho_1 / rho_2) * (alpha / omega);
            scaled(omega, v);           // v *= omega
            sub(v, p);                  // p -= v; so  p = p - omega v

            scaled(beta, p);            // p *= beta
            add(r, p);                  // p += r; so  p = r + beta p
            }

        p_direct(diag_precond, p, phat);// phat = direct_product(diag_precond, p);
        mult(A, phat, v);               //  v = A phat;
        applyMask(ld,v);
        alpha=rho_1/dot(v, rt);         // alpha = rho_1 /(v'*rtilde);
        s.assign(r.begin(), r.end());   // s = r
        scaled_add(v, -alpha, s);       // s += -alpha v; so s = r -alpha v

        if (iter.finished(norm(s)))
            {
            scaled_add(phat, alpha, x); // x += alpha * phat
            break;
            }
        else if ((iter.status == ITER_OVERFLOW)||(iter.status == CANNOT_CONVERGE))
            {break;}

        p_direct(diag_precond, s, shat);// shat = direct_product(diag_precond, s);
        mult(A, shat, t);               //  t = A shat;
        applyMask(ld,t);
        omega = dot(t, s)/dot(t,t);    // omega = (t'* s) / (t'*t);
        scaled_add(phat, alpha, x);    // x += alpha * phat;
        scaled_add(shat, omega, x);    // x += omega * shat;

        scaled(omega, t);              // t *= omega
        r.assign(s.begin(), s.end());  // r = s
        sub(t, r);                     // r -= t; so r = s - omega t

        rho_2 = rho_1;
        ++iter;
        }
    add(xd, x);                // x += xd
    }

/** directional stabilized biconjugate gradient (mask is ld) with diagonal preconditioner, returns residu
 * iter is an iteration object
 * the linear system to solve is A x = rhs
 * ld is a mask: a list of indices that will be zeroed when using applyMask(ld, vector)
 */
template<typename T>
T bicg_dir(iteration<T> &iter, SparseMatrix& A, std::vector<T> & x, const std::vector<T> & rhs, const std::vector<int>& ld)
    {
    T rho_1(0.0), rho_2(0.0), alpha(0.0), beta(0.0), omega(0.0);
    const size_t DIM = x.size();
    if (rhs.size()!=DIM){std::cout << "rhs size mismatch" << std::endl; exit(1);}

    std::vector<T> p(DIM), phat(DIM), shat(DIM), r(DIM), rt(DIM), s(DIM), t(DIM), v(DIM), diag_precond(DIM), b(rhs);

    A.build_diag_precond<T>(diag_precond);
    applyMask(ld,b);
    applyMask(ld,diag_precond);

    iter.set_rhsnorm(norm(b));

    r.assign(b.begin(),b.end());// r = b;
    mult(A, x, v);         // v = A x;
    sub(v, r);             // r -= v; donc r = b - A x;
    applyMask(ld,r);
    rt.assign(r.begin(),r.end()); // rt = r;
    p.assign(r.begin(),r.end()); // p = r;
    while (!iter.finished(norm(r)))
        {
        rho_1 = dot(rt,r);

        if(iter.get_iteration() > 0)
            {
            if((rho_2 == 0)||(omega == 0))
                {
                iter.status = CANNOT_CONVERGE;
                break;
                }
            beta = (rho_1 / rho_2) * (alpha / omega);
            scaled(omega, v);           // v *= omega
            sub(v, p);                  // p -= v; so  p = p - omega v

            scaled(beta, p);            // p *= beta
            add(r, p);                  // p += r; so  p = r + beta p
            }

        p_direct(diag_precond, p, phat);// phat = direct_product(diag_precond, p);
        mult(A, phat, v);               //  v = A phat;
        applyMask(ld,v);
        alpha=rho_1/dot(v, rt);         // alpha = rho_1 /(v'*rtilde);
        s.assign(r.begin(), r.end());   // s = r
        scaled_add(v, -alpha, s);       // s += -alpha v; so s = r -alpha v

        if (iter.finished(norm(s)))
            {
            scaled_add(phat, alpha, x); // x += alpha * phat
            break;
            }
        else if ((iter.status == ITER_OVERFLOW)||(iter.status == CANNOT_CONVERGE))
            {break;}

        p_direct(diag_precond, s, shat);// shat = direct_product(diag_precond, s);
        mult(A, shat, t);               //  t = A shat;
        applyMask(ld,t);
        omega = dot(t, s)/dot(t,t);    // omega = (t'* s) / (t'*t);
        scaled_add(phat, alpha, x);    // x += alpha * phat;
        scaled_add(shat, omega, x);    // x += omega * shat;

        scaled(omega, t);              // t *= omega
        r.assign(s.begin(), s.end());  // r = s
        sub(t, r);                     // r -= t; so r = s - omega t

        rho_2 = rho_1;
        ++iter;
        }
    return iter.get_res()/iter.get_rhsnorm();
    }
} // end namespace alg
#endif
