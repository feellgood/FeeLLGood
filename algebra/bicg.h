#include "algebra.h"

#ifndef BICG_H
#define BICG_H

namespace algebra
{
template<bool MASK, typename T>
T generic_bicg(iteration &iter, r_sparseMat& A, std::vector<T> & x,
               const std::vector<T> & rhs, const std::vector<T>& xd, const std::vector<int>& ld)
    {
    T rho_1(0.0), rho_2(0.0), alpha(0.0), beta(0.0), omega(0.0);
    const size_t DIM = x.size();
    if (rhs.size()!=DIM){std::cout << "rhs size mismatch" << std::endl; exit(1);}

    std::vector<T> p(DIM), phat(DIM), shat(DIM), r(DIM), rt(DIM), s(DIM), t(DIM), v(DIM), diag_precond(DIM), b(DIM);
    b.assign(rhs.begin(),rhs.end());// b = rhs;

    for(unsigned int i=0;i<diag_precond.size();i++)
        { diag_precond[i] = 1.0/A(i,i); }
    if(MASK)
        {
        mult(A, xd, v);
        sub(v, b);      // b -= A xd
        applyMask(ld,b);
        applyMask(ld,diag_precond);
        }
    iter.set_rhsnorm(norm(b));

    r.assign(b.begin(),b.end());// r = b;
    mult(A, x, v);         // v = A x;
    sub(v, r);             // r -= v; donc r = b - A x;

    if(MASK)
        { applyMask(ld,r); }

    rt.assign(r.begin(),r.end()); // copy(r, rt);
    p.assign(r.begin(),r.end()); // copy(r, p );

    while (!iter.finished(norm(r)))
        {
        rho_1 = dot(rt,r);
        if (!iter.first())
            {
            beta = (rho_1 / rho_2) * (alpha / omega);
            scaled(omega, v);           // v *= omega
            sub(v, p);                  // p -= v; so  p = p - omega v

            scaled(beta, p);            // p *= beta
            add(r, p);                  // p += r; so  p = r + beta p
            }

        p_direct(diag_precond, p, phat);// phat = M p;
        mult(A, phat, v);               //  v = A phat;
        if (MASK)
            { applyMask(ld,v); }

        alpha=rho_1/dot(v, rt);         // alpha = rho_1 /(v'*rtilde);
        s.assign(r.begin(), r.end());   // s = r
        scaled_add(v, -alpha, s);       // s += -alpha v; so s = r -alpha v

        if (iter.finished(norm(s)))
            {
            scaled_add(phat, alpha, x); // x += alpha * phat
            break;
            }

        p_direct(diag_precond, s, shat);// shat = M s;
        mult(A, shat, t);               //  t = A shat;
        if(MASK)
            { applyMask(ld,t); }

        omega = dot(t, s)/dot(t,t);    // omega = (t'* s) / (t'*t);
        scaled_add(phat, alpha, x);    // x += alpha * phat;
        scaled_add(shat, omega, x);    // x += omega * shat;

        scaled(omega, t);              // t *= omega
        r.assign(s.begin(), s.end());  // r = s
        sub(t, r);                     // r -= t; so r = s - omega t

        rho_2 = rho_1;
        ++iter;
        }
    if(MASK)
        { add(xd, x); }                // x += xd
    return iter.get_res()/iter.get_rhsnorm();
    }
} // end namespace alg
#endif
