#include "algebra.h"

#ifndef CG_H
#define CG_H

namespace algebra
{
/** conjugate gradient with diagonal preconditioner, returns residu. */
template<typename T>
T cg(iteration<T> &iter, r_sparseMat& A, std::vector<T> & x, const std::vector<T> & rhs)
    {
    T rho, rho_1(0.0);
    const size_t DIM = x.size();
    if (rhs.size()!=DIM){std::cout << "rhs size mismatch" << std::endl; exit(1);}

    std::vector<T> p(DIM),q(DIM),r(DIM),z(DIM),diag_precond(DIM), b(DIM);
    b.assign(rhs.begin(),rhs.end());// b = rhs;

    A.build_diag_precond<T>(diag_precond);
    iter.set_rhsnorm(norm(b));
    r.assign(b.begin(),b.end());          // r = b;
    std::vector<T> v_temp(x.size());
    mult(A,x,v_temp);                     // v_temp = A x;
    sub(v_temp,r);                        // r -= v_temp; donc r = b - A x;

    p_direct(diag_precond,r,z);           //mult(P, r, z);
    rho = dot(z,r);                       //rho = vect_sp(z, r);
    p.assign(z.begin(),z.end());          //copy(z, p);
    while (!iter.finished(norm(r)))
        {
        if (!iter.first())
            {
            p_direct(diag_precond,r,z);   //mult(P, r, z);
            rho = dot(z,r);
            scaled(rho/rho_1,p);          // p *= (rho/rho1)
            add(z,p);                     // p += z, so p = z + (rho/rho_1)*p
            }
        mult(A, p, q);

        T q_dot_p = dot(q,p);
        if (q_dot_p == 0.0) { std::cout<< "cg cannot converge: q orthogonal to p."; exit(1); }
        T a=rho/q_dot_p;
        scaled_add(p, +a, x);             // add(scaled(p, +a), x);
        scaled_add(q, -a, r);             // add(scaled(q, -a), r);
        rho_1 = rho;
        ++iter;
        }
    return norm(r)/norm(b);
    }

/** conjugate gradient with diagonal preconditioner with Dirichlet conditions (through masking technique), returns residu */
template<typename T>
T cg_dir(iteration<T> &iter, r_sparseMat& A, std::vector<T> & x, const std::vector<T> & rhs, const std::vector<T> & xd, const std::vector<int>& ld)
    {
    T rho, rho_1(0.0);
    const size_t DIM = x.size();
    if (rhs.size()!=DIM){std::cout << "rhs size mismatch" << std::endl; exit(1);}

    std::vector<T> p(DIM),q(DIM),r(DIM),z(DIM),diag_precond(DIM), b(DIM);
    b.assign(rhs.begin(),rhs.end());// b = rhs;

    A.build_diag_precond<T>(diag_precond);
    mult(A, xd, z);
    sub(z, b);                        // b -= A xd
    applyMask(ld,b);
    applyMask(ld,diag_precond);

    iter.set_rhsnorm(norm(b));
    r.assign(b.begin(),b.end());          // r = b;
    std::vector<T> v_temp(x.size());
    mult(A,x,v_temp);                     // v_temp = A x;
    sub(v_temp,r);                        // r -= v_temp; donc r = b - A x;

    applyMask(ld,r);

    p_direct(diag_precond,r,z);           //mult(P, r, z);
    rho = dot(z,r);                       //rho = vect_sp(z, r);
    p.assign(z.begin(),z.end());          //copy(z, p);
    while (!iter.finished(norm(r)))
        {
        if (!iter.first())
            {
            p_direct(diag_precond,r,z);   //mult(P, r, z);
            rho = dot(z,r);
            scaled(rho/rho_1,p);          // p *= (rho/rho1)
            add(z,p);                     // p += z, so p = z + (rho/rho_1)*p
            }
        mult(A, p, q);

        applyMask(ld,q);

        T q_dot_p = dot(q,p);
        if (q_dot_p == 0.0) { std::cout<< "cg cannot converge: q orthogonal to p."; exit(1); }
        T a=rho/q_dot_p;
        scaled_add(p, +a, x);             // add(scaled(p, +a), x);
        scaled_add(q, -a, r);             // add(scaled(q, -a), r);
        rho_1 = rho;
        ++iter;
        }
    add(xd, x);                   // x += xd
    return norm(r)/norm(b);
    }

} // end namespace alg
#endif
