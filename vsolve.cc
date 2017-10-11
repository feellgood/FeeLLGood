#include "fem.h"

int vsolve(Fem &fem, long nt)
{
boost::timer time;
const int  VERBOSE = 0;
const int NOD = fem.NOD;
const int TET = fem.TET;
const int FAC = fem.FAC;

const int MAXITER = 500;
const int REFRESH_PRC = 20;

base_projection(fem);   // definit plan tangent

write_matrix Kw(2*NOD, 2*NOD);
write_vector Lw(2*NOD);

/* changement de referentiel */
/* bcarvello, 2017: doesn't this belong in evolution.cc ? */
fem.DW_vz += fem.DW_dir*v_moy(fem, 2)*fem.lz/2.;
IF_VERBOSE(fem){
cout << boost::format("%5t average velocity %30T.") <<flush;
cout << fem.DW_vz << endl;

//cout << " matrix size " << 2*NOD << endl;
cout << boost::format("%5t assembling %30T.");
}

time.restart();

for (int t=0; t<TET; t++){
    Tet &tet = fem.tet[t];
    gmm::dense_matrix <double> K(3*Tet::N,3*Tet::N), Kp(2*Tet::N,2*Tet::N);
    vector <double> L(3*Tet::N), Lp(2*Tet::N);
//    gmm::clear(K), gmm::clear(Kp);
    integrales(fem, tet, K, L);     
    projection(fem, tet, K, L, Kp, Lp);
    assemblage(fem, tet, Kp, Lp, Kw, Lw);
    }

for (int t=0; t<FAC; t++){
    Fac &fac = fem.fac[t];
    gmm::dense_matrix <double> K(3*Fac::N,3*Fac::N), Kp(2*Fac::N,2*Fac::N);
    vector <double> L(3*Fac::N), Lp(2*Fac::N);
//    gmm::clear(K), gmm::clear(Kp);
    integrales(fem, fac, K, L);     
    projection(fem, fac, K, L, Kp, Lp);
    assemblage(fem, fac, Kp, Lp, Kw, Lw);
    }

IF_VERBOSE(fem) cout << time.elapsed() << endl;

read_matrix Kr(2*NOD,2*NOD);    gmm::copy(Kw, Kr);
read_vector Lr(2*NOD);          gmm::copy(Lw, Lr);
write_vector Xw(2*NOD);

gmm::iteration bicg_iter(1e-6);
gmm::iteration gmr_iter(1e-6);
bicg_iter.set_maxiter(MAXITER);
gmr_iter.set_maxiter(MAXITER);
bicg_iter.set_noisy(VERBOSE);
gmr_iter.set_noisy(VERBOSE);

time.restart();

//   gmm::identity_matrix Precond; 
//   gmm::diagonal_precond <read_matrix> Precond(Kr);
//   gmm::mr_approx_inverse_precond <read_matrix> Precond(Kr, 10, 1e-30);
//   gmm::ilu_precond <read_matrix> Precond(Kr);
//   gmm::ilutp_precond < read_matrix > Precond(Kr,10,1e-30);

if (!nt) {
    IF_VERBOSE(fem){
    cout << boost::format("%5t computing prc %30T.");
    fflush(NULL);
    }
    fem.prc = new gmm::diagonal_precond <read_matrix> (Kr);
//    fem.prc = new gmm::ilutp_precond < read_matrix > (Kr,1,1e-30);
    IF_VERBOSE(fem) cout << time.elapsed() << endl;
    }
else if (!(nt % REFRESH_PRC)) {
    delete fem.prc;
    IF_VERBOSE(fem){
    cout << boost::format("%5t computing prc %30T.");
    fflush(NULL);
    }
    fem.prc = new gmm::diagonal_precond <read_matrix> (Kr);
//    fem.prc = new gmm::ilutp_precond < read_matrix > (Kr,1,1e-30);
    IF_VERBOSE(fem) cout << time.elapsed() << endl;
    } 

time.restart();
gmm::bicgstab(Kr, Xw, Lr, *fem.prc, bicg_iter);

if (!(bicg_iter.converged())) {
    IF_VERBOSE(fem) cout << boost::format("%5t  bicg FAILED in %d %30T. ") % bicg_iter.get_iteration() 
         << time.elapsed() << endl;
    gmm::diagonal_precond <read_matrix>  gmr_prc (Kr);
    time.restart();
    gmm::clear(Xw);
    gmm::gmres(Kr, Xw, Lr, gmr_prc, 50, gmr_iter);
    if (!(gmr_iter.converged())) {
    IF_VERBOSE(fem) cout << boost::format("%5t gmres FAILED in ")
         << gmr_iter.get_iteration() << " !! ......... "
         << time.elapsed() << endl;
    return 1;
    }
    else {
        IF_VERBOSE(fem) cout << boost::format("%5t v-solve in ") << gmr_iter.get_iteration() 
             << " (gmres) ............. " << time.elapsed() << endl;
        }
    }
else {
    IF_VERBOSE(fem) cout << boost::format("%5t v-solve in ") << bicg_iter.get_iteration()
         << " (bicg,prc-" << (nt % REFRESH_PRC) << ")" 
         << " .......... " << time.elapsed() << endl;
    }

read_vector Xr(2*NOD);    gmm::copy(Xw, Xr);

double v2max = 0.0;

for (int i=0; i<NOD; i++) {
    double vp,vq,v2;
    vp = Xr[i]; vq = Xr[NOD+i];
    v2 = vp*vp + vq*vq;
    if (v2>v2max)
        v2max = v2;

    Node &node = fem.node[i];
    double u2 = 0.0;
    for (int d=0; d<3; d++) {
        node.v[d]  = vp*node.ep[d] + vq*node.eq[d];
        node.u[d] = node.u0[d] + fem.dt * node.v[d];
        u2+= sq (node.u[d]);
        }
    for (int d=0; d<3; d++)
        node.u[d]/= sqrt(u2);
    }

fem.vmax = sqrt(v2max);
//printf("vmax = %g \n", fem.vmax);

return 0;
}
