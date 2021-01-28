#include <boost/progress.hpp>

#include "fem.h"

#include "config.h"

//#include "gmm/gmm.h"
//#include "gmm/gmm_MUMPS_interface.h"
//#include "gmm/gmm_precond_ilu.h"
#include "gmm/gmm_iter.h"
#include "gmm/gmm_solver_bicgstab.h"

#include "mesh.h"

#include "tetra.h"
#include "facette.h"
#include "tiny.h"

class electrostatSolver {
public:
    inline electrostatSolver(const std::shared_ptr<Nodes::Node[]> _p_node /**< [in] pointer to the nodes */,const int max_iter): refNode(_p_node), MAXITER(max_iter) {}

private:
    const std::shared_ptr<Nodes::Node[]> refNode;/**< direct access to the Nodes */
    
    const int MAXITER;//5000
    
    inline double getSigma(const int reg) const {return 0;}
    inline double getJn(const int reg) const {return 0;}
    inline double getV(const int reg) const {return 0;}
    
    template <class T>
    void assembling(T const& obj, gmm::dense_matrix <double> &Ke, std::vector <double> &Le, write_matrix &K, write_vector &L)
    {
    const int N = obj.getN();
    for (int ie=0; ie<N; ie++)
        {
        int i= obj.ind[ie];             
        for (int je=0; je<N; je++)
            {
            int j= obj.ind[je];
            K(i, j)+= Ke(ie, je);
            }
        L[i]+= Le[ie];
        }
    }

    void integrales(Tetra::Tet &tet, gmm::dense_matrix <double> &AE)
    { //sigma is the region conductivity
    double sigma = getSigma(tet.getRegion());

    for (int npi=0; npi<Tetra::NPI; npi++)
        {
        double w, dai_dx, dai_dy, dai_dz, daj_dx, daj_dy, daj_dz;
        w = tet.weight[npi];

        for (int ie=0; ie<Tetra::N; ie++)
            {
            dai_dx= tet.dadx[ie][npi];  dai_dy= tet.dady[ie][npi];  dai_dz= tet.dadz[ie][npi];

            for (int je=0; je<Tetra::N; je++)
                {
                daj_dx= tet.dadx[je][npi];  daj_dy= tet.dady[je][npi];  daj_dz= tet.dadz[je][npi];
                AE(ie, je) += sigma*(dai_dx*daj_dx + dai_dy*daj_dy + dai_dz*daj_dz)*w;		
                }
            }
        }
    }

void integrales(Facette::Fac &fac, std::vector <double> &BE)
{
double Jn = getJn(fac.getRegion());

for (int npi=0; npi<Facette::NPI; npi++)
    {
    double w =fac.weight(npi);
    for (int ie=0; ie<Facette::N; ie++) { BE[ie] -= Facette::a[ie][npi]*Jn *w; }
    }
}

int solve(mesh &msh)
{
    const int NOD = msh.getNbNodes();
    const int TET = msh.getNbTets();
    const int FAC = msh.getNbFacs();
    
boost::timer time;
const int  VERBOSE = 0;

write_matrix Kw(NOD, NOD);
write_vector Lw(NOD);
write_vector Xw(NOD);

std::cout << "assembling..." << std::endl;
time.restart();

for (int ne=0; ne<TET; ne++){
    Tetra::Tet &tet = msh.tet[ne];
    gmm::dense_matrix <double> K(Tetra::N, Tetra::N);
    std::vector <double> L(Tetra::N);
    integrales(tet, K);
    assembling<Tetra::Tet>(tet, K, L, Kw, Lw);
    }

for (int ne=0; ne<FAC; ne++){
    Facette::Fac &fac = msh.fac[ne];
    gmm::dense_matrix <double> K(Facette::N, Facette::N);
    std::vector <double> L(Facette::N);
    integrales(fac, L);     
    assembling<Facette::Fac>(fac, K, L, Kw, Lw);
    }

std::cout << time.elapsed() << std::endl;

time.restart();

read_matrix  Kr(NOD, NOD);    gmm::copy(Kw, Kr);

// conditions de Dirichlet
for (int ne=0; ne<FAC; ne++)
    {
    Facette::Fac &fac = msh.fac[ne];
    const int reg = fac.getRegion();
    
    double V = getV(reg);
       for (int ie=0; ie<Facette::N; ie++)
            {
           int i= fac.ind[ie];
           gmm::linalg_traits<gmm::row_matrix<write_vector > >::const_sub_row_type row = mat_const_row(Kw, i);
           for (write_vector::const_iterator it=vect_const_begin(row); it!=vect_const_end(row); ++it)
               Kw(i, it->first)=0;
           Kw(i, i) =  1;
           Lw[i]    =  V;  
           Xw[i]    =  V;
           }
    }

/* equilibrage des lignes */
for (int i=0; i<NOD; i++){
    gmm::linalg_traits<gmm::row_matrix<write_vector > >::const_sub_row_type row = mat_const_row(Kw, i);
    double norme=gmm::vect_norminf(row);
    Lw[i]/=norme;
    for (write_vector::const_iterator it=vect_const_begin(row); it!=vect_const_end(row); ++it)
        Kw(i, it->first)/=norme;
    }

gmm::copy(Kw, Kr);
read_vector  Lr(NOD);        gmm::copy(Lw, Lr);

std::cout << time.elapsed() << std::endl;

std::cout << "preconditionning ..." << std::endl;
time.restart();

gmm::diagonal_precond <read_matrix> prc(Kr);

//gmm::ilu_precond <read_matrix> prc(Kr);
std::cout << time.elapsed() << std::endl;

std::cout << "solving ..." << std::endl;
time.restart();

gmm::iteration iter(1e-8);
iter.set_maxiter(MAXITER);
iter.set_noisy(VERBOSE);

gmm::bicgstab(Kr, Xw, Lr, prc, iter);

read_vector Xr(NOD); gmm::copy(Xw, Xr);

std::cout << time.elapsed() << std::endl;

for (int i=0; i<NOD; i++) { refNode[i].V = Xr[i]; }

return 0;
}
    
    
}; //end class electrostatSolver
