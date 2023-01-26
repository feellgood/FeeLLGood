/** \file electrostatSolver.h
  \brief solver for electrostatic problem when STT is required
  header containing electrostatSolver class. It uses biconjugate stabilized gradient with diagonal preconditioner. The solver is only called once to compute voltages V for each nodes of the mesh, when STT computation is involved.
 */

#include <map>
#include <iostream>

#include "gmm/gmm_iter.h"
#include "gmm/gmm_solver_bicgstab.h"
#include "gmm/gmm_solver_cg.h"


#include "fem.h"
#include "config.h"
#include "mesh.h"

#include "tetra.h"
#include "facette.h"
#include "tiny.h"

#include "spinTransferTorque.h"

    /** assemble the matrix K from tet and Ke inputs */
    inline void assembling_mat(Tetra::Tet const& tet, gmm::dense_matrix <double> const &Ke, write_matrix &K)
    {
    for (int ie=0; ie<Tetra::N; ie++)
        { for (int je=0; je<Tetra::N; je++) { K(tet.ind[ie], tet.ind[je]) += Ke(ie, je); } }
    }
    
    /** assemble the vector L from fac and Le inputs */
    inline void assembling_vect(Facette::Fac const& fac, std::vector <double> const &Le, write_vector &L)
    { for (int ie=0; ie<Facette::N; ie++) { L[ fac.ind[ie] ] += Le[ie]; } }
    
    
    /** compute side problem (electrostatic potential on the nodes) integrales for matrix coefficients,input from tet ; sigma is the region conductivity */
    void integrales(Tetra::Tet const& tet, double sigma, gmm::dense_matrix <double> &AE)
    {
    for (int npi=0; npi<Tetra::NPI; npi++)
            {
            double w = tet.weight[npi];

            for (int ie=0; ie<Tetra::N; ie++)
                {
                double dai_dx= tet.dadx[ie][npi];
                double dai_dy= tet.dady[ie][npi];
                double dai_dz= tet.dadz[ie][npi];

                for (int je=0; je<Tetra::N; je++)
                    {
                    double daj_dx= tet.dadx[je][npi];
                    double daj_dy= tet.dady[je][npi];
                    double daj_dz= tet.dadz[je][npi];
                    AE(ie, je) += sigma*(dai_dx*daj_dx + dai_dy*daj_dy + dai_dz*daj_dz)*w;		
                    }
                }
            }
    }
    
    /** compute integrales for vector coefficients, input from facette */
void integrales(Facette::Fac const& fac,double pot_val, std::vector <double> &BE)
    {
    for (int npi=0; npi<Facette::NPI; npi++)
        for (int ie=0; ie<Facette::N; ie++)
            { BE[ie] -= Facette::a[ie][npi]*pot_val*fac.weight(npi); }
    }

/** \class electrostatSolver
this class is containing both data and a solver to compute potential from dirichlet boundary conditions problem for the current density flowing in the sample.
*/
class electrostatSolver {
public:
    /** constructor */
    inline electrostatSolver(Mesh::mesh const& _msh /**< [in] reference to the mesh */,
                             STT const& _p_stt /**< all spin transfer torque parameters */,
                             const double _tol /**< [in] tolerance for solvers */,
                             const bool v /**< [in] verbose bool */,
                             const int max_iter /**< [in] maximum number of iteration */ ):
                             msh(_msh), p_stt(_p_stt), verbose(v), MAXITER(max_iter), precision(PRECISION_STT)
                             {
                             ksi = Pt::sq(p_stt.lJ/p_stt.lsf);
                             D0 = 2.0*p_stt.sigma/(Pt::sq(CHARGE_ELECTRON)*p_stt.N0);
                             pf = Pt::sq(p_stt.lJ)/(D0*(1.+ksi*ksi)) * BOHRS_MUB*p_stt.beta/CHARGE_ELECTRON;

                             if(verbose)
                             	{ std::cout << "Dirichlet boundary conditions..." << std::endl; infos(); }
                             bool has_converged = solve(_tol);
                             if(has_converged)
                             	{	
                             	if (p_stt.V_file) 
                             		{
                             		std::string fileName = "V.sol";
                             		if(verbose)
                                        { std::cout << "writing electrostatic potential solutions to file " << fileName << std::endl; }
                                    bool iznogood = msh.savesol(precision,fileName,"##columns: index\tV\n",V);
                                    if (verbose && iznogood)
                                        { std::cout << "file "<< fileName <<" status : " << iznogood <<std::endl; }
                             		}
                             	std::for_each(msh.tet.begin(),msh.tet.end(), [this]( Tetra::Tet const& tet )
                             		{
                                    std::array<Pt::pt3D,Tetra::NPI> _gradV;
                                    calc_gradV(tet,_gradV);
									gradV.push_back(_gradV);

									std::array<Pt::pt3D,Tetra::NPI> _Hm;
									calc_Hm(tet,_gradV,_Hm);
                             		Hm.push_back(_Hm);
                             		} );
                             	
                             	}
                             else
                                { std::cerr << "Solver (STT) has not converged" << std::endl; exit(1); }
                             }

/** computes the gradient(V) for tetra tet */
void calc_gradV(Tetra::Tet const& tet, std::array<Pt::pt3D,Tetra::NPI> & _gradV) const
    {
    for (int npi=0; npi<Tetra::NPI; npi++) 
    	{ 
    	double vx(0),vy(0),vz(0);
      	for (int i=0; i<Tetra::N; i++)
       	{
        	const double V_tet = V[ tet.ind[i] ];
        	vx += V_tet*tet.dadx[i][npi];
        	vy += V_tet*tet.dady[i][npi];
        	vz += V_tet*tet.dadz[i][npi];
        	}
    	_gradV[npi] = Pt::pt3D(vx,vy,vz);
    	}
    }

/** computes Hm contributions for each npi for tetrahedron tet */
void calc_Hm(Tetra::Tet const& tet, std::array<Pt::pt3D,Tetra::NPI> const& _gradV, std::array<Pt::pt3D,Tetra::NPI> & _Hm ) const
{
Pt::pt3D p_g[Tetra::NPI];
tet.interpolation(Nodes::get_p,p_g);

for (int npi=0; npi<Tetra::NPI; npi++)
	{ _Hm[npi] = -p_stt.sigma*_gradV[npi]*p_g[npi]; }
}

/** affect extraField function and extraCoeffs_BE function for all the tetrahedrons */
void prepareExtras(void)
	{
	std::for_each(msh.tet.begin(),msh.tet.end(), [this]( Tetra::Tet & tet )
		{
		const int _idx = tet.idx;
		tet.extraField = [this,_idx]( int npi, Pt::pt3D & _Hm ) { _Hm = this->Hm[_idx][npi]; } ;

		tet.extraCoeffs_BE = [this,&tet,_idx](int npi,double Js,Pt::pt3D &U,Pt::pt3D &dUdx,Pt::pt3D &dUdy,Pt::pt3D &dUdz, Pt::pt3D (&BE)[Tetra::N] )
			{
			const double prefactor = D0/Pt::sq(p_stt.lJ)/(gamma0*nu0*Js);

	        Pt::pt3D const& _gV = gradV[_idx][npi];

            Pt::pt3D j_grad_u = -p_stt.sigma*Pt::pt3D(Pt::pScal(_gV,Pt::pt3D(dUdx(Pt::IDX_X),dUdy(Pt::IDX_X),dUdz(Pt::IDX_X)) ),
                                 Pt::pScal(_gV,Pt::pt3D(dUdx(Pt::IDX_Y),dUdy(Pt::IDX_Y),dUdz(Pt::IDX_Y)) ),
                                 Pt::pScal(_gV,Pt::pt3D(dUdx(Pt::IDX_Z),dUdy(Pt::IDX_Z),dUdz(Pt::IDX_Z)) ));

            Pt::pt3D m = pf*(ksi*j_grad_u+U*j_grad_u);

            for (int i=0; i<Tetra::N; i++)
                { BE[i] += tet.weight[npi]*Tetra::a[i][npi]*(Hm[_idx][npi] + prefactor*m); }
			};
		});
	}

private:
    /** ksi is in Thiaville notations beta_DW */
    double ksi;

	/** density of states */
    double D0;

    /** a prefactor for BE coefficient coefficients*/
    double pf;
    
    /** mesh object to store nodes, fac, tet, and others geometrical values related to the mesh ( const ref ) */
    Mesh::mesh msh;

    /** spin transfer torque parameters */
    STT p_stt;

    /** if verbose set to true, some printing are sent to terminal */
    const bool verbose;

    /** maximum number of iteration for biconjugate stabilized gradient */
    const unsigned int MAXITER; //fixed to 5000 in ref code

    /** number of digits in the optional output file */
    const int precision;

    /** electrostatic potential values for boundary conditions, V.size() is the size of the vector of nodes */ 
    std::vector<double> V;

    /** table of the gradients of the potential, gradV.size() is the number of tetra */
    std::vector< std::array<Pt::pt3D,Tetra::NPI> > gradV;

    /** table of the Hm vectors (contribution of the STT to the tet::integrales) ; Hm.size() is the number of tetra */
    std::vector< std::array<Pt::pt3D,Tetra::NPI> > Hm;

    /** basic informations on boundary conditions */
    inline void infos(void) 
        {
        std::cout << "sigma: " << p_stt.sigma << std::endl;
        
        std::for_each(p_stt.boundaryCond.begin(),p_stt.boundaryCond.end(),
                      [](std::pair<std::string,double> const& p)
                      { std::cout << "regName: " << p.first << "\tV :" << p.second << std::endl; } );
        }
        
    /** fill matrix and vector to solve potential values on each node */
    void prepareData(write_matrix &Kw, write_vector & Lw)
        {
        const double sigma = p_stt.sigma;

        std::for_each(msh.tet.begin(),msh.tet.end(), [this,sigma,&Kw]( Tetra::Tet const& tet )
            {
            gmm::dense_matrix <double> K(Tetra::N, Tetra::N);
            integrales(tet,sigma, K);
            assembling_mat(tet, K, Kw);
            });

        double pot_val = 0; // we initialize pot_val to the average of the potentials set by the boundary conditions (does not seem to change convergence speed whatever value it is ..)
        std::for_each( p_stt.boundaryCond.begin(),p_stt.boundaryCond.end(), [&pot_val](auto const& it) { pot_val += it.second;} );
        pot_val /= p_stt.boundaryCond.size();

        std::for_each(msh.fac.begin(),msh.fac.end(), [this,pot_val,&Lw]( Facette::Fac const& fac )
		    {
		    std::vector <double> L(Facette::N);
		    integrales(fac, pot_val, L);
            assembling_vect(fac, L, Lw);
		    });
        }

/** solver, using biconjugate stabilized gradient, with diagonal preconditionner and Dirichlet boundary conditions */
int solve(const double iter_tol)
{
const int NOD = msh.getNbNodes();

write_matrix Kw(NOD, NOD);
write_vector Lw(NOD);
write_vector Xw(NOD);

prepareData(Kw,Lw);

read_matrix  Kr(NOD, NOD);    gmm::copy(Kw, Kr);

std::for_each( p_stt.boundaryCond.begin(),p_stt.boundaryCond.end(), [this,&Kw,&Lw,&Xw](auto const& it) 
	{ 
	double V = it.second;
	std::string name = it.first;
	
	auto surf_it = std::find_if(msh.s.begin(),msh.s.end(), [name] (Mesh::Surf const& _s) { return (_s.getName() == name); }  );
	
	if(surf_it != msh.s.end())
		{
		std::for_each( surf_it->elem.begin(),surf_it->elem.end(), [V,&Kw,&Lw,&Xw](Mesh::Triangle const& tri)
			{
			for (int ie=0; ie<Facette::N ; ie++) // should be Triangle::N
				{
				const int i= tri.ind[ie];
				std::for_each(mat_row(Kw, i).begin(),mat_row(Kw, i).end(), [](std::pair<const long unsigned int, double> & _it) { _it.second = 0.0; } );
            			Kw(i, i) =  1e9;
            			Lw[i]    =  V*1e9;  
            			Xw[i]    =  V;
				}
			});
		}
	} );

if(verbose)
    { std::cout << "line weighting..." << std::endl; }

for (int i=0; i<NOD; i++)
    {
    double norme=gmm::vect_norminf(mat_row(Kw, i));
    Lw[i]/=norme;
    std::for_each(mat_row(Kw, i).begin(),mat_row(Kw, i).end(), [norme](std::pair<const long unsigned int, double> & it) { it.second/=norme; } );
    }

gmm::copy(Kw, Kr);
read_vector  Lr(NOD);
gmm::copy(Lw, Lr);

if(verbose)
    { std::cout << "solving ..." << std::endl; }

gmm::iteration iter(iter_tol);
iter.set_maxiter(MAXITER);
iter.set_noisy(verbose);

gmm::bicgstab(Kr, Xw, Lr, gmm::diagonal_precond <read_matrix>(Kr), iter);

V.resize(NOD);
gmm::copy(Xw, V);
return (iter.get_iteration() < MAXITER);
}

}; //end class electrostatSolver
