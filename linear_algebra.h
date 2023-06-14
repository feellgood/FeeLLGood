#ifndef linear_algebra_h
#define linear_algebra_h

/** \file linear_algebra.h
\brief secondary header, it grabs altogether the linear algebra by the solver to apply fem method
<br> It encapsulates the calls to GMM, the assemblage and projection of the matrix for all elements
<br> projection and matrix assembly is multithreaded for tetrahedron, monothread for facette
*/
#include <random>
#include <execution>

#include "gmm/gmm_iter.h"
#include "gmm/gmm_precond_diagonal.h"
#include "gmm/gmm_solver_bicgstab.h"
#include "gmm/gmm_solver_gmres.h"

#include "config.h"

#include "facette.h"
#include "feellgoodSettings.h"
#include "mesh.h"
#include "node.h"
#include "pt3D.h"
#include "tetra.h"

/** \class LinAlgebra
convenient class to grab altogether some part of the calculations involved using gmm solver at each
timestep
*/
class LinAlgebra
    {
public:
    /** constructor */
    inline LinAlgebra(Settings &s /**< [in] */, Mesh::mesh &my_msh /**< [in] */)
        : NOD(my_msh.getNbNodes()), refMsh(&my_msh), settings(s)
        {
        base_projection();
        }

    /** destructor */
    ~LinAlgebra()
        {
        if (prc != nullptr)
            {
            delete prc;
            }
        }

    /** pointer to diagonal preconditioner  */
    gmm::diagonal_precond<read_matrix> *prc = nullptr;

    /** computes inner data structures of tetraedrons and triangular facettes (K matrices and L
     * vectors) */
    void prepareElements(Pt::pt3D const &Hext /**< [in] applied field */,
                         timing const &t_prm /**< [in] */);

    /**  solver, uses bicgstab and gmres, sparse matrix and vector are filled with multiThreading */
    int solver(timing const &t_prm /**< [in] */, long nt /**< [in] */);

    /** setter for DW_dz */
    inline void set_DW_vz(double vz /**< [in] */) { DW_vz = vz; }

    /** getter for v_max */
    inline double get_v_max(void) { return v_max; }

private:
    const int NOD; /**< total number of nodes, also an offset for filling sparseMatrix, initialized
                      by constructor */
    Mesh::mesh *refMsh; /**< direct access to the mesh */

    double DW_vz;             /**< speed of the domain wall */
    const Settings &settings; /**< settings */
    double v_max;             /**< maximum speed */

    long prc_time_step = -1; /**< time step when prc was built */

    /** update nodes.u and v values and find the maximum speed value  */
    void updateNodes(std::vector<double> const &X, const double dt);

    /** computes the local vector basis {ep,eq} in the tangeant plane for projection on the elements
     */
    void base_projection();

    /** template to assemble the big sparse matrix K with T tetra or facette */
    template<class T>
    void assemblage_mat(T &elem, write_matrix &K) const
        {
        const int N = elem.getN();
        for (int i = 0; i < N; i++)
            {
            int i_ = elem.ind[i];

            for (int j = 0; j < N; j++)
                {
                int j_ = elem.ind[j];
                K(NOD + i_, j_) += elem.Kp[i][j];
                K(NOD + i_, NOD + j_) += elem.Kp[i][N + j];
                K(i_, j_) += elem.Kp[N + i][j];
                K(i_, NOD + j_) += elem.Kp[N + i][N + j];
                }
            }
        }

    /** template to assemble the big vector L with T tetra or facette */
    template<class T>
    void assemblage_vect(T &elem, std::vector<double> &L) const
        {
        const int N = elem.getN();
        for (int i = 0; i < N; i++)
            {
            const int i_ = elem.ind[i];
            L[NOD + i_] += elem.Lp[i];
            L[i_] += elem.Lp[N + i];
            }
        }

    /** template to insert coeff in sparse matrix K_TH and vector L_TH, T is Tetra or Facette */
    template<class T>
    void insertCoeff(std::vector<T> &container, write_matrix &K_TH, std::vector<double> &L_TH)
        {
        std::for_each(container.begin(), container.end(),
                      [this,&K_TH, &L_TH](T &my_elem)
                      {
                          assemblage_mat<T>(my_elem,K_TH);
                          assemblage_vect<T>(my_elem,L_TH);
                      });
        }

    };  // fin class linAlgebra

#endif
