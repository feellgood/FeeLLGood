#ifndef fem_h
#define fem_h

/** \file fem.h
\brief principal header, contains the struct fem
This file is called by many source files in feellgood.
It does also contains the definition of many constants for the solver, and for scalfmm
*/
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <sstream>
#include <utility>
#include <vector>

#include <sys/times.h>
#include <time.h>
#include <unistd.h>

#include "mesh.h"
#include "fmm_demag.h"
#include "spinAccumulationSolver.h"
#include "linear_algebra.h"
#include "chronometer.h"
#include "ANN.h"
#include "feellgoodSettings.h"

/** number of energy terms in energy array E */
const int NB_ENERGY_TERMS = 4;

/** correspondance nature/index of the energy terms in energy array E */
enum ENERGY_TYPE
    {
    EXCHANGE = 0,
    ANISOTROPY = 1,
    DEMAG = 2,
    ZEEMAN = 3
    };

/** \class Fem
class container to grab altogether all parameters of a simulation, including mesh geometry,
containers for the mesh
*/
class Fem
    {
public:
    /** constructor: call mesh constructor, initialize pts,kdtree and many inner variables */
    inline Fem(Settings &mySets, timing &t_prm) : msh(mySets)
        {
        vmax = 0.0;
        std::fill(E.begin(),E.end(),0);
        Etot0 = INFINITY;  // avoid "WARNING: energy increased" on first time step
        Etot = 0.0;

        recenter_mem = false;
        if (mySets.recenter)
            {
            if (mySets.verbose)
                {
                std::cout << "Approximate nearest neighbors:\n";
                }
            pts = annAllocPts(msh.getNbNodes(), Nodes::DIM);
            if (!pts)
                {
                std::cout << "ANN memory error while allocating points" << std::endl;
                SYSTEM_ERROR;
                }
            else if (mySets.verbose)
                {
                std::cout << "  points allocated\n";
                }

            for (int i = 0; i < msh.getNbNodes(); i++)
                {
                this->pts[i][0] = msh.getNode_p(i).x();
                this->pts[i][1] = msh.getNode_p(i).y();
                this->pts[i][2] = msh.getNode_p(i).z();
                }

            if (mySets.verbose)
                {
                std::cout << "  building kd_tree\n";
                }

            kdtree = new ANNkd_tree(pts, msh.getNbNodes(), Nodes::DIM);
            if (!kdtree)
                {
                std::cout << "ANN memory error while allocating kd_tree" << std::endl;
                SYSTEM_ERROR;
                }
            recenter_mem = true;

            if (mySets.verbose)
                {
                std::cout << "  kd_tree allocated\n";
                }
            }
        else if (mySets.verbose)
            {
            std::cout << "No recentering.\n";
            }

        if (mySets.restoreFileName == "")
            {
            msh.init_distrib(mySets);
            }
        else
            {
            t_prm.set_t(msh.readSol(mySets.verbose, mySets.restoreFileName));
            }

        /* This potentially overrides the initial time set above by msh.readSol(). */
        if (!isnan(mySets.initial_time))
            {
            t_prm.set_t(mySets.initial_time);
            }

        if (mySets.recenter)
            {
            direction(Nodes::IDX_Z);
            } /* DW propagation direction for recentering */
        }

    /** destructor */
    ~Fem()
        {
        if (recenter_mem)
            {
            annDeallocPts(pts);
            delete kdtree;
            }
        }

    /** maximum speed of magnetization */
    double vmax;

    /** current iteration energies */
    std::array<double,NB_ENERGY_TERMS> E;

    /** direction of the domain wall */
    double DW_dir;

    /** initial total energy */
    double Etot0;
    
    /** total energy */
    double Etot;

    /** mesh object to store nodes, fac, tet, and others geometrical values related to the mesh */
    Mesh::mesh msh;

    /** computes all the energies */
    void energy(double const t /**< [in] time in second, used to compute zeeman contribution if
                                  applied field is time dependant */
                ,Settings &settings /**< [in] */);

    /**
    time evolution : one step in time
    */
    inline void evolution(void)
        {
        msh.evolution();
        Etot0 = Etot;
        }

    /** saving function for a solution */
    void saver(Settings &settings /**< [in] */, timing const &t_prm /**< [in] */,
               std::ofstream &fout /**< [out] */, const int nt /**< [in] */,
               std::vector<Eigen::Vector3d> &s /**< [in] */ ) const;

    /** recentering algorithm for the study of the motion of a domain wall.

    if \f$ D_i>0 \f$				        if  \f$ D_i<0 \f$

    <----------------|------->		------->|<----------------	\f$ m_i = < u_i > < 0 \f$

    or					or

    ---------------->|<-------		<-------|---------------->	\f$ m_i = <u_i> > 0 \f$
    */

    bool recenter(double thres /**< [in] threshold parameter */,
                  char recentering_direction /**< [in] X|Y|Z */);

    /** performs the time integration of the LLG problem, with spin accumulation if any */
    int time_integration(Settings &settings /**< [in] */,
                         LinAlgebra &linAlg /**< [in] */,
                         spinAcc &spinAcc_solver /**< [in] */,
                         scal_fmm::fmm &myFMM /**< [in] */,
                         timing &t_prm, int &nt /**< [out] number of time steps performed */);

private:
    bool recenter_mem;  /**< flag to know if kdtree and pts are allocated */
    ANNkd_tree *kdtree; /**< ANN kdtree to find efficiently the closest set of nodes to a physical
                           point in the mesh  */
    ANNpointArray pts;  /**< container for the building of the kdtree (handled by ANN library) */

    /** find direction of motion of DW */
    void direction(enum Nodes::index idx_dir /**< [in] */);

    /** compute demagnetizing field, energies, and prepare for next time step quantitites at time t */
    void compute_all(Settings &settings, spinAcc &spinAcc_solver /**< [in] */,
                     scal_fmm::fmm &myFMM, const double t)
        {
        bool success(false);
        chronometer fmm_counter(2);
        myFMM.calc_demag(msh);
        if (settings.spin_acc && spinAcc_solver.compute())
            { success = true; }
        if (settings.verbose)
            {
            std::cout << "magnetostatics done in " << fmm_counter.millis() << std::endl;
            if (settings.spin_acc)
                {
                if (success)
                    { std::cout << "spin diffusion solved.\n"; }
                else
                    { std::cout << "spin diffusion solver failed.\n"; }
                }
            }
        energy(t, settings);
        evolution();
        }
    };  // end class

#endif
