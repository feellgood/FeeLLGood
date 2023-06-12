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

#include "ANN.h"
#include "feellgoodSettings.h"

/** \class Fem
class container to grab altogether all parameters of a simulation, including mesh geometry,
containers for the mesh
*/
class Fem
    {
public:
    /** constructor: call mesh constructor, initialize pts,kdtree and many inner variables */
    inline Fem(Settings const &mySets, timing &t_prm) : msh(mySets)
        {
        vmax = 0.0;

        E_exch0 = E_exch = 0.0;
        E_aniso0 = E_aniso = 0.0;
        E_demag0 = E_demag = 0.0;
        E_zeeman0 = E_zeeman = 0.0;
        Etot0 = INFINITY;  // avoid "WARNING: energy increased" on first time step
        Etot = 0.0;

        recenter_mem = false;
        if (mySets.recenter)
            {
            if (mySets.verbose)
                {
                std::cout << "Approximate nearest neighbors:\n";
                }
            pts = annAllocPts(msh.getNbNodes(), Pt::DIM);
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
                this->pts[i][0] = msh.getNode(i).p.x();
                this->pts[i][1] = msh.getNode(i).p.y();
                this->pts[i][2] = msh.getNode(i).p.z();
                }

            if (mySets.verbose)
                {
                std::cout << "  building kd_tree\n";
                }

            kdtree = new ANNkd_tree(pts, msh.getNbNodes(), Pt::DIM);
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

        if (mySets.recenter)
            {
            direction(Pt::IDX_Z);
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

    double vmax; /**< maximum speed of magnetization */

    double E_exch0;   /**< previous iteration exchange energy  */
    double E_aniso0;  /**< previous iteration anisotropy energy  */
    double E_demag0;  /**< previous iteration demagnetizing energy  */
    double E_zeeman0; /**< previous iteration zeeman energy  */

    double E_exch;   /**< exchange energy */
    double E_aniso;  /**< anisotropy energy  */
    double E_demag;  /**< demagnetizing energy  */
    double E_zeeman; /**< zeeman energy  */

    double DW_dir; /**< direction of the domain wall */

    double Etot0; /**< initial total energy */
    double Etot;  /**< total energy */

    /** mesh object to store nodes, fac, tet, and others geometrical values related to the mesh */
    Mesh::mesh msh;

    /** computes all the energies */
    void energy(double const t /**< [in] time in second, used to compute zeeman contribution if
                                  applied field is time dependant */
                ,
                Settings &settings /**< [in] */);

    /**
    time evolution : one step in time
    */
    inline void evolution(void)
        {
        msh.evolution();
        E_exch0 = E_exch;
        E_aniso0 = E_aniso;
        E_demag0 = E_demag;
        E_zeeman0 = E_zeeman;
        Etot0 = Etot;
        }

    /** saving function for a solution */
    void saver(Settings &settings /**< [in] */, timing const &t_prm /**< [in] */,
               std::ofstream &fout /**< [out] */, const int nt /**< [in] */) const;

    /** recentering algorithm for the study of the motion of a domain wall.

    if \f$ D_i>0 \f$				        if  \f$ D_i<0 \f$

    <----------------|------->		------->|<----------------	\f$ m_i = < u_i > < 0 \f$

    or					or

    ---------------->|<-------		<-------|---------------->	\f$ m_i = <u_i> > 0 \f$
    */

    bool recenter(double thres /**< [in] threshold parameter */,
                  char recentering_direction /**< [in] X|Y|Z */);

private:
    bool recenter_mem;  /**< flag to know if kdtree and pts are allocated */
    ANNkd_tree *kdtree; /**< ANN kdtree to find efficiently the closest set of nodes to a physical
                           point in the mesh  */
    ANNpointArray pts;  /**< container for the building of the kdtree (handled by ANN library) */

    /** find direction of motion of DW */
    void direction(enum Pt::index idx_dir /**< [in] */);

    /** zeroing of all energies */
    inline void zeroEnergy(void)
        {
        E_exch = 0.0;
        E_aniso = 0.0;
        E_demag = 0.0;
        E_zeeman = 0.0;
        Etot = 0.0;
        }

    /** computes the sum of all energies */
    inline void calc_Etot(void) { Etot = E_exch + E_aniso + E_demag + E_zeeman; }

    };  // end class

#endif
