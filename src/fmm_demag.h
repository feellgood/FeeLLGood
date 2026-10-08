#ifndef FMM_DEMAG_H
#define FMM_DEMAG_H

/** \file fmm_demag.h
\brief this header is the interface to scalfmm 3. Its purpose is to prepare a source tree and a
target tree for the application of the fast multipole algorithm, and to compute the scalar magnetic
potential and the demagnetizing field.
scalfmm 3 headers define non inline functions, hence they must be included in a single translation
unit: they are only included by fmm_demag.cpp, and class fmm hides them behind a pointer to its
implementation.
*/

#include <memory>
#include <vector>

#include "mesh.h"

/** \namespace scal_fmm
to grab altogether the templates and functions using scalfmm for the computation of the demag field
*/

namespace scal_fmm
    {
/** \class fmm
to initialize the trees and the operators for the computation of the demagnetizing field, and launch
the computation easily with calc_demag public member
*/
class fmm
    {
public:
    /** constructor, initialize memory for trees, operators, sources corrections, initialize all
     * sources and targets
     */
    fmm(Mesh::mesh &msh /**< [in] */,
        std::vector<Tetra::prm> & prmTet /**< [in] */,
        std::vector<Triangle::prm> & prmTri /**< [in] */,
        const int ScalfmmNbThreads /**< [in] */,
        const int order /**< [in] order of the interpolation polynomials of the far field */,
        const int treeHeight /**< [in] height of the trees, 0 means automatic */,
        const int groupSize /**< [in] number of leaves and cells per group, 0 means automatic */,
        const bool verbose /**< [in] if true, print the durations of the steps */);

    /** destructor, defined where the implementation is complete */
    ~fmm();

    /**
    Compute the demagnetizing field. Include the second order corrections if FIRST_ORDER=OFF
    (which is the default).
    */
    void calc_demag(Mesh::mesh &msh /**< [in] */);

    /** height of the trees, root included (useful when it is computed automatically) */
    int treeHeight() const;

    /** group size of the source tree (useful when it is computed automatically) */
    int sourceGroupSize() const;

    /** group size of the target tree (useful when it is computed automatically) */
    int targetGroupSize() const;

private:
    /** implementation, see fmm_demag.cpp */
    struct impl;

    /** pointer to implementation */
    std::unique_ptr<impl> pImpl;
    };  // end class fmm

    }  // namespace scal_fmm
#endif
