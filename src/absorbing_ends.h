#ifndef ABSORBING_ENDS_H
#define ABSORBING_ENDS_H

/** \file absorbing_ends.h
\brief absorbing ends, to absorb the spin waves reaching the ends of a wire instead of reflecting
them (settings: absorbing_ends, see dev-notes/spin_waves_absorbing_damping_factor.md).
*/

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "node.h"

/** \struct AbsorbingEnds
Along an axis, over a length \f$ L \f$ at each end of the union of the selected volume regions (not
at the ends of each of them), the damping goes from the physical value \f$ \alpha_0 \f$ of the region
to \f$ \alpha_{\max} \f$:
\f$ \alpha = \alpha_0 + (\alpha_{\max}-\alpha_0) \sin^2\left(\frac{\pi}{2}\frac{s}{L}\right) \f$,
where \f$ s \f$ is the depth in the zone (0 at its inner boundary, \f$ L \f$ at the end of the
union), so that \f$ d\alpha/ds = 0 \f$ at both boundaries of the zone. The damping is evaluated at
the Gauss points of the LLG assembly.
*/
struct AbsorbingEnds
    {
    /** if false, the damping is the one of the regions everywhere */
    bool enable = false;

    /** axis of the wire */
    Nodes::index axis = Nodes::IDX_Z;

    /** length of each absorbing zone, in meters */
    double length = 0;

    /** damping at the two ends of the union of the selected regions */
    double alpha_max = 1;

    /** names of the selected volume regions (empty: all the volume regions) */
    std::vector<std::string> regions;

    /** lower bound of the union of the selected regions along the axis, in meters (set from the
     * mesh) */
    double lo = 0;

    /** upper bound of the union of the selected regions along the axis, in meters (set from the
     * mesh) */
    double hi = 0;

    /** damping at the coordinate x (along the axis, in meters) of a region of damping alpha0 */
    double alpha(const double alpha0, const double x) const
        {
        if (!enable)
            { return alpha0; }
        const double s = std::max(x - (hi - length), (lo + length) - x);  // depth in the nearest zone
        if (s <= 0)
            { return alpha0; }
        const double sn = std::sin(M_PI / 2 * std::min(s / length, 1.0));
        return alpha0 + (alpha_max - alpha0) * sn * sn;
        }
    };

#endif
