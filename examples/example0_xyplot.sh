#!/bin/bash

# This script regenerates the figures of the user note
# user-notes/fg-xyplot.md with fg-xyplot, from the output of example0.py
# (a uniformly magnetized ellipsoid relaxing in a field along x).
#
# Usage, from any directory:
#     python3 example0.py      # in examples/, once
#     ./example0_xyplot.sh     # writes user-notes/figures/fg-xyplot_*.png

set -e

examples_dir=$(dirname "$(readlink -f "$0")")
xyplot="$examples_dir/../tools/fg-xyplot"
figures="$examples_dir/../user-notes/figures"

cd "$examples_dir"
for f in ellipsoid.msh ellipsoid_iter0.sol ellipsoid_iter100.sol ellipsoid_iter200.sol; do
    if [ ! -f "$f" ]; then
        echo "$f not found: run 'python3 example0.py' in $examples_dir first" >&2
        exit 1
    fi
done
mkdir -p "$figures"

all_sol="ellipsoid_iter0.sol ellipsoid_iter100.sol ellipsoid_iter200.sol"

# Magnetization of all the nodes against z.
"$xyplot" ellipsoid.msh ellipsoid_iter200.sol -x z -y mx my mz \
    -o "$figures/fg-xyplot_m_z_nodes.png"

# Profiles averaged over 8 bins of z, at three times.
"$xyplot" ellipsoid.msh $all_sol -x z -y mx my -b 8 \
    -o "$figures/fg-xyplot_mxmy_z_times.png"

# Scalar potential near the axis: filter in nm (units of the mesh), z in m.
"$xyplot" ellipsoid.msh ellipsoid_iter200.sol -x z -y phi -f 'hypot(x, y) < 1.5' \
    -u 1e-9 --xlabel 'z (m)' --ylabel phi \
    --title 'scalar potential at less than 1.5 nm from the z axis' \
    -o "$figures/fg-xyplot_phi_axis.png"

# my against mx: one cloud of points per time step.
"$xyplot" ellipsoid.msh $all_sol -x mx -y my \
    -o "$figures/fg-xyplot_my_mx.png"

# Radial component of m against the azimuthal angle, near the surface.
"$xyplot" ellipsoid.msh ellipsoid_iter200.sol -x 'arctan2(y, x)' \
    -y '(mx*x + my*y)/hypot(x, y)' -f 'hypot(x, y) > 3' -o "$figures/fg-xyplot_mr_angle.png"
