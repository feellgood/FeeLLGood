#!/usr/bin/env python3
"""
Validation of the tunnel barrier (key `tunnel` of the surface regions, src/tunnelJunction.h) with a
non-uniform magnetization: the current concentrates in the domain parallel to the magnetization of
F1. Checked against a 2D reference written independently of feeLLGood.

Sample: slab F2 (two_domains.geo, Lx x Ly x t), current along z, tunnel barrier on the inlet
(z = 0), outlet at V = 0 (z = t). The magnetization of F1 is p = z; the one of F2 is a Néel wall
centred at x = 0, m = (sin theta, 0, cos theta), theta = 2 atan(exp(x/w)): parallel to p for
x << -w, antiparallel for x >> w. The conductance per unit area of the barrier follows the cosine
law, G_T(m) = (1 + c p.m)/((1 + c) RA), c = TMR/(2 + TMR), and the current leaving F2 through the
barrier is jn = G_T(m) (V - V1) (Robin condition).

Reference: the problem does not depend on y; div(sigma grad V) = 0 is solved in the plane (x, z) by
second-order finite volumes on a fine grid (Neumann on x = +-Lx/2, V = 0 on z = t, Robin on z = 0),
with the exact m(x).

Two cases:
- local: RA = 1e-12 Ohm m^2, much larger than the resistance t/sigma of F2: jn(x) is the one of a
  column, -V1/(RA(theta(x)) + t/sigma); far from the wall, jn_P/jn_AP is close to 1 + TMR;
- coupled: RA = 2e-15 Ohm m^2, comparable to t/sigma = 5e-15: the current spreads laterally in F2
  over about t around the wall (7.7 % from the local model of the columns), only the 2D reference
  describes it.

Usage: two_domains.py [path to feellgood], from any directory; the mesh is built with gmsh if
missing. Exit status 0 if the errors are below the tolerance in both cases.
"""

import json
import os
import subprocess
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

TOL = 1e-2                           # relative tolerance on V and jn at the inlet

LX, T = 100e-9, 20e-9                # as in two_domains.geo (m)
W = 5e-9                             # width parameter of the wall (m)
F2 = {"sigma": 4e6, "N0": 1.344e47, "P": 0.5, "lsf": 10e-9, "lsd": 1e-9}
BARRIER = {"TMR": 1.5, "V": -0.1}    # V1 < 0: electrons flowing from F1 into the sample
CASES = {"local": 1e-12, "coupled": 2e-15}   # RA, Ohm m^2


def theta(x):
    """angle between m and p = z across the wall"""
    return 2 * np.arctan(np.exp(x / W))


def conductance(RA, cos):
    """conductance per unit area of the barrier, cosine law"""
    c = BARRIER["TMR"] / (2 + BARRIER["TMR"])
    return (1 + c * cos) / ((1 + c) * RA)


def reference(RA, dx=0.2e-9):
    """potential V(x, z) of the 2D problem on a grid of step dx, rows along x, columns along z"""
    nx, nz = round(LX / dx) + 1, round(T / dx) + 1
    x = np.linspace(-LX / 2, LX / 2, nx)
    G = conductance(RA, np.cos(theta(x)))
    s, V1 = F2["sigma"], BARRIER["V"]
    idx = lambda i, k: i * nz + k
    rows, cols, vals = [], [], []
    rhs = np.zeros(nx * nz)

    def add(r, c, v):
        rows.append(r); cols.append(c); vals.append(v)

    for i in range(nx):
        for k in range(nz):
            r = idx(i, k)
            if k == nz - 1:                       # outlet, V = 0
                add(r, r, 1.0)
                continue
            # finite volume of node (i, k): half volumes on the boundaries, flux balance
            hx = dx if 0 < i < nx - 1 else dx / 2
            hz = dx if k > 0 else dx / 2
            for j in (i - 1, i + 1):              # lateral fluxes (none through x = +-Lx/2)
                if 0 <= j < nx:
                    add(r, idx(j, k), s * hz / dx)
                    add(r, r, -s * hz / dx)
            add(r, idx(i, k + 1), s * hx / dx)    # flux from above
            add(r, r, -s * hx / dx)
            if k > 0:
                add(r, idx(i, k - 1), s * hx / dx)
                add(r, r, -s * hx / dx)
            else:                                 # inlet: jn = G (V - V1) leaves the volume
                add(r, r, -G[i] * hx)
                rhs[r] -= G[i] * hx * V1
    A = sp.csr_matrix((vals, (rows, cols)), shape=(nx * nz, nx * nz))
    return x, spla.spsolve(A, rhs).reshape(nx, nz)


def settings(name, RA):
    return {
        "outputs": {"file_basename": name, "final_time": 1e-13, "evol_time_step": 1e-12,
                    "mag_config_every": 10},
        "mesh": {
            "filename": "two_domains.msh",
            "length_unit": 1e-9,
            "volume_regions": {"F2": {
                "Js": 1.0, "Ae": 1.3e-11, "alpha_LLG": 0.05, "K": 0,
                "sigma": F2["sigma"], "dens_state": F2["N0"], "P": F2["P"],
                "l_sd": F2["lsd"], "l_sf": F2["lsf"]}},
            "surface_regions": {
                "outlet": {"V": 0.0, "s": [0, 0, 0]},
                "inlet": {"tunnel": {"p": [0, 0, 1], "RA": RA, "TMR": BARRIER["TMR"],
                                     "V": BARRIER["V"]}}}},
        # x in m (the coordinates given to the function are in m, whatever length_unit)
        "initial_magnetization": "function(x, y, z) { const th = 2*Math.atan(Math.exp(x/%g)); "
                                 "return [Math.sin(th), 0, Math.cos(th)]; }" % W,
        "spin_accumulation": {"enable": True},
        "time_integration": {"max(du)": 0.02, "min(dt)": 5e-20, "max(dt)": 5e-13}}


def read_mesh_nodes(fileName):
    """node coordinates of a gmsh 4.1 text mesh, in file order (the order of the .sol files)"""
    lines = open(fileName).read().split('\n')
    i = lines.index('$Nodes')
    nb_blocks = int(lines[i + 1].split()[0])
    k, xyz = i + 2, []
    for _ in range(nb_blocks):
        n = int(lines[k].split()[3])
        k += 1 + n                   # skip the node tags
        xyz += [[float(v) for v in lines[k + j].split()] for j in range(n)]
        k += n
    return np.array(xyz)


def read_sol(fileName):
    """columns (mx, my, mz) and V of a .sol file"""
    rows = np.array([[float(v) for v in l.split()] for l in open(fileName)
                     if l.strip() and not l.startswith('#')])
    return rows[:, 1:4], rows[:, 8]


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    feellgood = sys.argv[1] if len(sys.argv) > 1 else "../../src/feellgood"
    if not os.path.exists("two_domains.msh"):
        subprocess.run(["gmsh", "-3", "two_domains.geo", "-o", "two_domains.msh"], check=True,
                       stdout=subprocess.DEVNULL)
    xyz = read_mesh_nodes("two_domains.msh") * 1e-9
    inlet = np.abs(xyz[:, 2]) < 1e-12
    x = xyz[inlet, 0]
    V1, TMR = BARRIER["V"], BARRIER["TMR"]

    ok = True
    for name, RA in CASES.items():
        subprocess.run([feellgood, "-"], input=json.dumps(settings(name, RA)), text=True,
                       check=True, stdout=subprocess.DEVNULL)
        m, V = read_sol(name + "_iter0.sol")
        m_err = np.abs(m[inlet, 2] - np.cos(theta(x))).max()
        # current leaving F2 through the barrier, from the potential and the magnetization
        jn = conductance(RA, m[inlet, 2]) * (V[inlet] - V1)
        xr, Vr = reference(RA)
        Vr0 = np.interp(x, xr, Vr[:, 0])
        jnr = np.interp(x, xr, conductance(RA, np.cos(theta(xr))) * (Vr[:, 0] - V1))
        err_V = np.abs(V[inlet] - Vr0).max() / np.abs(Vr0).max()
        err_jn = np.abs(jn - jnr).max() / np.abs(jnr).max()
        # far from the wall (|x| > Lx/2 - 5 w): columns in series with the barrier
        col = lambda cos: -V1 / (1 / conductance(RA, cos) + T / F2["sigma"])
        left, right = x < -LX / 2 + 5 * W, x > LX / 2 - 5 * W
        ratio = jn[left].mean() / jn[right].mean()
        ratio_col = col(1) / col(-1)
        # deviation from the local model of the columns, jn(x) = col(cos theta(x))
        dev_local = np.abs(jnr - col(np.cos(theta(x)))).max() / np.abs(jnr).max()
        print(f"{name:7s} RA = {RA:.0e} Ohm m^2: max error / max at the inlet: V {err_V:.2e}, "
              f"jn {err_jn:.2e} (m {m_err:.1e})")
        print(f"        jn_P/jn_AP far from the wall: feeLLGood {ratio:.4f}, columns "
              f"{ratio_col:.4f} (1 + TMR = {1 + TMR}); reference vs local model of the columns: "
              f"{dev_local:.1e}")
        ok = ok and err_V < TOL and err_jn < TOL
    print("OK" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
