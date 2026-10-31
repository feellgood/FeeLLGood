#!/usr/bin/env python3
"""
1D validation of the tunnel barrier (key `tunnel` of the surface regions, src/tunnelJunction.h) on a
wire with a uniform magnetization: potential and longitudinal spin accumulation, against the
analytic solution.

Wire F2 (wire.geo, length L along z), tunnel barrier on the inlet (z = 0), outlet at V = 0 and
s = 0 (z = L). The magnetization of F2 is m = z, the one of F1 is p = (sin theta, 0, cos theta),
theta = 0, 90, 180 deg.

- Potential: the barrier, of conductance per unit area G_T = (1 + c cos theta)/((1 + c) RA),
  c = TMR/(2 + TMR) (cosine law), is in series with the wire: the current density leaving the wire
  through the barrier is jn = -V1/(1/G_T + L/sigma) and the potential at the inlet is
  V(0) = -jn L/sigma.
- Spin: the spin current entering the wire is (mu_B/e) jn P_T,
  P_T = (P1 p + P2 m)/(1 + c cos theta), P1 P2 = c. Along m, it is the problem of
  validation/valet_fert_1d with a given injected polarization P_l = P_T.m: in the variables of Valet-Fert (spin accumulation phi in volts,
  s = mu_B N0 e phi/2), j_s = beta j + sigma (1 - beta^2) phi', phi'' = phi/l_sf^2, phi(L) = 0 and
  j_s(0) = P_l j, j = -jn, so that
      phi(z) = (P_l - beta) jn l_sf sinh((L - z)/l_sf) / (sigma (1 - beta^2) cosh(L/l_sf)).
  The transverse part (theta = 90 deg) is checked by thin_layer_angle.py.

Usage: wire_1d.py [path to feellgood], from any directory; the mesh is built with gmsh if missing.
Exit status 0 if the errors are below the tolerance at every angle.
"""

import json
import os
import subprocess
import sys

import numpy as np

E = 1.602176634e-19                  # C
MU_B = 9.2740100783e-24              # J/T

TOL = 1e-2                           # relative tolerance

L = 100e-9                           # length of the wire (m), as in wire.geo
F2 = {"sigma": 4e6, "N0": 1.344e47, "beta": 0.3, "lsf": 10e-9, "lsd": 1e-9}
BARRIER = {"RA": 1e-12, "TMR": 1.5, "P1": 0.7, "V": -0.1}   # V < 0: electrons from F1 into F2
ANGLES = [0, 90, 180]


def reference(theta, z):
    """potential at the inlet (V), normal current density (A/m^2) and spin accumulation s_z (A/m)
    at the positions z (m), for the angle theta (deg) between p and m"""
    cos = np.cos(np.radians(theta))
    c = BARRIER["TMR"] / (2 + BARRIER["TMR"])
    P1 = BARRIER["P1"]
    P2 = c / P1
    G = (1 + c * cos) / ((1 + c) * BARRIER["RA"])
    jn = -BARRIER["V"] / (1 / G + L / F2["sigma"])
    V0 = -jn * L / F2["sigma"]
    Pl = (P1 * cos + P2) / (1 + c * cos)
    b, l = F2["beta"], F2["lsf"]
    phi = (Pl - b) * jn * l * np.sinh((L - z) / l) / (F2["sigma"] * (1 - b**2) * np.cosh(L / l))
    return V0, jn, MU_B * F2["N0"] * E * phi / 2


def settings(theta, name):
    p = [float(np.sin(np.radians(theta))), 0.0, float(np.cos(np.radians(theta)))]
    return {
        "outputs": {"file_basename": name, "final_time": 1e-13, "evol_time_step": 1e-12,
                    "mag_config_every": 10},
        "mesh": {
            "filename": "wire.msh",
            "length_unit": 1e-9,
            "volume_regions": {"wire": {
                "Js": 1.0, "Ae": 1.3e-11, "alpha_LLG": 0.05, "K": 0,
                "sigma": F2["sigma"], "dens_state": F2["N0"], "P": F2["beta"],
                "l_sd": F2["lsd"], "l_sf": F2["lsf"]}},
            "surface_regions": {
                "outlet": {"V": 0.0, "s": [0, 0, 0]},
                "inlet": {"tunnel": {"p": p, "RA": BARRIER["RA"], "TMR": BARRIER["TMR"],
                                     "P1": BARRIER["P1"], "V": BARRIER["V"]}}}},
        "initial_magnetization": [0, 0, 1],
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
    """columns (sx, sy, sz) and V of a .sol file"""
    rows = np.array([[float(v) for v in l.split()] for l in open(fileName)
                     if l.strip() and not l.startswith('#')])
    return rows[:, 5:8], rows[:, 8]


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    feellgood = sys.argv[1] if len(sys.argv) > 1 else "../../src/feellgood"
    if not os.path.exists("wire.msh"):
        subprocess.run(["gmsh", "-3", "wire.geo", "-o", "wire.msh"], check=True,
                       stdout=subprocess.DEVNULL)
    z = read_mesh_nodes("wire.msh")[:, 2] * 1e-9
    inlet = np.abs(z) < 1e-12

    ok = True
    for theta in ANGLES:
        name = f"wire_theta{theta}"
        subprocess.run([feellgood, "-"], input=json.dumps(settings(theta, name)), text=True,
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        s, V = read_sol(name + "_iter0.sol")
        V0, jn, sz = reference(theta, z)
        err_V = np.abs(V[inlet] - V0).max() / abs(V0)
        err_s = np.abs(s[:, 2] - sz).max() / np.abs(sz).max()
        print(f"theta {theta:3d}: jn = {jn:.4e} A/m^2, V(inlet) feeLLGood {V[inlet].mean():.6e} V, "
              f"ref {V0:.6e} V, error {err_V:.1e}; s_z(inlet) feeLLGood {s[inlet, 2].mean():8.4f}"
              f" A/m, ref {sz[inlet].mean():8.4f} A/m, max error / max|s_z| {err_s:.2e}")
        ok = ok and err_V < TOL and err_s < TOL
    print("OK" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
