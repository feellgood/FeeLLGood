#!/usr/bin/env python3
"""
1D collinear validation of the continuity of the spin accumulation s between two meshed volume
regions of different densities of states N0 (dens_state).

Wire: non-magnetic spacer N (0 < z < tN, meshed) / free layer F2 (tN < z < L) / reservoir s = 0.
The spin is injected on z = 0 through a polarizer stack without N/F2 interface (key `slonczewski`,
the spacer being meshed): F1 (ideal reservoir) / interface F1/N. The magnetizations are along z,
F1 parallel (P) or antiparallel (AP) to F2.

feeLLGood takes s continuous between the regions, as is the spin current. With s = mu_B N0 mu_s/2
(validation/valet_fert_1d), this imposes N0_N mu_s,N = N0_F mu_s,F on the interface N/F2: mu_s is
continuous only if the densities of states are equal. The reference solves this model in the
variables of Valet-Fert (phi = mu_s/e in volts, charge current density j along +z):
    bulk:      j_s = beta j + sigma (1 - beta^2) dphi/dz,    phi'' = phi/lsf^2
    inlet:     j_s = gamma j + phi(0)/AR                      (F1 reservoir, phi = 0)
    z = tN:    N0_N phi_N = N0_F phi_F (s continuous), j_s continuous
    z = L:     phi = 0
s = mu_B N0 e phi/2 in each region. The script also gives the solution with mu_s continuous
(phi_N = phi_F), to show the difference between the two models.

Usage: two_regions.py [path to feellgood], from any directory; the mesh is built with gmsh if
missing. Exit status 0 if the relative error is below the tolerance in both configurations.
"""

import json
import os
import subprocess
import sys

import numpy as np

E = 1.602176634e-19                  # C
MU_B = 9.2740100783e-24              # J/T

TOL = 1e-2                           # relative tolerance on the spin accumulation

L, TN = 300e-9, 40e-9                # as in two_regions.geo
N = {"sigma": 6e7, "N0": 0.6e47, "lsf": 50e-9}                     # spacer (Cu-like)
F2 = {"sigma": 4e6, "N0": 1.344e47, "beta": 0.5, "lsf": 10e-9, "lsd": 1e-9}
ITF1 = {"AR": 0.5e-15, "gamma": 0.75, "Gr": 0.55e15}              # F1/N
JN = 1e12                            # A/m^2 on the inlet (outward normal -z): electrons from F1 to F2


def reference(sign_F1, z, s_continuous=True):
    """spin accumulation s_z (A/m) at the positions z (m), F1 along sign_F1 z"""
    j = -JN                          # charge current density along +z
    g1 = sign_F1 * ITF1["gamma"]
    bF = F2["beta"]
    lN, lF = N["lsf"], F2["lsf"]
    cN, cF = N["sigma"], F2["sigma"] * (1 - bF**2)
    # unknowns X = (a, b, c, d):
    # phi_N(z) = a exp(z/lN) + b exp(-z/lN) on [0, tN]
    # phi_F(z) = c exp(-(z - tN)/lF) + d exp((z - L)/lF) on [tN, L]
    eN, eF = np.exp(TN / lN), np.exp(-(L - TN) / lF)
    phiN_0, dphiN_0 = np.array([1, 1, 0, 0]), np.array([1 / lN, -1 / lN, 0, 0])
    phiN_t, dphiN_t = np.array([eN, 1 / eN, 0, 0]), np.array([eN / lN, -1 / (eN * lN), 0, 0])
    phiF_t, dphiF_t = np.array([0, 0, 1, eF]), np.array([0, 0, -1 / lF, eF / lF])
    phiF_L = np.array([0, 0, np.exp(-(L - TN) / lF), 1])
    rN, rF = (N["N0"], F2["N0"]) if s_continuous else (1.0, 1.0)
    A = np.array([
        cN * dphiN_0 - phiN_0 / ITF1["AR"],          # bulk N = interface F1/N at z = 0
        rN * phiN_t - rF * phiF_t,                   # rule on the interface N/F2
        cN * dphiN_t - cF * dphiF_t,                 # continuity of the spin current
        phiF_L])                                     # reservoir at z = L
    rhs = np.array([g1 * j, 0, bF * j, 0])
    a, b, c, d = np.linalg.solve(A, rhs)
    phi = np.where(z < TN, a * np.exp(z / lN) + b * np.exp(-z / lN),
                   c * np.exp(-(z - TN) / lF) + d * np.exp((z - L) / lF))
    N0 = np.where(z < TN, N["N0"], F2["N0"])
    return MU_B * N0 * E * phi / 2


def settings(sign_F1, name):
    stack = {"p": [0, 0, sign_F1],
             "F1N": {"AR": ITF1["AR"], "gamma": ITF1["gamma"], "Gr": ITF1["Gr"]}}
    return {
        "outputs": {"file_basename": name, "final_time": 1e-13, "evol_time_step": 1e-12,
                    "mag_config_every": 10},
        "mesh": {
            "filename": "two_regions.msh",
            "length_unit": 1e-9,
            "volume_regions": {
                "N": {"Js": 0, "sigma": N["sigma"], "dens_state": N["N0"], "P": 0,
                      "l_sd": 1e-9, "l_sf": N["lsf"]},
                "F2": {"Js": 1.0, "Ae": 1.3e-11, "alpha_LLG": 0.05, "K": 0,
                       "sigma": F2["sigma"], "dens_state": F2["N0"], "P": F2["beta"],
                       "l_sd": F2["lsd"], "l_sf": F2["lsf"]}},
            "surface_regions": {
                "outlet": {"V": 0.0, "s": [0, 0, 0]},
                "inlet": {"jn": JN, "slonczewski": stack}}},
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


def read_s(fileName):
    """spin accumulation columns (sx, sy, sz) of a .sol file"""
    rows = [l.split() for l in open(fileName) if l.strip() and not l.startswith('#')]
    return np.array([[float(v) for v in r[5:8]] for r in rows])


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    feellgood = sys.argv[1] if len(sys.argv) > 1 else "../../src/feellgood"
    if not os.path.exists("two_regions.msh"):
        subprocess.run(["gmsh", "-3", "two_regions.geo", "-o", "two_regions.msh"], check=True,
                       stdout=subprocess.DEVNULL)
    z = read_mesh_nodes("two_regions.msh")[:, 2] * 1e-9
    z = np.where(np.abs(z - TN) < 1e-12, TN, z)      # nodes of the interface: s continuous there

    ok = True
    print(f"dens_state: N {N['N0']:.3e}, F2 {F2['N0']:.3e} J^-1 m^-3")
    for sign_F1, name in [(1, "P"), (-1, "AP")]:
        subprocess.run([feellgood, "-"], input=json.dumps(settings(sign_F1, name)), text=True,
                       check=True, stdout=subprocess.DEVNULL)
        s = read_s(name + "_iter0.sol")
        s_ref = reference(sign_F1, z)
        s_mu = reference(sign_F1, z, s_continuous=False)
        err = np.abs(s[:, 2] - s_ref).max() / np.abs(s_ref).max()
        dev = np.abs(s_mu - s_ref).max() / np.abs(s_ref).max()
        itf = np.abs(z - TN) < 1e-12
        print(f"{name:2s}: s_z(N/F2) feeLLGood = {s[itf, 2].mean():9.3f} A/m, "
              f"model s continuous = {s_ref[itf].mean():9.3f} A/m; max error / max|s| = "
              f"{err:.2e}; model mu_s continuous: max deviation / max|s| = {dev:.2e}")
        ok = ok and err < TOL
    print("OK" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
