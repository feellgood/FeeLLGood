#!/usr/bin/env python3
"""
1D collinear validation of the spin injection through a polarizer stack (key `slonczewski` of the
surface regions, src/slonczewskiInterface.h), against the analytic solution of Valet and Fert.

Stack: F1 (bulk, thickness t1, reservoir on its far side) / interface F1/N / spacer N of zero
thickness / interface N/F2 / F2 (the meshed wire, length L, reservoir s = 0 at its end). The
magnetizations are along z, F1 parallel (P) or antiparallel (AP) to F2: the problem is 1D and
collinear, the mixing conductances play no role.

Valet-Fert, per channel of moment +-z (charge current densities j+-, conductivities
sigma+- = sigma (1 +- beta)/2, interface resistances r+- = 2 AR (1 -+ gamma)), with the spin
accumulation phi = (mu+ - mu-)/(2e) in volts:
    bulk:      j+ - j- = beta j + sigma (1 - beta^2) dphi/dx,    phi'' = phi/lsf^2
    interface: j+ - j- = gamma j + (phi_right - phi_left)/AR
The moment flux is -(mu_B/e)(j+ - j-), and the spin accumulation of feeLLGood is
s = mu_B N0 e phi/2 (in A/m), so that the bulk equations are the ones of feeLLGood with
beta = beta' = P, D0 = 2 sigma/(e^2 N0), and the interface conductances, in the units of the spin
diffusion problem, are
    Gpar = 2/(e^2 N0 AR)    [m/s]
which is the conversion checked here (done by feeLLGood, src/slonczewskiInterface.h: the settings
give AR and gamma). N0 is the density of states of F2: the stack being not meshed, all its
conductances (interfaces and bulk of F1) are converted with it, the density of states of F1 does
not enter.

Usage: valet_fert.py [path to feellgood], from any directory; the mesh is built with gmsh if
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

# F2: the meshed wire (lengths in m)
L = 300e-9
F2 = {"sigma": 4e6, "N0": 1.344e47, "beta": 0.3, "lsf": 10e-9, "lsd": 1e-9}
# polarizer stack
F1 = {"sigma": 2e6, "beta": 0.5, "lsf": 20e-9, "t": 15e-9}
ITF1 = {"AR": 0.5e-15, "gamma": 0.75}   # F1/N, AR in Ohm m^2
ITF2 = {"AR": 0.3e-15, "gamma": 0.6}    # N/F2
JN = 1e12                            # A/m^2, normal current density on the inlet (outward
                                     # normal -z): electrons flowing from F1 to F2


def valet_fert(sign_F1, z):
    """spin accumulation s_z (A/m) at the positions z (m) in F2, F1 along sign_F1 z"""
    j = -JN                          # charge current density along +z
    b1, g1 = sign_F1 * F1["beta"], sign_F1 * ITF1["gamma"]
    b2, g2 = F2["beta"], ITF2["gamma"]
    l1, l2, t1 = F1["lsf"], F2["lsf"], F1["t"]
    c1 = F1["sigma"] * (1 - b1**2)
    c2 = F2["sigma"] * (1 - b2**2)
    # unknowns X = (a1, d1, phiN, a2, d2):
    # phi1(x) = a1 exp(x/l1) + d1 exp(-(x + t1)/l1) on [-t1, 0]
    # phi2(x) = a2 exp(-x/l2) + d2 exp((x - L)/l2) on [0, L]
    e1, e2 = np.exp(-t1 / l1), np.exp(-L / l2)
    phi1_0 = np.array([1, e1, 0, 0, 0])                    # phi1(0)
    dphi1_0 = np.array([1 / l1, -e1 / l1, 0, 0, 0])        # phi1'(0)
    phi2_0 = np.array([0, 0, 0, 1, e2])                    # phi2(0)
    dphi2_0 = np.array([0, 0, 0, -1 / l2, e2 / l2])        # phi2'(0)
    phiN = np.array([0, 0, 1, 0, 0])
    q1 = (phiN - phi1_0) / ITF1["AR"]                      # interface F1/N, j+ - j- minus gamma j
    q2 = (phi2_0 - phiN) / ITF2["AR"]                      # interface N/F2
    A = np.array([
        [e1, 1, 0, 0, 0],                                  # phi1(-t1) = 0
        [0, 0, 0, e2, 1],                                  # phi2(L) = 0
        c1 * dphi1_0 - q1,                                 # bulk F1 = interface F1/N
        q1 - q2,                                           # through the spacer
        c2 * dphi2_0 - q2])                                # interface N/F2 = bulk F2
    rhs = np.array([0, 0, (g1 - b1) * j, (g2 - g1) * j, (g2 - b2) * j])
    a1, d1, pN, a2, d2 = np.linalg.solve(A, rhs)
    phi = a2 * np.exp(-z / l2) + d2 * np.exp((z - L) / l2)
    return MU_B * F2["N0"] * E * phi / 2


def settings(sign_F1, name):
    stack = {
        "p": [0, 0, sign_F1],
        "F1": {"sigma": F1["sigma"], "P": F1["beta"], "l_sf": F1["lsf"], "thickness": F1["t"]},
        "F1N": {"AR": ITF1["AR"], "gamma": ITF1["gamma"], "Gr": 5e14},
        "NF2": {"AR": ITF2["AR"], "gamma": ITF2["gamma"], "Gr": 5e14}}
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
    if not os.path.exists("wire.msh"):
        subprocess.run(["gmsh", "-3", "wire.geo", "-o", "wire.msh"], check=True,
                       stdout=subprocess.DEVNULL)
    z = read_mesh_nodes("wire.msh")[:, 2] * 1e-9

    ok = True
    for sign_F1, name in [(1, "P"), (-1, "AP")]:
        subprocess.run([feellgood, "-"], input=json.dumps(settings(sign_F1, name)), text=True,
                       check=True, stdout=subprocess.DEVNULL)
        s = read_s(name + "_iter0.sol")
        s_ref = valet_fert(sign_F1, z)
        err = np.abs(s[:, 2] - s_ref).max() / np.abs(s_ref).max()
        transverse = np.abs(s[:, :2]).max() / np.abs(s_ref).max()
        inlet = np.abs(z) < 1e-12
        print(f"{name:2s}: s_z(inlet) feeLLGood = {s[inlet, 2].mean():9.3f} A/m, "
              f"Valet-Fert = {s_ref[inlet].mean():9.3f} A/m; "
              f"max error / max|s| = {err:.2e}, transverse = {transverse:.1e}")
        ok = ok and err < TOL and transverse < TOL
    print("OK" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
