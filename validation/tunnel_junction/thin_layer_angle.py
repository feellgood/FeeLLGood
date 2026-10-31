#!/usr/bin/env python3
"""
Angular validation of the tunnel barrier (key `tunnel` of the surface regions,
src/tunnelJunction.h) on a thin free layer: non-collinear 1D case, where the transverse spin current
produces the spin-transfer torque, checked against a reference written independently of feeLLGood.

Free layer F2 (thin_layer.geo, thickness t, magnetization m = z), tunnel barrier on the inlet
(z = 0), outlet at V = 0 and s = 0 (z = t). The magnetization of F1 is
p = (sin theta, 0, cos theta), theta = 0, 30, ..., 180 deg.

Reference:
- current: the barrier, of conductance per unit area G_T = (1 + c cos theta)/((1 + c) RA),
  c = TMR/(2 + TMR) (cosine law), in series with F2: jn = -V1/(1/G_T + t/sigma) (electrons
  entering F2 for jn > 0);
- spin current entering F2: q = (mu_B/e) jn P_T, P_T = (P1 p + P2 m)/(1 + c cos theta),
  P1 P2 = c; its transverse part, (mu_B/e) jn P1 sin theta/(1 + c cos theta), is along x;
- F2: the equations of feeLLGood (beta = beta' = P, D0 = 2 sigma/(e^2 N0)), 1D, solved exactly by
  the eigenmodes of D0 (I - P^2 m m^T) s'' = [D0 (1 - P^2)/l_sf^2] s + (D0/l_sd^2) s x m, with
  s(t) = 0 and the moment flux q at z = 0: (mu_B/e) P jn m - D0 (I - P^2 m m^T) s'(0) = q.

The torque on F2 is proportional to m x s: s_y gives the damping-like torque (along
m x (m x p), that is -x) and s_x the field-like torque (along m x p, that is y). Their integrals
over the thickness are proportional to the transverse spin current entering F2, so that the
damping-like torque follows the form of Slonczewski, eps(theta) sin(theta) with
eps = C/((Lambda^2 + 1) + (Lambda^2 - 1) cos(theta)): per unit current (int s_y dz/jn),
Lambda^2 = (1 + c)/(1 - c) = 1 + TMR exactly; at fixed voltage (int s_y dz), Lambda^2 is close to
1, the barrier dominating the resistance.

Usage: thin_layer_angle.py [path to feellgood], from any directory; the mesh is built with gmsh if
missing. Exit status 0 if the errors are below the tolerance at every angle.
"""

import json
import os
import subprocess
import sys

import numpy as np

E = 1.602176634e-19                  # C
MU_B = 9.2740100783e-24              # J/T

TOL = 2e-2                           # relative tolerance

T2 = 3e-9                            # thickness of F2 (m), as in thin_layer.geo
F2 = {"sigma": 4e6, "N0": 1.344e47, "P": 0.5, "lsf": 5e-9, "lsd": 1e-9}
BARRIER = {"RA": 1e-12, "TMR": 1.5, "P1": 0.7, "V": -0.1}   # V < 0: electrons from F1 into F2
ANGLES = [0, 30, 60, 90, 120, 150, 180]


def current(theta):
    """normal current density leaving F2 through the barrier (A/m^2)"""
    c = BARRIER["TMR"] / (2 + BARRIER["TMR"])
    G = (1 + c * np.cos(np.radians(theta))) / ((1 + c) * BARRIER["RA"])
    return -BARRIER["V"] / (1 / G + T2 / F2["sigma"])


def reference(theta, z):
    """spin accumulation s (A/m), rows at the positions z (m) of F2, for the angle theta (deg)"""
    p = np.array([np.sin(np.radians(theta)), 0, np.cos(np.radians(theta))])
    m = np.array([0.0, 0.0, 1.0])
    c = BARRIER["TMR"] / (2 + BARRIER["TMR"])
    P1 = BARRIER["P1"]
    P2 = c / P1
    jn = current(theta)
    q = MU_B / E * jn * (P1 * p + P2 * m) / (1 + c * p.dot(m))
    D0 = 2 * F2["sigma"] / (E**2 * F2["N0"])
    P = F2["P"]
    I3 = np.eye(3)
    C = np.array([[0, m[2], -m[1]], [-m[2], 0, m[0]], [m[1], -m[0], 0]])  # C s = s x m
    A = D0 * (I3 - P**2 * np.outer(m, m))
    M = D0 * (1 - P**2) / F2["lsf"]**2 * I3 + D0 / F2["lsd"]**2 * C
    lam, V = np.linalg.eig(np.linalg.solve(A, M))
    kappa = np.sqrt(lam.astype(complex))
    # s(z) = V diag(sinh(kappa (T2 - z))) cf, s(T2) = 0; A s'(0) = (mu_B/e) jn P m - q
    dS0 = V @ np.diag(-kappa * np.cosh(kappa * T2))
    cf = np.linalg.solve(A @ dS0, MU_B / E * jn * P * m - q)
    s = np.array([(V @ (np.sinh(kappa * (T2 - zz)) * cf)) for zz in z])
    return s.real


def settings(theta, name):
    p = [float(np.sin(np.radians(theta))), 0.0, float(np.cos(np.radians(theta)))]
    return {
        "outputs": {"file_basename": name, "final_time": 1e-13, "evol_time_step": 1e-12,
                    "mag_config_every": 10},
        "mesh": {
            "filename": "thin_layer.msh",
            "length_unit": 1e-9,
            "volume_regions": {"F2": {
                "Js": 1.0, "Ae": 1.3e-11, "alpha_LLG": 0.05, "K": 0,
                "sigma": F2["sigma"], "dens_state": F2["N0"], "P": F2["P"],
                "l_sd": F2["lsd"], "l_sf": F2["lsf"]}},
            "surface_regions": {
                "outlet": {"V": 0.0, "s": [0, 0, 0]},
                "inlet": {"tunnel": {"p": p, "RA": BARRIER["RA"], "TMR": BARRIER["TMR"],
                                     "P1": BARRIER["P1"], "V": BARRIER["V"]}}}},
        "initial_magnetization": [0, 0, 1],
        "spin_accumulation": {"enable": True},
        "time_integration": {"max(du)": 0.02, "min(dt)": 5e-20, "max(dt)": 5e-13}}


def read_mesh(fileName):
    """node coordinates (file order) and tetrahedra (0-based indices in that order) of a gmsh 4.1
    text mesh"""
    lines = open(fileName).read().split('\n')
    i = lines.index('$Nodes')
    nb_blocks = int(lines[i + 1].split()[0])
    k, tags, xyz = i + 2, [], []
    for _ in range(nb_blocks):
        n = int(lines[k].split()[3])
        tags += [int(lines[k + 1 + j]) for j in range(n)]
        k += 1 + n
        xyz += [[float(v) for v in lines[k + j].split()] for j in range(n)]
        k += n
    idx = {tag: i for i, tag in enumerate(tags)}
    i = lines.index('$Elements')
    nb_blocks = int(lines[i + 1].split()[0])
    k, tets = i + 2, []
    for _ in range(nb_blocks):
        dim, _, etype, n = (int(v) for v in lines[k].split())
        if etype == 4:
            tets += [[idx[int(v)] for v in lines[k + 1 + j].split()[1:]] for j in range(n)]
        k += 1 + n
    return np.array(xyz) * 1e-9, np.array(tets)


def read_s(fileName):
    """spin accumulation columns (sx, sy, sz) of a .sol file"""
    rows = [l.split() for l in open(fileName) if l.strip() and not l.startswith('#')]
    return np.array([[float(v) for v in r[5:8]] for r in rows])


def thickness_integral(xyz, tets, s):
    """integral over the thickness of s (P1, exact), the volume integral divided by the section"""
    vol = np.abs(np.linalg.det(xyz[tets[:, 1:]] - xyz[tets[:, :1]])) / 6
    total = (vol[:, None] * s[tets].mean(axis=1)).sum(axis=0)
    return total / (vol.sum() / T2)


def fit_lambda2(thetas, I):
    """Lambda^2 of the form of Slonczewski fitted to the damping-like torque I(theta) (0 < theta <
    180 deg): sin(theta)/I = alpha + beta cos(theta), and the relative deviation of the fit"""
    th = np.radians(thetas)
    y = np.sin(th) / I
    beta, alpha = np.polyfit(np.cos(th), y, 1)
    resid = np.abs(np.polyval([beta, alpha], np.cos(th)) - y).max() / np.abs(y).max()
    return (alpha + beta) / (alpha - beta), resid


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    feellgood = sys.argv[1] if len(sys.argv) > 1 else "../../src/feellgood"
    if not os.path.exists("thin_layer.msh"):
        subprocess.run(["gmsh", "-3", "thin_layer.geo", "-o", "thin_layer.msh"], check=True,
                       stdout=subprocess.DEVNULL)
    xyz, tets = read_mesh("thin_layer.msh")

    ok = True
    rows = []
    print(" theta | jn [A/m^2] | damping-like, int s_y dz [A]: feeLLGood  ref | "
          "field-like, int s_x dz [A]: feeLLGood  ref | max err s")
    for theta in ANGLES:
        name = f"thin_layer_theta{theta}"
        subprocess.run([feellgood, "-"], input=json.dumps(settings(theta, name)), text=True,
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        s = read_s(name + "_iter0.sol")
        s_ref = reference(theta, xyz[:, 2])
        err = np.abs(s - s_ref).max() / np.abs(s_ref).max()
        I, I_ref = thickness_integral(xyz, tets, s), thickness_integral(xyz, tets, s_ref)
        print(f" {theta:5d} | {current(theta):.4e} | {I[1]:29.4e} {I_ref[1]:11.4e} | "
              f"{I[0]:27.4e} {I_ref[0]:11.4e} | {err:.2e}")
        ok = ok and err < TOL
        if 0 < theta < 180:
            rows.append((theta, I[1], current(theta)))

    th = np.array([r[0] for r in rows])
    Iy = np.array([r[1] for r in rows])
    jn = np.array([r[2] for r in rows])
    L2_V, dev_V = fit_lambda2(th, Iy)
    L2_I, dev_I = fit_lambda2(th, Iy / jn)
    print(f"damping-like torque (feeLLGood), form of Slonczewski: at fixed voltage "
          f"Lambda^2 = {L2_V:.4f} (deviation {dev_V:.1e}); per unit current Lambda^2 = {L2_I:.4f} "
          f"(deviation {dev_I:.1e}), 1 + TMR = {1 + BARRIER['TMR']}")
    print("OK" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
