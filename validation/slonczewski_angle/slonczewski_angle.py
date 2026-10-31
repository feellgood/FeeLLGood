#!/usr/bin/env python3
"""
Angular validation of the spin injection through a polarizer stack (key `slonczewski` of the
surface regions, src/slonczewskiInterface.h) on a thin free layer F2: non-collinear 1D case, checked
against a reference written independently of feeLLGood.

Stack: F1 (ideal reservoir) / interface F1/N / spacer N of zero thickness / interface N/F2 / F2 (the
meshed thin layer, thickness t2, magnetization m = z) / reservoir s = 0. The magnetization of F1 is
p = (sin theta, 0, cos theta).

Reference:
- interfaces, circuit theory of Brataas et al. (Phys. Rep. 427, 157, 2006), per unit area, in the
  variables of Valet-Fert (charge current density j along +z, spin accumulation phi = mu_s/e in
  volts, mu_s = (mu+ - mu-)/2, + being the majority spin, along the magnetization); the spin
  current from L to R through an interface of magnetization u is
      j_s = gamma j u + (1/AR) u (u.(phi_R - phi_L)) + 2 Re G_up,down (phi_R - phi_L)_perp
            + 2 Im G_up,down (phi_R - phi_L) x u,
  the moment flux being -(mu_B/e) j_s. This is the equation of Brataas et al. for the transverse
  spin current from N into F, 2 Re G_up,down m x (V_S x m) + 2 Im G_up,down V_S x m, with
  G_up,down = (e^2/h) [M - sum r_up r_down^*] (up: majority spin, along m) and V_S = -phi
  (mu = -e V); the sign of the term in Im G_up,down also follows from the scattering of a single
  channel: the transverse polarization absorbed is 1 - r_up^* r_down = conj(G_up,down);
  through N/F2, two models of the transverse spin current (key transverse of NF2): sink (default),
  it depends on the spin accumulation of the spacer only (circuit theory, phi_F replaced by
  m (m.phi_F) in N/F2), series, it depends on phi_F - phi_N;
- F2, equations of feeLLGood (beta = beta' = P, D0 = 2 sigma/(e^2 N0)), 1D, solved exactly by the
  eigenmodes of D0 (I - P^2 m m^T) s'' = [D0 (1 - P^2)/l_sf^2] s + (D0/l_sd^2) s x m;
- s = mu_B N0 e phi/2 (validation/valet_fert_1d).

The torque on F2 is proportional to m x s: with m = z, s_y gives the damping-like torque (along
m x (m x p), that is x) and s_x the field-like torque (along m x p, that is y). Their integrals
over the thickness are compared, and the damping-like one is compared to the form of Slonczewski,
eps(theta) sin(theta) with eps = C/((Lambda^2 + 1) + (Lambda^2 - 1) cos(theta)).

Usage: slonczewski_angle.py [path to feellgood], from any directory; the mesh is built with gmsh
if missing. Exit status 0 if the errors are below the tolerance at every angle.
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
# F1/N and N/F2, Gr + i Gi = G_up,down (S/m^2); Im G_up,down is large here, so that its sign matters
ITF1 = {"AR": 0.5e-15, "gamma": 0.75, "Gr": 0.55e15, "Gi": 0.15e15}
ITF2 = {"AR": 0.5e-15, "gamma": 0.75, "Gr": 0.55e15, "Gi": 0.25e15}
IM_SIGN = 1                          # sign of the term in Im G_up,down of the reference (-1: checks
                                     # that the test fails with the opposite sign)
JN = 1e12                            # A/m^2 on the inlet (outward normal -z): electrons from F1 to F2
ANGLES = [0, 30, 60, 90, 120, 150, 180]


def reference(theta, z, model="sink"):
    """spin accumulation s (A/m), rows at the positions z (m) of F2, for the angle theta (deg) and
    the model of the transverse spin current through N/F2 (sink or series)"""
    p = np.array([np.sin(np.radians(theta)), 0, np.cos(np.radians(theta))])
    m = np.array([0.0, 0.0, 1.0])
    j = -JN                                        # charge current density along +z
    D0 = 2 * F2["sigma"] / (E**2 * F2["N0"])
    P = F2["P"]
    I3 = np.eye(3)
    C = np.array([[0, m[2], -m[1]], [-m[2], 0, m[0]], [m[1], -m[0], 0]])  # C s = s x m
    A = D0 * (I3 - P**2 * np.outer(m, m))
    M = D0 * (1 - P**2) / F2["lsf"]**2 * I3 + D0 / F2["lsd"]**2 * C
    lam, V = np.linalg.eig(np.linalg.solve(A, M))
    kappa = np.sqrt(lam.astype(complex))
    # s(z) = V diag(sinh(kappa (T2 - z))) c, s(T2) = 0
    S0 = V @ np.diag(np.sinh(kappa * T2))                  # s(0) = S0 c
    dS0 = V @ np.diag(-kappa * np.cosh(kappa * T2))        # s'(0) = dS0 c
    to_phi = 2 / (MU_B * F2["N0"] * E)                     # phi = to_phi s

    def contact(u, itf):
        """j_s = gamma j u + K (phi_R - phi_L)"""
        Cu = np.array([[0, u[2], -u[1]], [-u[2], 0, u[0]], [u[1], -u[0], 0]])  # Cu w = w x u
        K = (np.outer(u, u) / itf["AR"] + 2 * itf["Gr"] * (I3 - np.outer(u, u))
             + IM_SIGN * 2 * itf["Gi"] * Cu)
        return itf["gamma"] * j * u, K

    b1, K1 = contact(p, ITF1)                    # L = F1 (phi = 0), R = N
    b2, K2 = contact(m, ITF2)                    # L = N, R = F2
    # unknowns X = (c, phiN); j_s1 = b1 + K1 phiN, j_s2 = b2 + K2 (Pm to_phi S0 c - phiN),
    # Pm = m m^T for the sink (only the longitudinal part of phi_F acts on N/F2)
    Pm = np.outer(m, m) if model == "sink" else I3
    Z = np.zeros((3, 3))
    js1 = np.hstack([Z, K1])
    js2 = np.hstack([K2 @ Pm @ (to_phi * S0), -K2])
    # moment flux in F2 at z = 0 (A includes D0): -A s'(0) - (mu_B/e) P j m = -(mu_B/e) j_s2
    Qx = np.hstack([-A @ dS0, Z])
    lhs = np.vstack([js1 - js2, Qx + (MU_B / E) * js2])
    rhs = np.concatenate([b2 - b1, (MU_B / E) * P * j * m - (MU_B / E) * b2])
    c = np.linalg.solve(lhs, rhs)[:3]
    s = np.array([(V @ (np.sinh(kappa * (T2 - zz)) * c)) for zz in z])
    return s.real


def settings(theta, name, model):
    p = [float(np.sin(np.radians(theta))), 0.0, float(np.cos(np.radians(theta)))]
    itf = lambda d: {"AR": d["AR"], "gamma": d["gamma"], "Gr": d["Gr"], "Gi": d["Gi"]}
    NF2 = dict(itf(ITF2), transverse=model)
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
                "inlet": {"jn": JN, "slonczewski": {"p": p, "F1N": itf(ITF1), "NF2": NF2}}}},
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


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    feellgood = sys.argv[1] if len(sys.argv) > 1 else "../../src/feellgood"
    if not os.path.exists("thin_layer.msh"):
        subprocess.run(["gmsh", "-3", "thin_layer.geo", "-o", "thin_layer.msh"], check=True,
                       stdout=subprocess.DEVNULL)
    xyz, tets = read_mesh("thin_layer.msh")

    ok = True
    for model in ["sink", "series"]:
        rows = []
        print(f"transverse spin current through N/F2: {model}")
        print(" theta | damping-like, int s_y dz [A]: feeLLGood  ref | "
              "field-like, int s_x dz [A]: feeLLGood  ref | max err s")
        for theta in ANGLES:
            name = f"{model}_theta{theta}"
            subprocess.run([feellgood, "-"], input=json.dumps(settings(theta, name, model)),
                           text=True, check=True, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
            s = read_s(name + "_iter0.sol")
            s_ref = reference(theta, xyz[:, 2], model)
            err = np.abs(s - s_ref).max() / np.abs(s_ref).max()
            I, I_ref = thickness_integral(xyz, tets, s), thickness_integral(xyz, tets, s_ref)
            print(f" {theta:5d} | {I[1]:29.4e} {I_ref[1]:11.4e} | {I[0]:27.4e} {I_ref[0]:11.4e} | "
                  f"{err:.2e}")
            ok = ok and err < TOL
            rows.append((theta, I_ref[1]))

        # damping-like torque: sin(theta)/int s_y = alpha + beta cos(theta), Slonczewski's form
        th = np.radians([r[0] for r in rows if 0 < r[0] < 180])
        y = np.sin(th) / np.array([r[1] for r in rows if 0 < r[0] < 180])
        beta, alpha = np.polyfit(np.cos(th), y, 1)
        resid = np.abs(np.polyval([beta, alpha], np.cos(th)) - y).max() / np.abs(y).max()
        # alpha = (Lambda^2 + 1)/C, beta = (Lambda^2 - 1)/C
        print(f"damping-like torque (reference): Lambda^2 = {(alpha + beta)/(alpha - beta):.4f}, "
              f"C = {2/(alpha - beta):.4e} A, deviation from the form of Slonczewski {resid:.1e}")
    print("OK" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
