#!/usr/bin/env python3
"""
Macrospin estimate of the spin-transfer switching of the example
pillar_slonczewski_current_density_ramp1ns, from the parameters of settings.yml, and comparison
with the simulation (pillar.evol).

The polarizer stack is converted as feeLLGood does (src/slonczewskiInterface.h): k = 2/(e^2 N0),
Gpar = k/AR, Gr + i Gi = 2 k G_up,down, bulk of F1 folded into F1/N (renormalizeF1), spacer
eliminated (effectiveBC). For a thin free layer whose spin accumulation is neglected (s_F = 0),
the spin current entering F2 is q = (mu_B/e) jn Peff(theta), theta being the angle between m and
p; its transverse part is absorbed by F2. The damping-like torque per unit area is
(hbar/2e) jn eta(theta) sin(theta), with eta(theta) sin(theta) the component of Peff along
m x (m x p) (up to the sign, Peff . (p - (p.m) m)/|...|).

Critical current density of the instability of the antiparallel state (m = -p, AP -> P), for a
perpendicular free layer of thickness t (macrospin, linearized LLG):
    jc = (2e/hbar) alpha mu0 Ms H_K,eff t / eta(pi),
    H_K,eff = 2K/(mu0 Ms) - (Nz - Nx) Ms,
Nz, Nx = (1 - Nz)/2 being the demagnetizing factors of the disk (radius R, thickness t):
    Nz = (2R/t) int_0^inf J1(x)^2 (1 - exp(-x t/R))/x^2 dx.
The field-like part of Peff and the Oersted field are neglected in jc.

Usage: stack_torque.py [--evol pillar.evol] [--no_plot]
"""

import argparse

import numpy as np
from scipy.integrate import quad
from scipy.special import j1

E = 1.602176634e-19                  # C
HBAR = 1.054571817e-34               # J s
MU0 = 1.25663706212e-6               # T m/A

# settings.yml
R, T = 20e-9, 4e-9                   # pillar.geo
MS, ALPHA, K = 1.0e6, 0.02, 8.0e5
N0 = 1.344e47
P_VEC = np.array([0.0, 0.0, -1.0])
F1 = {"sigma": 2.0e6, "P": 0.5, "lsf": 5.0e-9, "t": 10.0e-9}
F1N = {"AR": 0.5e-15, "gamma": 0.75, "Gr": 0.55e15, "Gi": 0.0}
NF2 = {"AR": 0.5e-15, "gamma": 0.75, "Gr": 0.55e15, "Gi": 0.0}


def conductance(u, Gpar, Gr, Gi):
    """3x3 conductance matrix of an interface of magnetization u (slonczewskiInterface.h)"""
    ux = np.array([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]])  # ux w = u x w
    uuT = np.outer(u, u)
    return Gpar * uuT + Gr * (np.eye(3) - uuT) - Gi * ux


def stack():
    """polarization and conductances of F1/N (renormalized by F1) and N/F2, in m/s"""
    k = 2 / (E**2 * N0)
    gF = k * F1["sigma"] * (1 - F1["P"]**2) / (F1["lsf"] * np.tanh(F1["t"] / F1["lsf"]))
    G1par = k / F1N["AR"]
    P1 = (F1["P"] * G1par + F1N["gamma"] * gF) / (gF + G1par)
    G1par = gF * G1par / (gF + G1par)
    return (P1, G1par, 2 * k * F1N["Gr"], 2 * k * F1N["Gi"],
            NF2["gamma"], k / NF2["AR"], 2 * k * NF2["Gr"], 2 * k * NF2["Gi"])


def p_eff(m):
    """effective polarization vector for the magnetization m of F2 (effectiveBC)"""
    P1, G1par, G1r, G1i, P2, G2par, G2r, G2i = stack()
    G1 = conductance(P_VEC, G1par, G1r, G1i)
    G2 = conductance(m, G2par, G2r, G2i)
    A = G1 + G2
    return G2 @ np.linalg.solve(A, P1 * P_VEC) + G1 @ np.linalg.solve(A, P2 * m)


def eta(theta):
    """damping-like and field-like efficiencies, Peff along m x (m x p) and m x p divided by
    sin(theta), theta being the angle between m and p"""
    m = np.array([np.sin(theta), 0.0, np.cos(theta)])
    p = np.array([0.0, 0.0, 1.0])               # p along z for the angle; then rotate the result
    # the stack has p = P_VEC: express m in the frame of P_VEC (rotation by pi around x if p = -z)
    if P_VEC[2] < 0:
        m = np.array([m[0], -m[1], -m[2]])
        p = P_VEC
    q = p_eff(m)
    dl = np.cross(m, np.cross(m, p))            # direction m x (m x p), norm sin(theta)
    fl = np.cross(m, p)
    s = np.linalg.norm(dl)
    return -q @ dl / s**2, q @ fl / s**2       # sign: torque towards p is positive


def demag_Nz():
    tau = T / R
    f = lambda x: j1(x)**2 * (1 - np.exp(-tau * x)) / x**2
    val = sum(quad(f, a, a + 50, limit=200)[0] for a in range(0, 2000, 50))
    return 2 / tau * val


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evol", help="feeLLGood .evol file to compare with")
    parser.add_argument("--no_plot", action="store_true")
    args = parser.parse_args()

    P1, G1par, G1r, G1i, P2, G2par, G2r, G2i = stack()
    print(f"stack (m/s): F1/N renormalized by F1: P = {P1:.4f}, Gpar = {G1par:.4e}, "
          f"Gr = {G1r:.4e}; N/F2: P = {P2:.4f}, Gpar = {G2par:.4e}, Gr = {G2r:.4e}")
    thetas = np.radians([1, 30, 60, 90, 120, 150, 179])
    print(" theta (deg) | eta damping-like | eta field-like")
    for th in thetas:
        dl, fl = eta(th)
        print(f" {np.degrees(th):11.0f} | {dl:16.4f} | {fl:14.4f}")

    Nz = demag_Nz()
    Nx = (1 - Nz) / 2
    HK = 2 * K / (MU0 * MS) - (Nz - Nx) * MS
    eta_AP = eta(np.radians(179.9))[0]
    jc = 2 * E / HBAR * ALPHA * MU0 * MS * HK * T / eta_AP
    print(f"demagnetizing factors of the disk: Nz = {Nz:.4f}, Nx = Ny = {Nx:.4f}")
    print(f"mu0 H_K,eff = {MU0 * HK:.4f} T, eta(pi) = {eta_AP:.4f}")
    print(f"critical current density AP -> P (macrospin): jc = {jc:.3e} A/m^2")

    if args.evol:
        cols = None
        rows = []
        for line in open(args.evol):
            if line.startswith('#'):
                if 'columns' in line or line.lstrip('#').strip().startswith('t'):
                    cols = line.lstrip('#').replace('columns:', '').split()
                continue
            if line.strip():
                rows.append([float(v) for v in line.split()])
        data = np.array(rows)
        if cols is None or len(cols) != data.shape[1]:
            cols = ["t", "dt", "max_dm", "<Mx>", "<My>", "<Mz>", "E_ex", "E_aniso", "E_demag",
                    "E_zeeman", "E_tot", "current_density"]
        c = {n: i for i, n in enumerate(cols)}
        t, mz, jn = data[:, c["t"]], data[:, c["<Mz>"]], data[:, c["current_density"]]
        below = np.nonzero(mz < 0)[0]
        if below.size:
            i = below[0]
            print(f"simulation: <Mz> crosses 0 at t = {t[i]:.4e} s, current_density = "
                  f"{jn[i]:.3e} A/m^2 ({jn[i]/jc:.2f} jc)")
        else:
            print(f"simulation: no switching up to t = {t[-1]:.4e} s, current_density = "
                  f"{jn[-1]:.3e} A/m^2 ({jn[-1]/jc:.2f} jc)")
        if not args.no_plot:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(7, 4))
            ax.plot(jn, data[:, c["<Mx>"]], label="<Mx>")
            ax.plot(jn, data[:, c["<My>"]], label="<My>")
            ax.plot(jn, mz, label="<Mz>")
            ax.axvline(jc, color="k", ls="--", label="jc (macrospin)")
            ax.set_xlabel("current_density [A/m^2]")
            ax.set_ylabel("average magnetization")
            ax.legend()
            fig.tight_layout()
            fig.savefig("pillar_switching.png", dpi=120)
            print("figure: pillar_switching.png")


if __name__ == "__main__":
    main()
