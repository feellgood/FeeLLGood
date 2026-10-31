#!/usr/bin/env python3
"""
Macrospin estimate of the spin-transfer switching of the example pillar_tunnel_voltage_ramp1ns, from
the parameters of settings.yml, and comparison with the simulation (pillar.evol, pillar_iter*.sol).

The tunnel barrier is the one of feeLLGood (src/tunnelJunction.h): conductance per unit area
G_T = (1 + c cos theta)/((1 + c) RA), c = TMR/(2 + TMR), theta being the angle between m and p, and
spin current entering F2 (mu_B/e) jn P_T, P_T = (P1 p + P2 m)/(1 + c cos theta), P1 P2 = c. Its
transverse part, absorbed by F2, gives the damping-like torque per unit area
(hbar/2e) jn eta(theta) sin(theta), eta(theta) = P1/(1 + c cos theta): per unit current, the torque
is larger near the antiparallel state (Lambda^2 = 1 + TMR in the form of Slonczewski); at fixed
voltage, jn = G_T V, it is proportional to sin(theta).

Critical current density of the instability of the antiparallel state (m = -p, AP -> P), for a
perpendicular free layer of thickness t (macrospin, linearized LLG):
    jc = (2e/hbar) alpha mu0 Ms H_K,eff t / eta(pi),
    H_K,eff = 2K/(mu0 Ms) - (Nz - Nx) Ms,
Nz, Nx = (1 - Nz)/2 being the demagnetizing factors of the disk (radius R, thickness t), and
critical voltage Vc = jc (RA_AP + t/sigma), RA_AP = RA (1 + TMR), the free layer being in series
with the barrier. The spin accumulation of F2 and the field-like torque are neglected.

With the simulation, the script gives the time and the voltage of the switching (<Mz> crossing 0,
pillar.evol) and, from the .sol files, the mean current density through the barrier
jn = G_T(m)(V - V1), averaged over the nodes of the bottom face, and the resistance-area product
-V1/jn: it jumps from RA_AP to RA_P at the switching.

Usage: junction_torque.py [--evol pillar.evol] [--no_plot]
"""

import argparse
import glob
import re

import numpy as np
from scipy.integrate import quad
from scipy.special import j1

E = 1.602176634e-19                  # C
HBAR = 1.054571817e-34               # J s
MU0 = 1.25663706212e-6               # T m/A

# settings.yml
R, T = 20e-9, 4e-9                   # pillar.geo
MS, ALPHA, K = 1.0e6, 0.02, 8.0e5
SIGMA = 2.0e6
P_VEC = np.array([0.0, 0.0, -1.0])
RA, TMR = 0.5e-12, 1.5
C = TMR / (2 + TMR)
P1 = np.sqrt(C)                      # default of feeLLGood, symmetric junction
VOLTAGE = lambda t: -0.8 * t / 1e-9  # voltage(t) of settings.yml


def eta(theta):
    """damping-like efficiency per unit current, theta being the angle between m and p"""
    return P1 / (1 + C * np.cos(theta))


def demag_Nz():
    tau = T / R
    f = lambda x: j1(x)**2 * (1 - np.exp(-tau * x)) / x**2
    val = sum(quad(f, a, a + 50, limit=200)[0] for a in range(0, 2000, 50))
    return 2 / tau * val


def read_evol(fileName):
    """columns of a .evol file, as a dictionary"""
    cols, rows = None, []
    for line in open(fileName):
        if line.startswith('#'):
            if 'columns' in line:
                cols = line.lstrip('#').replace('columns:', '').split()
            continue
        if line.strip():
            rows.append([float(v) for v in line.split()])
    data = np.array(rows)
    return {n: data[:, i] for i, n in enumerate(cols)}


def barrier_current(fileName):
    """time and mean current density through the barrier (A/m^2) of a .sol file, averaged over
    the nodes of the bottom face (z = -t/2)"""
    t = float(next(l for l in open(fileName) if l.startswith('## time:')).split()[-1])
    rows = np.array([[float(v) for v in l.split()] for l in open(fileName)
                     if l.strip() and not l.startswith('#')])
    xyz = read_nodes("pillar.msh")
    bottom = np.abs(xyz[:, 2] + T / 2) < 1e-12
    m, V = rows[bottom, 1:4], rows[bottom, 8]
    G = (1 + C * (m @ P_VEC)) / ((1 + C) * RA)
    return t, (G * (V - VOLTAGE(t))).mean()


def read_nodes(fileName):
    """node coordinates (m) of a gmsh 4.1 text mesh, in file order (the order of the .sol files)"""
    lines = open(fileName).read().split('\n')
    i = lines.index('$Nodes')
    nb_blocks = int(lines[i + 1].split()[0])
    k, xyz = i + 2, []
    for _ in range(nb_blocks):
        n = int(lines[k].split()[3])
        k += 1 + n
        xyz += [[float(v) for v in lines[k + j].split()] for j in range(n)]
        k += n
    return np.array(xyz) * 1e-9


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evol", help="feeLLGood .evol file to compare with")
    parser.add_argument("--no_plot", action="store_true")
    args = parser.parse_args()

    print(f"barrier: RA_P = {RA:.2e} Ohm m^2, RA_AP = {RA * (1 + TMR):.2e} Ohm m^2, "
          f"c = TMR/(2 + TMR) = {C:.4f}, P1 = P2 = {P1:.4f}")
    print(" theta (deg) | eta damping-like (per unit current)")
    for th in [0, 30, 60, 90, 120, 150, 180]:
        print(f" {th:11d} | {eta(np.radians(th)):.4f}")

    Nz = demag_Nz()
    Nx = (1 - Nz) / 2
    HK = 2 * K / (MU0 * MS) - (Nz - Nx) * MS
    jc = 2 * E / HBAR * ALPHA * MU0 * MS * HK * T / eta(np.pi)
    Vc = jc * (RA * (1 + TMR) + T / SIGMA)
    print(f"demagnetizing factors of the disk: Nz = {Nz:.4f}, Nx = Ny = {Nx:.4f}")
    print(f"mu0 H_K,eff = {MU0 * HK:.4f} T, eta(pi) = {eta(np.pi):.4f}")
    print(f"critical current density AP -> P (macrospin): jc = {jc:.3e} A/m^2, "
          f"critical voltage |Vc| = {Vc:.3f} V")
    print(f"Oersted field at the edge, mu0 jc R/2 = {MU0 * jc * R / 2 * 1e3:.2f} mT")

    if args.evol:
        d = read_evol(args.evol)
        t, mz, v = d["t"], d["<Mz>"], d["voltage"]
        below = np.nonzero(mz < 0)[0]
        if below.size:
            i = below[0]
            print(f"simulation: <Mz> crosses 0 at t = {t[i]:.4e} s, voltage = {v[i]:.3f} V "
                  f"({abs(v[i]) / Vc:.2f} |Vc|)")
        else:
            print(f"simulation: no switching up to t = {t[-1]:.4e} s, voltage = {v[-1]:.3f} V")
        sols = sorted(glob.glob("pillar_iter*.sol"),
                      key=lambda f: int(re.search(r"iter(\d+)", f).group(1)))
        cur = np.array([barrier_current(f) for f in sols if VOLTAGE(barrier_current(f)[0]) != 0])
        if cur.size:
            print("   t (s)    | voltage (V) | jn (A/m^2) | RA = -V1/jn (Ohm m^2)")
            for tt, jn in cur:
                print(f" {tt:.3e} | {VOLTAGE(tt):11.3f} | {jn:.3e} | {-VOLTAGE(tt) / jn:.3e}")
        if not args.no_plot:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, (ax, ax2) = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
            ax.plot(-v, d["<Mx>"], label="<Mx>")
            ax.plot(-v, d["<My>"], label="<My>")
            ax.plot(-v, mz, label="<Mz>")
            ax.axvline(Vc, color="k", ls="--", label="|Vc| (macrospin)")
            ax.set_ylabel("average magnetization")
            ax.legend()
            if cur.size:
                ax2.plot(-VOLTAGE(cur[:, 0]), -VOLTAGE(cur[:, 0]) / cur[:, 1] * 1e12, "o-")
            ax2.axhline(RA * (1 + TMR) * 1e12, color="C3", ls=":", label="RA_AP")
            ax2.axhline(RA * 1e12, color="C2", ls=":", label="RA_P")
            ax2.axvline(Vc, color="k", ls="--")
            ax2.set_xlabel("-voltage (V)")
            ax2.set_ylabel(r"RA ($\Omega\,\mu$m$^2$)")
            ax2.legend()
            fig.tight_layout()
            fig.savefig("pillar_switching.png", dpi=120)
            print("figure: pillar_switching.png")


if __name__ == "__main__":
    main()
