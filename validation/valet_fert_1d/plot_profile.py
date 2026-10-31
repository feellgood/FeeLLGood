#!/usr/bin/env python3
"""
Figure valet_fert_1d_profile.pdf of the documents valet_fert_1d.tex and valet_fert_1d_french.tex:
spin accumulation s_z along the wire, Valet-Fert solution (lines) against feeLLGood at 250
randomly chosen nodes (circles), configurations P and AP.

Usage: plot_profile.py, from any directory, after valet_fert.py (it reads wire.msh, P_iter0.sol
and AP_iter0.sol).
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import valet_fert as vf

Z_MAX = 60e-9                        # m, end of the plotted range
NB_NODES = 250                       # number of nodes plotted per configuration


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    z = vf.read_mesh_nodes("wire.msh")[:, 2] * 1e-9
    zz = np.linspace(0, vf.L, 2000)
    nodes = np.where(z < Z_MAX)[0]
    nodes = np.random.default_rng(0).choice(nodes, NB_NODES, replace=False)

    fig, ax = plt.subplots(figsize=(6, 3.6))
    for sign_F1, name, color in [(1, "P", "C0"), (-1, "AP", "C3")]:
        s = vf.read_s(name + "_iter0.sol")
        ax.plot(zz * 1e9, vf.valet_fert(sign_F1, zz), color, lw=1.5, label=f"{name}: Valet–Fert")
        ax.plot(z[nodes] * 1e9, s[nodes, 2], "o", ms=3.5, mfc="none", color=color,
                label=f"{name}: feeLLGood")
    ax.set_xlim(0, Z_MAX * 1e9)
    ax.set_xlabel("$z$ (nm)")
    ax.set_ylabel("$s_z$ (A/m)")
    ax.axhline(0, color="k", lw=0.5)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig("valet_fert_1d_profile.pdf")


if __name__ == "__main__":
    main()
