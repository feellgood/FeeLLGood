#!/usr/bin/env python3
"""
Figures of the documents slonczewski_angle.tex and slonczewski_angle_french.tex:
- slonczewski_angle_torque.pdf: damping-like and field-like accumulations integrated over the
  thickness of F2 against theta, reference (lines) and feeLLGood (symbols), models sink and series,
  with the fit of the damping-like torque to the form of Slonczewski (sink);
- slonczewski_angle_profile.pdf: components of s across F2 at theta = 90 deg (sink), reference
  (lines) and feeLLGood at 200 randomly chosen nodes (symbols).

Usage: plot_figures.py, from any directory, after slonczewski_angle.py (it reads thin_layer.msh and
the files <model>_theta<theta>_iter0.sol).
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import slonczewski_angle as sa

NB_NODES = 200                       # number of nodes plotted in the profile


def integral(z, s):
    """integral over the thickness of the columns of s, trapezoidal rule"""
    return ((s[1:] + s[:-1]) / 2 * np.diff(z)[:, None]).sum(axis=0)


def torque_figure(xyz, tets):
    z = np.linspace(0, sa.T2, 601)
    thetas = np.linspace(0, 180, 91)
    fig, ax = plt.subplots(figsize=(6, 3.8))
    for model, color in [("sink", "C0"), ("series", "C1")]:
        I_ref = np.array([integral(z, sa.reference(t, z, model)) for t in thetas])
        ax.plot(thetas, I_ref[:, 1] * 1e9, color, lw=1.5, label=f"damping-like, {model}")
        ax.plot(thetas, I_ref[:, 0] * 1e9, color, lw=1.0, ls="--", label=f"field-like, {model}")
        I = np.array([sa.thickness_integral(xyz, tets, sa.read_s(f"{model}_theta{t}_iter0.sol"))
                      for t in sa.ANGLES])
        ax.plot(sa.ANGLES, I[:, 1] * 1e9, "o", mfc="none", color=color)
        ax.plot(sa.ANGLES, I[:, 0] * 1e9, "s", ms=4, mfc="none", color=color)
        if model == "sink":
            # sin(theta)/int s_y = alpha + beta cos(theta), fitted at the angles of
            # slonczewski_angle.py
            th = np.radians([t for t in sa.ANGLES if 0 < t < 180])
            I_fit = np.array([integral(z, sa.reference(t, z, model))[1]
                              for t in np.degrees(th)])
            beta, alpha = np.polyfit(np.cos(th), np.sin(th) / I_fit, 1)
            Lambda2 = (alpha + beta) / (alpha - beta)
            th = np.radians(thetas)
            ax.plot(thetas, np.sin(th) / (alpha + beta * np.cos(th)) * 1e9, "k:", lw=1,
                    label=f"Slonczewski's form, sink ($\\Lambda^2$ = {Lambda2:.2f})")
    ax.set_xlim(0, 180)
    ax.set_ylim(top=12.5)                # room for the legend
    ax.set_xticks(range(0, 181, 30))
    ax.set_xlabel(r"$\theta$ (°)")
    ax.set_ylabel(r"$\int s\,dz$ (nA)")
    ax.legend(fontsize=7.5, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig("slonczewski_angle_torque.pdf")


def profile_figure(xyz):
    z = np.linspace(0, sa.T2, 601)
    s_ref = sa.reference(90, z, "sink")
    s = sa.read_s("sink_theta90_iter0.sol")
    nodes = np.random.default_rng(0).choice(len(xyz), NB_NODES, replace=False)
    fig, ax = plt.subplots(figsize=(6, 3.4))
    for k, (name, color) in enumerate([("$s_x$", "C3"), ("$s_y$", "C0"), ("$s_z$", "C2")]):
        ax.plot(z * 1e9, s_ref[:, k], color, lw=1.5, label=name)
        ax.plot(xyz[nodes, 2] * 1e9, s[nodes, k], "o", ms=3, mfc="none", color=color)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlim(0, sa.T2 * 1e9)
    ax.set_xlabel("$z$ (nm)")
    ax.set_ylabel("$s$ (A/m)")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig("slonczewski_angle_profile.pdf")


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    xyz, tets = sa.read_mesh("thin_layer.msh")
    torque_figure(xyz, tets)
    profile_figure(xyz)


if __name__ == "__main__":
    main()
