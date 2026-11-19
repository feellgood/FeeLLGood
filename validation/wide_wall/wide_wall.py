#!/usr/bin/env python3
"""
Adiabatic and non-adiabatic spin-transfer torques on a domain wall, against the local limit of
Zhang-Li and its first non-local correction (tools/fg-torques, dev-notes/adiabatic_torques.tex).

Wire along z (wire.geo, 4 x 4 x 300 nm, h = 1 nm), uniform current along the wire (jn on the
inlet z = 0, V = 0 and s = 0 on the outlet z = L), and a wall of width delta at its center,
imposed as the initial magnetization:
    theta(z) = 2 atan(exp((z - z0)/delta)),    m = (sin theta, 0, cos theta).
The magnetization is not relaxed: the spin accumulation of the iteration 0 is computed with this
profile, and fg-torques fits the Zhang-Li form T = eta [-(u.grad) m + beta m x (u.grad) m] to the
torque -m x s/(tau_sd Ms) (least squares, weights |(u.grad) m|^2).

Local limit (texture wide compared with l_sd), with xi = tau_sd/tau_sf = (1 - P P') l_sd^2/l_sf^2:
    beta = xi,    eta = 1/(1 + xi^2).
First non-local correction: the transverse diffusion D0 s'' of the spin accumulation of the local
limit, s0 = -K u theta' [xi e_theta + y], K = tau_sd Ms/(1 + xi^2), e_theta = dm/dtheta, y = m x
e_theta, adds lsd^2 s0'' to the source of the next order. Its dominant part, -K u theta''' y,
gives the component -lsd^2 K u theta''' along e_theta, hence
    beta_loc(z) = xi - lsd^2 theta'''/theta' = xi - (lsd^2/delta^2)(1 - 2 sech^2((z - z0)/delta)),
and, averaged with the weight theta'^2 (sech^2) of the fit,
    beta_fit = xi + lsd^2/(3 delta^2).
The correction is absolute: the same for any P'. It is 1.3 % of xi for delta = 25 nm and P' = 0,
5.3 % for 12.5 nm, twice as much with P' = P.

The script checks, for delta = 12.5 and 25 nm and P' = 0 and P:
  - beta_fit against xi + lsd^2/(3 delta^2), the error being relative to the non-local correction
    lsd^2/(3 delta^2) (tolerance TOL_CORR): the test separates the first-order theory from the
    local limit, the error of which would be 1 in these units; the next order, about
    (lsd/delta)^4/xi, is the largest error (1.6 % for delta = 12.5 nm, P' = P);
  - beta_loc at the center of the wall (nodal value of fg-torques) against xi + lsd^2/delta^2, the
    error being relative to lsd^2/delta^2 (tolerance TOL_CORR);
  - eta_fit against 1/(1 + xi^2), tolerance TOL_ETA (discretization error, as (h/delta)^2).

Usage: wide_wall.py [path to feellgood], from any directory; the mesh is built with gmsh if
missing. Exit status 0 if all the checks pass.
"""

import json
import os
import subprocess
import sys

import numpy as np

TOL_CORR = 5e-2                      # tolerance on the non-local correction of beta
TOL_ETA = 2e-3                       # tolerance on eta_fit

L, Z0 = 300e-9, 150e-9               # wire length (wire.geo), center of the wall (m)
F = {"sigma": 5e6, "N0": 1.344e47, "P": 0.7, "lsf": 5e-9, "lsd": 1e-9}
JN = 1e12                            # A/m^2 on the inlet
DELTAS = [12.5e-9, 25e-9]            # wall widths (m): z0 = 6 delta at least from the ends


def settings(name, delta, p_prime):
    # the coordinates given to the function are in meters; Math.cosh/tanh do not exist (ES5)
    profile = (f"function(x, y, z) {{ var e = Math.exp((z - {Z0!r}) / {delta!r}); "
               "return [2 / (e + 1 / e), 0, (1 / e - e) / (e + 1 / e)]; }")
    return {
        "outputs": {"file_basename": name, "final_time": 1e-13, "evol_time_step": 1e-12,
                    "mag_config_every": 10},
        "mesh": {
            "filename": "wire.msh",
            "length_unit": 1e-9,
            "volume_regions": {"wire": {
                "Js": 1.0, "Ae": 1.3e-11, "alpha_LLG": 0.05, "K": 0,
                "sigma": F["sigma"], "dens_state": F["N0"], "P": F["P"], "P_prime": p_prime,
                "l_sd": F["lsd"], "l_sf": F["lsf"]}},
            "surface_regions": {
                "inlet": {"jn": JN},
                "outlet": {"V": 0.0, "s": [0, 0, 0]}}},
        "initial_magnetization": profile,
        "spin_accumulation": {"enable": True},
        "time_integration": {"max(du)": 0.02, "min(dt)": 5e-20, "max(dt)": 5e-13}}


def read_columns(file_name):
    """columns of a file with a '## columns:' line, as a dict name -> array of strings"""
    names, rows = None, []
    for line in open(file_name):
        if line.startswith('## columns:'):
            names = line.split(':', 1)[1].split()
        elif line.strip() and not line.startswith('#'):
            rows.append(line.split())
    return {n: np.array([r[i] for r in rows]) for i, n in enumerate(names)}


def read_mesh_z(file_name):
    """z of the nodes of a gmsh 4.1 text mesh, in file order (the order of the .sol files)"""
    lines = open(file_name).read().split('\n')
    i = lines.index('$Nodes')
    nb_blocks = int(lines[i + 1].split()[0])
    k, z = i + 2, []
    for _ in range(nb_blocks):
        n = int(lines[k].split()[3])
        k += 1 + n
        z += [float(lines[k + j].split()[2]) for j in range(n)]
        k += n
    return np.array(z) * 1e-9


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    feellgood = sys.argv[1] if len(sys.argv) > 1 else "../../src/feellgood"
    fg_torques = os.path.join("..", "..", "tools", "fg-torques")
    if not os.path.exists("wire.msh"):
        subprocess.run(["gmsh", "-3", "wire.geo", "-o", "wire.msh"], check=True,
                       stdout=subprocess.DEVNULL)
    z = read_mesh_z("wire.msh")
    lsd2 = F["lsd"]**2

    ok = True
    print("errors of beta relative to the non-local correction, of eta absolute")
    print("delta    P' | xi       correction/xi | beta_fit  first order  error   "
          "| beta(z0)  first order  error   | eta_fit  1/(1+xi^2) error")
    for delta in DELTAS:
        for p_prime in [0.0, F["P"]]:
            name = f"wall_{delta * 1e9:g}nm_Pprime{p_prime:g}"
            with open(name + ".json", "w") as f:
                json.dump(settings(name, delta, p_prime), f)
            subprocess.run([feellgood, name + ".json"], check=True, stdout=subprocess.DEVNULL)
            subprocess.run([sys.executable, fg_torques, name + ".json", name + "_iter0.sol",
                            "-o", "torques.txt"], check=True, stdout=subprocess.DEVNULL)
            table = read_columns("torques.txt")
            eta, beta = float(table["eta"][0]), float(table["beta"][0])
            xi = float(table["beta_local"][0])
            xi_ref = (1 - F["P"] * p_prime) * lsd2 / F["lsf"]**2
            nodal = read_columns(name + "_iter0_torques.sol")
            center = np.abs(z - Z0) < 0.25e-9
            beta_center = np.mean(nodal["beta"][center].astype(float))

            beta_ref = xi + lsd2 / (3 * delta**2)
            center_ref = xi + lsd2 / delta**2
            eta_ref = 1 / (1 + xi**2)
            err_beta = abs(beta - beta_ref) / (lsd2 / (3 * delta**2))
            err_center = abs(beta_center - center_ref) / (lsd2 / delta**2)
            err_eta = abs(eta - eta_ref)
            correction = lsd2 / (3 * delta**2) / xi
            print(f"{delta * 1e9:4.1f} nm {p_prime:3.1f} | {xi:.6f} {correction:9.1e}"
                  f"     | {beta:.6f}  {beta_ref:.6f}"
                  f"     {err_beta:.1e} | {beta_center:.6f}  {center_ref:.6f}     {err_center:.1e}"
                  f" | {eta:.5f}  {eta_ref:.5f}    {err_eta:.1e}")
            ok = ok and abs(xi / xi_ref - 1) < 1e-12
            ok = ok and err_beta < TOL_CORR and err_center < TOL_CORR and err_eta < TOL_ETA
    os.remove("torques.txt")
    print("OK" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
