"""
Results-section sweep of the anode base case (sec:res-anode-sweep, user 17 Oct 2026), radial loss at R = 3.5 mm:
the electrode quantities as the bulk temperature (T_bulk = 13-17 kK at |J_app| = 4e6 A/m^2) and the
current density (|J_app| = 2e6-1e7 A/m^2 at T_bulk = 15 000 K) are varied.

    python anode_sweep.py [-o outdir]

Runs: PlasmeX runs/Tedge1000/anode_sweep_R3.5mm/T<T>_J<J> (launch.sh anode_sweep_R there), each restarted
from the base case runs/Tedge1000/anode_T15000_J4e6_L0.5mm_R3.5mm with only arc_Te / arc_Japp changed
(prob.use_radial_loss = 4, R = 3.5 mm, T_edge = 1000 K, L = 0.5 mm, W = 4.65 eV, T_wall = 1000 K, steady_tol 1e-7;
T15000_J2e6 1e-5). The earlier no-loss sweep (runs/Tedge1000/anode_sweep, base _noloss) is no longer used. Writes, each as
.pdf and .png:
    anode_sweep     left column against T_bulk, right column against |J_app|; rows as base_sweep.py:
                    q_w, j_i/j at the surface, phi_s and phi(L), sheath width x_s
and anode_sweep.csv (base_sweep.quantities, plus n_e at x = L and T_e at 300 um from the electrode).
The base case is drawn as a filled marker.
"""
import argparse
import csv
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from base_case import HERE, C_X, load, rc, sci_tex, style
from base_sweep import RUNS, ROWS, KSCALE, LEG, COL, last_plotfile, quantities

TB, JB = 15000.0, 4e6
BASE = "Tedge1000/anode_T15000_J4e6_L0.5mm_R3.5mm"
CASES = ([(T, JB, f"Tedge1000/anode_sweep_R3.5mm/T{T}_J4e6") for T in (13000, 14000, 16000, 17000)]
         + [(TB, float(J), f"Tedge1000/anode_sweep_R3.5mm/T15000_J{J}") for J in ("2e6", "6e6", "8e6", "1e7")])


def collect():
    rows = []
    for T, J, sub in [(TB, JB, BASE)] + CASES:
        pf = last_plotfile(os.path.join(RUNS, sub))
        if pf is None:
            print(f"  T={T:.0f} J={J:.1e}: no steady plotfile in {sub}, skipped")
            continue
        s = load(pf)
        q = quantities(s)
        q.update(T_bulk=T, J_app=J, run=sub, base=(sub == BASE), ne_L=s["ne"][-1], Te_300=float(
            np.interp(3e-4, s["d"], s["Te"])))
        flag = "" if q["steady"] < 2e-3 else "   <-- NOT STEADY?"
        print(f"  T={T:7.0f} J={J:.2e}  {sub}/{q['plotfile']}  t={q['t']:.3e}"
              f"  steady {q['steady']:.1e}  j spread {q['j_spread']:.1e}{flag}")
        rows.append(q)
    return rows


def fig_sweep(rows, stem):
    rc()
    fig, axs = plt.subplots(len(ROWS), 2, figsize=(7.0, 7.1), sharex="col",
                            gridspec_kw=dict(hspace=0.32, wspace=0.28))
    tags = "abcdefghijklmnop"
    for c, var in enumerate(("T_bulk", "J_app")):
        sel = [r for r in rows if (r["J_app"] == JB if var == "T_bulk" else r["T_bulk"] == TB)]
        sel.sort(key=lambda r: r[var])
        x = np.array([r[var] for r in sel]) / (1e3 if var == "T_bulk" else 1e6)
        for r_i, (keys, ylab) in enumerate(ROWS):
            ax = axs[r_i, c]
            for k in keys:
                y = np.array([r[k] for r in sel]) / KSCALE.get(k, 1.0)
                col = COL.get(k, C_X)
                ax.plot(x, y, "-", color=col, lw=1.0, label=LEG.get(k))
                ax.plot(x, y, "o", ms=3.8, mfc="white", mec=col, mew=0.9)
                b = np.array([r["base"] for r in sel])
                ax.plot(x[b], y[b], "o", ms=3.8, color=col)
            if len(keys) > 1:
                ax.legend(loc="best", frameon=False, handlelength=1.6)
            if r_i == 1:
                ax.axhline(0, color="#b0b0b0", lw=0.5)
            ax.set_ylabel(ylab)
            style(ax, xlog=False)
            ax.text(0.0, 1.03, f"({tags[2 * r_i + c]})", transform=ax.transAxes, ha="left",
                    va="bottom", fontsize=8.5)
        if var == "T_bulk":
            axs[0, c].set_title(f"$|J_{{\\rm app}}|={sci_tex(JB)}$ A m$^{{-2}}$", fontsize=9, pad=14)
            axs[-1, c].set_xlabel(r"$T_{\rm bulk}$ (10$^3$ K)")
        else:
            axs[0, c].set_title(f"$T_{{\\rm bulk}}={TB:.0f}$ K", fontsize=9, pad=14)
            axs[-1, c].set_xlabel(r"$|J_{\rm app}|$ (10$^6$ A m$^{-2}$)")
    for ext in ("pdf", "png"):
        fig.savefig(f"{stem}.{ext}", bbox_inches="tight", dpi=300)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("-o", "--outdir", default=HERE)
    a = ap.parse_args(argv)
    rows = collect()
    os.makedirs(a.outdir, exist_ok=True)
    keys = ["T_bulk", "J_app", "q_w", "ji_j", "phi_s", "x_s", "phi_L", "Te_w", "ne_w", "ni_w", "E_w",
            "ne_L", "Te_300", "steady", "j_spread", "t", "run", "plotfile"]
    with open(os.path.join(a.outdir, "anode_sweep.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for r in sorted(rows, key=lambda r: (r["J_app"], r["T_bulk"])):
            w.writerow([f"{r[k]:.5g}" if isinstance(r[k], float) else r[k] for k in keys])
    print("\n" + " ".join(f"{k:>9s}" for k in keys[:13]))
    for r in sorted(rows, key=lambda r: (r["J_app"], r["T_bulk"])):
        print(" ".join(f"{r[k]:9.4g}" if isinstance(r[k], float) else f"{str(r[k]):>9s}" for k in keys[:13]))
    fig_sweep(rows, os.path.join(a.outdir, "anode_sweep"))
    print(f"\nwrote anode_sweep (.pdf/.png) and anode_sweep.csv in {a.outdir}")


if __name__ == "__main__":
    main()
