#!/usr/bin/env python3
"""
Steady anode-layer profiles from PlasmeX runs that differ only in prob.rxn_rates.

    python3 compare_steady.py hoffert_lien=<plt> almeida=<plt> saha=<plt> [--out anode_steady]

Each plotfile is read on its AMR leaf cells with PlasmeX exe/bohm.py::sample (loaded live through
figures/validation_fig/plasmex.py). The Saha density n_S(T_e, T_h, p) is almeida_arc_bc.saha_ne,
with the same g+/g1 = 6, I = 15.75 eV as the code's analytic rates. Writes <out>.{pdf,png} and
prints the values quoted in reaction_rates.md section 5.
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "figures", "validation_fig"))
from plasmex import exe_module  # noqa: E402

B = exe_module("bohm")
A = exe_module("almeida_arc_bc")

COL = {"hoffert_lien": "#2a78d6", "almeida": "#eb6834", "saha": "#1baf7a"}
LAB = {"hoffert_lien": "Hoffert & Lien", "almeida": "Benilov & Naidis (Almeida)",
       "saha": "Hinnov & Hirschberg (MK / Pfender)"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def load(fn):
    d = B.sample(fn, want_derived=False)
    knobs, _ = B.read_knobs(d["dir"], None, fn)
    p = knobs["prob.pres"]
    d["dist"] = d["hi"] - d["x"]                                  # distance from the anode
    d["nS"] = np.array([A.saha_ne(te, th, p) for te, th in zip(d["ETEMP"], d["TGAS"])])
    return d


def style(ax):
    ax.grid(True, which="major", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.set_xscale("log")
    ax.set_xlim(1e-7, 5e-4)
    ax.set_xlabel("distance from anode (m)", color=INK)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+", help="model=plotfile")
    ap.add_argument("--out", default=os.path.join(HERE, "anode_steady"))
    a = ap.parse_args()
    runs = [r.split("=", 1) for r in a.runs]

    fig, axs = plt.subplots(2, 2, figsize=(10, 7.4), constrained_layout=True)
    (ax1, ax2), (ax3, ax4) = axs
    print(f"{'model':14s} {'n_e wall':>10s} {'n_e/n_S wall':>12s} {'T_e wall':>9s} {'T_h wall':>9s} "
          f"{'n_e 10um':>10s} {'n_e/n_S 10um':>12s} {'n_e/n_S 100um':>13s} {'phi_wall-phi_arc':>16s}")
    for m, fn in runs:
        d = load(fn)
        c = COL.get(m, INK)
        x = d["dist"]
        ax1.loglog(x, d["NE"], color=c, lw=2, label=LAB.get(m, m))
        ax2.loglog(x, d["NE"] / d["nS"], color=c, lw=2)
        ax3.semilogx(x, d["ETEMP"], color=c, lw=2)
        ax3.semilogx(x, d["TGAS"], color=c, lw=1.5, ls=":")
        ax4.semilogx(x, d["PHI"] - d["PHI"][0], color=c, lw=2)
        r = d["NE"] / d["nS"]
        i10, i100 = np.argmin(abs(x - 1e-5)), np.argmin(abs(x - 1e-4))
        print(f"{m:14s} {d['NE'][-1]:10.3e} {r[-1]:12.3e} {d['ETEMP'][-1]:9.0f} {d['TGAS'][-1]:9.0f} "
              f"{d['NE'][i10]:10.3e} {r[i10]:12.3e} {r[i100]:13.3e} {d['PHI'][-1] - d['PHI'][0]:16.3f}")

    ax1.set_ylabel(r"$n_e$ (m$^{-3}$)")
    ax1.set_title("(a) Electron density", loc="left", fontsize=10, color=INK)
    ax1.legend(fontsize=8, frameon=False, loc="lower right", labelcolor=INK)
    ax2.axhline(1.0, color=INK2, lw=0.8)
    ax2.set_ylabel(r"$n_e/n_S(T_e,T_h)$")
    ax2.set_title("(b) Departure from Saha at the local $T_e$", loc="left", fontsize=10, color=INK)
    ax3.set_ylabel("temperature (K)")
    ax3.set_title(r"(c) $T_e$ (solid), $T_h$ (dotted)", loc="left", fontsize=10, color=INK)
    ax4.set_ylabel(r"$\phi-\phi_{\rm arc}$ (V)")
    ax4.set_title("(d) Potential", loc="left", fontsize=10, color=INK)
    for ax in axs.flat:
        style(ax)
        ax.yaxis.label.set_color(INK)
    for ext in ("pdf", "png"):
        fig.savefig(f"{a.out}.{ext}", dpi=200)


if __name__ == "__main__":
    main()
