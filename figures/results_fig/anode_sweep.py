"""
Results-section sweep of the anode base case (sec:res-anode-sweep, user 17 Oct 2026), radial loss at R = sqrt(200 A/(pi |J_app|)) (23 Oct 2026; R = 3.5 mm before):
the anode layer as the bulk temperature (T_bulk = 13-17 kK at |J_app| = 4e6 A/m^2) and the
current density (|J_app| = 2e6-1e7 A/m^2 at T_bulk = 15 000 K) are varied.

    python anode_sweep.py [-o outdir]

Runs: PlasmeX runs/wannier_rel/anode_sweep_R3.5mm/T<T>_J<J> (launch.sh sweep there), each restarted
from its 17 Oct state (runs/Tedge1000/anode_sweep_R3.5mm) with plasmex.wannier_drift = 1 added
(prob.use_radial_loss = 4, R = 3.5 mm, T_edge = 1000 K, L = 0.5 mm, W = 4.65 eV, T_wall = 1000 K, steady_tol 1e-7;
T15000_J2e6 1e-5). Writes, each as .pdf and .png:
    anode_sweep     left column the T_bulk series, right column the |J_app| series (18 Oct 2026, user: show
                    which profiles follow T_bulk and which the current; j_i instead of j_i/j; compact):
                    (a)-(f) profiles against the distance from the anode (base_case's log-linear axis, anode
                    on the right): T_e and T_h, n_e, phi with the sheath edge (x_s, phi_s) marked on each
                    curve; (g)-(n) q_w, the ion current density j_i into the anode, the sheath width x_s and the sheath
                    drop phi_s against the parameter (x_s, phi_s added 18 Oct 2026, user; J = 2e6 replaced by 3e6).
                    Colour = the case, one ordinal ramp per series (orange T_bulk, blue |J_app|; validated
                    with the dataviz validate_palette.js --ordinal); the base case is drawn thicker.
and anode_sweep.csv (base_sweep.quantities, plus j_i at the surface, n_e at x = L and T_e at 300 um).
"""
import argparse
import csv
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from base_case import (HERE, C_X, away, load, mark_electrode, orient, rc, sci_tex, sheath_edge, side,
                       style, xlabel)
from base_sweep import RUNS, last_plotfile, quantities

TB, JB = 15000.0, 4e6
# 23 Oct 2026 (user): the column radius carries I = 200 A in every case, R = sqrt(I/(pi |J_app|)) (4.0 mm at 4e6;
# 4.6-2.5 mm over the J series), runs/I200 (launch.sh there), each restarted from the wannier_rel state below.
BASE = "I200/anode_T15000_J4e6"
CASES = ([(T, JB, f"I200/anode_sweep/T{T}_J4e6") for T in (13000, 14000, 16000, 17000)]
         + [(TB, float(J), f"I200/anode_sweep/T15000_J{J}") for J in ("3e6", "6e6", "8e6", "1e7")])
# Before 23 Oct (R = 3.5 mm everywhere): wannier_rel/anode_T15000_J4e6_L0.5mm_R3.5mm, wannier_rel/anode_sweep_R3.5mm/*
# 18 Oct 2026 (user): the J = 2e6 case (T15000_J2e6, converged only to steady_tol 1e-5) is no longer in the paper;
# J = 3e6 (T15000_J3e6, launch.sh J3e6 in runs/wannier_rel) was added instead.


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
            np.interp(3e-4, s["d"], s["Te"])), ji=q["ji_j"] * abs(np.median(s["j"])), prof=s)
        flag = "" if q["steady"] < 2e-3 else "   <-- NOT STEADY?"
        print(f"  T={T:7.0f} J={J:.2e}  {sub}/{q['plotfile']}  t={q['t']:.3e}"
              f"  steady {q['steady']:.1e}  j spread {q['j_spread']:.1e}{flag}")
        rows.append(q)
    return rows


# One ordinal ramp per series, light -> dark with the parameter (dataviz validate_palette.js --ordinal:
# monotone lightness, adjacent gaps >= 0.06, light end >= 2:1 on white, one hue). Blue: the palette's
# sequential steps 250-650; orange: the same five lightness steps on the hue of its slot-2 orange.
RAMP = {"T_bulk": ["#f39a78", "#ea7a4f", "#db5d2a", "#b3481c", "#8a3918"],
        "J_app": ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]}


def _series(rows, var):
    sel = [r for r in rows if (r["J_app"] == JB if var == "T_bulk" else r["T_bulk"] == TB)]
    return sorted(sel, key=lambda r: r[var])


def _case_label(r, var):
    return (f"{r['T_bulk'] / 1e3:g} kK" if var == "T_bulk"
            else rf"${sci_tex(r['J_app'])}$").replace(r"\times", r"\!\times\!")


def fig_sweep(rows, stem):
    rc()
    fig = plt.figure(figsize=(7.0, 9.3))
    outer = fig.add_gridspec(2, 1, height_ratios=[3.3, 2.1], hspace=0.17)
    top = outer[0].subgridspec(3, 2, hspace=0.13, wspace=0.28)
    bot = outer[1].subgridspec(2, 4, wspace=0.62, hspace=0.62)
    prof = np.empty((3, 2), dtype=object)
    for c in range(2):
        for r in range(3):
            prof[r, c] = fig.add_subplot(top[r, c], sharex=prof[0, c] if r else None)
    scal = [fig.add_subplot(bot[r, k]) for r in range(2) for k in range(4)]
    tag = "Anode"
    for c, var in enumerate(("T_bulk", "J_app")):
        sel = _series(rows, var)
        cols = RAMP[var][:len(sel)]
        for r, col in zip(sel, cols):
            p = r["prof"]
            d = p["d"]
            lw = 1.7 if r["base"] else 1.0
            lab = _case_label(r, var)
            prof[0, c].plot(d, p["Te"] / 1e3, color=col, lw=lw, label=lab)
            prof[0, c].plot(d, p["Th"] / 1e3, color=col, lw=lw, ls=(0, (4, 2)))
            prof[1, c].plot(d, p["ne"], color=col, lw=lw)
            prof[2, c].plot(d, p["phi"], color=col, lw=lw)
            x_s, phi_s, _ = sheath_edge(d, p["ne"], p["ni"], p["phi"])
            prof[2, c].plot([x_s], [phi_s], "o", ms=3.6, color=col, mec="white", mew=0.5, zorder=4)
        ax = prof[0, c]
        ax.set_ylim(0, 18.5)
        ax.set_ylabel(r"Temperature (10$^3$ K)")
        title = (f"$J_{{\\rm app}}={sci_tex(JB)}$ A m$^{{-2}}$" if var == "T_bulk"
                 else f"$T_{{\\rm bulk}}={TB:.0f}$ K")
        ax.set_title(title, fontsize=9, pad=12)
        h, l = ax.get_legend_handles_labels()       # case key in the empty lower left of the n_e panel
        leg = prof[1, c].legend(h, l, frameon=False, loc="lower left", handlelength=1.6, fontsize=6.8,
                        title=(r"$T_{\rm bulk}$" if var == "T_bulk" else r"$J_{\rm app}$ (A m$^{-2}$)"),
                        title_fontsize=6.8)
        leg._legend_box.align = "left"
        # line-style key (T_e solid, T_h dashed) in neutral ink, so it never reads as a case
        ax.text(0.03, 0.06, r"solid $T_e$, dashed $T_h$", transform=ax.transAxes, fontsize=6.8,
                color="#444444", ha="left", va="bottom")
        prof[1, c].set_yscale("log")
        prof[1, c].set_ylim(5e19, 5e23)
        prof[1, c].set_ylabel(r"$n_e$ (m$^{-3}$)")
        prof[2, c].axhline(0, color="#b0b0b0", lw=0.5)
        prof[2, c].set_ylim(-0.2, 6.6)
        prof[2, c].set_ylabel(r"$\phi$ (V)")
        prof[2, c].text(0.03, 0.06, r"$\bullet$ sheath edge $(x_{\rm s},\phi_{\rm s})$",
                        transform=prof[2, c].transAxes, fontsize=6.8, color="#444444", va="bottom")
        prof[2, c].set_xlabel(xlabel(tag))
        for r_i in range(3):
            a = prof[r_i, c]
            style(a)
            if r_i < 2:
                plt.setp(a.get_xticklabels(), visible=False)
            mark_electrode(a, side(tag), label=tag.lower() if r_i == 0 else None)
            a.text(0.0 if away(tag) == "left" else 1.0, 1.015, f"({'abcdef'[2 * r_i + c]})",
                   transform=a.transAxes, fontsize=9, va="bottom", ha=away(tag))
        orient(prof[0, c], tag)                   # sharex: the whole column follows

        # electrode quantities against the series parameter, points in the profiles' colours
        x = np.array([r[var] for r in sel]) / (1e3 if var == "T_bulk" else 1e6)
        xl = r"$T_{\rm bulk}$ (10$^3$ K)" if var == "T_bulk" else r"$J_{\rm app}$ (10$^6$ A m$^{-2}$)"
        for k, (key, scale, ylab) in enumerate((("q_w", 1.0, r"$q_w$ (MW m$^{-2}$)"),
                                                 ("ji", 1e5, r"$j_i$ (10$^5$ A m$^{-2}$)"),
                                                 ("x_s", 1e-6, r"$x_{\rm s}$ ($\mu$m)"),
                                                 ("phi_s", 1.0, r"$\phi_{\rm s}$ (V)"))):
            a = scal[4 * (k // 2) + 2 * c + k % 2]
            y = np.array([r[key] for r in sel]) / scale
            a.plot(x, y, "-", color="#9a9a9a", lw=0.8, zorder=1)
            for xi, yi, col, r in zip(x, y, cols, sel):
                a.plot([xi], [yi], "o", ms=4.6 if r["base"] else 3.8, color=col, mec="white",
                       mew=0.5, zorder=3)
            a.set_ylim(0, 1.12 * y.max())
            a.set_ylabel(ylab, labelpad=2)
            a.set_xlabel(xl)
            style(a, xlog=False)
            a.text(0.0, 1.03, f"({'ghijklmn'[4 * (k // 2) + 2 * c + k % 2]})", transform=a.transAxes, ha="left",
                   va="bottom", fontsize=9)
    fig.align_ylabels(prof[:, 0])
    fig.align_ylabels(prof[:, 1])
    for ext in ("pdf", "png"):
        fig.savefig(f"{stem}.{ext}", bbox_inches="tight", dpi=300)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("-o", "--outdir", default=HERE)
    a = ap.parse_args(argv)
    rows = collect()
    os.makedirs(a.outdir, exist_ok=True)
    keys = ["T_bulk", "J_app", "q_w", "ji_j", "ji", "phi_s", "x_s", "phi_L", "Te_w", "ne_w", "ni_w", "E_w",
            "ne_L", "Te_300", "steady", "j_spread", "t", "run", "plotfile"]
    with open(os.path.join(a.outdir, "anode_sweep.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for r in sorted(rows, key=lambda r: (r["J_app"], r["T_bulk"])):
            w.writerow([f"{r[k]:.5g}" if isinstance(r[k], float) else r[k] for k in keys])
    print("\n" + " ".join(f"{k:>9s}" for k in keys[:14]))
    for r in sorted(rows, key=lambda r: (r["J_app"], r["T_bulk"])):
        print(" ".join(f"{r[k]:9.4g}" if isinstance(r[k], float) else f"{str(r[k]):>9s}" for k in keys[:14]))
    fig_sweep(rows, os.path.join(a.outdir, "anode_sweep"))
    print(f"\nwrote anode_sweep (.pdf/.png) and anode_sweep.csv in {a.outdir}")


if __name__ == "__main__":
    main()
