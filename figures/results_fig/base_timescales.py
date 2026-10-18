"""
Results-section relaxation times: the response of the electrode quantities to a +2 % step of the
applied current density, for the cathode case of table tab:base and for the anode case and the anode
sweep (sec:res-anode-sweep).

    python base_timescales.py [-o outdir]

Runs: PlasmeX runs/Tedge1000R_timescales (17 Oct 2026; launch.sh and pipeline.sh there), the
runs/base_timescales procedure on the runs/Tedge1000 states, anode at R = 3.5 mm (anode_T15000_J4e6_L0.5mm_R3.5mm
and anode_sweep_R3.5mm; its cathode/ links to the record in runs/Tedge1000_timescales, the no-loss anode set). Each steady state is first re-converged with
T_wall held at its applied value (steady/), then |J_app| is raised by 2 % and the ion current density
j_i and the heat flux q_w through the electrode surface are logged every time step (hist.txt). Records:
    cathode  e1e-4  the inertia of the heavy-gas momentum equation x 1e-4 (axial acoustics removed;
                    same steady state)
    anode    e1     physical sound speed (the acoustics play no part at the anode)
Measures, per quantity f (j_i, q_w):
    response  (f/f_0 - 1)/0.02: the change in units of the applied relative step
    t_s       the time after which f stays within 1 % of the applied step of its final value, i.e.
              |f - f_inf| < 2e-4 |f_0| (f_inf = mean of the record's last 5 %). The 1 % of a quantity's
              OWN change (t_99 before 17 Oct) is ill-defined here: the net change of the anode j_i is
              -0.5 .. +0.7 % and changes sign within the current sweep.
    tau       the slow exponential mode, f = f_inf + A exp(-t/tau) least-squares fitted to the signed
              record over the window FIT (anode 5-150 us, cathode 0.3-4 us); left out where |A| is
              below 1 % of the step (t_s's tolerance) or 10x the fit residual (anode 13 kK: |A| = 2e-5).

Writes base_timescales.{pdf,png}:
    (a) cathode response, (b) anode response, (c) anode t_s and tau against |J_app| at 15 kK,
    (d) against T_bulk at 4e6
and base_timescales.csv (one row per case), and prints the numbers quoted in the text
(base_timescales_report.txt).
"""
import argparse
import csv
import glob
import os
import re
import sys

import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from base_case import HERE, C_E, C_I, rc, sci_tex, style

TS = os.path.join(os.environ.get("PLASMEX_ROOT", os.path.expanduser("~/work/PlasmeX")), "runs",
                  "Tedge1000R_timescales")
sys.path.insert(0, TS)
import analyse  # noqa: E402

STEP = 0.02
BASE = {"cathode": (21000.0, 1e8), "anode": (15000.0, 4e6)}
REF = {"cathode": "C_J1e8", "anode": "A15_J4e6"}
REC = {"cathode": "e1e-4", "anode": "e1"}
FIT = {"cathode": (3e-7, 4e-6), "anode": (5e-6, 1.5e-4)}
OBS = (("J_i", "JI_wallface", r"$j_i$", C_E, "o"), ("q_wall", "q_wall", r"$q_w$", C_I, "s"))
GREY = "#8a8a8a"


def settle(t, q):
    n = max(5, int(0.05 * q.size))
    qe = q[-n:].mean()
    y = np.abs(q - qe) / abs(q[0])
    env = np.maximum.accumulate(y[::-1])[::-1]
    k = np.flatnonzero(env < 0.01 * STEP)
    return (t[k[0]] if k.size else np.nan), qe


def slow_mode(t, q, win):
    m = (t > win[0]) & (t < win[1])
    y = q / q[0] - 1
    f = lambda tt, yi, A, tau: yi + A * np.exp(-tt / tau)  # noqa: E731
    p, _ = curve_fit(f, t[m], y[m], p0=(y[-1], y[m][0] - y[-1], 0.3 * win[1]), maxfev=20000)
    res = np.std(y[m] - f(t[m], *p))
    # a slow mode smaller than t_s's tolerance (1 % of the step) is not reported
    return (p[2] if abs(p[1]) > max(10 * res, 0.01 * STEP) else np.nan), p[1]


def rows():
    out = []
    for el in ("cathode", "anode"):
        for cd in sorted(glob.glob(os.path.join(TS, el, "*"))):
            rec = os.path.join(cd, REC[el])
            if not os.path.exists(os.path.join(rec, "hist.txt")):
                continue
            how = open(os.path.join(cd, "steady", "HOW_STARTED.txt")).read()
            T = float(re.search(r"prob\.arc_Te=(\S+)", how).group(1))
            J = abs(float(re.search(r"prob\.arc_Japp=(\S+)", how).group(1)))
            done = "finalized" in open(os.path.join(rec, "run.log")).read()
            h = analyse.load(rec)
            r = dict(electrode=el, case=os.path.basename(cd), T_arc=T, J=J, rec=REC[el], done=done, t_end=h["t"][-1])
            for k, col, *_ in OBS:
                q = h[col]
                ts, qe = settle(h["t"], q)
                tau, A = slow_mode(h["t"], q, FIT[el])
                r.update({f"{k}_net": qe / q[0] - 1, f"{k}_ts": ts, f"{k}_tau": tau, f"{k}_A": A})
            out.append(r)
    return out


def response(ax, el):
    h = analyse.load(os.path.join(TS, el, REF[el], REC[el]))
    t = h["t"][1:]
    for _, col, olab, c, _ in OBS:
        q = h[col]
        ax.plot(t, (q[1:] / q[0] - 1) / STEP, "-", color=c, lw=1.2, label=olab)
    ax.set_xscale("log")
    ax.set_xlim(1e-12, t[-1])
    ax.axhline(0, color=GREY, lw=0.5)
    ax.set_xlabel(r"$t$ (s)")
    ax.set_ylabel(r"$(f/f_0-1)\,/\,(\Delta J/J)$")


def sweep(rs, var):
    Tb, Jb = BASE["anode"]
    held = {"J": Jb, "T_arc": Tb}
    sel = [r for r in rs if r["electrode"] == "anode"
           and all(abs(r[k] / v - 1) < 1e-3 for k, v in held.items() if k != var)]
    return sorted(sel, key=lambda r: r[var])


def ts_panel(ax, rs, var, scale, xlab):
    sel = sweep(rs, var)
    x = np.array([r[var] for r in sel]) / scale
    base = np.array([r["case"] == REF["anode"] for r in sel])
    for k, _, olab, c, mk in OBS:
        y = np.array([r[f"{k}_ts"] for r in sel]) * 1e6
        ax.plot(x, y, "-", color=c, lw=1.0)
        ax.plot(x, y, mk, ms=3.8, mfc="white", mec=c, mew=0.9)
        ax.plot(x[base], y[base], mk, ms=3.8, color=c)
    tau = np.array([r["J_i_tau"] for r in sel]) * 1e6
    ax.plot(x, tau, "--", color="k", lw=0.9)
    ax.plot(x, tau, "D", ms=3.0, mfc="white", mec="k", mew=0.8)
    ax.plot(x[base], tau[base], "D", ms=3.0, color="k")
    ax.set_ylim(bottom=0)
    ax.set_xlabel(xlab)
    ax.set_ylabel(r"$t_{\rm s}$, $\tau$ ($\mu$s)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("-o", "--outdir", default=HERE)
    a = ap.parse_args(argv)
    rs = rows()
    keys = list(rs[0])
    with open(os.path.join(a.outdir, "base_timescales.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(keys)
        for r in rs:
            w.writerow([f"{r[k]:.5g}" if isinstance(r[k], float) else r[k] for k in keys])
    rc()
    fig, axs = plt.subplots(2, 2, figsize=(7.0, 4.9), gridspec_kw=dict(hspace=0.55, wspace=0.3))
    for c, el in enumerate(("cathode", "anode")):
        response(axs[0, c], el)
        Tb, Jb = BASE[el]
        axs[0, c].set_title(f"{el.capitalize()}, $T_{{\\rm bulk}}={Tb:.0f}$ K, "
                            + (f"$J_{{\\rm app}}={sci_tex(Jb)}$ A m$^{{-2}}$" if el == "cathode"
                               else f"$|J_{{\\rm app}}|={sci_tex(Jb)}$ A m$^{{-2}}$"),
                            fontsize=9, pad=14)
    ts_panel(axs[1, 0], rs, "J", 1e6, r"$|J_{\rm app}|$ (10$^6$ A m$^{-2}$)")
    ts_panel(axs[1, 1], rs, "T_arc", 1e3, r"$T_{\rm bulk}$ (10$^3$ K)")
    axs[1, 0].set_title(r"Anode, $T_{\rm bulk}=15\,000$ K", fontsize=9, pad=14)
    axs[1, 1].set_title(r"Anode, $|J_{\rm app}|=4\times10^{6}$ A m$^{-2}$", fontsize=9, pad=14)
    for i, ax in enumerate(axs.flat):
        style(ax, xlog=False)
        ax.text(0.0, 1.03, f"({'abcd'[i]})", transform=ax.transAxes, ha="left", va="bottom", fontsize=8.5)
    h = [Line2D([], [], color=c, marker=mk, ms=3.8, mfc="white", mec=c, mew=0.9, lw=1.0, label=olab)
         for _, _, olab, c, mk in OBS]
    h.append(Line2D([], [], color="k", ls="--", marker="D", ms=3.0, mfc="white", mew=0.8, lw=0.9,
                    label=r"$\tau$ (slow mode)"))
    fig.legend(handles=h, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.02))
    # numbers quoted in the text
    for r in rs:
        print(f"{r['electrode']:7s} {r['case']:12s} T {r['T_arc']:.0f} J {r['J']:.3g} {'' if r['done'] else 'PARTIAL '}"
              f"t_end {r['t_end'] * 1e6:.1f} us | " + "  ".join(
                  f"{k}: net {r[k + '_net'] * 100:+.3f}% t_s {r[k + '_ts'] * 1e6:.2f} us tau {r[k + '_tau'] * 1e6:.3f} us "
                  f"A {r[k + '_A'] * 100:+.4f}%" for k, *_ in OBS))
    for el in ("cathode", "anode"):
        h = analyse.load(os.path.join(TS, el, REF[el], REC[el]))
        for k, col, *_ in OBS:
            q = h[col]
            print(f"{el:7s} base {k:6s} (f/f0-1)/step " + " ".join(
                f"{tt:.0e}:{np.interp(tt, h['t'], (q / q[0] - 1) / STEP):+.3f}"
                for tt in (1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4)))
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(a.outdir, f"base_timescales.{ext}"), bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"wrote base_timescales (.pdf/.png) and base_timescales.csv in {a.outdir}")


if __name__ == "__main__":
    main()
