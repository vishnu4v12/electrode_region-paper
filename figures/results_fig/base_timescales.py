"""
Results-section relaxation times: the response of the electrode quantities to a +2 % step of the
applied current density, and (18 Oct 2026, user) to a +2 % step of the bulk temperature at fixed current,
for the cathode case of table tab:base and for the anode case and the anode sweep (sec:res-anode-sweep).

    python base_timescales.py [-o outdir]

Runs: PlasmeX runs/Tedge1000R_timescales (17 Oct 2026; launch.sh and pipeline.sh there), the
runs/base_timescales procedure on the runs/Tedge1000 states, anode at R = 3.5 mm (anode_T15000_J4e6_L0.5mm_R3.5mm
and anode_sweep_R3.5mm; its cathode/ links to the record in runs/Tedge1000_timescales, the no-loss anode set). Each steady state is first re-converged with
T_wall held at its applied value (steady/), then |J_app| is raised by 2 % and the ion current density
j_i and the heat flux q_w through the electrode surface are logged every time step (hist.txt). Records:
    cathode  e1e-4  the inertia of the heavy-gas momentum equation x 1e-4 (axial acoustics removed;
                    same steady state)
    anode    e1e-4  the same, since 18 Oct 2026 (user: acoustics off in every case). The physical-sound records
                    e1 / e1_T2 give D within 0.1 % in every anode case.
Measures, per quantity f (j_i, q_w):
    response  (f/f_0 - 1)/0.02: the change in units of the applied relative step
    D         the quasi-steady error coefficient, D = int_0^inf (g_inf - g) dt with g = (f/f_0 - 1)/0.02
              (18 Oct 2026, replaces t_s in the paper). For a linear layer and a current that varies slowly
              compared with tau, f - f_qs = -D f_0 dlnJ/dt (first moment of the step response), so the layer is
              quasi-steady to a relative accuracy eps when J/|dJ/dt| >= |D|/eps.
    t_s       (csv and report only) the time after which f stays within 1 % of the applied step of its final value, i.e.
              |f - f_inf| < 2e-4 |f_0| (f_inf = mean over the last 5 % of the record's duration). The 1 % of a quantity's
              OWN change (t_99 before 17 Oct) is ill-defined here: the net change of the anode j_i is
              -0.5 .. +0.7 % and changes sign within the current sweep.
    tau       the slow exponential mode, f = f_inf + A exp(-t/tau) least-squares fitted to the signed
              record over the window FIT (anode 5-150 us, cathode 0.3-4 us); left out where |A| is
              below 1 % of the step (t_s's tolerance) or 10x the fit residual (anode 13 kK: |A| = 2e-5).

Temperature step (18 Oct 2026): record <REC>_T2, arc_Te x 1.02 at the same arc_Japp (launch.sh KICK=1 KICKT=1.02);
the same measures with g = (f/f_0 - 1)/(dT/T), columns suffixed _T in the csv. The table cases also have
<REC>_Tm2 (arc_Te x 0.98) and <REC>_Jm2 (arc_Japp x 0.98), the linearity checks; the report prints both
signs and their mean.

Writes base_timescales.{pdf,png} (18 Oct 2026, user: two panels, no tau, the cathode in the text only):
    (a) anode t_1% = |D|/0.01 against |J_app| at 15 kK, (b) against T_bulk at 4e6; solid: current step,
    dashed: temperature step
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
import matplotlib.ticker
from matplotlib.lines import Line2D

from base_case import HERE, C_E, C_I, rc, sci_tex, style

TS = os.path.join(os.environ.get("PLASMEX_ROOT", os.path.expanduser("~/work/PlasmeX")), "runs",
                  "I200_timescales")   # 23 Oct 2026: R for I = 200 A (was wannier_rel_timescales, R 0.5 / 3.5 mm)
sys.path.insert(0, TS)
import analyse  # noqa: E402

STEP = 0.02
BASE = {"cathode": (21000.0, 1e8), "anode": (15000.0, 4e6)}
REF = {"cathode": "C_J1e8", "anode": "A15_J4e6"}
REC = {"cathode": "e1e-4", "anode": "e1e-4"}   # 18 Oct 2026 (user): acoustics off at both electrodes (anode was e1)
FIT = {"cathode": (3e-7, 4e-6), "anode": (5e-6, 1.5e-4)}
OBS = (("J_i", "JI_wallface", r"$j_i$", C_E, "o"), ("q_wall", "q_wall", r"$q_w$", C_I, "s"))
GREY = "#8a8a8a"
# 18 Oct 2026 (user): the 2e6 A/m^2 anode case is no longer in the paper (sec:res-anode-sweep, sec:res-time)
EXCLUDE = {"A15_J2e6"}
# the +-2 % linearity pair, where it is not REC: the anode -2 % T step aborts at heavy_inertia 1e-4 whatever imp_ttol
# or min_time_step (failed_e1e-4_Tm2_* there), so that pair is run at 1e-3 (D_T of the +2 % step within 0.03 % of 1e-4)
LIN_REC = {("anode", "T"): "e1e-3"}


def final(t, q):
    """The final value: the mean over the last 5 % of the record's DURATION (not of its samples, which
    crowd into the fast transient when dt is small)."""
    return q[t >= 0.95 * t[-1]].mean()


def settle(t, q):
    qe = final(t, q)
    y = np.abs(q - qe) / abs(q[0])
    env = np.maximum.accumulate(y[::-1])[::-1]
    k = np.flatnonzero(env < 0.01 * STEP)
    return (t[k[0]] if k.size else np.nan), qe


# the step kinds: (key, column suffix, record suffix, signed relative step)
STEPS = (("J", "", "", STEP), ("T", "_T", "_T2", STEP))


def lag(t, q, step=STEP):
    """Quasi-steady error coefficient D [s]: the area between the final value and the step response g."""
    g = (q / q[0] - 1) / step
    return np.trapezoid(final(t, g) - g, t)


def slow_mode(t, q, win):
    m = (t > win[0]) & (t < win[1])
    if m.sum() < 5:   # a record still running
        return np.nan, np.nan
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
            if os.path.basename(cd) in EXCLUDE:
                continue
            if not os.path.exists(os.path.join(cd, REC[el], "hist.txt")):
                continue
            how = open(os.path.join(cd, "steady", "HOW_STARTED.txt")).read()
            T = float(re.search(r"prob\.arc_Te=(\S+)", how).group(1))
            J = abs(float(re.search(r"prob\.arc_Japp=(\S+)", how).group(1)))
            r = dict(electrode=el, case=os.path.basename(cd), T_arc=T, J=J, rec=REC[el])
            for _, sfx, rsfx, step in STEPS:
                rec = os.path.join(cd, REC[el] + rsfx)
                if not os.path.exists(os.path.join(rec, "hist.txt")):
                    continue
                h = analyse.load(rec)
                r[f"done{sfx}"] = "finalized" in open(os.path.join(rec, "run.log")).read()
                r[f"t_end{sfx}"] = h["t"][-1]
                for k, col, *_ in OBS:
                    q = h[col]
                    ts, qe = settle(h["t"], q)
                    tau, A = slow_mode(h["t"], q, FIT[el])
                    r.update({f"{k}_net{sfx}": qe / q[0] - 1, f"{k}_D{sfx}": lag(h["t"], q, step),
                              f"{k}_ts{sfx}": ts, f"{k}_tau{sfx}": tau, f"{k}_A{sfx}": A})
            out.append(r)
    return out


def sweep(rs, var):
    Tb, Jb = BASE["anode"]
    held = {"J": Jb, "T_arc": Tb}
    sel = [r for r in rs if r["electrode"] == "anode"
           and all(abs(r[k] / v - 1) < 1e-3 for k, v in held.items() if k != var)]
    return sorted(sel, key=lambda r: r[var])


EPS = 0.01   # the accuracy of the quasi-steady approximation: t_J >= |D|/EPS


def ts_panel(ax, rs, var, scale, xlab):
    sel = sweep(rs, var)
    x = np.array([r[var] for r in sel]) / scale
    base = np.array([r["case"] == REF["anode"] for r in sel])
    for (_, sfx, _, _), ls in zip(STEPS, ("-", "--")):
        for k, _, olab, c, mk in OBS:
            y = np.abs([r.get(f"{k}_D{sfx}", np.nan) for r in sel]) / EPS * 1e6
            ax.plot(x, y, ls, color=c, lw=1.0)
            ax.plot(x, y, mk, ms=3.8, mfc="white", mec=c, mew=0.9)
            ax.plot(x[base], y[base], mk, ms=3.8, color=c)
    ax.set_yscale("log")
    ax.set_ylim(0.3, 3e4)
    ax.set_xlabel(xlab)
    ax.set_ylabel(r"$t_{1\%}$ ($\mu$s)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("-o", "--outdir", default=HERE)
    a = ap.parse_args(argv)
    rs = rows()
    keys = list(dict.fromkeys(k for r in rs for k in r))
    with open(os.path.join(a.outdir, "base_timescales.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(keys)
        for r in rs:
            w.writerow([f"{r[k]:.5g}" if isinstance(r.get(k), float) else r.get(k, "") for k in keys])
    rc()
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.6), gridspec_kw=dict(wspace=0.3))
    ts_panel(axs[0], rs, "J", 1e6, r"$J_{\rm app}$ (10$^6$ A m$^{-2}$)")
    ts_panel(axs[1], rs, "T_arc", 1e3, r"$T_{\rm bulk}$ (10$^3$ K)")
    axs[0].set_title(r"Anode, $T_{\rm bulk}=15\,000$ K", fontsize=9, pad=14)
    axs[1].set_title(r"Anode, $J_{\rm app}=4\times10^{6}$ A m$^{-2}$", fontsize=9, pad=14)
    for i, ax in enumerate(axs.flat):
        style(ax, xlog=False)
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        ax.text(0.0, 1.03, f"({'ab'[i]})", transform=ax.transAxes, ha="left", va="bottom", fontsize=8.5)
    h = [Line2D([], [], color=c, marker=mk, ms=3.8, mfc="white", mec=c, mew=0.9, lw=1.0, label=olab)
         for _, _, olab, c, mk in OBS]
    h += [Line2D([], [], color="k", ls=ls, lw=1.0, label=lab)
          for ls, lab in (("-", r"step of $J_{\rm app}$"), ("--", r"step of $T_{\rm bulk}$"))]
    fig.legend(handles=h, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.08))
    # numbers quoted in the text
    for sk, sfx, _, _ in STEPS:
        for r in rs:
            if f"t_end{sfx}" not in r:
                print(f"{r['electrode']:7s} {r['case']:12s} {sk}-step: no record")
                continue
            print(f"{r['electrode']:7s} {r['case']:12s} {sk}-step T {r['T_arc']:.0f} J {r['J']:.3g} "
                  f"{'' if r['done' + sfx] else 'PARTIAL '}t_end {r['t_end' + sfx] * 1e6:.1f} us | " + "  ".join(
                      f"{k}: net {r[k + '_net' + sfx] * 100:+.3f}% g_inf {r[k + '_net' + sfx] / STEP:+.3f} "
                      f"D {r[k + '_D' + sfx] * 1e6:+.3f} us D/g_inf {r[k + '_D' + sfx] / r[k + '_net' + sfx] * STEP * 1e6:.2f} us "
                      f"t_s {r[k + '_ts' + sfx] * 1e6:.2f} us tau {r[k + '_tau' + sfx] * 1e6:.3f} us "
                      f"t_1% {abs(r[k + '_D' + sfx]) / EPS * 1e6:.1f} us" for k, *_ in OBS))
    # linearity: +2 % against -2 % steps of the current and of the temperature, for the table cases
    for sk, (rp, rm) in (("J", ("", "_Jm2")), ("T", ("_T2", "_Tm2"))):
        for el in ("cathode", "anode"):
            for k, col, *_ in OBS:
                d = {}
                for rs_, step in ((rp, STEP), (rm, -STEP)):
                    rec = os.path.join(TS, el, REF[el], LIN_REC.get((el, sk), REC[el]) + rs_)
                    if os.path.exists(os.path.join(rec, "run.log")) and "finalized" in open(os.path.join(rec, "run.log")).read():
                        h = analyse.load(rec)
                        q = h[col]
                        d[rs_] = ((final(h["t"], q) / q[0] - 1) / step, lag(h["t"], q, step))
                if len(d) == 2:
                    (gp, dp), (gm, dm) = d[rp], d[rm]
                    print(f"{el:7s} base {k:6s} {sk}-step linearity: +2% g_inf {gp:+.3f} D {dp * 1e6:+.3f} us | "
                          f"-2% g_inf {gm:+.3f} D {dm * 1e6:+.3f} us | mean D {(dp + dm) / 2 * 1e6:+.3f} us, "
                          f"half-difference {abs(dp - dm) / abs(dp + dm) * 100:.1f} %")
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
