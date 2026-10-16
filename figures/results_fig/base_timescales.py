"""
Results-section relaxation times: the response of the electrode quantities of the base cases and
sweeps to a +2 % step of the applied current density, with the axial acoustics removed.

    python base_timescales.py [-o outdir]

Runs: PlasmeX runs/base_timescales (launch.sh there; README.md has the procedure). Each steady state
of runs/base_case_* and runs/base_sweep is first re-converged with T_wall and the radial-loss
coefficient held at their values (steady/), then |J_app| is raised by 2 % and the ion current
density j_i and the heat flux q_w through the electrode surface are logged every time step.
The paper reports the records with the acoustics removed:
    cathode  e1e-4  the inertia of the heavy-gas momentum equation x 1e-4 (sound x 100; same steady
                    state; 1e-2, 1e-3 and 1e-6 agree on the base case to 3 %)
    anode    e1     physical sound speed: the acoustics play no part at the anode (base case:
                    tau 11.21 us at x 1, 11.23 us at x 1e-6), and the sweep was run at x 1 only
The 1-D axial acoustics are left out because the real arc relieves a pressure disturbance radially;
the physical-sound records (cathode e1) stay in the csv. The L-series records are not used.
The fit is runs/base_timescales/analyse.py (imported, so the numbers here are the ones in its
results.csv): t_99 (analyse's t_ss) = time after which the quantity stays within 1 % of its total
change, tau = fit of A exp(-t/tau) to the envelope of the remaining deviation between 1e-2 and 1e-4
of the change.

Writes base_timescales.{pdf,png}, cathode left, anode right:
    (a, b) normalised response of j_i and q_w of the base case
    (c, d) t_99 against |J_app|          (e, f) t_99 against T_bulk
and base_timescales.csv (one row per case and record, every record), and prints the numbers quoted
in the text (base_timescales_report.txt).
"""
import argparse
import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from base_case import HERE, C_E, C_I, rc, sci_tex, style

TS = os.path.join(os.environ.get("PLASMEX_ROOT", os.path.expanduser("~/work/PlasmeX")), "runs",
                  "base_timescales")
sys.path.insert(0, TS)
import analyse  # noqa: E402

# reference states (cathode 21 kK / 1e8 since 15 Oct 2026, anode 14 kK / 1e7 since 16 Oct) and their run directories
BASE = {"cathode": (21000.0, 1e8, 5e-4), "anode": (14000.0, 1e7, 5e-4)}
REF = {"cathode": "C_J1e8", "anode": "A14_J1e7"}
OBS = (("J_i", r"$j_i$", C_E, "o"), ("q_wall", r"$q_w$", C_I, "s"))
# records used, in order of preference; e1e-4_tt1e-4 is e1e-4 at imp_ttol 1e-4
# anode: e1_W465 = q_wall with the copper work function 4.65 eV (the decks have 2.63; W enters only q_wall there:
# the 14 kK / 1e7 state is bit-identical at 2.63 and 4.65). Only the runs/sweep_anode14 cases have this record.
REC = {"cathode": ("e1e-4", "e1e-4_tt1e-4"), "anode": ("e1_W465",)}
GREY = "#8a8a8a"


def rows():
    out = []
    for el, name, cd, plt_ in analyse.cases():
        for r in analyse.case_rows(el, name, cd, plt_):
            out.append(r)
    return out


def usable(r):
    return r["done"] or r["qinf_from"] not in ("own", "own (unfinished)")


def pick(rs, el, var):
    """The cases of one sweep (var = 'J' or 'T_arc'), the others held at the base values; per case the
    first usable record of REC[el] that has a tau."""
    Tb, Jb, Lb = BASE[el]
    held = {"L": Lb, "J": Jb, "T_arc": Tb}
    sel = []
    for case in sorted({r["case"] for r in rs if r["electrode"] == el}):
        cr = [r for r in rs if r["case"] == case and usable(r) and np.isfinite(r["J_i_tau"])
              and all(abs(r[k] / v - 1) < 1e-3 for k, v in held.items() if k != var)]
        for rec in REC[el]:
            m = [r for r in cr if r["rec"] == rec]
            if m:
                sel.append(m[0]); break
    return sorted(sel, key=lambda r: r[var])


def response(ax, el):
    cd = os.path.join(TS, el, REF[el])
    h = analyse.load(os.path.join(cd, REC[el][0]))
    t = h["t"][1:]
    for k, olab, col, _ in OBS:
        q = h[analyse.OBS[k]]
        n = max(5, int(0.05 * q.size))
        y = (q - q[0]) / (q[-n:].mean() - q[0])
        ax.plot(t, y[1:], "-", color=col, lw=1.2, label=olab)
    ax.set_xscale("log")
    ax.set_xlim(1e-12, t[-1])
    ax.axhline(1, color=GREY, lw=0.5); ax.axhline(0, color=GREY, lw=0.5)
    ax.set_xlabel(r"$t$ (s)")
    ax.set_ylabel(r"$(f-f_0)/(f_\infty-f_0)$")


def t99_panel(ax, rs, el, var, scale, xlab):
    sel = pick(rs, el, var)
    x = np.array([r[var] for r in sel]) / scale
    base = np.array([abs(r["J"] / BASE[el][1] - 1) < 1e-3 and abs(r["T_arc"] / BASE[el][0] - 1) < 1e-3
                     for r in sel])
    for k, olab, col, mk in OBS:
        y = np.array([r[f"{k}_tss"] for r in sel]) * 1e6
        ax.plot(x, y, "-", color=col, lw=1.0)
        ax.plot(x, y, mk, ms=3.8, mfc="white", mec=col, mew=0.9)
        ax.plot(x[base], y[base], mk, ms=3.8, color=col)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(xlab)
    ax.set_ylabel(r"$t_{99}$ ($\mu$s)")


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
    fig, axs = plt.subplots(3, 2, figsize=(7.0, 7.0), gridspec_kw=dict(hspace=0.5, wspace=0.28))
    tags = "abcdef"
    for c, el in enumerate(("cathode", "anode")):
        response(axs[0, c], el)
        t99_panel(axs[1, c], rs, el, "J", 1e6, r"$|J_{\rm app}|$ (10$^6$ A m$^{-2}$)")
        t99_panel(axs[2, c], rs, el, "T_arc", 1e3, r"$T_{\rm bulk}$ (10$^3$ K)")
        for r_i in range(3):
            style(axs[r_i, c], xlog=False)
            axs[r_i, c].text(0.0, 1.03, f"({tags[2 * r_i + c]})", transform=axs[r_i, c].transAxes,
                             ha="left", va="bottom", fontsize=8.5)
        Tb, Jb, _ = BASE[el]
        axs[0, c].set_title(f"{el.capitalize()}, $T_{{\\rm bulk}}={Tb:.0f}$ K, "
                            + (f"$J_{{\\rm app}}={sci_tex(Jb)}$ A m$^{{-2}}$" if el == "cathode"
                               else f"$|J_{{\\rm app}}|={sci_tex(Jb)}$ A m$^{{-2}}$"),
                            fontsize=9, pad=14)
    h = [Line2D([], [], color=col, marker=mk, ms=3.8, mfc="white", mec=col, mew=0.9, lw=1.0, label=olab)
         for _, olab, col, mk in OBS]
    fig.legend(handles=h, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.965))
    # numbers quoted in the text
    for el in ("cathode", "anode"):
        for var in ("J", "T_arc"):
            print(f"{el:7s} {var:5s} " + "  ".join(
                f"{r[var]:.3g}: t99 {r['J_i_tss'] * 1e6:.2f}/{r['q_wall_tss'] * 1e6:.2f} us,"
                f" tau {r['J_i_tau'] * 1e6:.3f}/{r['q_wall_tau'] * 1e6:.3f} us ({r['rec']})"
                for r in pick(rs, el, var)))
    for el in ("cathode", "anode"):
        cd = os.path.join(TS, el, REF[el])
        h = analyse.load(os.path.join(cd, REC[el][0]))
        for k, *_ in OBS:
            q = h[analyse.OBS[k]]
            n = max(5, int(0.05 * q.size))
            y = (q - q[0]) / (q[-n:].mean() - q[0])
            print(f"{el:7s} base {k:6s} response " + " ".join(
                f"{tt:.0e}:{np.interp(tt, h['t'], y):.3f}" for tt in (1e-11, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 5e-5)))
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(a.outdir, f"base_timescales.{ext}"), bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"wrote base_timescales (.pdf/.png) and base_timescales.csv in {a.outdir}")


if __name__ == "__main__":
    main()
