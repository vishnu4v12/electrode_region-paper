"""
Results-section sweep figures: the electrode quantities of the base cases as the bulk-plasma state
(T_bulk at fixed J_app) and the current density (J_app at fixed T_bulk) are varied.

    python base_sweep.py [-o outdir]

Runs (PlasmeX runs/base_sweep/<electrode>/<case>, launched by runs/base_sweep/launch.sh, plus the
two base cases): the last plotfile of each run directory is read. Writes, each as .pdf and .png:
    base_sweep_T   quantities against T_bulk (anode at 1e7 A/m^2, cathode at 1e8 A/m^2)
    base_sweep_J   quantities against |J_app| (anode at 14 000 K, cathode at 21 000 K)
cathode in the left column, anode in the right (the paper's two-case convention), rows:
    q_w                 net heat flux into the electrode, eq (24) (plasmex_info q_wall.total)
    j_i/j               ion fraction of the current at the surface, + = ions into the electrode
    phi_s, phi(L)       potential at the sheath edge (|n_i - n_e|/n_i = 1 %) and at x = L,
                        relative to the electrode
    x_s                 sheath width: the distance of the sheath edge from the surface
                        (base_case.sheath_edge, the interpolated 1 % crossing)
and base_sweep.csv with these, T_e at the surface, the cathode T_wall, the emission fraction j_em/j,
the surface field and densities.
The base case of each sweep is drawn as a filled marker.

Every number comes from base_case.load (bohm's leaf-cell sampler and the plotfile's plasmex_info),
the same route as table tab:base. The steady check of each run is printed, as in base_case.py.
"""
import argparse
import csv
import glob
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from base_case import HERE, K_B, E_CH, C_E, C_X, load, rc, sci_tex, sheath_edge, style

RUNS = os.path.join(os.environ.get("PLASMEX_ROOT", os.path.expanduser("~/work/PlasmeX")), "runs")
# Reference states: cathode 21 000 K / 1e8 (user, 15 Oct 2026), anode 14 000 K / 1e7 (user, 16 Oct 2026).
# Since 16 Oct 2026 (user) the column radius is ASSUMED, R = 0.5 mm at the cathode and 3.5 mm at the anode, in
# every case: runs/fixedR (run_case.sh, ladder.sh there), each restarted from the calibrated-R steady state of the
# same conditions. Cathode J >= 1.2e8 by a ladder in Lambda from the calibrated value (<case>/ladder/).
# Before: anode runs/sweep_anode14, cathode T runs/sweep_ref, cathode J runs/base_sweep (calibrated R).
BASE = {"Cathode": (21000.0, 1e8, "fixedR/cathode_T21000_J1e8"), "Anode": (14000.0, 1e7, "fixedR/anode_T14000_J1e7")}
CASES = {   # electrode -> [(T_bulk [K], |J_app| [A/m^2], run directory under RUNS)]
    "Anode":   [(T, 1e7, f"fixedR/anode_T{T}_J1e7") for T in (11000, 12000, 13000)]
             + [(14000.0, float(J), f"fixedR/anode_T14000_J{J}") for J in ("7e6", "8e6", "9e6", "2e7")],
    "Cathode": [(T, 1e8, f"fixedR/cathode_T{T}_J1e8") for T in (19000, 20000, 22000, 23000)]
             + [(21000.0, float(J), f"fixedR/cathode_T21000_J{J}")
                for J in ("8e7", "9e7", "1.2e8", "1.248e+08", "1.29792e+08", "1.34984e+08", "1.40383e+08",
                          "1.45998e+08", "1.48e+08", "1.5e+08")],
}


def last_plotfile(d):
    """The plotfile the steady-state trigger wrote (run.log's "Stopping: steady state"), or, for a
    base case (no such line when it ran to stop_time), the last one. None if neither. A case's
    tight/ (the restart at plasmex.steady_tol = 1e-5, which closes the cathode heat balance) wins."""
    if os.path.isdir(os.path.join(d, "tight")):
        pf = last_plotfile(os.path.join(d, "tight"))
        if pf is not None:
            return pf
    pfs = sorted(glob.glob(os.path.join(d, "plt[0-9]*")))
    log = os.path.join(d, "run.log")
    txt = open(log).read() if os.path.isfile(log) else ""
    if "Stopping: steady state" in txt:
        return os.path.join(d, txt.split("Stopping: steady state")[1].split("wrote ")[1].split(" ")[0])
    if os.path.basename(d).startswith("base_case_") and pfs:
        return pfs[-1]
    return None


def quantities(s):
    """The plotted electrode quantities of one steady state (base_case.load's dict)."""
    info = s["info"]
    g = lambda k: float(info.get(k, np.nan))
    d, ne, ni, phi = s["d"], s["ne"], s["ni"], s["phi"]
    x_s, phi_s, _ = sheath_edge(d, ne, ni, phi)        # as table tab:base
    bd = s["budget"]
    resid = np.nanmax(np.abs(bd["DEDT"][1:-1]))
    big = max(np.nanmax(np.abs(bd[k])) for k in ("JOULE", "ELASTIC", "ION"))
    j = abs(np.median(s["j"]))
    return dict(
        q_w=g("q_wall.total") / 1e6, ji_j=s["ji_j"][0],
        phi_s=phi_s, x_s=x_s, phi_L=phi[-1], Te_w=s["Te"][0], T_wall=g("prob.wall_temp"),
        jem_j=E_CH * g("q_wall.e_flux_emitted") / j, E_w=g("q_wall.E_face"), ne_w=ne[0], ni_w=ni[0],
        steady=resid / big, j_spread=np.ptp(s["j"]) / j, t=s["t"], plotfile=s["name"])


def collect():
    rows = []
    for el, cases in CASES.items():
        Tb, Jb, bdir = BASE[el]
        for T, J, sub in [(Tb, Jb, bdir)] + cases:
            pf = last_plotfile(os.path.join(RUNS, sub))
            if pf is None:
                print(f"  {el} T={T:.0f} J={J:.1e}: no plotfile in {sub}, skipped")
                continue
            q = quantities(load(pf))
            q.update(electrode=el, T_bulk=T, J_app=J, run=sub, base=(sub == bdir))
            flag = "" if q["steady"] < 2e-3 else "   <-- NOT STEADY?"
            print(f"  {el:7s} T={T:7.0f} J={J:.2e}  {sub}/{q['plotfile']}  t={q['t']:.3e}"
                  f"  steady {q['steady']:.1e}  j spread {q['j_spread']:.1e}{flag}")
            rows.append(q)
    return rows


ROWS = [   # key(s), y label
    (("q_w",),             r"$q_w$ (MW m$^{-2}$)"),
    (("ji_j",),            r"$j_i/j$ at the surface"),
    (("phi_s", "phi_L"),   r"Potential (V)"),
    (("x_s",),             r"Sheath width $x_{\rm s}$ ($\mu$m)"),
]
KSCALE = {"x_s": 1e-6}
LEG = {"phi_s": r"$\phi_{\rm s}$", "phi_L": r"$\phi(L)$"}
COL = {"phi_s": C_X, "phi_L": C_E}


def fig_sweep(rows, var, stem):
    rc()
    fig, axs = plt.subplots(len(ROWS), 2, figsize=(7.0, 7.1), sharex="col",
                            gridspec_kw=dict(hspace=0.32, wspace=0.28))
    tags = "abcdefghijklmnop"
    for c, el in enumerate(("Cathode", "Anode")):
        Tb, Jb, _ = BASE[el]
        sel = [r for r in rows if r["electrode"] == el
               and (r["J_app"] == Jb if var == "T_bulk" else r["T_bulk"] == Tb)]
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
            # A nearly constant quantity (the cathode against T_bulk) is drawn on a span of at least 10 % of its
            # value, so that autoscaling does not turn a 0.03 % change into a full-height trend.
            ys = np.concatenate([np.array([r[k] for r in sel]) / KSCALE.get(k, 1.0) for k in keys])
            span, mid = 0.1 * np.max(np.abs(ys)), 0.5 * (ys.max() + ys.min())
            if ys.max() - ys.min() < span:
                ax.set_ylim(mid - 0.6 * span, mid + 0.6 * span)
            if len(keys) > 1:
                ax.legend(loc="best", frameon=False, handlelength=1.6)
            if r_i == 1:
                ax.axhline(0, color="#b0b0b0", lw=0.5)
            ax.set_ylabel(ylab)
            style(ax, xlog=False)
            ax.text(0.0, 1.03, f"({tags[2 * r_i + c]})", transform=ax.transAxes, ha="left",
                    va="bottom", fontsize=8.5)
        if var == "T_bulk":
            axs[0, c].set_title(f"{el}, $|J_{{\\rm app}}|={sci_tex(Jb)}$ A m$^{{-2}}$"
                                if el == "Anode" else
                                f"{el}, $J_{{\\rm app}}={sci_tex(Jb)}$ A m$^{{-2}}$", fontsize=9, pad=14)
            axs[-1, c].set_xlabel(r"$T_{\rm bulk}$ (10$^3$ K)")
        else:
            axs[0, c].set_title(f"{el}, $T_{{\\rm bulk}}={Tb:.0f}$ K", fontsize=9, pad=14)
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
    keys = ["electrode", "T_bulk", "J_app", "q_w", "ji_j", "phi_s", "x_s", "phi_L", "Te_w",
            "T_wall", "jem_j", "E_w", "ne_w", "ni_w", "steady", "j_spread", "t", "run", "plotfile"]
    with open(os.path.join(a.outdir, "base_sweep.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for r in sorted(rows, key=lambda r: (r["electrode"], r["J_app"], r["T_bulk"])):
            w.writerow([f"{r[k]:.5g}" if isinstance(r[k], float) else r[k] for k in keys])
    print("\n" + " ".join(f"{k:>9s}" for k in keys[:13]))
    for r in sorted(rows, key=lambda r: (r["electrode"], r["J_app"], r["T_bulk"])):
        print(" ".join(f"{r[k]:9.4g}" if isinstance(r[k], float) else f"{str(r[k]):>9s}" for k in keys[:13]))
    fig_sweep(rows, "T_bulk", os.path.join(a.outdir, "base_sweep_T"))
    fig_sweep(rows, "J_app", os.path.join(a.outdir, "base_sweep_J"))
    print(f"\nwrote base_sweep_T, base_sweep_J (.pdf/.png) and base_sweep.csv in {a.outdir}")


if __name__ == "__main__":
    main()
