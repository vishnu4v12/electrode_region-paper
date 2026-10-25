"""
PlasmeX 1-D cathode plotfile(s) vs Almeida, Benilov & Naidis (2008) J. Phys. D 41 245201 at
jc = 1e7 A/m^2 (Ar, p0 = 1 bar, Tw = 3500 K): Fig. 5(b) and curve b of Fig. 6, in one figure.

    python compare.py plt19500 [plt...] [-o outdir]

One figure per plotfile, compareJ1e7_<plotfile>.{pdf,png}, three stacked panels against the
distance from the cathode:
    (a) n_e, n_i            Fig. 5(b)
    (b) T_e, T_h            Fig. 5(b)
    (c) j_i/j (left axis)   Fig. 5(b)   and   -E (right axis)   Fig. 6 curve b
A terminal table also lists the potential phi (Fig. 6 b), which is not plotted.

ORIENTATION AND SIGNS. The paper's x is the distance from the cathode, cathode on the LEFT. The
PlasmeX cathode is the x_hi face, so the simulation is plotted against d = x_hi - x. The
paper's E is along its x, i.e. E_paper = -EFLDX, negative in the sheath (the field points to
the cathode). Panel (c) plots -E_paper = EFLDX, the field magnitude towards the cathode, so
that it and j_i/j are both positive and share one zero line. PHI is zero on the cathode face,
as the paper's phi is.

The paper curves are read from paper_data/fig5b_extracted.csv and fig6b_extracted.csv
(written by paper_data/extract.py). The simulation is read off the AMR leaf cells by
bohm.sample, loaded from the live PlasmeX tree (../plasmex.py); j_i/j = JIX / |JEX + JIX|, the ion
current density INTO the cathode over |J_app| (the manuscript's convention since 21 Oct 2026). The
boundary cell holds the FACE value of JEX/JIX, so the wall j_i/j is the flux the solver applied.
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, NullFormatter, MultipleLocator

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from plasmex import exe_module  # noqa: E402
from electrode_mark import mark_electrode  # noqa: E402

DATA = os.path.join(HERE, "paper_data")
XLIM = (1e-7, 1e-2)                          # Almeida et al. Figs 5, 6 start at 1e-7 m

# Colour = particle population (electrons / heavy particles); line style = source.
C_E, C_H, C_J, C_F = "#2a78d6", "#eb6834", "#2b2b2b", "#1a9e77"
SIM = dict(lw=1.5, ls="-")
REF = dict(lw=1.3, ls=(0, (4.5, 2.2)))

# The paper's run conditions; a plotfile that departs from them is flagged.
JC, P0, T_WALL = 1e7, 1.0e5, 3500.0
TABLE1_JIW = 0.216                           # Almeida et al. Table 1, wall j_i/j for 5(b)

N_YLIM = (1e19, 1e24)
T_YLIM, T_MAJOR = (0, 21), 5                 # kK
JI_YLIM, JI_MAJOR = (-0.05, 0.35), 0.1
E_UNIT, E_TOP, E_MAJOR = 1e7, 2.5, 1.0       # -E axis; its zero is j_i/j's zero
LABELS = dict(ne=(1.45e-7, 1.4e19), ni=(1.45e-7, 4.5e21), Te=(1.45e-7, 19.3), Th=(1.45e-7, 4.3),
              ji=(2.0e-6, 0.226), E=(9.0e-7, 0.7))
SIM_LABEL, REF_LABEL = "Present model", "Almeida et al (2008)"


def load_csv(name):
    """{quantity: (x_m, value)} from one of the paper_data CSVs."""
    path = os.path.join(DATA, name)
    if not os.path.isfile(path):
        raise SystemExit(f"missing {path} -- run paper_data/extract.py")
    cols = {}
    with open(path) as f:
        for line in f:
            if line.startswith(("#", "quantity")):
                continue
            q, _, x, v = line.strip().split(",")
            cols.setdefault(q, ([], []))
            cols[q][0].append(float(x))
            cols[q][1].append(float(v))
    return {q: (np.array(x), np.array(v)) for q, (x, v) in cols.items()}


def load_sim(plt_file):
    b = exe_module("bohm")
    d = b.sample(plt_file, ("NE", "ARP", "ETEMP", "TGAS", "JEX", "JIX", "PHI", "EFLDX"),
                 want_derived=False)
    dist = d["hi"] - d["x"]
    o = np.argsort(dist)
    j = d["JEX"] + d["JIX"]
    return {"d": dist[o], "ne": d["NE"][o], "ni": d["ARP"][o], "Te": d["ETEMP"][o],
            "Th": d["TGAS"][o], "ji_j": (d["JIX"] / np.abs(j))[o], "j": j[o],
            "E": -d["EFLDX"][o], "phi": d["PHI"][o],           # paper orientation
            "L": d["hi"] - d["lo"], "t": float(d["ds"].current_time)}


def style_axes(ax):
    ax.set_xscale("log")
    ax.set_xlim(*XLIM)
    ax.xaxis.set_major_locator(LogLocator(base=10, numticks=12))
    ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=12))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.tick_params(which="major", length=3.5, width=0.6)
    ax.tick_params(which="minor", length=2.0, width=0.5)
    for sp in ax.spines.values():
        sp.set_linewidth(0.6)


def figure(p, s, out_stem):
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["STIX", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix", "font.size": 8.5, "axes.labelsize": 8.5,
        "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
    })
    fig, axs = plt.subplots(3, 1, figsize=(3.5, 6.6), sharex=True,
                            gridspec_kw=dict(height_ratios=[1.15, 1.0, 1.0], hspace=0.07))
    for ax in axs:
        style_axes(ax)
        # the region beyond the PlasmeX domain (arc boundary at d = L)
        ax.axvspan(s["L"], XLIM[1], color="#000000", alpha=0.045, lw=0)
        ax.axvline(s["L"], color="#8a8a8a", lw=0.6, ls=":")
    lab = LABELS

    # ---- (a) densities ----
    ax = axs[0]
    ax.plot(*p["ne"], color=C_E, **REF)
    ax.plot(*p["ni"], color=C_H, **REF)
    ax.plot(s["d"], s["ne"], color=C_E, **SIM)
    ax.plot(s["d"], s["ni"], color=C_H, **SIM)
    ax.set_yscale("log")
    ax.set_ylim(*N_YLIM)
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_ylabel(r"$n_\mathrm{e},\ n_\mathrm{i}$ (m$^{-3}$)")
    ax.text(*lab["ne"], r"$n_\mathrm{e}$", color=C_E, fontsize=9)
    ax.text(*lab["ni"], r"$n_\mathrm{i}$", color=C_H, fontsize=9)
    ax.text(s["L"] * 1.25, N_YLIM[0] * 3, r"$x > L$", fontsize=7, color="#6a6a6a",
            va="bottom")

    # ---- (b) temperatures ----
    ax = axs[1]
    ax.plot(p["Te"][0], p["Te"][1] / 1e3, color=C_E, **REF)
    ax.plot(p["Th"][0], p["Th"][1] / 1e3, color=C_H, **REF)
    ax.plot(s["d"], s["Te"] / 1e3, color=C_E, **SIM)
    ax.plot(s["d"], s["Th"] / 1e3, color=C_H, **SIM)
    ax.yaxis.set_major_locator(MultipleLocator(T_MAJOR))
    ax.yaxis.set_minor_locator(MultipleLocator(T_MAJOR / 5))
    ax.set_ylim(*T_YLIM)
    ax.set_ylabel(r"$T_\mathrm{e},\ T_\mathrm{h}$ (10$^3$ K)")
    ax.text(*lab["Te"], r"$T_\mathrm{e}$", color=C_E, fontsize=9, va="bottom")
    ax.text(*lab["Th"], r"$T_\mathrm{h}$", color=C_H, fontsize=9, va="bottom")
    handles = [Line2D([], [], color="#2b2b2b", **SIM), Line2D([], [], color="#2b2b2b", **REF)]
    ax.legend(handles, [SIM_LABEL, REF_LABEL], frameon=False, handlelength=2.6,
              borderaxespad=0.4, loc="center left", bbox_to_anchor=(0.0, 0.5))

    # ---- (c) j_i/j (left) and -E (right), one zero line ----
    ax = axs[2]
    ax.tick_params(which="both", right=False)
    ax.axhline(0.0, color="#b0b0b0", lw=0.5)
    ax.plot(*p["ji_j"], color=C_J, **REF)
    ax.plot(s["d"], s["ji_j"], color=C_J, **SIM)
    ji_lo, ji_hi = JI_YLIM
    ax.set_ylim(ji_lo, ji_hi)
    ax.yaxis.set_major_locator(MultipleLocator(JI_MAJOR))
    ax.yaxis.set_minor_locator(MultipleLocator(JI_MAJOR / 4))
    ax.set_ylabel(r"$j_\mathrm{i}/j$")
    ax.set_xlabel(r"Distance from cathode surface, $x$ (m)")
    ax.text(*lab["ji"], r"$j_\mathrm{i}/j$", color=C_J, fontsize=9, va="bottom")

    fx = ax.twinx()
    fx.plot(p["E"][0], -p["E"][1] / E_UNIT, color=C_F, **REF)
    fx.plot(s["d"], -s["E"] / E_UNIT, color=C_F, **SIM)
    fx.set_ylim(E_TOP * ji_lo / ji_hi, E_TOP)
    fx.yaxis.set_major_locator(MultipleLocator(E_MAJOR))
    fx.yaxis.set_minor_locator(MultipleLocator(E_MAJOR / 4))
    fx.tick_params(which="both", direction="in", colors=C_F)
    fx.tick_params(which="major", length=3.5, width=0.6)
    fx.tick_params(which="minor", length=2.0, width=0.5)
    fx.spines["right"].set_color(C_F)
    fx.set_ylabel(rf"$-E$ (10$^{{{int(round(np.log10(E_UNIT)))}}}$ V m$^{{-1}}$)", color=C_F)
    fx.text(*lab["E"], r"$-E$", color=C_F, fontsize=9, va="bottom")

    # the cathode surface: hatched strip + thick frame line on the left of every panel
    for i, ax in enumerate(axs):
        mark_electrode(ax, "left", label="cathode" if i == 0 else None,
                       spine_axes=(fx,) if i == 2 else ())
    for ax, tag in zip(axs, "abc"):
        ax.text(0.975, 0.93, f"({tag})", transform=ax.transAxes, ha="right", va="top",
                fontsize=9)
    # no title: the run conditions and the source figures go in the paper's caption
    fig.align_ylabels(axs)
    for ext in ("pdf", "png"):
        fig.savefig(f"{out_stem}.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)


def report(p, s):
    def at(c, xq):
        return np.interp(np.log(xq), np.log(c[0]), c[1])
    sim = {k: (s["d"], s[k]) for k in ("ne", "ni", "Te", "Th", "ji_j", "E", "phi")}
    cols = [("ne", 1, "9.2e"), ("ni", 1, "9.2e"), ("Te", 1e3, "6.2f"), ("Th", 1e3, "6.2f"),
            ("ji_j", 1, "8.3f"), ("E", E_UNIT, "7.3f"), ("phi", 1, "7.2f")]
    print(f"{'x (m)':>8s} | " + " | ".join(
        f"{k + ' p':>{len(format(0, f)) }s} {k + ' s':>{len(format(0, f))}s}"
        for k, _, f in cols))
    print(f"          (T in kK, E in 1e{int(round(np.log10(E_UNIT)))} V/m, paper "
          "orientation; phi in V)")
    for xq in (1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4):
        if xq > s["L"] * 1.0001:
            continue
        row = []
        for k, sc, f in cols:
            row.append(f"{format(at(p[k], xq) / sc, f)} {format(at(sim[k], xq) / sc, f)}")
        print(f"{xq:8.1e} | " + " | ".join(row))
    (xs,), (ps,) = p["phi_max"]
    i = np.argmax(s["phi"])
    print(f"phi max (near-cathode voltage drop): paper {ps:.2f} V at x = {xs:.2e} m (square);"
          f"  PlasmeX {s['phi'][i]:.2f} V at x = {s['d'][i]:.2e} m")
    print(f"wall j_iw/j_c: PlasmeX {s['ji_j'][0]:.3f}, Table 1 {TABLE1_JIW:.3f}")


def sanity(s, name):
    """Flag a plotfile whose conditions are not the paper's -- the comparison would be moot."""
    jc = np.median(s["j"])
    msgs = []
    if abs(jc / JC - 1) > 0.05:
        msgs.append(f"j = {jc:.3e} A/m^2, Fig. 5(b) is jc = {JC:.1e}")
    if abs(s["Th"][0] - T_WALL) > 100:
        msgs.append(f"wall T_h = {s['Th'][0]:.0f} K, the paper has Tw = {T_WALL:.0f} K")
    if np.ptp(s["j"]) > 0.01 * jc:
        msgs.append(f"j varies by {np.ptp(s['j']) / jc:.1%} across the domain (not steady?)")
    for m in msgs:
        print(f"WARNING [{name}]: {m}")


def main(argv=None):
    ap = argparse.ArgumentParser(description="Compare PlasmeX plotfile(s) with Almeida et al. "
                                 "(2008), jc = 1e7 A/m^2: Fig. 5(b) and Fig. 6 curve b.")
    ap.add_argument("plotfiles", nargs="+", help="PlasmeX plotfile directories (1-D)")
    ap.add_argument("-o", "--outdir", default=".", help="where the figures go (default: .)")
    a = ap.parse_args(argv)

    p = {**load_csv("fig5b_extracted.csv"), **load_csv("fig6b_extracted.csv")}
    os.makedirs(a.outdir, exist_ok=True)

    for pf in a.plotfiles:
        pf = pf.rstrip("/")
        if not os.path.isfile(os.path.join(pf, "Header")):
            print(f"skipping {pf}: not a plotfile (no Header)")
            continue
        name = os.path.basename(os.path.abspath(pf))
        s = load_sim(pf)
        print(f"\n=== {name}  t = {s['t']:.4e} s   j = {np.median(s['j']):.4e} A/m^2   "
              f"domain {s['L'] * 1e6:.0f} um   vs jc = {JC:.2g} (Fig. 5(b), Fig. 6 b)")
        sanity(s, name)
        report(p, s)
        stem = os.path.join(a.outdir, f"compareJ1e7_{name}")
        figure(p, s, stem)
        print(f"wrote {stem}.pdf / .png")


if __name__ == "__main__":
    main()
