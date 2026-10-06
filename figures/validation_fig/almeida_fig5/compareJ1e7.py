"""
PlasmeX 1-D cathode plotfile(s) vs Almeida, Benilov & Naidis (2008) J. Phys. D 41 245201 at
jc = 1e7 A/m^2 (Ar, p0 = 1 bar, Tw = 3500 K): Fig. 5(b) and curve b of Fig. 6, in one figure.

    python compareJ1e7.py ~/work/PlasmeX/exe/cathode_plt05000 [plt...] [-o outdir] [--reextract]

One figure per plotfile, compareJ1e7_<plotfile>.{pdf,png}, three stacked panels against the
distance from the cathode:
    (a) n_e, n_i            Fig. 5(b)
    (b) T_e, T_h            Fig. 5(b)
    (c) j_i/j (left axis)   Fig. 5(b)   and   -E (right axis)   Fig. 6 curve b
A terminal table also lists the potential phi (Fig. 6 b), which is not plotted.

The engine is driven by a CASE dict (the paper panel, where its CSVs live, the axes and label
positions); ../almeida_J7.8e7/compareJ7.8e7.py runs the same engine on Fig. 5(c) / Fig. 6 c.

ORIENTATION AND SIGNS. The paper's x is the distance from the cathode, cathode on the LEFT. The
PlasmeX cathode is the x_hi face, so the simulation is plotted against d = x_hi - x. The
paper's E is along its x, i.e. E_paper = -EFLDX, negative in the sheath (the field points to
the cathode). Panel (c) plots -E_paper = EFLDX, the field magnitude towards the cathode, so
that it and j_i/j are both positive and share one zero line. PHI is zero on the cathode face,
as the paper's phi is.

The paper curves come from fig5b_extracted.csv and fig6b_extracted.csv (fig5_extract.py,
fig6_extract.py; run automatically when a CSV is missing). The simulation is read off the AMR
leaf cells by bohm.sample, loaded from the live PlasmeX tree (../plasmex.py).
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter, MultipleLocator, FuncFormatter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
import fig5_extract  # noqa: E402
import fig6_extract  # noqa: E402
from fig5_compare import XLIM, C_E, C_H, C_J, SIM, REF, style_axes, sanity  # noqa: E402
from plasmex import exe_module  # noqa: E402

C_F = "#1a9e77"                               # the electric field

CASE = dict(
    fig="b", jc=1e7, tag="J1e7", table1_jiw=0.216,   # fig5_compare.sanity reads fig and jc
    datadir=HERE,                                   # where fig5b/fig6b_extracted.csv live
    n_ylim=(1e19, 1e24),
    T_log=False, T_ylim=(0, 21), T_major=5,          # kK
    ji_ylim=(-0.05, 0.35), ji_major=0.1,
    E_unit=1e7, E_top=2.5, E_major=1.0,              # -E axis; its zero is j_i/j's zero
    labels=dict(ne=(1.4e-8, 5.5e19), ni=(1.4e-8, 4.5e21), Te=(1.4e-8, 18.8), Th=(1.4e-8, 4.3),
                ji=(2.0e-6, 0.226), E=(1.2e-8, 2.2)),
    wall_note_at=(0.975, 0.62),                      # axes fraction, right-aligned, top
)


def load_sim(plt_file):
    b = exe_module("bohm")
    d = b.sample(plt_file, ("NE", "ARP", "ETEMP", "TGAS", "JEX", "JIX", "PHI", "EFLDX"),
                 want_derived=False)
    dist = d["hi"] - d["x"]
    o = np.argsort(dist)
    j = d["JEX"] + d["JIX"]
    return {"d": dist[o], "ne": d["NE"][o], "ni": d["ARP"][o], "Te": d["ETEMP"][o],
            "Th": d["TGAS"][o], "ji_j": (d["JIX"] / j)[o], "j": j[o],
            "E": -d["EFLDX"][o], "phi": d["PHI"][o],           # paper orientation
            "L": d["hi"] - d["lo"], "t": float(d["ds"].current_time)}


def figure(case, p, s, out_stem):
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["STIX", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix", "font.size": 8.5, "axes.labelsize": 8.5,
        "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
    })
    fig, axs = plt.subplots(3, 1, figsize=(3.5, 6.6), sharex=True,
                            gridspec_kw=dict(height_ratios=[1.15, 1.0, 1.0], hspace=0.07))
    for ax in axs:
        style_axes(ax)
        ax.axvspan(s["L"], XLIM[1], color="#000000", alpha=0.045, lw=0)
        ax.axvline(s["L"], color="#8a8a8a", lw=0.6, ls=":")

    # ---- (a) densities ----
    ax = axs[0]
    ax.plot(*p["ne"], color=C_E, **REF)
    ax.plot(*p["ni"], color=C_H, **REF)
    ax.plot(s["d"], s["ne"], color=C_E, **SIM)
    ax.plot(s["d"], s["ni"], color=C_H, **SIM)
    ax.set_yscale("log")
    ax.set_ylim(*case["n_ylim"])
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_ylabel(r"$n_\mathrm{e},\ n_\mathrm{i}$ (m$^{-3}$)")
    lab = case["labels"]
    ax.text(*lab["ne"], r"$n_\mathrm{e}$", color=C_E, fontsize=9)
    ax.text(*lab["ni"], r"$n_\mathrm{i}$", color=C_H, fontsize=9)
    ax.text(s["L"] * 1.25, case["n_ylim"][0] * 3, "outside\nPlasmeX\ndomain", fontsize=6.5,
            color="#6a6a6a", va="bottom", linespacing=1.0)

    # ---- (b) temperatures ----
    ax = axs[1]
    ax.plot(p["Te"][0], p["Te"][1] / 1e3, color=C_E, **REF)
    ax.plot(p["Th"][0], p["Th"][1] / 1e3, color=C_H, **REF)
    ax.plot(s["d"], s["Te"] / 1e3, color=C_E, **SIM)
    ax.plot(s["d"], s["Th"] / 1e3, color=C_H, **SIM)
    if case["T_log"]:
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.yaxis.set_minor_formatter(NullFormatter())
    else:
        ax.yaxis.set_major_locator(MultipleLocator(case["T_major"]))
        ax.yaxis.set_minor_locator(MultipleLocator(case["T_major"] / 5))
    ax.set_ylim(*case["T_ylim"])
    ax.set_ylabel(r"$T_\mathrm{e},\ T_\mathrm{h}$ (10$^3$ K)")
    ax.text(*lab["Te"], r"$T_\mathrm{e}$", color=C_E, fontsize=9, va="bottom")
    ax.text(*lab["Th"], r"$T_\mathrm{h}$", color=C_H, fontsize=9, va="bottom")
    handles = [Line2D([], [], color="#2b2b2b", **SIM), Line2D([], [], color="#2b2b2b", **REF)]
    ax.legend(handles, ["PlasmeX", "Almeida et al. (2008)"], frameon=False, handlelength=2.6,
              borderaxespad=0.4, loc="center left", bbox_to_anchor=(0.0, 0.5))

    # ---- (c) j_i/j (left) and -E (right), one zero line ----
    ax = axs[2]
    ax.tick_params(which="both", right=False)
    ax.axhline(0.0, color="#b0b0b0", lw=0.5)
    ax.plot(*p["ji_j"], color=C_J, **REF)
    ax.plot(s["d"], s["ji_j"], color=C_J, **SIM)
    ji_lo, ji_hi = case["ji_ylim"]
    ax.set_ylim(ji_lo, ji_hi)
    ax.yaxis.set_major_locator(MultipleLocator(case["ji_major"]))
    ax.yaxis.set_minor_locator(MultipleLocator(case["ji_major"] / 4))
    ax.set_ylabel(r"$j_\mathrm{i}/j$")
    ax.set_xlabel(r"Distance from cathode surface, $x$ (m)")
    ax.text(*lab["ji"], r"$j_\mathrm{i}/j$", color=C_J, fontsize=9, va="bottom")
    ax.text(*case["wall_note_at"],
            rf"wall $j_\mathrm{{iw}}/j_\mathrm{{c}}$" "\n" rf"PlasmeX {s['ji_j'][0]:.3g}" "\n"
            rf"Table 1 {case['table1_jiw']:.3g}",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.8, color="#3a3a3a",
            linespacing=1.3)

    fx = ax.twinx()
    u, top = case["E_unit"], case["E_top"]
    fx.plot(p["E"][0], -p["E"][1] / u, color=C_F, **REF)
    fx.plot(s["d"], -s["E"] / u, color=C_F, **SIM)
    fx.set_ylim(top * ji_lo / ji_hi, top)
    fx.yaxis.set_major_locator(MultipleLocator(case["E_major"]))
    fx.yaxis.set_minor_locator(MultipleLocator(case["E_major"] / 4))
    fx.tick_params(which="both", direction="in", colors=C_F)
    fx.tick_params(which="major", length=3.5, width=0.6)
    fx.tick_params(which="minor", length=2.0, width=0.5)
    fx.spines["right"].set_color(C_F)
    fx.set_ylabel(rf"$-E$ (10$^{{{int(round(np.log10(u)))}}}$ V m$^{{-1}}$)", color=C_F)
    fx.text(*lab["E"], r"$-E$", color=C_F, fontsize=9, va="bottom")

    for ax, tag in zip(axs, "abc"):
        ax.text(0.975, 0.93, f"({tag})", transform=ax.transAxes, ha="right", va="top",
                fontsize=9)
    mant, k = f"{np.median(s['j']):.2e}".split("e")
    axs[0].set_title(rf"PlasmeX vs Almeida et al. (2008), Fig. 5({case['fig']}) and "
                     rf"Fig. 6 ({case['fig']})" "\n"
                     rf"Ar, $p_0$ = 1 bar, $T_\mathrm{{w}}$ = 3500 K, "
                     rf"$j_\mathrm{{c}}$ = {mant}$\times$10$^{{{int(k)}}}$ A m$^{{-2}}$",
                     fontsize=8.2, pad=4, linespacing=1.35)
    fig.align_ylabels(axs)
    for ext in ("pdf", "png"):
        fig.savefig(f"{out_stem}.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)


def report(case, p, s):
    def at(c, xq):
        return np.interp(np.log(xq), np.log(c[0]), c[1])
    sim = {k: (s["d"], s[k]) for k in ("ne", "ni", "Te", "Th", "ji_j", "E", "phi")}
    cols = [("ne", 1, "9.2e"), ("ni", 1, "9.2e"), ("Te", 1e3, "6.2f"), ("Th", 1e3, "6.2f"),
            ("ji_j", 1, "8.3f"), ("E", case["E_unit"], "7.3f"), ("phi", 1, "7.2f")]
    print(f"{'x (m)':>8s} | " + " | ".join(
        f"{k + ' p':>{len(format(0, f)) }s} {k + ' s':>{len(format(0, f))}s}"
        for k, _, f in cols))
    print(f"          (T in kK, E in 1e{int(round(np.log10(case['E_unit'])))} V/m, paper "
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


def main(case=CASE, argv=None):
    fg = case["fig"]
    ap = argparse.ArgumentParser(description="Compare PlasmeX plotfile(s) with Almeida et al. "
                                 f"(2008), jc = {case['jc']:.2g} A/m^2: Fig. 5({fg}) and "
                                 f"Fig. 6 curve {fg}.")
    ap.add_argument("plotfiles", nargs="+", help="PlasmeX plotfile directories (1-D)")
    ap.add_argument("-o", "--outdir", default=".", help="where the figures go (default: .)")
    ap.add_argument("--pdf", default=None, help="path to the Almeida 2008 PDF (only needed "
                    "if an extracted CSV is missing or --reextract is given)")
    ap.add_argument("--reextract", action="store_true",
                    help="re-extract the paper curves from the PDF (writes the CSVs and the "
                         "verification figures)")
    a = ap.parse_args(argv)

    pdf = a.pdf and os.path.abspath(a.pdf)
    dd = case["datadir"]
    if a.reextract:
        p = {**fig5_extract.extract(fg, pdf, outdir=dd), **fig6_extract.extract(fg, pdf, outdir=dd)}
    else:
        p = {**fig5_extract.load(fg, pdf, outdir=dd), **fig6_extract.load(fg, pdf, outdir=dd)}
    os.makedirs(a.outdir, exist_ok=True)

    for pf in a.plotfiles:
        pf = pf.rstrip("/")
        if not os.path.isfile(os.path.join(pf, "Header")):
            print(f"skipping {pf}: not a plotfile (no Header)")
            continue
        name = os.path.basename(os.path.abspath(pf))
        s = load_sim(pf)
        print(f"\n=== {name}  t = {s['t']:.4e} s   j = {np.median(s['j']):.4e} A/m^2   "
              f"domain {s['L'] * 1e6:.0f} um   vs jc = {case['jc']:.2g} (Fig. 5({fg}), "
              f"Fig. 6 {fg})")
        sanity(case, s, name)
        report(case, p, s)
        stem = os.path.join(a.outdir, f"compare{case['tag']}_{name}")
        figure(case, p, s, stem)
        print(f"wrote {stem}.pdf / .png")


if __name__ == "__main__":
    main()
