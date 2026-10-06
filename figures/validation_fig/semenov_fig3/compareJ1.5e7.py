"""
PlasmeX 1-D ANODE plotfile(s) vs Semenov, Krivtsun & Reisgen (2016) J. Phys. D 49 105204 at
j = 1500 A/cm^2 = 1.5e7 A/m^2 (Ar, 1 atm, anode surface T_s = 1000 K, unified model):
Fig. 3(b), and the 1500 A/cm^2 curves of Fig. 4 (j_i/j) and Fig. 5 (E), in one figure.

    python compareJ1.5e7.py <plotfile> [plt...] [-o outdir] [--reextract] [--pdf paper.pdf]

One figure per plotfile, compareJ1.5e7_<plotfile>.{pdf,png}, three stacked panels against the
distance from the anode:
    (a) n_e, n_i, n_sa                         Fig. 3(b)
    (b) T_e, T_h                               Fig. 3(b), Table 2
    (c) j_i/j (left axis)  Fig. 4   and   E (right axis)  Fig. 5, one shared zero line

ORIENTATION AND SIGNS. The paper's x is the distance from the anode, anode on the LEFT, and
j > 0. The PlasmeX anode is the x_hi face, so the simulation is plotted against d = x_hi - x,
and both quantities map without any further sign choice: E_paper = -EFLDX, and j_i/j =
JIX / (JEX + JIX) (flipping x flips both currents, so the ratio is unchanged; ions flowing into
the anode give j_i/j < 0 in both). Near the anode both are negative, so they are plotted with
their own signs and share the zero line. n_sa for PlasmeX is fig3_compare.n_saha(T_e, T_h),
the same definition as the paper's (see fig3_compare.py).

The paper curves come from fig3b_extracted.csv (fig3_extract.py) and fig45_j1500_extracted.csv
(fig45_extract.py), each run automatically when missing. The simulation is read off the AMR
leaf cells by bohm.sample from the live PlasmeX tree (../plasmex.py), via
../almeida_fig5/compareJ1e7.load_sim.
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter, MultipleLocator

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "almeida_fig5"))
import fig3_extract  # noqa: E402
import fig45_extract  # noqa: E402
from fig3_compare import n_saha, sanity, report, C_S  # noqa: E402
from fig5_compare import XLIM, C_E, C_H, C_J, SIM, REF, style_axes  # noqa: E402
from compareJ1e7 import load_sim, C_F  # noqa: E402

CASE = dict(
    fig="b", j=1.5e7, j_Acm2=1500,             # fig3_compare.sanity reads these and L_min
    L_min=2e-3,                                # uniform plasma from x ~ 1-2 mm (Fig. 3b)
    n_ylim=(1e19, 1e24),
    T_ylim=(0, 19), T_major=4,                 # kK
    ji_ylim=(-0.025, 0.005), ji_major=0.005,
    E_unit=1e6, E_bottom=-2.0, E_major=0.5,    # E axis; its zero is j_i/j's zero
    labels=dict(ne=(1.4e-8, 4.0e20), ni=(1.4e-8, 1.4e21), nsa=(1.6e-6, 1.3e20),
                Te=(1.4e-8, 5.9), Th=(1.4e-8, 1.3), ji=(1.2e-8, -0.0047), E=(1.2e-8, -1.2)),
)


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
        if s["L"] < XLIM[1]:
            ax.axvspan(s["L"], XLIM[1], color="#000000", alpha=0.045, lw=0)
            ax.axvline(s["L"], color="#8a8a8a", lw=0.6, ls=":")
    lab = case["labels"]

    # ---- (a) densities ----
    ax = axs[0]
    for k, c in (("ne", C_E), ("ni", C_H), ("nsa", C_S)):
        ax.plot(*p[k], color=c, **REF)
        ax.plot(s["d"], s[k], color=c, **SIM)
    ax.set_yscale("log")
    ax.set_ylim(*case["n_ylim"])
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_ylabel(r"$n_\mathrm{e},\ n_\mathrm{i},\ n_\mathrm{sa}$ (m$^{-3}$)")
    ax.text(*lab["ne"], r"$n_\mathrm{e}$", color=C_E, fontsize=9)
    ax.text(*lab["ni"], r"$n_\mathrm{i}$", color=C_H, fontsize=9)
    ax.text(*lab["nsa"], r"$n_\mathrm{sa}$", color=C_S, fontsize=9)
    if s["L"] < XLIM[1] / 3:
        ax.text(s["L"] * 1.25, case["n_ylim"][0] * 2, "outside\nPlasmeX\ndomain",
                fontsize=6.5, color="#6a6a6a", va="bottom", linespacing=1.0)
    handles = [Line2D([], [], color="#2b2b2b", **SIM), Line2D([], [], color="#2b2b2b", **REF)]
    ax.legend(handles, ["PlasmeX", "Semenov et al. (2016)"], frameon=False, handlelength=2.6,
              borderaxespad=0.4, loc="upper left", bbox_to_anchor=(0.0, 0.88))

    # ---- (b) temperatures ----
    ax = axs[1]
    ax.plot(p["Te"][0], p["Te"][1] / 1e3, color=C_E, **REF)
    ax.plot(p["Th"][0], p["Th"][1] / 1e3, color=C_H, **REF)
    ax.plot(s["d"], s["Te"] / 1e3, color=C_E, **SIM)
    ax.plot(s["d"], s["Th"] / 1e3, color=C_H, **SIM)
    ax.set_ylim(*case["T_ylim"])
    ax.yaxis.set_major_locator(MultipleLocator(case["T_major"]))
    ax.yaxis.set_minor_locator(MultipleLocator(case["T_major"] / 4))
    ax.set_ylabel(r"$T_\mathrm{e},\ T_\mathrm{h}$ (10$^3$ K)")
    ax.text(*lab["Te"], r"$T_\mathrm{e}$", color=C_E, fontsize=9, va="bottom")
    ax.text(*lab["Th"], r"$T_\mathrm{h}$", color=C_H, fontsize=9, va="bottom")
    tw, tb = fig3_extract.TABLE2[case["j_Acm2"]]
    ax.text(0.03, 0.55,
            rf"$T_\mathrm{{e}}$ at the anode: {s['Te'][0] / 1e3:.2f} kK PlasmeX, "
            rf"{tw / 1e3:.2f} kK ref." "\n"
            rf"$T_\mathrm{{e}}$ at $x = L$: {s['Te'][-1] / 1e3:.2f} kK PlasmeX, "
            rf"{tb / 1e3:.2f} kK ref. (Table 2)",
            transform=ax.transAxes, fontsize=6.6, color="#3a3a3a", linespacing=1.3)

    # ---- (c) j_i/j (left) and E (right), one zero line ----
    ax = axs[2]
    ax.tick_params(which="both", right=False)
    ax.axhline(0.0, color="#b0b0b0", lw=0.5)
    ax.plot(*p["ji_j"], color=C_J, **REF)
    ax.plot(s["d"], s["ji_j"], color=C_J, **SIM)
    ji_lo, ji_hi = case["ji_ylim"]
    ax.set_ylim(ji_lo, ji_hi)
    ax.yaxis.set_major_locator(MultipleLocator(case["ji_major"]))
    ax.yaxis.set_minor_locator(MultipleLocator(case["ji_major"] / 5))
    ax.set_ylabel(r"$j_\mathrm{i}/j$")
    ax.set_xlabel(r"Distance from anode surface, $x$ (m)")
    ax.text(*lab["ji"], r"$j_\mathrm{i}/j$", color=C_J, fontsize=9, va="bottom")

    fx = ax.twinx()
    u, bot = case["E_unit"], case["E_bottom"]
    fx.plot(p["E"][0], p["E"][1] / u, color=C_F, **REF)
    fx.plot(s["d"], s["E"] / u, color=C_F, **SIM)
    fx.set_ylim(bot, bot * ji_hi / ji_lo)
    fx.yaxis.set_major_locator(MultipleLocator(case["E_major"]))
    fx.yaxis.set_minor_locator(MultipleLocator(case["E_major"] / 5))
    fx.tick_params(which="both", direction="in", colors=C_F)
    fx.tick_params(which="major", length=3.5, width=0.6)
    fx.tick_params(which="minor", length=2.0, width=0.5)
    fx.spines["right"].set_color(C_F)
    fx.set_ylabel(rf"$E$ (10$^{{{int(round(np.log10(u)))}}}$ V m$^{{-1}}$)", color=C_F)
    fx.text(*lab["E"], r"$E$", color=C_F, fontsize=9, va="bottom")

    for ax, tag in zip(axs, "abc"):
        ax.text(0.975 if tag == "c" else 0.03, 0.93 if tag != "c" else 0.07, f"({tag})",
                transform=ax.transAxes, ha="right" if tag == "c" else "left",
                va="top" if tag != "c" else "bottom", fontsize=9)
    axs[0].set_title(r"PlasmeX vs Semenov et al. (2016), Fig. 3(b), 4 and 5" "\n"
                     rf"Ar, $p_0$ = 1 atm, $T_\mathrm{{s}}$ = 1000 K, "
                     rf"$j$ = {case['j_Acm2']} A cm$^{{-2}}$",
                     fontsize=8.2, pad=4, linespacing=1.35)
    fig.align_ylabels(axs)
    for ext in ("pdf", "png"):
        fig.savefig(f"{out_stem}.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)


def report_j_e(case, p, s):
    def at(c, xq):
        return np.interp(np.log(xq), np.log(c[0]), c[1])
    u = case["E_unit"]
    print(f"{'x (m)':>9s} | {'ji/j p':>9s} {'ji/j s':>9s} | {'E p':>8s} {'E s':>8s}   "
          f"(E in 1e{int(round(np.log10(u)))} V/m, paper orientation)")
    for xq in (1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 3e-4, 1e-3, 3e-3):
        if xq > s["L"] * 1.0001:
            continue
        e_p = at(p["E"], xq) / u if xq <= p["E"][0].max() else np.nan
        print(f"{xq:9.1e} | {at(p['ji_j'], xq):9.5f} {at((s['d'], s['ji_j']), xq):9.5f} | "
              f"{e_p:8.4f} {at((s['d'], s['E']), xq) / u:8.4f}")
    i, k = np.argmin(p["ji_j"][1]), np.argmin(s["ji_j"])
    print(f"j_i/j minimum: paper {p['ji_j'][1][i]:.5f} at x = {p['ji_j'][0][i]:.2e} m;  "
          f"PlasmeX {s['ji_j'][k]:.5f} at x = {s['d'][k]:.2e} m")


def main(case=CASE, argv=None):
    ap = argparse.ArgumentParser(
        description="Compare PlasmeX anode plotfile(s) with Semenov et al. (2016) at "
                    "j = 1500 A/cm^2: Fig. 3(b), Fig. 4 (j_i/j) and Fig. 5 (E).")
    ap.add_argument("plotfiles", nargs="+", help="PlasmeX plotfile directories (1-D, anode)")
    ap.add_argument("-o", "--outdir", default=".", help="where the figures go (default: .)")
    ap.add_argument("--pdf", default=None, help="path to the Semenov 2016 PDF (only needed "
                    "if an extracted CSV is missing or --reextract is given)")
    ap.add_argument("--reextract", action="store_true",
                    help="re-extract the paper curves from the PDF (writes the CSVs and the "
                         "verification figures)")
    a = ap.parse_args(argv)

    pdf = a.pdf and os.path.abspath(a.pdf)
    if a.reextract:
        p = {**fig3_extract.extract(case["fig"], pdf), **fig45_extract.extract(pdf)}
    else:
        p = {**fig3_extract.load(case["fig"], pdf), **fig45_extract.load(pdf)}
    os.makedirs(a.outdir, exist_ok=True)

    for pf in a.plotfiles:
        pf = pf.rstrip("/")
        if not os.path.isfile(os.path.join(pf, "Header")):
            print(f"skipping {pf}: not a plotfile (no Header)")
            continue
        name = os.path.basename(os.path.abspath(pf))
        s = load_sim(pf)
        s["nsa"] = n_saha(s["Te"], s["Th"])
        print(f"\n=== {name}  t = {s['t']:.4e} s   j = {np.median(s['j']):.4e} A/m^2   "
              f"domain {s['L'] * 1e3:.3f} mm   vs Semenov j = 1500 A/cm^2")
        sanity(case, s, name)
        report(p, s)
        report_j_e(case, p, s)
        stem = os.path.join(a.outdir, f"compareJ1.5e7_{name}")
        figure(case, p, s, stem)
        print(f"wrote {stem}.pdf / .png")


if __name__ == "__main__":
    main()
