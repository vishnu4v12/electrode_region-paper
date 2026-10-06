"""
Shared engine of compare3a.py / compare3b.py: a PlasmeX 1-D ANODE plotfile against
Semenov, Krivtsun & Reisgen (2016) J. Phys. D 49 105204, Fig. 3(a) or 3(b).

One journal-style figure per plotfile, two stacked panels on a common axis:
(a) n_e, n_i, n_sa   (b) T_e, T_h.

ORIENTATION. The paper's x is the distance from the anode surface, anode on the LEFT. In
PlasmeX the electrode is the x_hi face, so the simulation is plotted against d = x_hi - x,
exactly as ../almeida_fig5/fig5_compare.py does for the cathode (whose sampler, axis style
and line conventions are reused here).

n_sa is "the charged particle number density calculated by means of the Saha equation and
quasi-neutrality condition" (the paper's caption). For PlasmeX it is computed here from the
local T_e, T_h at p0 = 1 atm with a two-temperature Saha law (ionization governed by T_e)
and Dalton's law p = n_e k T_e + (n_i + n_a) k T_h, with the Ar+ internal partition function
4 + 2 exp(-2062 K / T_e) and Z_Ar = 1. Applied to the paper's OWN T_e, T_h curves this
reproduces the paper's n_sa to 1-4 % over the whole layer (a constant g+/g1 = 6 is 1-9 %
high), so the two n_sa curves are defined the same way.

The paper curves come from fig3{a,b}_extracted.csv (fig3_extract.py, run automatically
when the CSV is missing).
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
from fig5_compare import load_sim, style_axes, XLIM, SIM, REF, C_E, C_H  # noqa: E402

C_S = "#3f9a3f"                       # n_sa
P0, T_WALL = 101325.0, 1000.0         # the conditions every Fig. 3 panel shares

KB, ME, HP, QE = 1.380649e-23, 9.1093837015e-31, 6.62607015e-34, 1.602176634e-19
I_AR = 15.7596                        # [eV]


def n_saha(Te, Th, p=P0):
    """Charged-particle density from Saha (at T_e) + quasi-neutrality + Dalton (see doc)."""
    Te, Th = np.asarray(Te, float), np.asarray(Th, float)
    g = 4.0 + 2.0 * np.exp(-2062.0 / Te)
    S = 2.0 * g * (2 * np.pi * ME * KB * Te / HP**2) ** 1.5 * np.exp(-I_AR * QE / (KB * Te))
    a, c = (Te + Th) / Th, p / (KB * Th)        # n^2 = S (c - a n)
    return 0.5 * (-S * a + np.sqrt((S * a) ** 2 + 4 * S * c))


def figure(case, p, s, out_stem):
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["STIX", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix", "font.size": 8.5, "axes.labelsize": 8.5,
        "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
    })
    fig, axs = plt.subplots(2, 1, figsize=(3.5, 4.9), sharex=True,
                            gridspec_kw=dict(height_ratios=[1.15, 1.0], hspace=0.07))
    for ax in axs:
        style_axes(ax)
        if s["L"] < XLIM[1]:
            ax.axvspan(s["L"], XLIM[1], color="#000000", alpha=0.045, lw=0)
            ax.axvline(s["L"], color="#8a8a8a", lw=0.6, ls=":")
    lab = case["labels"]

    # ---- (a) densities ----
    ax = axs[0]
    ax.plot(*p["ne"], color=C_E, **REF)
    ax.plot(*p["ni"], color=C_H, **REF)
    ax.plot(*p["nsa"], color=C_S, **REF)
    ax.plot(s["d"], s["ne"], color=C_E, **SIM)
    ax.plot(s["d"], s["ni"], color=C_H, **SIM)
    ax.plot(s["d"], s["nsa"], color=C_S, **SIM)
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
    ax.set_xlabel(r"Distance from anode surface, $x$ (m)")
    tw, tb = fig3_extract.TABLE2[case["j_Acm2"]]
    ax.text(*case["wall_note_at"],
            rf"$T_\mathrm{{e}}$ at the anode: {s['Te'][0] / 1e3:.2f} kK PlasmeX, "
            rf"{tw / 1e3:.2f} kK ref." "\n"
            rf"$T_\mathrm{{e}}$ at $x = L$: {s['Te'][-1] / 1e3:.2f} kK PlasmeX, "
            rf"{tb / 1e3:.2f} kK ref. (Table 2)",
            transform=ax.transAxes, fontsize=6.6, color="#3a3a3a", linespacing=1.3)

    for ax, tag in zip(axs, "ab"):
        ax.text(0.03, 0.93, f"({tag})", transform=ax.transAxes, ha="left", va="top",
                fontsize=9)

    handles = [Line2D([], [], color="#2b2b2b", **SIM), Line2D([], [], color="#2b2b2b", **REF)]
    lg_ax, lg_kw = case["legend"]
    axs[lg_ax].legend(handles, ["PlasmeX", "Semenov et al. (2016)"], frameon=False,
                      handlelength=2.6, borderaxespad=0.4, **lg_kw)
    axs[0].set_title(rf"PlasmeX vs Semenov et al. (2016), Fig. 3({case['fig']})" "\n"
                     rf"Ar, $p_0$ = 1 atm, $T_\mathrm{{s}}$ = 1000 K, "
                     rf"$j$ = {case['j_Acm2']} A cm$^{{-2}}$",
                     fontsize=8.2, pad=4, linespacing=1.35)

    fig.align_ylabels(axs)
    for ext in ("pdf", "png"):
        fig.savefig(f"{out_stem}.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)


def report(p, s):
    def at(x, y, xq):
        return np.interp(np.log(xq), np.log(x), y)
    print(f"{'x (m)':>9s} | {'ne paper':>9s} {'ne sim':>9s} | {'ni paper':>9s} {'ni sim':>9s} | "
          f"{'nsa p':>9s} {'nsa s':>9s} | {'Te p':>6s} {'Te s':>6s} | {'Th p':>6s} {'Th s':>6s}")
    for xq in (1e-7, 1e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2):
        if xq > s["L"] * 1.0001:
            continue
        nsa_p = at(*p["nsa"], xq) if xq >= p["nsa"][0].min() else np.nan
        print(f"{xq:9.1e} | {at(*p['ne'], xq):9.2e} {at(s['d'], s['ne'], xq):9.2e} | "
              f"{at(*p['ni'], xq):9.2e} {at(s['d'], s['ni'], xq):9.2e} | "
              f"{nsa_p:9.2e} {at(s['d'], s['nsa'], xq):9.2e} | "
              f"{at(*p['Te'], xq) / 1e3:6.2f} {at(s['d'], s['Te'], xq) / 1e3:6.2f} | "
              f"{at(*p['Th'], xq) / 1e3:6.2f} {at(s['d'], s['Th'], xq) / 1e3:6.2f}")


def sanity(case, s, name):
    """Flag a plotfile whose conditions are not the figure's -- the comparison would be moot."""
    j = np.median(np.abs(s["j"]))
    msgs = []
    if abs(j / case["j"] - 1) > 0.05:
        msgs.append(f"|j| = {j:.3e} A/m^2, Fig. 3({case['fig']}) is j = {case['j']:.2e}")
    if abs(s["Th"][0] - T_WALL) > 50:
        msgs.append(f"wall T_h = {s['Th'][0]:.0f} K, the paper has T_s = {T_WALL:.0f} K")
    tb = fig3_extract.TABLE2[case["j_Acm2"]][1]
    if abs(s["Te"][-1] / tb - 1) > 0.02:
        msgs.append(f"T_e at x = L is {s['Te'][-1]:.0f} K, the paper's undisturbed plasma is "
                    f"{tb:.0f} K (Table 2) -- set prob.arc_Te (semenov_arc_bc.py --deck)")
    if s["L"] < case["L_min"]:
        msgs.append(f"domain L = {s['L'] * 1e3:.2f} mm: the paper's plasma is uniform only "
                    f"beyond ~{case['L_min'] * 1e3:.0f} mm, so the arc-side Dirichlet state "
                    f"sits inside the layer")
    if np.ptp(s["j"]) > 0.01 * j:
        msgs.append(f"j varies by {np.ptp(s['j']) / j:.1%} across the domain (not steady?)")
    for m in msgs:
        print(f"WARNING [{name}]: {m}")


def main(case, argv=None):
    ap = argparse.ArgumentParser(
        description=f"Compare PlasmeX anode plotfile(s) with Semenov et al. (2016) "
                    f"Fig. 3({case['fig']}).")
    ap.add_argument("plotfiles", nargs="+", help="PlasmeX plotfile directories (1-D, anode)")
    ap.add_argument("-o", "--outdir", default=".", help="where the figures go (default: .)")
    ap.add_argument("--pdf", default=None, help="path to the Semenov 2016 PDF (only needed "
                    "if the extracted CSV is missing or --reextract is given)")
    ap.add_argument("--reextract", action="store_true",
                    help="re-extract the paper curves from the PDF (writes the CSV and the "
                         "verification figure)")
    a = ap.parse_args(argv)

    fig = case["fig"]
    pdf = a.pdf and os.path.abspath(a.pdf)
    p = fig3_extract.extract(fig, pdf) if a.reextract else fig3_extract.load(fig, pdf)
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
              f"domain {s['L'] * 1e3:.3f} mm   vs Fig. 3({fig})")
        sanity(case, s, name)
        report(p, s)
        stem = os.path.join(a.outdir, f"compare3{fig}_{name}")
        figure(case, p, s, stem)
        print(f"wrote {stem}.pdf / .png")
