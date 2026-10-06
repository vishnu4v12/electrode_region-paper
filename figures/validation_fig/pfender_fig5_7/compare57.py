"""
PlasmeX 1-D anode plotfile(s) vs Dinulescu & Pfender (1980) J. Appl. Phys. 51 3149,
Figs. 5 and 7: Ar, 100 kPa, j = 900 A/cm^2 (9e6 A/m^2), anode T_w = 400 K; at the arc end
x = L (0.77 mm) T = T_e = 12 400 K with dT/dx = 2.754e6 K/m (not a uniform plasma).

    python compare57.py plt... [-o outdir] [--logx] [--reextract] [--pdf paper.pdf]

One figure per plotfile, compare57_<plotfile>.{pdf,png}, two stacked panels on a common axis:
(a) n_e, n_sa (Fig. 7)   (b) T_e, T_h (Fig. 5), against the distance from the anode
(PlasmeX anode = x_hi, plotted as d = x_hi - x, as ../semenov_fig3/fig3_compare.py does).

Axes: by default the paper's own (linear x in mm, n in 1e22 m^-3, T in kK), where the
sheath is invisible; --logx uses the log-log axes of the Semenov comparison instead.

n_sa is the paper's n_e' ("n_e is computed from the Saha equation at the temperature T_e",
Sec. III): Saha at T_e with p = (2 n_e + n_a) k T_e, the Ar+ partition function and p = 100
kPa. On the paper's own T_e this reproduces its n_e' to -5 ... +17 % (fig57_extract.py);
a two-temperature Saha (Semenov's n_sa) is up to 3x off near the anode.

The paper curves come from fig{5,7}_extracted.csv (fig57_extract.py, run automatically when
a CSV is missing).
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator, NullFormatter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "semenov_fig3"))
sys.path.insert(0, os.path.join(HERE, "..", "almeida_fig5"))
import fig57_extract as fx  # noqa: E402
from fig3_compare import n_saha, C_S  # noqa: E402
from fig5_compare import load_sim, style_axes, XLIM, SIM, REF, C_E, C_H  # noqa: E402

L_REF = 0.769e-3                      # the paper's layer thickness (both figures)


def figure(p, s, out_stem, logx=False):
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["STIX", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix", "font.size": 8.5, "axes.labelsize": 8.5,
        "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
    })
    fig, axs = plt.subplots(2, 1, figsize=(3.5, 4.9), sharex=True,
                            gridspec_kw=dict(height_ratios=[1.0, 1.0], hspace=0.07))
    # x in m on log axes, mm on linear ones
    xs = 1.0 if logx else 1e3
    x_hi = max(1.0, min(2.0, 1.05 * s["L"] * 1e3))     # mm, linear axes only
    for ax in axs:
        if logx:
            style_axes(ax)
        else:
            ax.set_xlim(0, x_hi)
            ax.xaxis.set_major_locator(MultipleLocator(0.2))
            ax.xaxis.set_minor_locator(MultipleLocator(0.05))
            ax.tick_params(which="both", direction="in", top=True, right=True)
            ax.tick_params(which="major", length=3.5, width=0.6)
            ax.tick_params(which="minor", length=2.0, width=0.5)
            for sp in ax.spines.values():
                sp.set_linewidth(0.6)
        ax.axvline(L_REF * xs, color="#8a8a8a", lw=0.6, ls=":")
        if s["L"] * xs < ax.get_xlim()[1]:
            ax.axvspan(s["L"] * xs, ax.get_xlim()[1], color="#000000", alpha=0.045, lw=0)

    # ---- (a) densities (Fig. 7) ----
    ax = axs[0]
    nu = 1.0 if logx else 1e22
    ax.plot(p["ne"][0] * xs, p["ne"][1] / nu, color=C_E, **REF)
    ax.plot(p["nsa"][0] * xs, p["nsa"][1] / nu, color=C_S, **REF)
    ax.plot(s["d"] * xs, s["ne"] / nu, color=C_E, **SIM)
    ax.plot(s["d"] * xs, s["nsa"] / nu, color=C_S, **SIM)
    if logx:
        ax.set_yscale("log")
        ax.set_ylim(1e19, 1e24)
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.set_ylabel(r"$n_\mathrm{e},\ n_\mathrm{sa}$ (m$^{-3}$)")
    else:
        top = max(p["ne"][1].max(), s["ne"].max(), s["nsa"].max()) / nu
        ax.set_ylim(0, 2 * np.ceil(max(9.0, 1.05 * top) / 2))
        ax.yaxis.set_major_locator(MultipleLocator(2))
        ax.yaxis.set_minor_locator(MultipleLocator(0.5))
        ax.set_ylabel(r"$n_\mathrm{e},\ n_\mathrm{sa}$ (10$^{22}$ m$^{-3}$)")

    # ---- (b) temperatures (Fig. 5) ----
    ax = axs[1]
    ax.plot(p["Te"][0] * xs, p["Te"][1] / 1e3, color=C_E, **REF)
    ax.plot(p["Th"][0] * xs, p["Th"][1] / 1e3, color=C_H, **REF)
    ax.plot(s["d"] * xs, s["Te"] / 1e3, color=C_E, **SIM)
    ax.plot(s["d"] * xs, s["Th"] / 1e3, color=C_H, **SIM)
    ax.set_ylim(0, 18)                                # Fig. 5 runs to 1.8e4 K
    ax.yaxis.set_major_locator(MultipleLocator(3))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.set_ylabel(r"$T_\mathrm{e},\ T_\mathrm{h}$ (10$^3$ K)")
    ax.set_xlabel(r"Distance from anode surface, $x$ " + ("(m)" if logx else "(mm)"))
    ax.text(0.97, 0.95 if logx else 0.05,
            rf"$T_\mathrm{{e}}$ at the anode: {s['Te'][0] / 1e3:.2f} kK PlasmeX, "
            rf"{fx.TEXT['Te_anode'] / 1e3:.2f} kK ref." "\n"
            rf"$T_\mathrm{{e}}$ at $x = L$: {s['Te'][-1] / 1e3:.2f} kK PlasmeX, "
            rf"{fx.TEXT['T_arc'] / 1e3:.2f} kK ref.",
            transform=ax.transAxes, ha="right", va="top" if logx else "bottom",
            fontsize=6.6, color="#3a3a3a",
            linespacing=1.3)

    for ax, tag in zip(axs, "ab"):
        ax.text(0.03, 0.93, f"({tag})", transform=ax.transAxes, ha="left", va="top",
                fontsize=9)

    # line style = source, colour = quantity: two legends
    src = [Line2D([], [], color="#2b2b2b", **SIM), Line2D([], [], color="#2b2b2b", **REF)]
    lg = axs[0].legend(src, ["PlasmeX", "Dinulescu & Pfender (1980)"], frameon=False,
                       handlelength=2.6, borderaxespad=0.4, loc="upper left",
                       bbox_to_anchor=(0.0, 0.88))
    axs[0].add_artist(lg)
    qty = [Line2D([], [], color=C_E, lw=1.5), Line2D([], [], color=C_S, lw=1.5)]
    axs[0].legend(qty, [r"$n_\mathrm{e}$", r"$n_\mathrm{sa}$ (Saha at $T_\mathrm{e}$)"],
                  frameon=False, handlelength=1.6, loc="lower right")
    qty = [Line2D([], [], color=C_E, lw=1.5), Line2D([], [], color=C_H, lw=1.5)]
    axs[1].legend(qty, [r"$T_\mathrm{e}$", r"$T_\mathrm{h}$"], frameon=False,
                  handlelength=1.6, loc="upper left", bbox_to_anchor=(0.0, 0.85))
    axs[0].set_title(r"PlasmeX vs Dinulescu & Pfender (1980), Figs. 5 and 7" "\n"
                     r"Ar, $p_0$ = 100 kPa, $T_\mathrm{w}$ = 400 K, "
                     rf"$j$ = {fx.J_ACM2} A cm$^{{-2}}$", fontsize=8.2, pad=4,
                     linespacing=1.35)

    fig.align_ylabels(axs)
    for ext in ("pdf", "png"):
        fig.savefig(f"{out_stem}.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)


def report(p, s):
    def at(x, y, xq):
        return np.interp(xq, x, y)
    print(f"{'x (mm)':>7s} | {'ne paper':>9s} {'ne sim':>9s} | {'nsa p':>9s} {'nsa s':>9s} | "
          f"{'Te p':>6s} {'Te s':>6s} | {'Th p':>6s} {'Th s':>6s}")
    for xq in (0.01e-3, 0.05e-3, 0.1e-3, 0.2e-3, 0.3e-3, 0.4e-3, 0.5e-3, 0.6e-3, 0.7e-3,
               0.76e-3):
        if xq > s["L"] * 1.0001:
            continue
        print(f"{xq * 1e3:7.2f} | {at(*p['ne'], xq):9.2e} {at(s['d'], s['ne'], xq):9.2e} | "
              f"{at(*p['nsa'], xq):9.2e} {at(s['d'], s['nsa'], xq):9.2e} | "
              f"{at(*p['Te'], xq) / 1e3:6.2f} {at(s['d'], s['Te'], xq) / 1e3:6.2f} | "
              f"{at(*p['Th'], xq) / 1e3:6.2f} {at(s['d'], s['Th'], xq) / 1e3:6.2f}")


def sanity(s, name):
    """Flag a plotfile whose conditions are not the paper's -- the comparison would be moot."""
    j = np.median(np.abs(s["j"]))
    msgs = []
    if abs(j / fx.J - 1) > 0.05:
        msgs.append(f"|j| = {j:.3e} A/m^2, Dinulescu & Pfender is j = {fx.J:.2e}")
    if abs(s["Th"][0] - fx.T_WALL) > 50:
        msgs.append(f"wall T_h = {s['Th'][0]:.0f} K, the paper has T_w = {fx.T_WALL:.0f} K")
    if abs(s["Te"][-1] / fx.TEXT["T_arc"] - 1) > 0.02:
        msgs.append(f"T_e at x = L is {s['Te'][-1]:.0f} K, the paper's arc end is "
                    f"{fx.TEXT['T_arc']:.0f} K -- set prob.arc_Te")
    if abs(s["L"] / L_REF - 1) > 0.05:
        msgs.append(f"domain L = {s['L'] * 1e3:.3f} mm, the paper's layer ends at "
                    f"{L_REF * 1e3:.3f} mm with T = T_e = 12 400 K and a gradient")
    if np.ptp(s["j"]) > 0.01 * j:
        msgs.append(f"j varies by {np.ptp(s['j']) / j:.1%} across the domain (not steady?)")
    for m in msgs:
        print(f"WARNING [{name}]: {m}")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Compare PlasmeX anode plotfile(s) with Dinulescu & Pfender (1980) "
                    "Figs. 5 and 7.")
    ap.add_argument("plotfiles", nargs="+", help="PlasmeX plotfile directories (1-D, anode)")
    ap.add_argument("-o", "--outdir", default=".", help="where the figures go (default: .)")
    ap.add_argument("--logx", action="store_true",
                    help="log-log axes (as the Semenov comparison) instead of the paper's")
    ap.add_argument("--pdf", default=None, help="path to the Dinulescu & Pfender PDF (only "
                    "needed if an extracted CSV is missing or --reextract is given)")
    ap.add_argument("--reextract", action="store_true",
                    help="re-digitise the paper curves (writes the CSVs and the verification "
                         "figures)")
    a = ap.parse_args(argv)

    pdf = a.pdf and os.path.abspath(a.pdf)
    p = fx.extract(pdf) if a.reextract else fx.load(pdf)
    os.makedirs(a.outdir, exist_ok=True)

    for pf in a.plotfiles:
        pf = pf.rstrip("/")
        if not os.path.isfile(os.path.join(pf, "Header")):
            print(f"skipping {pf}: not a plotfile (no Header)")
            continue
        name = os.path.basename(os.path.abspath(pf))
        s = load_sim(pf)
        s["nsa"] = n_saha(s["Te"], s["Te"], fx.P0)
        print(f"\n=== {name}  t = {s['t']:.4e} s   j = {np.median(s['j']):.4e} A/m^2   "
              f"domain {s['L'] * 1e3:.3f} mm   vs Dinulescu & Pfender Figs. 5, 7")
        sanity(s, name)
        report(p, s)
        stem = os.path.join(a.outdir, f"compare57_{name}{'_logx' if a.logx else ''}")
        figure(p, s, stem, a.logx)
        print(f"wrote {stem}.pdf / .png")


if __name__ == "__main__":
    main()
