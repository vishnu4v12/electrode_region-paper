"""
Two published argon anode layers side by side (reference against reference):

    Dinulescu & Pfender (1980) Figs. 5, 7   j =  900 A/cm^2, p = 100 kPa, T_w =  400 K
    Semenov et al.     (2016) Fig. 3(b)     j = 1500 A/cm^2, p = 1 atm,   T_s = 1000 K

The current densities (and wall temperatures) differ, so this is not a like-for-like check:
it shows how the two models' layers differ in structure -- Pfender's quasi-neutral boundary
layer (sheath as a boundary condition, arc-end state with a gradient at x = L = 0.77 mm)
against Semenov's unified layer (sheath resolved, uniform plasma beyond ~1-2 mm).

    python compare_pfender_semenov3b.py [plotfile] [-o outdir]

compare_pfender_semenov3b.{pdf,png}: rows (a) n_e, n_sa (b) T_e, T_h; columns: the Pfender
figures' linear axes (x in mm) and the Semenov figure's log axes (x in m, sheath visible).
An optional PlasmeX plotfile is drawn over both (solid), with the name as a suffix.

n_sa is each paper's own: Pfender's n_e' is Saha at T_e (one temperature); Semenov's n_sa is
Saha at T_e with Dalton's law at T_e, T_h (two temperatures). A PlasmeX n_sa is drawn with
Semenov's definition. The curves come from fig{5,7}_extracted.csv (fig57_extract.py) and
../semenov_fig3/fig3b_extracted.csv (fig3_extract.py).
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
import fig3_extract  # noqa: E402
from fig3_compare import n_saha, C_S  # noqa: E402
from fig5_compare import load_sim, style_axes, XLIM, SIM, REF, C_E, C_H  # noqa: E402

SEM = REF                                            # Semenov: dashed, as in compare3b
PF = dict(lw=1.3, ls=(0, (6.0, 1.6, 1.3, 1.6)))      # Pfender: dash-dot
X_LIN = 2.5                                          # mm; Semenov uniform beyond ~1-2 mm


def style_lin(ax):
    ax.set_xlim(0, X_LIN)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.tick_params(which="major", length=3.5, width=0.6)
    ax.tick_params(which="minor", length=2.0, width=0.5)
    for sp in ax.spines.values():
        sp.set_linewidth(0.6)


def figure(pf, sm, s, out_stem):
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["STIX", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix", "font.size": 8.5, "axes.labelsize": 8.5,
        "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
    })
    fig, axs = plt.subplots(2, 2, figsize=(7.0, 5.0), sharex="col",
                            gridspec_kw=dict(hspace=0.07, wspace=0.26))
    sources = [(pf, PF), (sm, SEM)] + ([(s, SIM)] if s else [])

    for col, logx in enumerate((False, True)):
        xs, nu = (1.0, 1.0) if logx else (1e3, 1e22)
        for ax in axs[:, col]:
            style_axes(ax) if logx else style_lin(ax)
        # ---- densities ----
        ax = axs[0, col]
        for d, st in sources:
            ax.plot(d["ne"][0] * xs, d["ne"][1] / nu, color=C_E, **st)
            ax.plot(d["nsa"][0] * xs, d["nsa"][1] / nu, color=C_S, **st)
        if logx:
            ax.set_yscale("log")
            ax.set_ylim(1e19, 1e24)
            ax.yaxis.set_minor_formatter(NullFormatter())
            ax.set_ylabel(r"$n_\mathrm{e},\ n_\mathrm{sa}$ (m$^{-3}$)")
        else:
            top = max(np.max(d["ne"][1][d["ne"][0] * xs <= X_LIN]) for d, _ in sources) / nu
            ax.set_ylim(0, 2 * np.ceil(1.05 * top / 2))
            ax.yaxis.set_major_locator(MultipleLocator(2))
            ax.yaxis.set_minor_locator(MultipleLocator(0.5))
            ax.set_ylabel(r"$n_\mathrm{e},\ n_\mathrm{sa}$ (10$^{22}$ m$^{-3}$)")
        # ---- temperatures ----
        ax = axs[1, col]
        for d, st in sources:
            ax.plot(d["Te"][0] * xs, d["Te"][1] / 1e3, color=C_E, **st)
            ax.plot(d["Th"][0] * xs, d["Th"][1] / 1e3, color=C_H, **st)
        ax.set_ylim(0, 18)
        ax.yaxis.set_major_locator(MultipleLocator(3))
        ax.yaxis.set_minor_locator(MultipleLocator(1))
        ax.set_ylabel(r"$T_\mathrm{e},\ T_\mathrm{h}$ (10$^3$ K)")
        ax.set_xlabel(r"Distance from anode surface, $x$ " + ("(m)" if logx else "(mm)"))

    for ax, tag in zip(axs.flat, "abcd"):
        ax.text(0.03, 0.95, f"({tag})", transform=ax.transAxes, ha="left", va="top",
                fontsize=9)
    axs[0, 0].set_title("linear axes (as Dinulescu & Pfender)", fontsize=8.2, pad=3)
    axs[0, 1].set_title("log axes (as Semenov et al.)", fontsize=8.2, pad=3)

    # colour = quantity (in-panel legends), line style = source (figure legend)
    axs[0, 0].legend([Line2D([], [], color=C_E, lw=1.5), Line2D([], [], color=C_S, lw=1.5)],
                     [r"$n_\mathrm{e}$", r"$n_\mathrm{sa}$ (Saha)"], frameon=False,
                     handlelength=1.6, loc="lower right")
    axs[1, 0].legend([Line2D([], [], color=C_E, lw=1.5), Line2D([], [], color=C_H, lw=1.5)],
                     [r"$T_\mathrm{e}$", r"$T_\mathrm{h}$"], frameon=False,
                     handlelength=1.6, loc="lower right")
    handles = [Line2D([], [], color="#2b2b2b", **PF), Line2D([], [], color="#2b2b2b", **SEM)]
    labels = [r"Dinulescu & Pfender (1980), Figs. 5, 7: $j$ = 900 A cm$^{-2}$, "
              r"100 kPa, $T_\mathrm{w}$ = 400 K",
              r"Semenov et al. (2016), Fig. 3(b): $j$ = 1500 A cm$^{-2}$, "
              r"1 atm, $T_\mathrm{s}$ = 1000 K"]
    if s:
        handles.append(Line2D([], [], color="#2b2b2b", **SIM))
        labels.append(rf"PlasmeX {s['name']}: $|j|$ = {np.median(np.abs(s['j'])) / 1e4:.0f} "
                      rf"A cm$^{{-2}}$, $T_\mathrm{{h}}$(wall) = {s['Th'][1][0]:.0f} K")
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.93),
               ncol=1, frameon=False, handlelength=3.0)
    fig.align_ylabels(axs[:, 0])
    fig.align_ylabels(axs[:, 1])
    for ext in ("pdf", "png"):
        fig.savefig(f"{out_stem}.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)


def summary(name, d):
    """Key numbers of one layer: T_e at the wall and at the far end, where T_h catches up
    with T_e (|T_e - T_h| < 2 % of T_e), and n_e at 10 and 100 um."""
    xe, te = d["Te"]
    th = np.interp(xe, *d["Th"])
    neq = np.abs(te - th) > 0.02 * te
    x_eq = xe[np.flatnonzero(neq)[-1] + 1] if neq.any() and neq[-1] == 0 else np.nan
    ne = lambda xq: np.interp(xq, *d["ne"]) if xq >= d["ne"][0][0] else np.nan
    print(f"{name:32s} {te[0] / 1e3:6.2f} {te[-1] / 1e3:6.2f} {x_eq * 1e3:8.3f} "
          f"{ne(1e-5):9.2e} {ne(1e-4):9.2e}   (first point x = {xe[0] * 1e6:.2g} um)")


def main(argv=None):
    ap = argparse.ArgumentParser(description="Dinulescu & Pfender (1980) Figs. 5/7 against "
                                 "Semenov et al. (2016) Fig. 3(b), optionally with PlasmeX.")
    ap.add_argument("plotfile", nargs="?", help="optional PlasmeX plotfile (1-D, anode)")
    ap.add_argument("-o", "--outdir", default=HERE, help="where the figure goes "
                    "(default: this folder)")
    a = ap.parse_args(argv)

    pf = fx.load()
    sm = fig3_extract.load("b")
    stem = os.path.join(a.outdir, "compare_pfender_semenov3b")
    s = None
    if a.plotfile:
        s = load_sim(a.plotfile.rstrip("/"))
        s["name"] = os.path.basename(os.path.abspath(a.plotfile.rstrip("/")))
        s["nsa"] = n_saha(s["Te"], s["Th"])
        for q in ("ne", "nsa", "Te", "Th"):
            s[q] = (s["d"], s[q])
        stem += f"_{s['name']}"

    print(f"{'':32s} {'Te(0)':>6s} {'Te(L)':>6s} {'x_eq/mm':>8s} {'ne(10um)':>9s} "
          f"{'ne(100um)':>9s}   [T in kK; x_eq: T_h within 2 % of T_e beyond]")
    summary("Dinulescu & Pfender, 900 A/cm2", pf)
    summary("Semenov Fig. 3(b), 1500 A/cm2", sm)
    if s:
        summary(f"PlasmeX {s['name']}", s)
    os.makedirs(a.outdir, exist_ok=True)
    figure(pf, sm, s, stem)
    print(f"wrote {stem}.pdf / .png")


if __name__ == "__main__":
    main()
