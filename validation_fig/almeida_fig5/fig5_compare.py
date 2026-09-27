"""
Shared engine of compare5a.py / compare5b.py: a PlasmeX 1-D cathode plotfile against
Almeida, Benilov & Naidis (2008) J. Phys. D 41 245201, Fig. 5(a) or 5(b).

One journal-style figure per plotfile, three stacked panels on a common axis:
(a) n_e, n_i   (b) T_e, T_h   (c) j_i/j.

ORIENTATION. The paper's x is the distance from the cathode surface, cathode on the LEFT.
In PlasmeX the cathode is the x_hi face (e_bc_hi = electrode, T_gas(x_hi) = wall_temp), so
the simulation is plotted against d = x_hi - x and both share the paper's orientation.

The paper curves come from fig5{a,b}_extracted.csv (fig5_extract.py, run automatically when
the CSV is missing). The simulation is read off the AMR leaf cells by bohm.sample, as gif.py
does; j_i/j = JIX / (JEX + JIX). The boundary cell holds the FACE value of JEX/JIX
(BoundaryFaceOverwrite), so the wall j_i/j is the flux the solver actually applied.
"""
import argparse
import importlib.util
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import (LogLocator, NullFormatter, MultipleLocator, FixedLocator,
                               FuncFormatter)

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fig5_extract  # noqa: E402

XLIM = (1e-8, 1e-2)

# Colour = particle population (electrons / heavy particles); line style = source.
C_E, C_H, C_J = "#2a78d6", "#eb6834", "#2b2b2b"
SIM = dict(lw=1.5, ls="-")
REF = dict(lw=1.3, ls=(0, (4.5, 2.2)))

# The run conditions every Fig. 5 panel shares; a plotfile that departs from them is flagged.
P0, T_WALL = 1.0e5, 3500.0


def load_bohm(plt_file):
    """Import bohm.py by path -- beside the plotfile first, then here, then one level up."""
    for d in (os.path.dirname(os.path.abspath(plt_file)), HERE, os.path.join(HERE, "..")):
        path = os.path.join(d, "bohm.py")
        if os.path.isfile(path):
            name = "bohm_" + re.sub(r"\W+", "_", os.path.abspath(d))
            if name in sys.modules:
                return sys.modules[name]
            spec = importlib.util.spec_from_file_location(name, path)
            mod = importlib.util.module_from_spec(spec)
            sys.modules[name] = mod
            spec.loader.exec_module(mod)
            return mod
    raise SystemExit(f"cannot find bohm.py beside {plt_file}, in {HERE} or its parent")


def load_sim(plt_file):
    b = load_bohm(plt_file)
    d = b.sample(plt_file, ("NE", "ARP", "ETEMP", "TGAS", "JEX", "JIX"), want_derived=False)
    dist = d["hi"] - d["x"]
    o = np.argsort(dist)
    j = d["JEX"] + d["JIX"]
    return {"d": dist[o], "ne": d["NE"][o], "ni": d["ARP"][o], "Te": d["ETEMP"][o],
            "Th": d["TGAS"][o], "ji_j": (d["JIX"] / j)[o], "j": j[o],
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


def _sci(v, _pos=None):
    """Symlog tick label: 0, 10^k, -10^k."""
    if v == 0:
        return "0"
    k = int(round(np.log10(abs(v))))
    return rf"${'-' if v < 0 else ''}10^{{{k}}}$"


def figure(case, p, s, out_stem):
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["STIX", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix", "font.size": 8.5, "axes.labelsize": 8.5,
        "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
    })
    fig, axs = plt.subplots(3, 1, figsize=(3.5, 6.4), sharex=True,
                            gridspec_kw=dict(height_ratios=[1.15, 1.0, 0.85], hspace=0.07))
    for ax in axs:
        style_axes(ax)
        # the region beyond the PlasmeX domain (arc boundary at d = L)
        ax.axvspan(s["L"], XLIM[1], color="#000000", alpha=0.045, lw=0)
        ax.axvline(s["L"], color="#8a8a8a", lw=0.6, ls=":")
    lab = case["labels"]

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
    ax.set_ylim(*case["T_ylim"])
    ax.yaxis.set_major_locator(MultipleLocator(case["T_major"]))
    ax.yaxis.set_minor_locator(MultipleLocator(case["T_major"] / 5))
    ax.set_ylabel(r"$T_\mathrm{e},\ T_\mathrm{h}$ (10$^3$ K)")
    ax.text(*lab["Te"], r"$T_\mathrm{e}$", color=C_E, fontsize=9, va="bottom")
    ax.text(*lab["Th"], r"$T_\mathrm{h}$", color=C_H, fontsize=9, va="bottom")

    # ---- (c) ion fraction of the current ----
    ax = axs[2]
    ax.axhline(0.0, color="#b0b0b0", lw=0.5)
    ax.plot(*p["ji_j"], color=C_J, **REF)
    ax.plot(s["d"], s["ji_j"], color=C_J, **SIM)
    if case.get("ji_symlog"):
        # The two curves differ by orders of magnitude and the paper's goes negative: a log
        # scale above linthresh, linear through zero.
        lt = case["ji_symlog"]
        ax.set_yscale("symlog", linthresh=lt, linscale=0.6)
        ax.set_ylim(*case["ji_ylim"])
        top = int(np.floor(np.log10(case["ji_ylim"][1])))
        majors = [-lt, 0.0] + [10.0 ** k for k in range(int(np.log10(lt)), top + 1)]
        ax.yaxis.set_major_locator(FixedLocator(majors))
        ax.yaxis.set_major_formatter(FuncFormatter(_sci))
        minors = [m * 10.0 ** k for k in range(int(np.log10(lt)), top + 1)
                  for m in range(2, 10) if m * 10.0 ** k < case["ji_ylim"][1]]
        ax.yaxis.set_minor_locator(FixedLocator(minors))
        ax.set_ylabel(r"$j_\mathrm{i}/j$ (symlog)")
    else:
        ax.set_ylim(*case["ji_ylim"])
        ax.yaxis.set_major_locator(MultipleLocator(case["ji_major"]))
        ax.yaxis.set_minor_locator(MultipleLocator(case["ji_major"] / 4))
        ax.set_ylabel(r"$j_\mathrm{i}/j$")
    ax.set_xlabel(r"Distance from cathode surface, $x$ (m)")
    ax.text(*case["wall_note_at"],
            rf"wall value $j_\mathrm{{iw}}/j_\mathrm{{c}}$:" "\n"
            rf"{s['ji_j'][0]:.3g} PlasmeX,  {case['table1_jiw']:.3g} ref. (Table 1)",
            transform=ax.transAxes, fontsize=6.8, color="#3a3a3a", linespacing=1.3)

    for ax, tag in zip(axs, "abc"):
        ax.text(0.975, 0.93, f"({tag})", transform=ax.transAxes, ha="right", va="top",
                fontsize=9)

    handles = [Line2D([], [], color="#2b2b2b", **SIM), Line2D([], [], color="#2b2b2b", **REF)]
    lg_ax, lg_kw = case["legend"]
    axs[lg_ax].legend(handles, ["PlasmeX", f"Almeida et al. (2008)"], frameon=False,
                      handlelength=2.6, borderaxespad=0.4, **lg_kw)
    mant, k = f"{np.median(s['j']):.2e}".split("e")      # rounds before splitting: 1.00e+06
    axs[0].set_title(rf"PlasmeX vs Almeida et al. (2008), Fig. 5({case['fig']})" "\n"
                     rf"Ar, $p_0$ = 1 bar, $T_\mathrm{{w}}$ = 3500 K, "
                     rf"$j_\mathrm{{c}}$ = {mant}$\times$10$^{{{int(k)}}}$ A m$^{{-2}}$",
                     fontsize=8.2, pad=4, linespacing=1.35)

    fig.align_ylabels(axs)
    for ext in ("pdf", "png"):
        fig.savefig(f"{out_stem}.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)


def report(p, s):
    def at(x, y, xq):
        return np.interp(np.log(xq), np.log(x), y)
    print(f"{'x (m)':>9s} | {'ne paper':>9s} {'ne sim':>9s} | {'ni paper':>9s} {'ni sim':>9s} | "
          f"{'Te p':>6s} {'Te s':>6s} | {'Th p':>6s} {'Th s':>6s} | {'ji/j p':>9s} {'ji/j s':>9s}")
    for xq in (1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 2e-4):
        if xq > s["L"] * 1.0001:
            continue
        print(f"{xq:9.1e} | {at(*p['ne'], xq):9.2e} {at(s['d'], s['ne'], xq):9.2e} | "
              f"{at(*p['ni'], xq):9.2e} {at(s['d'], s['ni'], xq):9.2e} | "
              f"{at(*p['Te'], xq) / 1e3:6.2f} {at(s['d'], s['Te'], xq) / 1e3:6.2f} | "
              f"{at(*p['Th'], xq) / 1e3:6.2f} {at(s['d'], s['Th'], xq) / 1e3:6.2f} | "
              f"{at(*p['ji_j'], xq):9.2e} {at(s['d'], s['ji_j'], xq):9.2e}")


def sanity(case, s, name):
    """Flag a plotfile whose conditions are not the figure's -- the comparison would be moot."""
    jc = np.median(s["j"])
    msgs = []
    if abs(jc / case["jc"] - 1) > 0.05:
        msgs.append(f"j = {jc:.3e} A/m^2, Fig. 5({case['fig']}) is jc = {case['jc']:.1e}")
    if abs(s["Th"][0] - T_WALL) > 100:
        msgs.append(f"wall T_h = {s['Th'][0]:.0f} K, the paper has Tw = {T_WALL:.0f} K")
    if np.ptp(s["j"]) > 0.01 * jc:
        msgs.append(f"j varies by {np.ptp(s['j']) / jc:.1%} across the domain (not steady?)")
    for m in msgs:
        print(f"WARNING [{name}]: {m}")


def main(case, argv=None):
    ap = argparse.ArgumentParser(
        description=f"Compare PlasmeX plotfile(s) with Almeida et al. (2008) Fig. 5({case['fig']}).")
    ap.add_argument("plotfiles", nargs="+", help="PlasmeX plotfile directories (1-D)")
    ap.add_argument("-o", "--outdir", default=".", help="where the figures go (default: .)")
    ap.add_argument("--pdf", default=None, help="path to the Almeida 2008 PDF (only needed "
                    "if the extracted CSV is missing or --reextract is given)")
    ap.add_argument("--reextract", action="store_true",
                    help="re-extract the paper curves from the PDF (writes the CSV and the "
                         "verification figure)")
    a = ap.parse_args(argv)

    fig = case["fig"]
    if a.reextract:
        p = fig5_extract.extract(fig, a.pdf and os.path.abspath(a.pdf))
    else:
        p = fig5_extract.load(fig, a.pdf and os.path.abspath(a.pdf))
    os.makedirs(a.outdir, exist_ok=True)

    for pf in a.plotfiles:
        pf = pf.rstrip("/")
        if not os.path.isfile(os.path.join(pf, "Header")):
            print(f"skipping {pf}: not a plotfile (no Header)")
            continue
        name = os.path.basename(os.path.abspath(pf))
        s = load_sim(pf)
        print(f"\n=== {name}  t = {s['t']:.4e} s   j = {np.median(s['j']):.4e} A/m^2   "
              f"domain {s['L'] * 1e6:.0f} um   vs Fig. 5({fig})")
        sanity(case, s, name)
        report(p, s)
        stem = os.path.join(a.outdir, f"compare5{fig}_{name}")
        figure(case, p, s, stem)
        print(f"wrote {stem}.pdf / .png")
