"""
PlasmeX 1-D ANODE plotfile(s) vs Semenov, Krivtsun & Reisgen (2016) J. Phys. D 49 105204 at
j = 1500 A/cm^2 = 1.5e7 A/m^2 (Ar, 1 atm, anode surface T_s = 1000 K, unified model):
Fig. 3(b), and the 1500 A/cm^2 curves of Fig. 4 (j_i/j) and Fig. 5 (E), in one figure.

    python compare.py plt73000 [plt...] [-o outdir]

One figure per plotfile, compareJ1.5e7_<plotfile>.{pdf,png}, three stacked panels against the
distance from the anode:
    (a) n_e, n_i, n_sa                         Fig. 3(b)
    (b) T_e, T_h                               Fig. 3(b), Table 2
    (c) j_i/j (left axis)  Fig. 4   and   E (right axis)  Fig. 5, one shared zero line

ORIENTATION AND SIGNS. The paper's x is the distance from the anode, anode on the LEFT, and
j > 0. The PlasmeX anode is the x_hi face, so the simulation is plotted against d = x_hi - x.
The distance axis is MIRRORED (d increasing to the left), so the anode sits on the RIGHT, as in
the results figures (cathodes on the left, anodes on the right). Both quantities map without any further sign choice: E_paper = -EFLDX, and j_i/j =
JIX / (JEX + JIX) (flipping x flips both currents, so the ratio is unchanged; ions flowing into
the anode give j_i/j < 0 in both). Near the anode both are negative, so they are plotted with
their own signs and share the zero line.

n_sa is "the charged particle number density calculated by means of the Saha equation and
quasi-neutrality condition" (the paper's caption). For PlasmeX it is computed here from the
local T_e, T_h at p0 = 1 atm with a two-temperature Saha law (ionization governed by T_e)
and Dalton's law p = n_e k T_e + (n_i + n_a) k T_h, with the Ar+ internal partition function
4 + 2 exp(-2062 K / T_e) and Z_Ar = 1. Applied to the paper's OWN T_e, T_h curves this
reproduces the paper's n_sa to 1-4 % over the whole layer (a constant g+/g1 = 6 is 1-9 %
high), so the two n_sa curves are defined the same way.

The paper curves are read from paper_data/fig3b_extracted.csv and fig45_j1500_extracted.csv
(written by paper_data/extract.py). The simulation is read off the AMR leaf cells by
bohm.sample, loaded from the live PlasmeX tree (../plasmex.py).
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
XLIM = (1e-8, 1e-2)

# Colour = particle population (electrons / heavy particles); line style = source.
C_E, C_H, C_J, C_F, C_S = "#2a78d6", "#eb6834", "#2b2b2b", "#1a9e77", "#3f9a3f"
SIM = dict(lw=1.5, ls="-")
REF = dict(lw=1.3, ls=(0, (4.5, 2.2)))

# The paper's run conditions; a plotfile that departs from them is flagged.
J, J_ACM2, P0, T_WALL = 1.5e7, 1500, 101325.0, 1000.0
TE_WALL, TE_BULK = 5.59e3, 16.82e3           # Table 2, 1500 A/cm^2: T_e at the anode, undisturbed
L_MIN = 2e-3                                 # uniform plasma from x ~ 1-2 mm (Fig. 3b)

N_YLIM = (1e19, 1e24)
T_YLIM, T_MAJOR = (0, 19), 4                 # kK
JI_YLIM, JI_MAJOR = (-0.025, 0.005), 0.005
E_UNIT, E_BOTTOM, E_MAJOR = 1e6, -2.0, 0.5   # E axis; its zero is j_i/j's zero
LABELS = dict(ne=(1.8e-8, 4.0e20), ni=(1.8e-8, 1.4e21), nsa=(1.6e-6, 1.3e20),
              Te=(1.8e-8, 5.9), Th=(1.8e-8, 1.3), ji=(1.8e-8, -0.0047), E=(1.8e-8, -1.2))
SIM_LABEL, REF_LABEL = "Present model", "Semenov et al (2016)"

KB, ME, HP, QE = 1.380649e-23, 9.1093837015e-31, 6.62607015e-34, 1.602176634e-19
I_AR = 15.7596                               # [eV]


def n_saha(Te, Th, p=P0):
    """Charged-particle density from Saha (at T_e) + quasi-neutrality + Dalton (see doc)."""
    Te, Th = np.asarray(Te, float), np.asarray(Th, float)
    g = 4.0 + 2.0 * np.exp(-2062.0 / Te)
    S = 2.0 * g * (2 * np.pi * ME * KB * Te / HP**2) ** 1.5 * np.exp(-I_AR * QE / (KB * Te))
    a, c = (Te + Th) / Th, p / (KB * Th)        # n^2 = S (c - a n)
    return 0.5 * (-S * a + np.sqrt((S * a) ** 2 + 4 * S * c))


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
    d = b.sample(plt_file, ("NE", "ARP", "ETEMP", "TGAS", "JEX", "JIX", "EFLDX"),
                 want_derived=False)
    dist = d["hi"] - d["x"]
    o = np.argsort(dist)
    j = d["JEX"] + d["JIX"]
    s = {"d": dist[o], "ne": d["NE"][o], "ni": d["ARP"][o], "Te": d["ETEMP"][o],
         "Th": d["TGAS"][o], "ji_j": (d["JIX"] / j)[o], "j": j[o],
         "E": -d["EFLDX"][o],                                   # paper orientation
         "L": d["hi"] - d["lo"], "t": float(d["ds"].current_time)}
    s["nsa"] = n_saha(s["Te"], s["Th"])
    return s


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
    axs[0].invert_xaxis()                     # anode on the right (sharex: all three panels)
    for ax in axs:
        if s["L"] < XLIM[1]:
            # the region beyond the PlasmeX domain (arc boundary at d = L)
            ax.axvspan(s["L"], XLIM[1], color="#000000", alpha=0.045, lw=0)
            ax.axvline(s["L"], color="#8a8a8a", lw=0.6, ls=":")
    lab = LABELS

    # ---- (a) densities ----
    ax = axs[0]
    for k, c in (("ne", C_E), ("ni", C_H), ("nsa", C_S)):
        ax.plot(*p[k], color=c, **REF)
        ax.plot(s["d"], s[k], color=c, **SIM)
    ax.set_yscale("log")
    ax.set_ylim(*N_YLIM)
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_ylabel(r"$n_\mathrm{e},\ n_\mathrm{i},\ n_\mathrm{sa}$ (m$^{-3}$)")
    ax.text(*lab["ne"], r"$n_\mathrm{e}$", color=C_E, fontsize=9, ha="right")
    ax.text(*lab["ni"], r"$n_\mathrm{i}$", color=C_H, fontsize=9, ha="right")
    ax.text(*lab["nsa"], r"$n_\mathrm{sa}$", color=C_S, fontsize=9, ha="right")
    if s["L"] < XLIM[1] / 3:
        ax.text(s["L"] * 1.25, N_YLIM[0] * 2, r"$x > L$", fontsize=7, color="#6a6a6a",
                va="bottom", ha="right")

    # ---- (b) temperatures ----
    ax = axs[1]
    ax.plot(p["Te"][0], p["Te"][1] / 1e3, color=C_E, **REF)
    ax.plot(p["Th"][0], p["Th"][1] / 1e3, color=C_H, **REF)
    ax.plot(s["d"], s["Te"] / 1e3, color=C_E, **SIM)
    ax.plot(s["d"], s["Th"] / 1e3, color=C_H, **SIM)
    ax.set_ylim(*T_YLIM)
    ax.yaxis.set_major_locator(MultipleLocator(T_MAJOR))
    ax.yaxis.set_minor_locator(MultipleLocator(T_MAJOR / 4))
    ax.set_ylabel(r"$T_\mathrm{e},\ T_\mathrm{h}$ (10$^3$ K)")
    ax.text(*lab["Te"], r"$T_\mathrm{e}$", color=C_E, fontsize=9, va="bottom", ha="right")
    ax.text(*lab["Th"], r"$T_\mathrm{h}$", color=C_H, fontsize=9, va="bottom", ha="right")
    handles = [Line2D([], [], color="#2b2b2b", **SIM), Line2D([], [], color="#2b2b2b", **REF)]
    ax.legend(handles, [SIM_LABEL, REF_LABEL], frameon=False, handlelength=2.6,
              borderaxespad=0.6, loc="lower left")

    # ---- (c) j_i/j (left) and E (right), one zero line ----
    ax = axs[2]
    ax.tick_params(which="both", right=False)
    ax.axhline(0.0, color="#b0b0b0", lw=0.5)
    ax.plot(*p["ji_j"], color=C_J, **REF)
    ax.plot(s["d"], s["ji_j"], color=C_J, **SIM)
    ji_lo, ji_hi = JI_YLIM
    ax.set_ylim(ji_lo, ji_hi)
    ax.yaxis.set_major_locator(MultipleLocator(JI_MAJOR))
    ax.yaxis.set_minor_locator(MultipleLocator(JI_MAJOR / 5))
    ax.set_ylabel(r"$j_\mathrm{i}/j$")
    ax.set_xlabel(r"Distance from anode surface, $x$ (m)")
    ax.text(*lab["ji"], r"$j_\mathrm{i}/j$", color=C_J, fontsize=9, va="bottom", ha="right")

    fx = ax.twinx()
    fx.plot(p["E"][0], p["E"][1] / E_UNIT, color=C_F, **REF)
    fx.plot(s["d"], s["E"] / E_UNIT, color=C_F, **SIM)
    fx.set_ylim(E_BOTTOM, E_BOTTOM * ji_hi / ji_lo)
    fx.yaxis.set_major_locator(MultipleLocator(E_MAJOR))
    fx.yaxis.set_minor_locator(MultipleLocator(E_MAJOR / 5))
    fx.tick_params(which="both", direction="in", colors=C_F)
    fx.tick_params(which="major", length=3.5, width=0.6)
    fx.tick_params(which="minor", length=2.0, width=0.5)
    fx.spines["right"].set_color(C_F)
    fx.set_ylabel(rf"$E$ (10$^{{{int(round(np.log10(E_UNIT)))}}}$ V m$^{{-1}}$)", color=C_F)
    fx.text(*lab["E"], r"$E$", color=C_F, fontsize=9, va="bottom", ha="right")

    # the anode surface: hatched strip + thick frame line on the right of every panel
    for i, ax in enumerate(axs):
        mark_electrode(ax, "right", label="anode" if i == 0 else None,
                       spine_axes=(fx,) if i == 2 else ())
    for ax, tag in zip(axs, "abc"):
        ax.text(0.955, 0.93 if tag != "c" else 0.07, f"({tag})", transform=ax.transAxes,
                ha="right", va="top" if tag != "c" else "bottom", fontsize=9)
    # no title: the run conditions and the source figures go in the paper's caption
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

    u = E_UNIT
    print(f"{'x (m)':>9s} | {'ji/j p':>9s} {'ji/j s':>9s} | {'E p':>8s} {'E s':>8s}   "
          f"(E in 1e{int(round(np.log10(u)))} V/m, paper orientation)")
    for xq in (1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 3e-4, 1e-3, 3e-3):
        if xq > s["L"] * 1.0001:
            continue
        e_p = at(*p["E"], xq) / u if xq <= p["E"][0].max() else np.nan
        print(f"{xq:9.1e} | {at(*p['ji_j'], xq):9.5f} {at(s['d'], s['ji_j'], xq):9.5f} | "
              f"{e_p:8.4f} {at(s['d'], s['E'], xq) / u:8.4f}")
    i, k = np.argmin(p["ji_j"][1]), np.argmin(s["ji_j"])
    print(f"T_e at the anode: PlasmeX {s['Te'][0]:.0f} K, Table 2 {TE_WALL:.0f} K")
    print(f"j_i/j minimum: paper {p['ji_j'][1][i]:.5f} at x = {p['ji_j'][0][i]:.2e} m;  "
          f"PlasmeX {s['ji_j'][k]:.5f} at x = {s['d'][k]:.2e} m")


def sanity(s, name):
    """Flag a plotfile whose conditions are not the paper's -- the comparison would be moot."""
    j = np.median(np.abs(s["j"]))
    msgs = []
    if abs(j / J - 1) > 0.05:
        msgs.append(f"|j| = {j:.3e} A/m^2, Fig. 3(b) is j = {J:.2e}")
    if abs(s["Th"][0] - T_WALL) > 50:
        msgs.append(f"wall T_h = {s['Th'][0]:.0f} K, the paper has T_s = {T_WALL:.0f} K")
    if abs(s["Te"][-1] / TE_BULK - 1) > 0.02:
        msgs.append(f"T_e at x = L is {s['Te'][-1]:.0f} K, the paper's undisturbed plasma is "
                    f"{TE_BULK:.0f} K (Table 2) -- check prob.arc_Te")
    if s["L"] < L_MIN:
        msgs.append(f"domain L = {s['L'] * 1e3:.2f} mm: the paper's plasma is uniform only "
                    f"beyond ~{L_MIN * 1e3:.0f} mm, so the arc-side Dirichlet state "
                    f"sits inside the layer")
    if np.ptp(s["j"]) > 0.01 * j:
        msgs.append(f"j varies by {np.ptp(s['j']) / j:.1%} across the domain (not steady?)")
    for m in msgs:
        print(f"WARNING [{name}]: {m}")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Compare PlasmeX anode plotfile(s) with Semenov et al. (2016) at "
                    "j = 1500 A/cm^2: Fig. 3(b), Fig. 4 (j_i/j) and Fig. 5 (E).")
    ap.add_argument("plotfiles", nargs="+", help="PlasmeX plotfile directories (1-D, anode)")
    ap.add_argument("-o", "--outdir", default=".", help="where the figures go (default: .)")
    a = ap.parse_args(argv)

    p = {**load_csv("fig3b_extracted.csv"), **load_csv("fig45_j1500_extracted.csv")}
    os.makedirs(a.outdir, exist_ok=True)

    for pf in a.plotfiles:
        pf = pf.rstrip("/")
        if not os.path.isfile(os.path.join(pf, "Header")):
            print(f"skipping {pf}: not a plotfile (no Header)")
            continue
        name = os.path.basename(os.path.abspath(pf))
        s = load_sim(pf)
        print(f"\n=== {name}  t = {s['t']:.4e} s   j = {np.median(s['j']):.4e} A/m^2   "
              f"domain {s['L'] * 1e3:.3f} mm   vs Semenov j = 1500 A/cm^2")
        sanity(s, name)
        report(p, s)
        stem = os.path.join(a.outdir, f"compareJ1.5e7_{name}")
        figure(p, s, stem)
        print(f"wrote {stem}.pdf / .png")


if __name__ == "__main__":
    main()
