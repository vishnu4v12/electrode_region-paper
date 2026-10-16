"""
Results-section figures for the two base cases (PlasmeX runs/base_case_anode, base_case_cathode).

    python base_case.py <anode plotfile> <cathode plotfile> [-o outdir]

Writes, each as .pdf and .png:
    base_qwall        the terms of the electrode heat flux q_w, main.tex eq (24) (eq:bc-qwall),
                      the code's own values from each plotfile's plasmex_info
    base_profiles     n_e, n_i, n_n / T_e, T_h, T_eff / phi / j_i/j against the distance from
                      the electrode, cathode (left column) and anode (right column)
    base_een_budget   every steady term of the electron energy equation, main.tex eq (11)
                      (eq:eenergy), against the distance from the electrode: the divergence of
                      the whole electron energy flux h_e, eq (12), Joule heating, elastic
                      exchange, ionization, radiation and radial loss (drawn only when nonzero)
DISTANCE AXIS. The axis coordinate is ln x + x/X_C (X_C = 50 um): logarithmic near the
electrode, so the sheath is resolved, and smoothly becoming linear towards the arc end (90 %
linear at 500 um), so the profile there is seen undistorted. Ticks in um.

ORIENTATION. Cathode panels on the left, anode panels on the right. x is the distance from the
electrode on a log axis. The cathode panels have the cathode on the LEFT (x increasing to the
right); the anode panels are mirrored, with the anode on the RIGHT (x increasing to the left).
The electrode edge of every panel is marked by electrode_mark.mark_electrode (hatched strip,
thick frame line, a label on the top panel); panel tags go on the side away from the electrode.
T_eff (the Wannier ion temperature) is drawn only where n_n > n_i: the closure is ion-atom only.
The script also prints a table of the electrode quantities (T_wall, sheath edge and layer voltage, current
partition, the eq (24) heat-flux terms) for the paper's table, and the maximum ion impact energy (the mean
energy at the sheath edge plus e phi_s, for an ion crossing the sheath without a collision).

The plotfiles are read off the AMR leaf cells by bohm.sample, and every budget term comes from
bohm.electron_energy_flux / electron_energy_budget, loaded from the live PlasmeX tree
(../validation_fig/plasmex.py). Nothing is rebuilt here that bohm already has.

DIVERGENCE OF h_e. The budget plots bohm's NDIVH, the total -div h_e on each leaf cell's own
two faces (the wall face taking the BC's q_wall). It is NaN in the arc cell, whose face flux the
plotfile does not carry, and in the cells next to a coarse-fine interface, where the two-cell
grad T_e is not the operator's closure gradient.
"""
import argparse
import os
import sys

import numpy as np
from scipy.special import lambertw
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import (FixedLocator, FuncFormatter, NullFormatter,
                               SymmetricalLogLocator)

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "validation_fig"))
from plasmex import EXE, exe_module  # noqa: E402
from electrode_mark import mark_electrode  # noqa: E402

FIELDS = ("NE", "ARP", "AR", "PHI", "ETEMP", "TGAS", "EFLDX", "JEX", "JIX", "UNX", "HVELX",
          "EENRXN", "EJH", "ERXN")
XLIM = (1e-9, 5e-4)
# Log-linear distance axis: the axis coordinate is u = ln(x/x0) + (x - x0)/X_C. Its slope
# 1/x + 1/X_C makes it logarithmic for x << X_C (the sheath) and linear for x >> X_C (the arc end),
# with a smooth change in between: no joint, so no corner or curvature jump in any curve.
X_C = 5e-5
K_B, E_CH = 1.380649e-23, 1.602176634e-19
M_AR = 39.948 * 1.66053906660e-27

# Population colours as in the validation figures (electrons blue, heavy particles orange).
C_E, C_I, C_N, C_X = "#2a78d6", "#eb6834", "#1baf7a", "#2b2b2b"

# eq (11) terms: key -> (label, colour, line style). -div h_e is dashed, the volume sources
# solid. Colours: the validated categorical order of the dataviz palette.
BUDGET = [
    ("NDIVH",   r"$-\nabla\!\cdot\!\mathbf{h}_e$",                   "#2a78d6", "--"),
    ("JOULE",   r"$-e\Gamma_e\!\cdot\!E$",                           "#eda100", "-"),
    ("ELASTIC", "elastic exchange",                                  "#e87ba4", "-"),
    ("ION",     r"$-\varepsilon_{\rm ion}S$",                        "#008300", "-"),
    ("RAD",     r"$-w_{\rm rad}$",                                   "#4a3aa7", "-"),
    ("RADL",    r"$-Q_e$",                                           "#e34948", "-"),
]


def _fwd(x):
    x = np.clip(np.asarray(x, dtype=float), XLIM[0] * 1e-3, None)
    return np.log(x / XLIM[0]) + (x - XLIM[0]) / X_C


def _inv(u):
    # ln x + x/X_C = v  ->  x = X_C W(e^v / X_C), W the Lambert function (principal branch)
    v = np.clip(np.asarray(u, dtype=float), -50.0, 600.0) + np.log(XLIM[0]) + XLIM[0] / X_C
    return X_C * np.real(lambertw(np.exp(v) / X_C))


def sci_tex(v):
    """Mathtext for a positive value: 10^{n} when the mantissa is 1, else m\\times10^{n} (9e6, 1.5e8)."""
    e = int(np.floor(np.log10(v) + 1e-12))
    m = v / 10.0 ** e
    return rf"10^{{{e}}}" if abs(m - 1) < 1e-9 else rf"{m:g}\times10^{{{e}}}"


def _um_label(x, _pos=None):
    um = x * 1e6
    if um < 1.0:
        return rf"$10^{{{int(round(np.log10(um)))}}}$"
    return f"{um:.0f}"


def split_axis(ax):
    """The log-linear distance axis (u = ln x + x/X_C), ticks in um."""
    ax.set_xscale("function", functions=(_fwd, _inv))
    ax.set_xlim(*XLIM)
    dec = 10.0 ** np.arange(np.ceil(np.log10(XLIM[0])), -4.5)          # 1 nm ... 10 um
    lin = np.arange(1e-4, XLIM[1] * (1 + 1e-9), 1e-4)                   # 100 ... 500 um
    ax.xaxis.set_major_locator(FixedLocator(np.unique(np.r_[dec, lin])))
    minor = [m * d for d in dec for m in range(2, 10)]                    # up to 90 um
    minor += list(np.arange(1.5e-4, XLIM[1], 1e-4))
    ax.xaxis.set_minor_locator(FixedLocator(minor))
    ax.xaxis.set_major_formatter(FuncFormatter(_um_label))
    ax.xaxis.set_minor_formatter(NullFormatter())


def style(ax, xlog=True):
    if xlog:
        split_axis(ax)
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.tick_params(which="major", length=3.5, width=0.6)
    ax.tick_params(which="minor", length=2.0, width=0.5)
    for sp in ax.spines.values():
        sp.set_linewidth(0.6)


def orient(ax, tag):
    """Anode on the right: mirror the anode panels' distance axis (cathode stays on the left)."""
    if tag == "Anode":
        ax.invert_xaxis()


def side(tag):
    """Screen side of the electrode: 'left' for the cathode, 'right' for the (mirrored) anode."""
    return "right" if tag == "Anode" else "left"


def away(tag):
    return "left" if side(tag) == "right" else "right"


def xlabel(tag):
    return rf"Distance from the {tag.lower()}, $x$ ($\mu$m)"


def rc():
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["STIX", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix", "font.size": 8.5, "axes.labelsize": 8.5,
        "legend.fontsize": 7.0, "xtick.labelsize": 8, "ytick.labelsize": 8,
    })


def load(pf):
    """Leaf-cell state, the eq (11)/(12) terms and the wall quantities, ordered by distance."""
    b = exe_module("bohm")
    d = b.sample(pf, FIELDS, want_derived=False)
    if d["missing"]:
        raise SystemExit(f"{pf}: missing {d['missing']}")
    knobs, _ = b.read_knobs(d["dir"], None, pf)
    q = b.electron_energy_flux(d, knobs, EXE)
    if q is None:
        raise SystemExit(f"{pf}: cannot rebuild h_e (transport tables in the run directory?)")
    bud = b.electron_energy_budget(d, knobs, q, EXE)
    parts = {k: np.asarray(bud[k], dtype=float)
             for k in ("NDIVH", "JOULE", "ELASTIC", "ION", "RAD", "RADL", "DEDT")}

    t_eff = b.ion_transport(d["TGAS"], d["AR"], np.abs(d["EFLDX"]))[0]
    dist = d["hi"] - d["x"]
    o = np.argsort(dist)
    j = d["JEX"] + d["JIX"]
    s = {"d": dist[o], "ne": d["NE"][o], "ni": d["ARP"][o], "nn": d["AR"][o],
         "Te": d["ETEMP"][o], "Th": d["TGAS"][o], "Teff": t_eff[o], "phi": d["PHI"][o],
         "E": d["EFLDX"][o], "ji_j": (d["JIX"] / j)[o], "j": j[o], "dxc": d["dxc"][o],
         "budget": {k: v[o] for k, v in parts.items()}, "rxn_err": bud["RXN_ERR"],
         "info": b.read_run_info(pf), "t": float(d["ds"].current_time),
         "L": d["hi"] - d["lo"], "name": os.path.basename(os.path.abspath(pf))}
    return s


def sheath_edge(d, ne, ni, phi):
    """Sheath edge x_s and its potential phi_s: where |n_i - n_e|/n_i falls to 1 % for the last time,
    log-linearly interpolated between the outermost cell above 1 % and the next one (the cells there
    are 15 nm at the cathode and 61 nm at the anode, too coarse to quote the cell centre). Also
    returns that outermost cell's index. (The anode has an electron-rich layer at the wall and an
    ion-rich one out to ~1 um; the first quasi-neutral cell from the wall falls between the two.)"""
    sep = np.abs(ni - ne) / ni
    i = int(np.nonzero(sep > 0.01)[0].max())
    f = (np.log(0.01) - np.log(sep[i])) / (np.log(sep[i + 1]) - np.log(sep[i]))
    return d[i] + f * (d[i + 1] - d[i]), phi[i] + f * (phi[i + 1] - phi[i]), i


P0 = 101325.0   # [Pa], prob.pres


def n_eq_Th(Th, p=P0):
    """Equilibrium charged-particle density n_eq(T_h) of the model's OWN reactions (prob.rxn_rates =
    almeida, eq S: k_ion n_e n_n = k_rec n_e^2 n_i) at T_e = T_h = T, with quasi-neutrality and Dalton's
    law p = (2 n + n_n) k T: exactly the function that writes exe/ar_eq.dat (gen_arc_eq_table.n_eq), i.e.
    the arc-face density at x = L. Below ~3 kK n_eq is below any plotted density; NaN where the root
    finder's bracket (n >= 1e-30 p/kT) holds no root."""
    g = exe_module("gen_arc_eq_table")
    out = np.full(np.shape(Th), np.nan)
    for k, T in enumerate(np.asarray(Th, float)):
        try:
            out[k] = g.n_eq(float(T), p, None, None, "almeida")
        except ValueError:
            pass
    return out


def fig_profiles(runs, stem):
    rc()
    fig, axs = plt.subplots(4, 2, figsize=(7.0, 7.6), sharex="col",
                            gridspec_kw=dict(height_ratios=[1.15, 1, 0.9, 0.9],
                                             hspace=0.17, wspace=0.28))
    for col, (tag, s) in enumerate(runs):
        d = s["d"]
        ax = axs[0, col]
        ax.plot(d, s["ne"], color=C_E, lw=1.4, label=r"$n_e$")
        ax.plot(d, s["ni"], color=C_I, lw=1.4, ls=(0, (4, 2)), label=r"$n_i$")
        neq = n_eq_Th(s["Th"])
        ax.plot(d, neq, color=C_N, lw=1.2, ls="-.", label=r"$n_{\rm eq}(T_h)$")
        ax.set_yscale("log")
        lo = 0.3 * min(s["ne"].min(), s["ni"].min())
        ax.set_ylim(lo, 3 * max(np.nanmax(neq), s["ne"].max()))
        ax.set_ylabel(r"Number density (m$^{-3}$)")
        ax.set_title(tag, fontsize=9, pad=10)
        ax.legend(frameon=False, loc="best", handlelength=2.2)

        ax = axs[1, col]
        ax.plot(d, s["Te"] / 1e3, color=C_E, lw=1.4, label=r"$T_e$")
        ax.plot(d, s["Th"] / 1e3, color=C_I, lw=1.4, label=r"$T_h$")
        ax.set_ylim(0, 1.08 * max(s["Te"].max(), s["Th"].max()) / 1e3)
        ax.set_ylabel(r"Temperature (10$^3$ K)")
        ax.legend(frameon=False, loc="best", handlelength=2.2)

        ax = axs[2, col]
        ax.axhline(0, color="#b0b0b0", lw=0.5)
        ax.plot(d, s["phi"], color=C_X, lw=1.4)
        ax.set_ylabel(r"$\phi$ (V)")

        ax = axs[3, col]
        ax.axhline(0, color="#b0b0b0", lw=0.5)
        ax.plot(d, s["ji_j"], color=C_X, lw=1.4)
        ax.set_ylabel(r"$j_i/j$")
        ax.set_xlabel(xlabel(tag))
    for (r, c), ax in np.ndenumerate(axs):
        style(ax)
        # tags just above each frame, so they never sit on a curve or a legend
        tag = runs[c][0]
        mark_electrode(ax, side(tag), label=tag.lower() if r == 0 else None)
        ax.text(0.0 if away(tag) == "left" else 1.0, 1.015, f"({'abcdefgh'[2 * r + c]})",
                transform=ax.transAxes, fontsize=9, va="bottom", ha=away(tag))
    for col, (tag, _) in enumerate(runs):
        orient(axs[0, col], tag)                  # sharex="col": the whole column follows
    fig.align_ylabels(axs[:, 0])
    fig.align_ylabels(axs[:, 1])
    for ext in ("pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)


# Linear band of the budget's symlog axis, per panel. The cathode's -div h_e is the small residual
# of a column balanced by the radial loss and changes sign three times (~20, ~70, ~496 um); with the
# default band (1e9) every crossing is drawn as a vertical wall through several decades.
LINTHRESH = {"Cathode": 1e11}


def fig_budget(runs, stem):
    rc()
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.6), gridspec_kw=dict(wspace=0.25))
    for ax, (tag, s) in zip(axs, runs):
        d, bd = s["d"], s["budget"]
        big = max(np.nanmax(np.abs(bd[k])) for k, *_ in BUDGET)
        lin = LINTHRESH.get(tag, 10 ** np.floor(np.log10(big) - 6))   # default: six decades of log range
        ax.axhline(0, color="#b0b0b0", lw=0.5)
        for k, lab, c, ls in BUDGET:
            v = bd[k]
            ok = np.isfinite(v)
            if not ok.any() or np.nanmax(np.abs(v)) == 0.0:
                continue
            ax.plot(d[ok], v[ok], color=c, ls=ls, lw=1.3, label=lab)
        ax.set_yscale("symlog", linthresh=lin, linscale=0.6)
        top = 10 ** np.ceil(np.log10(big))
        ax.set_ylim(-top, top)
        ax.yaxis.set_major_locator(SymmetricalLogLocator(base=10, linthresh=lin,
                                                         subs=[1.0]))
        ax.set_ylabel(r"Electron energy source (W m$^{-3}$)")
        ax.set_xlabel(xlabel(tag))
        ax.set_title(tag, fontsize=9)
        style(ax)
        orient(ax, tag)
        mark_electrode(ax, side(tag), label=tag.lower())
        # every other decade label, or the symlog axis is unreadable
        labs = []
        for t in ax.get_yticks():
            if t == 0:
                labs.append("0")
            else:
                e = int(round(np.log10(abs(t))))
                keep = e % 2 == 0 and abs(t) > lin
                labs.append((f"$-10^{{{e}}}$" if t < 0 else f"$10^{{{e}}}$") if keep else "")
        ax.set_yticks(ax.get_yticks())
        ax.set_yticklabels(labs)
    leg = {}                                      # every term drawn in either panel
    for ax in axs:
        for hh, ll in zip(*ax.get_legend_handles_labels()):
            leg.setdefault(ll, hh)
    h, l = list(leg.values()), list(leg.keys())
    fig.legend(h, l, loc="lower center", ncol=min(len(l), 6), frameon=False,
               bbox_to_anchor=(0.5, -0.1),
               handlelength=2.6, columnspacing=1.4)
    for ax, t, (tag, _) in zip(axs, "ab", runs):
        ax.text(0.97 if away(tag) == "right" else 0.03, 0.95, f"({t})", transform=ax.transAxes,
                ha=away(tag), va="top", fontsize=9)
    for ext in ("pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)


QW_TERMS = [
    ("electron_kinetic",       r"absorbed e$^-$, $2k_BT_e$"),
    ("electron_work_function", r"absorbed e$^-$, $W$"),
    ("emission_kinetic",       r"emitted e$^-$, $-2k_BT_{\rm wall}$"),
    ("emission_work_function", r"emitted e$^-$, $-W$"),
    ("ion_recombination",      r"ions, $\varepsilon_{\rm ion}-W$"),
    ("gas_conduction",         r"gas conduction"),
    ("total",                  r"total $q_w$"),
]


def fig_heat(runs, stem):
    """The eq (24) terms of the electrode heat flux, one horizontal bar each (+ = into the wall)."""
    rc()
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.5), sharey=True, gridspec_kw=dict(wspace=0.08))
    y = np.arange(len(QW_TERMS))[::-1]
    for ax, (tag, s) in zip(axs, runs):
        v = np.array([float(s["info"].get("q_wall." + k, 0.0)) / 1e6 for k, _ in QW_TERMS])
        col = [C_X if k == "total" else (C_I if x >= 0 else C_E)
               for (k, _), x in zip(QW_TERMS, v)]
        ax.barh(y, v, height=0.62, color=col, edgecolor="none")
        ax.axvline(0, color="#8a8a8a", lw=0.6)
        span = np.max(np.abs(v))
        ax.set_xlim(-1.45 * span if v.min() < 0 else -0.08 * span, 1.32 * span)
        for yy, x in zip(y, v):
            if x == 0.0:
                continue
            ax.text(x + (0.02 * span if x >= 0 else -0.02 * span), yy, f"{x:.1f}".replace("-", "\u2212"),
                    va="center", ha="left" if x >= 0 else "right", fontsize=7, color="#3a3a3a")
        ax.set_yticks(y)
        ax.set_yticklabels([l for _, l in QW_TERMS])
        ax.set_xlabel(r"Heat flux into the electrode (MW m$^{-2}$)")
        ax.set_title(f"{tag}, $T_{{\\rm wall}}={float(s['info']['prob.wall_temp']):.0f}$ K",
                     fontsize=9)
        style(ax, xlog=False)
        ax.tick_params(axis="y", which="both", left=False, right=False)
    for ext in ("pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)


def report(tag, s):
    info = s["info"]
    g = lambda k: float(info.get(k, np.nan))
    d, phi = s["d"], s["phi"]
    i = np.argmax(np.abs(phi))
    bd = s["budget"]
    resid = np.nanmax(np.abs(bd["DEDT"][1:-1]))
    big = max(np.nanmax(np.abs(bd[k])) for k in ("JOULE", "ELASTIC", "ION"))
    print(f"\n=== {tag}: {s['name']}  t = {s['t']:.4e} s  j = {np.median(s['j']):.4e} A/m^2")
    print(f"  steady check: max|dE/dt residual| / max|source| = {resid / big:.2e};"
          f"  rebuild vs EENRXN {s['rxn_err']:.1e};  j spread {np.ptp(s['j']) / abs(np.median(s['j'])):.2e}")
    print(f"  T_wall {g('prob.wall_temp'):.0f} K   (balance {int(g('prob.wall_temp_balance'))})")
    print(f"  wall cell: n_e {s['ne'][0]:.3e}  n_i {s['ni'][0]:.3e}  T_e {s['Te'][0]:.0f} K"
          f"  T_h {s['Th'][0]:.0f} K  T_eff {s['Teff'][0]:.0f} K  E {s['E'][0]:.3e} V/m")
    print(f"  phi: at L {phi[-1]:.2f} V;  extreme {phi[i]:.2f} V at x = {d[i]:.2e} m")
    x_s, phi_s, i_s = sheath_edge(d, s["ne"], s["ni"], phi)
    print(f"  sheath edge (|n_i - n_e|/n_i = 1 %): x_s {x_s:.3e} m  phi_s {phi_s:.3f} V")
    # Maximum ion impact energy: an ion that enters the sheath with the local mean energy there
    # (Wannier: random 3/2 k T_eff plus directed 1/2 M v_d^2, i.e. 3/2 k T_h + M v_d^2) and crosses
    # the sheath without a collision, gaining e phi_s.
    v_d2 = 3.0 * K_B * (s["Teff"][i_s] - s["Th"][i_s]) / M_AR
    eps_s = (1.5 * K_B * s["Teff"][i_s] + 0.5 * M_AR * v_d2) / E_CH
    print(f"  ion mean energy at the sheath edge {eps_s:.2f} eV;  max impact energy"
          f" {eps_s + phi_s:.2f} eV;  x ion flux -> {(eps_s + phi_s) * E_CH * g('q_wall.ion_flux') / 1e6:.1f} MW/m^2")
    print(f"  j_i/j at wall {s['ji_j'][0]:.4f};  min {s['ji_j'].min():.4f}  max {s['ji_j'].max():.4f}")
    for k in ("electron_kinetic", "electron_work_function", "emission_kinetic",
              "emission_work_function", "ion_recombination", "gas_conduction", "total"):
        print(f"  q_wall.{k:24s} {g('q_wall.' + k) / 1e6:10.3f} MW/m^2")
    for k in ("e_flux_absorbed", "e_flux_emitted", "ion_flux", "ion_flux_maxwell",
              "ion_flux_drift", "E_face"):
        print(f"  q_wall.{k:24s} {g('q_wall.' + k):10.4e}")
    # where each term peaks
    for k, *_ in BUDGET:
        v = bd[k]
        if np.isfinite(v).any():
            m = np.nanargmax(np.abs(v))
            print(f"  {k:8s} peak {v[m]:+.3e} W/m^3 at x = {d[m]:.2e} m")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("anode")
    ap.add_argument("cathode")
    ap.add_argument("-o", "--outdir", default=HERE)
    a = ap.parse_args(argv)
    runs = [("Cathode", load(a.cathode.rstrip("/"))), ("Anode", load(a.anode.rstrip("/")))]
    for tag, s in runs:
        report(tag, s)
    os.makedirs(a.outdir, exist_ok=True)
    fig_profiles(runs, os.path.join(a.outdir, "base_profiles"))
    fig_budget(runs, os.path.join(a.outdir, "base_een_budget"))
    fig_heat(runs, os.path.join(a.outdir, "base_qwall"))
    print(f"\nwrote base_profiles, base_een_budget, base_qwall (.pdf/.png) in {a.outdir}")


if __name__ == "__main__":
    main()
