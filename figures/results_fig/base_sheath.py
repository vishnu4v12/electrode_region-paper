"""
Results-section figures on the momentum balance and the collisionality of the near-electrode layers
of the two base cases (PlasmeX runs/base_case_cathode, base_case_anode).

    python base_sheath.py <anode plotfile> <cathode plotfile> [-o outdir]

Writes, each as .pdf and .png, and prints the numbers quoted in the text (base_sheath_report.txt):
    base_momentum   every term of the electron and ion momentum balances, main.tex eqs (2), (3)
                    (eq:mom-e, eq:mom-i), against the distance from the electrode: electric force,
                    pressure gradient, the two collisional drags combined, the thermal (Soret) force,
                    and the species inertia that the model neglects, evaluated a posteriori
    base_sheath     (a, b) mean free paths against the distance from the electrode;
                    (c, d) the electron Boltzmann relation; (e, f) the ion Mach number u_i/u_B

TERMS (per unit volume, x-component, positive = towards the electrode), on the AMR leaf cells
(bohm.sample), with the solver's own coefficients (bohm.electron_transport, ion_transport,
electron_thermal_coeffs; src/aider.H CalcTransportCoeffs / CalcExtraVelocities):
  electron   electric  -e n_e E
             pressure  -d(n_e k T_e)/dx
             drag      -(e n_e/|mu_e|) [w_en (u_e - u_n) + w_ei (u_e - u_i)],  w_en = coul_fac, w_ei = 1 - w_en
             thermal   -C_e n_e k dT_e/dx
             inertia   -m_e n_e u_e du_e/dx                     (not in the model)
  ion        electric  +e n_i E
             pressure  -d(n_i k T_h)/dx
             drag      -(e n_i/mu_in) [(u_i - u_n) + r (u_i - u_e)],  r = nu_ie/nu_in
             thermal   +A_i n_e k dT_e/dx                       (C_i of the paper)
             mixture   -(n_i/n_h) [e(n_i - n_e) E - dp/dx]       (Y_i F; zero at steady state, printed)
             inertia   -M n_i u_i du_i/dx                       (not in the model)
u_i = JIX/(e n_i), u_e = -JEX/(e n_e), u_n = UNX. The sum of the retained terms is printed as a check
of the reconstruction (cell-centred here, face-based in the solver). The wall cell is dropped: JIX/JEX
there hold the FACE value while n is the cell value.

MEAN FREE PATHS as in table tab:scales (reaction_rates/scale_table.py), with the local state:
  lambda_e = 1/(n_n Qbar_en(T_e) + n_i Qbar_ei(T_e, n_e)),
  lambda_i = 1/(n_n Qbar_in(T_h) + n_i Qbar_ii),  Qbar_ii = Qbar_ei at T_h,
  lambda_in = 1/(n_n Qbar_in(T_h)), the ion-atom part alone (the ion friction; report only),
  lambda_Th = 1/((n_i + n_n) Qbar_in(T_h)), the heavy-particle energy-exchange length of tab:scales,
  L_Th = T_h/|dT_h/dx|, the length over which T_h changes (single heavy temperature: lambda_Th << L_Th),
  lambda_De = sqrt(eps0 k T_e/(n_e e^2)).
N = int_0^x_s dx/lambda is the number of mean free paths across the sheath.

BOLTZMANN: n_B(x) = n_e(x_s) exp(-int_{x_s}^x e E/(k T_e) dx'), the Boltzmann distribution with the local
T_e, started at the sheath edge (outermost cell with |n_i - n_e|/n_i > 1 %; x_s itself is the interpolated
1 % crossing, base_case.sheath_edge, as table tab:base).
BOHM: u_B = sqrt(k (T_e + T_h)/M); the ion Mach number is u_i/u_B with u_i towards the electrode.
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter, SymmetricalLogLocator

from base_case import HERE, C_E, C_I, C_N, C_X, rc, orient, side, away, _um_label, sheath_edge
sys.path.insert(0, os.path.join(HERE, "..", "validation_fig"))
from plasmex import EXE, exe_module  # noqa: E402
from electrode_mark import mark_electrode  # noqa: E402

b = exe_module("bohm")
A = exe_module("almeida_arc_bc")
E, KB, M, ME, EPS0 = b.ECHARGE, b.K_B, b.ARMASS, b.EMASS, b.EPS0
FIELDS = ("NE", "ARP", "AR", "PHI", "ETEMP", "TGAS", "EFLDX", "JEX", "JIX", "UNX", "HVELX")
XLIM = (3e-9, 1e-4)          # near-electrode window [m]; the first interior cell is at 5.7 nm

# force terms: key -> (label, colour, line style, width); the neglected inertia in black
TERMS = [
    ("E",   "electric force",           "#2a78d6", "-",  1.4),
    ("p",   "pressure gradient",        "#eb6834", "-",  1.4),
    ("dr",  "drag (both partners)",     "#1baf7a", "-",  1.4),
    ("T",   "thermal force",            "#eda100", "--", 1.3),
    ("in",  "inertia (neglected)",      "#2b2b2b", ":",  1.6),
]


def grad(x, f):
    return np.array([b.cell_gradient({"x": x}, f, i)[0] for i in range(len(x))])


def load(pf):
    d = b.sample(pf, FIELDS, want_derived=False)
    if d["missing"]:
        raise SystemExit(f"{pf}: missing {d['missing']}")
    x = d["x"]
    ne, ni, nn = d["NE"], d["ARP"], d["AR"]
    Te, Th, Ex, phi = d["ETEMP"], d["TGAS"], d["EFLDX"], d["PHI"]
    ui = d["JIX"] / (E * ni)
    ue = -d["JEX"] / (E * ne)
    un = d["UNX"]

    et = b.electron_transport(Te, nn, ne, ni, exe_dir=EXE)
    mu_e, cf = et["mu"], et["fac"]
    _, mu_in, _ = b.ion_transport(Th, nn, np.abs(Ex))
    r = (1.0 - cf) * mu_in * ne / (np.abs(mu_e) * ni)
    c_e, a_i, _ = b.electron_thermal_coeffs(ni, ne, et["te_ev"], nn, et["D"], True, 1, EXE)
    dTe = grad(x, Te)
    nh = ni + nn
    p = ne * KB * Te + nh * KB * Th

    t = {"x": x, "d": d["hi"] - x, "ne": ne, "ni": ni, "nn": nn, "Te": Te, "Th": Th, "E": Ex,
         "phi": phi, "ui": ui, "ue": ue, "un": un, "u": d["HVELX"]}
    t["e"] = {"E": -E * ne * Ex,
              "p": -grad(x, ne * KB * Te),
              "dr": -(E * ne / np.abs(mu_e)) * (cf * (ue - un) + (1.0 - cf) * (ue - ui)),
              "T": -c_e * ne * KB * dTe,
              "in": -ME * ne * ue * grad(x, ue)}
    t["i"] = {"E": E * ni * Ex,
              "p": -grad(x, ni * KB * Th),
              "dr": -(E * ni / mu_in) * ((ui - un) + r * (ui - ue)),
              "T": a_i * ne * KB * dTe,
              "in": -M * ni * ui * grad(x, ui),
              "mix": -(ni / nh) * (E * (ni - ne) * Ex - grad(x, p))}

    # mean free paths (table tab:scales definitions, local state)
    qen = np.array([A.qbar_ea(T)[0] for T in Te])
    qei = A.qbar_ei(Te, ne)
    qin = A.qbar_ia(Th)
    qii = A.qbar_ei(Th, ne)
    t["lam_e"] = 1.0 / (nn * qen + ni * qei)
    t["lam_i"] = 1.0 / (nn * qin + ni * qii)
    t["lam_in"] = 1.0 / (nn * qin)
    t["lam_Th"] = 1.0 / ((ni + nn) * qin)
    t["L_Th"] = Th / np.maximum(np.abs(grad(x, Th)), 1e-300)
    t["Teff"] = b.ion_transport(Th, nn, np.abs(Ex))[0]
    t["lam_D"] = np.sqrt(EPS0 * KB * Te / (ne * E * E))

    # Bohm velocity and the ion Mach number
    t["uB"] = np.sqrt(KB * (Te + Th) / M)
    t["uB0"] = np.sqrt(KB * Te / M)

    # drop the wall cell (face-valued J against cell-valued n), order by distance
    keep = np.argsort(t["d"])[1:]
    for k in list(t):
        if isinstance(t[k], np.ndarray):
            t[k] = t[k][keep]
        elif isinstance(t[k], dict):
            t[k] = {kk: vv[keep] for kk, vv in t[k].items()}
    # sheath edge (table tab:base, wall cell included): the interpolated 1 % crossing x_s; the
    # "at x_s" values are those of its cell, the outermost one with |n_i - n_e|/n_i > 1 %
    dd = d["hi"] - x
    o = np.argsort(dd)
    x_s, _, i_c = sheath_edge(dd[o], ne[o], ni[o], phi[o])
    t["x_s"] = float(x_s)
    # Boltzmann density from the sheath edge, with the local T_e
    dist = t["d"]
    i_s = int(np.argmin(np.abs(dist - dd[o][i_c])))
    f = E * t["E"] / (KB * t["Te"])          # d ln n_B / dx (lab x) = -e E/(k T_e); x = hi - dist
    # integrate in dist: d ln n_B / d(dist) = +e E/(k T_e)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(dist))])
    t["nB"] = t["ne"][i_s] * np.exp(cum - cum[i_s])
    t["i_s"] = i_s
    t["name"] = os.path.basename(os.path.abspath(pf))
    return t


def report(tag, t):
    dist, x_s, i_s = t["d"], t["x_s"], t["i_s"]
    print(f"\n=== {tag}: {t['name']}   sheath edge x_s = {x_s:.3e} m")
    for sp, nm in (("i", "ion"), ("e", "electron")):
        F = t[sp]
        tot = F["E"] + F["p"] + F["dr"] + F["T"] + (F.get("mix", 0.0))
        big = np.max(np.abs(np.vstack([F[k] for k in ("E", "p", "dr", "T")])), axis=0)
        print(f"  {nm}: |sum of retained terms| / max|term|: median {np.median(np.abs(tot) / big):.1e}, "
              f"max within {XLIM[1]:.0e} m {np.max((np.abs(tot) / big)[dist <= XLIM[1]]):.1e}")
        if sp == "i":
            print(f"  ion mixture term max |Y_i F|/max|term| = {np.max(np.abs(F['mix']) / big):.1e}")
        for xq in (2e-9, 1e-8, 1e-7, 3e-7, x_s, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4):
            k = int(np.argmin(np.abs(dist - xq)))
            print(f"    x={dist[k]:.2e}: " + "  ".join(f"{key} {F[key][k]:+.2e}" for key in ("E", "p", "dr", "T", "in"))
                  + f"   |in|/|dr| {abs(F['in'][k] / F['dr'][k]):.2e}  |T|/|dr| {abs(F['T'][k] / F['dr'][k]):.2e}")
        rat = np.abs(F["in"]) / np.abs(F["dr"])
        for lim in (1.0, 0.1, 0.01):
            w = dist[rat > lim]
            print(f"    |inertia|/|drag| > {lim}: within {w.max() if w.size else 0:.2e} m;"
                  f" max {rat.max():.2f} at {dist[np.argmax(rat)]:.2e} m")
    # mean free paths
    for k in ("lam_e", "lam_i", "lam_in", "lam_Th", "L_Th", "lam_D"):
        v = t[k]
        print(f"  {k:7s} wall-adjacent {v[0]:.2e}  at x_s {v[i_s]:.2e}  at 10 um "
              f"{v[np.argmin(np.abs(dist - 1e-5))]:.2e}  at 100 um {v[np.argmin(np.abs(dist - 1e-4))]:.2e}")
    for k in ("lam_e", "lam_i", "lam_in", "lam_Th"):
        m = dist <= x_s
        n = np.trapezoid(1.0 / t[k][m], dist[m]) + dist[m][0] / t[k][m][0]
        print(f"  mean free paths across the sheath, {k}: N = {n:.2f}")
    # where lambda = distance
    kn = t["lam_Th"] / t["L_Th"]
    for xq in (2e-9, 1e-8, 1e-7, x_s, 1e-6, 1e-5, 1e-4, 3e-4):
        k = int(np.argmin(np.abs(dist - xq)))
        print(f"    x={dist[k]:.2e}: lam_Th {t['lam_Th'][k]:.2e}  L_Th {t['L_Th'][k]:.2e}  Kn_Th {kn[k]:.2e}"
              f"  T_h {t['Th'][k]:.0f}  T_eff {t['Teff'][k]:.0f}  (T_eff-T_h)/T_h {(t['Teff'][k]-t['Th'][k])/t['Th'][k]:.3f}")
    for lim in (1.0, 0.1, 0.01):
        w = dist[kn > lim]
        print(f"    lam_Th/L_Th > {lim}: within {w.max() if w.size else 0:.2e} m; max {kn.max():.2e} at {dist[np.argmax(kn)]:.2e} m")
    rt = (t["Teff"] - t["Th"]) / t["Th"]
    for lim in (1.0, 0.1):
        w = dist[rt > lim]
        print(f"    (T_eff - T_h)/T_h > {lim}: within {w.max() if w.size else 0:.2e} m")
    for k in ("lam_e", "lam_i", "lam_in", "lam_Th"):
        w = dist[t[k] > dist]
        print(f"  {k} > x within {w.max() if w.size else 0:.2e} m")
    # Boltzmann
    m = dist <= x_s
    dev = t["ne"] / t["nB"] - 1.0
    print(f"  Boltzmann n_e/n_B - 1 in the sheath: wall-adjacent {dev[0]:+.3f}, max |.| {np.max(np.abs(dev[m])):.3f}")
    for xq in (1e-6, 3e-6, 1e-5, 3e-5, 1e-4):
        k = int(np.argmin(np.abs(dist - xq)))
        print(f"    x={dist[k]:.2e}: n_e/n_B - 1 = {dev[k]:+.3f}")
    eb = np.abs(t["e"]["E"] + t["e"]["p"]) / np.abs(t["e"]["E"])
    print(f"  |eE n_e + grad p_e|/|e n_e E|: wall-adjacent {eb[0]:.3f}, at x_s {eb[i_s]:.3f}, "
          f"max in sheath {np.max(eb[m]):.3f}, at 10 um {eb[np.argmin(np.abs(dist - 1e-5))]:.3f}")
    # Bohm
    Mi = t["ui"] / t["uB"]
    print(f"  u_i: wall-adjacent {t['ui'][0]:.0f} m/s, at x_s {t['ui'][i_s]:.0f} m/s;  u_B at x_s {t['uB'][i_s]:.0f}"
          f" (sqrt(kTe/M) {t['uB0'][i_s]:.0f}) m/s")
    print(f"  Mach u_i/u_B: wall-adjacent {Mi[0]:.3f}, at x_s {Mi[i_s]:.3f}, max {Mi.max():.3f}")
    w = dist[Mi >= 1.0]
    print(f"  u_i >= u_B within {w.max() if w.size else 0:.2e} m")
    print(f"  heavy-gas velocity max |u| {np.max(np.abs(t['u'])):.2e} m/s; max |u_n| {np.max(np.abs(t['un'])):.2e}")


def near_axis(ax):
    ax.set_xscale("log")
    ax.set_xlim(*XLIM)
    ax.xaxis.set_major_locator(FixedLocator(10.0 ** np.arange(-9, -3.5)))
    ax.xaxis.set_minor_locator(FixedLocator([m * 10.0 ** e for e in range(-9, -4) for m in range(2, 10)]))
    ax.xaxis.set_major_formatter(FuncFormatter(_um_label))
    ax.xaxis.set_minor_formatter(NullFormatter())


def style(ax):
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.tick_params(which="major", length=3.5, width=0.6)
    ax.tick_params(which="minor", length=2.0, width=0.5)
    for sp in ax.spines.values():
        sp.set_linewidth(0.6)


def xlabel(tag):
    return rf"Distance from the {tag.lower()}, $x$ ($\mu$m)"


def mark_sheath_edge(ax, t):
    ax.axvline(t["x_s"], color="#8a8a8a", lw=0.7, ls=(0, (3, 2)), zorder=0.6)


def symlog_labels(ax, lin):
    labs = []
    for v in ax.get_yticks():
        if v == 0:
            labs.append("0")
        else:
            e = int(round(np.log10(abs(v))))
            keep = e % 2 == 0 and abs(v) > lin
            labs.append((f"$-10^{{{e}}}$" if v < 0 else f"$10^{{{e}}}$") if keep else "")
    ax.set_yticks(ax.get_yticks())
    ax.set_yticklabels(labs)


def fig_momentum(runs, stem):
    rc()
    fig, axs = plt.subplots(2, 2, figsize=(7.0, 5.8), sharex="col",
                            gridspec_kw=dict(hspace=0.16, wspace=0.28))
    for c, (tag, t) in enumerate(runs):
        for r, (sp, nm) in enumerate((("i", "Ion"), ("e", "Electron"))):
            ax = axs[r, c]
            F = t[sp]
            big = max(np.max(np.abs(F[k])) for k, *_ in TERMS)
            lin = 10 ** np.floor(np.log10(big) - 6)
            ax.axhline(0, color="#b0b0b0", lw=0.5)
            for k, lab, col, ls, lw in TERMS:
                ax.plot(t["d"], F[k], color=col, ls=ls, lw=lw, label=lab)
            ax.set_yscale("symlog", linthresh=lin, linscale=0.5)
            top = 10 ** np.ceil(np.log10(2.0 * big))        # a curve never rides the frame
            ax.set_ylim(-top, top)
            ax.yaxis.set_major_locator(SymmetricalLogLocator(base=10, linthresh=lin, subs=[1.0]))
            near_axis(ax)
            style(ax)
            symlog_labels(ax, lin)
            mark_sheath_edge(ax, t)
            ax.set_ylabel(f"{nm} momentum source (N m$^{{-3}}$)")
            mark_electrode(ax, side(tag), label=tag.lower() if r == 0 else None)
            ax.text(0.0 if away(tag) == "left" else 1.0, 1.015, f"({'abcd'[2 * r + c]})",
                    transform=ax.transAxes, fontsize=9, va="bottom", ha=away(tag))
        axs[0, c].set_title(tag, fontsize=9, pad=10)
        axs[1, c].set_xlabel(xlabel(tag))
        orient(axs[0, c], tag)
    for c in range(2):
        fig.align_ylabels(axs[:, c])
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 0.035),
               handlelength=2.6, columnspacing=1.3)
    for ext in ("pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)


def fig_sheath(runs, stem):
    rc()
    fig, axs = plt.subplots(3, 2, figsize=(7.0, 7.4), sharex="col",
                            gridspec_kw=dict(height_ratios=[1.2, 1, 0.9], hspace=0.16, wspace=0.28))
    for c, (tag, t) in enumerate(runs):
        d = t["d"]
        ax = axs[0, c]
        ax.plot(d, t["lam_e"], color=C_E, lw=1.4, label=r"$\lambda_e$")
        ax.plot(d, t["lam_i"], color=C_I, lw=1.4, label=r"$\lambda_i$")
        ax.plot(d, t["lam_Th"], color=C_I, lw=1.0, ls="--", label=r"$\lambda_{T_h}$")
        ax.plot(d, t["lam_D"], color=C_N, lw=1.4, ls="-.", label=r"$\lambda_{De}$")
        ax.plot(XLIM, XLIM, color="#8a8a8a", lw=0.7, ls=":", label=r"$x$")
        ax.set_yscale("log")
        ax.set_ylim(3e-10, 3e-5)
        ax.set_ylabel("Length (m)")
        ax.set_title(tag, fontsize=9, pad=10)

        ax = axs[1, c]
        ax.plot(d, t["ne"], color=C_E, lw=1.4, label=r"$n_e$")
        ax.plot(d, t["nB"], color=C_X, lw=1.1, ls="--", label=r"Boltzmann relation, $n_B$")
        ax.plot(d, t["ni"], color=C_I, lw=1.0, ls=(0, (4, 2)), label=r"$n_i$")
        ax.set_yscale("log")
        lo = 0.5 * min(t["ne"].min(), t["nB"][d < XLIM[1]].min())
        hi = 2.0 * max(t["ne"][d < XLIM[1]].max(), t["ni"][d < XLIM[1]].max())
        ax.set_ylim(lo, hi)
        ax.set_ylabel(r"Number density (m$^{-3}$)")

        ax = axs[2, c]
        ax.axhline(1.0, color="#8a8a8a", lw=0.6)
        ax.plot(d, t["ui"] / t["uB"], color=C_I, lw=1.4, label=r"$u_i/u_B$")
        ax.set_yscale("log")
        ax.set_ylim(3e-3, 6.0)
        ax.set_ylabel(r"Ion Mach number $u_i/u_B$")
        ax.set_xlabel(xlabel(tag))
        ax.text(0.5, 0.95, r"$u_B=\sqrt{k_B(T_e+T_h)/M}$",
                transform=ax.transAxes, ha="center", va="top", fontsize=7.5, color="#3a3a3a")
        for r in range(3):
            near_axis(axs[r, c])
            style(axs[r, c])
            mark_sheath_edge(axs[r, c], t)
            mark_electrode(axs[r, c], side(tag), label=tag.lower() if r == 0 else None)
            axs[r, c].text(0.0 if away(tag) == "left" else 1.0, 1.015, f"({'abcdef'[2 * r + c]})",
                           transform=axs[r, c].transAxes, fontsize=9, va="bottom", ha=away(tag))
        orient(axs[0, c], tag)
    for r in range(2):
        axs[r, 0].legend(frameon=False, loc="lower right", handlelength=2.4, fontsize=7)
    for c in range(2):
        fig.align_ylabels(axs[:, c])
    for ext in ("pdf", "png"):
        fig.savefig(f"{stem}.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("anode")
    ap.add_argument("cathode")
    ap.add_argument("-o", "--outdir", default=HERE)
    a = ap.parse_args(argv)
    runs = [("Cathode", load(a.cathode.rstrip("/"))), ("Anode", load(a.anode.rstrip("/")))]
    for tag, t in runs:
        report(tag, t)
    os.makedirs(a.outdir, exist_ok=True)
    fig_momentum(runs, os.path.join(a.outdir, "base_momentum"))
    fig_sheath(runs, os.path.join(a.outdir, "base_sheath"))
    print(f"\nwrote base_momentum, base_sheath (.pdf/.png) in {a.outdir}")


if __name__ == "__main__":
    main()
