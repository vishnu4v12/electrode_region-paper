#!/usr/bin/env python3
"""
Compare the three argon ionization/recombination rate models used for electrode regions:

  HL   Hoffert & Lien (1967)      PlasmeX: exe/hoff_f.dat, exe/hoff_r.dat (k_i, k_r vs T_e[eV])
  BN   Benilov & Naidis (1998)    Almeida et al. (2008) App. A, Semenov et al. (2016);
                                  implemented in PlasmeX exe/almeida_arc_bc.py (rates, saha)
  HH   Hinnov & Hirschberg (1962) Mitchner & Kruger Ch. IX eqs (4.13), (4.20), (6.8);
                                  used by Dinulescu & Pfender (1980) eq (8)

All three are written in the common form  n_dot = k_i n_e n_a - k_r n_e^2 n_i.
For HH:  MK (4.13)  n_dot = alpha [n_1 (n_e^2/n_1)* - n_e^2],  alpha = 5.6e-27 n_e (kT_e)^-4.5 cm^3/s,
so k_r = 5.6e-39 (kT_e[eV])^-4.5 m^6/s and k_i = k_r (n_e^2/n_a)_Saha.

Nothing from PlasmeX is copied: tables and BN formulas are loaded from the live tree through
figures/validation_fig/plasmex.py ($PLASMEX_ROOT, else ~/work/PlasmeX).

Writes rates_comparison.{pdf,png} and prints the numbers quoted in reaction_rates.md.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "figures", "validation_fig"))
from plasmex import exe_module, EXE  # noqa: E402

A = exe_module("almeida_arc_bc")     # BN rates, Saha, binary diffusion coefficients
EV_K, KB = A.EV_K, A.KB
P0 = 101325.0                        # [Pa]

# ---------------------------------------------------------------- the three models
_f = np.loadtxt(os.path.join(EXE, "hoff_f.dat"))
_r = np.loadtxt(os.path.join(EXE, "hoff_r.dat"))


def hl(Te):
    """k_i [m^3/s], k_r [m^6/s] exactly as exe/table.H reads them (linear in T_e[eV] on log k)."""
    x = Te / EV_K
    ki = np.exp(np.interp(x, _f[:, 0], np.log(_f[:, 1])))
    kr = np.exp(np.interp(x, _r[:, 0], np.log(_r[:, 1])))
    return ki, kr


def bn(Te, ne=None):
    """Almeida App. A: k_i = k_dir + beta k_step, k_r = (k_dir + k_step)/S. ne=None -> beta = 1."""
    k_dir, _, k_step = A.rates(Te)
    beta = 1.0 if ne is None else ne / (ne + A.NE0)
    return k_dir + beta * k_step, (k_dir + k_step) / A.saha(Te)


def hh(Te):
    """MK (4.20a): alpha/n_e = 5.6e-27 cm^6/s (kT_e/eV)^-4.5; k_i by detailed balance (4.13)."""
    kr = 5.6e-39 * (Te / EV_K) ** -4.5
    return kr * A.saha(Te), kr


# ---------------------------------------------------------------- LTE relaxation length
def lte_lengths(T, model):
    """Ionization-recombination relaxation at the LTE (Saha, 1 atm, T_e = T_h = T) composition.
    tau = 1/(2 k_i n_a): the linearisation of k_i n_e n_a - k_r n_e^3 about equilibrium at fixed n_a
    (for BN the small d(beta)/dn_e term is dropped).  lambda_r = sqrt(D_amb tau), with the
    ambipolar D_amb = (1 + T_e/T_h) D_ia = 2 D_ia and D_ia from Almeida (45) with BN's
    charge-exchange Qbar_ia -- the same functions almeida_arc_bc.py uses."""
    ne = np.array([A.saha_ne(t, t, P0) for t in T])
    na = P0 / (KB * T) - 2 * ne
    ki = model(T, ne)[0] if model is bn else model(T)[0]
    tau = 1.0 / (2 * ki * na)
    n = na + 2 * ne
    D_ia = A.nD(T, A.M_AR / 2, A.qbar_ia(T)) / n
    return tau, np.sqrt(2 * D_ia * tau), ne


# ---------------------------------------------------------------- numbers for the .md
def report():
    print("hoff tables vs Saha (g+/g1 = 6, I = 15.75 eV): k_f/(k_r S)")
    for x in (1.0, 2.0, 4.0, 8.0):
        Te = x * EV_K
        ki, kr = hl(Te)
        print(f"  T_e = {x:4.1f} eV   {ki / (kr * A.saha(Te)):.4f}")

    # least-squares fit of HL k_i to C vbar (eps* + 2 kT_e) exp(-eps*/kT_e)
    x = _f[_f[:, 0] > 0.6, 0]
    vbar = np.sqrt(8 * A.E_CH * x / (np.pi * A.ME))
    best = min(((np.std(np.log(_f[_f[:, 0] > 0.6, 1] / (vbar * (e + 2 * x) * np.exp(-e / x)))), e)
                for e in np.arange(11.40, 11.90, 0.005)))
    e = best[1]
    C = np.exp(np.mean(np.log(_f[_f[:, 0] > 0.6, 1] / (vbar * (e + 2 * x) * np.exp(-e / x)))))
    print(f"HL k_i fit: eps* = {e:.3f} eV, C = {C:.3e} m^2/eV, rms log residual {best[0]:.3e}")

    print("\nrate coefficients (beta = 1 for BN)")
    print("   T[K]   T_e[eV]    k_i HL     k_i BN     k_i HH    | k_r HL     k_r BN     k_r HH   | BN/HL  HH/HL")
    for T in (6000, 8000, 10000, 12000, 15000, 20000, 25000):
        a, b, c = hl(T), bn(T), hh(T)
        print(f"  {T:6d}  {T / EV_K:6.3f}   {a[0]:.2e}  {b[0]:.2e}  {c[0]:.2e} | "
              f"{a[1]:.2e}  {b[1]:.2e}  {c[1]:.2e} | {b[0] / a[0]:5.2f}  {c[0] / a[0]:7.1f}")

    print("\nBN radiative-escape factor: k_i(n_e)/k_i(beta=1)")
    for T in (8000, 10000, 12000):
        print(f"  T = {T}: " + "  ".join(f"n_e={n:.0e}: {bn(T, n)[0] / bn(T)[0]:.2f}"
                                         for n in (1e19, 1e20, 1e21, 1e22)))

    print("\nDinulescu & Pfender (1980) Table I: their alpha vs alpha = k_r n_e from each model")
    print("  T_e[K]   n_e[m^-3]   alpha_DP    alpha_HH    alpha_BN    alpha_HL")
    for Te, ne, a_dp in ((12000, 6.2e22, 2.6e-17), (12500, 8.1e22, 2.6e-17),
                         (13000, 1.1e23, 2.6e-17), (14000, 1.4e23, 2.7e-17)):
        print(f"  {Te:6d}   {ne:.1e}    {a_dp:.1e}     {hh(Te)[1] * ne:.2e}   "
              f"{bn(Te)[1] * ne:.2e}   {hl(Te)[1] * ne:.2e}")

    print("\nLTE relaxation at 1 atm: tau [s], lambda_r [mm]")
    T = np.array([6000., 8000., 10000., 12000., 15000.])
    rows = {m.__name__: lte_lengths(T, m) for m in (hl, bn, hh)}
    print("   T[K]   n_e,Saha    " + "   ".join(f"tau {k:2s}     lam {k:2s} " for k in rows))
    for i, t in enumerate(T):
        print(f"  {t:6.0f}  {rows['hl'][2][i]:.2e}  " +
              "  ".join(f"{rows[k][0][i]:.2e}  {rows[k][1][i] * 1e3:.2e}" for k in rows))


# ---------------------------------------------------------------- figure
# categorical slots 1-3 of the reference palette, fixed per model (identity follows the model)
COL = {"hl": "#2a78d6", "bn": "#eb6834", "hh": "#1baf7a"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
LAB = {"hl": "Hoffert & Lien (PlasmeX)", "bn": "Benilov & Naidis (Almeida)",
       "hh": "Hinnov & Hirschberg (MK / Pfender)"}


def style(ax):
    ax.grid(True, which="major", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def ev_axis(ax):
    top = ax.secondary_xaxis("top", functions=(lambda t: t / EV_K, lambda e: e * EV_K))
    top.set_xlabel(r"$kT_e$ (eV)", fontsize=8, color=INK2)
    top.tick_params(colors=INK2, labelsize=7)
    for s in ("top",):
        top.spines[s].set_color(INK2)


def plot():
    T = np.linspace(5000, 30000, 600)
    ne_low = 1e20
    fig, axs = plt.subplots(2, 2, figsize=(10, 7.8), constrained_layout=True)
    (a, b), (c, d) = axs

    for ax, j, ylab, title in ((a, 0, r"$k_i$ (m$^3$ s$^{-1}$)", "(a) Ionization"),
                               (b, 1, r"$k_r$ (m$^6$ s$^{-1}$)", "(b) Three-body recombination")):
        ax.semilogy(T, hl(T)[j], color=COL["hl"], lw=2, label=LAB["hl"])
        ax.semilogy(T, bn(T)[j], color=COL["bn"], lw=2, label=LAB["bn"] + r", $n_e\gg n_e^{(0)}$")
        if j == 0:
            ax.semilogy(T, bn(T, ne_low)[0], color=COL["bn"], lw=2, ls="--",
                        label=LAB["bn"] + r", $n_e=10^{20}$ m$^{-3}$")
        ax.semilogy(T, hh(T)[j], color=COL["hh"], lw=2, label=LAB["hh"])
        ax.set_xlim(5000, 30000)
        ax.set_xlabel(r"$T_e$ (K)")
        ax.set_ylabel(ylab)
        ax.set_title(title, loc="left", fontsize=10, color=INK)
        style(ax)
        ev_axis(ax)
    a.set_ylim(1e-26, 1e-13)
    b.set_ylim(1e-44, 1e-36)
    a.legend(fontsize=7.5, frameon=False, loc="lower right", labelcolor=INK)

    # (c) ratio to the PlasmeX rates -- the same for k_i and k_r when both obey Saha
    r_bn = bn(T)[0] / hl(T)[0]
    r_hh = hh(T)[0] / hl(T)[0]
    c.semilogy(T, r_bn, color=COL["bn"], lw=2)
    c.semilogy(T, bn(T, ne_low)[0] / hl(T)[0], color=COL["bn"], lw=2, ls="--")
    c.semilogy(T, r_hh, color=COL["hh"], lw=2)
    c.axhline(1.0, color=COL["hl"], lw=2)
    i = np.searchsorted(T, 16000)
    c.text(T[i], 0.85, "Hoffert & Lien", color=INK2, fontsize=8, ha="center", va="top")
    c.text(T[i], r_bn[i] * 1.15, r"BN, $n_e\gg n_e^{(0)}$", color=INK2, fontsize=8, ha="center", va="bottom")
    c.text(T[i], bn(T[i], ne_low)[0] / hl(T[i])[0] * 0.8, r"BN, $n_e=10^{20}$ m$^{-3}$", color=INK2,
           fontsize=8, ha="center", va="top")
    c.text(T[i], r_hh[i] * 1.3, "Hinnov & Hirschberg", color=INK2, fontsize=8, ha="center", va="bottom")
    c.set_xlim(5000, 30000)
    c.set_ylim(0.05, 3e3)
    c.set_xlabel(r"$T_e$ (K)")
    c.set_ylabel(r"$k/k_{\rm HL}$  (same for $k_i$ and $k_r$)")
    c.set_title("(c) Ratio to the PlasmeX (Hoffert & Lien) rates", loc="left", fontsize=10, color=INK)
    style(c)
    ev_axis(c)

    # (d) relaxation length at the LTE composition
    TT = np.linspace(6000, 20000, 300)
    for k, m in (("hl", hl), ("bn", bn), ("hh", hh)):
        d.semilogy(TT, lte_lengths(TT, m)[1] * 1e3, color=COL[k], lw=2, label=LAB[k])
    d.axhspan(0.1, 1.0, color="#ececea", zorder=0)
    d.text(15500, 0.55, "anode boundary layer, 0.1-1 mm (DP 1980)", color=INK2, fontsize=8,
           ha="center", va="center")
    d.set_xlim(6000, 20000)
    d.set_xlabel(r"$T = T_e = T_h$ (K), Saha composition at 1 atm")
    d.set_ylabel(r"$\lambda_r=\sqrt{D_{\rm amb}\,\tau_r}$ (mm)")
    d.set_title(r"(d) Ionization-recombination length, $\tau_r = 1/(2k_in_a)$", loc="left",
                fontsize=10, color=INK)
    d.legend(fontsize=7.5, frameon=False, loc="upper right", labelcolor=INK)
    style(d)

    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, f"rates_comparison.{ext}"), dpi=200)


if __name__ == "__main__":
    report()
    plot()
