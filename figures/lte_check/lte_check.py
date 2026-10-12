"""
Accuracy of the two-limit bridging of Appendix A, checked at LTE.

At T_e = T_h = T and the Saha composition of argon at p = 1 atm, the transport coefficients of
the model are evaluated exactly as in Appendix A and compared with the LTE argon properties in
Ar_LTE_prop.dat (Boulos, Fauchais & Pfender tables). The two figures are the journal versions
used in Appendix A:

  sigma_LTE.pdf   electrical conductivity
                    sigma = n_e e^2 / (m_e nubar^M_e),  nubar^M_e = nubar^M_en + nubar^E_ei/1.975
                  plus the two single-limit curves (electron-neutral only, Spitzer-Harm only).

  kappa_LTE.pdf   thermal conductivity
                    kappa_e (lambda_en, lambda_ei joined in series)  +  kappa_h (eq. kappa-h)
                    + kappa_R, the reactive conductivity of Mutation++ (lambda.dat, column 3).
                  The model has no kappa_R term of its own; the LTE value is a total conductivity.

Printed only (not plotted): the model's own implied kappa_R, from the ambipolar flux of ion
pairs (derivation in reactive_conductivity.md):
    kappa_R = (eps_ion + 5/2 k T) (n_n/n_h) 2 mu_in (k/e) d(n_e T)/dT

The composition is PlasmeX's LTE fit (argon_lte, via bohm.n_neutral_lte); it reproduces the n_e
column of PlasmeX/exe/LTE_COND.dat to 0.1% (printed as n_e/n_ref).

Nothing from PlasmeX is copied: bohm.py and ar_e_lorentz.dat are read from the live tree
(validation_fig/plasmex.py).   Run:  python3 lte_check.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "validation_fig"))
from plasmex import exe_module, EXE  # noqa: E402

b = exe_module("bohm")
KB, E, ME = b.K_B, b.ECHARGE, b.EMASS

P = 101325.0
LTE_FILE = os.path.join(HERE, "Ar_LTE_prop.dat")
MUTPP_LAMBDA = os.path.join(HERE, "lambda.dat")     # Mutation++: T, total, reactive

# reference palette (dataviz skill, light mode); the LTE data is the neutral ink
INK, BLUE, ORANGE, AQUA, MAGENTA = "#1f1f1e", "#2a78d6", "#eb6834", "#1baf7a", "#e87ba4"
# journal style: serif text and Computer Modern maths to match the LaTeX body, boxed axes with
# inward ticks, no grid; sized to be included at about half the text width without scaling
FIGSIZE = (3.4, 2.6)
plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"], "mathtext.fontset": "cm",
                     "font.size": 8, "axes.linewidth": 0.6, "lines.linewidth": 1.1,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
                     "legend.frameon": False, "legend.fontsize": 7.5, "legend.handlelength": 2.4})


def load_lte():
    """T, kappa, sigma from the tab-separated table; rows with empty columns are dropped."""
    rows = []
    with open(LTE_FILE, encoding="utf-8", errors="replace") as f:
        next(f)
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) > 6 and c[5].strip() and c[6].strip():
                rows.append((float(c[0]), float(c[5]), float(c[6])))
    a = np.array(rows)
    return a[:, 0], a[:, 1], a[:, 2]


def composition(T):
    """Saha at p, T_e = T_h = T: (n_e = n_i, n_n), from the code's own LTE fit."""
    nn = b.n_neutral_lte(T, P)
    ne = 0.5 * (P / (KB * T) - nn)
    return ne, nn


def lorentz(col, Te_ev):
    x, y = b.load_table("ar_e_lorentz.dat", EXE, col=col)
    return np.interp(Te_ev, x, y)


def model(T):
    ne, nn = composition(T)
    ni = ne
    te = T * KB / E
    lnl = b.coulomb_log(ne, te)

    # collision frequencies, eqs. (nu-M), (spitzer), (nu-ei)
    nuM_en = nn * E**2 / (ME * lorentz(2, te))            # = n_e e^2/(m_e sigma_en)
    nuE_ei = b.NUEI_C * ni * lnl / te**1.5
    nuM_ei = nuE_ei / b.SIGMA_SH
    sig = ne * E**2 / (ME * (nuM_en + nuM_ei))
    sig_en = ne * E**2 / (ME * nuM_en)
    sig_ei = ne * E**2 / (ME * nuM_ei)

    # electron thermal conductivity, eq. (bridge)
    lam_en = (ne / nn) * lorentz(5, te)
    lam_ei = b.LAMBDA_SH * ne * KB**2 * T / (ME * nuE_ei)
    k_e = 1.0 / (1.0 / lam_en + 1.0 / lam_ei)
    k_h = b.k_heavy(nn, ni, T, lnl)

    # reactive part carried by the ambipolar flux (see module docstring)
    h = 1.0
    ne_p, _ = composition(T + h)
    ne_m, _ = composition(T - h)
    dneT = ((T + h) * ne_p - (T - h) * ne_m) / (2 * h)
    mu_in = b.ion_mobility0(T, nn)
    k_r = (b.E_ION * E + 2.5 * KB * T) * (nn / (nn + ni)) * 2 * mu_in * (KB / E) * dneT

    return dict(ne=ne, nn=nn, sig=sig, sig_en=sig_en, sig_ei=sig_ei,
                k_e=k_e, k_h=k_h, k_r=k_r)


def journal_axes(ax):
    ax.tick_params(direction="in", which="both", top=True, right=True)
    ax.set_xlim(5000, 24000)
    ax.set_xticks([5000, 10000, 15000, 20000])
    ax.set_xlabel(r"$T$ (K)")


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, f"{name}.{ext}"), bbox_inches="tight", dpi=300)
    plt.close(fig)


def fig_sigma(T, sig_lte, m, Tmin):
    k = T >= Tmin
    d = k & (T % 500 == 0)                       # data markers every 500 K
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.semilogy(T[d], sig_lte[d], "o", ms=3.5, mfc="none", mew=0.8, color=INK, label="LTE data")
    ax.semilogy(T[k], m["sig"][k], color=BLUE, label=r"$\bar\nu^M_{en}+\bar\nu^M_{ei}$")
    ax.semilogy(T[k], m["sig_en"][k], "--", color=ORANGE, label=r"$\bar\nu^M_{en}$ only")
    ax.semilogy(T[k], m["sig_ei"][k], ":", color=AQUA, lw=1.4, label=r"$\bar\nu^M_{ei}$ only")
    ax.set_ylim(1e0, 1e5)
    ax.set_ylabel(r"$\sigma$ (S m$^{-1}$)")
    journal_axes(ax)
    ax.legend(loc="lower right")
    save(fig, "sigma_LTE")


def fig_kappa(T, k_lte, m, Tmin):
    k = T >= Tmin
    d = k & (T % 500 == 0)
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.plot(T[d], k_lte[d], "o", ms=3.5, mfc="none", mew=0.8, color=INK, label="LTE data")
    ax.plot(T[k], m["k_tot"][k], color=BLUE, label=r"$\kappa_e+\kappa_h+\kappa_R$")
    ax.plot(T[k], m["k_e"][k], "-.", color=MAGENTA, label=r"$\kappa_e$")
    ax.plot(T[k], m["k_h"][k], ":", color=AQUA, lw=1.4, label=r"$\kappa_h$")
    ax.plot(T[k], m["k_r_mpp"][k], "--", color=ORANGE, label=r"$\kappa_R$")
    ax.set_ylim(0, 4)
    ax.set_ylabel(r"$\kappa$ (W m$^{-1}$ K$^{-1}$)")
    journal_axes(ax)
    ax.legend(loc="upper left")
    save(fig, "kappa_LTE")


def main():
    T, k_lte, sig_lte = load_lte()
    with np.errstate(divide="ignore", invalid="ignore"):    # n_e = 0 below ~2000 K
        m = model(T)
    mp = np.loadtxt(MUTPP_LAMBDA, skiprows=1)
    m["k_r_mpp"] = np.interp(T, mp[:, 0], mp[:, 2])
    m["k_tot"] = m["k_e"] + m["k_h"] + m["k_r_mpp"]
    Tmin = 5000.0
    fig_sigma(T, sig_lte, m, Tmin)
    fig_kappa(T, k_lte, m, Tmin)
    ref = np.loadtxt(os.path.join(EXE, "LTE_COND.dat"))
    Tn, ne_ref = ref[:, 0], ref[:, 1]
    ne, _ = composition(Tn)

    print("Figure: kappa_e + kappa_h + kappa_R(Mutation++) / Boulos:",
          " ".join(f"{Tq/1e3:.0f}kK {m['k_tot'][np.argmin(abs(T - Tq))] / k_lte[np.argmin(abs(T - Tq))]:.3f}"
                   for Tq in (6000, 8000, 10000, 12000, 14000, 16000, 18000, 20000, 24000)))
    print("Diagnostics (k_R = the ambipolar estimate, reactive_conductivity.md):")
    print("   T[K]   n_e/n_ref  sigma/LTE  (k_e+k_h)/LTE  (k_e+k_h+k_R)/LTE    k_e    k_h    k_R   k_LTE"
          "  k_M++  kR_M++  k_R/kR_M++")
    for Tq in (6000, 7000, 8000, 10000, 12000, 14000, 15000, 16000, 18000, 20000, 22000, 24000):
        i = np.argmin(abs(T - Tq)); j = np.argmin(abs(Tn - Tq))
        tr = m["k_e"][i] + m["k_h"][i]
        print(f"{T[i]:7.0f}  {ne[j]/ne_ref[j]:9.3f}  {m['sig'][i]/sig_lte[i]:9.3f}  {tr/k_lte[i]:13.3f}"
              f"  {(tr + m['k_r'][i])/k_lte[i]:17.3f}  {m['k_e'][i]:6.3f} {m['k_h'][i]:6.3f}"
              f" {m['k_r'][i]:6.3f} {k_lte[i]:6.3f}  {np.interp(T[i], mp[:, 0], mp[:, 1]):5.3f}"
              f"  {np.interp(T[i], mp[:, 0], mp[:, 2]):6.3f}  {m['k_r'][i] / np.interp(T[i], mp[:, 0], mp[:, 2]):9.3f}")


if __name__ == "__main__":
    main()
