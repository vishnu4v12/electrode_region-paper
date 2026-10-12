"""
Ar+ mobility in Ar of Appendix A, with and without the two-temperature (Wannier) correction,
compared with the Monte Carlo reduced mobility of Nikitovic et al. (Sci. Sintering 58 (2026) 409, Fig. 7).

Model (eqs. Qia, mu_i0, wannier, mu-in), evaluated with PlasmeX's own bohm.ion_transport:
    mu_i^0 n_n depends on T_h only, and v_d, T_eff depend on E/n_n only, so the reduced mobility
    K_0 = mu n_n / N_0 (N_0 = 2.6868e25 m^-3) is a function of E/N at fixed gas temperature.
      "without"  K_0 = mu_i^0 n_n / N_0                     (constant)
      "with"     K_0 = mu_in  n_n / N_0 = mu_i^0 sqrt(T_h/T_eff) n_n / N_0

Reference: the Monte Carlo (bulk) K_0(E/N) of Nikitovic et al. at T = 300 K, digitized from their
Fig. 7 into nikitovic2026_fig7_MC.dat (see its header). Their MC agrees with the drift-tube data
(Helm & Elford, Ellis et al., Hegerberg et al., Phelps, Basurto et al., Dashdorj et al.) shown in
that figure. The model curves are drawn at T_h = 300 K.

Nothing from PlasmeX is copied (validation_fig/plasmex.py).   Run:  python3 ion_mobility.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "validation_fig"))
from plasmex import exe_module  # noqa: E402

b = exe_module("bohm")
N0 = 2.6867811e25                       # Loschmidt constant, 273.15 K and 101.325 kPa
TD = 1e-21                              # 1 Td in V m^2
DATA = os.path.join(HERE, "nikitovic2026_fig7_MC.dat")
T_MC = 300.0                            # gas temperature of the MC data

# same journal style and palette as figures/lte_check/lte_check.py
INK, BLUE, ORANGE, AQUA, MAGENTA = "#1f1f1e", "#2a78d6", "#eb6834", "#1baf7a", "#e87ba4"
MARKERS = ["o", "s", "^", "D", "v"]
FIGSIZE = (3.4, 2.6)
plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"], "mathtext.fontset": "cm",
                     "font.size": 8, "axes.linewidth": 0.6, "lines.linewidth": 1.1,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
                     "legend.frameon": False, "legend.fontsize": 7.5, "legend.handlelength": 2.4})


def model_k0(en_td, t_gas):
    """K0 [cm^2/Vs] without and with the Wannier correction, at E/N [Td] and T_h = t_gas."""
    n = N0                                           # any density: K0 depends on E/N only
    mu0 = b.ion_mobility0(t_gas, n)
    _, mu_in, _ = b.ion_transport(t_gas, n, en_td * TD * n)
    return np.full_like(en_td, mu0 * n / N0 * 1e4), mu_in * n / N0 * 1e4


def main():
    en_mc, k_mc = np.loadtxt(DATA, unpack=True)
    en = np.logspace(0, 4, 300)
    k_no, k_w = model_k0(en, T_MC)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.semilogx(en_mc, k_mc, "-", color=INK, label="Monte Carlo, Nikitović et al.")
    ax.semilogx(en, k_no, "--", color=ORANGE, label=r"$\mu_i^0$ (single temperature)")
    ax.semilogx(en, k_w, color=BLUE, label=r"$\mu_{in}$ (two-temperature)")
    ax.tick_params(direction="in", which="both", top=True, right=True)
    ax.set_xlim(en[0], en[-1])
    ax.set_ylim(0, 1.8)
    ax.set_xlabel(r"$E/N$ (Td)")
    ax.set_ylabel(r"$K_0$ (cm$^2$ V$^{-1}$ s$^{-1}$)")
    ax.legend(loc="lower left")
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, f"mobility_Ar.{ext}"), bbox_inches="tight", dpi=300)
    plt.close(fig)

    print(f"model at T_h = {T_MC:.0f} K: K0(single temperature) = {k_no[0]:.3f} cm^2/Vs")
    print("  E/N [Td]   K0 MC   K0 two-temperature   K0 single temperature")
    for q in (1, 10, 30, 100, 300, 1000, 3000, 7000):
        lq = np.log10(q)
        print(f"  {q:7d}   {np.interp(lq, np.log10(en_mc), k_mc):.3f}   "
              f"{np.interp(q, en, k_w):.3f}               {k_no[0]:.3f}")


if __name__ == "__main__":
    main()
