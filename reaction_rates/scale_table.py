#!/usr/bin/env python3
"""
LTE length and time scales quoted in Table tab:scales (sections/model.tex).

At each T the plasma is in LTE at p = 1 atm (T_e = T_h = T). The composition is Ar, Ar+ and electrons
from the Saha equation, quasi-neutrality n_e = n_i and Dalton's law, with internal partition functions
summed over the low-lying levels below (NIST ASD energies, no lowering of the ionization potential).
  tau_rxn    = 1/(2 k_rec n_e^2 + 2 k_ion n_e)           linearised ionization-recombination, fixed p, T
  lambda_rxn = sqrt(D_a tau_rxn),  D_a = (1 + T_e/T_h) D_in = 2 D_in   (D_in: Almeida (45), n_h = n_n + n_i)
  lambda_e   = 1/(n_n Qbar_en + n_i Qbar_ei)              Qbar_en: Almeida (46) on exe/ar_Q_en.dat; Qbar_ei: (47)
  lambda_T   = lambda_e / sqrt(2 m_e/M)
  lambda_Th  = 1/[(n_i + n_n) Qbar_in]                    ion-neutral temperature equilibration;
                                                          Qbar_in: Benilov & Naidis resonant charge exchange
  lambda_i   = 1/(n_n Qbar_in + n_i Qbar_ii)              ion mean free path; Qbar_ii = (47) at T_h
  lambda_De  = sqrt(eps0 k T/(n_e e^2))
Rates are Hoffert & Lien (exe/hoff_f.dat, hoff_r.dat), as in the model. Everything is loaded from
the live PlasmeX tree through figures/validation_fig/plasmex.py.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "figures", "validation_fig"))
from plasmex import exe_module, EXE  # noqa: E402

A = exe_module("almeida_arc_bc")
KB, E, ME, M, EPS0, EV_K, H = A.KB, A.E_CH, A.ME, A.M_AR, A.EPS0, A.EV_K, A.H
P0 = 101325.0

# ---------------------------------------------------------------- Saha (Ar, Ar+)
I1 = 15.7596                                  # [eV] ionization energy of Ar (NIST)
# (energy [eV], statistical weight) of the low-lying levels (NIST ASD, grouped by configuration)
LEV_AR = [(0.0, 1), (11.548, 5), (11.624, 3), (11.723, 1), (11.828, 3),   # 3p6, 3p5 4s
          (13.17, 36), (14.0, 60), (14.15, 12), (14.6, 36)]               # 4p, 3d, 5s, 5p
LEV_ARP = [(0.0, 4), (0.178, 2), (13.48, 2)]                              # 3p5 2P, 3s3p6 2S


def Z(levels, T):
    kT = T / EV_K
    return sum(g * np.exp(-e / kT) for e, g in levels)


def composition(T, p=P0):
    """n_e = n_i, n_n in LTE: n_e^2/n_n = S, n_n + 2 n_e = p/kT"""
    S = 2 * Z(LEV_ARP, T) / Z(LEV_AR, T) * (2 * np.pi * ME * KB * T / H**2) ** 1.5 \
        * np.exp(-I1 * EV_K / T)
    N = p / (KB * T)
    ne = 0.5 * (-2 * S + np.sqrt(4 * S**2 + 4 * S * N))
    return ne, N - 2 * ne


# ---------------------------------------------------------------- rates and scales
_f = np.loadtxt(os.path.join(EXE, "hoff_f.dat"))
_r = np.loadtxt(os.path.join(EXE, "hoff_r.dat"))


def hl(T):
    x = T / EV_K
    return (np.exp(np.interp(x, _f[:, 0], np.log(_f[:, 1]))),
            np.exp(np.interp(x, _r[:, 0], np.log(_r[:, 1]))))


def scales(T):
    ne, nn = composition(T)
    ni = ne
    nh = nn + ni
    ki, kr = hl(T)
    tau = 1.0 / (2 * kr * ne**2 + 2 * ki * ne)
    Da = 2 * A.nD(T, M / 2, A.qbar_ia(T)) / nh
    Qen = A.qbar_ea(T)[0]
    Qei = A.qbar_ei(T, ne)
    Qin = A.qbar_ia(T)
    lam_e = 1.0 / (nn * Qen + ni * Qei)
    Qii = A.qbar_ei(T, ne)                     # (47) with T_h = T_e = T
    return dict(T=T, ne=ne, nn=nn, lnL=A.ln_lambda(T, ne),
                bal=ki * nn / (kr * ne**2), tau=tau, Da=Da,
                lDe=np.sqrt(EPS0 * KB * T / (ne * E**2)), le=lam_e, li=1.0 / (nn * Qin + ni * Qii),
                lTh=1.0 / ((ni + nn) * Qin),
                lrxn=np.sqrt(Da * tau), lT=lam_e / np.sqrt(2 * ME / M))


if __name__ == "__main__":
    keys = ["T", "ne", "nn", "lnL", "bal", "Da",
            "lDe", "le", "li", "lTh", "lrxn", "lT", "tau"]
    print(" ".join(f"{k:>9s}" for k in keys))
    for T in (12e3, 14e3, 18e3, 20e3, 22e3):
        s = scales(T)
        print(" ".join(f"{s[k]:9.3g}" for k in keys))
