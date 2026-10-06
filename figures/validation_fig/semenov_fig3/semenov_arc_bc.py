#!/usr/bin/env python3
"""
Arc-side boundary state of Semenov, Krivtsun & Reisgen (2016), J. Phys. D 49 105204,
eqs (19)-(24), for the ARGON anode layer -- the counterpart of PlasmeX/exe/almeida_arc_bc.py
(Almeida, Benilov & Naidis 2008, eqs 26-30), whose transport functions it reuses.

At x = L Semenov drops every gradient term and closes the uniform plasma by

    (19)  k_i n_a n_e = k_r n_e^2 n_i          with k_r = k_i (n_a/n_e^2)_Saha  ->  EXACT Saha
    (20)  p = n_a k T_h + n_i k T_h + n_e k T_e,  p = 1 atm
    (21)  n_i = n_e
    (22)  J_e = -(j/e)(delta_i + delta_e)/(gamma_e + gamma_i)
    (23)  E = (j/e^2)(delta_i gamma_e - delta_e gamma_i)/(n_i gamma_e + n_e gamma_i)
    (24)  T_h = T_e,  j E = w_rad(T_e)

(22)-(23) are the two momentum balances (4) with (3), i.e. the same Ohm's law that
almeida_arc_bc.ohm() solves (Almeida eq 27). What differs from Almeida's closure:

  * ionization balance: Semenov defines k_r through the Saha ratio (his section 2.2), so (19)
    is Saha whatever k_i is -- no radiative-escape factor beta = n_e/(n_e + n_e0) as in
    Almeida's (26). Saha weights: G = g+/g_a = 6 (Almeida; Benilov & Naidis 1998 Table I)
    by default; --partition uses the Ar+ internal partition function 4 + 2 exp(-2062 K/T).
  * T_h = T_e exactly (24); Almeida iterates the small T_e - T_h of his (29).
  * w_rad = 2.6e25 p T_e^-2.52 exp(-1.69e5/T_e) [W/m^3], p in bar -- only Almeida's (58),
    which he blends with (60) below n_e0 = 1e21 m^-3 (irrelevant here: n_e > 1e23).
  * ion-atom momentum-transfer cross section: Semenov's own fit (section 2.2)
        Qbar_ia^(1,1) = 1.28 - 0.194 ln(T_h)^0.9   [1e-18 m^2, T_h in kK], 0.5-25 kK,
    read as 1.28 - 0.194 * 0.9 ln(T_h) (the other reading, (ln T_h)^0.9, is not real below
    1 kK, inside the stated range; the two differ by < 1 % in Qbar for 10-20 kK).
    Almeida uses 2 [a - b ln(2 k T_h/eV)]^2; they differ by ~2 % at 14-17 kK. It barely
    matters: the ions carry ~1e-3 of the current.
  * Qbar_ea: Semenov uses Biagi-v8.9 (LXCat); here PlasmeX's ar_Q_en.dat, as in
    almeida_arc_bc.py. At the ionization degrees of Table 2 (~0.5) sigma is dominated by
    electron-ion collisions, so this is a secondary effect (see --sens).

Usage:
    python3 semenov_arc_bc.py                  # Table 2 check, all five current densities
    python3 semenov_arc_bc.py --deck 700 1500  # anode_inputs lines (the Fig. 3 cases)
    python3 semenov_arc_bc.py --j 1.2e7        # any |j| [A/m^2]
    python3 semenov_arc_bc.py --sens           # sensitivity of T to the closure choices

PlasmeX: loaded from the live source tree, $PLASMEX_ROOT or ~/work/PlasmeX (../plasmex.py).
"""
import argparse
import os
import sys

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, ".."))
from plasmex import exe_module  # noqa: E402

abc = exe_module("almeida_arc_bc")

KB, ME, HP, EV_K = abc.KB, abc.ME, abc.H, abc.EV_K
P0 = 101325.0

# Semenov et al. (2016) Table 2: j [A/cm^2] -> (T_e at the anode surface, T_e undisturbed) [K]
TABLE2 = {500: (4.69e3, 13.53e3), 700: (4.93e3, 14.37e3), 1000: (5.20e3, 15.40e3),
          1500: (5.59e3, 16.82e3), 2000: (5.99e3, 18.04e3)}
# domain length for the deck: the paper's plasma is uniform only beyond x ~ 1-3 mm
L_OF = {500: 5e-3, 700: 5e-3, 1000: 4e-3, 1500: 3e-3, 2000: 3e-3}
FIG_OF = {700: "a", 1500: "b"}                 # Table 2 rows that have a Fig. 3 panel


def qbar_ia_semenov(Th):
    """Semenov section 2.2 fit [m^2]; Th in K."""
    return (1.28 - 0.194 * 0.9 * np.log(Th / 1e3)) * 1e-18


def saha(T, partition=False):
    """(n_e n_i/n_a)_Saha at T_e = T [m^-3]."""
    g = 4.0 + 2.0 * np.exp(-2062.0 / T) if partition else abc.G_PLUS / abc.G_1
    return 2 * g * (2 * np.pi * ME * KB * T / HP**2)**1.5 * np.exp(-abc.I_ION * EV_K / T)


def composition(T, p=P0, partition=False):
    """(19)-(21) at T_h = T_e = T: n_e^2 = S n_a, p = (2 n_e + n_a) k T. Returns n_e, n_a."""
    S = saha(T, partition)
    c = p / (KB * T)                           # n_e^2 = S (c - 2 n_e)
    ne = -S + np.sqrt(S * S + S * c)
    return ne, c - 2 * ne


def w_rad(T, p=P0):
    """Semenov / Almeida (58), p in bar."""
    return 2.6e25 * (p / 1e5) / T**2.52 * np.exp(-1.69e5 / T)


def ohm(T, ne, na, j, qia="semenov"):
    """(22)-(23): E [V/m] and the fluxes, via almeida_arc_bc.ohm with Semenov's Qbar_ia."""
    saved = abc.qbar_ia
    if qia == "semenov":
        abc.qbar_ia = qbar_ia_semenov
    try:
        return abc.ohm(T, T, ne, na, j)
    finally:
        abc.qbar_ia = saved


def solve(j, p=P0, partition=False, qia="semenov"):
    """The uniform-plasma state for current density |j| [A/m^2] at pressure p [Pa]."""
    j = abs(j)

    def F(T):
        ne, na = composition(T, p, partition)
        return np.log(j * ohm(T, ne, na, j, qia)[0] / w_rad(T, p))

    T = brentq(F, 5e3, 4e4, xtol=1e-6)
    ne, na = composition(T, p, partition)
    E, v_i, v_e, _, _ = ohm(T, ne, na, j, qia)
    return dict(j=j, p=p, T=T, ne=ne, na=na, E=E, sigma=j / E, ji_frac=abc.E_CH * ne * v_i / j,
                Je=ne * v_e, w=w_rad(T, p), ion_deg=ne / (ne + na))


def fig3_far_field(j_Acm2):
    """n_e and T_e at x = 8 mm on the extracted Fig. 3 curve, or None."""
    if j_Acm2 not in FIG_OF:
        return None
    import fig3_extract
    d = fig3_extract.load(FIG_OF[j_Acm2])
    at = lambda q: float(np.interp(np.log(8e-3), np.log(d[q][0]), d[q][1]))
    return at("ne"), at("Te")


def table2(partition=False):
    print(f"Semenov et al. (2016) eqs (19)-(24), argon, p = {P0:g} Pa"
          f"{', Ar+ partition function' if partition else ', g+/g_a = 6'}.")
    print("T 'paper' = Table 2, T_e in the undisturbed plasma; n_e 'Fig. 3' = extracted "
          "curve at x = 8 mm.\n")
    print(f"{'j[A/cm2]':>8} {'T paper':>8} {'T here':>8} {'dT/T':>7}   {'n_e here':>9} "
          f"{'n_e Fig3':>9} {'dn/n':>7}   {'E[V/m]':>7} {'sigma[S/m]':>10} {'x_ion':>6} "
          f"{'j_i/j':>8}")
    for jA, (_, Tp) in TABLE2.items():
        r = solve(jA * 1e4, partition=partition)
        ff = fig3_far_field(jA)
        nf = (f"{ff[0]:9.3e} {r['ne'] / ff[0] - 1:+7.1%}" if ff else f"{'-':>9} {'-':>7}")
        print(f"{jA:8d} {Tp:8.0f} {r['T']:8.0f} {(r['T'] - Tp) / Tp:+7.2%}   {r['ne']:9.3e} "
              f"{nf}   {r['E']:7.0f} {r['sigma']:10.4g} {r['ion_deg']:6.3f} {r['ji_frac']:8.1e}")


def sensitivity():
    """How T responds to each closure choice, at the two Fig. 3 current densities."""
    print("Sensitivity of the undisturbed T to the closure (Table 2 T in brackets):\n")
    cases = [("this script (Semenov closure)", {}),
             ("Ar+ partition function instead of G = 6", dict(partition=True)),
             ("Almeida's Qbar_ia instead of Semenov's", dict(qia="almeida"))]
    for jA in (700, 1500):
        print(f"j = {jA} A/cm^2 [{TABLE2[jA][1]:.0f} K]")
        for label, kw in cases:
            r = solve(jA * 1e4, **kw)
            print(f"  {label:42s} T = {r['T']:8.1f} K  ({r['T'] / TABLE2[jA][1] - 1:+.2%})")
        ra = abc.solve(jA * 1e4, P0)
        print(f"  {'almeida_arc_bc.solve (beta, w_rad blend, (29))':42s} T = {ra['Te']:8.1f} K"
              f"  ({ra['Te'] / TABLE2[jA][1] - 1:+.2%})")
        # T that would be needed to match Table 2 exactly, and what E the paper implies
        T = TABLE2[jA][1]
        print(f"  at T = {T:.0f} K: w_rad/j = {w_rad(T) / (jA * 1e4):.0f} V/m is the E Semenov's "
              f"column must have; the Ohm's law here gives "
              f"{ohm(T, *composition(T), jA * 1e4)[0]:.0f} V/m\n")


def deck(j_Acm2, partition=False, use_table2=True):
    """anode_inputs lines. use_table2: set T from Table 2 (the paper's own bulk) rather than
    from the radiation balance solved here; n_e and E are then evaluated at that T."""
    j = j_Acm2 * 1e4
    r = solve(j, partition=partition)
    T = TABLE2[j_Acm2][1] if use_table2 else r["T"]
    ne, na = composition(T, P0, partition)
    E_ohm = ohm(T, ne, na, j)[0]
    # At the paper's own T its column satisfies (24) exactly, so its field is w_rad(T)/j --
    # no conductivity needed. At the solved T the two coincide by construction.
    E = w_rad(T) / j if use_table2 else E_ohm
    print(f"# --- Semenov et al. (2016), j = {j_Acm2} A/cm^2 = {j:.3g} A/m^2; T from "
          f"{'Table 2' if use_table2 else 'jE = w_rad solved here'} ---")
    print(f"geometry.prob_hi       = {L_OF[j_Acm2]:g}  1e-4  1e-4   # >= 3 x the layer; "
          f"check U_a does not change when doubled")
    print("prob.pres              = 101325")
    print("prob.wall_temp         = 1000")
    print("prob.wall_temp_balance = 0")
    print("prob.electrode_emission = 0")
    print(f"prob.arc_Japp          = {-j:.4g}             # anode: same sign as anode_inputs")
    print(f"prob.arc_Te            = {T:.0f}")
    print(f"prob.arc_e             = {ne:.4e}")
    print(f"prob.arc_arp           = {ne:.4e}")
    print(f"prob.arc_efldx         = {-E:.4g}             # same sign as arc_Japp (j = sigma E)")
    print(f"# implied n_Ar = {na:.4e} m^-3; the radiation balance here gives T = {r['T']:.0f} K "
          f"({r['T'] / TABLE2[j_Acm2][1] - 1:+.2%} vs Table 2)")
    if use_table2:
        print(f"# arc_efldx = w_rad(T)/j (the paper's eq. 24 at its T); the Ohm's law here "
              f"(Almeida App. A coefficients) gives {E_ohm:.4g} V/m ({E_ohm / E - 1:+.1%})")
    ff = fig3_far_field(j_Acm2)
    if ff:
        print(f"# Fig. 3({FIG_OF[j_Acm2]}) at x = 8 mm: n_e = {ff[0]:.4e} m^-3 "
              f"({ne / ff[0] - 1:+.1%} vs arc_e above), T_e = {ff[1]:.0f} K")
    print()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--j", type=float, nargs="+", help="|j| [A/m^2]: print the state")
    ap.add_argument("--deck", type=int, nargs="*", metavar="J_Acm2",
                    help="anode_inputs lines for Table 2 rows (default 700 1500)")
    ap.add_argument("--solved-T", action="store_true",
                    help="with --deck: use the T solved here instead of Table 2")
    ap.add_argument("--partition", action="store_true",
                    help="Saha with the Ar+ partition function instead of G = 6")
    ap.add_argument("--sens", action="store_true", help="sensitivity to the closure choices")
    a = ap.parse_args()
    if a.sens:
        sensitivity()
    elif a.deck is not None:
        for jA in a.deck or [700, 1500]:
            if jA not in TABLE2:
                ap.error(f"{jA} A/cm^2 is not a Table 2 row ({sorted(TABLE2)})")
            deck(jA, a.partition, not a.solved_T)
    elif a.j:
        for jj in a.j:
            r = solve(jj, partition=a.partition)
            print(f"|j| = {r['j']:.4g} A/m^2: T = {r['T']:.1f} K, n_e = {r['ne']:.4e}, "
                  f"n_a = {r['na']:.4e} m^-3, E = {r['E']:.1f} V/m, sigma = {r['sigma']:.4g} S/m")
    else:
        table2(a.partition)


if __name__ == "__main__":
    main()
