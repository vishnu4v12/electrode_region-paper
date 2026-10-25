"""
Numbers quoted in sections/results.tex (sec:res-structure, sec:res-energy, sec:res-heat) that
base_case.py does not print: profile values at fixed distances, domain and sheath integrals of the
electron energy terms, and the components of the electron energy flux.

    python text_numbers.py <anode plotfile> <cathode plotfile>   (> text_numbers_report.txt)

Same leaf-cell route as base_case.py (base_case.load); the energy-flux components are
bohm.electron_energy_flux (ENTH, COND, TD on the interior faces; TOT; q_wall at the wall face).
Distances are from the electrode; fluxes are reported positive TOWARDS the electrode.
"""
import sys

import numpy as np

import base_case as bc
from base_case import K_B, E_CH

b = bc.exe_module("bohm")


def at(d, f, x):
    return float(np.interp(np.log(x), np.log(d), f))


def flux_parts(pf, s):
    dd = b.sample(pf, bc.FIELDS, want_derived=False)
    knobs, _ = b.read_knobs(dd["dir"], None, pf)
    q = b.electron_energy_flux(dd, knobs, bc.EXE)
    hi = dd["hi"]
    df = hi - q["xf"]
    o = np.argsort(df)
    parts = {k: np.asarray(q[k])[o] for k in ("TOT", "ENTH", "COND", "TD")}
    return df[o], parts, float(q["q_wall"])


def integ(s, term, lo=0.0, hi=np.inf):
    d, dx = s["d"], s["dxc"]
    v = s["budget"][term]
    m = (d >= lo) & (d <= hi) & np.isfinite(v)
    return float(np.sum(v[m] * dx[m]))


def shares(v, step=1.0):
    """Percentages of sum(v), rounded to `step` by largest remainder so that they sum to exactly 100."""
    p = 100.0 * np.asarray(v) / np.sum(v) / step
    f = np.floor(p)
    k = int(round(100.0 / step - f.sum()))
    f[np.argsort(f - p)[:k]] += 1
    return f * step


def report(tag, pf):
    s = bc.load(pf)
    d, ne, ni, nn, Te, Th, phi, E = (s[k] for k in ("d", "ne", "ni", "nn", "Te", "Th", "phi", "E"))
    x_s, phi_s, _ = bc.sheath_edge(d, ne, ni, phi)
    info = s["info"]
    g = lambda k: float(info.get(k, np.nan))
    print(f"\n===== {tag}: {pf}")
    print(f"x_s {x_s:.3e}  phi_s {phi_s:.3f}  phi(L) {phi[-1]:.3f}  T_wall {g('prob.wall_temp'):.0f}")
    print(f"surface cell: n_e {ne[0]:.3e} n_i {ni[0]:.3e} (n_e/n_i-1 {ne[0]/ni[0]-1:+.3f}, n_i/n_e {ni[0]/ne[0]:.2f})"
          f"  n_n {nn[0]:.3e}  T_e {Te[0]:.0f}  T_h {Th[0]:.0f}  T_eff {s['Teff'][0]:.0f}  E_face {g('q_wall.E_face'):.3e}")
    print(f"max n_i/n_e in sheath {np.max((ni/ne)[d <= x_s]):.2f}")
    for x in (1e-7, 1e-6, x_s, 1e-5, 1.4e-5, 1e-4, 3e-4, 5e-4 * 0.999):
        print(f"  x={x:.2e}: phi {at(d, phi, x):.3f}  n_e {np.exp(at(d, np.log(ne), x)):.3e}"
              f"  n_n {np.exp(at(d, np.log(nn), x)):.3e}  T_e {at(d, Te, x):.0f}  T_h {at(d, Th, x):.0f}"
              f"  T_e-T_h {at(d, Te - Th, x):.0f}  E {at(d, E, x):.3e}  j_i/j {at(d, s['ji_j'], x):+.4f}")
    print(f"  bulk n_e {ne[-1]:.3e}; bulk/surface n_e {ne[-1]/ne[0]:.0f};  T_e range {Te.min():.0f}-{Te.max():.0f}")
    print(f"  column field at L (last 10 um mean) {np.mean(E[d > d[-1] - 1e-5]):.3e} V/m")
    gem, gab, gi = g("q_wall.e_flux_emitted"), g("q_wall.e_flux_absorbed"), g("q_wall.ion_flux")
    J = abs(np.median(s["j"]))
    print(f"  Gamma_em {gem:.3e} -> j_em {E_CH*gem:.3e} A/m^2 = {E_CH*gem/J:.3f} J;  Gamma_abs {gab:.3e};  Gamma_i {gi:.3e}")
    print(f"  j_i/j surface {s['ji_j'][0]:.4f}")

    # electron energy budget
    MW = 1e-6
    for term in ("JOULE", "ELASTIC", "ION", "RAD", "RADL"):
        print(f"  int {term:7s}: domain {MW*integ(s, term):+9.2f}  sheath {MW*integ(s, term, 0, x_s):+9.2f}"
              f"  x<10um {MW*integ(s, term, 0, 1e-5):+9.2f}  x<100um {MW*integ(s, term, 0, 1e-4):+9.2f}"
              f"  x>100um {MW*integ(s, term, 1e-4):+9.2f}  MW/m^2")
    i = int(np.nanargmax(np.abs(s["budget"]["JOULE"])))
    print(f"  Joule peak {s['budget']['JOULE'][i]:.3e} W/m^3 at {d[i]:.2e}")
    df, P, qw = flux_parts(pf, s)
    # flux towards the electrode: lab +x is towards the electrode (electrode at hi), so TOT etc. as is
    print(f"  h_e at the wall face (q_wall eeng) {MW*qw:+.2f} MW/m^2")
    for x in (1e-7, x_s, 1e-6, 1e-5, 1e-4, df[-1]):
        k = int(np.argmin(np.abs(df - x)))
        print(f"  h_e at x={df[k]:.2e}: TOT {MW*P['TOT'][k]:+9.2f}  ENTH {MW*P['ENTH'][k]:+9.2f}"
              f"  COND {MW*P['COND'][k]:+9.2f}  TD {MW*P['TD'][k]:+9.2f}  MW/m^2 (towards electrode)")
    # Shares of the domain energy budget quoted in sec:res-energy. Each side is normalised by its own
    # total and rounded to 1 % by largest remainder, so the quoted percentages sum to exactly 100.
    h_L, h_w = P["TOT"][-1], qw                  # towards the electrode, at x ~ L and at the wall face
    I = {t: integ(s, t) for t in ("JOULE", "ELASTIC", "ION", "RAD", "RADL")}
    if tag == "Cathode":
        sides = {"Joule heat goes to": {"h_e into the column at L": -h_L, "h_e into the cathode": h_w,
                                        "radial loss": -I["RADL"], "ionization": -I["ION"],
                                        "radiation": -I["RAD"], "elastic": -I["ELASTIC"]}}
        gain = I["JOULE"]
    else:
        sides = {"energy supplied by": {"-div h_e (h_L - h_w)": h_L - h_w, "recombination": I["ION"]},
                 "energy taken by": {"Joule loss": -I["JOULE"], "elastic": -I["ELASTIC"],
                                     "radial loss": -I["RADL"], "radiation": -I["RAD"]}}
        gain = h_L - h_w + I["ION"]
    for name, terms in sides.items():
        tot = sum(terms.values())
        pct = shares(list(terms.values()))
        print(f"  {name} (total {MW*tot:.2f} MW/m^2, {100*(tot/gain - 1):+.3f} % vs the gains): "
              + ", ".join(f"{k} {p:.0f} % ({100 * v / tot:.2f})" for (k, v), p in zip(terms.items(), pct)) + f"  [sum {sum(pct):.0f}]")
    print(f"  h_e reaching the surface / h_e entering at L: {h_w / h_L:.3f}")
    # conduction sign change
    c = P["COND"]
    sc = np.nonzero(np.diff(np.sign(c)))[0]
    print(f"  COND sign changes at x = " + ", ".join(f"{df[k]:.2e}" for k in sc[:6]))
    S = s["budget"]["ION"]
    neg = d[(S > 0) & np.isfinite(S)]
    print(f"  ION term > 0 (recombination heating) at x in [{neg.min() if neg.size else np.nan:.2e}, {neg.max() if neg.size else np.nan:.2e}]"
          f"  fraction of cells {neg.size / d.size:.2f}")
    for k in ("electron_kinetic", "electron_work_function", "emission_kinetic", "emission_work_function",
              "ion_recombination", "gas_conduction", "total"):
        print(f"  q_wall.{k:24s} {MW*g('q_wall.' + k):9.3f}")
    heat = sum(g("q_wall." + k) for k in ("electron_kinetic", "electron_work_function", "ion_recombination", "gas_conduction"))
    cool = g("q_wall.emission_kinetic") + g("q_wall.emission_work_function")
    print(f"  heating {MW*heat:.1f}  cooling {MW*cool:.1f}  absorbed electrons {MW*(g('q_wall.electron_kinetic')+g('q_wall.electron_work_function')):.1f}"
          f"  ion share of heating {g('q_wall.ion_recombination')/heat:.3f}  electron share of q_w {(g('q_wall.electron_kinetic')+g('q_wall.electron_work_function'))/g('q_wall.total'):.3f}")


if __name__ == "__main__":
    report("Anode", sys.argv[1])
    report("Cathode", sys.argv[2])
