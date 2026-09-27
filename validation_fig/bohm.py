"""
1-D sheath plots from a PlasmeX plotfile.

Four stacked panels against the axial coordinate:

    1  n_e, n_i and the Boltzmann prediction for n_e          (log y)
    2  the potential
    3  Te, T_gas and the ion effective temperature T_eff
    4  the ion fluid velocity against the Bohm speed

plus a terminal report: the boundary-cell (edge) values, the sheath fall, the Bohm velocity,
the E field at the arc edge, the drift/diffusion split of the ion and electron flux there,
the min/max of every plotted field inside the plotted window, and the heat flux to the
electrode.

    python bohm.py plt00200
    python bohm.py 'plt*'

Everything else is configured in the USER CONFIG block below. This module is also the shared
library for transport.py and gif.py -- the sampler, the LTE T_gas inverse, the Wannier ion
closure and the electron transport rebuild live here once, so the three scripts cannot drift
apart from each other or from the solver.

--------------------------------------------------------------------------------------------
EDGE VALUES ARE READ, NEVER EXTRAPOLATED
--------------------------------------------------------------------------------------------
Every "at the arc" / "at the electrode" number below is the value stored in the FIRST or LAST
LEAF CELL of the plotfile, verbatim. Nothing is extrapolated to a face.

That is the right thing to do because the solver has already put the face value where it
matters. `BoundaryFaceOverwrite` (src/aider.H) overwrites the boundary cell of the three
conduction-current fields with the domain FACE value, so JEX/JIX/JX in the wall cell ARE the
flux the solver applied across the wall, not a cell average. Every other field is an ordinary
cell-centred quantity there:

    JX, JEX, JIX          boundary cell = domain FACE value        (FACE_CELL_FIELDS)
    everything else       boundary cell = ordinary cell value

The report tags each edge quantity `[FACE]` or `[cell]` so it is never ambiguous which one is
being quoted. FACE_CELL_FIELDS below is the single place that list is written down; if
src/aider.H gains or drops a BoundaryFaceOverwrite call, change it there.

A consequence worth stating: bc.H's Schottky barrier lowering reads the boundary FACE field
`efld_face`, but EFLDX is not face-overwritten, so the wall-cell EFLDX this script quotes is
a cell average and the Jem it produces is NOT bc.H's Jem. It is reported as the cell-value
estimate it is rather than papered over with an extrapolation.

--------------------------------------------------------------------------------------------
LEAF CELLS, ANY REFINEMENT RATIO
--------------------------------------------------------------------------------------------
sample() returns the AMR LEAF cells -- every cell the solver actually solved on, at its own
width, in x order -- and never a uniform fine resampling. So amr.ref_ratio is read off the
grids themselves and may be any value, or a different value per level; nothing here computes
ref_ratio**(finest-level), an assumption that was wrong and silently so on this deck (which
ships ref_ratio = 4).

It is also what makes the scripts fast: the shipped 4-level deck has ~640 leaf cells against
32768 fine samples, so every gradient, every reconstruction and every matplotlib call below
does ~50x less work, and a coarse cell is one point rather than 64 identical copies.

Only two of the plotted quantities are not in the plotfile:

    T_gas   the neutral density AR inverted through the LTE (Saha) map, which is what set AR
            in the first place -- argon_lte::n_neutral in exe/chemistry.H. NOT p/(k_B T),
            which is the TOTAL density and reads 152x high at a 21000 K arc.
    T_eff   the ion effective temperature of the Wannier two-temperature closure,
            k T_eff = k T_gas + (M/3) v_d^2, as src/aider.H::CalcTransportCoeffs forms it.
            It is an intermediate there and never reaches the plotfile.

1-D ONLY, by assertion. Use gif_2d.py for r-z.
"""

# ============================================================
# USER CONFIG
# ============================================================

# ---- x axis ------------------------------------------------------------------------
# Plot window in metres. None takes the domain edge from the plotfile.
PROB_LO = None
PROB_HI = None

# Log x axis. A plain log10(x) would collapse the sheath into the last decade -- the
# opposite of useful -- so this plots against the DISTANCE FROM THE ELECTRODE,
# w = x_domain_hi - x, on a log scale, with the axis inverted so the electrode stays on
# the right and the arc on the left exactly as in the linear view. PROB_LO/PROB_HI are
# always given in real x either way.
XLOG = False

# ---- deck ---------------------------------------------------------------------------
# The input deck to read prob.pres / prob.wall_temp / prob.electrode_emission from, so the
# heat-flux block cannot quietly disagree with the run it is describing. Searched next to the
# plotfile first, then next to this script. Set to None to use the literals below instead.
DECK_FILE = "inputs"

# Fallbacks, used only when DECK_FILE is None or does not carry the parameter.
PRESSURE = 101325.0      # prob.pres [Pa]; rebuilds T_gas (and hence T_eff) from AR
WALL_TEMP = 3500.0       # prob.wall_temp [K]
EMISSION = True          # prob.electrode_emission

# ---- sheath edge -------------------------------------------------------------------
# Charge-separation threshold defining the sheath edge, walking in from the wall:
#     |n_i - n_e| / n_i  >  QN_EPS   ->  inside the sheath
QN_EPS = 0.01

# Where the ARC (bulk) reference is taken, as a fraction of the domain measured from the
# arc. 0.0 = the arc boundary cell itself, which is where the deck imposes the arc state and
# is as far from the sheath as the domain goes.
ARC_REF_FRAC = 0.0

# ---- electron transport tables ------------------------------------------------------
# Directory holding ar_e_lorentz.dat. Needed only for the arc-side drift/diffusion split
# and for transport.py -- the plotfile carries no mobility. None searches the plotfile's
# directory and then this script's.
TABLE_DIR = None

# ---- electrode heat flux -- these must match exe/bc.H, not be tuned ------------------
WORK_FN = 2.63           # prob.work_function [eV] fallback, the code default (exe/prob.H)
RICHARDSON_A = 3.2e4     # prob.richardson_A [A m^-2 K^-2] fallback, the code default
NUEI_C  = 2.9063e-12     # Mitchner & Kruger (VIII 5.1d) constant, Te in eV: src/coulomb.H
E_ION = 15.7596          # Ar first ionization potential [eV], argon_lte::EION_J

# ---- output ------------------------------------------------------------------------
OUT_PREFIX = "bohm_"
DPI = 150


# ============================================================
# PHYSICAL CONSTANTS  (must match exe/chemistry.H)
# ============================================================

K_B     = 1.380649e-23
ECHARGE = 1.602176634e-19
ARMASS  = 6.63e-26
EPS0    = 8.8541878128e-12
EMASS   = 9.10938356e-31

# Saha constants, mirroring argon_lte in exe/chemistry.H exactly.
SAHA_PREFAC = 4.829365966618981e21     # 2 (2 pi m_e k_B / h^2)^{3/2}
EION_J      = 2.524966288118640e-18    # Ar first ionization potential, 15.7596 eV
TFIT_LO     = 1000.0
TFIT_HI     = 25000.0
# g(u) = 2 Q_Ar+ / Q_Ar, u = T/1e4, highest order first
G_COEFF = [2.74894191130579e-02, -4.25604572769701e-01,  2.82670401706654e+00,
          -1.04950971422059e+01,  2.40084948769187e+01, -3.55184984989637e+01,
           3.46901606404912e+01, -2.23756159163732e+01,  9.38724110129288e+00,
           3.50219221604479e+00]

# aider.H's reference density for the low-field reduced ion mobility.
N0_REF = 2.686763e25

# plasmex_coulomb::TE_FLOOR_EV (src/coulomb.H): below this the correction is switched off.
# It guards the divisions by Te; it is not a physical cutoff.
TE_FLOOR_EV = 1.0e-3

# 1-D orientation: the electrode is at x-hi, so ELEC_SGN = +1 (src/var_def.H). Every
# "towards the wall" quantity below is sgn-oriented with this.
SGN = 1.0

# Fields whose BOUNDARY CELL was overwritten with the domain FACE value by
# BoundaryFaceOverwrite (src/aider.H). See the module docstring. Everything not listed here
# holds an ordinary cell-centred value in the boundary cell.
FACE_CELL_FIELDS = frozenset(("JX", "JEX", "JIX", "JY", "JEY", "JIY", "UNX", "UNY"))

# The fields every consumer needs. sample() reads exactly what the caller asks for and
# nothing else -- reading all 16 components of a plotfile to plot four of them is the single
# largest avoidable cost in these scripts, and gif.py in particular names only four.
BASE_FIELDS = ("NE", "ARP", "AR", "PHI", "ETEMP", "TGAS", "EFLDX")

# What sample() reads when the caller names nothing: everything this module's own analyse() /
# electrode_heat_flux() / report() can use. Anything absent from the plotfile is skipped, so
# older plotfiles that predate JEX/JIX/EEN still load. steady.py calls sample(fn) bare and
# relies on this covering ION_VEL.
DEFAULT_FIELDS = BASE_FIELDS + ("ION_VEL", "JEX", "JIX", "JX", "EEN", "UNX")


# ============================================================
# IMPORTS
# ============================================================

import argparse
import glob
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yt

yt.set_log_level("error")


# ============================================================
# HELPERS
# ============================================================

def extract_number(fname):
    m = re.search(r'(\d+)$', os.path.basename(fname))
    return int(m.group(1)) if m else -1


def face_tag(name):
    """`[FACE]` or `[cell]` -- what the boundary cell of this field actually holds."""
    return "[FACE]" if name in FACE_CELL_FIELDS else "[cell]"


def search_dirs(start):
    """The two places a data file may live: beside the plotfile, then beside this script."""
    here = os.path.dirname(os.path.abspath(__file__))
    return [start, here] if start and start != here else [here]


_DECK_CACHE = {}


def deck_values(exe_dir, deck=None):
    """
    Every `name = value` of the input deck `deck` (default DECK_FILE), as a dict of the FIRST
    whitespace token of each value, plus the path it was read from -- or ({}, None) when the
    deck is not found (or `deck` is None). A minimal ParmParse reader: `#` comments, last
    assignment wins. Searched next to the plotfile first, then next to this script. Cached.
    """
    deck = DECK_FILE if deck is None else deck
    if deck is None:
        return {}, None

    key = (exe_dir, deck)
    if key in _DECK_CACHE:
        return _DECK_CACHE[key]

    path = next((p for p in (os.path.join(d, deck) for d in search_dirs(exe_dir))
                 if os.path.isfile(p)), None)
    vals = {}
    if path is not None:
        with open(path) as fh:
            for line in fh:
                line = line.split("#", 1)[0]
                if "=" not in line:
                    continue
                k, v = line.split("=", 1)
                v = v.split()
                if v:
                    vals[k.strip()] = v[0]
    _DECK_CACHE[key] = (vals, path)
    return vals, path


def read_deck(exe_dir):
    """
    Pull prob.pres / prob.wall_temp / prob.electrode_emission out of the input deck.

    The point is that the heat-flux block cannot silently describe a different run from the
    one that wrote the plotfile -- the shipped deck has moved wall_temp and electrode_emission
    more than once, and a stale literal here is wrong output, not a missing feature.

    Returns (pressure, wall_temp, emission, source) and falls back to the USER CONFIG
    literals for anything the deck does not carry.
    """
    if DECK_FILE is None:
        return PRESSURE, WALL_TEMP, EMISSION, "USER CONFIG (DECK_FILE = None)"

    vals, path = deck_values(exe_dir)
    if path is None:
        return PRESSURE, WALL_TEMP, EMISSION, f"USER CONFIG ({DECK_FILE} not found)"

    def get(name, default, cast=float):
        try:
            return cast(vals[name])
        except (KeyError, ValueError):
            return default

    return (get("prob.pres", PRESSURE),
            get("prob.wall_temp", WALL_TEMP),
            bool(get("prob.electrode_emission", EMISSION, lambda s: int(float(s)))),
            path)


# The run-time knobs electron_energy_flux() needs, with the CODE defaults (src/PlasmeX.cpp /
# CLAUDE.md: all three physics switches default OFF in the code, ON in the decks). A deck that
# does not set one ran with the code default, so that -- not the shipped-deck value -- is the
# right fallback.
ENERGY_FLUX_KNOBS = {
    "plasmex.thermal_diffusion":   0,
    "plasmex.electron_kappa":      0,
    "plasmex.stefan_maxwell_drag": 0,
    "prob.wall_temp":              WALL_TEMP,
    "prob.electrode_emission":     int(EMISSION),
    "prob.work_function":          WORK_FN,
    "prob.richardson_A":           RICHARDSON_A,
}


def read_knobs(exe_dir, deck):
    """
    The ENERGY_FLUX_KNOBS values from `deck`, as floats, plus the path read (None when the
    deck was not found, in which case every value is its code default).
    """
    vals, path = deck_values(exe_dir, deck)
    out = {}
    for name, default in ENERGY_FLUX_KNOBS.items():
        try:
            out[name] = float(vals[name])
        except (KeyError, ValueError):
            out[name] = float(default)
    return out, path


def crossing(x, f, level, i_from, i_to):
    """
    First crossing of f(x) == level walking from index i_from towards i_to, linearly
    interpolated. Locating a crossing INSIDE the domain is interpolation between two solved
    cells; it is not the boundary extrapolation this module refuses to do.

    Vectorised: the sign-change test is one numpy pass rather than a Python loop over cells.
    Returns (x_cross, i_lo) or (None, None).
    """
    g = np.asarray(f, dtype=float) - level
    up = i_to > i_from
    lo, hi = (i_from, i_to) if up else (i_to, i_from)
    a, b = g[lo:hi], g[lo + 1:hi + 1]

    # An EXACT zero counts only on the side of the pair the walk arrives from -- the far node
    # of pair k is the near node of pair k+-1 and is tested there. Getting this backwards
    # reports the crossing one cell early on the inward walk from the wall, which is the only
    # direction either caller uses.
    hit = np.nonzero(((a if up else b) == 0.0) | (a * b < 0.0))[0]
    if hit.size == 0:
        return None, None

    k = int(hit[0] if up else hit[-1]) + lo          # walk order: i_from towards i_to
    near = k if up else k + 1
    if g[near] == 0.0:
        return float(x[near]), k
    w = g[k] / (g[k] - g[k + 1])
    return float(x[k] + w * (x[k + 1] - x[k])), k


def at(x, xq, f):
    """Linear interpolation of f onto the single interior point xq (x increasing)."""
    if xq is None:
        return np.nan
    return float(np.interp(xq, x, f))


def n_neutral_lte(t, pres):
    """
    Equilibrium NEUTRAL argon density [m^-3]. Vectorised port of argon_lte::n_neutral in
    exe/chemistry.H -- same coefficients, same closed form, same endpoint clamping on g.

    This is NOT p/(k_B T): that is the total number density, and at 21000 K argon is ~50%
    ionized, so it runs 152x high.
    """
    t = np.asarray(t, dtype=float)
    u = np.clip(t, TFIT_LO, TFIT_HI) * 1.0e-4
    g = np.polyval(G_COEFF, u)
    n_tot = pres / (K_B * t)
    arg = EION_J / (K_B * t)
    # exp(-arg) flushes to zero below ~258 K; there S = 0 and n_Ar = n_tot exactly.
    s_saha = np.where(arg > 700.0, 0.0,
                      g * SAHA_PREFAC * t * np.sqrt(t) * np.exp(-np.minimum(arg, 700.0)))
    q = np.sqrt(s_saha * s_saha + s_saha * n_tot)
    return np.where(arg > 700.0, n_tot,
                    s_saha * n_tot * n_tot / np.maximum((s_saha + q) ** 2, 1e-300))


_NAR_INV_CACHE = {}


def t_gas_lte(n_ar, pres):
    """
    Inverse of n_neutral_lte, by the same construction probinit() uses for ProbParm::nar_inv
    (exe/prob.H): 4096 nodes geometric in T over [250, 60000] K, piecewise linear in
    log(n_Ar), clamped at the ends. Built from n_neutral_lte itself, so this script and the
    solver agree by construction rather than by two independent fits.
    """
    key = float(pres)
    if key not in _NAR_INV_CACHE:
        t_nodes = np.geomspace(250.0, 60000.0, 4096)[::-1]      # n_Ar increasing
        _NAR_INV_CACHE[key] = (np.log(n_neutral_lte(t_nodes, pres)), t_nodes)
    x_nodes, t_nodes = _NAR_INV_CACHE[key]
    return np.interp(np.log(np.maximum(n_ar, 1e-300)), x_nodes, t_nodes)


def gas_temperature(d, pres):
    """
    The gas temperature for a sampled plotfile.

    TGAS is a SOLVED state component since 7 Sep 2026 (src/solve_neutrals.cpp), so read it.
    The LTE inverse below is kept only for plotfiles written before that, where the gas
    temperature really was implicit in AR -- and it is WRONG for a current plotfile, because
    AR is now set by Dalton's law rather than by Saha, so inverting it through the LTE curve
    answers a different question.

    Returns None if neither is available.
    """
    if "TGAS" in d:
        return d["TGAS"]
    if "AR" in d:
        print("NOTE: this plotfile has no TGAS component, so the gas temperature is being "
              "reconstructed by inverting AR through the LTE curve. That is only right for "
              "plotfiles written before the gas temperature became a solved field.")
        return t_gas_lte(d["AR"], pres)
    return None


def ion_transport(t_gas, ndens, emag):
    """
    Ar+ effective temperature, mobility and diffusion coefficient, mirroring the closed-form
    Wannier two-temperature solution in src/aider.H::CalcTransportCoeffs exactly (same
    reference density N0, same quartic solved in closed form for the drift velocity squared).

    T_eff is the reason this is here: it is an intermediate in the solver and never reaches
    the plotfile, so panel 3 has to reconstruct it.

        v_d = mu(T_eff) E ,  mu = mu_0 sqrt(T_gas/T_eff) ,  k T_eff = k T_gas + (M/3) v_d^2
    """
    ndens = np.maximum(np.asarray(ndens, dtype=float), 1e-300)
    mu0_i = 1.85e-4 * N0_REF / ndens
    a_w = (mu0_i * np.asarray(emag, dtype=float)) ** 2
    b_w = ARMASS / (3.0 * K_B * np.asarray(t_gas, dtype=float))
    # vd2 -> a_w as a_w*b_w -> 0, so the E -> 0 limit is exact and finite
    vd2 = (np.sqrt(1.0 + 4.0 * a_w * b_w) - 1.0) / (2.0 * b_w)
    t_eff = t_gas + ARMASS * vd2 / (3.0 * K_B)
    mu_i = mu0_i * np.sqrt(t_gas / t_eff)
    d_i = mu_i * t_eff * (K_B / ECHARGE)                # D = mu k T_eff / e
    return t_eff, mu_i, d_i


_TABLE_CACHE = {}


def load_table(name, exe_dir=None, col=1):
    """
    Load column `col` of one of the whitespace-separated tables as (x, y). ar_e_lorentz.dat
    has six (Te, nubar_ea/n_a, sigma, phi, lambda', lambda), read per column exactly as
    exe/prob.H calls define() once per column on it.

    Same contract as exe/table.H's LinearTableInterpolator with Scale::Linear: strictly
    increasing x, linear in y, and CLAMPED to the endpoints outside the tabulated range --
    which np.interp does by construction, so the two agree without a special case. Returns
    None if the file is not found, so a caller can degrade instead of dying.
    """
    key = (name, col)
    if key in _TABLE_CACHE:
        return _TABLE_CACHE[key]

    dirs = [TABLE_DIR] if TABLE_DIR else search_dirs(exe_dir or os.getcwd())
    tab = None
    for dd in dirs:
        fn = os.path.join(dd, name)
        if os.path.isfile(fn):
            # np.loadtxt skips '#' comment lines by default, which is the same rule
            # exe/table.H's define() applies -- ar_e_lorentz.dat carries a header.
            a = np.loadtxt(fn)
            tab = (np.asarray(a[:, 0], dtype=float), np.asarray(a[:, col], dtype=float))
            break
    _TABLE_CACHE[key] = tab
    return tab


def k_heavy(n_a, n_i, t_h, lnl):
    """
    Heavy-particle thermal conductivity [W/(m K)], a port of plasmex_heavy::HeavyConductivity
    (src/heavy_transport.H): Almeida, Benilov & Naidis (2008) eqs. (54)-(57) with the argon fits
        Qaa = 1.12e-18 T^-0.2,  Qia = 3.6e-18 T^-0.3,  Qii = e^4 lnL/(36 pi eps0^2 (kT)^2),
        k_a = 75 k v n_a/(64 (n_a Qaa + n_i Qia)),  k_i = 75 k v n_i/(64 (n_a Qia + n_i Qii)),
        v = sqrt(pi k T/M).
    Replaced the Ar.dat table (LTE composition) on 20 Sep 2026. If the header changes, change
    this too -- there is no way for the script to read it.
    """
    n_a = np.asarray(n_a, dtype=float); n_i = np.asarray(n_i, dtype=float)
    t_h = np.maximum(np.asarray(t_h, dtype=float), 1.0)
    qaa = 1.12e-18 * t_h**-0.2
    qia = 3.6e-18 * t_h**-0.3
    kt = K_B * t_h
    qii = ECHARGE**4 * lnl / (36.0 * np.pi * EPS0**2 * kt * kt)
    pref = (75.0 / 64.0) * K_B * np.sqrt(np.pi * kt / ARMASS)
    return pref * n_a / (n_a * qaa + n_i * qia) + pref * n_i / (n_a * qia + n_i * qii)


def coulomb_log(ne, te_ev):
    """plasmex_coulomb::CoulombLog: max(2, ln(lam_D/b90)), lam_D on n_e."""
    tes = np.maximum(np.asarray(te_ev, dtype=float), TE_FLOOR_EV)
    nes = np.maximum(np.asarray(ne, dtype=float), 1.0)
    lam_d = np.sqrt(EPS0 * tes / (nes * ECHARGE))
    b90 = ECHARGE / (12.0 * np.pi * EPS0 * tes)
    return np.maximum(2.0, np.log(lam_d / b90))


def k_gas_at(n_a, n_i, ne, te_k, t_face):
    """
    Heavy-gas thermal conductivity [W/(m K)] at the electrode FACE, mirroring exe/bc.H's
    q_wall: HeavyConductivity on the arithmetic face means -- which for T_g under a Dirichlet
    wall is 0.5*(T_gas + T_wall), and for the densities and the lnL inputs is 0.5*(ghost +
    interior) of sborder. The script has no ghosts, so it uses the wall-cell values for those
    (n_Ar is extrapolated there anyway; n_i/n_e/Te differ by the half-cell gradient).
    """
    te_ev = np.asarray(te_k, dtype=float) * K_B / ECHARGE
    return float(k_heavy(n_a, n_i, t_face, coulomb_log(ne, te_ev)))


def coulomb_factor(ni, ne, te_ev, mu_en):
    """
    Electron-ion Coulomb correction to the neutral-only transport, a vectorised port of
    plasmex_coulomb::TransportFactor (src/coulomb.H):

        nu_en = e / (m_e mu_en)
        lnL   = max(2, log(lam_D/b90)),  lam_D = sqrt(eps0 Te/(n_e e)),  b90 = e/(12 pi eps0 Te)
        nu_ei = NUEI_C n_i lnL Te^-3/2              Mitchner & Kruger (5.1d): the ION density
        ->  nu_en / (nu_en + 0.5 nu_ei)

    NOTE the 0.5. That is what src/coulomb.H currently carries; CLAUDE.md still records the
    line as nu_en/(nu_en+nu_ei), and an earlier version carried the Spitzer 0.582. This
    mirrors the TREE. If src/coulomb.H changes, change it here too -- there is no way for
    these scripts to read it.

    mu_en is the MAGNITUDE of the neutral-only mobility sigma/(n_e e). It multiplies
    BOTH mu_e and D_e, which is why the Einstein ratio D/mu is untouched by it.

    All three of the solver's guards are kept, so the factor is exactly 1 wherever the
    correction is switched off -- which is what makes "with" and "without" coincide in the
    cold, unionized part of the domain rather than diverge on a division by zero.
    """
    ni = np.asarray(ni, dtype=float)
    ne = np.asarray(ne, dtype=float)
    te_ev = np.asarray(te_ev, dtype=float)
    mu_en = np.asarray(mu_en, dtype=float)

    ok = (ne > 0.0) & (ni > 0.0) & (te_ev > TE_FLOOR_EV) & (mu_en > 0.0)
    tes = np.where(ok, te_ev, 1.0)
    nes = np.where(ok, ne, 1.0)
    nis = np.where(ok, ni, 1.0)
    mus = np.where(ok, mu_en, 1.0)

    nu_en = ECHARGE / (EMASS * mus)
    lam_d = np.sqrt(EPS0 * tes / (nes * ECHARGE))
    b90 = ECHARGE / (12.0 * np.pi * EPS0 * tes)
    lnl = np.maximum(2.0, np.log(lam_d / b90))
    nu_ei = NUEI_C * nis * lnl / (tes * np.sqrt(tes))

    return np.where(ok, nu_en / (nu_en + 0.5 * nu_ei), 1.0)


def electron_transport(te_k, n_ar, ne, ni, exe_dir=None):
    """
    Electron mobility and diffusion coefficient [m^2/(V s)], [m^2/s], mirroring
    src/aider.H::CalcTransportCoeffs: the reduced tables divided by the neutral density, then
    multiplied by the Coulomb factor.

    mu_e is NEGATIVE -- sigma = n_e e mu_en is tabulated positive and the sign of the charge
    is put back here, as in aider.H -- so v_drift = mu_e E already carries the right direction
    and must not be given another minus.

    Returns a dict, or None when the tables are not on disk. The UNCORRECTED pair is returned
    alongside the corrected one because transport.py plots both, and `clamped` marks where Te
    left the table and np.interp held an endpoint (exe/table.H does the same, silently).
    """
    te_ev = np.asarray(te_k, dtype=float) * K_B / ECHARGE
    nd = np.maximum(np.asarray(n_ar, dtype=float), 1e-300)

    # Mitchner & Kruger's Lorentz-gas conductivity sigma (VIII 2.14a) per unit n_e/n_a, column 2
    # of ar_e_lorentz.dat: mu_en = -sigma/(n_e e), and D = |mu| kT/e with Te clamped to the
    # table's range as the lookup is -- src/aider.H, electron block.
    sig = load_table("ar_e_lorentz.dat", exe_dir, col=2)
    if sig is None:
        return None
    te_cl = np.clip(te_ev, sig[0][0], sig[0][-1])
    mu_neut = -np.interp(te_ev, sig[0], sig[1]) / (nd * ECHARGE)    # negative
    d_neut = -mu_neut * te_cl
    rng = sig
    fac = coulomb_factor(ni, ne, te_ev, np.abs(mu_neut))
    return dict(mu=mu_neut * fac, D=d_neut * fac, mu_neut=mu_neut, D_neut=d_neut,
                fac=fac, te_ev=te_ev,
                clamped=(te_ev < rng[0][0]) | (te_ev > rng[0][-1]),
                te_range=(float(rng[0][0]), float(rng[0][-1])))


# Spitzer-Harm coefficients, Mitchner & Kruger (VIII 5.1a-c) -- src/coulomb.H SIGMA_SH,
# LAMBDA_SH, PHI_SH. Mirrors the tree; if coulomb.H changes, change these too.
SIGMA_SH  = 1.975
LAMBDA_SH = 3.203
PHI_SH    = 1.389


def electron_thermal_coeffs(ni, ne, te_ev, n_ar, d_e, td_on, kappa_model, exe_dir=None):
    """
    Vectorised port of plasmex_coulomb::ElectronThermalCoeffs (src/coulomb.H): returns
    (C_e, A_i, K) -- the thermal-diffusion coefficient C_e = (e/k) phi/sigma, its ion part A_i
    (Almeida A-(49)), and the conductivity factor K = lambda/(n_e k D_e).

        P^-1   = nubar_ei / nubar_ea,   nubar_ea = n_a * col 1 of ar_e_lorentz.dat
        phi/s  = (P^-1 (phi/s)_SH + (phi/s)_ea) / (1 + P^-1)                  A-(48)
        1/lam  = 1/lam_ea + 1/lam_ei,   lam_ea = (n_e/n_a) col 5,
                 lam_ei = LAMBDA_SH n_e k^2 T_e/(m_e nubar_ei)                A-(50)

    td_on False gives C_e = A_i = 0; kappa_model != 1 gives K = 5/2 -- the code's defaults.
    `d_e` is the Coulomb-corrected electron diffusivity (electron_transport()["D"]), which is
    what aider.H hands to the function. Mirrors the TREE; if coulomb.H changes, change it here.
    Returns None when ar_e_lorentz.dat cannot be found.
    """
    tabs = [load_table("ar_e_lorentz.dat", exe_dir, col=c) for c in (1, 2, 3, 5)]
    if any(t is None for t in tabs):
        return None
    te_ev = np.asarray(te_ev, dtype=float)
    ni = np.asarray(ni, dtype=float)
    ne = np.asarray(ne, dtype=float)
    nn = np.maximum(np.asarray(n_ar, dtype=float), 1e-300)
    d_e = np.asarray(d_e, dtype=float)
    nub1, sig1, phi1, lam1 = (np.interp(te_ev, t[0], t[1]) for t in tabs)

    # NuEI with its guards: zero without ions/electrons or below the Te floor.
    ok = (ne > 0.0) & (ni > 0.0) & (te_ev > TE_FLOOR_EV)
    tes = np.where(ok, te_ev, 1.0)
    nu_ei = np.where(ok, NUEI_C * ni * coulomb_log(np.where(ok, ne, 1.0), tes) / (tes * np.sqrt(tes)), 0.0)

    zero = np.zeros_like(te_ev)
    c_e, a_i, kfac = zero.copy(), zero.copy(), np.full_like(te_ev, 2.5)
    if td_on:
        inv_p = nu_ei / (nn * nub1)
        ps_ei = (PHI_SH / SIGMA_SH) * (K_B / ECHARGE)
        ps_e = (inv_p * ps_ei + phi1 / sig1) / (1.0 + inv_p)
        c_e = (ECHARGE / K_B) * ps_e
        a_i = (ECHARGE / K_B) * inv_p * ps_ei / (1.0 + inv_p)
    if kappa_model == 1:
        te_k = np.maximum(te_ev, TE_FLOOR_EV) * ECHARGE / K_B
        inv_lam_ei_ne = EMASS * nu_ei / (LAMBDA_SH * K_B * K_B * te_k)
        kfac = 1.0 / (K_B * d_e * (1.0 / (lam1 / nn) + inv_lam_ei_ne))
    return c_e, a_i, kfac


def emission_flux(e_wall, t_wall, emission, work_fn=WORK_FN, rich_a=RICHARDSON_A):
    """
    Thermionic emission at the electrode, exe/bc.H's species_bc() literal for literal:
    Richardson-Schottky with A = rich_a [A/(m^2 K^2)], work function work_fn [eV] (the deck's
    prob.richardson_A / prob.work_function, defaulting to the code's values), and the barrier
    lowered by the field when it points into the wall (sgn*E > 0).

    Returns (jem, jem_e, jem_app, jem_e_app): the particle flux [m^-2 s^-1, NEGATIVE = into
    the plasma] and its energy flux (2 k T_wall each), always computed, and the pair actually
    APPLIED (zero unless `emission`). bc.H evaluates the lowering on the boundary FACE field,
    which no plotfile carries; callers pass the wall cell's EFLDX and must say so.
    """
    wf = work_fn * ECHARGE
    corrector = np.sqrt(ECHARGE ** 3 * max(SGN * float(e_wall), 0.0) / (4.0 * np.pi * EPS0))
    jem = -rich_a * t_wall * t_wall * np.exp(-(wf - corrector) / (K_B * t_wall)) / ECHARGE
    jem_e = jem * 2.0 * K_B * t_wall
    return jem, jem_e, (jem if emission else 0.0), (jem_e if emission else 0.0)


def electron_energy_flux(d, knobs, exe_dir=None):
    """
    The TOTAL electron energy flux [W/m^2, coordinate-oriented: + = towards x-hi, the
    electrode] -- the flux whose divergence the EEN source terms (EENRXN: Joule, inelastic,
    elastic exchange) balance -- rebuilt on the INTERIOR faces between adjacent leaf cells.
    Only its explicit part ever reaches a plotfile (DIVEENEXPLFLX, and only its divergence);
    the implicit conduction half does not. The continuum form of what src/aider.H assembles
    as mob(EEN), diff(EEN) and xvel(EEN) (EEN note in CalcTransportCoeffs, CalcExtraVelocities):

        Q = (5/2 + C_e) k T_e Gamma_e                         "enthalpy" (+ Onsager)
          - K n_e k D_e grad T_e                              "conduction"
          - k T_e n_e [A_i (u_i - u) + f_drag A_a (u_n - u)]  "heavy" (A-(15)), A_a = C_e - A_i

    Gamma_e = -JEX/e is the solver's own electron flux (drift, diffusion, thermal, drag and
    frame terms included), so it is NOT rebuilt; u_i = JIX/(e ARP), u_n = UNX, u = HVELX, each
    0 when absent. The face value of every input is the arithmetic mean of the two cells and
    grad T_e is the two-cell difference over the centre spacing -- non-uniform across a
    coarse-fine interface, where it is not the operator's closure gradient, so expect a
    blip there. Everything is time-n plotfile state, whereas the solver lags u_i and u_n one
    step and uses T_e^n in its coefficients: agreement to ~1-2% is what to expect.

    And the flux the SOURCES imply, anchored on the wall BC and integrated inward:

        Q_src(x_f) = Q_wall - sum_{leaf cells right of x_f} EENRXN * dx,
        Q_wall     = 2 k T_e n_e cbar/4 + jem_e_app        (exe/bc.H: eeng_maxw + Jem_E)

    This is the STEADY-STATE flux -- d(EEN)/dt is dropped -- so Q and Q_src coincide exactly
    when the run is steady and the budget closes, and their gap is otherwise the transient.
    Jem uses the wall CELL field (see emission_flux()).

    `knobs` is read_knobs()'s dict. Returns a dict of numpy arrays: xf (interior faces), TOT,
    ENTH, COND, HEAVY and HEAT = COND + HEAVY (the heat flux proper) on xf; xf_src (interior faces + the wall face), SRC on xf_src; plus
    q_wall and jem_app. None when the transport tables cannot be found or JEX/EENRXN are absent.
    """
    if "JEX" not in d or "EENRXN" not in d:
        return None
    x = d["x"]
    ne, ni, na, te = d["NE"], d["ARP"], d["AR"], d["ETEMP"]
    te_ev = te * K_B / ECHARGE

    et = electron_transport(te, na, ne, ni, exe_dir)
    if et is None:
        return None
    th = electron_thermal_coeffs(ni, ne, te_ev, na, et["D"],
                                 bool(knobs["plasmex.thermal_diffusion"]),
                                 int(knobs["plasmex.electron_kappa"]), exe_dir)
    if th is None:
        return None
    c_e, a_i, kfac = th
    f_drag = 1.0 if knobs["plasmex.stefan_maxwell_drag"] else 0.0

    zeros = np.zeros_like(x)
    u_i = d["JIX"] / (ECHARGE * np.maximum(ni, 1e-300)) if "JIX" in d else zeros
    u_n = d.get("UNX", zeros)
    u_h = d.get("HVELX", zeros)

    def face(f):
        f = np.asarray(f, dtype=float)
        return 0.5 * (f[1:] + f[:-1])

    xf = face(x)
    dxf = x[1:] - x[:-1]
    gam = face(-d["JEX"] / ECHARGE)
    te_f, ne_f, c_f = face(te), face(ne), face(c_e)
    enth = (2.5 + c_f) * K_B * te_f * gam
    cond = -face(kfac * ne * K_B * et["D"]) * (te[1:] - te[:-1]) / dxf
    heavy = -K_B * te_f * ne_f * (face(a_i) * face(u_i - u_h)
                                  + f_drag * face(c_e - a_i) * face(u_n - u_h))

    # ---- the source-integrated flux, from the wall face inward ----
    iw = len(x) - 1
    t_wall = knobs["prob.wall_temp"]
    _, _, jem_app, jem_e_app = emission_flux(d["EFLDX"][iw] if "EFLDX" in d else 0.0,
                                             t_wall, bool(knobs["prob.electrode_emission"]),
                                             knobs["prob.work_function"], knobs["prob.richardson_A"])
    cbar = np.sqrt(8.0 * K_B * te[iw] / (np.pi * EMASS))
    q_wall = SGN * (2.0 * K_B * te[iw] * ne[iw] * cbar / 4.0 + jem_e_app)
    s_w = np.asarray(d["EENRXN"], dtype=float) * d["dxc"]
    # tail[i] = sum over cells i..end; the face between cells i-1 and i sees tail[i].
    tail = np.cumsum(s_w[::-1])[::-1]
    src = np.r_[q_wall - tail[1:], q_wall]
    xf_src = np.r_[xf, d["hi"]]

    return dict(xf=xf, TOT=enth + cond + heavy, ENTH=enth, COND=cond, HEAVY=heavy,
                HEAT=cond + heavy,
                xf_src=xf_src, SRC=src, q_wall=q_wall, jem_app=jem_app)


def cell_gradient(d, f, i):
    """
    d(f)/dx at leaf cell i, from that cell and its two leaf neighbours on the non-uniform
    3-point stencil (exact for a quadratic, and the ordinary central difference when the
    three cells are the same size). One-sided at the ends of the domain.

    Because sample() returns LEAF CELLS, cell i-1 and i+1 are genuinely the adjacent cells at
    whatever level they sit on, for any refinement ratio. The old fine-resampled layout needed
    a stride table to find them and got it wrong when ref_ratio was not 2 -- returning a
    bit-identically zero gradient, which for electrons flips the sign of the total flux.

    Returns (dfdx, x_centre).
    """
    x = d["x"]
    f = np.asarray(f, dtype=float)
    n = len(x)

    if n < 2:
        return np.nan, float(x[i])
    if i == 0:
        return float((f[1] - f[0]) / (x[1] - x[0])), float(x[0])
    if i == n - 1:
        return float((f[i] - f[i - 1]) / (x[i] - x[i - 1])), float(x[i])

    hl, hr = x[i] - x[i - 1], x[i + 1] - x[i]
    g = (-hr / (hl * (hl + hr)) * f[i - 1]
         + (hr - hl) / (hl * hr) * f[i]
         + hl / (hr * (hl + hr)) * f[i + 1])
    return float(g), float(x[i])


def drift_diffusion(d, i):
    """
    Split the flux of each charged species at leaf cell i into its drift and diffusion halves,
    sgn-oriented so positive means "towards the electrode":

        Gamma_drift = mu n E        Gamma_diff = -D grad(n)        Gamma = the sum

    and the velocities that follow, v = Gamma/n. Both are reported; they differ by a factor n
    and are not interchangeable.

    Everything is rebuilt at the CELL CENTRE from cell-centred plotfile fields, so this is not
    bit-for-bit what the solver formed on the faces (mu and D there come from a two-cell
    average, and the drift flux is upwinded rather than centred). The report prints the
    solver's own total beside it as the check -- ION_VEL for the ions, JEX for the electrons.

    IN THE SHEATH THE ELECTRON RECONSTRUCTION IS NOT TRUSTWORTHY, and the failure is
    intrinsic, not a bug: drift and diffusion each hugely exceed their sum (Boltzmann
    balance), so the net is a small difference of two large near-cancelling numbers and a
    percent-level error in grad(n_e) becomes an order-unity error in the total. The IONS are
    unaffected because their drift dominates. Read the electron total in the bulk; in the
    sheath read JEX itself.

    THE ION CHECK CARRIES A ONE-STEP LAG. mu_i is field-dependent (Wannier) and src/aider.H
    builds it from the PREVIOUS step's cell-centred E, while this rebuilds it from the E in
    the plotfile. On a transient that alone is worth ~10%; the two coincide at steady state.

    mu_e is negative and mu_i positive, so v_drift already carries the species' sign; only
    grad(n) needs the explicit minus.
    """
    ni = float(d["ARP"][i])
    ne = float(d["NE"][i])
    e_x = float(d["EFLDX"][i])

    dni, x0 = cell_gradient(d, d["ARP"], i)
    dne, _ = cell_gradient(d, d["NE"], i)

    _, mu_i, d_i = ion_transport(float(d["T_GAS"][i]), float(d["AR"][i]), abs(e_x))

    out = {"x": x0, "e_field": SGN * e_x, "ne": ne, "ni": ni,
           "grad_ni": dni, "grad_ne": dne,
           "boundary_cell": (i == 0) or (i == len(d["x"]) - 1)}

    out["ion"] = _channel(float(mu_i), float(d_i), e_x, ni, dni)
    # The solver's own cell-centred ion fluid velocity, drift plus diffusion (ION_VEL, written
    # by DiagnoseIonCurrent from the two axial faces). ION_VEL is NOT face-overwritten, so in a
    # boundary cell it is a two-face average that includes the BC-set domain face -- a 10-30%
    # gap there is the boundary, not an error in the reconstruction.
    out["ion"]["solver"] = SGN * float(d["ION_VEL"][i]) if "ION_VEL" in d else np.nan
    out["ion"]["solver_src"] = f"ION_VEL {face_tag('ION_VEL')}"

    et = electron_transport(float(d["ETEMP"][i]), float(d["AR"][i]), ne, float(d["ARP"][i]), d["dir"])
    if et is None:
        out["electron"] = None
        out["e_off"] = ("ar_e_lorentz.dat not found -- set TABLE_DIR to the directory "
                        "the run used")
        return out

    out["electron"] = _channel(float(et["mu"]), float(et["D"]), e_x, ne, dne)
    out["e_off"] = None

    # Gamma_e = JEX / (-e). An EXACT zero is not a measurement: JEX is an ordinary state
    # component, so it is always in the header and carries whatever the last diagnostic step
    # left there (a t = 0 plotfile carries zero). MEASURED, it is also exactly zero in the ARC
    # boundary cell of every plotfile checked -- BoundaryFaceOverwrite puts the x-lo domain
    # FACE value there and the flux arrays are not written on that face -- so at the default
    # ARC_REF_FRAC = 0 the check is simply unavailable. Say so rather than reporting zero.
    jex = float(d["JEX"][i]) if "JEX" in d else None
    if jex is None or jex == 0.0 or ne <= 0.0:
        out["electron"]["solver"] = np.nan
        out["electron"]["solver_src"] = (
            "no JEX in this plotfile" if jex is None else
            "JEX is exactly zero in this cell (the x-lo domain face carries no current in "
            "the plotfile; also true of a t = 0 plotfile)")
    else:
        out["electron"]["solver"] = SGN * jex / (-ECHARGE) / ne
        out["electron"]["solver_src"] = f"JEX/(-e n_e) {face_tag('JEX')}"

    return out


def _channel(mu, dcoef, e_x, n, dndx):
    """
    One species' drift/diffusion decomposition, as FLUXES and as the velocities that follow:

        Gamma_drift = mu n E        Gamma_diff = -D grad(n)        Gamma = the sum

    The velocity is that divided by the local density, which is the ordinary fluid-velocity
    definition and the form ION_VEL is written in. So v_diff = -D grad(n)/n, NOT -D grad(n):
    the two differ by a factor n and are not the same units.
    """
    g_dr = SGN * mu * n * e_x
    g_df = SGN * (-dcoef * dndx)
    inv_n = (1.0 / n) if n > 0.0 else np.nan
    return {"mu": mu, "D": dcoef, "n": n,
            "g_drift": g_dr, "g_diff": g_df, "g_total": g_dr + g_df,
            "v_drift": g_dr * inv_n, "v_diff": g_df * inv_n,
            "v_total": (g_dr + g_df) * inv_n}


# ============================================================
# SAMPLE ONE PLOTFILE ONTO ITS AMR LEAF CELLS
# ============================================================

def sample(fn, fields=None, want_derived=True):
    """
    Return the AMR LEAF cells of a 1-D plotfile in x order -- the cells the solver actually
    solved on, each at its own width, with no resampling and no interpolation.

    A cell at level L is a leaf when no level-(L+1) grid covers it. That test is done in x,
    against the finer grids' own edges, so amr.ref_ratio never appears: any ratio works, and
    a per-level list of different ratios works. Contrast the two things this replaces:

      * yt.LineBuffer (the old gif.py) does not sample at cell centres at all -- measured
        against a cell-centred linspace on the production deck it drifts from 0.5*dx of
        offset at the arc to 1.5*dx at the electrode, a cell and a half of misregistration
        exactly where the sheath is -- and it re-traces the ray once per field.
      * a uniform fine resampling (the old bohm.py) replicates each coarse cell across the
        fine samples it covers, so consecutive samples are identical, a naive gradient over
        them is zero, and finding the true neighbour needs a per-level stride that is wrong
        unless ref_ratio is 2.

    `fields` is the list of plotfile components to read; None reads DEFAULT_FIELDS. Whatever
    is asked for but absent from the header is skipped silently -- older plotfiles predate
    JEX/JIX/EEN and the callers degrade rather than raise -- and the names that were missing
    come back in d["missing"] so a report can say so.

    Keys: x, dxc (per-cell width), dx (finest width), level, lo, hi, ncell, finest, ds, file,
    dir, plus one array per field and, with want_derived, T_GAS and T_EFF.
    """
    ds = yt.load(fn)

    if int(ds.dimensionality) != 1:
        raise RuntimeError(f"{fn} is {ds.dimensionality}-D; these scripts are 1-D only "
                           f"(use gif_2d.py for r-z).")

    lo = float(ds.domain_left_edge.d[0])
    hi = float(ds.domain_right_edge.d[0])

    have = {f[1] for f in ds.field_list if f[0] == "boxlib"}
    want = DEFAULT_FIELDS if fields is None else fields
    names = [n for n in want if n in have]
    missing = [n for n in want if n not in have]

    by_level = {}
    for g in ds.index.grids:
        by_level.setdefault(int(g.Level), []).append(g)
    finest = max(by_level)

    xs, ws, lv = [], [], []
    cols = {n: [] for n in names}

    for level in sorted(by_level):
        # x-extent of every grid one level finer; a cell whose centre falls inside one of
        # them is covered and is not a leaf.
        child = [(float(g.LeftEdge.d[0]), float(g.RightEdge.d[0]))
                 for g in by_level.get(level + 1, ())]

        for g in by_level[level]:
            gl = float(g.LeftEdge.d[0])
            w = float(np.asarray(g.dds, dtype=float)[0])
            n = int(np.asarray(g.ActiveDimensions, dtype=int)[0])

            xc = gl + (np.arange(n) + 0.5) * w
            keep = np.ones(n, dtype=bool)
            for a, b in child:
                keep &= ~((xc > a) & (xc < b))
            if not keep.any():
                continue

            xs.append(xc[keep])
            ws.append(np.full(int(keep.sum()), w))
            lv.append(np.full(int(keep.sum()), level, dtype=int))
            # yt pads a 1-D AMReX dataset out to 3-D with a single transverse cell, so the
            # axial line is the whole block; ravel is the axial run.
            for name in names:
                cols[name].append(np.asarray(g[("boxlib", name)]).ravel()[keep])

    x = np.concatenate(xs)
    order = np.argsort(x, kind="stable")

    d = {name: np.concatenate(cols[name])[order] for name in names}
    d["x"] = x[order]
    d["dxc"] = np.concatenate(ws)[order]
    d["level"] = np.concatenate(lv)[order]

    # The leaves must tile the domain exactly: AMReX grids are properly nested and a fine
    # patch replaces a whole number of coarse cells. A gap or an overlap here means the leaf
    # test above missed something, and every edge value below would then be quoting the wrong
    # cell -- so it is checked rather than assumed.
    span = float(d["dxc"].sum())
    if abs(span - (hi - lo)) > 1e-9 * (hi - lo):
        raise RuntimeError(f"{fn}: leaf cells span {span:.6e} m but the domain is "
                           f"{hi - lo:.6e} m -- the AMR leaf walk is inconsistent")

    d.update({"dx": float(d["dxc"].min()), "lo": lo, "hi": hi, "ncell": len(d["x"]),
              "finest": finest, "ds": ds, "file": fn,
              "dir": os.path.dirname(os.path.abspath(fn)), "missing": missing})

    if want_derived:
        tg = gas_temperature(d, d.get("pres", PRESSURE))
        if tg is not None:
            d["T_GAS"] = tg
            if "EFLDX" in d and "AR" in d:
                d["T_EFF"] = ion_transport(d["T_GAS"], d["AR"], np.abs(d["EFLDX"]))[0]

    return d


def set_window(d):
    """
    Apply PROB_LO / PROB_HI, the plot window in metres. None takes the domain edge out of the
    plotfile, so the default is always the whole domain of whatever was loaded.
    """
    xlo = d["lo"] if PROB_LO is None else float(PROB_LO)
    xhi = d["hi"] if PROB_HI is None else float(PROB_HI)
    if xhi <= xlo:
        raise SystemExit(f"PROB_HI ({xhi}) must be greater than PROB_LO ({xlo})")
    d["xlim"] = (xlo, xhi)
    d["mask"] = (d["x"] >= xlo) & (d["x"] <= xhi)
    return d


# ============================================================
# ANALYSIS
# ============================================================

def analyse(d):
    """
    Locate the sheath edge and the Bohm point, and build the Boltzmann curve.

    Everything is oriented with SGN = +1 (electrode at x-hi), so a positive velocity means
    "towards the electrode".
    """
    x = d["x"]
    ne, ni, phi = d["NE"], d["ARP"], d["PHI"]

    i_wall = len(x) - 1
    i_bulk = 0

    te_v = d["ETEMP"] * K_B / ECHARGE                # electron temperature in volts

    # ION_VEL is the full ion fluid velocity, drift plus diffusion, straight off the solver's
    # own flux arrays. The Bohm criterion is about the total velocity, so it is not decomposed.
    u_flux = SGN * d["ION_VEL"] if "ION_VEL" in d else np.full_like(x, np.nan)
    cs_cold = np.sqrt(ECHARGE * te_v / ARMASS)       # the cold-ion Bohm speed

    # ---- the sheath edge: quasi-neutrality breaks, walking in from the wall ----
    sep = np.abs(ni - ne) / np.maximum(ni, 1e-300)
    x_edge, _ = crossing(x, sep, QN_EPS, i_wall, i_bulk)

    # ---- the Bohm point: where the ion fluid velocity first reaches c_s ----
    x_bohm, _ = crossing(x, u_flux - cs_cold, 0.0, i_wall, i_bulk)

    # ---- the arc (bulk) reference ----
    i_ref = i_bulk if ARC_REF_FRAC <= 0.0 else int(
        min(max(ARC_REF_FRAC * len(x), 0), len(x) - 1))
    arc = dict(i=i_ref, x=float(x[i_ref]), ne=float(ne[i_ref]), ni=float(ni[i_ref]),
               n_ar=float(d["AR"][i_ref]), te_v=float(te_v[i_ref]), phi=float(phi[i_ref]),
               t_gas=float(d["T_GAS"][i_ref]), e_field=float(SGN * d["EFLDX"][i_ref]),
               vel=drift_diffusion(d, i_ref))

    # ---- the Boltzmann prediction, anchored at the sheath edge ----
    # n_e should follow n_e,edge exp((Phi - Phi_edge)/Te_edge) across a retarding sheath.
    # Te_edge, not the local Te, is the right normaliser: a reflective retarding barrier is
    # energy conserving, so by Liouville it maps a Maxwellian onto a Maxwellian at the SAME
    # Te, and Te(x) = const is the collisionless prediction.
    edge = None
    ne_boltz = np.full_like(ne, np.nan)
    if x_edge is not None:
        edge = dict(x=float(x_edge), ne=at(x, x_edge, ne), ni=at(x, x_edge, ni),
                    te_v=at(x, x_edge, te_v), phi=at(x, x_edge, phi),
                    cs=at(x, x_edge, cs_cold), u=at(x, x_edge, u_flux))
        ne_boltz = edge["ne"] * np.exp((phi - edge["phi"]) / edge["te_v"])

    # Panel 4 and the min/max table read these off d like any plotfile field. u_flux is
    # sgn-oriented, so positive means "towards the electrode" -- the direction the Bohm
    # criterion is about, and the reason the sign is kept rather than taking |u_i|.
    d["NE_BOLTZ"], d["U_ION"], d["CS"] = ne_boltz, u_flux, cs_cold

    return dict(x=x, i_wall=i_wall, i_bulk=i_bulk, te_v=te_v, u_flux=u_flux, cs_cold=cs_cold,
                x_edge=x_edge, x_bohm=x_bohm, edge=edge, arc=arc,
                sheath_fall=np.nan if edge is None else edge["phi"] - float(phi[i_wall]))


# ============================================================
# ELECTRODE HEAT FLUX
# ============================================================

def electrode_heat_flux(d, a):
    """
    Power delivered to the electrode, channel by channel [W/m^2], MEASURED off the plotfile's
    WALL CELL. Positive = into the wall.

        q = <electron kinetic>               energy carried by the electron flux
          + <electron flux> * work_fn        condensation of arriving electrons
          + i_wall * (E_ion - work_fn)       ion neutralisation at the surface
          + 2 k(T_face) (T_gas - T_wall)/dx  heavy-gas conduction across the half cell,
                                             k_hp (Almeida 54-57) at T_face = (T_gas+T_wall)/2

    Both charged fluxes are read out of the wall cell, not rebuilt from exe/bc.H:

        Gamma_e = sgn * JEX / (-e)      [FACE]
        Gamma_i = sgn * JIX / ( e)      [FACE]

    JEX and JIX are face-overwritten there, so those ARE the fluxes the solver applied across
    the wall. That matters for the ions in particular: bc.H's own i_wall expression is DEAD
    whenever the deck declares the heavy group Dirichlet at the electrode (heavy_bc_hi = 1,
    which exe/inputs ships), so reconstructing it would report a wall the run did not have.
    And at implicit_level = 3 the electron flux crossing the face is not the n_e c_bar/4 the
    BC asked for either -- ExplicitDriftFlx is re-run after the species solve.

    THE SCHOTTKY TERM IS THE ONE PLACE THIS CANNOT MATCH bc.H. bc.H lowers the barrier with
    the boundary FACE field efld_face, but EFLDX is NOT in FACE_CELL_FIELDS -- the wall cell
    holds the ordinary face-to-cell average -- and this script does not extrapolate. So Jem
    here is the CELL-VALUE estimate, reported as such. E rises steeply into the wall and the
    lowering is exp(sqrt(E)), so expect a real difference, not a rounding one.

    The energy SPLIT is a model in either case: there is no JEEN diagnostic, so only the
    number flux is measured. Arriving electrons are given 2 k Te each and emitted ones
    2 k T_wall each, which is the Maxwellian half-flux identity
    eeng_maxw/e_maxw = (4/3)(EEN/n_e) = 2 k Te, exact in bc.H by construction.
    """
    iw = a["i_wall"]
    t_wall, emission = d["wall_temp"], d["emission"]
    wf = d["work_fn"] * ECHARGE

    te_k = float(d["ETEMP"][iw])
    t_gas = float(d["T_GAS"][iw])

    # sgn*E > 0 means the field points out of the domain into the wall, the orientation that
    # accelerates an emitted electron back into the plasma; a repelling field gives zero
    # lowering rather than a negative sqrt argument.
    e_cell = float(d["EFLDX"][iw])

    # jem/jem_e are always COMPUTED -- they are reported either way -- but only APPLIED when
    # the deck asks for emission. Every consumer must use the applied pair, or emission = 0
    # silently does nothing. exe/bc.H carries the same split under the same names.
    jem, jem_e, jem_app, jem_e_app = emission_flux(e_cell, t_wall, emission,
                                                   d["work_fn"], d["rich_a"])

    # The wall cell's OWN width, straight off its grid, so this is the same dx[dir] bc.H
    # divides by at any refinement ratio. The 2 is the half-cell distance to the wall.
    # Conductivity at the WALL FACE, k_hp(n_Ar, n_i, 0.5*(T_gas + T_wall)) -- the same face
    # temperature exe/bc.H uses, not k(T_gas) of the cell centre.
    k_face = k_gas_at(float(d["AR"][iw]), float(d["ARP"][iw]), float(d["NE"][iw]),
                      float(d["ETEMP"][iw]), 0.5 * (t_gas + t_wall))
    q_cond = 2.0 * k_face * (t_gas - t_wall) / float(d["dxc"][iw])

    # ---- ion channel ----
    if "JIX" in d:
        i_wall_flux = SGN * float(d["JIX"][iw]) / ECHARGE
        i_src = f"JIX {face_tag('JIX')}"
    elif "ION_VEL" in d:
        # ION_VEL is a two-face average, not the boundary face, so this is only an estimate.
        i_wall_flux = SGN * float(d["ARP"][iw]) * float(d["ION_VEL"][iw])
        i_src = f"ARP*ION_VEL {face_tag('ION_VEL')} (no JIX in this plotfile)"
    else:
        i_wall_flux, i_src = np.nan, "unavailable (no JIX, no ION_VEL)"
    q_ion = i_wall_flux * (E_ION * ECHARGE - wf)

    # ---- electron channel ----
    # An EXACT zero means the split currents were never written for this plotfile, not that no
    # electrons reached the wall. JEX/JIX are ordinary state components, so they are always in
    # the header -- they just hold whatever the last diagnostic step left there, and t = 0
    # leaves zero. Guard that rather than reporting it as a measurement.
    e_solved = SGN * float(d["JEX"][iw]) / (-ECHARGE) if "JEX" in d else None
    if e_solved is None:
        off = "no JEX in this plotfile"
    elif e_solved == 0.0:
        off = ("JEX is exactly zero in the wall cell -- the split currents were never "
               "written for this plotfile (t = 0?)")
    else:
        off = None

    terms = total = None
    if off is None:
        # Split the net flux back into arriving and emitted: they carry different energies per
        # particle (2 k Te vs 2 k T_wall), so the net alone cannot be turned into an energy
        # flux. jem is negative and is already zero when emission is off.
        gamma_arr = e_solved - jem_app
        terms = [
            ("electron kinetic  (2 k Te * flux)", 2.0 * K_B * te_k * gamma_arr + jem_e_app),
            ("electron work fn  (flux * phi_w)", e_solved * wf),
            ("ion neutralisation", q_ion),
            ("gas conduction", q_cond),
        ]
        total = sum(v for _, v in terms)

    return dict(terms=terms, total=total, off=off, jem=jem, jem_app=jem_app,
                e_solved=e_solved, i_wall_flux=i_wall_flux, i_src=i_src, t_gas=t_gas,
                te_k=te_k, ne=float(d["NE"][iw]), ni=float(d["ARP"][iw]),
                e_cell=e_cell, dx_wall=float(d["dxc"][iw]),
                k_face=k_face)


# ============================================================
# REPORT
# ============================================================

# (key, label) -- the min/max table and the edge table walk this list.
PLOTTED = [
    ("NE",       "n_e            [m^-3]"),
    ("ARP",      "n_i            [m^-3]"),
    ("NE_BOLTZ", "n_e Boltzmann  [m^-3]"),
    ("PHI",      "Phi            [V]"),
    ("EFLDX",    "E_x            [V/m]"),
    ("ETEMP",    "Te             [K]"),
    ("T_GAS",    "T_gas          [K]"),
    ("T_EFF",    "T_eff          [K]"),
    ("U_ION",    "u_i            [m/s]"),
    ("CS",       "c_s (Bohm)     [m/s]"),
    ("JEX",      "J_e            [A/m^2]"),
    ("JIX",      "J_i            [A/m^2]"),
    ("JX",       "J_total        [A/m^2]"),
]


def edge_table(d):
    """
    The boundary-CELL values at both domain ends, tagged with what each one actually is.

    This is the block the whole no-extrapolation policy exists for. Every number is the value
    stored in the first / last leaf cell, read verbatim. `[FACE]` marks the fields whose
    boundary cell was overwritten with the domain face value by BoundaryFaceOverwrite, i.e.
    the ones where the cell value IS the face value; `[cell]` marks an ordinary cell-centred
    value, which for a face quantity is an average over the cell and NOT what crossed the
    boundary.
    """
    ia, iw = 0, len(d["x"]) - 1
    print("  EDGE VALUES -- boundary CELL values, read verbatim, never extrapolated")
    print(f"    arc  cell 0      x = {d['x'][ia]:.6e} m   width {d['dxc'][ia]:.4e} m   "
          f"level {d['level'][ia]}")
    print(f"    wall cell {iw:<6d} x = {d['x'][iw]:.6e} m   width {d['dxc'][iw]:.4e} m   "
          f"level {d['level'][iw]}")
    print(f"    {'':24s} {'':6s} {'arc (x-lo)':>16s} {'electrode (x-hi)':>18s}")
    for key, label in PLOTTED:
        if key not in d:
            continue
        v = np.asarray(d[key], dtype=float)
        print(f"    {label:24s} {face_tag(key):6s} {v[ia]:16.6e} {v[iw]:18.6e}")
    print("    [FACE] = boundary cell overwritten with the domain FACE value "
          "(BoundaryFaceOverwrite, src/aider.H)")
    print("    [cell] = ordinary cell-centred value; NOT the face value, and not "
          "extrapolated to one")
    print()


def report(d, a, q):
    x, iw = a["x"], a["i_wall"]
    xw = x[iw]
    xlo, xhi = d["xlim"]
    edge = a["edge"]

    print()
    print("=" * 78)
    print(f"  {d['file']}   t = {float(d['ds'].current_time):.6e} s")
    print(f"  1-D, finest level {d['finest']}, {d['ncell']} leaf cells, "
          f"dx {d['dxc'].max():.3e} .. {d['dx']:.3e} m")
    print(f"  deck: {d['deck_src']}   pres = {d['pres']:.6g} Pa   "
          f"wall_temp = {d['wall_temp']:.6g} K   emission = {int(d['emission'])}")
    print("=" * 78)
    print(f"  electrode at x-hi, wall cell centre x = {xw:.6e} m")
    print(f"  plot window  x = [{xlo:.6e}, {xhi:.6e}] m"
          f"{'   (log x, distance from the electrode)' if XLOG else ''}")
    if d["missing"]:
        print(f"  NOT IN THIS PLOTFILE: {', '.join(d['missing'])}")
    print()

    edge_table(d)

    # ---------------------------------------------------------------- sheath
    print("  SHEATH")
    if edge is None:
        print(f"    sheath edge  |n_i-n_e|/n_i = {QN_EPS:g}      NOT FOUND")
    else:
        print(f"    sheath edge  |n_i-n_e|/n_i = {QN_EPS:g}      "
              f"x = {edge['x']:.6e} m   width = {abs(xw - edge['x']):.3e} m")
    if a["x_bohm"] is None:
        print("    Bohm point   u_i = c_s                    NOT FOUND")
    else:
        print(f"    Bohm point   u_i = c_s                    x = {a['x_bohm']:.6e} m   "
              f"{abs(xw - a['x_bohm']):.3e} m from the electrode")

    if edge is not None:
        print(f"    SHEATH FALL      Phi(edge) - Phi(wall cell) = {a['sheath_fall']:.6f} V"
              f"   = {a['sheath_fall'] / edge['te_v']:.3f} Te,edge")
        print(f"    BOHM VELOCITY    c_s at the sheath edge = {edge['cs']:.6e} m/s"
              f"   (Te,edge = {edge['te_v']:.4f} V)")
        print(f"                     u_i there              = {edge['u']:.6e} m/s"
              f"   -> u_i/c_s = {edge['u'] / edge['cs']:.4f}")
    print()

    # ---------------------------------------------------------------- arc
    arc = a["arc"]
    v = arc["vel"]
    print(f"  ARC EDGE  (cell {arc['i']}, x = {arc['x']:.6e} m, "
          f"{abs(xw - arc['x']):.3e} m from the electrode"
          f"{'; this IS the boundary cell' if v['boundary_cell'] else ''})")
    print(f"    E   = {arc['e_field']:.6e} V/m   {face_tag('EFLDX')}   "
          f"(positive = towards the electrode)")
    print(f"    n_e = {arc['ne']:.4e} m^-3   Te = {arc['te_v']:.4f} V "
          f"({arc['te_v'] * ECHARGE / K_B:.1f} K)   T_gas = {arc['t_gas']:.1f} K")
    print()

    print("  DRIFT / DIFFUSION SPLIT THERE   (positive = towards the electrode)")
    print(f"    grad(n_i) = {v['grad_ni']:.6e}   grad(n_e) = {v['grad_ne']:.6e} m^-4")
    print("    (mu and D are rebuilt as src/aider.H does; neither reaches the plotfile. The")
    print("     solver's own total is the check -- see drift_diffusion() for the two reasons")
    print("     it legitimately differs: the ion mu carries a one-step lag in E, and ION_VEL")
    print("     in a boundary cell averages in the BC-set domain face.)")
    for tag, label in (("ion", "Ar+"), ("electron", "electron")):
        c = v[tag]
        if c is None:
            print(f"    {label:8s}  NOT AVAILABLE: {v['e_off']}")
            continue
        print(f"    {label:8s}  mu = {c['mu']:.6e} m^2/(V s)   D = {c['D']:.6e} m^2/s")
        print(f"      {'':26s} {'flux [m^-2 s^-1]':>16s}  {'vel = flux/n [m/s]':>20s}")
        print(f"      {'drift      mu n E':26s} {c['g_drift']:16.6e}  {c['v_drift']:20.6e}")
        print(f"      {'diffusion  -D grad(n)':26s} {c['g_diff']:16.6e}  {c['v_diff']:20.6e}")
        print(f"      {'TOTAL':26s} {c['g_total']:16.6e}  {c['v_total']:20.6e}")
        sv = c["solver"]
        if np.isfinite(sv):
            rel = (sv - c["v_total"]) / sv if sv != 0.0 else np.nan
            print(f"      {'solver ' + c['solver_src']:26s} {sv * c['n']:16.6e}"
                  f"  {sv:20.6e}   ({100.0 * rel:+.2f} %)")
        else:
            print(f"      solver check NOT AVAILABLE ({c['solver_src']})")
    print()

    # ---------------------------------------------------------------- min / max
    print("  MIN / MAX OF THE PLOTTED FIELDS, over the plot window")
    m = d["mask"]
    xs_all = x[m]
    for key, label in PLOTTED:
        if key not in d:
            continue
        vv = np.asarray(d[key], dtype=float)[m]
        good = np.isfinite(vv)
        if not good.any():
            print(f"    {label:24s}  (no finite samples in the window)")
            continue
        vv, xs = vv[good], xs_all[good]
        imn, imx = int(np.argmin(vv)), int(np.argmax(vv))
        print(f"    {label:24s}  min {vv[imn]:14.6e} at x = {xs[imn]:.6e}"
              f"      max {vv[imx]:14.6e} at x = {xs[imx]:.6e}")
    print()

    # ---------------------------------------------------------------- heat flux
    print(f"  HEAT FLUX TO THE ELECTRODE   [W/m^2]   "
          f"(emission {'ON' if d['emission'] else 'OFF'})")
    print(f"    wall cell: n_e = {q['ne']:.4e}   n_i = {q['ni']:.4e} m^-3   "
          f"Te = {q['te_k']:.1f} K   T_gas = {q['t_gas']:.1f} K   dx = {q['dx_wall']:.4e} m")
    print(f"    Jem = {q['jem']:.6e} m^-2 s^-1"
          f"{'' if d['emission'] else '  (NOT applied)'}")
    print(f"      from E = {q['e_cell']:.4e} V/m {face_tag('EFLDX')} -- bc.H uses the "
          f"boundary FACE field here,")
    print( "      which the plotfile does not carry, so this is the cell-value estimate; "
           "see electrode_heat_flux().")
    print(f"    ion flux from {q['i_src']}: {q['i_wall_flux']:.6e} m^-2 s^-1")
    if q["terms"] is None:
        print(f"    NOT AVAILABLE: {q['off']}")
    else:
        print(f"    electron flux from JEX {face_tag('JEX')}: "
              f"{q['e_solved']:.6e} m^-2 s^-1")
        for label, val in q["terms"]:
            share = 100.0 * val / q["total"] if q["total"] != 0.0 else np.nan
            print(f"      {label:36s} {val:14.6e}   {share:6.1f} %")
        print(f"      {'TOTAL':36s} {q['total']:14.6e}")
    print()


# ============================================================
# PLOT HELPERS  (shared with transport.py and gif.py)
# ============================================================

def autoscale(ax, mask, series, logy, symlog=False):
    """
    Set the y limits from the data actually inside the plotted x window, so a zoomed
    PROB_LO/PROB_HI does not leave the panel on the full-domain range.

    `symlog` is for a panel whose data changes sign: log would drop the negative branch and
    linear would hide everything but the spike. Padding there is multiplicative in magnitude,
    so it works across zero.
    """
    vals = []
    for s in series:
        v = np.asarray(s, dtype=float)[mask]
        v = v[np.isfinite(v)]
        if logy:
            v = v[v > 0.0]
        if v.size:
            vals.append(v)
    if not vals:
        return
    v = np.concatenate(vals)
    vlo, vhi = float(v.min()), float(v.max())
    if symlog:
        f = 1.3
        ax.set_ylim(vlo * f if vlo < 0.0 else vlo / f, vhi * f if vhi > 0.0 else vhi / f)
    elif logy:
        f = 10.0 ** 0.05
        ax.set_ylim(vlo / f, vhi * f)
    else:
        pad = 0.05 * (vhi - vlo) if vhi > vlo else max(abs(vhi) * 0.05, 1e-30)
        ax.set_ylim(vlo - pad, vhi + pad)


def xaxis(d, xlog=None):
    """
    The x-axis transform shared by all three scripts.

    Returns (tx, xlim, xlabel). In log mode everything is plotted against the distance from
    the electrode FACE, which is strictly positive, with the axis inverted so the electrode
    stays on the right exactly as in the linear view.
    """
    xlog = XLOG if xlog is None else xlog
    dx, hi = d["dx"], d["hi"]
    xlo, xhi = d["xlim"]

    def tx(v):
        return (hi - np.asarray(v, dtype=float)) if xlog else np.asarray(v, dtype=float)

    if xlog:
        w_lo = max(hi - xhi, 0.4 * dx)
        w_hi = max(hi - xlo, w_lo * 1.001)
        return tx, (w_hi, w_lo), r"distance from the electrode  $x_{wall}-x$  [m]"
    return tx, (xlo, xhi), "x [m]"


def shade_amr(axes, d, tx, xlog=None):
    """
    Shade the refined regions, one span per contiguous run of EQUAL level.

    Not one per level >= L: nested spans stack their alphas, and at max_level = 6 that sums
    past 1 and paints the sheath -- the one region worth looking at -- solid grey. The spans
    are built from each leaf cell's OWN width, so they are exact at any refinement ratio.
    """
    xlog = XLOG if xlog is None else xlog
    lvl, x, w = d["level"], d["x"], d["dxc"]
    lmax = int(lvl.max())
    if lmax == 0:
        return
    edges = np.nonzero(np.diff(lvl))[0] + 1
    for i0, i1 in zip(np.r_[0, edges], np.r_[edges, len(x)]):
        level = int(lvl[i0])
        if level == 0:
            continue
        xa = float(tx(x[i0] - 0.5 * w[i0]))
        xb = float(tx(x[i1 - 1] + 0.5 * w[i1 - 1]))
        if xlog:
            xa, xb = max(xa, 0.2 * d["dx"]), max(xb, 0.2 * d["dx"])
        for ax in axes:
            ax.axvspan(min(xa, xb), max(xa, xb), alpha=0.30 * level / lmax,
                       color="0.55", lw=0, zorder=0)


# ============================================================
# PLOT
# ============================================================

def make_plot(d, a, out_name):
    x = d["x"]
    tx, xlim, xlabel = xaxis(d)
    xp, mask = tx(x), d["mask"]

    fig, axes = plt.subplots(4, 1, figsize=(10, 14), sharex=True)

    # ---- panel 1: densities + the Boltzmann prediction ----
    ax = axes[0]
    ax.plot(xp, d["NE"], lw=2, color="#2a78d6", label=r"$n_e$")
    ax.plot(xp, d["ARP"], lw=2, color="#d95f02", label=r"$n_i$  (Ar$^+$)")
    if a["edge"] is not None:
        ax.plot(xp, d["NE_BOLTZ"], lw=1.6, ls="--", color="#1b7837",
                label=r"$n_e^{\rm Boltz}=n_{e,\rm edge}e^{(\Phi-\Phi_{\rm edge})/T_{e,\rm edge}}$")
    ax.set_yscale("log")
    ax.set_ylabel(r"density [m$^{-3}$]")
    autoscale(ax, mask, [d["NE"], d["ARP"], d["NE_BOLTZ"]], True)

    # ---- panel 2: potential ----
    ax = axes[1]
    ax.plot(xp, d["PHI"], lw=2, color="#7570b3", label=r"$\Phi$")
    ax.set_ylabel("potential [V]")
    autoscale(ax, mask, [d["PHI"]], False)

    # ---- panel 3: temperatures ----
    ax = axes[2]
    ax.plot(xp, d["ETEMP"], lw=2, color="#d62728", label=r"$T_e$")
    ax.plot(xp, d["T_GAS"], lw=2, color="#2ca02c", label=r"$T_{gas}$ (solved)")
    ax.plot(xp, d["T_EFF"], lw=2, color="#8c564b", label=r"$T_{eff}$ (Wannier, ion)")
    ax.set_ylabel("temperature [K]")
    autoscale(ax, mask, [d["ETEMP"], d["T_GAS"], d["T_EFF"]], False)

    # ---- panel 4: ion velocity against the Bohm speed ----
    # u_i is the FULL ion fluid velocity, drift plus diffusion, off ION_VEL. It is
    # sgn-oriented, so positive = towards the electrode, and it is NOT plotted as |u_i|: the
    # sign is the physically meaningful part, and the crossing of the two curves is the Bohm
    # point marked on every panel. symlog because u_i changes sign and spans decades either
    # way; the linear threshold is the Bohm speed itself.
    ax = axes[3]
    ax.plot(xp, d["U_ION"], lw=2, color="#e7298a",
            label=r"$u_i$  (ION_VEL, +ve towards wall)")
    ax.plot(xp, d["CS"], lw=2, color="#1b9e77", label=r"$c_s=\sqrt{eT_e/M}$  (Bohm speed)")
    ax.axhline(0.0, color="0.5", lw=0.8, zorder=0)
    cs_w = np.asarray(d["CS"], dtype=float)[mask]
    cs_w = cs_w[np.isfinite(cs_w) & (cs_w > 0.0)]
    ax.set_yscale("symlog", linthresh=(float(cs_w.min()) if cs_w.size else 1.0))
    ax.set_ylabel("velocity [m/s]")
    ax.set_xlabel(xlabel)
    autoscale(ax, mask, [d["U_ION"], d["CS"]], False, symlog=True)

    # ---- the two markers, on every panel ----
    for k, ax in enumerate(axes):
        if a["x_edge"] is not None:
            ax.axvline(float(tx(a["x_edge"])), color="k", ls="--", lw=1.5,
                       label=("sheath edge" if k == 0 else None))
        if a["x_bohm"] is not None:
            ax.axvline(float(tx(a["x_bohm"])), color="crimson", ls="-.", lw=1.5,
                       label=("Bohm point  $u_i=c_s$" if k == 0 else None))

    shade_amr(axes, d, tx)

    for ax in axes:
        if XLOG:
            ax.set_xscale("log")
        ax.set_xlim(*xlim)
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(fontsize=9, loc="best")

    fall = a["sheath_fall"]
    title = (f"{os.path.basename(d['file'])}   "
             f"t = {float(d['ds'].current_time):.4e} s   "
             f"finest level {d['finest']}   {d['ncell']} leaf cells")
    if np.isfinite(fall):
        title += f"   sheath fall = {fall:.4f} V"
    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0, 1, 0.975])
    fig.savefig(out_name, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_name}")


# ============================================================
# MAIN
# ============================================================

def main():
    p = argparse.ArgumentParser(description="1-D sheath plots from a PlasmeX plotfile")
    p.add_argument("pattern", help="plotfile, or a glob like 'plt*'")
    args = p.parse_args()

    files = sorted([f for f in glob.glob(args.pattern) if os.path.isdir(f)],
                   key=extract_number)
    if not files:
        raise SystemExit(f"no plotfiles match {args.pattern!r}")

    for fn in files:
        d = sample(fn, want_derived=False)
        d["pres"], d["wall_temp"], d["emission"], d["deck_src"] = read_deck(d["dir"])
        knobs, _ = read_knobs(d["dir"], None)
        d["work_fn"], d["rich_a"] = knobs["prob.work_function"], knobs["prob.richardson_A"]
        d["T_GAS"] = gas_temperature(d, d["pres"])
        d["T_EFF"] = ion_transport(d["T_GAS"], d["AR"], np.abs(d["EFLDX"]))[0]
        set_window(d)

        if "ION_VEL" in d and np.all(d["ION_VEL"] == 0.0):
            print(f"WARNING: {fn} has an all-zero ION_VEL -- was it written on a step where "
                  f"the diagnostics ran? (see the `diagnostics` gate in src/Advance.cpp)")

        a = analyse(d)
        q = electrode_heat_flux(d, a)
        report(d, a, q)
        make_plot(d, a, OUT_PREFIX + os.path.basename(fn) + ".png")


if __name__ == "__main__":
    main()
