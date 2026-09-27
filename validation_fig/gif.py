"""
Animate a 1-D PlasmeX plotfile series into a GIF.

    python gif.py 'plt*'            every plotfile matching the glob
    python gif.py 'plt*' 5          every 5th one
    python gif.py 'plt*' 1 out.gif  explicit output name

One stacked figure per frame: the groups configured in PLOT_GROUPS below, plus an AMR-level
panel showing where the mesh is refined.

Data comes straight off the AMR LEAF CELLS via bohm.sample() -- see bohm.py's module docstring
for what that means and why. Three things follow that matter here:

  * ANY REFINEMENT RATIO. The leaf walk reads each grid's own dds, so amr.ref_ratio may be
    2, 4, or a different value per level. Nothing computes ref_ratio**(finest-level).

  * FAST. This used to sample with yt.LineBuffer, which re-traces the ray once PER FIELD and
    resamples onto a uniform fine grid: 8.3 s per frame on the shipped 4-level deck, plotting
    32768 points of which 32124 were duplicates of a coarse cell. The leaf walk reads each
    grid once and hands matplotlib ~640 real cells. Only the fields PLOT_GROUPS actually names
    are read out of the plotfile.

  * NO EXTRAPOLATION AT THE EDGES. The first and last points of every curve are the boundary
    LEAF CELL values, verbatim. For JX / JEX / JIX that cell was overwritten by the solver
    with the domain FACE value (BoundaryFaceOverwrite, src/aider.H), so those curves end on
    the flux that actually crossed the boundary; every other field ends on an ordinary
    cell-centred value. The legend tags each curve `[FACE]` or `[cell]` and the per-frame
    console line prints both edge values, so which one is being read is never a guess.

ELECTRON ENERGY FLUX PANEL (PLOT_EEN_FLUX). The plotfile carries only the EXPLICIT part of the
electron energy flux, and only its divergence (DIVEENEXPLFLX); the implicit conduction half never
reaches disk. This panel rebuilds the TOTAL flux on the interior faces (bohm.electron_energy_flux:
enthalpy + Onsager (5/2 + C_e) k T_e Gamma_e with the solver's own Gamma_e = -JEX/e, conduction
-lambda grad T_e, and the heavy-velocity half of Almeida A-(15)), and overlays the flux the
SOURCES imply -- the exe/bc.H wall value minus the integral of EENRXN from the wall. At steady
state the two coincide (div Q = S); their gap is d(EEN)/dt. The console prints that gap.

1-D ONLY, by assertion in bohm.sample. yt's cylindrical coordinate handler has no
pixelize_line, so the old LineBuffer path raised a MESSAGELESS NotImplementedError on every
r-z plotfile and produced a GIF of completely empty axes with nothing obviously wrong. Use
gif_2d.py for r-z.
"""

# ============================================================
# USER CONFIG
# ============================================================

# ---- x axis ------------------------------------------------------------------------
# Plot window in metres. None takes the domain edge from the plotfile.
PROB_LO = None
PROB_HI = None

# Log x axis, in DISTANCE FROM THE ELECTRODE (a plain log10(x) would collapse the sheath into
# the last decade). The axis is inverted so the electrode stays on the right.
XLOG = False

# ---- what to plot ------------------------------------------------------------------
# One panel per group. "vars" are plotfile component names; only these are read.
#   mode        "individual" (one curve each) or "sum" (one curve, the sum)
#   logy        log y axis
#   ylabel      y axis label
#   gas_temp    also draw T_gas -- the solved TGAS component (see bohm.gas_temperature)
PLOT_GROUPS = [
    {"vars": ["NE", "ARP"], "mode": "individual", "logy": True,  "ylabel": "Density [m$^{-3}$]"},
    {"vars": ["AR"], "mode": "individual", "logy": False,  "ylabel": "Density [m$^{-3}$]"},
    # VION is not a plotfile field -- it's derived below as JIX/(e*ARP), the ion charge
    # being +e, so u_i = J_i/(e n_i) in the same coordinate-component convention as UNX/HVELX.
    #{"vars": ["UNX", "HVELX", "VION"], "mode": "individual", "logy": False,  "ylabel": "velocity [m/s]"},
    {"vars": [ "HVELX"], "mode": "individual", "logy": False,  "ylabel": "velocity [m/s]"},
    {"vars": ["PRES"], "mode": "individual", "logy": False,  "ylabel": "Pressure [Pa]$]"},
    {"vars": ["PHI"],       "mode": "individual", "logy": False, "ylabel": "Potential [V]"},
    {"vars": ["JEX", "JIX"],"mode": "individual", "logy": False, "ylabel": "Current [A/m$^2$]"},
    {"vars": ["ETEMP","TGAS"],     "mode": "individual", "logy": False, "ylabel": "Temperature [K]"},
]

# ---- electron energy flux panel -----------------------------------------------------
# The TOTAL electron energy flux [W/m^2, + = towards the electrode, x-hi] -- the flux whose
# divergence the electron-energy sources (EENRXN: Joule, inelastic, elastic exchange) balance --
# rebuilt on the interior faces by bohm.electron_energy_flux. One extra panel above the AMR one.
PLOT_EEN_FLUX = True
# Which curves:
#   "total"       enthalpy + conduction + heavy, the whole flux
#   "enthalpy"    (5/2 + C_e) k T_e Gamma_e, Gamma_e = -JEX/e (the solver's electron flux)
#   "heat"        the heat-transfer terms lumped: "conduction" + "heavy"
#   "conduction"  -lambda grad T_e, lambda = K n_e k D_e (M&K/Almeida, src/coulomb.H)
#   "heavy"       -k T_e n_e [A_i (u_i - u) + A_a (u_n - u)], the heavy-velocity Onsager half
#   "source"      Q_wall - integral of EENRXN from the wall: the flux the sources imply at
#                 STEADY STATE (d EEN/dt dropped). It lies on "total" when the run is steady.
EEN_FLUX_PARTS = ["total", "enthalpy", "heat", "source"]
# The deck the run used, read for plasmex.thermal_diffusion / electron_kappa /
# stefan_maxwell_drag and prob.wall_temp / electrode_emission (beside the plotfiles first, then
# beside this script). A knob the deck does not set takes the CODE default (off).
DECK_FILE = "cathode_inputs"
# The flux changes sign near the wall (the Te maximum); symlog shows both branches.
EEN_FLUX_SYMLOG = False

# ---- y axis across the animation ----------------------------------------------------
# True  -> every frame shares one y range per panel, taken from a first pass over all the
#          files, so the animation shows the solution moving rather than the axes moving.
# False -> each frame autoscales to itself.
# The first pass costs one extra read of each plotfile and is usually worth it.
FIXED_YLIM = True

# ---- output ------------------------------------------------------------------------
GIF_NAME = "animation.gif"
FRAME_DIR = "frames"
FRAME_DURATION = 0.2     # seconds per frame
DPI = 110


# ============================================================
# IMPORTS
# ============================================================

import glob
import importlib.util
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import imageio.v2 as imageio


HERE = os.path.dirname(os.path.abspath(__file__))

# Fields computed here rather than read off the plotfile, and the label to use for them
# instead of "NAME [face_tag]" (b.face_tag doesn't know about them).
DERIVED_LABELS = {"VION": "ion velocity  J_i/(e n_i)"}

# The electron energy flux panel: EEN_FLUX_PARTS name -> (key in bohm.electron_energy_flux's
# dict, legend label, line style).
EEN_CURVES = {
    "total":      ("TOT",   r"total  $Q_e$",                                   dict(lw=2.4, color="k")),
    "enthalpy":   ("ENTH",  r"$(5/2 + C_e)\,kT_e\,\Gamma_e$",                   dict(lw=1.6)),
    "heat":       ("HEAT",  r"heat  $-\lambda_e \nabla T_e - kT_e n_e[A_i(u_i-u)+A_a(u_n-u)]$", dict(lw=1.6)),
    "conduction": ("COND",  r"$-\lambda_e \nabla T_e$",                          dict(lw=1.6)),
    "heavy":      ("HEAVY", r"$-kT_e n_e[A_i(u_i-u)+A_a(u_n-u)]$",              dict(lw=1.2)),
    "source":     ("SRC",   r"$Q_{wall} - \int$ EENRXN dx  (steady state)",     dict(lw=1.6, ls="--", color="#d62728")),
}
# The plotfile components the flux needs (JIX/UNX/HVELX are optional: 0 when absent).
EEN_FIELDS = ("NE", "ARP", "AR", "ETEMP", "JEX", "JIX", "UNX", "HVELX", "EENRXN", "EFLDX")


def een_series(d, part):
    """(x, values) of one EEN_FLUX_PARTS curve, or None -- faces, not cells."""
    q = d.get("QE")
    if q is None or part not in EEN_CURVES:
        return None
    key = EEN_CURVES[part][0]
    return (q["xf_src"] if key == "SRC" else q["xf"]), q[key]


def een_mask(d, xs):
    """The plot-window mask for a face array (bohm.set_window only masks the cells)."""
    lo, hi = d["xlim"]
    return (xs >= lo) & (xs <= hi)


def load_bohm(exe_dir):
    """Import bohm.py by path -- beside the plotfiles first, then beside this script."""
    for d in (exe_dir, HERE):
        path = os.path.join(d, "bohm.py")
        if os.path.isfile(path):
            modname = "bohm_" + re.sub(r"\W+", "_", os.path.abspath(d))
            if modname in sys.modules:
                return sys.modules[modname]
            spec = importlib.util.spec_from_file_location(modname, path)
            mod = importlib.util.module_from_spec(spec)
            sys.modules[modname] = mod
            spec.loader.exec_module(mod)
            return mod
    raise SystemExit(f"cannot find bohm.py in {exe_dir} or {HERE}")


# ============================================================
# ONE FRAME
# ============================================================

def frame(b, d, ylims, out_name):
    """
    Draw one frame. `ylims` is a list of (lo, hi) per group, or None to autoscale per frame.
    """
    x = d["x"]
    tx, xlim, xlabel = b.xaxis(d, XLOG)
    xp, mask = tx(x), d["mask"]

    een = PLOT_EEN_FLUX and d.get("QE") is not None
    nplots = len(PLOT_GROUPS) + (1 if een else 0) + 1
    fig, axes = plt.subplots(nplots, 1, figsize=(9, 3.1 * nplots), sharex=True, squeeze=False)
    axes = axes[:, 0]

    for k, (ax, group) in enumerate(zip(axes[:-1], PLOT_GROUPS)):
        series = []

        if group.get("mode", "individual") == "sum":
            present = [v for v in group["vars"] if v in d]
            if present:
                tot = np.sum([np.asarray(d[v], dtype=float) for v in present], axis=0)
                ax.plot(xp, tot, lw=2, label=group.get("sum_label", " + ".join(present)))
                series.append(tot)
        else:
            for name in group["vars"]:
                if name not in d:
                    continue
                v = np.asarray(d[name], dtype=float)
                label = DERIVED_LABELS.get(name, f"{name} {b.face_tag(name)}")
                ax.plot(xp, v, lw=2, label=label)
                series.append(v)

        if group.get("gas_temp") and "T_GAS" in d:
            ax.plot(xp, d["T_GAS"], lw=1.6, ls="--", color="#2ca02c",
                    label="T_gas (LTE from AR) [cell]")
            series.append(d["T_GAS"])

        if group.get("logy"):
            ax.set_yscale("log")
        ax.set_ylabel(group.get("ylabel", ""))
        if ylims is not None and ylims[k] is not None:
            ax.set_ylim(*ylims[k])
        elif series:
            b.autoscale(ax, mask, series, bool(group.get("logy")))

    # ---- electron energy flux panel ----
    # Faces, not cells: each curve has its own x array and window mask.
    if een:
        ax = axes[len(PLOT_GROUPS)]
        series, masks = [], []
        for part in EEN_FLUX_PARTS:
            got = een_series(d, part)
            if got is None:
                continue
            xs, v = got
            ax.plot(tx(xs), v, label=EEN_CURVES[part][1], **EEN_CURVES[part][2])
            series.append(v)
            masks.append(een_mask(d, xs))
        ax.axhline(0.0, color="#888888", lw=0.8)
        if EEN_FLUX_SYMLOG:
            vmax = max((np.nanmax(np.abs(v[m])) for v, m in zip(series, masks) if m.any()), default=1.0)
            ax.set_yscale("symlog", linthresh=1e-3 * vmax)
        ax.set_ylabel("e energy flux [W/m$^2$]\n(+ = towards electrode)")
        k_een = len(PLOT_GROUPS)
        if ylims is not None and ylims[k_een] is not None:
            ax.set_ylim(*ylims[k_een])
        elif series:
            vals = np.concatenate([v[m] for v, m in zip(series, masks)])
            b.autoscale(ax, np.ones(vals.size, dtype=bool), [vals], False, symlog=EEN_FLUX_SYMLOG)

    # ---- AMR level panel ----
    # Drawn as a step over each leaf cell's OWN edges, so it is exact at any refinement ratio.
    ax = axes[-1]
    lvl, w = d["level"], d["dxc"]
    # Cell EDGES, not centres. In XLOG the last edge is the electrode face itself, where the
    # transformed coordinate is exactly 0 and a log axis has nothing to draw, so it is pulled
    # in to the same 0.2*dx floor shade_amr uses.
    xe = tx(np.r_[x - 0.5 * w, x[-1] + 0.5 * w[-1]])
    if XLOG:
        xe = np.maximum(xe, 0.2 * d["dx"])
    ax.step(xe, np.r_[lvl, lvl[-1]], where="post", lw=2, color="#444444")
    ax.set_yticks(np.arange(int(lvl.max()) + 1))
    ax.set_ylim(-0.25, int(lvl.max()) + 0.5)
    ax.set_ylabel("AMR level")
    ax.set_xlabel(xlabel)

    b.shade_amr(axes[:-1], d, tx, XLOG)

    for ax in axes:
        if XLOG:
            ax.set_xscale("log")
        ax.set_xlim(*xlim)
        ax.grid(True, which="both", alpha=0.3)
        if ax is not axes[-1]:
            ax.legend(fontsize=8.5, loc="best")

    fig.suptitle(f"{os.path.basename(d['file'])}   "
                 f"t = {float(d['ds'].current_time):.4e} s   "
                 f"finest level {d['finest']}   {d['ncell']} leaf cells")
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(out_name, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


# ============================================================
# MAIN
# ============================================================

def main():
    if len(sys.argv) < 2:
        raise SystemExit("Usage: python gif.py 'plt*' [stride] [output.gif]")

    pattern = sys.argv[1]
    stride = int(sys.argv[2]) if len(sys.argv) >= 3 else 1
    gif_name = sys.argv[3] if len(sys.argv) >= 4 else GIF_NAME

    files = [f for f in glob.glob(pattern) if os.path.isdir(f)]
    if not files:
        raise SystemExit(f"no plotfiles match {pattern!r}")

    b = load_bohm(os.path.dirname(os.path.abspath(files[0])))
    b.PROB_LO, b.PROB_HI, b.XLOG = PROB_LO, PROB_HI, XLOG
    files = sorted(files, key=b.extract_number)[::stride]

    # Only the components the panels name, plus AR when a gas temperature is wanted.
    wanted = [v for g in PLOT_GROUPS for v in g["vars"] if v not in DERIVED_LABELS]
    if any(g.get("gas_temp") for g in PLOT_GROUPS):
        # TGAS is the solved field; AR is still wanted so the LTE fallback in
        # bohm.gas_temperature() can serve a plotfile written before TGAS existed.
        wanted.extend(("TGAS", "AR"))
    if any("VION" in g["vars"] for g in PLOT_GROUPS):
        # VION = JIX/(e*ARP); JIX is FACE-overwritten at the boundary cell (BoundaryFaceOverwrite,
        # src/aider.H), ARP is the ordinary cell density -- same mixed convention bohm.py's own
        # electrode_heat_flux() uses for the ion channel.
        wanted.extend(("JIX", "ARP"))
    if PLOT_EEN_FLUX:
        wanted.extend(EEN_FIELDS)
    fields = tuple(dict.fromkeys(wanted))

    print(f"{len(files)} frame(s), reading {', '.join(fields)}")

    def read(fn):
        d = b.sample(fn, fields, want_derived=False)
        tg = b.gas_temperature(d, b.read_deck(d["dir"])[0])
        if tg is not None:
            d["T_GAS"] = tg
        if "JIX" in d and "ARP" in d:
            d["VION"] = d["JIX"] / (b.ECHARGE * d["ARP"])
        d = b.set_window(d)
        if PLOT_EEN_FLUX:
            knobs, path = b.read_knobs(d["dir"], DECK_FILE)
            if path is None and not warned:
                print(f"WARNING: {DECK_FILE} not found beside the plotfiles or gif.py -- the "
                      f"energy flux uses the CODE defaults (thermal_diffusion, electron_kappa, "
                      f"stefan_maxwell_drag all OFF), which the shipped decks do not run with")
                warned.append(True)
            d["QE"] = b.electron_energy_flux(d, knobs)
            if d["QE"] is None and not warned:
                print("WARNING: no electron energy flux panel -- ar_e_lorentz.dat not found or "
                      "JEX/EENRXN absent from the plotfile")
                warned.append(True)
        return d

    warned = []

    # ---- pass 1: common y limits ----
    # One read per file, then the datasets are dropped; holding every frame's arrays would
    # scale memory with the length of the series for no benefit.
    ylims = None
    if FIXED_YLIM and len(files) > 1:
        acc = [[np.inf, -np.inf] for _ in PLOT_GROUPS]
        acc_een = [np.inf, -np.inf]
        for fn in files:
            d = read(fn)
            m = d["mask"]
            for part in (EEN_FLUX_PARTS if PLOT_EEN_FLUX else ()):
                got = een_series(d, part)
                if got is None:
                    continue
                v = got[1][een_mask(d, got[0])]
                v = v[np.isfinite(v)]
                if v.size:
                    acc_een[0] = min(acc_een[0], float(v.min()))
                    acc_een[1] = max(acc_een[1], float(v.max()))
            for k, group in enumerate(PLOT_GROUPS):
                names = list(group["vars"]) + (["T_GAS"] if group.get("gas_temp") else [])
                for name in names:
                    if name not in d:
                        continue
                    v = np.asarray(d[name], dtype=float)[m]
                    v = v[np.isfinite(v)]
                    if group.get("logy"):
                        v = v[v > 0.0]
                    if v.size:
                        acc[k][0] = min(acc[k][0], float(v.min()))
                        acc[k][1] = max(acc[k][1], float(v.max()))
        ylims = []
        for k, group in enumerate(PLOT_GROUPS):
            vlo, vhi = acc[k]
            if not np.isfinite(vlo) or not np.isfinite(vhi):
                ylims.append(None)
            elif group.get("logy"):
                f = 10.0 ** 0.05
                ylims.append((vlo / f, vhi * f))
            else:
                pad = 0.05 * (vhi - vlo) if vhi > vlo else max(abs(vhi) * 0.05, 1e-30)
                ylims.append((vlo - pad, vhi + pad))
        # The flux panel's slot, after the groups (frame() indexes it as len(PLOT_GROUPS)).
        vlo, vhi = acc_een
        if not (np.isfinite(vlo) and np.isfinite(vhi)):
            ylims.append(None)
        elif EEN_FLUX_SYMLOG:
            ylims.append((vlo * 1.3 if vlo < 0.0 else vlo / 1.3, vhi * 1.3 if vhi > 0.0 else vhi / 1.3))
        else:
            pad = 0.05 * (vhi - vlo) if vhi > vlo else max(abs(vhi) * 0.05, 1e-30)
            ylims.append((vlo - pad, vhi + pad))

    # ---- pass 2: the frames ----
    os.makedirs(FRAME_DIR, exist_ok=True)
    names = []
    for i, fn in enumerate(files):
        d = read(fn)
        out = os.path.join(FRAME_DIR, f"frame_{i:04d}.png")
        frame(b, d, ylims, out)
        names.append(out)

        # Edge values, printed so a frame's boundary numbers can be read without opening the
        # plot -- and tagged, because the ones marked [FACE] are the flux the solver applied
        # across the boundary while the rest are ordinary cell averages.
        edges = "  ".join(
            f"{v}{'' if v in DERIVED_LABELS else b.face_tag(v)} "
            f"{float(d[v][0]):.3e}/{float(d[v][-1]):.3e}"
            for g in PLOT_GROUPS for v in g["vars"] if v in d)
        print(f"  [{i + 1}/{len(files)}] {os.path.basename(fn)}  "
              f"t = {float(d['ds'].current_time):.4e} s  lvl<={d['finest']}  "
              f"{d['ncell']} cells   arc/wall: {edges}")
        q = d.get("QE")
        if q is not None:
            # Closure of div Q = S over the window, normalised by the largest |Q_src| there --
            # NOT pointwise, which blows up where the flux crosses zero near the wall.
            mq = een_mask(d, q["xf"])
            if mq.any():
                gap = np.abs(q["TOT"] - q["SRC"][:-1])[mq]
                scale = float(np.abs(q["SRC"][:-1][mq]).max())
                j = int(np.argmax(gap))
                print(f"      e energy flux: Q_wall {q['q_wall']:.4e} W/m^2   "
                      f"max|Q_e - Q_src|/max|Q_src| = {gap[j] / scale:.3e} "
                      f"at x = {q['xf'][mq][j]:.6e} m")

    # ---- the GIF ----
    # Every frame is drawn at the same figsize and DPI with bbox_inches="tight", which can
    # still differ by a pixel when a tick label grows. imageio needs identical shapes, so the
    # frames are cropped to the smallest common size rather than resampled -- a crop of a
    # couple of pixels of margin is invisible, a LANCZOS resize of a line plot is not.
    imgs = [imageio.imread(f) for f in names]
    h = min(im.shape[0] for im in imgs)
    w = min(im.shape[1] for im in imgs)
    imgs = [im[:h, :w] for im in imgs]

    imageio.mimsave(gif_name, imgs, duration=FRAME_DURATION, loop=0)
    print(f"wrote {gif_name}  ({len(imgs)} frames, {w}x{h})")


if __name__ == "__main__":
    main()
