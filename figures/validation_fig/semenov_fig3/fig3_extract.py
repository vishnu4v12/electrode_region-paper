"""
Extract the curves of Semenov, Krivtsun & Reisgen (2016) J. Phys. D 49 105204, figure 3(a)
or 3(b) (Ar, 1 atm, anode surface T_s = 1 kK; j = 700 / 1500 A cm^-2; unified model) straight
from the PDF's VECTOR paths -- the same method as ../almeida_fig5/fig5_extract.py.

`mutool trace` dumps every stroked path on page 9 (printed page 8) in PDF points. Each curve
is ONE continuous polyline (the dashing is applied at render time), identified by its CMYK
colour and dash pattern:

    solid      red  n_e      blue n_i         (left axis, log, m^-3)
    dashed     red  T_e      blue T_h         (right axis, linear, kK)
    dash-dot   green n_sa (Saha + quasi-neutrality)

Each axis is calibrated by least squares on the figure's own tick marks (majors and minors);
the residuals are printed, and the extracted curves are checked against the paper's Table 2
(T_e at the anode surface and in the undisturbed plasma) and against the LTE far field.

    python fig3_extract.py a        -> fig3a_extracted.csv, fig3a_verification.{pdf,png}
    python fig3_extract.py b        -> fig3b_extracted.csv, fig3b_verification.{pdf,png}

fig3_compare.py calls extract() itself when the CSV is missing.

x in the paper is the distance from the anode surface (anode on the LEFT).
"""
import os
import subprocess
import sys
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "almeida_fig5"))
from pdftrace import parse  # noqa: E402

PDF_NAME = "Semenov_(Almeida Anode).pdf"
PAGE = 9                                    # printed page 8

UNITS = {"ne": "m^-3", "ni": "m^-3", "nsa": "m^-3", "Te": "K", "Th": "K"}
AXIS_OF = {"ne": "n", "ni": "n", "nsa": "n", "Te": "T", "Th": "T"}

# CMYK, as the PDF writes them
BLACK = (0.0, 0.0, 0.0, 1.0)
RED = (0.0, 0.994, 1.0, 0.0)
BLUE = (0.884, 0.769, 0.0, 0.0)
GREEN = (0.859, 0.174, 1.0, 0.056)
SOLID, DASHED, DASHDOT = "", "6.569 2.986", ".597 2.389 4.18 2.389"


def _decades(lo, hi):
    """Tick values of a log axis from 10^lo to 10^hi: 1..9 x 10^k, then 10^hi."""
    return np.array([m * 10.0 ** k for k in range(lo, hi) for m in range(1, 10)] + [10.0 ** hi])


_AXES = {
    "x": dict(side="bottom", log=True, values=_decades(-8, -2)),
    # labels 10^20 ... 10^23; minor ticks only between them (the frame runs a little beyond)
    "n": dict(side="left", log=True, values=_decades(20, 23)),
    # labels 1, 5, ..., 21 kK (long ticks), short ticks every 2 kK
    "T": dict(side="right_out", log=False, values=np.arange(1, 22, 2) * 1e3),
}
# name, colour, dash -- in the order the PDF draws them (checked, not just trusted)
_CURVES = [("Te", RED, DASHED), ("Th", BLUE, DASHED), ("nsa", GREEN, DASHDOT),
           ("ne", RED, SOLID), ("ni", BLUE, SOLID)]

PANELS = {
    "a": dict(caption="Ar, 1 atm, T_s = 1 kK, j = 700 A/cm^2 (7e6 A/m^2)", j=7e6,
              frame_x=(350.476, 509.923), frame_y=(71.527, 214.551),
              axes=_AXES, curves=_CURVES,
              table2={"Te_wall": 4.93e3, "Te_bulk": 14.37e3},
              crop=(318.0, 540.0, 60.0, 240.0)),
    "b": dict(caption="Ar, 1 atm, T_s = 1 kK, j = 1500 A/cm^2 (1.5e7 A/m^2)", j=1.5e7,
              frame_x=(350.476, 509.923), frame_y=(263.220, 406.245),
              axes=_AXES, curves=_CURVES,
              table2={"Te_wall": 5.59e3, "Te_bulk": 16.82e3},
              crop=(318.0, 540.0, 252.0, 432.0)),
}
# Table 2 of the paper, all rows: j [A/cm^2] -> (T_e at the anode surface, T_e undisturbed) [K]
TABLE2 = {500: (4.69e3, 13.53e3), 700: (4.93e3, 14.37e3), 1000: (5.20e3, 15.40e3),
          1500: (5.59e3, 16.82e3), 2000: (5.99e3, 18.04e3)}
# Table 3: j [A/cm^2] -> U_a, U_1, U_sh [V], L_a [mm], L_sh [um]
TABLE3 = {500: (2.13, 2.05, 0.33, 1.46, 1.564), 1000: (2.63, 2.60, 0.56, 0.84, 1.525),
          2000: (3.18, 2.35, 1.00, 0.40, 1.420)}


def find_pdf(extra=()):
    # the repository's References/ (any level up), beside this script, or in its parent
    ups = [os.path.join(HERE, *[".."] * k) for k in range(4)]
    for d in (*extra, *(os.path.join(u, "References") for u in ups), *ups[:2]):
        p = os.path.join(d, PDF_NAME)
        if os.path.isfile(p):
            return os.path.abspath(p)
    raise SystemExit(f"cannot find '{PDF_NAME}' in References/, beside {HERE} or its parent")


def near(a, b, tol=0.6):
    return abs(a - b) < tol


def fit(pix, val, log):
    """Least-squares page-coord -> value map; returns (to_value, to_page, residuals in pt)."""
    v = np.log10(val) if log else np.asarray(val, float)
    A = np.vstack([pix, np.ones_like(pix)]).T
    (m, c), *_ = np.linalg.lstsq(A, v, rcond=None)
    res_pt = (v - c) / m - pix
    f = (lambda p: 10.0 ** (m * p + c)) if log else (lambda p: m * p + c)
    inv = (lambda q: (np.log10(q) - c) / m) if log else (lambda q: (q - c) / m)
    return f, inv, res_pt


def tick_positions(paths, panel, side):
    """Page positions of every tick on one side of the frame, deduplicated (a major tick is
    drawn twice, long and short), sorted away from the axis origin."""
    (xl, xr), (yt, yb) = panel["frame_x"], panel["frame_y"]
    pos = set()
    for p in paths:
        if len(p["xy"]) != 2 or not np.allclose(p["color"], BLACK):
            continue
        (x0, y0), (x1, y1) = p["xy"]
        if not 0.8 < np.hypot(x1 - x0, y1 - y0) < 6.0:
            continue
        if side == "bottom" and near(x0, x1, .05) and near(min(y0, y1), yb, .1) \
                and xl - 1 < x0 < xr + 1:
            pos.add(round(x0, 2))
        elif near(y0, y1, .05) and yt - 1 < y0 < yb + 1:
            lo, hi = min(x0, x1), max(x0, x1)
            if side == "left" and near(hi, xl, .1) and lo < xl - 0.8:
                pos.add(round(y0, 2))
            elif side == "right_out" and near(lo, xr, .1) and hi > xr + 0.8:
                pos.add(round(y0, 2))
    return np.array(sorted(pos, reverse=(side != "bottom")))


def calibrate(paths, panel, verbose=True):
    cal = {}
    if verbose:
        print("axis calibration (max |residual| of every tick against the fitted map):")
    for name, ax in panel["axes"].items():
        pos = tick_positions(paths, panel, ax["side"])
        assert len(pos) == len(ax["values"]), \
            f"axis {name}: {len(pos)} ticks found, {len(ax['values'])} expected"
        f, inv, res = fit(pos, ax["values"], ax["log"])
        cal[name] = (f, inv)
        if verbose:
            print(f"  {name} ({'log' if ax['log'] else 'lin'})  {len(pos):3d} ticks   "
                  f"max |res| = {np.abs(res).max():.3f} pt")
    return cal


def pick_curves(paths, panel):
    """The five data curves inside the frame, in drawing order, checked against _CURVES.
    Points on the frame edge are the plot clipping the curve (n_sa below 10^20 m^-3 at
    the anode) and are dropped."""
    (xl, xr), (yt, yb) = panel["frame_x"], panel["frame_y"]
    found = [p for p in paths if len(p["xy"]) > 100
             and xl - 0.2 <= p["xy"][:, 0].min() and p["xy"][:, 0].max() <= xr + 0.2
             and yt - 0.2 <= p["xy"][:, 1].min() and p["xy"][:, 1].max() <= yb + 0.2]
    want = panel["curves"]
    assert len(found) == len(want), [(p["color"], p["dash"], len(p["xy"])) for p in found]
    out = {}
    for p, (name, col, dash) in zip(found, want):
        assert np.allclose(p["color"], col) and p["dash"] == dash, (name, p["color"], p["dash"])
        xy = p["xy"]
        xy = xy[(xy[:, 1] < yb - 0.01) & (xy[:, 1] > yt + 0.01)]
        xy = xy[np.argsort(xy[:, 0], kind="stable")]
        out[name] = thin(xy)
    return out


def thin(xy, step=0.08):
    """Drop points closer than `step` pt (along the path) to the last kept one. The source
    polylines carry ~15 000 points each at a few thousandths of a point apart; 0.08 pt is
    well below the line width (0.6 pt), so nothing visible is lost."""
    keep = [0]
    for i in range(1, len(xy)):
        if np.hypot(*(xy[i] - xy[keep[-1]])) >= step:
            keep.append(i)
    if keep[-1] != len(xy) - 1:
        keep.append(len(xy) - 1)
    return xy[keep]


def trace_page(pdf, page=PAGE):
    with tempfile.NamedTemporaryFile(suffix=".xml") as tmp:
        subprocess.run(["mutool", "trace", pdf, str(page)], check=True, stdout=tmp,
                       stderr=subprocess.DEVNULL)
        tmp.flush()
        return parse(tmp.name)[0]


def csv_path(fig):
    return os.path.join(HERE, f"fig3{fig}_extracted.csv")


def extract(fig, pdf=None, verify_fig=True):
    """Extract panel `fig` ('a' or 'b'), run the checks, write the CSV (and the
    verification figure). Returns {name: (x_m, value)}."""
    panel = PANELS[fig]
    pdf = pdf or find_pdf()
    paths = trace_page(pdf)
    print(f"--- Semenov et al. 2016, Fig. 3({fig}): {panel['caption']}")
    cal = calibrate(paths, panel)
    curves = pick_curves(paths, panel)
    data = {n: (cal["x"][0](xy[:, 0]), cal[AXIS_OF[n]][0](xy[:, 1])) for n, xy in curves.items()}
    checks(data, panel)

    path = csv_path(fig)
    with open(path, "w") as f:
        f.write(f"# Semenov, Krivtsun & Reisgen 2016, J. Phys. D 49 105204, Fig. 3({fig}): "
                f"{panel['caption']}\n")
        f.write("# extracted from the PDF vector paths (fig3_extract.py); "
                "x = distance from anode\n")
        f.write("quantity,unit,x_m,value\n")
        for name, (xs, vs) in data.items():
            for a, b in zip(xs, vs):
                f.write(f"{name},{UNITS[name]},{a:.6e},{b:.6e}\n")
    print(f"wrote {path}  ({sum(len(v[0]) for v in data.values())} points, {len(data)} curves)")
    if verify_fig:
        verification_figure(fig, pdf, data, cal)
    return data


def load(fig, pdf=None):
    """The extracted curves of Fig. 3(fig), extracting first if the CSV is not there."""
    path = csv_path(fig)
    if not os.path.isfile(path):
        return extract(fig, pdf)
    cols = {}
    with open(path) as f:
        for line in f:
            if line.startswith(("#", "quantity")):
                continue
            q, _, x, v = line.strip().split(",")
            cols.setdefault(q, ([], []))
            cols[q][0].append(float(x))
            cols[q][1].append(float(v))
    return {q: (np.array(x), np.array(v)) for q, (x, v) in cols.items()}


def checks(data, panel):
    t2 = panel["table2"]
    xe, te = data["Te"]
    xh, th = data["Th"]
    print("Table 2 vs extraction:")
    print(f"  T_e at the anode   paper {t2['Te_wall'] / 1e3:.2f} kK   extracted at "
          f"x = {xe[0]:.1e} m: {te[0] / 1e3:.2f} kK")
    print(f"  T_e undisturbed    paper {t2['Te_bulk'] / 1e3:.2f} kK   extracted at "
          f"x = {xe[-1]:.1e} m: {te[-1] / 1e3:.2f} kK")
    print(f"  T_h at the anode   T_s = 1.00 kK         extracted: {th[0] / 1e3:.2f} kK")
    far = {k: np.interp(np.log(5e-3), np.log(v[0]), v[1]) for k, v in data.items()}
    print(f"  far field (x = 5 mm): n_i/n_e - 1 = {far['ni'] / far['ne'] - 1:+.3f},  "
          f"n_sa/n_e - 1 = {far['nsa'] / far['ne'] - 1:+.3f},  T_h/T_e - 1 = "
          f"{far['Th'] / far['Te'] - 1:+.3f}")


def verification_figure(fig, pdf, data, cal):
    """(a) the original panel, rasterised from the PDF; (b) the same raster, greyed, with
    every extracted curve drawn back over it through the inverse of the fitted axis maps."""
    panel = PANELS[fig]
    dpi = 600
    x0, x1, y0, y1 = panel["crop"]
    s = dpi / 72.0
    with tempfile.TemporaryDirectory() as td:
        stem = os.path.join(td, "crop")
        subprocess.run(["pdftoppm", "-f", str(PAGE), "-l", str(PAGE), "-r", str(dpi), "-png",
                        "-singlefile", "-x", str(int(x0 * s)), "-y", str(int(y0 * s)),
                        "-W", str(int((x1 - x0) * s)), "-H", str(int((y1 - y0) * s)),
                        pdf, stem], check=True)
        img = plt.imread(stem + ".png")
    ext = [int(x0 * s) / s, int(x0 * s) / s + img.shape[1] / s,
           int(y0 * s) / s + img.shape[0] / s, int(y0 * s) / s]

    # Register the raster on the frame itself (pdftoppm and mutool may place the page by
    # different boxes): find the bottom and left frame lines in the image by an ink-weighted
    # centroid and shift the raster so they land on the traced frame.
    ink = 1.0 - img[..., :3].mean(axis=2)
    (fxl, fxr), (fyt, fyb) = panel["frame_x"], panel["frame_y"]

    def centroid(profile, guess, half=15):
        i = np.arange(max(guess - half, 0), guess + half)
        w = profile[i]
        return (i * w).sum() / w.sum() + 0.5

    # Probe where no curve meets the frame: the bottom edge at 15 % of the width (n_sa
    # reaches the bottom edge near x = 1e-5 m), the left edge at 5 % of the height (above the last
    # tick, 10^23).
    col = int((fxl + 0.15 * (fxr - fxl) - ext[0]) * s)
    row = int((fyt + 0.05 * (fyb - fyt) - ext[3]) * s)
    y_found = ext[3] + centroid(ink[:, col], int((fyb - ext[3]) * s)) / s
    x_found = ext[0] + centroid(ink[row, :], int((fxl - ext[0]) * s)) / s
    dx, dy = fxl - x_found, fyb - y_found
    print(f"raster registration: dx = {dx:+.2f} pt, dy = {dy:+.2f} pt "
          f"(page-box offset; the extracted data are unaffected)")
    ext = (ext[0] + dx, ext[1] + dx, ext[2] + dy, ext[3] + dy)

    style = {"ne": "#d62728", "ni": "#1f3fff", "nsa": "#2ca02c", "Te": "#ff7f0e",
             "Th": "#17becf"}
    plt.rcParams.update({"font.family": "serif", "font.serif": ["STIX", "STIXGeneral",
                         "DejaVu Serif"], "mathtext.fontset": "stix", "font.size": 8})
    f, axs = plt.subplots(1, 2, figsize=(7.0, 3.1))
    axs[0].imshow(img, extent=ext)
    axs[0].set_title(f"(a) Semenov et al. (2016), Fig. 3({fig}) — original", fontsize=8.5)
    grey = img[..., :3].mean(axis=2)
    axs[1].imshow(0.55 + 0.45 * grey, cmap="gray", vmin=0, vmax=1, extent=ext)
    for name, (xs, vs) in data.items():
        axs[1].plot(cal["x"][1](xs), cal[AXIS_OF[name]][1](vs), lw=0.9, color=style[name],
                    label={"nsa": "$n_{sa}$"}.get(name, name))
    axs[1].set_title("(b) extracted data re-plotted over the original (greyed)", fontsize=8.5)
    axs[1].legend(loc="lower left", bbox_to_anchor=(1.0, 0.0), fontsize=7, frameon=False,
                  handlelength=1.6)
    for ax in axs:
        ax.set_xlim(ext[0], ext[1])
        ax.set_ylim(ext[2], ext[3])
        ax.set_axis_off()
    f.tight_layout()
    for e in ("pdf", "png"):
        f.savefig(os.path.join(HERE, f"fig3{fig}_verification.{e}"), dpi=300,
                  bbox_inches="tight")
    plt.close(f)
    print(f"wrote fig3{fig}_verification.pdf/.png")


if __name__ == "__main__":
    for fg in sys.argv[1:] or ["a", "b"]:
        extract(fg)
