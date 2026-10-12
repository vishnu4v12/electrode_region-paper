"""
Extract the j = 1500 A/cm^2 (1.5e7 A/m^2) curves of Semenov, Krivtsun & Reisgen (2016)
J. Phys. D 49 105204 (Ar, 1 atm, anode surface T_s = 1 kK, unified model) straight from the
PDF's VECTOR paths:

    Fig. 3(b), page 9    n_e, n_i, n_sa, T_e, T_h          -> fig3b_extracted.csv
    Fig. 4,    page 10   j_i/j  (red curve; blue is 700)   -> fig45_j1500_extracted.csv
    Fig. 5,    page 10   E      (red curve)                -> fig45_j1500_extracted.csv

    python extract.py [--pdf paper.pdf]

Also writes fig3b_verification.{pdf,png} and fig45_j1500_verification.{pdf,png}: the original
panel, and the same raster greyed with every extracted curve drawn back over it.

`mutool trace` dumps every stroked path on a page in PDF points. Each curve is ONE continuous
polyline (the dashing is applied at render time), identified by its CMYK colour and dash
pattern. In Fig. 3:

    solid      red  n_e      blue n_i         (left axis, log, m^-3)
    dashed     red  T_e      blue T_h         (right axis, linear, kK)
    dash-dot   green n_sa (Saha + quasi-neutrality)

Each axis is calibrated by least squares on the figure's own tick marks (majors and minors);
the residuals are printed, and the extracted curves are checked against the paper's Table 2
(T_e at the anode surface and in the undisturbed plasma) and against the LTE far field.

x is the distance from the anode surface (anode on the LEFT). j_i/j and E carry the paper's
signs: x points away from the anode and j > 0, so ions flowing INTO the anode give j_i/j < 0,
and E < 0 near the anode (the paper: E > 0 in the uniform plasma, changing sign inside the
layer). The paper's j_e/j (Fig. 4, right axis) is j_i/j - 1 (J_i - J_e = j/e).

Needs mutool and pdftoppm. The PDF is looked for in the repository's References/.
"""
import argparse
import os
import subprocess
import sys
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from pdftrace import parse  # noqa: E402

PDF_NAME = "Semenov_(Almeida Anode).pdf"
CAPTION = "Ar, 1 atm, T_s = 1 kK, j = 1500 A/cm^2 (1.5e7 A/m^2)"

# CMYK, as the PDF writes them
BLACK = (0.0, 0.0, 0.0, 1.0)
RED = (0.0, 0.994, 1.0, 0.0)
BLUE = (0.884, 0.769, 0.0, 0.0)
GREEN = (0.859, 0.174, 1.0, 0.056)
SOLID, DASHED, DASHDOT = "", "6.569 2.986", ".597 2.389 4.18 2.389"


def _decades(lo, hi):
    """Tick values of a log axis from 10^lo to 10^hi: 1..9 x 10^k, then 10^hi."""
    return np.array([m * 10.0 ** k for k in range(lo, hi) for m in range(1, 10)] + [10.0 ** hi])


# Each panel: its page, its frame in page points (y measured from the TOP of the page), the
# tick values of each axis in the order the ticks are met going away from the axis origin,
# the crop for the verification raster and a probe point (fractions of the frame) where no
# curve or tick meets the frame, used to register the raster.
FIG3 = dict(
    page=9,                                  # printed page 8
    frame_x=(350.476, 509.923), frame_y=(263.220, 406.245),
    axes={
        "x": dict(side="bottom", log=True, values=_decades(-8, -2)),
        # labels 10^20 ... 10^23; minor ticks only between them (the frame runs a little beyond)
        "n": dict(side="left", log=True, values=_decades(20, 23)),
        # labels 1, 5, ..., 21 kK (long ticks), short ticks every 2 kK
        "T": dict(side="right_out", log=False, values=np.arange(1, 22, 2) * 1e3),
    },
    # name, colour, dash -- in the order the PDF draws them (checked, not just trusted)
    curves=[("Te", RED, DASHED), ("Th", BLUE, DASHED), ("nsa", GREEN, DASHDOT),
            ("ne", RED, SOLID), ("ni", BLUE, SOLID)],
    crop=(318.0, 540.0, 252.0, 432.0),
    # the bottom edge at 15 % of the width (n_sa reaches the bottom edge near x = 1e-5 m),
    # the left edge at 5 % of the height (above the last tick, 10^23)
    probe=(0.15, 0.05),
)
FIG3_UNITS = {"ne": "m^-3", "ni": "m^-3", "nsa": "m^-3", "Te": "K", "Th": "K"}
FIG3_AXIS = {"ne": "n", "ni": "n", "nsa": "n", "Te": "T", "Th": "T"}
TABLE2 = {"Te_wall": 5.59e3, "Te_bulk": 16.82e3}      # Table 2, 1500 A/cm^2

FIG45 = {
    "ji_j": dict(
        fig="Fig. 4", unit="-", page=10,     # printed page 9
        frame_x=(94.45, 247.52), frame_y=(197.65, 318.12),
        axes={
            "x": dict(side="bottom", log=True, values=_decades(-8, -2)),
            # labels 0, -0.01, -0.02; ticks every 0.005
            "y": dict(side="left", log=False, values=np.arange(-4, 1) * 0.005),
        },
        crop=(52.0, 262.0, 183.0, 340.0),
        probe=(0.469, 0.393),                # between ticks (mid-frame is on 1e-5 and -0.01)
    ),
    "E": dict(
        fig="Fig. 5", unit="V/m", page=10,
        frame_x=(356.4, 515.84), frame_y=(78.02, 221.05),
        axes={
            "x": dict(side="bottom", log=True, values=_decades(-8, -2)),
            # labels 0, -0.5, ..., -2.0 (10^6 V/m); ticks every 0.25
            "y": dict(side="left", log=False, values=np.arange(-8, 1) * 0.25e6),
        },
        crop=(318.0, 530.0, 62.0, 243.0),
        probe=(0.457, 0.5),                  # between ticks (mid-frame is on 1e-5)
    ),
}


# ---------------------------------------------------------------- PDF tracing and axes

def find_pdf():
    # the repository's References/ (any level up)
    for k in range(5):
        p = os.path.join(HERE, *[".."] * k, "References", PDF_NAME)
        if os.path.isfile(p):
            return os.path.abspath(p)
    raise SystemExit(f"cannot find '{PDF_NAME}' in References/ -- pass --pdf")


def trace_page(pdf, page):
    with tempfile.NamedTemporaryFile(suffix=".xml") as tmp:
        subprocess.run(["mutool", "trace", pdf, str(page)], check=True, stdout=tmp,
                       stderr=subprocess.DEVNULL)
        tmp.flush()
        return parse(tmp.name)[0]


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


def calibrate(paths, panel):
    cal = {}
    print("axis calibration (max |residual| of every tick against the fitted map):")
    for name, ax in panel["axes"].items():
        pos = tick_positions(paths, panel, ax["side"])
        assert len(pos) == len(ax["values"]), \
            f"axis {name}: {len(pos)} ticks found, {len(ax['values'])} expected"
        f, inv, res = fit(pos, ax["values"], ax["log"])
        cal[name] = (f, inv)
        print(f"  {name} ({'log' if ax['log'] else 'lin'})  {len(pos):3d} ticks   "
              f"max |res| = {np.abs(res).max():.3f} pt")
    return cal


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


def in_frame(p, panel):
    (xl, xr), (yt, yb) = panel["frame_x"], panel["frame_y"]
    return (len(p["xy"]) > 100
            and xl - 0.2 <= p["xy"][:, 0].min() and p["xy"][:, 0].max() <= xr + 0.2
            and yt - 0.2 <= p["xy"][:, 1].min() and p["xy"][:, 1].max() <= yb + 0.2)


def write_csv(path, header, units, data):
    with open(path, "w") as f:
        for h in header:
            f.write(f"# {h}\n")
        f.write("quantity,unit,x_m,value\n")
        for name, (xs, vs) in data.items():
            for a, b in zip(xs, vs):
                f.write(f"{name},{units[name]},{a:.6e},{b:.6e}\n")
    print(f"wrote {path}  ({sum(len(v[0]) for v in data.values())} points, {len(data)} curves)")


# ---------------------------------------------------------------- Fig. 3(b)

def pick_fig3(paths):
    """The five data curves inside the frame, in drawing order, checked against the panel's
    list. Points on the frame edge are the plot clipping the curve (n_sa below 10^20 m^-3 at
    the anode) and are dropped."""
    (xl, xr), (yt, yb) = FIG3["frame_x"], FIG3["frame_y"]
    found = [p for p in paths if in_frame(p, FIG3)]
    want = FIG3["curves"]
    assert len(found) == len(want), [(p["color"], p["dash"], len(p["xy"])) for p in found]
    out = {}
    for p, (name, col, dash) in zip(found, want):
        assert np.allclose(p["color"], col) and p["dash"] == dash, (name, p["color"], p["dash"])
        xy = p["xy"]
        xy = xy[(xy[:, 1] < yb - 0.01) & (xy[:, 1] > yt + 0.01)]
        xy = xy[np.argsort(xy[:, 0], kind="stable")]
        out[name] = thin(xy)
    return out


def extract_fig3(pdf):
    paths = trace_page(pdf, FIG3["page"])
    print(f"--- Semenov et al. 2016, Fig. 3(b): {CAPTION}")
    cal = calibrate(paths, FIG3)
    curves = pick_fig3(paths)
    data = {n: (cal["x"][0](xy[:, 0]), cal[FIG3_AXIS[n]][0](xy[:, 1]))
            for n, xy in curves.items()}

    xe, te = data["Te"]
    xh, th = data["Th"]
    print("Table 2 vs extraction:")
    print(f"  T_e at the anode   paper {TABLE2['Te_wall'] / 1e3:.2f} kK   extracted at "
          f"x = {xe[0]:.1e} m: {te[0] / 1e3:.2f} kK")
    print(f"  T_e undisturbed    paper {TABLE2['Te_bulk'] / 1e3:.2f} kK   extracted at "
          f"x = {xe[-1]:.1e} m: {te[-1] / 1e3:.2f} kK")
    print(f"  T_h at the anode   T_s = 1.00 kK         extracted: {th[0] / 1e3:.2f} kK")
    far = {k: np.interp(np.log(5e-3), np.log(v[0]), v[1]) for k, v in data.items()}
    print(f"  far field (x = 5 mm): n_i/n_e - 1 = {far['ni'] / far['ne'] - 1:+.3f},  "
          f"n_sa/n_e - 1 = {far['nsa'] / far['ne'] - 1:+.3f},  T_h/T_e - 1 = "
          f"{far['Th'] / far['Te'] - 1:+.3f}")

    write_csv(os.path.join(HERE, "fig3b_extracted.csv"),
              [f"Semenov, Krivtsun & Reisgen 2016, J. Phys. D 49 105204, Fig. 3(b): {CAPTION}",
               "extracted from the PDF vector paths (extract.py); x = distance from anode"],
              FIG3_UNITS, data)
    return data, cal


# ---------------------------------------------------------------- Figs. 4 and 5

def pick_red(paths, panel):
    """The red (1500 A/cm^2) polyline inside the frame, thinned, in x order."""
    found = [p for p in paths if in_frame(p, panel) and np.allclose(p["color"], RED)]
    assert len(found) == 1, [(p["color"], len(p["xy"])) for p in found]
    xy = found[0]["xy"]
    return thin(xy[np.argsort(xy[:, 0], kind="stable")])


def extract_fig45(pdf):
    paths = trace_page(pdf, 10)
    data, cals = {}, {}
    for name, panel in FIG45.items():
        print(f"--- Semenov et al. 2016, {panel['fig']}: {CAPTION}")
        cal = calibrate(paths, panel)
        xy = pick_red(paths, panel)
        data[name] = (cal["x"][0](xy[:, 0]), cal["y"][0](xy[:, 1]))
        cals[name] = cal

    x, j = data["ji_j"]
    i = np.argmin(j)
    print("checks:")
    print(f"  j_i/j   at the anode (x = {x[0]:.1e} m): {j[0]:.5f}   minimum {j[i]:.5f} at "
          f"x = {x[i]:.2e} m   at x = {x[-1]:.1e} m: {j[-1]:+.5f}")
    x, e = data["E"]
    k = np.nonzero(np.diff(np.sign(e)) > 0)[0]
    zero = f"{x[k[0]]:.2e} m" if len(k) else "none"
    print(f"  E       at the anode (x = {x[0]:.1e} m): {e[0] / 1e6:.3f}e6 V/m   sign change "
          f"(- to +) at x = {zero}   at x = {x[-1]:.1e} m: {e[-1]:+.3e} V/m")

    write_csv(os.path.join(HERE, "fig45_j1500_extracted.csv"),
              [f"Semenov, Krivtsun & Reisgen 2016, J. Phys. D 49 105204, Fig. 4 (ji_j) and "
               f"Fig. 5 (E), j = 1500 A/cm^2: {CAPTION}",
               "extracted from the PDF vector paths (extract.py); x = distance from "
               "anode; paper signs (j > 0, x away from the anode)"],
              {n: p["unit"] for n, p in FIG45.items()}, data)
    return data, cals


# ---------------------------------------------------------------- verification figures

def panel_raster(pdf, panel, dpi=600):
    """The panel's crop of the page as an image, and its extent in page points registered
    on the traced frame (so the inverse axis maps land on the raster)."""
    page = panel["page"]
    x0, x1, y0, y1 = panel["crop"]
    s = dpi / 72.0
    with tempfile.TemporaryDirectory() as td:
        stem = os.path.join(td, "crop")
        subprocess.run(["pdftoppm", "-f", str(page), "-l", str(page), "-r", str(dpi), "-png",
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
        return (i * w).sum() / w.sum() + 0.5            # pixel centre, image-edge origin

    px, py = panel["probe"]
    col = int((fxl + px * (fxr - fxl) - ext[0]) * s)
    row = int((fyt + py * (fyb - fyt) - ext[3]) * s)
    y_found = ext[3] + centroid(ink[:, col], int((fyb - ext[3]) * s)) / s
    x_found = ext[0] + centroid(ink[row, :], int((fxl - ext[0]) * s)) / s
    dx, dy = fxl - x_found, fyb - y_found
    print(f"raster registration: dx = {dx:+.2f} pt, dy = {dy:+.2f} pt "
          f"(page-box offset; the extracted data are unaffected)")
    return img, (ext[0] + dx, ext[1] + dx, ext[2] + dy, ext[3] + dy)


def _rc():
    plt.rcParams.update({"font.family": "serif", "font.serif": ["STIX", "STIXGeneral",
                         "DejaVu Serif"], "mathtext.fontset": "stix", "font.size": 8})


def _save(f, axs, ext_of, stem):
    for ax, ext in zip(axs.flat, ext_of):
        ax.set_xlim(ext[0], ext[1])
        ax.set_ylim(ext[2], ext[3])
        ax.set_axis_off()
    f.tight_layout()
    for e in ("pdf", "png"):
        f.savefig(os.path.join(HERE, f"{stem}.{e}"), dpi=300, bbox_inches="tight")
    plt.close(f)
    print(f"wrote {stem}.pdf/.png")


def verify_fig3(pdf, data, cal):
    """(a) the original panel; (b) the same raster, greyed, with every extracted curve drawn
    back over it through the inverse of the fitted axis maps."""
    img, ext = panel_raster(pdf, FIG3)
    style = {"ne": "#d62728", "ni": "#1f3fff", "nsa": "#2ca02c", "Te": "#ff7f0e",
             "Th": "#17becf"}
    _rc()
    f, axs = plt.subplots(1, 2, figsize=(7.0, 3.1))
    axs[0].imshow(img, extent=ext)
    axs[0].set_title("(a) Semenov et al. (2016), Fig. 3(b) — original", fontsize=8.5)
    axs[1].imshow(0.55 + 0.45 * img[..., :3].mean(axis=2), cmap="gray", vmin=0, vmax=1,
                  extent=ext)
    for name, (xs, vs) in data.items():
        axs[1].plot(cal["x"][1](xs), cal[FIG3_AXIS[name]][1](vs), lw=0.9, color=style[name],
                    label={"nsa": "$n_{sa}$"}.get(name, name))
    axs[1].set_title("(b) extracted data re-plotted over the original (greyed)", fontsize=8.5)
    axs[1].legend(loc="lower left", bbox_to_anchor=(1.0, 0.0), fontsize=7, frameon=False,
                  handlelength=1.6)
    _save(f, axs, (ext, ext), "fig3b_verification")


def verify_fig45(pdf, data, cals):
    """Per figure: the original panel, and the same raster greyed with the extracted curve
    drawn back over it."""
    _rc()
    f, axs = plt.subplots(2, 2, figsize=(7.0, 6.0))
    exts = []
    for row, (name, panel) in zip(axs, FIG45.items()):
        img, ext = panel_raster(pdf, panel)
        row[0].imshow(img, extent=ext)
        row[0].set_title(f"Semenov et al. (2016), {panel['fig']} — original", fontsize=8.5)
        row[1].imshow(0.55 + 0.45 * img[..., :3].mean(axis=2), cmap="gray", vmin=0, vmax=1,
                      extent=ext)
        cal = cals[name]
        row[1].plot(cal["x"][1](data[name][0]), cal["y"][1](data[name][1]), lw=0.9,
                    color="#d62728")
        row[1].set_title("1500 A cm$^{-2}$ extracted, over the original (greyed)",
                         fontsize=8.5)
        exts += [ext, ext]
    _save(f, axs, exts, "fig45_j1500_verification")


def main():
    ap = argparse.ArgumentParser(description="Extract Semenov et al. (2016) Fig. 3(b) and the "
                                 "1500 A/cm^2 curves of Figs. 4 and 5 from the PDF into CSVs "
                                 "beside this script.")
    ap.add_argument("--pdf", default=None, help=f"path to '{PDF_NAME}' (default: References/)")
    ap.add_argument("--no-verify", action="store_true", help="skip the verification figures")
    a = ap.parse_args()
    pdf = os.path.abspath(a.pdf) if a.pdf else find_pdf()
    d3, c3 = extract_fig3(pdf)
    d45, c45 = extract_fig45(pdf)
    if not a.no_verify:
        verify_fig3(pdf, d3, c3)
        verify_fig45(pdf, d45, c45)


if __name__ == "__main__":
    main()
