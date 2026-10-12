"""
Extract the jc = 1e7 A/m^2 curves of Almeida, Benilov & Naidis (2008) J. Phys. D 41 245201
(Ar, p0 = 1 bar, Tw = 3500 K) straight from the PDF's VECTOR paths, both on page 9:

    Fig. 5(b)        n_e, n_i, n_S, T_e, T_h, T_eq, j_i/j      -> fig5b_extracted.csv
    Fig. 6 curve b   E (solid), phi (dashed), and the two markers on phi
                                                                -> fig6b_extracted.csv

    python extract.py [--pdf paper.pdf]

Also writes fig5b_verification.{pdf,png} and fig6b_verification.{pdf,png}: the original panel,
and the same raster greyed with every extracted curve drawn back over it.

No pixel digitising: `mutool trace` dumps every stroked path on the page in PDF points, the
curves are picked out by colour/line style, and each axis is calibrated by least squares on
the figure's OWN tick marks (every major and minor tick). The tick residuals are printed, and
the extracted data are checked against the paper:
    Fig. 5(b)  Table 1 (wall j_i/j, T_e max), quasi-neutrality away from the sheath
    Fig. 6 b   the phi-max square against the curve's maximum, phi against -int E dx, and the
               1 % charge-separation circle against the Fig. 5(b) densities

Fig. 6 has three y axes: E of curves b and c on the LEFT in 10^8 V/m, E of curve a on the inner
RIGHT in 10^6 V/m, phi on the outer RIGHT in V. The three solid (blue) and the three dashed
(red) curves are drawn in the order a, b, c; each is chained from its pieces (a solid curve is
two strokes, a dashed one ~40 dash strokes) by end-to-start continuity.

x is the distance from the cathode surface (cathode on the LEFT). E is the paper's, i.e. along
x (away from the cathode): negative in the sheath. phi is relative to the cathode surface.

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

PDF_NAME = "Almeida_2008_J._Phys._D__Appl._Phys._41_245201.pdf"
PAGE = 9                                    # printed page 8

BLACK, RED, GREEN, ORANGE, BLUE = (0, 0, 0), (1, 0, 0), (0, .604, .2), (1, .396, .2), (0, 0, 1)


def _decades(lo, hi):
    """Tick values of a log axis from 10^lo to 10^hi: 1..9 x 10^k, then 10^hi."""
    return np.array([m * 10.0 ** k for k in range(lo, hi) for m in range(1, 10)] + [10.0 ** hi])


# Each panel: its frame in page points (y measured from the TOP of the page), the tick values
# of each axis in the order the ticks are met going away from the axis origin (bottom-up for
# y axes, left-right for x).
FIG5 = dict(
    caption="Ar, p0 = 1 bar, Tw = 3500 K, jc = 1e7 A/m^2",
    frame_x=(337.9, 474.7), frame_y=(80.3, 217.2),
    axes={
        "x": dict(side="bottom",    log=True,  values=_decades(-7, -2)),
        "n": dict(side="left",      log=True,  values=_decades(19, 24)),
        # j_i/j: majors 0 ... 0.3 (labels 0, 0.2, 0.3), minors every 1/30, one below 0
        "J": dict(side="right_in",  log=False, values=(np.arange(11) - 1) / 30.0),
        # T: log, 3 ... 10 and 20 ... 50 kK (the long tick is 10 kK)
        "T": dict(side="right_out", log=True,
                  values=np.array([3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50]) * 1e3),
    },
    # the curves in the order the PDF draws them with the colour each must have; dotted black
    # curves are told apart by that order
    curves=[  # name, rgb, solid?
        ("ne", BLACK, True), ("ni", RED, False), ("nS", BLACK, False),
        ("Te", GREEN, False), ("Th", ORANGE, False), ("Teq", BLACK, False),
        ("ji_j", BLUE, False),
    ],
    table1={"ji_w": 0.216, "Te_max": 18.5e3},
    crop=(300.0, 505.0, 70.0, 240.0),
)
FIG5_UNITS = {"ne": "m^-3", "ni": "m^-3", "nS": "m^-3", "Te": "K", "Th": "K", "Teq": "K",
              "ji_j": "-"}
FIG5_AXIS = {"ne": "n", "ni": "n", "nS": "n", "Te": "T", "Th": "T", "Teq": "T", "ji_j": "J"}

FIG6 = dict(
    caption="Ar, p0 = 1 bar, Tw = 3500 K, jc = 1e7 A/m^2 (curve b)",
    frame_x=(355.51, 507.35), frame_y=(453.87, 605.82),
    tick_lw_max=0.7,                         # this figure's frame and ticks are 0.65 pt
    axes={
        "x":   dict(side="bottom",    log=True,  values=_decades(-7, -2)),
        # E of b, c: labels -1.2, -0.8, -0.4, 0 (10^8 V/m), minors halfway
        "E8":  dict(side="left",      log=False, values=np.arange(-6, 1) * 0.2e8),
        # E of a: labels 0, 1, 2 (10^6 V/m), minors halfway
        "E6":  dict(side="right_in",  log=False, values=np.arange(5) * 0.5e6),
        # phi: labels 0, 4, ..., 16 V, minors halfway
        "phi": dict(side="right_out", log=False, values=np.arange(0, 17, 2.0)),
    },
    crop=(322.0, 528.0, 445.0, 625.0),
)
FIG6_INDEX = 1                               # curve b in the drawing order a, b, c
FIG6_UNITS = {"E": "V/m", "phi": "V", "phi_1pct": "V", "phi_max": "V"}


# ---------------------------------------------------------------- PDF tracing and axes

def find_pdf():
    # the repository's References/ (any level up)
    for k in range(5):
        p = os.path.join(HERE, *[".."] * k, "References", PDF_NAME)
        if os.path.isfile(p):
            return os.path.abspath(p)
    raise SystemExit(f"cannot find {PDF_NAME} in References/ -- pass --pdf")


def trace_page(pdf):
    with tempfile.NamedTemporaryFile(suffix=".xml") as tmp:
        subprocess.run(["mutool", "trace", pdf, str(PAGE)], check=True, stdout=tmp,
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
    """Page positions of every tick on one side of the frame, deduplicated, sorted away
    from the axis origin (bottom-up for y, left-right for x)."""
    (xl, xr), (yt, yb) = panel["frame_x"], panel["frame_y"]
    segs = [p["xy"] for p in paths if len(p["xy"]) == 2 and p["color"] == BLACK
            and p["lw"] < panel.get("tick_lw_max", 0.5)]
    pos = set()
    for s in segs:
        (x0, y0), (x1, y1) = s
        length = np.hypot(x1 - x0, y1 - y0)
        if not 1.0 < length < 6.0:
            continue
        if side == "bottom" and near(x0, x1, .05) and near(min(y0, y1), yb, .1) \
                and xl - 1 < x0 < xr + 1:
            pos.add(round(x0, 2))
        elif near(y0, y1, .05) and yt - 1 < y0 < yb + 1:
            lo, hi = min(x0, x1), max(x0, x1)
            if side == "left" and near(hi, xl, .1) and lo < xl - 1:
                pos.add(round(y0, 2))
            elif side == "right_in" and near(hi, xr, .1) and lo < xr - 1:
                pos.add(round(y0, 2))
            elif side == "right_out" and near(lo, xr, .1) and hi > xr + 1:
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
        print(f"  {name:3s} ({'log' if ax['log'] else 'lin'})  {len(pos):3d} ticks   "
              f"max |res| = {np.abs(res).max():.3f} pt")
    return cal


def write_csv(path, header, units, data):
    with open(path, "w") as f:
        for h in header:
            f.write(f"# {h}\n")
        f.write("quantity,unit,x_m,value\n")
        for name, (xs, vs) in data.items():
            for a, b in zip(xs, vs):
                f.write(f"{name},{units[name]},{a:.6e},{b:.6e}\n")
    print(f"wrote {path}  ({sum(len(v[0]) for v in data.values())} points, {len(data)} curves)")


# ---------------------------------------------------------------- Fig. 5(b)

def curve_runs(paths, panel):
    """Consecutive runs of same-style strokes inside the frame: [(rgb, solid, [xy...])]."""
    (xl, xr), (yt, yb) = panel["frame_x"], panel["frame_y"]
    runs = []
    for p in paths:
        xy = p["xy"]
        if not (xl - 0.2 <= xy[:, 0].min() and xy[:, 0].max() <= xr + 0.2
                and yt - 6 < xy[:, 1].min() and xy[:, 1].max() < yb + 0.2):
            continue
        if p["lw"] > 0.9:                    # the n_em arrow on the left axis
            continue
        if len(xy) == 2 and p["color"] == BLACK and \
                (near(xy[0, 0], xy[1, 0], .05) or near(xy[0, 1], xy[1, 1], .05)):
            continue                         # frame / tick marks
        key = (tuple(float(c) for c in p["color"]), len(xy) > 100)
        if runs and runs[-1][:2] == key:
            runs[-1][2].append(xy)
        else:
            runs.append((*key, [xy]))
    return runs


def pick_fig5(paths, panel):
    runs = curve_runs(paths, panel)
    want = panel["curves"]
    assert len(runs) == len(want), [(r[0], r[1], len(r[2])) for r in runs]
    out = {}
    for (rgb, solid, segs), (name, rgb_w, solid_w) in zip(runs, want):
        assert np.allclose(rgb, rgb_w) and solid == solid_w, (name, rgb, solid)
        xy = np.concatenate(segs)
        out[name] = xy[np.argsort(xy[:, 0])]
    return out


def extract_fig5(paths):
    panel = FIG5
    print(f"--- Almeida et al. 2008, Fig. 5(b): {panel['caption']}")
    cal = calibrate(paths, panel)
    curves = pick_fig5(paths, panel)
    data = {n: (cal["x"][0](xy[:, 0]), cal[FIG5_AXIS[n]][0](xy[:, 1]))
            for n, xy in curves.items()}

    t1 = panel["table1"]
    x, j = data["ji_j"]
    print("Table 1 vs extraction:")
    print(f"  j_iw/j_c   paper {t1['ji_w']:.4g}   extracted at x = {x[0]:.2e} m: {j[0]:.4g}")
    xt, te = data["Te"]
    print(f"  Te(max)    paper {t1['Te_max'] / 1e3:.1f} kK  extracted max: "
          f"{te.max() / 1e3:.2f} kK (at x = {xt[np.argmax(te)]:.2e} m)")
    x_ne, ne = data["ne"]
    x_ni, ni = data["ni"]
    x0 = 3e-6
    far = (x_ne > x0) & (x_ne < 5e-3)
    r = np.interp(np.log(x_ne[far]), np.log(x_ni), ni) / ne[far] - 1
    print(f"  quasi-neutrality, {x0 * 1e6:.0f} um < x < 5 mm: max |n_i/n_e - 1| = "
          f"{np.abs(r).max():.3f}")

    write_csv(os.path.join(HERE, "fig5b_extracted.csv"),
              [f"Almeida, Benilov & Naidis 2008, J. Phys. D 41 245201, Fig. 5(b): "
               f"{panel['caption']}",
               "extracted from the PDF vector paths (extract.py); x = distance from cathode"],
              FIG5_UNITS, data)
    return data, cal


# ---------------------------------------------------------------- Fig. 6 curve b

def chains(paths, rgb):
    """Curves of one colour, in drawing order: each stroke continues the previous one if it
    starts within 4 pt of where that one ended, else it opens a new chain. Only chains that
    start on the left frame are curves (the arrows pointing to the axes are not)."""
    (xl, xr), (yt, yb) = FIG6["frame_x"], FIG6["frame_y"]
    out = []
    for p in paths:
        xy = p["xy"]
        if p["color"] != rgb or p["lw"] < 0.5 or not (
                xl - 0.2 <= xy[:, 0].min() and xy[:, 0].max() <= xr + 0.2
                and yt - 0.2 <= xy[:, 1].min() and xy[:, 1].max() <= yb + 0.2):
            continue
        if out and np.hypot(*(xy[0] - out[-1][-1][-1])) < 4.0:
            out[-1].append(xy)
        else:
            out.append([xy])
    out = [np.concatenate(c) for c in out if near(c[0][0, 0], xl, 0.2)]
    return [c[np.argsort(c[:, 0], kind="stable")] for c in out]


def markers(paths, curve):
    """Centres of the marker outlines (thin black closed strokes) lying on `curve`:
    {'circle': (x, y), 'square': (x, y)} in page points."""
    found = {}
    for p in paths:
        xy = p["xy"]
        if p["color"] != BLACK or p["lw"] > 0.2 or len(xy) < 5:
            continue
        c = 0.5 * (xy.min(axis=0) + xy.max(axis=0))
        if abs(np.interp(c[0], curve[:, 0], curve[:, 1]) - c[1]) < 1.0:
            found["square" if len(xy) == 5 else "circle"] = tuple(c)
    assert set(found) == {"circle", "square"}, found
    return found


def extract_fig6(paths, fig5):
    print(f"--- Almeida et al. 2008, Fig. 6: {FIG6['caption']}")
    cal = calibrate(paths, FIG6)
    solid, dashed = chains(paths, BLUE), chains(paths, RED)
    assert len(solid) == 3 and len(dashed) == 3, (len(solid), len(dashed))
    e_xy, phi_xy = solid[FIG6_INDEX], dashed[FIG6_INDEX]
    mk = markers(paths, phi_xy)

    X, E, PHI = cal["x"][0], cal["E8"][0], cal["phi"][0]
    data = {"E": (X(e_xy[:, 0]), E(e_xy[:, 1])), "phi": (X(phi_xy[:, 0]), PHI(phi_xy[:, 1]))}
    for q, m in (("phi_1pct", "circle"), ("phi_max", "square")):
        data[q] = (X(np.array([mk[m][0]])), PHI(np.array([mk[m][1]])))

    xe, e = data["E"]
    xp, phi = data["phi"]
    i = np.argmax(phi)
    (xs,), (ps,) = data["phi_max"]
    print("checks:")
    print(f"  phi max     square at x = {xs:.3e} m, {ps:.2f} V;  curve max at "
          f"x = {xp[i]:.3e} m, {phi[i]:.2f} V")
    # phi(x) - phi(x0) = -int_x0^x E dx, over the sheath where E is resolved on 1e8 V/m
    print("  phi(x) - phi(1e-7 m) vs -int E dx (trapezoid on the extracted E):")
    for xq in (3e-7, 1e-6, 3e-6, 1e-5):
        m = xe <= xq
        xx, ee = np.r_[xe[m], xq], np.r_[e[m], np.interp(xq, xe, e)]
        integ = -np.sum(0.5 * (ee[1:] + ee[:-1]) * np.diff(xx))
        dphi = np.interp(xq, xp, phi) - phi[0]
        print(f"    x = {xq:.0e} m:  dphi = {dphi:7.2f} V   -int E dx = {integ:7.2f} V")
    # the circle: |n_i - n_e| / n_i = 1 %, from the Fig. 5(b) densities
    x5 = np.geomspace(1e-7, 1e-4, 4000)
    ne = np.exp(np.interp(np.log(x5), np.log(fig5["ne"][0]), np.log(fig5["ne"][1])))
    ni = np.exp(np.interp(np.log(x5), np.log(fig5["ni"][0]), np.log(fig5["ni"][1])))
    x1 = x5[np.argmax(np.abs(ni - ne) / ni < 0.01)]
    (xc,), (pc,) = data["phi_1pct"]
    print(f"  1 % charge separation   circle at x = {xc:.3e} m ({pc:.2f} V);  "
          f"Fig. 5(b) densities: x = {x1:.3e} m")

    write_csv(os.path.join(HERE, "fig6b_extracted.csv"),
              [f"Almeida, Benilov & Naidis 2008, J. Phys. D 41 245201, Fig. 6, curve b: "
               f"{FIG6['caption']}",
               "extracted from the PDF vector paths (extract.py); x = distance from "
               "cathode; E along x (paper sign); phi relative to the cathode, in V",
               "phi_1pct: circle, charge separation reaches 1%; phi_max: square, phi maximum"],
              FIG6_UNITS, data)
    return data, cal


# ---------------------------------------------------------------- verification figures

def panel_raster(pdf, panel, dpi=600):
    """The panel's crop of the page as an image, and its extent in page points registered
    on the traced frame (so the inverse axis maps land on the raster)."""
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

    # Register the raster on the frame itself. pdftoppm places the page by its MediaBox
    # (842 pt tall) while mutool trace uses the CropBox (841 pt), so the raster sits ~1 pt
    # low if the two are assumed to agree. Find the bottom and left frame lines in the image
    # (ink-weighted centroid) and shift the raster so they land on the traced frame.
    ink = 1.0 - img[..., :3].mean(axis=2)
    (fxl, fxr), (fyt, fyb) = panel["frame_x"], panel["frame_y"]

    def centroid(profile, guess, half=15):
        i = np.arange(max(guess - half, 0), guess + half)
        w = profile[i]
        return (i * w).sum() / w.sum() + 0.5            # pixel centre, image-edge origin

    # a column and a row through the plot, clear of tick marks and curve ends at the frame
    col = int((fxl + 0.5 * (fxr - fxl) - ext[0]) * s)
    row = int((fyt + 0.5 * (fyb - fyt) - ext[3]) * s)
    y_found = ext[3] + centroid(ink[:, col], int((fyb - ext[3]) * s)) / s
    x_found = ext[0] + centroid(ink[row, :], int((fxl - ext[0]) * s)) / s
    dx, dy = fxl - x_found, fyb - y_found
    print(f"raster registration: dx = {dx:+.2f} pt, dy = {dy:+.2f} pt "
          f"(MediaBox vs CropBox; the extracted data are unaffected)")
    return img, (ext[0] + dx, ext[1] + dx, ext[2] + dy, ext[3] + dy)


def verification_figure(pdf, panel, title, lines, stem, what="data"):
    """(a) the original panel, rasterised from the PDF; (b) the same raster, greyed, with
    `lines` [(page_x, page_y, plot kwargs)] drawn over it."""
    img, ext = panel_raster(pdf, panel)
    plt.rcParams.update({"font.family": "serif", "font.serif": ["STIX", "STIXGeneral",
                         "DejaVu Serif"], "mathtext.fontset": "stix", "font.size": 8})
    f, axs = plt.subplots(1, 2, figsize=(7.0, 3.1))
    axs[0].imshow(img, extent=ext)
    axs[0].set_title(f"(a) Almeida et al. (2008), {title} — original", fontsize=8.5)
    axs[1].imshow(0.55 + 0.45 * img[..., :3].mean(axis=2), cmap="gray", vmin=0, vmax=1,
                  extent=ext)
    for px, py, kw in lines:
        axs[1].plot(px, py, **kw)
    axs[1].set_title(f"(b) extracted {what} re-plotted over the original (greyed)", fontsize=8.5)
    axs[1].legend(loc="lower left", bbox_to_anchor=(1.0, 0.0), fontsize=7, frameon=False,
                  handlelength=1.6)
    for ax in axs:
        ax.set_xlim(ext[0], ext[1])
        ax.set_ylim(ext[2], ext[3])
        ax.set_axis_off()
    f.tight_layout()
    for e in ("pdf", "png"):
        f.savefig(os.path.join(HERE, f"{stem}.{e}"), dpi=300, bbox_inches="tight")
    plt.close(f)
    print(f"wrote {stem}.pdf/.png")


def verify_fig5(pdf, data, cal):
    style = {"ne": "#000000", "ni": "#d62728", "nS": "#7f7f7f", "Te": "#2ca02c",
             "Th": "#ff7f0e", "Teq": "#7f7f7f", "ji_j": "#1f3fff"}
    lines = [(cal["x"][1](xs), cal[FIG5_AXIS[n]][1](vs),
              dict(lw=0.9, color=style[n], label=n.replace("ji_j", "$j_i/j$")))
             for n, (xs, vs) in data.items()]
    verification_figure(pdf, FIG5, "Fig. 5(b)", lines, "fig5b_verification")


def verify_fig6(pdf, data, cal):
    X, E, PHI = cal["x"][1], cal["E8"][1], cal["phi"][1]
    lines = [(X(data["E"][0]), E(data["E"][1]), dict(lw=0.9, color="#1f3fff", label="E (b)")),
             (X(data["phi"][0]), PHI(data["phi"][1]),
              dict(lw=0.9, color="#d62728", label=r"$\varphi$ (b)"))]
    for q, mk, lab in (("phi_1pct", "o", "1 % sep."), ("phi_max", "s", r"$\varphi$ max")):
        lines.append((X(data[q][0]), PHI(data[q][1]),
                      dict(marker=mk, ls="none", ms=5, mfc="none", mec="#2ca02c", mew=0.9,
                           label=lab)))
    verification_figure(pdf, FIG6, "Fig. 6", lines, "fig6b_verification", "curve b")


def main():
    ap = argparse.ArgumentParser(description="Extract Almeida et al. (2008) Fig. 5(b) and "
                                 "Fig. 6 curve b from the PDF into CSVs beside this script.")
    ap.add_argument("--pdf", default=None, help=f"path to {PDF_NAME} (default: References/)")
    ap.add_argument("--no-verify", action="store_true", help="skip the verification figures")
    a = ap.parse_args()
    pdf = os.path.abspath(a.pdf) if a.pdf else find_pdf()
    paths = trace_page(pdf)
    d5, c5 = extract_fig5(paths)
    d6, c6 = extract_fig6(paths, d5)
    if not a.no_verify:
        verify_fig5(pdf, d5, c5)
        verify_fig6(pdf, d6, c6)


if __name__ == "__main__":
    main()
