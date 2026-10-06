"""
Extract the curves of Almeida, Benilov & Naidis (2008) J. Phys. D 41 245201, figure 5(a), 5(b)
or 5(c) (Ar, p0 = 1 bar, Tw = 3500 K; jc = 1e6 / 1e7 / 7.8e7 A/m^2) straight from the PDF's
VECTOR paths.

No pixel digitising: `mutool trace` dumps every stroked path on page 9 in PDF points, the
curves are picked out by colour/line style, and each axis is calibrated by least squares on
the figure's OWN tick marks (every major and minor tick). The tick residuals are printed, and
the extracted data are checked against the numbers in the paper's Table 1.

    python fig5_extract.py a        -> fig5a_extracted.csv, fig5a_verification.{pdf,png}
    python fig5_extract.py b        -> fig5b_extracted.csv, fig5b_verification.{pdf,png}

compare5a.py / compareJ1e7.py call extract() themselves when the CSV is missing. Fig. 5(c) is
extracted into ../almeida_J7.8e7/ (its extract.py passes outdir=), as any panel can be.

x in the paper is the distance from the cathode surface (cathode on the LEFT).
"""
import glob
import os
import subprocess
import sys
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pdftrace import parse

HERE = os.path.dirname(os.path.abspath(__file__))
PDF_NAME = "Almeida_2008_J._Phys._D__Appl._Phys._41_245201.pdf"
PAGE = 9                                    # printed page 8

UNITS = {"ne": "m^-3", "ni": "m^-3", "nS": "m^-3", "nb": "m^-3",
         "Te": "K", "Th": "K", "Teq": "K", "ji_j": "-"}
AXIS_OF = {"ne": "n", "ni": "n", "nS": "n", "nb": "n", "Te": "T", "Th": "T", "Teq": "T",
           "ji_j": "J"}

BLACK, RED, GREEN, ORANGE, BLUE = (0, 0, 0), (1, 0, 0), (0, .604, .2), (1, .396, .2), (0, 0, 1)


def _decades(lo, hi):
    """Tick values of a log axis from 10^lo to 10^hi: 1..9 x 10^k, then 10^hi."""
    return np.array([m * 10.0 ** k for k in range(lo, hi) for m in range(1, 10)] + [10.0 ** hi])


# Every panel: its frame in page points (y measured from the TOP of the page), the tick values
# of each axis in the order the ticks are met going away from the axis origin (bottom-up for
# y axes, left-right for x), and the curves in the order the PDF draws them with the colour
# each must have. Dotted black curves are told apart by that order; their labels in the PDF
# sit at the start of each run.
PANELS = {
    "a": dict(
        caption="Ar, p0 = 1 bar, Tw = 3500 K, jc = 1e6 A/m^2",
        frame_x=(119.2, 261.8), frame_y=(74.7, 217.2),
        axes={
            "x": dict(side="bottom",    log=True,  values=_decades(-7, -2)),
            "n": dict(side="left",      log=True,  values=_decades(18, 22)),
            # j_i/j: majors -2e-4 ... 6e-4 (labels "0", "0.0006"), minors halfway
            "J": dict(side="right_in",  log=False, values=np.arange(-2, 7) * 1e-4),
            # T: labels 4, 8, 12 (10^3 K), majors every 2 kK, no minors
            "T": dict(side="right_out", log=False, values=np.arange(4, 13, 2) * 1e3),
        },
        # the two dotted density curves start together at the top; n_b is the one that
        # ends ON n_e at x = 1e-2 m (the paper: "n_e = n_b ... imposed by the boundary
        # conditions at the upper boundary"), n_S the one that stays above it
        curves=[  # name, rgb, solid?
            ("ne", BLACK, True), ("ni", RED, False), ("nS", BLACK, False),
            ("Te", GREEN, False), ("Th", ORANGE, False), ("Teq", BLACK, False),
            ("ji_j", BLUE, False), ("nb", BLACK, False),
        ],
        table1={"ji_w": 4e-4},
        crop=(82.0, 292.0, 67.0, 240.0),
    ),
    "b": dict(
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
        curves=[  # name, rgb, solid?
            ("ne", BLACK, True), ("ni", RED, False), ("nS", BLACK, False),
            ("Te", GREEN, False), ("Th", ORANGE, False), ("Teq", BLACK, False),
            ("ji_j", BLUE, False),
        ],
        table1={"ji_w": 0.216, "Te_max": 18.5e3},
        crop=(300.0, 505.0, 70.0, 240.0),
    ),
    "c": dict(
        caption="Ar, p0 = 1 bar, Tw = 3500 K, jc = 7.8e7 A/m^2",
        frame_x=(226.09, 363.31), frame_y=(253.03, 390.25),
        axes={
            "x": dict(side="bottom",    log=True,  values=_decades(-7, -2)),
            # labels 10^17 ... 10^23; the frame top is the 10^24 tick
            "n": dict(side="left",      log=True,  values=_decades(17, 24)),
            # j_i/j: labels 0, 0.2, ..., 0.8, ticks every 0.1 from -0.1 (the bottom corner)
            "J": dict(side="right_in",  log=False, values=(np.arange(10) - 1) / 10.0),
            # T: log, labels 10 and 100 (10^3 K); ticks 3 ... 9, 10 ... 90, 100 ... 600 kK
            "T": dict(side="right_out", log=True,
                      values=np.r_[np.arange(3, 10), np.arange(10, 100, 10),
                                   np.arange(100, 700, 100)] * 1e3),
        },
        # n_e is two solid strokes (one run); otherwise drawn in the same order as 5(b)
        curves=[  # name, rgb, solid?
            ("ne", BLACK, True), ("ni", RED, False), ("nS", BLACK, False),
            ("Te", GREEN, False), ("Th", ORANGE, False), ("Teq", BLACK, False),
            ("ji_j", BLUE, False),
        ],
        table1={"ji_w": 0.760, "Te_max": 430e3},
        qn_from=1e-5,                        # the sheath here reaches ~5 um
        crop=(185.0, 400.0, 245.0, 412.0),
    ),
}


def find_pdf(extra=()):
    # beside this script, in its parent, or in the repository's References/ (any level up)
    ups = [os.path.join(HERE, *[".."] * k) for k in range(4)]
    for d in (*extra, *ups[:2], *(os.path.join(u, "References") for u in ups)):
        p = os.path.join(d, PDF_NAME)
        if os.path.isfile(p):
            return os.path.abspath(p)
    raise SystemExit(f"cannot find {PDF_NAME} beside {HERE}, in its parent or in References/")


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
            print(f"  {name:2s} ({'log' if ax['log'] else 'lin'})  {len(pos):3d} ticks   "
                  f"max |res| = {np.abs(res).max():.3f} pt")
    return cal


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


def pick_curves(paths, panel):
    runs = curve_runs(paths, panel)
    want = panel["curves"]
    assert len(runs) == len(want), [(r[0], r[1], len(r[2])) for r in runs]
    out = {}
    for (rgb, solid, segs), (name, rgb_w, solid_w) in zip(runs, want):
        assert np.allclose(rgb, rgb_w) and solid == solid_w, (name, rgb, solid)
        xy = np.concatenate(segs)
        out[name] = xy[np.argsort(xy[:, 0])]
    return out


def trace_page(pdf):
    with tempfile.NamedTemporaryFile(suffix=".xml") as tmp:
        subprocess.run(["mutool", "trace", pdf, str(PAGE)], check=True, stdout=tmp,
                       stderr=subprocess.DEVNULL)
        tmp.flush()
        return parse(tmp.name)[0]


def csv_path(fig, outdir=HERE):
    return os.path.join(outdir, f"fig5{fig}_extracted.csv")


def extract(fig, pdf=None, verify_fig=True, outdir=HERE):
    """Extract panel `fig` ('a', 'b' or 'c'), run the checks, write the CSV (and the
    verification figure) into `outdir`. Returns {name: (x_m, value)}."""
    panel = PANELS[fig]
    pdf = pdf or find_pdf()
    paths = trace_page(pdf)
    print(f"--- Almeida et al. 2008, Fig. 5({fig}): {panel['caption']}")
    cal = calibrate(paths, panel)
    curves = pick_curves(paths, panel)
    data = {n: (cal["x"][0](xy[:, 0]), cal[AXIS_OF[n]][0](xy[:, 1])) for n, xy in curves.items()}
    checks(data, panel)

    path = csv_path(fig, outdir)
    with open(path, "w") as f:
        f.write(f"# Almeida, Benilov & Naidis 2008, J. Phys. D 41 245201, Fig. 5({fig}): "
                f"{panel['caption']}\n")
        f.write("# extracted from the PDF vector paths (fig5_extract.py); "
                "x = distance from cathode\n")
        f.write("quantity,unit,x_m,value\n")
        for name, (xs, vs) in data.items():
            for a, b in zip(xs, vs):
                f.write(f"{name},{UNITS[name]},{a:.6e},{b:.6e}\n")
    print(f"wrote {path}  ({sum(len(v[0]) for v in data.values())} points, {len(data)} curves)")
    if verify_fig:
        verification_figure(fig, pdf, data, cal, outdir)
    return data


def load(fig, pdf=None, outdir=HERE):
    """The extracted curves of Fig. 5(fig), extracting first if the CSV is not there."""
    path = csv_path(fig, outdir)
    if not os.path.isfile(path):
        return extract(fig, pdf, outdir=outdir)
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
    t1 = panel["table1"]
    x, j = data["ji_j"]
    print("Table 1 vs extraction:")
    print(f"  j_iw/j_c   paper {t1['ji_w']:.4g}   extracted at x = {x[0]:.2e} m: {j[0]:.4g}")
    if "Te_max" in t1:
        xt, te = data["Te"]
        print(f"  Te(max)    paper {t1['Te_max'] / 1e3:.1f} kK  extracted max: "
              f"{te.max() / 1e3:.2f} kK (at x = {xt[np.argmax(te)]:.2e} m)")
    x_ne, ne = data["ne"]
    x_ni, ni = data["ni"]
    x0 = panel.get("qn_from", 3e-6)
    far = (x_ne > x0) & (x_ne < 5e-3)
    r = np.interp(np.log(x_ne[far]), np.log(x_ni), ni) / ne[far] - 1
    print(f"  quasi-neutrality, {x0 * 1e6:.0f} um < x < 5 mm: max |n_i/n_e - 1| = "
          f"{np.abs(r).max():.3f}")


def panel_raster(pdf, panel, dpi=600, page=PAGE):
    """The panel's crop of the page as an image, and its extent in page points registered
    on the traced frame (so the inverse axis maps land on the raster)."""
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
    px, py = panel.get("probe", (0.5, 0.5))
    col = int((fxl + px * (fxr - fxl) - ext[0]) * s)
    row = int((fyt + py * (fyb - fyt) - ext[3]) * s)
    y_found = ext[3] + centroid(ink[:, col], int((fyb - ext[3]) * s)) / s
    x_found = ext[0] + centroid(ink[row, :], int((fxl - ext[0]) * s)) / s
    dx, dy = fxl - x_found, fyb - y_found
    print(f"raster registration: dx = {dx:+.2f} pt, dy = {dy:+.2f} pt "
          f"(MediaBox vs CropBox; the extracted data are unaffected)")
    return img, (ext[0] + dx, ext[1] + dx, ext[2] + dy, ext[3] + dy)


def verification_figure(fig, pdf, data, cal, outdir=HERE):
    """(a) the original panel, rasterised from the PDF; (b) the same raster, greyed, with
    every extracted curve drawn back over it through the inverse of the fitted axis maps."""
    panel = PANELS[fig]
    img, ext = panel_raster(pdf, panel)
    style = {"ne": "#000000", "ni": "#d62728", "nS": "#7f7f7f", "nb": "#9467bd",
             "Te": "#2ca02c", "Th": "#ff7f0e", "Teq": "#7f7f7f", "ji_j": "#1f3fff"}
    plt.rcParams.update({"font.family": "serif", "font.serif": ["STIX", "STIXGeneral",
                         "DejaVu Serif"], "mathtext.fontset": "stix", "font.size": 8})
    f, axs = plt.subplots(1, 2, figsize=(7.0, 3.1))
    axs[0].imshow(img, extent=ext)
    axs[0].set_title(f"(a) Almeida et al. (2008), Fig. 5({fig}) — original", fontsize=8.5)
    grey = img[..., :3].mean(axis=2)
    axs[1].imshow(0.55 + 0.45 * grey, cmap="gray", vmin=0, vmax=1, extent=ext)
    for name, (xs, vs) in data.items():
        axs[1].plot(cal["x"][1](xs), cal[AXIS_OF[name]][1](vs), lw=0.9, color=style[name],
                    label=name.replace("ji_j", "$j_i/j$"))
    axs[1].set_title("(b) extracted data re-plotted over the original (greyed)", fontsize=8.5)
    axs[1].legend(loc="lower left", bbox_to_anchor=(1.0, 0.0), fontsize=7, frameon=False,
                  handlelength=1.6)
    for ax in axs:
        ax.set_xlim(ext[0], ext[1])
        ax.set_ylim(ext[2], ext[3])
        ax.set_axis_off()
    f.tight_layout()
    for e in ("pdf", "png"):
        f.savefig(os.path.join(outdir, f"fig5{fig}_verification.{e}"), dpi=300,
                  bbox_inches="tight")
    plt.close(f)
    print(f"wrote fig5{fig}_verification.pdf/.png")


if __name__ == "__main__":
    figs = sys.argv[1:] or ["a", "b"]
    for fg in figs:
        extract(fg)
