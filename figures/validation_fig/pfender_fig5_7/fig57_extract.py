"""
Digitise Dinulescu & Pfender (1980) J. Appl. Phys. 51 3149, "Analysis of the anode boundary
layer of high intensity arcs", Figs. 5 and 7 (Ar, 100 kPa, j = 900 A/cm^2, anode T_w = 400 K):

    Fig. 5   T_e, T (heavy particles)      x in 1e-4 m, T in 1e4 K
    Fig. 7   n_e (diffusion eq.), n_e' (Saha at T_e)    x in 1e-3 m, n in 1e22 m^-3

Unlike Almeida (2008) and Semenov (2016), the PDF is a 300 dpi bitmap scan (one JBIG2 image
per page), so there are no vector paths to read. The curves are digitised from the raster
of page 8 (printed page 3155):

  * AXES. The scan is rotated by ~0.3-0.5 deg, so each frame is fitted as two straight lines
    (left and bottom edges, ink centroids along the edge) and a pixel is mapped to oblique
    coordinates (u along the bottom edge, v along the left edge) from their intersection.
    u and v are then calibrated by least squares on the tick marks (x ticks just below the
    bottom edge, y ticks just outside the left edge); residuals are printed.
  * CURVES. Inside the frame the two curves form ONE connected component (they merge at the
    arc end); every text label is a separate, smaller component and is dropped. Column by
    column, the ink runs of that component are read as: two runs -> upper and lower curve;
    one run no thicker than a line -> both curves (merged); one thicker run -> the two
    curves overlap partly, upper = top edge + w/2, lower = bottom edge - w/2 (w = line
    width). Where a curve is steep (its run in a column longer than 2.5 w) the column
    reading is replaced by a row-by-row reading of the same curve.
  * CHECKS against the paper's text (Sec. V): T = T_e = 12 400 K and dT_e/dx = 2.754e6 K/m
    at the arc end x = L, T_e = 8430 K and T = 400 K at the anode; the two figures must give
    the same L; n_e' against Saha evaluated on the digitised T_e (Sec. III: "n_e is computed
    from the Saha equation at the temperature T_e").

    python fig57_extract.py          -> fig5_extracted.csv, fig7_extracted.csv,
                                        fig5_verification.{pdf,png}, fig7_verification.{pdf,png}

fig57_compare.py calls load() and extracts automatically when a CSV is missing.
x in the paper is the distance from the anode surface (anode on the LEFT); it is written
here in metres.
"""
import os
import subprocess
import sys
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import ndimage as ndi

HERE = os.path.dirname(os.path.abspath(__file__))
PDF_NAME = "Analysis of anode BL.pdf"
PAGE = 8                                     # printed page 3155
DPI = 300                                    # the scan's own resolution

# The conditions of every curve in both figures, and the values quoted in the text (Sec. V)
J, J_ACM2, P0, T_WALL = 9e6, 900, 1.0e5, 400.0
TEXT = dict(T_arc=12400.0, dTdx_arc=2.754e6, Te_anode=8430.0, T_anode=400.0)

UNITS = {"Te": "K", "Th": "K", "ne": "m^-3", "nsa": "m^-3"}

# Pixel boxes (page raster, row/col) and axes of the two panels. `box` only bounds the search
# for the frame lines; the frame itself is measured. Tick values run away from the origin.
FIGS = {
    "5": dict(caption="Species temperatures, Ar, 100 kPa, j = 900 A/cm^2",
              frame_guess=dict(left=331, right=1098, bottom=2935, top=2332),
              x_ticks=np.arange(0, 11) * 0.8e-4,            # 0.0 ... 8.0 (x 1e-4 m) labelled
              y_ticks=np.arange(0, 10) * 0.2e4,             # 0.0 ... 1.8 (x 1e4 K)
              curves=("Te", "Th"),                          # upper, lower
              crop=(2230, 3010, 60, 1240)),
    "7": dict(caption="Charge carrier number density, Ar, 100 kPa, j = 900 A/cm^2",
              frame_guess=dict(left=1430, right=2113, bottom=2943, top=2432),
              x_ticks=np.arange(0, 12) * 0.08e-3,           # 0.00 ... 0.88 (x 1e-3 m)
              y_ticks=np.arange(0, 10) * 1.0e22,            # 0.0 ... 9.0 (x 1e22 m^-3)
              curves=("ne", "nsa"),                         # upper, lower
              crop=(2330, 3020, 1270, 2420)),
}


def find_pdf(extra=()):
    refs = os.path.join(HERE, "..", "..", "References")
    for d in (*extra, refs, HERE, os.path.join(HERE, "..")):
        p = os.path.join(d, PDF_NAME)
        if os.path.isfile(p):
            return os.path.abspath(p)
    raise SystemExit(f"cannot find '{PDF_NAME}' in References/, beside {HERE} or its parent")


def page_ink(pdf):
    """The page's scan (the single embedded bitmap, at its native 300 dpi) as a bool ink mask."""
    with tempfile.TemporaryDirectory() as td:
        subprocess.run(["pdfimages", "-f", str(PAGE), "-l", str(PAGE), "-png", pdf,
                        os.path.join(td, "p")], check=True)
        pngs = sorted(f for f in os.listdir(td) if f.endswith(".png"))
        img = max((plt.imread(os.path.join(td, f)) for f in pngs), key=lambda a: a.size)
    if img.ndim == 3:
        img = img[..., :3].mean(axis=2)
    return img, img < 0.5


def _runs(a):
    """(start, stop) of every True run of a 1-D bool array."""
    d = np.diff(np.r_[0, a.astype(np.int8), 0])
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)))


# ---------------------------------------------------------------- axes ----------------------
def frame(ink, g, half=6):
    """Left and bottom frame edges as fitted lines: col = a + b row, row = c + d col."""
    rows = np.arange(g["top"] + 20, g["bottom"] - 20, 10)
    cols = np.arange(g["left"] + 20, g["right"] - 20, 10)

    def centroid_line(ts, guess, along_rows):
        pts = []
        for t in ts:
            k = np.arange(guess - half, guess + half + 1)
            w = ink[t, k] if along_rows else ink[k, t]
            if 2 <= w.sum() <= 6:                    # the line alone (no curve, no tick)
                pts.append((t, (k * w).sum() / w.sum()))
        t, c = np.array(pts).T
        fitc = np.polyfit(t, c, 1)
        keep = np.abs(np.polyval(fitc, t) - c) < 1.0  # drop the odd blemish
        return np.polyfit(t[keep], c[keep], 1)[::-1]

    left = centroid_line(rows, g["left"], True)      # col = a + b*row
    bottom = centroid_line(cols, g["bottom"], False)  # row = c + d*col
    return left, bottom


def oblique(left, bottom):
    """Pixel (col, row) -> (u, v): coordinates along the bottom and left edges, in pixels,
    measured from their intersection (the plot origin)."""
    a, b = left
    c, d = bottom
    r0 = (c + d * a) / (1 - d * b)                   # intersection
    c0 = a + b * r0
    ex = np.array([1.0, d]) / np.hypot(1.0, d)       # (col, row) along the bottom edge
    ey = np.array([b, 1.0]) / np.hypot(b, 1.0)       # along the left edge, pointing DOWN
    M = np.linalg.inv(np.column_stack([ex, -ey]))    # (dc, dr) = u ex - v ey

    def to_uv(col, row):
        uv = M @ np.vstack([np.asarray(col, float) - c0, np.asarray(row, float) - r0])
        return uv[0], uv[1]
    return to_uv, (c0, r0), ex, ey


def ticks(ink, left, bottom, g, fig):
    """Tick centres: x ticks in the row 5 px below the bottom edge, y ticks in the column 5 px
    left of the left edge. Each is a 2-4 px wide ink run."""
    a, b = left
    c, d = bottom
    xt, yt = [], []
    for col in range(g["left"] - 10, g["right"] + 10):
        row = int(round(c + d * col)) + 5
        if ink[row, col]:
            xt.append((col, row))
    for row in range(g["top"] - 10, g["bottom"] + 10):
        col = int(round(a + b * row)) - 5
        if ink[row, col]:
            yt.append((col, row))

    def group(pts, key):
        pts = sorted(pts, key=key)
        out, cur = [], [pts[0]]
        for p in pts[1:]:
            if key(p) - key(cur[-1]) <= 1:
                cur.append(p)
            else:
                out.append(np.mean(cur, axis=0))
                cur = [p]
        out.append(np.mean(cur, axis=0))
        return np.array(out)

    xt = group(xt, key=lambda p: p[0])
    yt = group(yt, key=lambda p: -p[1])              # away from the origin (upwards)
    n_x, n_y = len(FIGS[fig]["x_ticks"]), len(FIGS[fig]["y_ticks"])
    assert len(xt) == n_x and len(yt) == n_y, \
        f"Fig. {fig}: {len(xt)} x ticks (want {n_x}), {len(yt)} y ticks (want {n_y})"
    return xt, yt


def calibrate(ink, fig, verbose=True):
    g = FIGS[fig]["frame_guess"]
    left, bottom = frame(ink, g)
    to_uv, origin, ex, ey = oblique(left, bottom)
    xt, yt = ticks(ink, left, bottom, g, fig)
    u, _ = to_uv(xt[:, 0], xt[:, 1])
    _, v = to_uv(yt[:, 0], yt[:, 1])
    X, Y = FIGS[fig]["x_ticks"], FIGS[fig]["y_ticks"]
    px, py = np.polyfit(u, X, 1), np.polyfit(v, Y, 1)
    res_x = (X - np.polyval(px, u)) / px[0]
    res_y = (Y - np.polyval(py, v)) / py[0]
    if verbose:
        rot = np.degrees(np.arctan(bottom[1])), np.degrees(np.arctan(-left[1]))
        print(f"frame: bottom edge tilted {rot[0]:+.2f} deg, left edge {rot[1]:+.2f} deg "
              f"(scan rotation/shear)")
        print(f"axis calibration: x {len(u)} ticks, max |res| = {np.abs(res_x).max():.2f} px; "
              f"y {len(v)} ticks, max |res| = {np.abs(res_y).max():.2f} px "
              f"(1 px = {px[0]:.3g} m, {py[0]:.3g} {'K' if fig == '5' else 'm^-3'})")

    def to_data(col, row):
        uu, vv = to_uv(col, row)
        return np.polyval(px, uu), np.polyval(py, vv)

    def to_pix(x, y):
        uu, vv = (x - px[1]) / px[0], (y - py[1]) / py[0]
        c0, r0 = origin
        return c0 + uu * ex[0] - vv * ey[0], r0 + uu * ex[1] - vv * ey[1]
    return dict(to_data=to_data, to_pix=to_pix, left=left, bottom=bottom, g=g,
                ticks=(xt, yt))


# ---------------------------------------------------------------- curves --------------------
def curve_pixels(ink, cal, inset=4):
    """The one large component inside the frame: both curves, no text."""
    g = cal["g"]
    a, b = cal["left"]
    c, d = cal["bottom"]
    rr, cc = np.mgrid[0:ink.shape[0], 0:ink.shape[1]]
    inside = ((cc > a + b * rr + inset) & (rr < c + d * cc - inset)
              & (cc < g["right"] - inset) & (rr > g["top"] + inset))
    lab, n = ndi.label(ink & inside, structure=np.ones((3, 3)))
    sizes = ndi.sum(np.ones_like(lab), lab, index=np.arange(1, n + 1))
    k = 1 + int(np.argmax(sizes))
    assert sizes[k - 1] > 5 * np.sort(sizes)[-2], "the curves are not one dominant component"
    return lab == k


def line_width(mask):
    """Median length of the single, thin column runs: the pen width in pixels."""
    w = []
    for col in range(mask.shape[1]):
        r = _runs(mask[:, col])
        if len(r) == 2:
            w += [e - s for s, e in r]
    return float(np.median(w))


def trace(mask, cal):
    """(upper, lower) curves as (col, row) pixel arrays, x increasing."""
    w = line_width(mask)
    cols = np.flatnonzero(mask.any(axis=0))
    runs = {col: _runs(mask[:, col]) for col in cols}
    two = np.array([col for col in cols if len(runs[col]) >= 2])
    up, lo = [], []
    steep_up, steep_lo = [], []
    for col in cols:
        r = runs[col]
        if len(r) >= 2:
            (s1, e1), (s2, e2) = r[0], r[-1]
            up.append((col, 0.5 * (s1 + e1 - 1)))
            lo.append((col, 0.5 * (s2 + e2 - 1)))
            if e1 - s1 > 2.5 * w:
                steep_up.append(col)
            if e2 - s2 > 2.5 * w:
                steep_lo.append(col)
            continue
        s, e = r[0]
        # One run is both curves only where they have merged: the nearest two-run column
        # must show them less than 2.5 w apart. Otherwise the other curve is simply absent
        # from this column (e.g. at the frame edge) and the run goes to the nearer curve.
        near = two[np.argmin(np.abs(two - col))]
        (s1, e1), (s2, e2) = runs[near][0], runs[near][-1]
        c_up, c_lo = 0.5 * (s1 + e1 - 1), 0.5 * (s2 + e2 - 1)
        if c_lo - c_up > 2.5 * w:
            mid = 0.5 * (s + e - 1)
            (up if abs(mid - c_up) < abs(mid - c_lo) else lo).append((col, mid))
        elif e - s <= w + 2:
            up.append((col, 0.5 * (s + e - 1)))
            lo.append((col, 0.5 * (s + e - 1)))
        else:
            up.append((col, s + 0.5 * (w - 1)))
            lo.append((col, e - 0.5 * (w + 1)))
    up, lo = np.array(up), np.array(lo)
    return _resteep(mask, up, steep_up), _resteep(mask, lo, steep_lo), w


def _resteep(mask, pts, steep_cols):
    """Replace the column reading over the steep stretch by a row reading: in each row of
    that stretch, the centre of the run nearest the curve's column reading."""
    if not steep_cols:
        return pts
    c_lo, c_hi = min(steep_cols), max(steep_cols)
    band = pts[(pts[:, 0] >= c_lo) & (pts[:, 0] <= c_hi)]
    r_lo, r_hi = int(band[:, 1].min()) - 3, int(band[:, 1].max()) + 3
    rows = []
    for row in range(r_lo, r_hi + 1):
        r = [(s, e) for s, e in _runs(mask[row, c_lo:c_hi + 1])]
        if not r:
            continue
        # the curve's own run in this row is the one nearest its column reading
        ref = np.interp(row, *(band[np.argsort(band[:, 1])][:, ::-1].T))
        s, e = min(r, key=lambda se: abs(c_lo + 0.5 * (se[0] + se[1] - 1) - ref))
        rows.append((c_lo + 0.5 * (s + e - 1), row))
    keep = pts[(pts[:, 0] < c_lo) | (pts[:, 0] > c_hi)]
    out = np.vstack([keep, np.array(rows)])
    return out[np.argsort(out[:, 0], kind="stable")]


# ---------------------------------------------------------------- driver --------------------
def csv_path(fig):
    return os.path.join(HERE, f"fig{fig}_extracted.csv")


def extract_one(fig, img, ink, verify_fig=True):
    spec = FIGS[fig]
    print(f"--- Dinulescu & Pfender 1980, Fig. {fig}: {spec['caption']}")
    cal = calibrate(ink, fig)
    mask = curve_pixels(ink, cal)
    up, lo, w = trace(mask, cal)
    print(f"curves: pen width {w:.1f} px, {len(up)} / {len(lo)} points")
    data = {}
    for name, pts in zip(spec["curves"], (up, lo)):
        x, y = cal["to_data"](pts[:, 0], pts[:, 1])
        x = np.clip(x, 0.0, None)
        o = np.argsort(x, kind="stable")
        data[name] = (x[o], y[o])

    path = csv_path(fig)
    with open(path, "w") as f:
        f.write(f"# Dinulescu & Pfender 1980, J. Appl. Phys. 51 3149, Fig. {fig}: "
                f"{spec['caption']}, T_w = 400 K\n")
        f.write("# digitised from the 300 dpi page scan (fig57_extract.py); "
                "x = distance from anode\n")
        f.write("quantity,unit,x_m,value\n")
        for name, (xs, vs) in data.items():
            for a, b in zip(xs, vs):
                f.write(f"{name},{UNITS[name]},{a:.6e},{b:.6e}\n")
    print(f"wrote {path}  ({sum(len(v[0]) for v in data.values())} points, {len(data)} curves)")
    if verify_fig:
        verification_figure(fig, img, data, cal, (up, lo))
    return data


def extract(pdf=None, verify_fig=True):
    """Both figures; returns {name: (x_m, value)} with Te, Th, ne, nsa."""
    pdf = pdf or find_pdf()
    img, ink = page_ink(pdf)
    data = {}
    for fig in FIGS:
        data.update(extract_one(fig, img, ink, verify_fig))
    checks(data)
    return data


def load(pdf=None):
    """The digitised curves of Figs. 5 and 7, extracting first if a CSV is not there."""
    if not all(os.path.isfile(csv_path(f)) for f in FIGS):
        return extract(pdf)
    cols = {}
    for fig in FIGS:
        with open(csv_path(fig)) as f:
            for line in f:
                if line.startswith(("#", "quantity")):
                    continue
                q, _, x, v = line.strip().split(",")
                cols.setdefault(q, ([], []))
                cols[q][0].append(float(x))
                cols[q][1].append(float(v))
    return {q: (np.array(x), np.array(v)) for q, (x, v) in cols.items()}


def n_saha(Te, Th, p=P0):
    """Saha (at T_e) + quasi-neutrality + Dalton, as ../semenov_fig3/fig3_compare.py.
    n_saha(Te, Te) is the paper's n_e' (one temperature); n_saha(Te, Th) is Semenov's n_sa."""
    sys.path.insert(0, os.path.join(HERE, "..", "semenov_fig3"))
    from fig3_compare import n_saha as ns
    return ns(Te, Th, p)


def checks(d):
    xe, te = d["Te"]
    xh, th = d["Th"]
    print("Sec. V text vs digitisation:")
    print(f"  arc end x = L        Fig. 5: {xe.max() * 1e3:.3f} mm   Fig. 7: "
          f"{max(d['ne'][0].max(), d['nsa'][0].max()) * 1e3:.3f} mm")
    print(f"  T_e at x = L         text {TEXT['T_arc']:.0f} K   digitised {te[-1]:.0f} K")
    print(f"  T   at x = L         text {TEXT['T_arc']:.0f} K   digitised {th[-1]:.0f} K")
    print(f"  T_e at the anode     text {TEXT['Te_anode']:.0f} K    digitised {te[0]:.0f} K")
    print(f"  T   at the anode     text {TEXT['T_anode']:.0f} K     digitised {th[0]:.0f} K")
    # The drawn curve flattens towards its end (a quadratic fit gives an even smaller end
    # slope), so the figure is only loosely consistent with the quoted gradient.
    tail = xe > xe.max() - 0.08e-3
    slope = np.polyfit(xe[tail], te[tail], 1)[0]
    print(f"  dT_e/dx at x = L     text {TEXT['dTdx_arc']:.3e} K/m  digitised {slope:.3e} K/m "
          f"(chord over the last 0.08 mm)")
    print(f"  (the digitised anode values are the first points, at x = {xe[0] * 1e6:.1f} um "
          f"(T_e) and {xh[0] * 1e6:.1f} um (T): T rises steeply off the wall and its curve "
          f"merges into the frame line below that)")
    xs = d["nsa"][0]
    ok = (xs >= xe.min()) & (xs <= xe.max())
    te_s, th_s = np.interp(xs[ok], xe, te), np.interp(xs[ok], xh, th)
    for label, ns in (("one-T Saha at T_e (the paper's definition, Sec. III)", n_saha(te_s, te_s)),
                      ("two-T Saha(T_e, T) + Dalton (Semenov's n_sa)", n_saha(te_s, th_s))):
        rel = ns / d["nsa"][1][ok] - 1
        print(f"  n_e' vs {label}: {rel.min():+.0%} ... {rel.max():+.0%} "
              f"(median {np.median(rel):+.0%})")

def verification_figure(fig, img, data, cal, pix):
    """(a) the original panel from the scan; (b) the same, greyed, with every digitised curve
    drawn back over it through the inverse of the fitted axis map, and the ticks used."""
    r0, r1, c0, c1 = FIGS[fig]["crop"]
    crop = img[r0:r1, c0:c1]
    ext = [c0, c1, r1, r0]
    style = {"Te": "#d62728", "Th": "#1f3fff", "ne": "#d62728", "nsa": "#1f3fff"}
    names = {"Te": "$T_e$", "Th": "$T$", "ne": "$n_e$", "nsa": "$n_e'$ (Saha)"}
    plt.rcParams.update({"font.family": "serif", "font.serif": ["STIX", "STIXGeneral",
                         "DejaVu Serif"], "mathtext.fontset": "stix", "font.size": 8})
    f, axs = plt.subplots(1, 2, figsize=(7.0, 2.9))
    axs[0].imshow(crop, cmap="gray", vmin=0, vmax=1, extent=ext)
    axs[0].set_title(f"(a) Dinulescu & Pfender (1980), Fig. {fig} — scan", fontsize=8.5)
    axs[1].imshow(0.55 + 0.45 * crop, cmap="gray", vmin=0, vmax=1, extent=ext)
    for name, (xs, vs) in data.items():
        axs[1].plot(*cal["to_pix"](xs, vs), lw=0.8, color=style[name], label=names[name])
    xt, yt = cal["ticks"]
    axs[1].plot(xt[:, 0], xt[:, 1], "+", ms=4, mew=0.6, color="#2ca02c", label="ticks used")
    axs[1].plot(yt[:, 0], yt[:, 1], "+", ms=4, mew=0.6, color="#2ca02c")
    axs[1].set_title("(b) digitised data re-plotted over the scan (greyed)", fontsize=8.5)
    axs[1].legend(loc="lower left", bbox_to_anchor=(1.0, 0.0), fontsize=7, frameon=False,
                  handlelength=1.6)
    for ax in axs:
        ax.set_xlim(ext[0], ext[1])
        ax.set_ylim(ext[2], ext[3])
        ax.set_axis_off()
    f.tight_layout()
    for e in ("pdf", "png"):
        f.savefig(os.path.join(HERE, f"fig{fig}_verification.{e}"), dpi=300,
                  bbox_inches="tight")
    plt.close(f)
    print(f"wrote fig{fig}_verification.pdf/.png")


if __name__ == "__main__":
    extract(sys.argv[1] if len(sys.argv) > 1 else None)
