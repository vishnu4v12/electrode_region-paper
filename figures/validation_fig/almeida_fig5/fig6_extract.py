"""
Extract curve b or c of Almeida, Benilov & Naidis (2008) J. Phys. D 41 245201, figure 6 (Ar,
p0 = 1 bar, Tw = 3500 K; b: jc = 1e7 A/m^2, the case of Fig. 5(b); c: jc = 7.8e7 A/m^2, Fig.
5(c)): the electric field E(x) (solid) and the electrostatic potential phi(x) (dashed),
straight from the PDF's VECTOR paths, the same way fig5_extract.py does Fig. 5 (same page,
same helpers).

    python fig6_extract.py [b]      -> fig6b_extracted.csv, fig6b_verification.{pdf,png}
    python fig6_extract.py c        -> the same for c, here; ../almeida_J7.8e7/extract.py
                                       writes them into that folder instead (outdir=)

compareJ1e7.py and compareJ7.8e7.py call extract() themselves when the CSV is missing.

The figure has three y axes: E of curves b and c on the LEFT in 10^8 V/m, E of curve a on the
inner RIGHT in 10^6 V/m, phi on the outer RIGHT in V. Curve c is drawn as phi/20 and is scaled
back here: the CSV holds phi itself. Each axis is calibrated by least squares on its own tick
marks. The three solid (blue) and the three dashed (red) curves are drawn in the order a, b, c;
each is chained from its pieces (a solid curve is two strokes, a dashed one ~40 dash strokes)
by end-to-start continuity.

Also written: the two markers on the phi(x) curve -- the circle (charge separation reaches 1%)
and the square (phi maximum, the paper's near-cathode voltage drop). Checks: the square
against the maximum of the extracted phi, phi against -integral(E dx), and the circle against
the 1% charge separation of the matching Fig. 5 panel's densities (fig5{b,c}_extracted.csv,
read from the same folder).

x is the distance from the cathode surface; E is the paper's, i.e. along x (away from the
cathode): negative in the sheath, the field points to the cathode. phi is relative to the
cathode surface.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import fig5_extract as f5
from fig5_extract import BLACK, RED, BLUE, _decades

HERE = os.path.dirname(os.path.abspath(__file__))

PANEL = dict(
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
# curve -> its index in the drawing order, the factor phi is drawn divided by, the caption
CURVES = {
    "b": dict(index=1, phi_scale=1.0,
              caption="Ar, p0 = 1 bar, Tw = 3500 K, jc = 1e7 A/m^2 (curve b)"),
    "c": dict(index=2, phi_scale=20.0,
              caption="Ar, p0 = 1 bar, Tw = 3500 K, jc = 7.8e7 A/m^2 (curve c, drawn as phi/20)"),
}
UNITS = {"E": "V/m", "phi": "V", "phi_1pct": "V", "phi_max": "V"}


def chains(paths, rgb):
    """Curves of one colour, in drawing order: each stroke continues the previous one if it
    starts within 4 pt of where that one ended, else it opens a new chain. Only chains that
    start on the left frame are curves (the arrows pointing to the axes are not)."""
    (xl, xr), (yt, yb) = PANEL["frame_x"], PANEL["frame_y"]
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
    out = [np.concatenate(c) for c in out if f5.near(c[0][0, 0], xl, 0.2)]
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


def csv_path(curve="b", outdir=HERE):
    return os.path.join(outdir, f"fig6{curve}_extracted.csv")


def extract(curve="b", pdf=None, verify_fig=True, outdir=HERE):
    """Extract curve `curve` ('b' or 'c'), run the checks, write the CSV (and the verification
    figure) into `outdir`. Returns {name: (x_m, value)}; the markers are one-point curves."""
    cv = CURVES[curve]
    pdf = pdf or f5.find_pdf()
    paths = f5.trace_page(pdf)
    print(f"--- Almeida et al. 2008, Fig. 6: {cv['caption']}")
    cal = f5.calibrate(paths, PANEL)
    solid, dashed = chains(paths, BLUE), chains(paths, RED)
    assert len(solid) == 3 and len(dashed) == 3, (len(solid), len(dashed))
    e_xy, phi_xy = solid[cv["index"]], dashed[cv["index"]]     # drawing order a, b, c
    mk = markers(paths, phi_xy)

    X, E = cal["x"][0], cal["E8"][0]
    PHI = lambda y: cv["phi_scale"] * cal["phi"][0](y)        # noqa: E731
    data = {"E": (X(e_xy[:, 0]), E(e_xy[:, 1])), "phi": (X(phi_xy[:, 0]), PHI(phi_xy[:, 1]))}
    for q, m in (("phi_1pct", "circle"), ("phi_max", "square")):
        data[q] = (X(np.array([mk[m][0]])), PHI(np.array([mk[m][1]])))
    checks(data, curve, outdir)

    path = csv_path(curve, outdir)
    with open(path, "w") as f:
        f.write(f"# Almeida, Benilov & Naidis 2008, J. Phys. D 41 245201, Fig. 6, curve {curve}: "
                f"{cv['caption']}\n")
        f.write("# extracted from the PDF vector paths (fig6_extract.py); x = distance from "
                "cathode; E along x (paper sign); phi relative to the cathode, in V\n")
        f.write("# phi_1pct: circle, charge separation reaches 1%; phi_max: square, phi maximum\n")
        f.write("quantity,unit,x_m,value\n")
        for name, (xs, vs) in data.items():
            for a, b in zip(xs, vs):
                f.write(f"{name},{UNITS[name]},{a:.6e},{b:.6e}\n")
    print(f"wrote {path}  ({sum(len(v[0]) for v in data.values())} points)")
    if verify_fig:
        verification_figure(curve, pdf, data, cal, outdir)
    return data


def load(curve="b", pdf=None, outdir=HERE):
    """The extracted curve of Fig. 6, extracting first if the CSV is not there."""
    path = csv_path(curve, outdir)
    if not os.path.isfile(path):
        return extract(curve, pdf, outdir=outdir)
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


def checks(data, curve, outdir):
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
    # the circle: |n_i - n_e| / n_i = 1 %, from the densities of the matching Fig. 5 panel
    p5 = f5.load(curve, outdir=outdir)
    x5 = np.geomspace(1e-7, 1e-4, 4000)
    ne = np.exp(np.interp(np.log(x5), np.log(p5["ne"][0]), np.log(p5["ne"][1])))
    ni = np.exp(np.interp(np.log(x5), np.log(p5["ni"][0]), np.log(p5["ni"][1])))
    sep = np.abs(ni - ne) / ni
    x1 = x5[np.argmax(sep < 0.01)]
    (xc,), (pc,) = data["phi_1pct"]
    print(f"  1 % charge separation   circle at x = {xc:.3e} m ({pc:.2f} V);  "
          f"Fig. 5({curve}) densities: x = {x1:.3e} m")


def verification_figure(curve, pdf, data, cal, outdir=HERE):
    """(a) the original figure; (b) greyed, with the extracted curve drawn back over it."""
    img, ext = f5.panel_raster(pdf, PANEL)
    plt.rcParams.update({"font.family": "serif", "font.serif": ["STIX", "STIXGeneral",
                         "DejaVu Serif"], "mathtext.fontset": "stix", "font.size": 8})
    f, axs = plt.subplots(1, 2, figsize=(7.0, 3.1))
    axs[0].imshow(img, extent=ext)
    axs[0].set_title("(a) Almeida et al. (2008), Fig. 6 — original", fontsize=8.5)
    axs[1].imshow(0.55 + 0.45 * img[..., :3].mean(axis=2), cmap="gray", vmin=0, vmax=1,
                  extent=ext)
    sc = CURVES[curve]["phi_scale"]
    X, E = cal["x"][1], cal["E8"][1]
    PHI = lambda v: cal["phi"][1](v / sc)                     # noqa: E731
    axs[1].plot(X(data["E"][0]), E(data["E"][1]), lw=0.9, color="#1f3fff", label=f"E ({curve})")
    axs[1].plot(X(data["phi"][0]), PHI(data["phi"][1]), lw=0.9, color="#d62728",
                label=rf"$\varphi$ ({curve})" + (f" / {sc:.0f}" if sc != 1 else ""))
    for q, mk in (("phi_1pct", "o"), ("phi_max", "s")):
        axs[1].plot(X(data[q][0]), PHI(data[q][1]), mk, ms=5, mfc="none", mec="#2ca02c",
                    mew=0.9, label={"phi_1pct": "1 % sep.", "phi_max": r"$\varphi$ max"}[q])
    axs[1].set_title(f"(b) extracted curve {curve} re-plotted over the original (greyed)",
                     fontsize=8.5)
    axs[1].legend(loc="lower left", bbox_to_anchor=(1.0, 0.0), fontsize=7, frameon=False,
                  handlelength=1.6)
    for ax in axs:
        ax.set_xlim(ext[0], ext[1])
        ax.set_ylim(ext[2], ext[3])
        ax.set_axis_off()
    f.tight_layout()
    for e in ("pdf", "png"):
        f.savefig(os.path.join(outdir, f"fig6{curve}_verification.{e}"), dpi=300,
                  bbox_inches="tight")
    plt.close(f)
    print(f"wrote fig6{curve}_verification.pdf/.png")


if __name__ == "__main__":
    for cv in sys.argv[1:] or ["b"]:
        extract(cv)
