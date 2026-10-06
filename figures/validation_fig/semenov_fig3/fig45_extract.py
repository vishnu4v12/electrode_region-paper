"""
Extract the j = 1500 A/cm^2 (1.5e7 A/m^2) curves of Semenov, Krivtsun & Reisgen (2016)
J. Phys. D 49 105204 (Ar, 1 atm, anode surface T_s = 1 kK, unified model), the case of
Fig. 3(b), from two more figures on page 10 (printed page 9):

    Fig. 4   j_i/j   ion current density normalised to j   (red curve; blue is 700 A/cm^2)
    Fig. 5   E       electric field, 10^6 V/m               (red curve; also 500/1000/2000)

straight from the PDF's VECTOR paths -- the method of fig3_extract.py, whose helpers are used:
each curve is ONE continuous polyline picked by colour, each axis calibrated by least squares
on its own tick marks (residuals printed).

    python fig45_extract.py      -> fig45_j1500_extracted.csv, fig45_j1500_verification.{pdf,png}

compareJ1.5e7.py calls extract() itself when the CSV is missing.

x is the distance from the anode surface (anode on the LEFT). Both quantities carry the
paper's signs: x points away from the anode and j > 0, so ions flowing INTO the anode give
j_i/j < 0, and E < 0 near the anode (the paper: E > 0 in the uniform plasma, changing sign
inside the layer). The paper's j_e/j (Fig. 4, right axis) is j_i/j - 1 (J_i - J_e = j/e).
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "almeida_fig5"))
import fig3_extract as f3  # noqa: E402
import fig5_extract as f5  # noqa: E402  (panel_raster)
from fig3_extract import RED, _decades  # noqa: E402

PAGE = 10                                   # printed page 9
PANELS = {
    "ji_j": dict(
        fig="Fig. 4", unit="-",
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
        fig="Fig. 5", unit="V/m",
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
CAPTION = "Ar, 1 atm, T_s = 1 kK, j = 1500 A/cm^2 (1.5e7 A/m^2)"


def pick(paths, panel):
    """The red (1500 A/cm^2) polyline inside the frame, thinned, in x order."""
    (xl, xr), (yt, yb) = panel["frame_x"], panel["frame_y"]
    found = [p for p in paths if len(p["xy"]) > 100 and np.allclose(p["color"], RED)
             and xl - 0.2 <= p["xy"][:, 0].min() and p["xy"][:, 0].max() <= xr + 0.2
             and yt - 0.2 <= p["xy"][:, 1].min() and p["xy"][:, 1].max() <= yb + 0.2]
    assert len(found) == 1, [(p["color"], len(p["xy"])) for p in found]
    xy = found[0]["xy"]
    return f3.thin(xy[np.argsort(xy[:, 0], kind="stable")])


def csv_path():
    return os.path.join(HERE, "fig45_j1500_extracted.csv")


def extract(pdf=None, verify_fig=True):
    """Extract both curves, run the checks, write the CSV (and the verification figure).
    Returns {'ji_j': (x_m, value), 'E': (x_m, V/m)}."""
    pdf = pdf or f3.find_pdf()
    paths = f3.trace_page(pdf, PAGE)
    data, cals = {}, {}
    for name, panel in PANELS.items():
        print(f"--- Semenov et al. 2016, {panel['fig']}: {CAPTION}")
        cal = f3.calibrate(paths, panel)
        xy = pick(paths, panel)
        data[name] = (cal["x"][0](xy[:, 0]), cal["y"][0](xy[:, 1]))
        cals[name] = cal
    checks(data)

    with open(csv_path(), "w") as f:
        f.write(f"# Semenov, Krivtsun & Reisgen 2016, J. Phys. D 49 105204, Fig. 4 (ji_j) and "
                f"Fig. 5 (E), j = 1500 A/cm^2: {CAPTION}\n")
        f.write("# extracted from the PDF vector paths (fig45_extract.py); x = distance from "
                "anode; paper signs (j > 0, x away from the anode)\n")
        f.write("quantity,unit,x_m,value\n")
        for name, (xs, vs) in data.items():
            for a, b in zip(xs, vs):
                f.write(f"{name},{PANELS[name]['unit']},{a:.6e},{b:.6e}\n")
    print(f"wrote {csv_path()}  ({sum(len(v[0]) for v in data.values())} points)")
    if verify_fig:
        verification_figure(pdf, data, cals)
    return data


def load(pdf=None):
    """The extracted curves, extracting first if the CSV is not there."""
    if not os.path.isfile(csv_path()):
        return extract(pdf)
    cols = {}
    with open(csv_path()) as f:
        for line in f:
            if line.startswith(("#", "quantity")):
                continue
            q, _, x, v = line.strip().split(",")
            cols.setdefault(q, ([], []))
            cols[q][0].append(float(x))
            cols[q][1].append(float(v))
    return {q: (np.array(x), np.array(v)) for q, (x, v) in cols.items()}


def checks(data):
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


def verification_figure(pdf, data, cals):
    """Per figure: the original panel, and the same raster greyed with the extracted curve
    drawn back over it through the inverse of the fitted axis maps."""
    plt.rcParams.update({"font.family": "serif", "font.serif": ["STIX", "STIXGeneral",
                         "DejaVu Serif"], "mathtext.fontset": "stix", "font.size": 8})
    f, axs = plt.subplots(2, 2, figsize=(7.0, 6.0))
    for row, (name, panel) in zip(axs, PANELS.items()):
        img, ext = f5.panel_raster(pdf, panel, page=PAGE)
        row[0].imshow(img, extent=ext)
        row[0].set_title(f"Semenov et al. (2016), {panel['fig']} — original", fontsize=8.5)
        row[1].imshow(0.55 + 0.45 * img[..., :3].mean(axis=2), cmap="gray", vmin=0, vmax=1,
                      extent=ext)
        cal = cals[name]
        row[1].plot(cal["x"][1](data[name][0]), cal["y"][1](data[name][1]), lw=0.9,
                    color="#d62728")
        row[1].set_title("1500 A cm$^{-2}$ extracted, over the original (greyed)",
                         fontsize=8.5)
        for ax in row:
            ax.set_xlim(ext[0], ext[1])
            ax.set_ylim(ext[2], ext[3])
            ax.set_axis_off()
    f.tight_layout()
    for e in ("pdf", "png"):
        f.savefig(os.path.join(HERE, f"fig45_j1500_verification.{e}"), dpi=300,
                  bbox_inches="tight")
    plt.close(f)
    print("wrote fig45_j1500_verification.pdf/.png")


if __name__ == "__main__":
    extract(sys.argv[1] if len(sys.argv) > 1 else None)
