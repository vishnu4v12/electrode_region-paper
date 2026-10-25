"""
Schematic of the arc discharge and of the two near-electrode domains (figure fig:schematic).

    python arc_schematic.py [-o outdir]

(a) the whole arc: thoriated-tungsten cathode (top), copper anode (bottom), the LTE arc column between
    them, and the cathode and anode regions (thickness exaggerated) that the 1-D model resolves along the
    axis, from the electrode surface (x = 0) to the bulk-plasma boundary (x = L).
(b) the cathode domain and (c) the anode domain, not to scale: the space-charge sheath, the
    quasi-neutral non-equilibrium layer and the near-LTE plasma, the data imposed at each end, the
    particle fluxes at the surface, the radial loss, and qualitative T_e and T_h.
Figure conventions as in the results figures: the cathode on the left of its panel, the anode on the
right of its panel (mirrored axis).
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
C_E, C_I = "#2a78d6", "#eb6834"          # electrons blue, heavy particles / ions orange (results figures)
W_COL, CU_COL = "#8c8f94", "#b87333"     # tungsten grey, copper
PLASMA = "#f6b26b"
SHEATH = "#cfe0f5"
LAYER = "#fde3c4"
LTE = "#fff4e6"
INK = "#2b2b2b"


def rc():
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["STIX", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "stix", "font.size": 8.0,
    })


def arrow(ax, p0, p1, col, lw=1.1, ms=8, **kw):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=ms, color=col, lw=lw,
                                 shrinkA=0, shrinkB=0, **kw))


def brace(ax, x, y0, y1, text, side=1, **kw):
    """Vertical bracket at x from y0 to y1, label on `side` (+1 right, -1 left)."""
    d = 0.06 * side
    ax.plot([x, x + d, x + d, x], [y0, y0, y1, y1], color=INK, lw=0.7)
    ax.text(x + 2.4 * d, 0.5 * (y0 + y1), text, ha="left" if side > 0 else "right", va="center", **kw)


def panel_arc(ax):
    """(a) The arc: cathode rod on top, anode plate at the bottom, bell-shaped column between."""
    yc, ya = 7.3, 1.0                       # cathode tip, anode surface
    # column: narrow at the cathode spot, spreading towards the anode
    y = np.linspace(ya, yc, 200)
    s = (y - ya) / (yc - ya)
    r = 0.11 + 0.80 * (1 - s) ** 0.5 + 0.25 * (1 - s) ** 6   # bell: fast spread below the spot, flare at the anode
    for k, a in ((1.25, 0.18), (1.1, 0.30), (1.0, 0.55)):
        ax.add_patch(Polygon(np.c_[np.r_[k * r, -k * r[::-1]], np.r_[y, y[::-1]]], closed=True,
                             fc=PLASMA, ec="none", alpha=a))
    ax.add_patch(Polygon(np.c_[np.r_[0.45 * r, -0.45 * r[::-1]], np.r_[y, y[::-1]]], closed=True,
                         fc="#fbe7c6", ec="none", alpha=0.9))
    # cathode: rod + conical tip
    ax.add_patch(Polygon([(-0.32, 10.2), (0.32, 10.2), (0.32, 8.4), (0.06, yc), (-0.06, yc), (-0.32, 8.4)],
                         closed=True, fc=W_COL, ec=INK, lw=0.7))
    ax.text(0.42, 9.5, "cathode\n(W–ThO$_2$)", ha="left", va="center")
    # anode plate
    ax.add_patch(Rectangle((-1.55, 0.25), 3.1, ya - 0.25, fc=CU_COL, ec=INK, lw=0.7))
    ax.text(0.0, 0.62, "anode (Cu)", ha="center", va="center", color="white")
    # near-electrode regions (thickness exaggerated)
    hc, ha_ = 0.55, 0.55
    ax.add_patch(Rectangle((-0.32, yc - hc), 0.64, hc, fc="none", ec=C_E, lw=1.0, ls="--"))
    ax.add_patch(Rectangle((-1.25, ya), 2.5, ha_, fc="none", ec=C_E, lw=1.0, ls="--"))
    ax.text(-0.4, yc - 0.5 * hc, "cathode\nregion (b)", ha="right", va="center", color=C_E)
    ax.text(-1.3, ya + 0.5 * ha_ + 0.65, "anode region (c)", ha="left", va="bottom", color=C_E)
    # the 1-D domains along the axis
    arrow(ax, (0.0, yc), (0.0, yc - hc - 0.35), INK, lw=0.8, ms=6)
    ax.text(0.07, yc - hc - 0.3, "$x$", ha="left", va="center")
    arrow(ax, (0.9, ya), (0.9, ya + ha_ + 0.35), INK, lw=0.8, ms=6)
    ax.text(0.97, ya + ha_ + 0.3, "$x$", ha="left", va="center")
    ax.text(0.0, 4.2, "LTE arc\ncolumn", ha="center", va="center", fontsize=8.5)
    ax.text(0.0, 3.35, "(bulk model)", ha="center", va="center", fontsize=7)
    # current
    arrow(ax, (1.35, 3.0), (1.35, 6.0), INK, lw=0.9, ms=7)
    ax.text(1.42, 4.5, "$I$", ha="left", va="center")
    ax.set_xlim(-1.7, 1.75)
    ax.set_ylim(0.1, 10.3)
    ax.set_aspect("equal")
    ax.axis("off")


def profile(u, a, b, w, x0):
    return a + (b - a) / (1 + np.exp(-(u - x0) / w))


def panel_layer(ax, el):
    """(b) cathode / (c) anode domain on a not-to-scale coordinate u in [0, 1] from the surface to x = L.
    Drawn in u; the anode panel is mirrored on the page (anode on the right)."""
    us, ul = 0.11, 0.62                     # sheath edge, end of the non-equilibrium layer
    y0, y1 = 0.0, 1.0
    ax.add_patch(Rectangle((0, y0), us, y1 - y0, fc=SHEATH, ec="none"))
    ax.add_patch(Rectangle((us, y0), ul - us, y1 - y0, fc=LAYER, ec="none"))
    ax.add_patch(Rectangle((ul, y0), 1 - ul, y1 - y0, fc=LTE, ec="none"))
    for u in (us, ul):
        ax.plot([u, u], [y0, y1], color=INK, lw=0.5, ls=":")
    ax.plot([1, 1], [y0, y1], color=INK, lw=1.0, ls="--")
    elec = W_COL if el == "cathode" else CU_COL
    ax.add_patch(Rectangle((-0.07, y0 - 0.02), 0.07, y1 - y0 + 0.04, fc=elec, ec=INK, lw=0.7, clip_on=False))
    # region labels
    lab_y = 0.93
    ax.text(us / 2, lab_y, "sheath", ha="center", va="top", fontsize=7)
    ax.text(us / 2, lab_y - 0.1, r"$\lesssim1\,\mu$m", ha="center", va="top", fontsize=6.5)
    mid = 0.5 * (us + ul)
    if el == "cathode":
        ax.text(mid, lab_y, "ionization layer", ha="center", va="top", fontsize=7.5)
        ax.text(mid, lab_y - 0.1, r"$T_e\neq T_h$, $n_e\neq n_{\rm eq}$, $\sim10^{-4}$ m",
                ha="center", va="top", fontsize=6.5)
    else:
        ax.text(mid, lab_y, "quasi-neutral anode layer", ha="center", va="top", fontsize=7.5)
        ax.text(mid, lab_y - 0.1, r"$T_e\neq T_h$, $n_e>n_{\rm eq}$", ha="center", va="top", fontsize=6.5)
    ax.text(0.5 * (ul + 1), lab_y, "near-LTE plasma", ha="center", va="top", fontsize=7.5)
    ax.text(0.5 * (ul + 1), lab_y - 0.1, "radial loss $Q_e$, $Q_h$", ha="center", va="top", fontsize=6.5,
            color="#e34948")
    for u in (0.73, 0.89):
        arrow(ax, (u, 0.60), (u, 0.74), "#e34948", lw=0.8, ms=5)
    # qualitative temperatures (lower half of the strip)
    u = np.linspace(0.0, 1.0, 400)
    if el == "cathode":
        Te = profile(u, 0.53, 0.56, 0.08, 0.6) - 0.012 * np.exp(-((u - 0.3) / 0.06) ** 2)
        Th = 0.30 + (0.56 - 0.30) * (1 - np.exp(-u / 0.17))
    else:
        Te = 0.40 + (0.56 - 0.40) * (1 - np.exp(-(u + 0.02) / 0.4)) / (1 - np.exp(-1.02 / 0.4))
        Th = 0.24 + (0.56 - 0.24) * (1 - np.exp(-u / 0.28)) / (1 - np.exp(-1 / 0.28))
    ax.plot(u, Te, color=C_E, lw=1.2)
    ax.plot(u, Th, color=C_I, lw=1.2)
    k = 140 if el == "cathode" else 230
    ax.text(u[k], Te[k] + 0.035, "$T_e$", color=C_E, ha="center", va="bottom", fontsize=7.5)
    ax.text(u[k], Th[k] - 0.035, "$T_h$", color=C_I, ha="center", va="top", fontsize=7.5)
    # fluxes at the surface (in the sheath, below the temperature curves); labels point away from the electrode
    ha = "left" if el == "cathode" else "right"   # the anode panel is mirrored
    fl = ((C_I, "ions", True), (C_E, "emitted e$^-$", False)) if el == "cathode" else \
         ((C_E, "e$^-$", True), (C_I, "ions", True))
    for (col, lab, inward), yy in zip(fl, (0.2, 0.08)):
        p0, p1 = ((0.095, yy), (0.008, yy)) if inward else ((0.008, yy), (0.095, yy))
        arrow(ax, p0, p1, col, lw=1.0, ms=6)
        ax.text(0.105, yy, lab, color=col, ha=ha, va="center", fontsize=6.5)
    # boundary data
    surf = ("$x=0$: $T_{\\rm wall}$ (heat balance),\nemission, wall fluxes" if el == "cathode"
            else "$x=0$: $T_{\\rm wall}$ prescribed,\nwall fluxes")
    ax.text(0.0, -0.08, surf, ha="center" if False else ("left" if el == "cathode" else "right"),
            va="top", fontsize=6.5)
    ax.text(1.0, -0.08, "$x=L$: $T_{\\rm bulk}$, $n_{\\rm eq}(T_{\\rm bulk})$,\n$J_{\\rm app}$ (LTE bulk)",
            ha="right" if el == "cathode" else "left", va="top", fontsize=6.5)
    ax.set_xlim(-0.08, 1.02)
    ax.set_ylim(-0.3, 1.02)
    if el == "anode":
        ax.invert_xaxis()
    ax.axis("off")
    ax.set_title(("(b) cathode region" if el == "cathode" else "(c) anode region") + " (not to scale)",
                 fontsize=8.5, loc="left" if el == "cathode" else "right", pad=3)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("-o", "--outdir", default=HERE)
    a = ap.parse_args(argv)
    rc()
    fig = plt.figure(figsize=(7.0, 3.9))
    gs = fig.add_gridspec(2, 2, width_ratios=(0.9, 2.1), wspace=0.05, hspace=0.45)
    ax = fig.add_subplot(gs[:, 0])
    panel_arc(ax)
    ax.set_title("(a)", fontsize=8.5, loc="left", pad=3)
    panel_layer(fig.add_subplot(gs[0, 1]), "cathode")
    panel_layer(fig.add_subplot(gs[1, 1]), "anode")
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(a.outdir, f"arc_schematic.{ext}"), bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"wrote arc_schematic (.pdf/.png) in {a.outdir}")


if __name__ == "__main__":
    main()
