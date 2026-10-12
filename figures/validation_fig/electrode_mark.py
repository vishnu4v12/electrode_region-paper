"""
Marks the electrode edge of a profile plot, shared by the validation and results figures.

    from electrode_mark import mark_electrode
    mark_electrode(ax, "left", label="cathode")    # label: top panel only

The electrode side of the axes gets a hatched strip (the drawing convention for a solid wall)
and a thicker frame line. The strip sits BEHIND the curves (zorder 0.5), so it hides no data
even where the plotted range starts at the first cell. The optional label goes just above the
frame at the electrode corner.
"""
from matplotlib.patches import Rectangle

WIDTH = 0.03                       # strip width, fraction of the axes
SPINE_LW = 2.0
FILL, HATCH_C, TEXT_C = "#ececec", "#9a9a9a", "#3a3a3a"


def mark_electrode(ax, side, label=None, spine_axes=()):
    """side: 'left' or 'right' on screen (a mirrored axis puts the electrode on the right).
    spine_axes: twin axes whose spine on that side must be thickened as well."""
    x0 = 0.0 if side == "left" else 1.0 - WIDTH
    ax.add_patch(Rectangle((x0, 0.0), WIDTH, 1.0, transform=ax.transAxes, facecolor=FILL,
                           edgecolor=HATCH_C, hatch="//////", lw=0, zorder=0.5))
    for a in (ax, *spine_axes):
        a.spines[side].set_linewidth(SPINE_LW)
    if label:
        ax.text(0.0 if side == "left" else 1.0, 1.015, label, transform=ax.transAxes,
                ha=side, va="bottom", fontsize=8, fontstyle="italic", color=TEXT_C)
