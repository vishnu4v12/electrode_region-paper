"""
PlasmeX 1-D cathode plotfile(s) vs Almeida, Benilov & Naidis (2008), Fig. 5(a):
Ar, p0 = 1 bar, Tw = 3500 K, jc = 1e6 A/m^2.

    python compare5a.py plt130000000 [plt... ] [-o outdir] [--reextract] [--pdf paper.pdf]

One figure per plotfile, compare5a_<plotfile>.{pdf,png}: (a) n_e, n_i (b) T_e, T_h (c) j_i/j
against the distance from the cathode (PlasmeX cathode = x_hi, plotted as d = x_hi - x).
The engine and all the conventions are in fig5_compare.py; only the figure's axes are set here.

j_i/j is on a SYMLOG axis here: the paper's ion fraction is ~4e-4 and turns negative, so a
linear axis that holds a PlasmeX value orders of magnitude larger flattens the paper's curve
onto zero, and a log axis cannot show its negative tail.
"""
from fig5_compare import main

CASE = dict(
    fig="a",
    jc=1e6,
    table1_jiw=4e-4,                          # Almeida et al. Table 1, row 5(a)
    n_ylim=(1e18, 1e23),
    T_ylim=(0, 13), T_major=2,
    ji_ylim=(-3e-4, 0.5), ji_symlog=1e-4,     # symlog j_i/j, linear inside +-1e-4
    labels=dict(ne=(1.4e-8, 7e19), ni=(1.4e-8, 7e20), Te=(1.4e-8, 9.3), Th=(1.4e-8, 2.2)),
    legend=(1, dict(loc="lower right", bbox_to_anchor=(1.0, 0.02))),
    wall_note_at=(0.03, 0.56),
)

if __name__ == "__main__":
    main(CASE)
