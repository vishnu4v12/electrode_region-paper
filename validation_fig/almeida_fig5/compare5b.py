"""
PlasmeX 1-D cathode plotfile(s) vs Almeida, Benilov & Naidis (2008), Fig. 5(b):
Ar, p0 = 1 bar, Tw = 3500 K, jc = 1e7 A/m^2.

    python compare5b.py plt120000000 [plt... ] [-o outdir] [--reextract] [--pdf paper.pdf]

One figure per plotfile, compare5b_<plotfile>.{pdf,png}: (a) n_e, n_i (b) T_e, T_h (c) j_i/j
against the distance from the cathode (PlasmeX cathode = x_hi, plotted as d = x_hi - x).
The engine and all the conventions are in fig5_compare.py; only the figure's axes are set here.
"""
from fig5_compare import main

CASE = dict(
    fig="b",
    jc=1e7,
    table1_jiw=0.216,                         # Almeida et al. Table 1, row 5(b)
    n_ylim=(1e19, 1e24),
    T_ylim=(0, 21), T_major=5,
    ji_ylim=(-0.05, 0.35), ji_major=0.1,      # linear j_i/j
    labels=dict(ne=(1.4e-8, 5.5e19), ni=(1.4e-8, 4.5e21), Te=(1.4e-8, 18.8), Th=(1.4e-8, 4.3)),
    legend=(1, dict(loc="center left", bbox_to_anchor=(0.0, 0.5))),
    wall_note_at=(0.03, 0.2),
)

if __name__ == "__main__":
    main(CASE)
