"""
PlasmeX 1-D anode plotfile(s) vs Semenov, Krivtsun & Reisgen (2016), Fig. 3(b):
Ar, p0 = 1 atm, anode surface T_s = 1000 K, j = 1500 A/cm^2 (1.5e7 A/m^2), unified model.

    python compare3b.py plt... [-o outdir] [--reextract] [--pdf paper.pdf]

One figure per plotfile, compare3b_<plotfile>.{pdf,png}: (a) n_e, n_i, n_sa (b) T_e, T_h
against the distance from the anode (PlasmeX anode = x_hi, plotted as d = x_hi - x).
The engine and all the conventions are in fig3_compare.py; only the figure's axes are set here.
"""
from fig3_compare import main

CASE = dict(
    fig="b",
    j=1.5e7, j_Acm2=1500,
    L_min=2e-3,                               # uniform plasma from x ~ 1-2 mm (Fig. 3b)
    n_ylim=(1e19, 1e24),
    T_ylim=(0, 19), T_major=4,
    labels=dict(ne=(1.4e-8, 4.0e20), ni=(1.4e-8, 1.4e21), nsa=(1.6e-6, 1.3e20),
                Te=(1.4e-8, 5.9), Th=(1.4e-8, 1.3)),
    legend=(0, dict(loc="upper left", bbox_to_anchor=(0.0, 0.88))),
    wall_note_at=(0.03, 0.55),
)

if __name__ == "__main__":
    main(CASE)
