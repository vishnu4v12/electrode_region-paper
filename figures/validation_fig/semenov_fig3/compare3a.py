"""
PlasmeX 1-D anode plotfile(s) vs Semenov, Krivtsun & Reisgen (2016), Fig. 3(a):
Ar, p0 = 1 atm, anode surface T_s = 1000 K, j = 700 A/cm^2 (7e6 A/m^2), unified model.

    python compare3a.py plt... [-o outdir] [--reextract] [--pdf paper.pdf]

One figure per plotfile, compare3a_<plotfile>.{pdf,png}: (a) n_e, n_i, n_sa (b) T_e, T_h
against the distance from the anode (PlasmeX anode = x_hi, plotted as d = x_hi - x).
The engine and all the conventions are in fig3_compare.py; only the figure's axes are set here.
"""
from fig3_compare import main

CASE = dict(
    fig="a",
    j=7e6, j_Acm2=700,
    L_min=3e-3,                               # uniform plasma from x ~ 2-3 mm (Fig. 3a)
    n_ylim=(1e19, 1e24),
    T_ylim=(0, 16), T_major=4,
    labels=dict(ne=(1.4e-8, 2.2e20), ni=(1.4e-8, 7.5e20), nsa=(3.5e-6, 1.3e20),
                Te=(1.4e-8, 5.3), Th=(1.4e-8, 1.3)),
    legend=(0, dict(loc="upper left", bbox_to_anchor=(0.0, 0.88))),
    wall_note_at=(0.03, 0.55),
)

if __name__ == "__main__":
    main(CASE)
