"""
PlasmeX 1-D cathode plotfile(s) vs Almeida, Benilov & Naidis (2008) J. Phys. D 41 245201 at
jc = 7.8e7 A/m^2 (Ar, p0 = 1 bar, Tw = 3500 K): Fig. 5(c) and curve c of Fig. 6, in one figure.

    python compareJ7.8e7.py <plotfile> [plt...] [-o outdir] [--reextract] [--pdf paper.pdf]

One figure per plotfile, compareJ7.8e7_<plotfile>.{pdf,png}: (a) n_e, n_i (b) T_e, T_h (log,
as in the paper: T_e peaks at 430 kK in the sheath) (c) j_i/j and -E on a twin axis sharing
its zero. The engine, the orientation and the sign conventions are those of
../almeida_fig5/compareJ1e7.py; only the case is set here. The paper curves are this folder's
fig5c/fig6c_extracted.csv (extract.py, run automatically when missing).
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "almeida_fig5"))
from compareJ1e7 import main  # noqa: E402

CASE = dict(
    fig="c", jc=7.8e7, tag="J7.8e7", table1_jiw=0.760,
    datadir=HERE,
    n_ylim=(1e17, 1e24),
    T_log=True, T_ylim=(2.5, 700),                   # kK
    ji_ylim=(-0.1, 0.9), ji_major=0.2,
    E_unit=1e8, E_top=1.35, E_major=0.5,
    labels=dict(ne=(1.4e-8, 2.5e18), ni=(1.4e-8, 6e21), Te=(1.4e-8, 230), Th=(1.4e-8, 5.0),
                ji=(1.5e-6, 0.775), E=(1.2e-8, 1.2)),
    wall_note_at=(0.975, 0.62),
)

if __name__ == "__main__":
    main(CASE)
