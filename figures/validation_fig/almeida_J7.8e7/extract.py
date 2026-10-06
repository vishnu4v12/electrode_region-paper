"""
Extract the jc = 7.8e7 A/m^2 case of Almeida, Benilov & Naidis (2008) J. Phys. D 41 245201
(Ar, p0 = 1 bar, Tw = 3500 K) into THIS folder: Fig. 5(c) (n_e, n_i, n_S, T_e, T_h, T_eq,
j_i/j) and curve c of Fig. 6 (E and phi; drawn as phi/20, stored as phi).

    python extract.py   -> fig5c_extracted.csv, fig6c_extracted.csv and their
                           fig5c/fig6c_verification.{pdf,png}

The extraction code is ../almeida_fig5/fig5_extract.py and fig6_extract.py (same PDF page,
same method, same checks against Table 1 and between the two figures); only the output folder
differs. compareJ7.8e7.py runs this itself when a CSV is missing.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "almeida_fig5"))
import fig5_extract  # noqa: E402
import fig6_extract  # noqa: E402


def extract(pdf=None):
    return {**fig5_extract.extract("c", pdf, outdir=HERE),
            **fig6_extract.extract("c", pdf, outdir=HERE)}


def load(pdf=None):
    return {**fig5_extract.load("c", pdf, outdir=HERE),
            **fig6_extract.load("c", pdf, outdir=HERE)}


if __name__ == "__main__":
    extract(sys.argv[1] if len(sys.argv) > 1 else None)
