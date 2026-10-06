# Electrode-region paper

Draft manuscript: "A unified model for species fluxes, ion impact energy and heat
transfer at the electrode of an atmospheric-pressure arc discharge", for
*Journal of Physics D: Applied Physics* (IOP `iopjournal` class).
The solver is PlasmeX, in `../PlasmeX` (not part of this repository).

## Layout
| path | contents |
|---|---|
| `main.tex`, `sections/*.tex`, `references.bib` | the manuscript (`sections/appendix.tex` is not `\input` yet) |
| `iopjournal.cls`, `orcid.pdf` | needed by the build; keep them next to `main.tex` |
| `main.pdf` | last compiled manuscript |
| `response.md` | current Q&A and explanations (render the equations in VS Code) |
| `literature_review.md` | novelty positioning and the validation plan |
| `archive/` | earlier `response.md` Q&A, superseded but kept for reference |
| `References/` | PDFs of the key papers (Almeida 2008, Semenov 2016 anode, Dinulescu & Pfender 1980 anode, transport chapter) |
| `validation_fig/` | validation scripts and data (below) |
| `iop_template/` | IOP template zip and author guidelines |

## Validation figures
Both benchmarks read the published curves straight from the PDF's vector paths
(`mutool trace`, axes calibrated on the figure's own ticks), write them to CSV with a
verification figure (extracted curves drawn over the original), and overlay a PlasmeX
1-D plotfile against the distance from the electrode.

Nothing from PlasmeX is copied here. `validation_fig/plasmex.py` loads PlasmeX scripts and
data (`bohm.py` to sample plotfiles, `almeida_arc_bc.py`, the transport tables) from the live
source tree, `$PLASMEX_ROOT` or `~/work/PlasmeX`, so later PlasmeX changes are picked up.

- **Cathode**, `validation_fig/almeida_fig5/`: Almeida, Benilov & Naidis (2008) Fig. 5(a)/(b).
  `fig5_extract.py` -> `fig5{a,b}_extracted.csv`; `compare5a.py <plotfile>` for 5(a).
  Fig. 6 curve b (E and phi at jc = 1e7 A m^-2): `fig6_extract.py` -> `fig6b_extracted.csv`
  (checked: phi against -int E dx, the phi-max square, the 1 % charge-separation circle against
  Fig. 5(b)); `compareJ1e7.py <plotfile>` plots n_e, n_i / T_e, T_h / j_i/j with -E on a twin axis.
- **Cathode**, `validation_fig/almeida_J7.8e7/`: the jc = 7.8e7 A m^-2 case, Fig. 5(c) and Fig. 6
  curve c (drawn as phi/20, stored as phi). `extract.py` runs the two extractors above into this
  folder -> `fig5c_extracted.csv`, `fig6c_extracted.csv` (same checks, Table 1);
  `compareJ7.8e7.py <plotfile>` is `compareJ1e7.py` with this case's axes (T on a log scale).
- **Anode**, `validation_fig/semenov_fig3/`: Semenov, Krivtsun & Reisgen (2016) Fig. 3(a)/(b)
  (Ar, 1 atm, T_s = 1000 K, j = 700 / 1500 A cm^-2). `fig3_extract.py` ->
  `fig3{a,b}_extracted.csv` (checked against the paper's Table 2);
  `compare3a.py` / `compare3b.py <plotfile>` plot n_e, n_i, n_sa and T_e, T_h;
  `fig45_extract.py` -> `fig45_j1500_extracted.csv`: the 1500 A cm^-2 curves of Fig. 4 (j_i/j)
  and Fig. 5 (E); `compareJ1.5e7.py <plotfile>` plots n_e, n_i, n_sa / T_e, T_h / j_i/j with E
  on a twin axis;
  `semenov_arc_bc.py` solves the paper's arc-side column, eqs (19)-(24), checks it against
  Table 2 and prints the `anode_inputs` lines (`--deck`); it reuses
  `~/work/PlasmeX/exe/almeida_arc_bc.py`, the Almeida counterpart.
- **Anode**, `validation_fig/pfender_fig5_7/`: Dinulescu & Pfender (1980) Figs. 5 and 7
  (Ar, 100 kPa, T_w = 400 K, j = 900 A cm^-2). The PDF is a bitmap scan, so
  `fig57_extract.py` digitises the raster -> `fig{5,7}_extracted.csv` (T_e, T_h; n_e, n_e'
  from Saha at T_e), checked against the values in the paper's text;
  `compare57.py <plotfile> [--logx]` plots n_e, n_sa and T_e, T_h;
  `compare_pfender_semenov3b.py [plotfile]` plots Pfender against Semenov Fig. 3(b)
  (j = 900 vs 1500 A cm^-2), optionally with PlasmeX over both.

## Compiling
```
latexmk -pdf main.tex        # or pdflatex -> bibtex -> pdflatex x2 (the VS Code recipe)
```
Build intermediates (`*.aux`, `*.log`, `*.synctex.gz`, ...) are git-ignored.
