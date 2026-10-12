# Electrode-region paper

Draft manuscript: "A unified model for species fluxes, ion impact energy and heat
transfer at the electrode of an atmospheric-pressure arc discharge", for
*Journal of Physics D: Applied Physics* (IOP `iopjournal` class).
The solver is PlasmeX, in `../PlasmeX` (not part of this repository).

## Layout
| path | contents |
|---|---|
| `CLAUDE.md` | guidance for Claude Code; points to `../PlasmeX/CLAUDE.md` and the run directories |
| `main.tex`, `sections/*.tex`, `references.bib` | the manuscript (`sections/appendix.tex` is not `\input` yet) |
| `iopjournal.cls`, `orcid.pdf` | needed by the build; keep them next to `main.tex` |
| `main.pdf` | last compiled manuscript |
| `response.md` | current Q&A and explanations (render the equations in VS Code) |
| `literature_review.md` | novelty positioning and the validation plan |
| `archive/` | earlier `response.md` Q&A, superseded but kept for reference |
| `References/` | PDFs of the key papers (Almeida 2008, Semenov 2016 anode, Dinulescu & Pfender 1980 anode, transport chapter) |
| `figures/validation_fig/` | validation scripts and data (below) |
| `iop_template/` | IOP template zip and author guidelines |

## Validation figures
`figures/validation_fig/` holds one folder per benchmark. Each has a single `compare.py` that
takes PlasmeX 1-D plotfile(s) and writes the comparison figure, and a `paper_data/` folder with
the curves extracted from the paper (CSV), the script that extracts them (`extract.py`) and its
verification figures (the extracted curves drawn over the original). The extraction reads the
PDF's vector paths (`mutool trace`, axes calibrated on the figure's own ticks); it only needs
re-running if the extraction itself changes.

```
validation_fig/
  plasmex.py, pdftrace.py         shared: PlasmeX loader, `mutool trace` parser
  almeida_J1e7/
    compare.py <plotfile...>      -> compareJ1e7_<plotfile>.{pdf,png}
    plt19500/                     the plotfile behind the current figure
    paper_data/extract.py         -> fig5b_extracted.csv, fig6b_extracted.csv (+ verification)
  semenov_J1.5e7/
    compare.py <plotfile...>      -> compareJ1.5e7_<plotfile>.{pdf,png}
    plt73000/
    paper_data/extract.py         -> fig3b_extracted.csv, fig45_j1500_extracted.csv (+ verification)
```

Nothing from PlasmeX is copied here. `validation_fig/plasmex.py` loads PlasmeX scripts (`bohm.py`
to sample plotfiles) from the live source tree, `$PLASMEX_ROOT` or `~/work/PlasmeX`, so later
PlasmeX changes are picked up.

- **Cathode**, `almeida_J1e7/`: Almeida, Benilov & Naidis (2008), jc = 1e7 A m^-2 (Ar, 1 bar,
  T_w = 3500 K). Fig. 5(b) (checked against Table 1) and Fig. 6 curve b, E and phi (checked:
  phi against -int E dx, the phi-max square, the 1 % charge-separation circle against Fig. 5(b)).
  The figure: n_e, n_i / T_e, T_h / j_i/j with -E on a twin axis; the terminal table adds phi.
- **Anode**, `semenov_J1.5e7/`: Semenov, Krivtsun & Reisgen (2016), j = 1500 A cm^-2 (Ar, 1 atm,
  T_s = 1000 K). Fig. 3(b) (checked against Table 2) and the 1500 A cm^-2 curves of Fig. 4
  (j_i/j) and Fig. 5 (E). The figure: n_e, n_i, n_sa / T_e, T_h / j_i/j with E on a twin axis.

## Compiling
```
latexmk -pdf main.tex        # or pdflatex -> bibtex -> pdflatex x2 (the VS Code recipe)
```
Build intermediates (`*.aux`, `*.log`, `*.synctex.gz`, ...) are git-ignored.
