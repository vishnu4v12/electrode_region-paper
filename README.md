# Manuscript build notes

Draft manuscript for "A unified model for species fluxes, ion impact energy
and heat transfer at the electrode of an atmospheric-pressure arc
discharge," targeting *Journal of Physics D: Applied Physics* (IOP
`iopart` class, per the outline's reference-paper mix).

## Structure
- `main.tex` — preamble, title/abstract, `\input`s the section files below.
- `sections/intro.tex` — literature review + outline reconciliation (done).
- `sections/model.tex` — governing equations, scaling analysis, sheath
  collisionality (done).
- `sections/methods.tex` — unified model + PlasmeX numerical implementation
  (done).
- `sections/results.tex`, `sections/conclusion.tex` — **placeholders only**.
- `references.bib` — all citations used so far; two entries flagged `NOTE`
  need confirmation before submission.

## Compiling
Requires **`iopart.cls`** and **`iopart-num.bst`**, which IOP distributes
directly from the journal's own LaTeX author-guidelines/template page
(not on CTAN — confirmed while preparing this draft; CTAN only carries the
`iopart-num` bibliography-style companion, not the class itself). Place
both files next to `main.tex`, then:
```
latexmk -pdf main.tex
```
If you'd rather not chase down `iopart.cls` yet, the fastest way to proof
the content is to swap `\documentclass[12pt]{iopart}` for
`\documentclass[12pt]{article}` and comment out `\submitto{\jpd}`,
`\ead{...}` and the `iopart-num` bibliography style line (use `plain`
instead) — the body text/equations do not depend on the class.

## Provenance / how to keep this maintainable
The reference PDFs and PlasmeX's internal `doc/*.tex` derivations were read
once and condensed into `../notes/*.md` (not this folder) specifically so
future updates don't require re-reading the source PDFs:
- `../notes/refs_bulk_lte_models.md`
- `../notes/refs_sheath_electrode_models.md`
- `../notes/plasmex_model_notes.md`
- `../notes/outline_review.md` — the outline-vs-literature consistency
  check requested when this draft was created; re-read this first if the
  outline changes, since it explains several non-obvious mappings (e.g.
  why `Lowke2006` is cited under two different modelling categories, and
  why the outline's ~0.2 mm layer and PlasmeX's own ~165 nm sheath width
  are nested, not conflicting, numbers).

When the user updates `../Outline.pdf` or adds papers to
`../Reference Papers/`, update the relevant `../notes/*.md` file(s) first
(new/changed papers only, don't redo the whole batch), then revise the
`.tex` sections against the updated notes rather than re-deriving from the
PDFs directly.
