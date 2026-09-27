# Responses

*Started afresh 2026-09-27. The previous file, Questions 1–11, is archived verbatim in `saved_FAQ_2026-09-27.md`. `literature_review.md` §7.4 still cites "response.md Question 8/9"; those now live in the archive.*

---

## Question 1

How can I validate the **anode** region the same way I validate the cathode region against Fig. 5(b) of Almeida *et al.* (2008)? Search the literature (see also `literature_review.md`) and suggest any kind of validation that works for the anode region.

## Response

### 1.1 What the Fig. 5(b) validation actually is

`validation_fig/almeida_fig5/compare5b.py` is a **code-to-code benchmark** built from these parts:

1. The **same physics problem**: Ar, $p_0=1$ bar, a fixed $T_w=3500$ K, $j_c=10^7$ A m$^{-2}$, local LTE at the outer edge.
2. **Published profiles** $n_e,n_i,n_S,T_e,T_h,j_i/j$ plotted against the distance from the wall $d$ on a log axis. They are extracted from the vector paths of the IOP PDF (`fig5_extract.py`, `pdftrace.py`).
3. **One scalar check** from the paper's table: $j_i/j$ at the wall, Table 1.
4. An overlay of the PlasmeX plotfile, with $d=x_{hi}-x$.

The anode version needs the same four parts. The literature has **one paper that fits this almost exactly** (Semenov 2016). There are also two more code-to-code checks and a set of experimental checks.

---

### 1.2 Primary benchmark (the direct analogue of Fig. 5b): Semenov, Krivtsun & Reisgen (2016)

I. L. Semenov, I. V. Krivtsun, U. Reisgen, *J. Phys. D* **49** 105204 (2016), "Numerical study of the anode boundary layer in atmospheric pressure arc discharges".

It is a unified model: no thermal equilibrium, no ionization equilibrium and no quasi-neutrality are assumed, and the sheath is resolved. It runs for argon at the same $p$ as Almeida, so it is the anode counterpart of Almeida 2008.

| item | Semenov 2016 | current `PlasmeX/exe/anode_inputs` | action |
|---|---|---|---|
| gas, pressure | Ar, 1 atm | Ar, `prob.pres = 101325` | none |
| anode surface temperature | $T_s=1$ kK, fixed | `prob.wall_temp = 1000` | none |
| emission | cold anode, none | `prob.electrode_emission = 0` | none |
| current density | 500–2000 A cm$^{-2}$; **Fig. 3 profiles at 700 and 1500 A cm$^{-2}$**; Table 3 at 500, 1000, 2000 | `prob.arc_Japp = -1e7` (= 1000 A cm$^{-2}$, i.e. Table 3's middle row) | add runs at $7\times10^6$ and $1.5\times10^7$ A m$^{-2}$ for the profiles |
| domain length | ~1 mm; the layer thickness $L_a$ is **1.46 / 0.84 / 0.40 mm** at 500 / 1000 / 2000 A cm$^{-2}$ | `geometry.prob_hi = 2e-4` (**200 µm**) | **must increase** to ≳ 3 mm (see below) |
| outer (plasma-side) BC | ionization equilibrium, quasi-neutral, uniform; bulk $T$ = 13–18 kK depending on $j$ | Dirichlet: `arc_Te = 12000`, `arc_e = 6.46e22`, `arc_efldx = -700` | set it from the bulk state for each $j$ (below) |
| wall BCs | electron flux (their eq. 25), an ion-mobility-based field condition (eq. 28), an electron-energy condition (eq. 30) | `e_bc_hi = 3`, `ion_bc_hi = 3`, `een_bc_hi = 3` | read eqs. 25–30 and match the ion BC, as you did with `ion_bc_hi = 1` for Almeida |

**What to compare** (Semenov's figure and table numbers are as given on the IOP abstract page; check them against the PDF):

- **Fig. 3:** $n_e$, $n_i$, $n_{\rm Saha}$, $T_e$ and $T_h$ against the distance from the anode, at 700 and 1500 A cm$^{-2}$. This is the counterpart of Almeida's Fig. 5 panels (a) and (b).
- **Fig. 4:** the electron and ion current densities against distance. This is the counterpart of the $j_i/j$ panel, although at the anode $j_i/j$ is small and has the opposite sign to the cathode case.
- **Table 3** is the counterpart of Almeida's Table 1:

| $j$ (A cm$^{-2}$) | $U_a$ (V) | $U_1$ (V) | $U_{sh}$ (V) | $L_a$ (mm) | $L_{sh}$ (µm) |
|---|---|---|---|---|---|
| 500 | 2.13 | 2.05 | 0.33 | 1.46 | 1.564 |
| 1000 | 2.63 | 2.60 | 0.56 | 0.84 | 1.525 |
| 2000 | 3.18 | 2.35 | 1.00 | 0.40 | 1.420 |

Here $U_a$ is the anode voltage drop, $U_1$ the drop across the ionization layer, $U_{sh}$ the sheath drop, $L_a$ the thickness of the anode layer and $L_{sh}$ the sheath thickness.

**Settings that decide whether the comparison is fair:**

1. **Domain length.** The layer thickness $L_a$ is the distance over which $n_e$ relaxes to $n_{\rm Saha}(T)$. At 1000 A cm$^{-2}$ it is 0.84 mm, and at 500 A cm$^{-2}$ it is 1.46 mm. A 200 µm domain puts the Dirichlet LTE state *inside* the non-equilibrium layer, which clamps $n_e$ and $T_e$ and shifts $U_a$. Use $L\ge 3L_a$, e.g. `geometry.prob_hi = 3e-3` or `5e-3`, and check that $U_a$ does not change when $L$ is doubled. This is the domain-length issue already noted in `literature_review.md` §7.4, but it is now quantified by Semenov's $L_a$. Expect slower convergence: the heavy-gas relaxation time scales as $L^2/\alpha$.

2. **Bulk state.** At the outer edge Semenov imposes an equilibrium, uniform plasma:
$$T_e=T_h=T_\infty,\qquad n_e=n_i=n_{\rm Saha}(T_\infty,p_0),\qquad E_\infty=\frac{j}{\sigma(T_\infty)} ,$$
where $T_\infty$ comes from a column energy balance, and the abstract gives 13–18 kK over their $j$ range. `PlasmeX/exe/almeida_arc_bc.py` already solves the same local problem: ionization balance, Ohm's law, and $jE=w_{\rm rad}$ (Almeida eqs. 26–30). Run it at $j=7\times10^6$, $10^7$ and $1.5\times10^7$ A m$^{-2}$, **check the resulting $T_\infty$ against the values Semenov states** (text or Fig. 3 at $x=L$), and put $T_\infty$, $n_{\rm Saha}$ and $E_\infty$ into `arc_Te`, `arc_e`/`arc_arp` and `arc_efldx`. If Semenov's $T_\infty$ differs from the one the script gives, use Semenov's value, because the aim is the same problem, not the same column model.

3. **The definition of $U_a$.** The potential difference across the whole domain includes the ohmic drop of the equilibrium plasma. Subtract it:
$$U_a=\big[\varphi(x_{\rm edge})-\varphi(x_{hi})\big]-E_\infty\,(x_{hi}-x_{\rm edge}) ,$$
with $x_{\rm edge}$ in the uniform region. Semenov quotes $U_a>0$ for what is physically a *negative* anode fall (the anode sits below the adjacent plasma). Check which sign convention they use before comparing, and do the same for $U_{sh}$ and $L_a$: use their definitions of the sheath edge and the layer edge, e.g. a threshold on $n_i-n_e$ and on $n_e/n_{\rm Saha}$.

4. **Resolution.** $L_{sh}\approx1.5$ µm, which your `tag_electrode_depth` (1, 0.5, 0.2 µm) already resolves.

**Reusing the tools.** `fig5_compare.py` already handles "electrode at $x_{hi}$, plot against $d=x_{hi}-x$". Write a thin `compare_semenov3.py` alongside it with a `CASE` dict holding $j$, the axis limits and the Table 3 row. If the IOP PDF of Semenov 2016 contains **vector** curves, as Almeida's does, `pdftrace.py`/`fig5_extract.py` can pull them out directly. If they are raster, digitize them with WebPlotDigitizer and save the result in the same CSV format (`quantity,unit,x_m,value`).

**What counts as agreement.** Transport data differ between the codes: Semenov's are not the M&K/LXCat set PlasmeX uses. Aim for the same pattern as your cathode comparison: profiles on top of each other in the quasi-neutral part, the location of the $n_e/n_{\rm Saha}$ departure within ~20 %, $U_a$ within ~10 %, and the correct **trends with $j$**: $|U_a|$ increases, $L_a$ shrinks and $U_{sh}$ grows. The trend test is the most robust part, because it does not depend on small transport differences.

---

### 1.3 Second code-to-code check: Khrabry, Kaganovich, Nemchinsky & Khodak (2018), open access

A. Khrabry *et al.*, *Phys. Plasmas* **25** 013521 (2018), "Investigation of the short argon arc with hot anode. Part I", open access as [arXiv:1708.01301](https://arxiv.org/abs/1708.01301). I read the full text.

This paper is useful for three reasons.

- **It closes a triangle.** Their Fig. 2 benchmarks their code against **Almeida 2008 at the cathode**: Ar, 1 atm, $T_w=3500$ K, $j=5\times10^6$ A m$^{-2}$, and they write "the profiles … almost coincide". PlasmeX is already benchmarked against Almeida at the cathode. So PlasmeX and Khrabry agree through a common reference, and the Khrabry anode results become a fair target.
- **Fig. 3: near-anode profiles** at the same $p$ and $j=5\times10^6$ A m$^{-2}$, with a **hot, emitting anode**, $T_w=3900$ K. The potential is non-monotonic (a double layer): a negative field in the quasi-neutral region and a positive sheath drop that suppresses emission. This is the only published case I found that tests PlasmeX's anode **with `electrode_emission = 1`**. The anode case normally runs with emission off, so this is an extra test, not a replacement.
- **Figs. 13–14: cold-anode fall.** They give the sheath voltage and the voltage of the non-equilibrium layer as functions of $j$ over about $2\times10^6$–$2\times10^7$ A m$^{-2}$, for $T_w=1000$ K and $3000$ K, and state that the result hardly depends on $T_w$. For a cold anode, both voltages are negative and grow in magnitude with $j$ (a falling VAC of the anode layer). **This is the anode analogue of the cathode CVC.** A $j$-sweep of `anode_inputs` at $T_w=1000$ K (and one repeat at 3000 K) gives the comparison curve directly.

Caveats: their model covers the whole gap, with radiation and radial losses, so the bulk-edge state is not imposed as in PlasmeX. Take the plasma state at the edge of their anode layer from their Figs. 4–5 and use it as your Dirichlet state. Their figures are **raster JPEGs** (I checked with `pdfimages`), so the vector extractor will not work on them; use WebPlotDigitizer.

---

### 1.4 Scalar check that needs no profiles: anode heating voltage (Almeida, Cunha & Benilov 2017)

N. A. Almeida, M. D. Cunha, M. S. Benilov, *J. Phys. D* **50** 385203 (2017). This is the same group and the same unified model as the Fig. 5 benchmark. Its result is simple and very testable:

$$q_a \simeq j\,U_h,\qquad U_h=6.1\ \text{V for Ar at 1 bar}\ (5.9\text{ V at 10 bar},\ 6.7\text{ V at 100 bar}),$$

and $q_a$ is only weakly affected by $T_w$ (their Fig. 2 and Table 1).

**How to use it:** run `anode_inputs` over a small matrix of $j$ ($2\times10^6$ to $2\times10^7$ A m$^{-2}$) and $T_w$ (1000, 2000, 3000 K). Compute $q_w$ from eq. (`eq:bc-qwall`) and plot $q_w/j$ against $j$. Pass criteria: **(i)** $q_w$ is linear in $j$, **(ii)** $q_w/j\approx6$ V, **(iii)** $q_w/j$ hardly depends on $T_w$.

The comparison is like for like. Your $q_w$ contains the work-function term $W\cdot\tfrac14 n_e\bar c_e$, and so does their $q$: most of the 6.1 V is $W_{\rm W}\approx4.5$ V plus the electron enthalpy term $2kT_e/e$. At the anode, emission is negligible for $T_w\lesssim3000$ K. Above that the "emission off" setting no longer matches their model.

Almeida 2017 then uses $U_h$ in a thin-rod anode model and reports "good agreement" with several lamp experiments. So agreement on $U_h$ also ties PlasmeX, indirectly, to those experiments.

---

### 1.5 Further code-to-code option

**Santos, Almeida, Lisnyak, Gonnet & Benilov 2022**, *Phys. Plasmas* **29** 043503: the passive-anode regime (low/moderate $T_w$, a positive space-charge sheath with a negative drop) against the active one. In the passive regime, their near-anode profiles and $U_a$ at 1 atm Ar are a third independent check. The paper is behind a paywall, and I have not seen its figures, so check which $j$ and $T_w$ it shows before you plan runs.

---

### 1.6 Experimental validation of the anode region

No single anode experiment is a clean 1D test. The attachment is multidimensional, and the probe or spectroscopic resolution (~50–100 µm) is coarser than the sheath. So use the experiments for **integral quantities** (fall voltage, heat flux or power input, $T_e$ and $n_e$ at the edge of the layer), after setting the PlasmeX bulk state to the measured one.

| data | conditions | what to compare | how |
|---|---|---|---|
| **Sanders & Pfender 1984**, *J. Appl. Phys.* **55** 714 | free-burning Ar, 1 atm, water-cooled Cu anode, 100–250 A; probes through the anode | anode fall **−2.1 to −1.4 V**; anode heat flux | take $j$ at the axis and the bulk $T$ from their measurements; set $T_w\approx$ Cu surface ($\lesssim$ 1000 K); check the sign and magnitude of $U_a$ and $q_w$ |
| **Sanders, Etemadi, Hsu & Pfender 1982**, *J. Appl. Phys.* **53** 4136 | Ar, 1 atm, 100–300 A; spectroscopy of the whole anode region | $T_e$, $n_e$ against the distance from the anode; departure from LTE near the anode | overlay the PlasmeX $T_e(d)$, $n_e(d)$, averaged over the instrument resolution, on the measured points in the outer ~0.1–1 mm (this is the region Semenov's $L_a$ says is non-equilibrium) |
| **Redwitz *et al.* 2006**, "Arc attachment at HID anodes: measurements and interpretation", *J. Phys. D* **39** 2160 | W rod anodes in Ar at 0.1–1 MPa (HID model lamp) | anode fall; **$T_e$ and $n_e$ in the anodic boundary layer**, spectroscopic, high spatial resolution; power input ∝ current | the diffuse-attachment cases only (a negative fall in the centre of a diffuse attachment); the power input per ampere is a direct check of $U_h$ in §1.4 |
| **Dabringhausen *et al.* 2005**, *J. Phys. D* **38** 3128; **Nandelstädt *et al.* 2002** I–III | Ar 0.26 MPa, 1.5–6 A, W electrodes | electrically measured anode falls; pyrometric $T_w$ | the anode fall against $I$; $T_w$ as input |
| **Heberlein, Mentel & Pfender 2010**, *J. Phys. D* **43** 023001 (survey) | compilation | typical anode falls, boundary-layer $n_e$, $T_e$ and heat fluxes ($>10^7$ W m$^{-2}$) | a "reasonable range" plot: PlasmeX points inside the band of the measured values |

Also worth citing, as non-sheath-resolved models with experimental comparison: Amakawa, Jenista, Heberlein & Pfender 1998 (*J. Phys. D* **31** 2826; 200 A transferred Ar arc; negative anode fall and heat flux against the flow rate) and Yang & Heberlein 2007 (anode attachment modes, *PSST* **16** 529).

---

### 1.7 Recommended anode validation set for the paper

| # | run | reference | figure or table in the paper |
|---|---|---|---|
| A1 | `anode_inputs`, $L$ = 3–5 mm, $T_w=1000$ K, $j=7\times10^6$ and $1.5\times10^7$ A m$^{-2}$, bulk state from §1.2 (2) | **Semenov 2016 Fig. 3 (and 4)** | a 3-panel overlay like `compare5b`: $n_e,n_i,n_S$ / $T_e,T_h$ / $j_e,j_i$. **This is the anode Fig. 5b.** |
| A2 | the same at 500, 1000, 2000 A cm$^{-2}$ | **Semenov Table 3** | a table of $U_a$, $U_{sh}$, $L_a$, $L_{sh}$: PlasmeX against Semenov |
| A3 | $j$ sweep $2\times10^6$–$2\times10^7$, $T_w$ = 1000 and 3000 K | **Khrabry 2018 Fig. 14** | $U_{sh}(j)$ and $U_{\rm layer}(j)$, the anode "CVC" |
| A4 | the A3 runs, post-processed | **Almeida 2017** ($U_h$ = 6.1 V) | $q_w/j$ against $j$ for three $T_w$ |
| A5 (optional) | $T_w=3900$ K, `electrode_emission = 1`, $j=5\times10^6$ | Khrabry 2018 Fig. 3 | double-layer potential profile |
| A6 | bulk state from the measured $T$; Cu anode at $T_w\approx$ 500–1000 K | Sanders & Pfender 1984; Sanders *et al.* 1982; Redwitz 2006 | anode fall and $T_e(d)$ against experiment |
| A0 | before all of the above: $L$ = 1, 2, 4 mm and a halved $\Delta x$ at the wall | – | a convergence note: $U_a$ and $q_w$ unchanged |

A1 + A2 replace Fig. 5b one-to-one. A3 + A4 are cheap because they reuse the same runs, and they give two independent code-to-code checks from two other groups. A6 is the experimental anchor. Together they give the anode section the same weight as the cathode section, a code benchmark plus experiment, without needing any new physics in the code.

**Before running:** get the Semenov 2016 PDF and note (a) their definitions of $U_a$, $U_{sh}$ and $L_a$, (b) the bulk $T_\infty$ at each $j$, (c) their wall BCs, eqs. 25–30, and (d) whether the Fig. 3 curves are vector paths. Items (a)–(c) fix the deck; item (d) decides whether `fig5_extract.py` can be reused.

---

### Sources (Question 1)

- [Semenov, Krivtsun & Reisgen 2016, J. Phys. D 49 105204](https://iopscience.iop.org/article/10.1088/0022-3727/49/10/105204): abstract, Fig. 3/4 content, Table 3
- [Khrabry et al. 2018, Part I (arXiv 1708.01301, full text)](https://arxiv.org/abs/1708.01301)
- [Almeida, Cunha & Benilov 2017, J. Phys. D 50 385203](https://iopscience.iop.org/article/10.1088/1361-6463/aa7f97): $U_h$ = 6.1 V (Ar, 1 bar)
- [Santos et al. 2022, Phys. Plasmas 29 043503](https://pubs.aip.org/aip/pop/article-abstract/29/4/043503/2843564/Numerical-investigation-of-regimes-of-current)
- [Sanders & Pfender 1984, J. Appl. Phys. 55 714](https://pubs.aip.org/aip/jap/article-abstract/55/3/714/168524/Measurement-of-anode-falls-and-anode-heat-transfer)
- [Sanders, Etemadi, Hsu & Pfender 1982, J. Appl. Phys. 53 4136 (related search result)](https://link.springer.com/article/10.1007/BF00566213)
- [Arc attachment at HID anodes: measurements and interpretation (Redwitz et al. 2006)](https://www.osti.gov/etdeweb/biblio/20826050)
- [Dabringhausen et al. 2005, Different modes of arc attachment at HID cathodes](https://www.osti.gov/etdeweb/biblio/20665787)
- [Determination of HID electrode falls in a model lamp II: Langmuir-probe measurements](https://www.osti.gov/etdeweb/biblio/20310618)
- [Heberlein, Mentel & Pfender 2010, The anode region of electric arcs: a survey](https://iopscience.iop.org/article/10.1088/0022-3727/43/2/023001)
- [Amakawa et al. 1998, Anode-boundary-layer behaviour in a transferred high-intensity arc](https://iopscience.iop.org/article/10.1088/0022-3727/31/20/017)
- [Anode attachment modes and their formation in a high intensity argon arc (Yang & Heberlein 2007)](https://www.osti.gov/etdeweb/biblio/20906408)
- [Almeida et al. 2009, near-anode layers of very high-pressure arcs (Hg, Xe: not an Ar benchmark)](https://ui.adsabs.harvard.edu/abs/2009JPhD...42d5210A/abstract)

*Checked only at the abstract or search level:* Semenov's BC equation numbers and bulk temperatures; the figure contents of Santos 2022; the journal volume and pages of Redwitz 2006 and Sanders 1982 (confirm them before adding them to `references.bib`). I looked at Almeida's PhD thesis (2009): its anode chapter is Hg/Xe at 100 bar, so it gives no Ar anode profiles.
