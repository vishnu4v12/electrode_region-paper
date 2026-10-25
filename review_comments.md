# Referee report

**Manuscript:** "Modelling near electrode regions of an atmospheric-pressure arc discharge"
**Target journal:** J. Phys. D: Appl. Phys.
**Version reviewed:** working tree of 22 Oct 2026 (`main.tex`, `sections/*.tex` as compiled; `sections/methods.tex` is not `\input` and was not reviewed as part of the paper)
**Recommendation:** major revision

---

## 1. Summary

The manuscript presents a one-dimensional, transient, unified drift-diffusion–Poisson model of the near-cathode and near-anode regions of an atmospheric-pressure argon arc. It has three species, separate electron and heavy-particle energy equations, Stefan–Maxwell-type drag and thermal forces, a heavy-gas momentum equation, a volumetric radial heat loss, and (at the cathode) a surface heat balance with Schottky-enhanced thermionic emission. The outer boundary is placed at L = 0.5 mm, where an LTE state is imposed. The model is compared with the unified solutions of Almeida et al. (2008, cathode) and Semenov et al. (2016, anode). It is then applied to a cathode case (21 kK, 1e8 A/m²) and an anode case (15 kK, 4e6 A/m²). The results cover the layer structure, the electron energy budget, the momentum balance and collisionality of the sheaths, the electrode heat flux, an anode parameter sweep, and a linear-response estimate of when the layers can be treated as quasi-steady.

**Strengths.**
- The model is carefully formulated. The species momentum balances are written with explicit Newton's-third-law pairs (eq. newton3), and the electron transport is mapped consistently onto Mitchner & Kruger's coefficients (app. electrons).
- The validation against two independent published unified models agrees well on the integral quantities (j_i/j at the wall, sheath voltage, T_e at the surface).
- The a posteriori checks in sec:res-momentum are uncommon and useful: the neglected inertia evaluated from the solution, the mean free paths across the sheaths, the Boltzmann test, and λ_{T_h}/L_{T_h}.
- The quasi-steadiness analysis (eqs g → qs-error) is a clean and reusable way to turn a step response into a criterion for coupling to an LTE model.
- The anode sweep has a clear physical message: T_bulk controls j_i and gas conduction, while J controls the electron terms.

**Main concerns, in brief.** (i) The paper does not show that the results are independent of the domain length L, the radial-loss parameters (R, T_edge) or the mesh. (ii) The numerical method is under-described for a paper whose transient results reach the ps scale. (iii) Some headline conclusions (Bohm overestimate, "drift-diffusion adequate for electrons", cathode q_w) are more parameter- or definition-dependent than the text suggests. (iv) The cathode validation is not like-for-like in geometry.

---

## 2. Major comments

### M1. Dependence on the domain length L (and on what is imposed at x = L)
The whole study rests on imposing an LTE state at x = L = 0.5 mm. The paper's own results suggest the layers are not complete within that distance:
- **Anode.** n_e exceeds n_eq(T_h) by ×2.4 at 200 µm and ×1.2 at 300 µm, and "comes within 10 % of n_eq(T_h) only within 0.15 mm of x = L" (results §res-structure). T_e − T_h is still 1.3 kK at 100 µm. The base_profiles panels show finite T and n gradients at x = L. So the anode layer occupies essentially the whole domain, and the Dirichlet LTE state is imposed where the plasma is only just reaching equilibrium.
- **Cathode.** 47 % of the domain's Joule heat leaves through x = L as electron energy flux (§res-energy). Whatever receives this flux is not modelled, and T_h is still rising at x = L in fig. base-profiles(c).
- Table tab:scales gives λ_rxn ≈ λ_T ≈ 1.4–3.4 × 10⁻⁴ m at 11–14 kK, which is the same order as L. Section 2.1 itself says "the electrode affects the plasma over a region of thickness of order 10⁻⁴ m".

**Request:** repeat both base cases (at least) with L = 1 mm and 2 mm (same bulk T, J, R), and report the changes in q_w, j_i/j, φ_s, x_s and T_e,w. If the anode quantities move by more than a few per cent, either choose a longer L or state clearly that the results are conditional on L = 0.5 mm. In that case, explain how an LTE model would supply the gradients, not only the values, at x = L. The paper says gradients "may exist at this boundary", but with a Dirichlet condition they are outputs of the model, not inputs (see also the j_i/j-maximum discussion in §val-anode, which makes exactly this point).

### M2. Sensitivity to the radial-loss closure (R, T_edge, form of Q)
The radial loss (eq. radial-loss) is a modelling novelty claimed in the introduction, but its parameters are assumed:
- R = 0.5 mm (cathode) and 3.5 mm (anode) are "typical", and T_edge = 1000 K is not justified.
- With uniform κ, Λκ(T − T_edge) is the Elenbaas–Heller estimate. κ varies by more than an order of magnitude between the core and an edge at 1000 K, so evaluating κ at the local (core) temperature overestimates the loss. A Kirchhoff-potential form, Λ[S(T) − S(T_edge)] with S = ∫κ dT, is the standard fix and costs nothing.
- κ_R is an LTE quantity (Mutation++). It is applied in the electron channel, including in the non-equilibrium layer. The paper argues that the loss is negligible near the electrode, which is plausible, but it should be shown (for example, the share of the Q_e integral within 100 µm).
- Consistency with the stated arc current: J·πR² gives 79 A at the cathode (1e8 A/m², R = 0.5 mm) and 154 A at the anode (4e6 A/m², R = 3.5 mm). Both values are quoted as "typical of the axis of a 200 A arc" (results ¶1). Either the radii or the current densities do not correspond to 200 A, or the text should explain why they need not.

**Request:** a short sensitivity table (R ×0.5/×2, T_edge 300/1000/3000 K) for the electrode quantities of tab:base. Without it, the reader cannot tell which conclusions come from the near-electrode physics and which come from the closure.

### M3. Mesh convergence is not demonstrated
Section numerics gives the mesh (944 cells, finest 3.8 nm), and §validation states that the cells are smaller than λ_De by a factor of ≥ 1.4 (cathode) and 6 (anode). No convergence evidence is given. A factor of 1.4 in the cathode sheath is modest, and the electrode quantities are wall-cell values (n_e,w, T_e,w, E_w).

**Request:** a mesh-refinement table (finest Δx halved and doubled, or one level added and one removed) for the quantities of tab:base, with the observed order. One paragraph and one table in §numerics would suffice.

### M4. The time integration is not described, yet transient results are central
§numerics says only that "all equations … are solved together at each time step". There is nothing on the time discretisation, the nonlinear solver, the time-step control or the accuracy tolerance. But §res-time reports a q_w dip of 32 % "within 10 ps", D coefficients to two significant figures, and t_s values down to 0.13 µs.

**Request:** state the integrator (implicit scheme, order, Newton, Jacobian treatment), the error tolerance and the resulting Δt range during a step response. Also state any lagged terms (e.g., if the drag velocities or the Wannier drift are taken from the previous step). Show that D_X is converged with respect to the time tolerance, for example D at two tolerances for one table case. Also describe how the ε = 10⁻⁴ heavy-inertia scaling is implemented, and confirm that it does not affect the non-acoustic dynamics: the 0.3 % statement covers only the steady state.

### M5. Validation against Almeida et al.: the geometry difference is large at x = L
Almeida et al. solve in spherical symmetry about a hemisphere with r_c = 1 mm, so their j(x) = j_c r_c²/(r_c + x)². At x = L = 0.5 mm this is 0.44 j_c, and their Joule heating there is about 0.2 of the planar value. The present model keeps j = j_c throughout but imposes their state at x = L. The agreement is good, but the reader cannot tell whether it is robust or partly fortuitous. The cathode-sheath discrepancies (densities ×2.5 at 0.1 µm, E 2.3 vs 1.6 × 10⁷ V/m, x_s 1.1 vs 1.5 µm) are attributed only generally to "transport properties and ion boundary condition".

**Requests:**
- Either include the area factor in the 1-D equations (a quasi-1-D (r_c + x)² metric) for the validation case, or quantify its effect.
- Identify which model difference produces the sheath discrepancy, e.g., by switching the ion mobility to field-independent or the ion BC to Almeida's form. The near-wall field enters the Schottky lowering, so a 44 % field difference is not cosmetic.
- State which rate coefficients and transport were used for the Semenov case. The paper uses Almeida's rates throughout. If Semenov et al. used others, "reproduced as closely as possible" should say so.

### M6. The claim "drift-diffusion is adequate for the electrons in the sheaths"
The sheaths are 0.18 and 0.30 electron mean free paths thick. The Boltzmann test of fig. base-sheath(c,d) compares n_e with n_B built from the model's own φ and T_e. When the drag is small, a fluid balance −en E − ∇p_e ≈ 0 is Boltzmann by construction. So the test shows that the drag is small, not that the fluid description is valid. Two kinetic effects are not addressed:
1. **Absorbing-wall tail depletion.** In a collisionless retarding sheath, the returning half of the distribution is missing near the wall. The wall flux is then not ¼ n_e c̄_e with the full local density: the Hagelaar-type form has ½ n c̄ for the same reason. This affects n_e,w, the absorbed electron flux (408 MW/m² at the cathode) and the anode q_w, which is 66 % electronic.
2. **The emitted-electron beam at the cathode.** Emitted electrons (1.65 J) are accelerated through 3.7 V over ≈ 0.2 λ_e, so they form a beam, not part of a Maxwellian at T_e. The model deposits the sheath "Joule heating" (1.8 × 10¹⁵ W/m³ at the wall) into the local thermal population. The mean velocity u_e = Γ_e/n_e also cannot represent counter-streaming emitted and plasma electrons. For the same reason the "electron inertia exceeds drag by ×44" result is not physically meaningful.

**Request:** temper the conclusion ("remains adequate") to what is shown. Discuss the wall-flux factor and its effect on the electrode quantities: a run with the ½ n c̄ form would bound it. Argue, with λ_e versus the energy-deposition length, that the beam energy is thermalised within the ionization layer (as in Benilov-type models), so that the integral energy balance is unaffected even if the local profile is.

### M7. Bohm-criterion comparison depends on the sheath-edge definition
The sheath edge is the 1 % charge-separation point. At the cathode the ions reach u_B at 4.2 % charge separation, inside the "sheath". So the conclusion that "a Bohm flux would overestimate the ion flux by 35 %" (abstract, conclusion) depends directly on the 1 % threshold. In addition, the relevant criterion for a collisional presheath/sheath is the modified Bohm criterion (Godyak & Sternberg 1990; Riemann 1991; Valentini 1996), not the collisionless u_i ≥ u_B. That the collisionless criterion fails in a sheath of 1.9–9.4 ion mean free paths is expected.

**Request:**
- Give the ion flux ratio for two or three edge thresholds (e.g., 0.1, 1, 5 %), or define the edge where u_i = u_B and report the charge separation there.
- Compare with a collisional sheath criterion or with a matched collisional-sheath model.
- Rephrase the abstract and conclusion accordingly. The anode factor of 5 (u_i = 0.2 u_B) is robust. The cathode 35 % is not.

### M8. Neglected ion inertia at the cathode, and the ion impact energy
§res-momentum shows that the neglected ion inertia exceeds the drag by ×2 at the cathode surface and is above 10 % of the drag out to 1.9 µm. The text concludes that the local field–velocity relation "overestimates the ion velocity near the surface". The quoted ion impact energy, 5.1 eV from v_d = 5.0 km/s, is evaluated exactly there. The ion kinetic-energy flux it implies, 5.1 eV × 3.3 × 10²⁵ m⁻² s⁻¹ ≈ 27 MW/m², equals the computed cathode q_w (26.6 MW/m²). The paper says only that eq. bc-qwall assumes the ion kinetic energy is "returned to the neutrals".

**Request:**
- State the size of this omitted flux next to q_w.
- Explain whether it is double-counted with the heavy-gas conduction term (in the model the ion energy gained from the field is dissipated through Q_in into T_h and returns partly as κ_h ∇T_h at the wall).
- Bound the effect of ion inertia, at least on φ_s, E_w and q_w, for example by including the ion momentum equation in a sensitivity run.

### M9. The cathode heat flux and wall temperature are largely set by the assumed electrode model
q_w is 26.6 MW/m² out of 566 MW/m² heating and 540 MW/m² emission cooling. The text rightly says T_wall is fixed mainly by the emission balance. It follows that q_w ≈ (k_w/L_w)(T_wall − T₀) is essentially determined by the assumed k_w = 100 W/m/K, L_w = 1 cm and T₀ = 1000 K. Also:
- T_wall = 3669 K is within ~25 K of the melting point of tungsten (3695 K). This deserves a comment, because it suggests that the cathode conduction model or the operating point is at the edge of physical plausibility.
- The quasi-steady analysis holds T_wall fixed during the steps. For AC operation (§res-time, last paragraph), the cathode surface temperature has its own thermal time (ms–s), comparable to or longer than the 0.4 ms quoted. Because emission is exponential in T_wall, the AC conclusion should either include the wall thermal response or be restricted to the plasma side.

**Request:** report ∂q_w/∂(k_w/L_w) (or q_w for two L_w values), and qualify the conclusion bullet "set mainly by the electrode rather than by the plasma" to say it is set by the assumed conduction model.

### M10. Scope of the quasi-steadiness conclusions
- The response kernels are linearised about one operating point per electrode (plus the anode sweep). The AC statement, "quasi-steady except within about 0.4 ms of zero current", applies the cathode D_J computed at 1e8 A/m² down to currents near zero, where the cathode layer (emission, sheath, T_wall) is completely different. There is no cathode current sweep to support this.
- After a T_bulk step, the bulk densities jump instantly to n_eq(T_bulk). In an LTE model they would follow the bulk T on the reaction time. This is fine as a modelling choice, but it should be stated as one and its effect on D_T noted.
- The acoustic suppression is justified by "relieved radially in a real arc". A reference or a one-line estimate (radial acoustic time R/c compared with the axial one) would help.

---

## 3. Minor comments (technical)

1. **n_eq(T_h) as the reference density** (results §res-structure, fig. base-profiles). k_rec follows from k_ion by detailed balance at T_e (§model-governing), so the state the reactions relax towards is Saha at T_e (two-temperature), not n_eq(T_h). The text calls n_eq(T_h) "the equilibrium density that the reactions would give at the local heavy-particle temperature", which is misleading. Define n_eq precisely: which Saha, which pressure law. Either plot the T_e-based equilibrium as well or rephrase. The cathode panel also shows n_eq(T_h) > n_e around 50–100 µm (an ionization deficit relative to it), which the text does not mention.
2. **Scale hierarchy** (§2.1). "λ_rxn ≈ λ_T" holds only at 11–14 kK. At 18–23 kK, λ_T/λ_rxn is 2.6–8.4 (tab:scales). Write λ_rxn ≲ λ_T, or say where they are comparable.
3. **Radial loss and the introduction.** The introduction says the radial loss "can be important, especially in the near-cathode region". §model says it is "comparatively negligible" close to the electrodes, and the cathode budget puts 30 % of the Joule heat into the radial loss, "most of it at x > 100 µm". Make the introduction consistent: it matters in the column part of the domain, not in the layer.
4. **Electron wall energy.** Eq. bc-electron-energy uses 2k_BT_e per absorbed electron. This is consistent with the half-Maxwellian assumption in eq. bc-electron, but if M6 leads to a different wall-flux factor, the energy flux factor should change with it.
5. **Ion wall flux** (eq. bc-ion). The effusion term uses c̄_i at T_wall but the local n_i. Say whether n_i is the wall-cell value. The same factor-of-2 question as M6.1 applies.
6. **Heavy-gas temperature jump.** T_h(0) = T_wall is imposed with no accommodation/Knudsen jump. With λ_{T_h}/L_{T_h} up to 0.10 at the cathode surface, and gas conduction being 31 % of the anode q_w and 88 MW/m² at the cathode, a one-line estimate of the jump would be appropriate.
7. **Wannier mobility at T_h vs T_eff** (app. ions). The closure is 15–26 % below the Monte Carlo mobility at 300–3000 Td. E/N at the cathode surface is roughly 10³ Td (E = 1.9 × 10⁷ V/m, n ≈ 2 × 10²⁴ m⁻³). State the E/N range in the sheaths, and the expected effect of a ~20 % mobility error on j_i and φ_s.
8. **The mobility comparison is at 300 K only.** The application has T_h from 1 to 21 kK. Say whether eq. Qia has been checked against data at elevated T.
9. **κ_h, κ_i "as in Almeida et al."** Give the formulas or the equation numbers in Almeida 2008 so the paper is self-contained.
10. **LTE conductivity comparison** (app. lte-check). "Overestimates the conductivity at intermediate temperatures": quantify it (maximum % and at which T). Say whether 13–17 kK (the anode sweep) is affected, since σ sets the column field and the Joule term.
11. **Notation clashes.**
    - λ is used for mean free paths (λ_e, λ_i), lengths (λ_T, λ_rxn) and thermal conductivities (λ, λ_en, λ_ei in app. electrons).
    - Λ is both the Coulomb logarithm argument and the radial-loss factor 4/R².
    - C is the thermal-force coefficient (C_e, C_i, C_n) and also the electron speed (C, C̄_e).
    - x is both position and reduced energy (eq. Qbar-e, eq. nu-M).
    - T_w, T_s and T_wall all appear for the surface temperature (validation captions).
12. **"superscript m"** (model.tex:85) should be M.
13. **"944 cells in total"** (§numerics). Say "leaf cells". The full hierarchy has more. Also give the base-level cell size (977 nm), since only the refined sizes are listed.
14. **Steady-state criterion.** Give a physical residual of the converged states, e.g., the variation of the total current density across the domain and the closure of the global energy balance. This complements the ‖ΔU‖/‖U‖ criterion.
15. **Table tab:base.** Add the effective anode voltage q_w/J (8.9 V) and the cathode equivalent. These are the quantities usually compared with measurements, e.g., anode heat flux measurements (Nestor 1962; Sanders & Pfender 1984, already cited) and measured cathode falls. A short comparison with experiment, even order of magnitude, would strengthen the results section.
16. **Statement of novelty.** The introduction would benefit from an explicit list of what is new relative to Almeida (2008), Semenov (2016) and Krivtsun (2025): the near-LTE outer boundary at 0.5 mm, the radial loss, the heavy-gas momentum equation, the a posteriori collisionality analysis, the anode sweep, and the quasi-steadiness criterion. Each contribution should then be tied to a section.
17. **"In the steady states the heavy gas is at rest"** (§res-momentum). In 1-D with an impermeable wall this is exact (zero total heavy flux), so say so instead of calling it negligible.
18. **Section structure.** §2 is titled "Physical model and scale analysis" but also contains the boundary conditions and the numerical implementation, both as subsections. Either retitle §2 or make §numerics its own section.
19. **The abstract** says the bulk state is "supplied by measurements or by an LTE model". In this paper it is assumed typical values (results ¶1). Rephrase to "can be supplied".

---

## 4. Language and typesetting

| Location | Issue |
|---|---|
| `main.tex` title | "near electrode regions" → "near-electrode regions" |
| `intro.tex:5` | "Resolving these wide range of length scales" → "this wide range" |
| `model.tex:7` | "perturbed by the presence of relatively cold and absorbing nature of electrode surfaces" (ungrammatical) |
| `model.tex:46` | "thier" → "their"; "in the cathode region even though they are present but are shown to have" (double conjunction) |
| `model.tex:85` | "the superscript \(m\)" → \(M\) |
| `model.tex:191` | "We therefore solve for following qunatities" → "the following quantities" |
| `numerics.tex:4` | "smaller than either Section (\ref{…})" → "than either (Section~\ref{…})" |
| `validation.tex:6` | "1\,\%,the cells" (missing space) |
| `validation.tex:35` | "For anode we consider" → "For the anode, we consider" |
| `validation.tex:37` | "obatined" → "obtained"; double space before "sheath width" |
| `results.tex:4` | "200 amps" → "200 A"; "Bulk-plasma at the cathode is" → "The bulk plasma at the cathode is"; "At anode however" → "At the anode, however," |
| `results.tex:6` | "The LTE models and experiments … uses" → "use" |
| `results.tex:117` | "( Figure" → "(figure" (lower-case, as elsewhere) |
| `results.tex:158` | "with a increase" → "with an increase" |
| `results.tex:287` | "where as" → "whereas" |
| `conclusion.tex:15` | "Coupled to an LTE model of the bulk plasma, would enable us todetermine" → "Coupling to an LTE model of the bulk plasma would enable us to determine" |
| general | `\label{sec:scales involved}` contains a space. It works, but `sec:scales` is safer. |
| general | Units: write "µm" consistently with `\ \mu`m spacing; use "kK" or "K" consistently in prose. |

---

## 5. Before submission (housekeeping)

- **Author, affiliation and acknowledgements** are placeholders ("Author Name", "Department/Institution", `\ack{To be completed.}`).
- **Data/code availability statement.** IOP requires one. Now that the code name is gone from the text, decide what is said about code availability ("available on request" or a repository link). Referees will ask.
- **`sections/methods.tex`** is not `\input` but is in the source tree. It contains stale and contradictory text (2-D Hypre runs, "100× larger time step", a reference to a non-existent `sec:intro-electrode`, the LiangTrelles conductivity) and still names PlasmeX twice. Delete it or keep it out of the submission bundle.
- **LaTeX comments in `results.tex`** (lines ~12–17, 92, 142–191) record run paths and the code name ("PlasmeX runs/…", "plasmex.wannier_drift"). They do not appear in the PDF, but they do appear in an arXiv or journal source upload. Strip them, e.g. with `arxiv_latex_cleaner`, before uploading.
- Rebuild with the full pdflatex → bibtex → pdflatex ×2 chain after the edits. The current log has no undefined references, and bibtex reports 0 warnings.

---

## 6. Notes for the author (not for the referee report)

These point to material you already have that could answer some of the requests above. All of it predates the 17–18 Oct base cases (Tedge1000 / wannier_rel), so it would need rerunning on the current states before it can be quoted.
- **M3:** the mesh study `runs/mesh_base` (15 Oct) found first-order convergence, 1–2.3 % in the cathode sheath and < 0.4 % at the anode. It was on the previous base cases, and the table was removed from §numerics.
- **M1:** earlier anode runs (9–10 Oct) were not L-independent between 0.4 and 1 mm. Expect a referee to push on this most for the anode.
- **M8:** `runs/ion_inertia` / `doc/ion_mom.md` (17 Oct, before `wannier_drift = 1`) gave cathode φ_s 3.6 → 4.5 V, E_w +64 %, q_w −2.7 %. That is a direct answer to M8 if rerun on the current state.
