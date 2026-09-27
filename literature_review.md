# Literature review: what this paper can claim as new, and how to validate it

*Prepared 2026-09-26. Supersedes the novelty discussion in `response.md`, Question 7. That discussion relied on the old `todo.txt` (now deleted) and treated the heavy-gas momentum equation as a physics claim, which it is not (see §0).*

---

## 0. Two corrections that shape everything below

1. **The heavy-gas momentum equation is not a novelty.** In 1D, with an impermeable electrode on one side, the steady state has $\mathbf u=0$. Eq. (17) then reduces to $\nabla p=e(n_i-n_e)\mathbf E$, which is exactly the mechanical-equilibrium condition Almeida *et al.* (2008) impose (their eq. 9). In PlasmeX the momentum equation is a device for relaxing the transient acoustically. At steady state the model is physically the same as the mechanical-equilibrium model. The paper should say this in one sentence, and ideally show it numerically ($\mathbf u\to0$ and $\nabla p\to e(n_i-n_e)\mathbf E$ at convergence).
2. **Time dependence by itself is not a novelty either.** The Benilov group has published non-stationary unified 1D simulations: AC arc ignition on cold electrodes (Santos *et al.* 2021) and the development of the active-anode regime through a thermal instability (Santos *et al.* 2022).

---

## 1. What has already been done

| work | geometry / time | electrodes | sheath | electrode heat | conditions and outputs |
|---|---|---|---|---|---|
| **Almeida, Benilov & Naidis 2008**, J. Phys. D 41 245201 | 1D spherical, steady, Newton + Petukhov | cathode | resolved (unified), $n_i=0$ at the wall | $T_w$ **prescribed** | Hg at 100 bar, Ar at 1 bar; $T_w=3500$ K, $j_c=10^6$–$10^8$ A m$^{-2}$; profiles, CVC, $q_c(T_w)$; compared with Redwitz *et al.* 2005 |
| **Almeida *et al.* 2009** (cited as `almeida2009investigating`) | 1D, steady | anode | unified | – | very-high-pressure Hg |
| **Semenov, Krivtsun & Reisgen 2016**, J. Phys. D 49 105204 | 1D, steady | **anode** | unified; shows that the collisionless sheath model fails | – | **Ar, 1 atm, 500–2000 A cm$^{-2}$**; $U_a$ = 2.13 → 3.18 V; heat flux; agrees with existing theory and data |
| **Benilov, Almeida, Baeva *et al.* 2016**, J. Phys. D 49 215201 | 1D | cathode | compares unified, NLTE-sheath and 2T-ionization-layer-sheath | – | $10^5$–$10^7$ A m$^{-2}$: the approaches agree, sometimes "surprisingly well" |
| **Baeva, Benilov, Almeida & Uhrlandt 2016**, J. Phys. D 49 245205 | **2D** axisymmetric, steady (COMSOL) | both | **interfaces** (NLTE-sheath BCs), not resolved | cathode and anode bodies | free-burning Ar arc, 20–200 A; CVC |
| **Almeida, Cunha & Benilov 2017**, J. Phys. D 50 (aa7f97) | 1D near-anode, steady | anode | unified | rod model | Ar, Xe, Hg; **anode heating voltage** vs $j$, $T_w$, $p$; $q_a\propto j$, weakly dependent on $T_w$ |
| **Lisnyak, Cunha, Bauchire & Benilov 2017**, J. Phys. D 50 315203 | 2D LTE core + layer BCs | both | via the 1D near-cathode model; **anode: constant heating voltage** | yes | free-burning Ar arc, 20–200 A; arc voltage validated |
| **Khrabry, Kaganovich, Nemchinsky & Khodak 2018**, Phys. Plasmas 25 | 1D whole gap, steady | both (hot, ablating anode) | fluid and collisionless compared; same profiles | electrodes | short Ar arc; V–A characteristics |
| **Baeva, Loffhagen & Uhrlandt 2019**, PCPP | 1D whole gap, steady | both | unified | thermionic cathode conduction | TIG **micro**arc, 0.4 mm, $5\times10^3$–$2.3\times10^6$ A m$^{-2}$ |
| **Baeva *et al.* 2019**, PCPP (electron transport description) | 1D | both | unified | – | DC Ar microplasmas; effect of the electron transport model |
| **Liang & Trelles 2019**, PSST 28 115012 | 3D, NLTE | both | sheath folded in as an effective conductivity | cathode | TIG 100–200 A; field reversal at the anode |
| **Santos *et al.* 2021**, J. Phys. D 54 195202 | 1D whole gap, **transient** | both | unified | electrodes | AC arc ignition on cold W in Ar at 1 atm: displacement → ion → thermionic current |
| **Baeva *et al.* 2021**, J. Phys. D 54 025203 | 1D whole gap | both | unified | Cu electrodes | low-current short Ar/Cu-vapour arcs |
| **Santos *et al.* 2022**, Phys. Plasmas 29 043503 | 1D whole gap, **non-stationary** | anode focus | unified | yes | passive and active anode regimes; thermal instability |
| **Baeva 2024**, J. Phys. D 57 LT01 | 1D whole gap, quasi-static sequence | both | unified | Cu, 10 mm | air, 0.88 A, contact opening, multiple field reversals. States that "unified models … are currently mostly limited to one dimension" |
| **Santos *et al.* 2024**, J. Phys. D (ad5f3c) | 1D **reduced** near-cathode model for use as a BC | cathode | explicit sheath | coupled | ignition on cold refractory cathodes; designed for multi-D codes |
| **Gonzalez, Cayla, Freton & Teulet 2009**, J. Phys. D | 2D | cathode | sheath and presheath models (2T) | cathode body | free-burning Ar arc, 200 A |
| **Krivtsun** (anode processes in welding arcs) | 2D | anode | layer models | – | the anode fall is negative and non-uniform across the attachment; constricted vs diffuse modes |

A related but separate thread: a high-field correction to the ion mobility in the cathode sheath has been used in near-cathode models (one search result attributes to it a ~10× mobility reduction and a 2–3× thinner sheath). I could not pin down the source, so treat "field-dependent ion mobility in the sheath" as *not new* until it is checked.

### Consequences: what **not** to claim

- unified (sheath-resolved) modelling of the cathode layer: Almeida 2008;
- unified modelling of the anode layer in atmospheric argon: Semenov 2016; also Almeida 2017 and Santos 2022;
- unified modelling of both electrodes with the same equations: Khrabry 2018, Baeva 2019/2021/2024, Santos 2021/2022;
- time-dependent unified modelling: Santos 2021/2022;
- a collision-dominated sheath, or failure of the collisionless sheath: Almeida 2008, Semenov 2016;
- a cathode surface heat balance: Baeva 2019 (conduction in the cathode); Almeida 2008 already uses the same $k=100$ W m$^{-1}$ K$^{-1}$, 10 mm, 1000 K rod estimate that `cathode_inputs` uses;
- electrode boundary data for LTE codes: Lisnyak 2017, Almeida 2017 (anode heating voltage), Santos 2024.

---

## 2. Gaps that remain, and that PlasmeX can fill

They are ranked by strength of novelty times feasibility with the current code.

### N1. Two-dimensional (r–z) unified modelling of the near-electrode region, with the sheath resolved by AMR — **strongest**

- **Gap:** every sheath-*resolved* arc model found above is 1D. Baeva (2024) says explicitly that unified models are "mostly limited to one dimension". The 2D and 3D arc models treat the sheath as an interface (Baeva 2016), an effective conductivity (Liang & Trelles 2019), a 1D-model-derived BC (Lisnyak 2017, Santos 2024), or a separate sheath and presheath model (Gonzalez 2009).
- **What PlasmeX already has:** a single source tree that builds 1D or 2D r–z (`DIM=2`, `inputs_2d`); radial bulk-edge profiles $n_e(r)$, $T_e(r)$, $J_{\rm app}(r)$ read from an LTE arc solution (`arc_ne.dat`, `arc_Te.dat`, `arc_Japp.dat`); thermionic emission; AMR to sub-Debye resolution; and a Hypre potential solver for the anisotropic mesh.
- **Physics this would show that 1D cannot:** the radial variation of the fall voltage, heat flux and ion flux across the attachment; current spreading and constriction inside the non-equilibrium layer; and edge effects where the layer meets the cold gas. Krivtsun reports that the anode fall is negative and non-uniform across the attachment. That is a directly testable 2D prediction.
- **What is needed:** a converged 2D anode (and ideally cathode) case; a comparison of the 2D solution on the axis with the 1D solution; and radial profiles of $U$ and $q$. The cathode heat balance is 1D-only at present (`probinit` aborts in 2D), so a 2D cathode case needs a prescribed $T_w(r)$ or an extension of the balance.
- **Positioning risk:** low, provided you verify with a final targeted search (Baeva, Benilov, Semenov/Krivtsun, Khrabry, Trelles, 2023–2026) that nobody has published a 2D sheath-resolved arc electrode region in the meantime.

### N2. Energy budget of the collision-dominated sheath and the ion impact energy — **matches the title**

- **Gap:** the unified models do not contain an ion kinetic-energy term in the electrode heat flux. Almeida's eq. (32) has emission cooling, electron heating, ion neutralisation and heavy-particle conduction, and Almeida remarks that simplified collisionless models "have a term accounting for kinetic energy of incident ions" instead. Our $q_w$ (eq. 24) has the same structure. So in a collisional sheath, where does $e\,U_{\rm sheath}\,j_i$ go? It could go into gas heating ($Q_{in}$) and then conduction to the wall; into the column (Li & Benilov 2007 show that the power deposited in the near-cathode layer also goes to the column, mainly by electron enthalpy transport); or directly into the wall as ion kinetic energy. This is not quantified in the unified literature, as far as I found. It decides which ion term a reduced BC for LTE codes should contain.
- **What PlasmeX can report:** a full power budget of the sheath (electron Joule heating, ion Joule heating, $Q_{in}$, electron–heavy exchange, fluxes to the wall and to the column), and a mean ion impact energy at the wall from the field-heated ion temperature ($T_{\rm eff}$, eq. `eq:mu-in`) plus the drift, compared with the collisionless bound $eU_{\rm sheath}$.
- **Caution:** the model has $T_i=T_n=T_h$, so the ion impact energy is a diagnostic built from local quantities, not a transported quantity. Define it carefully in the paper, check that it stays below $eU_{\rm sheath}$, and state its limits. If it cannot be made robust, drop "ion impact energy" from the title.

### N3. A single framework for both electrodes, with a systematic "which physics matters" assessment — **moderate novelty, high usefulness**

- **Gap:** the papers above differ in transport closure (Chapman–Enskog/Stefan–Maxwell interpolation, drift-diffusion), in the ion wall BC ($n_i=0$ against flux-type), and in what they include (thermal diffusion, e–i friction on the ions, Coulomb corrections). Baeva 2019 examined the electron transport description for microplasmas, but I found no systematic sensitivity study for high-current electrode layers.
- **What PlasmeX can do:** the switches already exist: `thermal_diffusion`, `ion_ei_friction`, `electron_kappa`, `stefan_maxwell_drag`, the ion wall BC (`ion_bc_hi = 1` gives Almeida's $n_i=0$), and emission on/off. Produce a table of how $U_{\rm fall}$, $q_w$, $j_i/j$ and $T_w$ at both electrodes change as each term is removed. Reviewers value this, and it tells LTE modellers which effects their BCs must keep.
- The paper's transport section (M&K averages on LXCat cross sections, joined to Spitzer–Härm, with the $M$/$E$ distinction and Onsager-consistent thermal diffusion) is the natural basis for this. Present it as careful and consistent, not as new.

### N4. Electrode boundary data derived from a resolved layer, driven by an actual LTE bulk solution — **moderate**

- **Gap:** Lisnyak 2017 couples the LTE core to the cathode through a 1D model, but to the **anode through a constant heating voltage**. Almeida 2017 then computed the anode heating voltage. What is still missing is boundary data at *both* electrodes, from the same resolved-layer model, as functions of the bulk-edge state that an LTE solution actually delivers ($T_{\rm bulk}$, $n_e^{\rm bulk}$, $J$) and of $T_w$.
- **What is needed:** a coarse sweep that produces $U(J,T_{\rm bulk},T_w)$ and $q_w(J,T_{\rm bulk},T_w)$ with its partition, compared with Almeida 2017's anode heating voltage. This extends existing work rather than being "first".

### N5. Numerical method — **supporting, not a headline for J. Phys. D**

A time-marching, semi-implicit AMR solver (AMReX; 1D and 2D from one source; GPU-capable). It reaches steady states from arbitrary initial data without Newton continuation, and uses the heavy-gas momentum equation only to relax the acoustic transient. This is what makes N1 feasible, so describe it in Methods. It is a headline only in a methods journal (J. Comput. Phys., Comput. Phys. Commun.) or a companion paper.

---

## 3. Recommended positioning

**Option A: one paper.** 1D verification and validation at both electrodes (§4), plus the energy budget and ion impact energy (N2), plus the sensitivity table (N3), **plus one 2D anode demonstration** (N1). The 2D result is what makes the paper clearly different from the Benilov/Baeva/Semenov line.

**Option B: two papers.**
- *Paper 1 (this one):* "A sheath-resolving AMR model of the near-electrode regions of atmospheric-pressure argon arcs: verification, validation, and the energy budget of the collisional sheath". Contents: §4 validation, N2, N3, N4.
- *Paper 2:* "Two-dimensional unified modelling of arc attachment…" (N1).

**Draft contribution paragraph** (Option A; edit to match the final results):

> Unified models, in which a single set of transport equations together with Poisson's equation is solved up to the electrode surface, have established the structure of the near-electrode layers of high-pressure arcs [Almeida 2008; Semenov 2016; Khrabry 2018; Baeva 2019–2024; Santos 2021–2022]. All of them are one-dimensional, and multidimensional arc models still represent the sheath through interface conditions or effective properties [Baeva 2016; Lisnyak 2017; Liang & Trelles 2019]. The present work (i) introduces a time-marching, adaptively refined solver that resolves the space-charge sheath in one and two (axisymmetric) dimensions, with the bulk-edge state taken from an LTE arc solution; (ii) verifies it against the established one-dimensional unified results for the cathode and the anode of an atmospheric-pressure argon arc and validates it against probe and spectroscopic measurements; (iii) quantifies the energy budget of the collision-dominated sheath, including the energy with which ions reach the electrode, which is absent from unified models and assumed in simplified ones; (iv) assesses which transport effects control the fall voltage and the heat flux at each electrode; and (v) presents a first two-dimensional sheath-resolved solution of the near-anode region, showing the variation of the anode fall and heat flux across the arc attachment.

---

## 4. Validation plan

### V1. Verification (no external data)

| test | pass criterion |
|---|---|
| grid convergence with $\Delta x$ and $\Delta t$ varied **independently** | $U_{\rm fall}$, $q_w$ and the wall $n_e$, $E$ converge at the expected order |
| steady state of the heavy gas | $\mathbf u\to0$, and $\nabla p-e(n_i-n_e)\mathbf E\to0$ (Almeida's mechanical equilibrium) |
| conservation | total current uniform in $x$ at steady state; heavy-particle number conserved; the energy budget closes (sum of the terms in N2 equals the wall + column fluxes) |
| 2D vs 1D | the 2D solution on the axis tends to the 1D solution when the radial profiles are flat |

### V2. Code-to-code benchmarks (same equations, published numbers)

| benchmark | conditions to reproduce | compare | adjustments needed |
|---|---|---|---|
| **Almeida 2008, argon cathode** | $p=1$ bar, W, $T_w=3500$ K (fixed: `wall_temp_balance = 0`), $j_c=10^6,10^7,10^8$ A m$^{-2}$ | profiles (their Fig. 5), $U(j_c)$ (CVC), $q_c(T_w)$ (Fig. 12a), $j_i/j$, $j_{em}/j$ | set the ion wall BC to $n_i=0$ (`ion_bc_hi = 1`) for a like-for-like run, then show the effect of the flux BC. They use spherical $r_c=1$ mm; the difference matters only in the constriction zone. Their outer BC is local LTE with a radiation balance. `exe/almeida_arc_bc.py` already exists. |
| **Benilov *et al.* 2016** | cathode, $10^5$–$10^7$ A m$^{-2}$ | $U$, $q$ against three published approaches | – |
| **Semenov *et al.* 2016, argon anode** | 1 atm, $j=5\times10^6$–$2\times10^7$ A m$^{-2}$ | $U_a$ (2.13 → 3.18 V), layer thickness, heat flux | check their sign convention for $U_a$ and their wall BCs |
| **Almeida *et al.* 2017, anode heating voltage** | Ar, range of $j$ and $T_w$ | $q_a/j$ and its weak $T_w$ dependence | – |
| **Santos *et al.* 2022** | passive-anode regime | profiles, $U_a$ | – |

### V3. Experimental validation

| electrode | data | quantity |
|---|---|---|
| cathode | **Redwitz, Langenscheidt & Mentel 2005**, J. Phys. D 38 3143: high-resolution spectroscopy in front of W cathodes (and anodes), Ar, 0.1–0.3 MPa, 0.5–10 A | $n_e(x)$, $T_e(x)$. Almeida used 2.6 bar, 4 A, $T_w=3000$ K, $U=15$ V: he got $T_e$ to 10–20% but $n_e$ off by about 3×, which leaves room to improve. |
| cathode | **Nandelstädt, Luhmann, Mentel *et al.* 2002**, "Determination of HID electrode falls in a model lamp" I–III, J. Phys. D: pyrometry and Langmuir probes | $T_w$ and the cathode fall, for testing the heat balance |
| cathode | Dickson & von Engel 1967 (Proc. R. Soc. A 300 316); classic fall-voltage data | cathode fall about 11–15 V at 1 atm for ≥10 A |
| cathode | Zhou & Heberlein 1994; Haidar & Farmer | cathode surface temperature and near-cathode $T$ |
| anode | **Sanders & Pfender 1984**, J. Appl. Phys. 55 714: probes through a water-cooled Cu anode, 100–250 A | anode fall −2.1 to −1.4 V; $T_e$ at the anode 8800–9800 K; anode heat flux |
| anode | **Heberlein, Mentel & Pfender 2010**, J. Phys. D 43 023001 (survey) | compiled boundary-layer $n_e$, $T_e$, fall voltages and heat fluxes ($>10^7$ W m$^{-2}$) |
| anode | Tanaka & Ushio 1999 (`tanaka1999plasma`), Nestor 1962 | anode heat flux distribution |
| anode | Dinulescu & Pfender 1980 (theory) | fall −1.5 to −2.7 V |

### V4. System-level check (uses the LTE bulk)

Total arc voltage = $U_{\rm cathode}$ + $U_{\rm column}$ (from the LTE solution) + $U_{\rm anode}$, for the 200 A, 10 mm free-burning argon arc. Compare with Hsu *et al.* 1983 and the 20–200 A CVC used by Lisnyak 2017 and Baeva 2016. This is the check that shows the whole LTE-plus-resolved-layer approach works.

---

## 5. Suggested run matrix

| # | case | purpose |
|---|---|---|
| 1 | anode, 1D, Semenov conditions, 3 current densities | V2 |
| 2 | cathode, 1D, Almeida conditions, $T_w=3500$ K fixed, $n_i=0$ and flux BC | V2, and the effect of the ion BC |
| 3 | cathode, 1D, heat balance on, Redwitz/Nandelstädt conditions | V3, and $T_w$ prediction |
| 4 | anode, 1D, Sanders & Pfender conditions (100–250 A) | V3 |
| 5 | runs 1–4 with the physics switches toggled one at a time | N3 |
| 6 | energy-budget post-processing of runs 1–4 | N2 |
| 7 | anode, 2D r–z, bulk profiles from an LTE 200 A arc | N1 |
| 8 | $J$ sweep at both electrodes | N4, V4 |

---

## 6. References to add to `references.bib`

Not yet cited in the paper:
Semenov, Krivtsun & Reisgen 2016; Benilov *et al.* 2016 (the three-approach comparison; `Benilov2016` may already be this paper, so check); Baeva *et al.* 2016, 2019 (both), 2021, 2024; Almeida, Cunha & Benilov 2017; Lisnyak *et al.* 2017; Khrabry *et al.* 2018; Santos *et al.* 2021, 2022, 2024; Li & Benilov 2007; Gonzalez *et al.* 2009; Redwitz *et al.* 2005; Nandelstädt *et al.* 2002 (I–III); Heberlein *et al.* 2010; Dickson & von Engel 1967; Krivtsun (anode processes in welding arcs).

Several entries above were checked only at the abstract level. Confirm volume and page numbers from the publisher before adding them.

---

## 7. Positioning around the dynamics of the arc: "a transient LTE model needs a transient electrode patch"

*Added 2026-09-27 in answer to the question: an arc is itself dynamic (electrode length and shape, surface reactions and ablation all keep changing), so a transient LTE model needs a transient electrode-region model as a patch. What does the literature say, and what can PlasmeX do?*

### 7.1 The idea, and the objection a referee will raise

The idea is sound in spirit, but in the form "*a transient LTE model requires a transient electrode model*" a referee from the Benilov group will reject it, for this reason:

> The near-electrode layer relaxes in about a microsecond. Almost every bulk transient (flow, AC at 50 Hz, pulsing at kHz, erosion) is much slower. So the layer is **quasi-stationary**: a steady layer model, evaluated with the instantaneous bulk state and wall temperature, is enough. That is exactly what Lisnyak *et al.* (2017) and Santos *et al.* (2024) provide.

Santos 2024 says so directly: their near-cathode model is quasi-stationary, and they propose it as the boundary condition for multidimensional, time-dependent MHD simulations of arc ignition. So "the arc is transient, so the layer model must be transient" is not a sufficient argument. The version that survives review is narrower and **quantitative**: the layer model has to be transient *when, and only when*, one of the following holds.

1. **The event is faster than the layer's relaxation time.** Current zero in switchgear, ignition, re-ignition after polarity reversal, arc restrike in torches, spot formation, and high-frequency pulsing.
2. **The layer is quasi-stationary but the electrode is not.** The electrode's thermal time (ms to s) is comparable with an AC or pulse period. Then the layer state depends on the history of $T_w$, and a lookup table in the instantaneous current alone is wrong. Santos *et al.* 2021 show this at 50 Hz (§7.3).
3. **The steady problem has more than one solution.** Diffuse or spot cathode modes, passive or active anode, S-shaped CVCs. Which branch is realised depends on history, so a steady lookup table is ambiguous and a time-marching model decides it.
4. **The electrode surface itself drifts** (tip recession, gap growth, emitter depletion that changes the work function $\varphi$, metal vapour). These drifts are very slow (minutes to hours), so the layer does **not** need to be transient for them. It needs to be **parametric** in the surface state ($\varphi$, $T_w$, geometry, vapour fraction), and the erosion model needs the layer's *sensitivities* to that state.

So the defensible thesis is a **time-scale-resolved electrode boundary** rather than simply a transient one: a map that says which form of electrode coupling (quasi-stationary lookup, lookup plus electrode thermal inertia, or fully transient layer) is valid for a given arc transient, with the thresholds computed from a resolved model. PlasmeX is a time-marching, sheath-resolving solver, so it is well suited to computing those thresholds. As far as I found, nobody has published them.

### 7.2 Time-scale hierarchy (order-of-magnitude, Ar, 1 atm; my estimates, to be replaced by PlasmeX numbers)

| process | scale | time | how estimated |
|---|---|---|---|
| sheath formation (ion transit) | 1–10 µm, $v_i\sim10^3$–$10^4$ m/s | ≲ 10 ns | $d/v_i$ |
| ionization layer / presheath | 10–100 µm | ~ 1–10 µs | $d^2/D_a$; $1/(k_i n_a)$ |
| heavy-gas thermal layer | ~100 µm | ~ 1 µs | $d^2/\alpha_{\rm gas}$ |
| acoustic relaxation across 1 mm | 1 mm | ~ 1 µs | $L/c_s$ |
| electrode thermal response | 1 mm tip to 10 mm rod, $\alpha_W\approx3\times10^{-5}$ m² s⁻¹ | ~ 30 ms to 3 s | $L^2/\alpha_W$. Benilov 2020 quotes "fractions to tens of milliseconds" for cathode heating during ignition |
| erosion, tip reshaping, emitter depletion | – | minutes to hours | Nemchinsky 2012: ~1 µg s⁻¹ at 200 A; Qu *et al.* 2024 |

Bulk and operating transients the LTE model would impose:

| transient | typical time | regime (numbers in §7.1) |
|---|---|---|
| column flow, $L/u$ (10 mm at 100–300 m/s) | 30–100 µs | 1 barely; mostly quasi-stationary |
| torch restrike, anode-root motion | ~0.1–10 ms (75 Hz for a side-blown planar anode; kHz in spray torches) | 2 (anode thermal), 3 (reattachment) |
| AC HID lamps, AC TIG | 50 Hz to kHz | 2; 1 just after current zero |
| pulsed TIG | 1 Hz to ~20 kHz | 2; 1 at high frequency |
| switchgear current zero, $dI/dt\sim10$ A µs⁻¹ | µs | 1 |
| ignition on cold electrodes | µs to ms | 1 and 2 |
| erosion, gap growth, emitter depletion | min to h | 4 |

The table shows that "transient" is really three different problems. That structure is itself a useful contribution, and PlasmeX can fill in the numbers.

### 7.3 What the literature already shows

**Evidence that quasi-stationary coupling fails** (supports the idea):

- **Santos *et al.* 2021** (J. Phys. D 54 195202): a 1D unified AC simulation, $j=10^7\sin(100\pi t)$ A m$^{-2}$, electrode heat conduction non-stationary. Quote: "*the glow-to-arc transition in the AC case occurs in a way substantially different from the quasi-stationary glow-to-arc transition*". The electrodes start heating at significantly higher current densities than in the quasi-stationary case. **This is the strongest existing citation for the idea.** It covers only the first quarter-period after ignition, though, not the periodic steady AC state and not polarity reversal.
- **Bötticher, Graser & Kloss 2004** (J. Phys. D 37 55): in an HID model lamp, a current step gives diffuse → spot → diffuse cathode attachment. Their model is 3D heat conduction in the cathode with *nonlinear surface heating*, i.e. a quasi-stationary layer times a transient electrode. That is regime 2. It works because the layer is fast and the cathode is slow.
- **Flesch & Neiger 2003** (J. Phys. D 36 849), and the D2 automotive-lamp AC model: time-dependent plasma plus electrodes for AC HID lamps. The near-electrode non-LTE region is represented by a **non-LTE conductivity equation, not a resolved sheath**. They state that AC operation makes time-dependent electrode models "essential".
- **Phase-resolved AC HID measurements** (Mentel group: "boundary layers of ac-arcs at HID electrodes"; frequency dependence of electrode temperature; "from low to RF frequencies"): the electrode temperature rises in the anode phase and falls in the cathode phase, and a high cathode fall after current zero can trigger a diffuse-to-spot transition with a voltage peak. **These are ready-made validation data for a transient unified model** that nobody has used for one.
- **Benilov and co-workers on multiple solutions**: bifurcations (Benilov & Cunha 2003, PRE 68 056407); stability of steady solutions studied by evolving perturbations in time; spot–diffuse transitions with hysteresis; and **Almeida *et al.* 2017, Plasma Process. Polym. "Computing different modes on cathodes of DC glow and high-pressure arc discharges: time-dependent versus stationary solvers"** (abstract only; read it before claiming anything about solvers). Also the Appl. Phys. Rev. 11 031310 (2024) review of discharge modes and transitions between them.
- **Santos *et al.* 2022**: the active-anode regime develops through a thermal instability in a non-stationary simulation, i.e. regime 3 at the anode.
- **Trelles *et al.* 2007** (J. Phys. D 40 5937): with a 2T model, restrike of the torch anode attachment appears **without a separate reattachment model**, because of the non-equilibrium near the anode. The LTE torch models need an artificial reattachment criterion. This is a strong argument that bulk dynamics depend on the near-anode physics.
- **Baeva 2024**: the contact opening (gap growing with time) is computed as a **quasi-static sequence** of steady unified solutions, not as a transient. That is an open gap for a time-marching unified model.

**Evidence on the slowly changing electrode** (regime 4):

- **Nemchinsky 2012** (J. Phys. D 45 135201): evaporated W atoms are ionized almost immediately in the near-cathode layer and most return to the cathode. The **escape factor** is a few per cent, and the net erosion is ~1 µg s⁻¹ for a 200 A free-burning Ar arc. That calculation used an *ionization-layer* model with measured plasma data.
- **Qu *et al.* 2024** (Sci. Technol. Weld. Join.): experimental evolution of TIG cathode erosion and the resulting change in arc behaviour. Studies of Ce/La-doped tungsten report dopant migration to the tip, followed by volatilization, which changes $\varphi$ over time.
- **Javidi Shirvan *et al.* 2017** (Comput. Phys. Commun.): a coupling BC for an arc on a **non-homogeneous** cathode (spatially varying $\varphi$), with sheath, Knudsen layer and presheath; plus the review by Choquet *et al.*, "Gas tungsten arc models including the physics of the cathode layer: remaining issues" (Weld. World 2017).
- **Khrabry *et al.* 2018, 2020 (J. Appl. Phys. 128 123303)**, Kundrapu & Keidar: carbon-arc anode ablation, with low and high ablation regimes and a Knudsen layer. In the experiments the gap is held constant by moving the electrode as it ablates. **Musikhin, Nemchinsky, Gu, Koel & Raitses 2025** (arXiv 2506.11308): 2-colour pyrometry of a molten, ablating metal anode.
- **Baeva *et al.* 2021**: Cu vapour in low-current short arcs, treated in a unified model.

**Honest summary:** fully transient, sheath-resolved models exist only in 1D and only for ignition and anode instability (Santos 2021, 2022). AC periodic operation with **polarity reversal**, current-zero re-ignition, and the frequency at which quasi-stationary coupling breaks down have not been studied with a unified model, as far as I found. Erosion and ablation are handled with layer models that do not resolve the sheath.

### 7.4 What PlasmeX can do: concrete studies, ranked by cost

First, check what the current code supports (verified in `../PlasmeX`, 2026-09-27):

- `prob.wall_temp_balance` is **algebraic**: $q_w(T_w)=k(T_w-T_{\rm base})/L$, lagged one step, and "the cathode has no heat capacity in this model" (`src/aider.H`, `UpdateWallTemperature`). It is 1D only.
- `prob.arc_Japp` is a **constant** read by `arc_state()` (`exe/prob.H`). There is no current waveform.
- The heavy-gas momentum equation is there to relax the transient (§0). **Before any transient result is presented as physical, the heavy-gas transient must be shown to be time-accurate**: no artificial damping or speed-up, and $\Delta t$ convergence of the *transient response*, not just of the steady state. This is the main technical risk of the whole programme.

| # | study | code changes | answers | novelty |
|---|---|---|---|---|
| **T1** | **Relaxation-time map.** Apply small steps in $J$, in the bulk-edge state and in $T_w$, then fit the response of $U_{\rm fall}$, $q_w$ and $j_i/j$ to $\tau_{\rm layer}(j,T_w)$ at the cathode and the anode | $J_{\rm app}(t)$ (a few lines in `arc_state`) | the threshold for regime 1; checks §7.2 | high, cheap. I found no published relaxation times for the unified layer |
| **T2** | **Frequency sweep.** Drive $J(t)=\bar J+\Delta J\sin\omega t$ and plot the dynamic $U$–$J$ hysteresis loops against frequency. Find $\omega^\ast$, where the deviation from the quasi-stationary curve exceeds 5–10% | as T1 | the frequency limit of lookup-table BCs (pulsed TIG, HID) | high |
| **T3** | **Electrode thermal inertia.** Replace the algebraic balance with 1D transient heat conduction in the electrode (a small 1D solid grid, $\rho c\,\partial_tT=\partial_x(k\partial_xT)$, coupled through $q_w$) | moderate | regime 2; comparison with Bötticher 2004 and the Mentel AC data | medium on its own, high combined with T4 |
| **T4** | **AC polarity reversal.** One electrode, sinusoidal current: the cathode phase has emission on, the anode phase emission off (per the current BC logic, the switch must follow the sign of $J$). Run to a periodic steady state | T3 + sign-dependent BC | re-ignition after current zero, anode/cathode asymmetry of electrode heating, and whether the layer is quasi-stationary right after the zero crossing | **high**. Santos 2021 did only the first quarter-period |
| **T5** | **Surface-state sensitivities.** Sweep $\varphi$ (dopant depletion: pure W 4.55 eV, down to ~2.6–3 eV for doped or thoriated W), $T_w$ and the bulk-edge state. Report $\partial U/\partial\varphi$, $\partial q_w/\partial\varphi$, etc. | none (steady sweeps) | the input erosion and lifetime models need (regime 4) | medium, cheap, and it fits naturally with N4 |
| **T6** | **History dependence.** Reach the same final operating point by ramps from above and from below. Check for hysteresis (passive/active anode, S-shaped CVC) | $J_{\rm app}(t)$ | regime 3; whether a lookup table is single-valued | medium |
| **T7** | **Evaporated-metal escape factor.** Add W (or Cu) atoms and ions as species with an evaporation flux at the wall. Compute the escape factor *in the resolved sheath* and compare with Nemchinsky 2012 | new species plus cross sections (large) | net erosion rate from first principles | high, but a separate paper |
| **T8** | **Moving boundary** (gap growth, recession) | large (moving mesh or domain) | truly transient contact opening vs Baeva's quasi-static sequence | high, but late |

*T1 status (2026-09-27):* `PlasmeX/exe/timescales.py` fits the tail relaxation. The existing `anode_inputs` start-up run (0–8 µs) gives a single slow mode, $\tau=2.5$–$3.0$ µs, in every diagnostic (heavy-gas dominated); the charged-particle mode is faster than the 0.2 µs frame spacing. Not yet publishable: the layer is wider than the 200 µm domain (test $L$ = 200/400/800 µm, i.e. 2–8 times the electron–gas exchange length $\delta_T\approx100$ µm; see `response.md` Question 9), and step tests by restart with a changed `prob.arc_Japp` are still needed. See `response.md` Question 8.

**A deliverable that ties it together: a dynamic reduced BC.** From T1–T3, propose a low-order model an LTE code can use:
$$\tau_{\rm layer}(j,T_w)\,\frac{dU}{dt}+U=U_{\rm qs}(j,T_{\rm bulk},T_w,\varphi),\qquad \rho c\,\partial_t T=\nabla\cdot(k\nabla T)\ \text{in the electrode},$$
with $U_{\rm qs}$ and $q_{w,\rm qs}$ tabulated from steady PlasmeX runs (N4) and $\tau_{\rm layer}$ from T1. Also give the criterion $\omega\,\tau_{\rm layer}\ll1$ for when the lag can be dropped. This turns "the arc is dynamic" into a tool other groups can adopt, which is the kind of claim J. Phys. D and PSST referees accept.

### 7.5 How this fits the paper plan (§3)

- **This paper:** keep the steady validation core (§4, N2, N3). Add **T1** (and T2 if time allows) as a short section, "Dynamic response of the near-electrode layer". It justifies the time-marching solver physically rather than only numerically, and gives the paper a forward-looking hook. Say in the introduction that arcs are dynamic, then state the *quantitative* claim: *"we compute the relaxation time of the resolved layer and hence the range of validity of quasi-stationary electrode coupling for transient LTE simulations."* Avoid "a transient LTE model requires a transient electrode model"; Santos 2024 and Lisnyak 2017 contradict it as a general statement.
- **Paper 2:** T3 + T4 (AC with electrode thermal inertia and polarity reversal), validated against the Mentel-group phase-resolved HID data and compared with Santos 2021.
- **Paper 3:** T5 + T7 (+ T8): the electrode as a slowly changing boundary, from erosion rate and $\varphi$ drift to the arc response.
- The 2D line (N1) is independent and can go in either paper 2 or its own paper. Transient 2D (restrike, spot formation) is the long-term goal, and the only one PlasmeX's AMR is really needed for.

### 7.6 Other lines worth considering

1. **Uncertainty quantification.** Cross sections, the emission constants (Richardson $A$, $\varphi$, the Schottky treatment), the ion wall BC and the rod parameters all feed $U_{\rm fall}$ and $q_w$. A variance-based sensitivity study, i.e. which input dominates the uncertainty of $q_w$, is rare in this literature and cheap with a fast 1D code.
2. **Other gases**: N$_2$, air, CO$_2$ (switchgear, where current zero is regime 1). Baeva 2024 is in air. This needs molecular chemistry, which is a big step.
3. **Validation against transient experiments rather than steady ones.** The Mentel group's phase-resolved AC data (electrode temperature, electrode fall vs phase) are rarely used by unified modellers and would differentiate the paper.
4. **Surrogate or ML boundary model** trained on PlasmeX sweeps (N4 + T5), exported as a table or small network for LTE codes. Only worth doing after the physics is validated; on its own it is not a physics paper.

### 7.7 To verify before writing

- Read Almeida *et al.* 2017 (PPAP, time-dependent vs stationary solvers) in full. It may already contain relaxation-time or solver-comparison results that overlap with T1 and T6.
- Read Santos 2021 §4 and Santos 2022 for any stated relaxation times.
- Search again for "AC unified model arc polarity reversal", "current zero sheath-resolved atmospheric arc" and "relaxation time near-cathode layer" (2023–2026) before claiming T1 or T4 as first.
- The Qu 2024, Javidi Shirvan 2017, Flesch & Neiger 2003 and Almeida 2017 PPAP entries were checked only at the abstract or search-result level.

---

## Sources

- [Almeida, Benilov & Naidis 2008 — Unified modelling of near-cathode plasma layers](https://nano.pppl.gov/Publications/pubs/Almeida_UnifiedModeling.pdf) (local copy: `validation_fig/Almeida_2008_...pdf`)
- [Semenov, Krivtsun & Reisgen 2016 — Numerical study of the anode boundary layer in atmospheric pressure arc discharges](https://iopscience.iop.org/article/10.1088/0022-3727/49/10/105204)
- [Benilov et al. 2016 — Account of near-cathode sheath in numerical models of high-pressure arc discharges](https://iopscience.iop.org/article/10.1088/0022-3727/49/21/215201)
- [Baeva, Benilov, Almeida & Uhrlandt 2016 — Novel non-equilibrium modelling of a DC electric arc in argon](https://iopscience.iop.org/article/10.1088/0022-3727/49/24/245205)
- [Almeida, Cunha & Benilov 2017 — Computing anode heating voltage in high-pressure arc discharges](https://iopscience.iop.org/article/10.1088/1361-6463/aa7f97)
- [Lisnyak, Cunha, Bauchire & Benilov 2017 — Matching the LTE arc core with the electrodes](https://iopscience.iop.org/article/10.1088/1361-6463/aa76d3)
- [Khrabry et al. 2018 — Investigation of the short argon arc with hot anode, Part I](https://www.osti.gov/pages/biblio/1419791)
- [Baeva, Loffhagen & Uhrlandt 2019 — Unified non-equilibrium modelling of TIG microarcs](https://link.springer.com/article/10.1007/s11090-019-10020-x)
- [Baeva et al. 2019 — Fluid modelling of DC argon microplasmas: effects of the electron transport description](https://link.springer.com/article/10.1007/s11090-019-09994-5)
- [Liang & Trelles 2019 — 3D free-burning argon arc with a novel sheath coupling procedure](https://iopscience.iop.org/article/10.1088/1361-6595/ab4bb6)
- [Santos et al. 2021 — Numerical investigation of AC arc ignition on cold electrodes in atmospheric-pressure argon](https://digituma.uma.pt/bitstream/10400.13/4576/1/Numerical%20investigation%20of%20AC%20arc%20ignition%20on%20cold%20electrodes%20in%20atmospheric-pressure%20argon.pdf)
- [Baeva et al. 2021 — Unified modelling of low-current short-length arcs between copper electrodes](https://iopscience.iop.org/article/10.1088/1361-6463/abba5d)
- [Santos et al. 2022 — Regimes of current transfer to anodes of high-pressure arc discharges](https://pubs.aip.org/aip/pop/article/29/4/043503/2843564/Numerical-investigation-of-regimes-of-current)
- [Baeva 2024 — Reversal of the electric field and the anode fall in DC arcs in air during contact opening](https://iopscience.iop.org/article/10.1088/1361-6463/ad5c73)
- [Santos et al. 2024 — Near-cathode layer model for ignition of high-pressure arcs on cold refractory cathodes](https://iopscience.iop.org/article/10.1088/1361-6463/ad5f3c)
- [Li & Benilov 2007 — Effect of a near-cathode sheath on heat transfer in high-pressure arc plasmas](https://iopscience.iop.org/article/10.1088/0022-3727/40/7/024)
- [Gonzalez, Cayla, Freton & Teulet — Two-dimensional self-consistent modelling of the arc/cathode interaction](https://www.osti.gov/etdeweb/biblio/21457644)
- [Krivtsun — Anode processes in welding arcs](https://www.researchgate.net/publication/329970731_Anode_processes_in_welding_arcs)
- [Redwitz, Langenscheidt & Mentel 2005 — Spectroscopic investigation of the plasma boundary layer in front of HID electrodes](https://www.osti.gov/etdeweb/biblio/20665784)
- [Determination of HID electrode falls in a model lamp I–III](https://www.osti.gov/etdeweb/biblio/20310619)
- [Sanders & Pfender 1984 — Measurement of anode falls and anode heat transfer](https://pubs.aip.org/aip/jap/article-abstract/55/3/714/168524/Measurement-of-anode-falls-and-anode-heat-transfer)
- [Heberlein, Mentel & Pfender 2010 — The anode region of electric arcs: a survey](https://iopscience.iop.org/article/10.1088/0022-3727/43/2/023001)
- [Dickson & von Engel 1967 — Resolving the electrode fall spaces of electric arcs](https://royalsocietypublishing.org/doi/10.1098/rspa.1967.0172)
- [Khrabry et al. 2019 — Validated modeling of atmospheric-pressure anodic arc](https://arxiv.org/abs/1902.09991)

Added for §7 (transient positioning):

- [Santos et al. 2021 — AC arc ignition on cold electrodes (full text; the AC vs quasi-stationary comparison is in §5)](https://digituma.uma.pt/bitstream/10400.13/4576/1/Numerical%20investigation%20of%20AC%20arc%20ignition%20on%20cold%20electrodes%20in%20atmospheric-pressure%20argon.pdf)
- [Santos, Almeida, Benilova & Benilov 2024 — J. Phys. D 57 405202 (quasi-stationary layer model proposed as a BC for multi-D MHD)](https://iopscience.iop.org/article/10.1088/1361-6463/ad5f3c)
- [Benilov 2020 — Modeling the physics of interaction of high-pressure arcs with their electrodes: advances and challenges, J. Phys. D 53 013002](https://iopscience.iop.org/article/10.1088/1361-6463/ab47be)
- [Benilov 2008 — Understanding and modelling plasma–electrode interaction in high-pressure arc discharges: a review](https://iopscience.iop.org/article/10.1088/0022-3727/41/14/144001)
- [Bötticher, Graser & Kloss 2004 — Cathodic arc attachment in a HID model lamp during a current step, J. Phys. D 37 55](https://www.osti.gov/etdeweb/biblio/20458246)
- [Flesch & Neiger 2003 — Time-dependent simulation of plasma and electrodes in HID lamps with different electrode shapes, J. Phys. D 36 849](https://www.osti.gov/etdeweb/biblio/20372139)
- [Ac modelling of D2 automotive HID lamps including plasma and electrodes](https://www.osti.gov/etdeweb/biblio/20590905)
- [The boundary layers of ac-arcs at HID-electrodes: phase resolved electrical measurements and optical observations](https://www.osti.gov/etdeweb/biblio/20927522)
- [HID electrode temperature vs operating frequency, pyrometric measurements](https://www.osti.gov/etdeweb/biblio/21569773)
- [Development of the arc attachment at HID lamp electrodes from low to RF frequencies](https://www.osti.gov/etdeweb/biblio/21166413)
- [Benilov & Cunha 2003 — Bifurcation points in the theory of axially symmetric arc cathodes, PRE 68 056407](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.68.056407)
- [Stability of direct current transfer to thermionic cathodes II: numerical simulation](https://www.osti.gov/etdeweb/biblio/21009954)
- [Almeida et al. 2017 — Computing different modes on cathodes of DC glow and high-pressure arc discharges: time-dependent versus stationary solvers, PPAP](https://onlinelibrary.wiley.com/doi/10.1002/ppap.201600122)
- [Glows, arcs, ohmic discharges: an electrode-centered review on discharge modes and the transitions between them, Appl. Phys. Rev. 11 031310 (2024)](https://pubs.aip.org/aip/apr/article/11/3/031310/3303850/Glows-arcs-ohmic-discharges-An-electrode-centered)
- [Santos et al. 2022 — Numerical investigation of regimes of current transfer to anodes, Phys. Plasmas 29 043503](https://pubs.aip.org/aip/pop/article-abstract/29/4/043503/2843564/Numerical-investigation-of-regimes-of-current)
- [Trelles, Pfender & Heberlein 2007 — Non-equilibrium modelling of arc plasma torches, J. Phys. D 40 5937](https://iopscience.iop.org/article/10.1088/0022-3727/40/19/024)
- [Experimental study on the restrike mode of a DC arc anode attachment, PCPP](https://link.springer.com/article/10.1007/s11090-021-10196-1)
- [Nemchinsky 2012 — Cathode erosion in a high-pressure high-current arc: tungsten cathode in a free-burning argon arc, J. Phys. D 45 135201](https://www.osti.gov/etdeweb/biblio/22025007)
- [Qu et al. 2024 — Evolution rules for cathode erosion and corresponding arc behaviour in TIG welding, Sci. Technol. Weld. Join.](https://journals.sagepub.com/doi/full/10.1177/13621718231214948)
- [Javidi Shirvan et al. 2017 — Coupling boundary condition for high-intensity electric arc attached on a non-homogeneous refractory cathode, Comput. Phys. Commun.](https://www.sciencedirect.com/science/article/abs/pii/S0010465517302977)
- [Choquet et al. 2017 — Gas tungsten arc models including the physics of the cathode layer: remaining issues, Weld. World](https://link.springer.com/article/10.1007/s40194-017-0513-2)
- [Khrabry et al. 2020 — Analytical model of low and high ablation regimes in carbon arcs, J. Appl. Phys. 128 123303](https://pubs.aip.org/aip/jap/article/128/12/123303/1027201/Analytical-model-of-low-and-high-ablation-regimes)
- [Musikhin, Nemchinsky, Gu, Koel & Raitses 2025 — In situ studies of a molten metal anode ablation in a nearly atmospheric pressure DC arc (arXiv 2506.11308)](https://arxiv.org/abs/2506.11308)
