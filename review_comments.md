# Referee report

**Manuscript:** "A unified model for species fluxes, ion impact energy and heat transfer at the electrode of an atmospheric-pressure arc discharge"

**Journal:** Journal of Physics D: Applied Physics

**Version reviewed:** `main.pdf` built 14 Oct 2026 (sections `intro`, `model`, `bc`, `validation`, `results`, `conclusion`, `appendix`)

**Recommendation:** Major revision

---

## 1. Summary of the submission

The authors present a one-dimensional, time-dependent, unified drift-diffusion–Poisson model of the near-electrode layers of an atmospheric-pressure argon arc. The model has three species (e, Ar$^+$, Ar) and separate electron and heavy-particle energy equations. The species fluxes come from inertialess species momentum balances that keep the inter-species drag and the electron thermal force. A single-fluid heavy-gas momentum equation is solved, and the mesh is refined (AMR, AMReX) down to a fraction of the Debye length at the electrode. The outer boundary sits $L = 0.5$ mm from the electrode, where Dirichlet values from an LTE bulk state are imposed. A volumetric "radial loss" term, calibrated so that the imposed bulk state is a uniform column solution, stands in for lateral conduction. The cathode surface temperature follows from a simple 1-D conduction balance that includes thermionic emission with Schottky lowering.

The model is validated against the unified cathode-layer model of Almeida et al. (2008) and the anode-layer model of Semenov et al. (2016). It is then applied to one cathode case (21 kK, $1.2\times10^{8}\ \mathrm{A\,m^{-2}}$) and one anode case (12.4 kK, $9\times10^{6}\ \mathrm{A\,m^{-2}}$), with sweeps in bulk temperature and current density. The results cover the layer structure, the electron energy budget, an a-posteriori momentum and collisionality analysis of the sheaths, the electrode heat flux, and relaxation times after a 2 % current step.

The manuscript is carefully written, and the quantitative statements are mostly checkable from the tables. I verified several of them: the current partition at the cathode, the heat-flux decomposition, the 0.27 %/K emission sensitivity, the $1.9\times10^{8}\ \mathrm{A\,m^{-2}}$ emitted current density and the anode W/A figures. All are internally consistent. The a-posteriori momentum and Boltzmann/Bohm analysis (Section 5.3) is the most interesting part of the paper and is done better than in most comparable work.

However, the paper currently has four problems. It gives no description of the numerical method. It has no grid or domain-length sensitivity study. It reports several "electrode" quantities that depend on the truncation of the domain. And some of its headline conclusions are direct consequences of modelling choices (the 1-D cathode heat balance, the work function used at the anode) rather than findings. The novelty claimed in the introduction is also not exercised by the results. I detail these points below.

---

## 2. Major comments

### M1. There is no numerical-methods section

The introduction and the conclusion present AMR resolution of the sheath as a central feature, and the keywords name the code (PlasmeX). Yet the manuscript contains no description of:

- the spatial discretisation (finite volume? how are the drift–diffusion fluxes discretised — Scharfetter–Gummel, upwind?);
- the time integration (implicit? operator split? how the stiff dielectric-relaxation time and the electron energy equation are handled);
- the coupling between Poisson and the charged-species equations;
- the coarse–fine treatment at refinement interfaces;
- the refinement criterion (is the hierarchy actually adaptive, or a fixed geometric refinement towards the wall?);
- the linear and non-linear solvers;
- the treatment of the galvanostatic condition, equation (bc-current), as an "evolution equation for the boundary field".

None of the results can be reproduced from the paper as it stands. Please add a methods section. A `sections/methods.tex` file exists in the source tree but is not `\input` in `main.tex`; it may only need updating and including.

### M2. Grid convergence is not demonstrated

The only resolution statement is that the sheath cells are smaller than the local Debye length "by a factor of at least 1.5 at the cathode". That is barely resolved for a layer of 20–30 Debye lengths with an exponential electron density profile. It is also the region that sets $j_i/j$, $\phi_s$, the Schottky field and therefore $\Gamma_{\rm em}$ and $q_w$.

Please provide a grid-refinement study for at least one cathode and one anode case: one more and one fewer refinement level, or a refinement ratio of 2 instead of 4. Report the observed order and the estimated discretisation error in $\phi_s$, $E_w$, $j_i/j$ at the surface, $q_w$ and $T_{\rm wall}$. Quantities are currently quoted to three significant figures (e.g. 38.2 $\mathrm{MW\,m^{-2}}$, 9.96 V, $j_i/j = 0.218$ versus $0.216$). The precision of the validation claims cannot be judged without this.

### M3. Dependence on the domain length $L$ and on the outer boundary condition

The paper's main methodological claim is that the domain "need not extend far into the arc" because the bulk state is imposed at $x = L = 0.5$ mm. This claim is not tested, and several results suggest the outer boundary is influencing the solution:

1. **Anode, base case.** In figure 2(h), $j_i/j$ decreases monotonically all the way to $x = L$. Its extremum is *at the boundary*, so the ion current partition in the outer half of the domain is set by the Dirichlet data.
2. **Anode, validation.** The largest discrepancy (the $j_i/j$ minimum at $3.5\times10^{-4}$ m versus $4.3\times10^{-4}$ m) lies within 0.15 mm of the truncation. The explanation given — that gradients at $x = L$ "follow from each model's governing equations" — concedes that the solution near $L$ depends on where $L$ is placed.
3. **Scale analysis.** Table 1 gives $\lambda_T \approx 2$–$3.4\times10^{-4}$ m and $\lambda_{\rm rxn} \approx 1.8$–$2.6\times10^{-4}$ m at 11–12 kK, so at the anode $L$ is only about 1.5–2 relaxation lengths.
4. **Anode relaxation time.** The anode relaxation time of 11 µs, and $t_{99} \approx 40$–$75$ µs for $j_i$, look like diffusive times over the domain. If so, they scale with $L^2$ and are artefacts of the truncation.

Please repeat at least the two base cases with $L = 1$ mm and 2 mm (same bulk state and same $\Lambda$). Report how $\phi_s$, $j_i/j$ at the wall, $q_w$, $T_{\rm wall}$ and the relaxation times change. Without this, neither the methodological claim nor the anode results can be accepted.

A related point: a Dirichlet condition on *all* of $n_e$, $n_i$, $T_e$ and $T_h$ at $x = L$ reflects any disturbance that does not have exactly the imposed state at $L$. Please discuss whether a non-reflecting or Robin-type outer condition was considered, and why Dirichlet was chosen.

### M4. Reported voltages include the ohmic drop of the column inside the domain

Table 2 and figures 6–7 report $\phi(L)$, the "voltage across the domain". This includes the ohmic drop of the column section contained in the domain, so it is proportional to the arbitrary $L$:

- **Cathode:** the column field is $1.1\times10^{4}\ \mathrm{V\,m^{-1}}$, which accounts for roughly 5.5 V of the 13.7 V over 0.5 mm.
- **Anode:** figure 2(f) gives a column slope of about $2\times10^{3}\ \mathrm{V\,m^{-1}}$ beyond 100 µm, i.e. about 1 V of the 2.6 V.

The usual and transferable quantity is the near-electrode voltage drop, $U_{\rm ne} = \phi(L) - E_{\rm col}L$, the linear extrapolation of the column potential back to $x = 0$. Only that quantity can be passed to an LTE column model, which is the stated purpose of the work (conclusion, last paragraph). Please report it instead of, or alongside, $\phi(L)$. Revise the sweep discussion ("the voltage across the domain from 11.4 to 15.4 V"; "the anode lies below the bulk-plasma boundary by 1.6 to 3.6 V") accordingly.

### M5. The cathode heat flux conclusion follows from the electrode model, not the plasma

The paper's strongest stated conclusion is that "the heat flux into the cathode body is set by the electrode rather than by the plasma" ($q_w = 26.3$–$27.3\ \mathrm{MW\,m^{-2}}$ while the surface heating changes by 50 %). The cathode model is

$$
q_w = \frac{k_w}{L_w}\,(T_{\rm wall} - T_0),
$$

so a 1 $\mathrm{MW\,m^{-2}}$ change in $q_w$ needs a 100 K change in $T_{\rm wall}$. That change in $T_{\rm wall}$ changes emission cooling by about 30 %. So *any* plasma model coupled to this 1-D conduction law and Richardson–Schottky emission gives a nearly constant $q_w$. The result is a property of the assumed boundary condition, and the parameters are acknowledged to be "crude".

The paper should either:

- (a) present this explicitly as a consequence of the cathode model and show how $q_w$ and $T_{\rm wall}$ move when $k_w/L_w$ is varied (e.g. by $\times 0.3$ and $\times 3$), or
- (b) couple to a more realistic cathode thermal model (e.g. a 2-D axisymmetric rod with a non-uniform attachment).

Related points:

1. **Melting point.** The computed $T_{\rm wall} = 3698$ K, and $T_{\rm wall}$ up to 3727 K in the sweep, *exceed the melting point of tungsten (3695 K)*. $W = 2.63$ eV for thoriated tungsten is also hard to defend at these temperatures, where the thoria monolayer is not maintained. Please comment. Possibly the electrode parameters ($L_w$, $T_0$ or the attachment model) need revisiting.
2. **Radiative cooling is missing.** Equation (bc-qwall) has no thermal radiation from the cathode surface. At 3700 K, $\varepsilon_{\rm rad}\,\sigma_{\rm SB}T_{\rm wall}^4 \approx 3$–$4\ \mathrm{MW\,m^{-2}}$ for an emissivity of 0.3–0.35, which is 12–15 % of $q_w$. That is larger than several effects the paper discusses (e.g. "$q_w$ by about 1 %").
3. **Above $1.5\times10^{8}\ \mathrm{A\,m^{-2}}$.** The statement "Above $1.5\times10^{8}\ \mathrm{A\,m^{-2}}$ no steady state was reached with this cathode heat balance" needs to be explained. Is this a loss of existence of the steady solution (e.g. a turning point, which would be physically significant and connect to the known multiplicity of solutions of the plasma–cathode interaction, Benilov et al.), or a failure of the solver? A reader cannot tell.
4. **Emission exceeds the arc current.** $\Gamma_{\rm em}$ exceeds the arc current by a factor of 1.44–1.77 in every cathode case, with about 40 % of the emitted flux returned by plasma electrons. Please put this regime in context with the literature on thermionic arc cathodes, where $j_{\rm em}/j < 1$ is more common at these current densities.

### M6. The anode work function

The anode is described as water-cooled with $T_{\rm wall} = 1000$ K, yet $W = 2.63$ eV (thoriated tungsten) is used at the anode too ("The electrode material has $W = 2.63$ eV"). Water-cooled anodes are normally copper ($W \approx 4.5$–$4.65$ eV), and in GTAW the anode is the workpiece. The work-function term $W\lvert J_{\rm app}\rvert/e$ is the single largest term of the anode heat flux (62 %; "2.63 W per ampere"). With $W = 4.6$ eV, $q_w$ would rise from 38 $\mathrm{MW\,m^{-2}}$ to about 56 $\mathrm{MW\,m^{-2}}$, and the "4.1–4.4 W per ampere" result would become about 6 W per ampere.

Please use a work function appropriate to the anode material, or justify the choice. Also compare the resulting effective anode heat-transfer voltage with measurements (e.g. Nestor 1962; Sanders & Pfender 1984; Heberlein et al. 2010; Tanaka et al.). An experimental comparison would greatly strengthen the heat-transfer part of the paper, which the title emphasises.

### M7. The validation does not cover the quantities the paper is about

The title highlights species fluxes, ion impact energy and heat transfer. The validation compares profiles, $j_i/j$ and the sheath drop, but not the electrode heat flux or the near-electrode voltage. Almeida et al. (2008) report integral characteristics of the near-cathode layer, including the voltage drop and the energy flux to the cathode, and Semenov et al. discuss the anode heat flux. Please include $q_w$ (or the equivalent heat-transfer voltage) and the near-electrode voltage in the validation, using the same $W$ and $A_R$ as the reference works. State those values explicitly for the validation cases. The present text gives $W = 4.55$ eV for the cathode only and gives $A_R$ for neither validation case.

Further points on the validation:

1. **Test the explanation of the sheath discrepancy.** The factor-2.5 density difference and the 40 % field difference in the cathode sheath are attributed, "most likely", to ion transport at high $E/N$ and to the ion boundary condition. Both are easy to test: run with Almeida's ion mobility and/or their ion wall condition. A unified model's main selling point is the sheath, so this discrepancy should be explained rather than speculated about.
2. **Geometry and the outer region.** Almeida's spherical geometry makes $j$ differ by a factor of 2.25 at $x = L$, so Joule heating $j^2/\sigma$ in the outer part of the domain differs by a factor of about 5. The good temperature agreement in the outer half of the domain therefore reflects the shared Dirichlet data rather than the physics. Please say so, and restrict the validation claims to the region where the two problems are actually the same.
3. **Precision of the claims.** "The temperatures agree within 3–5 % and the charged-particle densities within 15 %" should be qualified by the discretisation error (M2) and by the extraction error of the digitised curves.

### M8. Ion impact energy: the title promises it, but the paper declines to compute it

Section 5.4 correctly argues that the local Wannier closure cannot give the impact energy in a sheath with only 1–2 collisions, and reports only an upper bound (5.6 eV, $\le 38\ \mathrm{MW\,m^{-2}}$). That is honest, but then the title ("…ion impact energy…") overstates the contribution. Either retitle the paper or provide an actual estimate. For example, a test-particle Monte Carlo of Ar$^+$ through the computed sheath field, with charge exchange and elastic collisions, is cheap and would answer the question directly.

Please also clarify the energy accounting:

- In the model, the energy the ions gain from the field, $e\,\boldsymbol{\Gamma}_i\cdot\mathbf{E}$, is deposited in the heavy gas through $Q_{in}$. Part of it reaches the wall through the gas-conduction term $-\kappa_h\nabla T_h\cdot\hat{\mathbf{n}}$, which is already 105 $\mathrm{MW\,m^{-2}}$ at the cathode. Adding the kinetic impact energy on top would partly double-count it. Please state how much of $e\,\boldsymbol{\Gamma}_i\cdot\mathbf{E}$ in the sheath is already contained in the wall conduction.
- The bound of 38 $\mathrm{MW\,m^{-2}}$ is comparable to $q_w$ itself. The statement that it would change $T_{\rm wall}$ by only about 20 K relies again on the stiffness of the emission cooling (see M5).

### M9. Physical validity of the 1-D model and of the radial-loss closure at the cathode

1. **Column radius versus domain length.** The calibrated cathode column radius is $R = 0.42$ mm, smaller than the 0.5 mm domain and comparable to $\lambda_T$ and $\lambda_{\rm rxn}$. A planar 1-D description of the near-cathode layer needs the layer thickness to be small compared with the radius of the attachment. Here the two are comparable, and the volumetric radial loss removes 574 $\mathrm{MW\,m^{-2}}$ — 36 % of the Joule heat in the domain — including inside the ionization layer. Please discuss the 1-D validity at the cathode explicitly. Quantify the radial loss within $x < 100$ µm compared with the other terms, instead of only stating that it is "negligible there".
2. **The loss form.** The Elenbaas-type form $Q = \Lambda\,\kappa(T)\,(T - T_{\rm edge})$ uses the local conductivity times the full temperature difference to the edge. This is a strong linearisation of the heat-flux (Kirchhoff) potential $S = \int_{T_{\rm edge}}^{T}\kappa\,dT'$. Please justify it or use $S$. The value of $T_{\rm edge}$ (300 K in the computations, as far as I can tell) is not stated in the paper.
3. **$\kappa_R$ in non-equilibrium regions.** $\kappa_R$ is an LTE quantity, yet it is evaluated at the local $T_e$ throughout the domain, including the strongly non-equilibrium anode layer where $T_e \ne T_h$ and $n_e \gg n_{\rm Saha}$. Please justify this or restrict $\kappa_R$ to the column part.

### M10. Pressure boundary conditions appear over-determined

Section 3.3 states that "the pressure at the electrode surface is taken to be atmospheric". At $x = L$, the neutral density follows from Dalton's law at $p_0$. In a steady state with $\mathbf{u} = 0$, equation (bulk-balance) gives

$$
p(L) - p(0) = \int_0^L e\,(n_i - n_e)\,E\,dx = \frac{\varepsilon_0}{2}\left(E_L^2 - E_0^2\right).
$$

At the cathode, $E_0 = 2.2\times10^{7}\ \mathrm{V\,m^{-1}}$, so this is about 2.1 kPa ($\approx 2\,\%$ of $p_0$). Both ends cannot then be at $p_0$. Please state precisely which pressure condition is imposed at the wall (or which is an output), and how the heavy-gas momentum equation is closed at both ends.

### M11. Steady-state criterion versus the relaxation times reported

The steady state is declared when "the relative change of every variable over 1 µs fell below $10^{-3}$". This criterion bounds a rate of change, not the distance from the steady state. How far it leaves the solution from its fixed point depends on the relaxation time of the case, and Section 5.5 shows that this time differs by a factor of about 40 between the two electrodes.

**1. A rate criterion leaves a distance of $\varepsilon\,\tau/t_{\rm ref}$.** Near a stable steady state $U^*$, the slowest eigenmode dominates the remaining departure $\delta = (U - U^*)/U^*$ once the faster modes have died out:

$$
\delta(t) \approx \delta_0\, e^{-t/\tau}, \qquad \left|\frac{d\delta}{dt}\right| = \frac{|\delta|}{\tau}.
$$

The test "relative change over $t_{\rm ref}$ below $\varepsilon$" is $|d\delta/dt|\,t_{\rm ref} < \varepsilon$. It is therefore first met when

$$
|\delta| \approx \varepsilon\,\frac{\tau}{t_{\rm ref}}.
$$

That is, the remaining error is the tolerance multiplied by $\tau/t_{\rm ref}$. The solver does not set this ratio; the physics of the case does. With $\varepsilon = 10^{-3}$ and $t_{\rm ref} = 1$ µs:

| | $\tau$ (Section 5.5) | distance left, $\varepsilon\,\tau/t_{\rm ref}$ |
|---|---|---|
| cathode, acoustics removed | 0.29 µs | $3\times10^{-4}$ |
| cathode, physical sound | $\approx 1$ µs | $10^{-3}$ |
| anode | 11 µs | $1.1\times10^{-2}$ |

The same tolerance therefore gives a cathode converged to $10^{-3}$ or better but an anode converged only to about 1 %. For the error left at the anode to match the stated $10^{-3}$, the tolerance would have to be $\varepsilon \lesssim 10^{-4}$. A safer choice is $10^{-5}$ or below, because $\tau$ grows towards lower bulk temperature (Section 5.5).

**2. Three features of the test can make it weaker than this estimate.** The paper does not specify any of them, and each matters:

- *Which quantities are tested, and in which norm.* The quantities the paper reports ($j_i$, $q_w$, $\phi_s$, $E_w$) are surface fluxes and derived quantities, not the transported variables. The surface ion current at the anode is $0.4\,\%$ of the total current, a small difference of large terms, so its relative drift can be much faster than that of any state variable. If the "relative change of every variable" is a norm over the whole domain or level, the bulk values dominate it, and a drift confined to the electrode layer, where all the reported quantities are evaluated, hardly registers. The estimate above then holds for each reported quantity separately, $\delta_q \approx \tau\,|\dot q/q|$, with that quantity's own rate. The variable test does not bound that rate.
- *How long the test must hold.* Requiring the criterion over several consecutive checks adds nothing if the checks span a time much shorter than $\tau$. Over a window $\Delta t_w$ the slow mode decays only by $e^{-\Delta t_w/\tau}$.
- *Which integrator is used.* For a time-accurate integration, the estimate $\delta \approx \varepsilon\tau/t_{\rm ref}$ applies as it stands. For a backward-Euler pseudo-transient iteration with $\Delta t \gg \tau$, each step multiplies $\delta$ by $\tau/(\tau + \Delta t)$. After a few steps past the trigger, the state is then converged far below the tolerance. The same tolerance thus means very different things for the two kinds of integration. The paper does not say which one was used for which case (see also M1).

**3. Why it matters here.** An error of 1 % is comparable to the last quoted digit of several anode results. $j_i/j = -0.0037$ is given to two significant figures, a rounding interval of $\pm 1.4\,\%$, and $\phi_s = 0.41$ V is given to $\pm 1.2\,\%$. The ion current, on which the anode discussion relies, is the slowest quantity of all: the paper's own $t_{99}$ is 48 µs for $j_i$, against 20 µs for $q_w$. The anode validation case (Semenov et al.) is at the same risk. Its differences from the reference ($j_i/j = -0.0068$ versus $-0.0065$, i.e. 4.5 %; $\phi_s$ within 4 %) are larger than a 1 % error, so the conclusions of that comparison probably stand. But an error of this size should be demonstrated, not assumed. The cathode results are not affected at the stated tolerance.

**4. The residual that is reported does not close the question.** The paper reports that the sum of the terms of the electron energy equation, $\partial\mathcal{E}/\partial t$, is below $1.2\times10^{-3}$ of its largest term. This too is a rate, scaled by a local process rather than by the distance from the solution, and it covers one equation out of six. The electron and ion continuity equations, Poisson's equation, the heavy-particle energy and momentum equations, and the cathode surface heat balance ($q_w$ against the conducted flux $k_w(T_{\rm wall} - T_0)/L_w$) are not reported. The heat balance is an algebraic condition, and a test on the transported variables does not watch it.

Please:
1. state the criterion precisely: the variables tested, the norm, the number of consecutive checks and their spacing, and the integrator and time step at which it was applied for each case;
2. demonstrate the remaining distance for the reported quantities. Either continue each run with the tolerance reduced by $10^{2}$–$10^{4}$ and report the change of $j_i$, $q_w$, $\phi_s$, $E_w$ and $T_{\rm wall}$, or give $\delta_q \approx \tau\,|\dot q/q|$ for each, using the $\tau$ of Section 5.5;
3. report a residual for every equation solved, including the cathode surface heat balance;
4. do the same for the two validation cases.

### M12. Neglecting species inertia

Section 5.3 shows that the neglected inertia exceeds the drag in the inner 0.36 µm of the cathode sheath (ions, by a factor of 2.4 at the surface). For electrons near the cathode it exceeds the drag by up to a factor of 54 and is a third of the electric force at the surface. The authors are commendably frank about this. But the quantities the paper is about — $j_i/j$ at the surface, $n_i$ at the surface, the Schottky field and hence $\Gamma_{\rm em}$ and $q_w$ — are all evaluated in exactly this region.

Please quantify the consequence. One option is a run in which the ion momentum equation keeps the inertia term $M n_i(\mathbf{u}_i\cdot\nabla)\mathbf{u}_i$ (the ion is the critical species, and this is a standard extension). At a minimum, give an estimate of how much the surface ion density and $E_w$ could change. Related: the emitted electrons form a nearly monoenergetic beam (cold, at $T_{\rm wall}$, accelerated by about 4 V) that dominates $n_e$ at the cathode surface. Treating them as part of a single Maxwellian electron fluid at $T_e \approx 19.5$ kK deserves discussion. Benilov's group and others treat the emitted electrons separately.

### M13. The claimed novelty is not demonstrated by the results

The introduction claims three advantages:

- (i) the outer boundary can carry non-zero gradients taken from an LTE model or measurements;
- (ii) AMR makes the approach "especially helpful in multidimensional simulations";
- (iii) the model is intended to supply boundary conditions to an LTE arc model.

None is exercised:

- every bulk state is a uniform, calibrated column solution;
- all calculations are 1-D, and the four-level refinement towards the wall could equally be a fixed graded mesh (M1);
- no coupling to an LTE model, and no extraction of a boundary-condition table, is shown.

Either reduce the claims to what the paper shows, or add a demonstration — e.g. one case with a bulk state and gradient taken from the authors' LTE arc model (Korsipati et al. 2026), and a table or fit of $U_{\rm ne}$, $j_i/j$ and $q_w$ versus $(T_{\rm bulk}, J_{\rm app})$ that an LTE solver could actually use.

---

## 3. Minor comments

### Physics and model

1. **Electron wall flux (eq. bc-electron).** The flux is pinned at $\tfrac14 n_e\bar c_e$, while the ion flux (eq. bc-ion) includes a drift term. The asymmetry should be justified. Note also that Hagelaar et al. (2000), whom the authors cite, use $\tfrac12 v_{\rm th}$ with a reflection factor $(1-r)/(1+r)$, not $\tfrac14$, when the drift term is added. Please state precisely how eq. (bc-ion) differs from Hagelaar's form and why.
2. **Ion thermal speed (eq. bc-ion).** $\bar c_i$ is evaluated at $T_{\rm wall}$. Near the cathode the ions have $T_{\rm eff} \approx 48$ kK, and the heavy gas at the wall is at $T_{\rm wall}$ only by assumption. Please justify this choice, or use the local ion temperature.
3. **Secondary electron emission** is neglected without an estimate. With $\gamma \approx 0.1$ for Ar$^+$ on W, $\gamma\Gamma_i \approx 4\times10^{24}\ \mathrm{m^{-2}\,s^{-1}}$, below 0.4 % of $\Gamma_{\rm em}$. Saying so in one sentence would close the point.
4. **Heavy-particle energy equation (eq. henergy).** It contains no diffusive enthalpy flux of the ions relative to the heavy frame, $\tfrac52 k_BT_h\,n_i(\mathbf{u}_i - \mathbf{u})$. This is non-negligible in the cathode layer, where the ions carry a substantial flux through a strong $T_h$ gradient. Please comment.
5. **Transport bridge (eq. bridge).** The overestimate of $\sigma$ at intermediate ionisation (Appendix B.4, figure 9) is described qualitatively only. Please give the maximum deviation and the temperature at which it occurs. That range (5–12 kK) covers the whole anode layer, so it bears directly on the anode results.
6. **Ion mobility.** Section B.2 notes that evaluating $\overline{Q}_{in}$ at $T_{\rm eff}$ instead of $T_h$ brings the mobility within 7 % of the Monte Carlo data (against 15–26 % now), and that the cost is an iteration. At the cathode surface, $E/N \approx 10^3$ Td, exactly where the present closure is 20–25 % low. One local scalar iteration per face is cheap. Please adopt the better closure, or show that the results are insensitive to it.
7. **Mobility comparison.** It is made at 300 K, while the sheath gas is at 1000–3700 K. Please comment on how well the closure carries over to elevated $T_h$.
8. **Bohm criterion (Section 5.3).** The comparison uses the sheath edge defined by $\lvert n_i - n_e\rvert/n_i = 1\,\%$. Even in a collisionless sheath the ion velocity at a 1 % charge-separation point is below $u_B$, because the transition is smooth (Riemann 1991; Franklin). Please show that the conclusion holds for other reasonable edge definitions (0.1 %, 10 %), or use the Riemann asymptotic sheath edge. The anode conclusion ($0.13\,u_B$) is robust; the cathode one ($0.82\,u_B$) may not be.
9. **Notation for $u_B$.** $u_B$ is defined as $\sqrt{k_B(T_e+T_i)/M}$ in the text and as $\sqrt{k_B(T_e+T_h)/M}$ in the caption of figure 5. Please make them consistent.
10. **Removing axial acoustics (Section 5.6).** Scaling the heavy-gas inertia by $\epsilon = 10^{-4}$ is a reasonable device, but it should be presented as a modelling assumption, not a numerical one. It also suppresses real axial pressure dynamics (e.g. the cathode jet). Separately, figure 8(a) still shows oscillations in $q_w$ between $10^{-8}$ and $10^{-7}$ s after the acoustics have been removed. Please identify what they are.
11. **Anode response time scaling.** For the anode response, please state whether $t_{99}$ scales with $L$ (see M3).
12. **Total-current statement.** "The total current density is uniform to within $10^{-8}$" — please say whether this is conduction plus displacement current, and whether it is relative to $J_{\rm app}$.

### Presentation and consistency

13. **Placeholders.** The abstract is empty. The author and affiliation are placeholders. The acknowledgements say "To be completed". There is no data/code availability statement, which IOP requires.
14. **Wrong cross-reference.** In the first sentence of Section 5.1, "as estimated in Section 3.2" (`sec:model-governing`) should point to the scale analysis (`sec:scales involved`).
15. **Table 1 caption.** "$\lambda_i \approx \lambda_e$" is not true at 11 and 12 kK (factors 1.9 and 1.5).
16. **Symbol clashes.** Several symbols carry two meanings:
    - $\Lambda$ is both the radial-loss coefficient and the Coulomb-logarithm argument ($\ln\Lambda$);
    - $\phi$ is the potential and $\varphi$ the Mitchner–Kruger thermal-diffusion coefficient;
    - $\lambda$ denotes mean free paths, the Debye length and the electron thermal conductivity (eq. mk-q);
    - $C$ is both the electron speed and the thermal-force coefficients $C_e$, $C_i$, $C_n$;
    - $\varepsilon$ / $\epsilon$ is used for $\varepsilon_0$, $\varepsilon_{\rm ion}$, the electron energy and the inertia scaling;
    - the wall temperature appears as $T_{\rm w}$, $T_{\rm wall}$ and $T_{\rm s}$;
    - the axis scale is $x_c$ in captions and $X_C$ elsewhere.
17. **Undefined at first use.** $n_h$ and $\rho$ are used before they are defined (eq. flux-i; Section 3.2).
18. **"Exactly $W\lvert J_{\rm app}\rvert/e$" (Section 5.5, anode).** The absorbed electron flux is $\lvert J_{\rm app}\rvert/e + \Gamma_i$, because the ions also flow in, so the work-function term is not *exactly* $W\lvert J_{\rm app}\rvert/e$ (the difference is 0.4 %).
19. **Cathode budget.** $543 + 153 = 696\ \mathrm{MW\,m^{-2}}$ against the quoted $690\ \mathrm{MW\,m^{-2}}$. The conduction part at $x = L$ ($\approx -6\ \mathrm{MW\,m^{-2}}$) should be stated for the numbers to close.
20. **Anode surface densities.** "$n_e$ exceeds $n_i$ by 9 %" at the anode surface, but table 2 gives $5.2/4.7$, i.e. 11 % — or 9.6 % if normalised by $n_e$. Please make the definition explicit.
21. **Radial-loss parameters.** State $T_{\rm edge}$, the values of $\Lambda$ (not only $R$), and that $R$ is quoted from $\Lambda = 4/R^2$.
22. **Validation electrode parameters.** State $A_R$ and $W$ explicitly for each validation case. Section 3.3 currently reads as if $W = 2.63$ eV and $A_R = 3\times10^{4}\ \mathrm{A\,m^{-2}\,K^{-2}}$ apply everywhere.
23. **Resolution wording.** "Smaller than the local Debye length by a factor of at least 1.5" is ambiguous. Give cells per Debye length and the finest cell size (15 nm?) explicitly.
24. **Equation (heavy-momentum)** ends with a full stop followed by "which governs…".
25. **Language.**
    - "1 cm far from the electrode" → "1 cm from the electrode";
    - "based on scales estimates" → "based on the scale estimates" (also a double space);
    - "which is what makes arc discharges modelling challenging" → "arc-discharge modelling";
    - "Nor do field reversals need to be anticipated" is a fragment.
26. **Figure 2(d).** $T_{\rm eff}$ is listed in the legend of the anode panel but is indistinguishable from $T_h$. Either drop it there or say so in the caption.
27. **Figure 8.** The time axes start at $10^{-12}$ s, so the initial ($f = f_0$) state is not visible. Consider adding $t = 0$ or a linear inset.
28. **Results run records.** The `% Run record` comments in `results.tex` are commented out and do not appear in the PDF. That is fine, but please make sure they are removed before submission.

### References

29. **Mitchner & Kruger (1973)** is entered as an `@article` with the publisher in the title field. It should be a `@book` (Wiley, New York).
30. **Zhdanov (2002)** is cited as *Plasma Phys. Control. Fusion* 44, 2283. That page is a book review. Cite the book itself (*Transport Processes in Multicomponent Plasma*, Taylor & Francis, 2002).
31. **Liang & Trelles** has no volume, pages or DOI (the bib note says to confirm).
32. **Krivtsun et al. (2025)** has "The Paton Welding Journal 08" in the journal field; add the volume and issue properly.
33. **LXCat** should name the specific database used for $Q^{(1)}$ (e.g. Phelps, Biagi or Hayashi), with its citation; the access date is not enough.
34. **Korsipati et al. (2026)** is missing pages/article number.
35. **Missing literature.** Please consider citing and discussing:
    - Benilov's reviews of plasma–cathode interaction (e.g. *J. Phys. D* 41 144001, 2008) on solution multiplicity;
    - Riemann (1991) on the sheath edge;
    - experimental anode heat-flux data (Nestor 1962; Sanders & Pfender 1984);
    - Baeva et al. on anode-layer modelling with emitted/separate electron groups.

---

## 4. Assessment

The modelling is detailed and the a-posteriori checks (momentum balance, Boltzmann relation, Bohm criterion) are valuable. But in its present form the paper:

- lacks the numerical description and the convergence evidence (grid and domain length) required for a computational paper;
- reports some electrode quantities that depend on the arbitrary domain length;
- draws its main heat-transfer conclusions from boundary-condition choices that are either crude (cathode conduction law, no surface radiation, $T_{\rm wall}$ above the melting point) or questionable (thoriated-W work function at a water-cooled anode).

Addressing M1–M7 is needed before the paper can be considered further. M8–M13 would substantially strengthen it.

---

## Appendix (for the authors only, not part of the referee report)

Things found while checking the manuscript against the PlasmeX runs and decks, which an external referee could not see:

- **Decks match the paper.** The base-case decks (`runs/base_case_{anode,cathode}/*_inputs`) use `rxn_rates = almeida`, $W = 2.63$ eV, $A_R = 3\times10^4$, and `radial_loss_Tedge = 300` K. $T_{\rm edge} = 300$ K is not stated anywhere in the paper (minor comment 21).
- **Validation parameters are not in the paper.** Both validation decks use $W = 4.55$ eV and $A_R = 1.26\times10^6$. The paper gives $W$ only for the cathode and $A_R$ for neither (M7, minor comment 22). Note that the anode validation also uses $W = 4.55$, while the anode base case uses 2.63 (M6).
- **M3 will be hard to answer from existing runs.** The project's own notes (PlasmeX memory, 9–11 Oct) record that the anode is *not* $L$-independent at 0.4–1 mm, and that the anode relaxation time is $L$-dependent. Answering M3 will need new runs, or a reframing of the anode results.
- **Stale README.** The paper `README.md` names `plt19500` / `plt73000` as the validation plotfiles, but the figures in the manuscript use `plt18500` / `plt38500`.
- **`methods.tex` is orphaned.** It exists but is not `\input`, so the PDF has no numerics section (M1). Its text still describes the older "fourth strategy" framing and should be checked against the present code (Newton/BDF2 vs `Advance`, `flux_scheme = 2`, the tridiagonal composite solve).
- **M11 checked against the runs (14 Oct 2026).** The mechanism is real and the estimate is quantitatively right, but at the base case's quoted precision it changes no number in the paper.
  - *The criterion in the code* (`PlasmeX/src/SteadyState.cpp`): $R_f = \lVert f^{n+1}-f^n\rVert_2 / \lVert f^{n+1}\rVert_2 \cdot t_{\rm ref}/\Delta t$, max over AMR levels, for PHI, NE, ARP, EEN and TGAS. The settings are `steady_tol = 1e-3`, `steady_ref_time = 1e-6`, and `steady_window = 4` checks `steady_check_int = 500` steps apart. AR, HVELX and $T_{\rm wall}$ are not tested. Each L2 norm is taken over a whole level, so on level 0 the bulk dominates it.
  - *Anode base case* `runs/base_case_anode/plt47000`: BDF2 from $t = 0$, $\Delta t = 10^{-9}$ s at the end. The trigger fired at $t = 46.6$ µs $\approx 4.2\,\tau$. The window was $4 \times 500 \times 1$ ns $= 2$ µs $= 0.18\,\tau$, so it added almost nothing. The prediction $\delta \approx \tau \times$ (last-µs rate in `hist.txt`), with $\tau = 11$ µs, matches the change measured when `runs/base_timescales` re-converged the state (backward Euler, `steady_tol = 1e-7`, $\Lambda$ pinned; `anode/A_base/steady/plt47500`):

    | probe | rate at the end (per µs) | predicted $\tau\times$rate | measured on re-convergence |
    |---|---|---|---|
    | ion current at the wall face | $-9.9\times10^{-4}$ | $-1.09\,\%$ | $-1.07\,\%$ |
    | $\phi$ in the wall cell | $-1.2\times10^{-3}$ | $-1.30\,\%$ | $-1.27\,\%$ |
    | $n_e$ at mid-domain | $-7.3\times10^{-4}$ | $-0.80\,\%$ | $-0.78\,\%$ |
    | $T_g$ at mid-domain | $-2.3\times10^{-4}$ | $-0.25\,\%$ | $-0.25\,\%$ |
    | $q_w$ | $-8.6\times10^{-5}$ | $-0.094\,\%$ | $-0.091\,\%$ |
    | electron current at the wall face | $-3.7\times10^{-6}$ | $-0.004\,\%$ | $-0.004\,\%$ |

    The single-mode picture holds to about 3 %. The ion current drifts 270 times faster (relatively) than the electron current, which is the "which quantity" point of M11.2. The table values do not change at the quoted precision: $j_i/j$ goes from $-0.003704$ to $-0.003665$ (both $-0.0037$) and $q_w$ from 38.19 to 38.15 MW m$^{-2}$ (both 38.2). The re-converged state is the one the relaxation-time runs already start from. Replacing the canonical `plt47000` with `plt47500` in the paper is the user's decision. In the anode sweep curves the base case is currently the least-converged point.
  - *Everything else is converged far below the tolerance*, because it ran backward Euler at large $\Delta t$ (the M11.2 integrator point). The anode sweeps ran at $\Delta t = 10^{-6}$ s to 360–440 µs, i.e. 30–40 $\tau$. The cathode base ran to 1 ms with a rate of about $10^{-5}$ per µs at the end. On re-convergence at `1e-7`, every one of these states moved by $\le 2\times10^{-6}$ (`runs/base_timescales/README.md`). The cathode base heat balance closes to $5.4\times10^{-4}$ (`plasmex_info`: $q_w$ 26.964 against conduction 26.978 MW m$^{-2}$); the cathode sweep `tight/` states close to $\le 9\times10^{-4}$.
  - *Unchecked: the validation runs.* `plt18500` (Almeida, $t = 18.0$ µs) and `plt38500` (Semenov, $t = 38.0$ µs) ran BDF2 at $\Delta t \le 10^{-9}$ s under the same trigger (`exe/input_files/*`: `imp_order = 2`, `max_time_step = 1e-9`). Their $\tau$ has not been measured and they have not been re-converged. The Semenov anode case ($\Lambda = 0$) is the one at risk. A backward-Euler restart at `steady_tol=1e-7` from their checkpoints would settle it.
  - *The manuscript has changed since this review.* The quoted sentence ("relative change of every variable over 1 µs ... $10^{-3}$") has been removed from `sections/results.tex`, which now gives only the energy residual and the current uniformity. That removes a criterion that sounded tighter than it was. But the paper now states no criterion at all (M1), and a referee will still ask.
  - *Correction to the first draft of M11.* It said a 1 % error was "larger than" the validation differences ($-0.0068$ vs $-0.0065$; $\phi_s$ within 4 %). It is not: those differences are 4.5 % and 4 %. M11.3 now says so.
- **Arithmetic checks passed.** I recomputed from the stated inputs:
  - the cathode current closure, $e(\Gamma_{\rm em} - \Gamma_{\rm abs} + \Gamma_i) = 1.206\times10^{8}\ \mathrm{A\,m^{-2}}$;
  - all six $q_w$ terms in figure 7;
  - the emission sensitivity (0.26 %/K), the $\pm 20$ K $T_{\rm wall}$ estimates and the 38 $\mathrm{MW\,m^{-2}}$ ion-energy bound;
  - the anode W/A values and the Debye-length scaling of the anode sheath width;
  - the Mason–Schamp constant ($1.11\times10^5$) and the Wannier closed form.

  All are consistent with the text.
