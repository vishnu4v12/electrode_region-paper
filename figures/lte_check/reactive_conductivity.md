# Reactive thermal conductivity used in `kappa_LTE.pdf`

The model has no reactive conductivity coefficient. In an LTE table, the reactive part
$\kappa_R$ describes heat carried by ion pairs: they diffuse down the temperature gradient
and recombine on the cold side, where they release the ionization energy. The model carries
the same heat through its own equations:

- the ambipolar ion flux $\boldsymbol\Gamma_i$;
- the ionization sink $\varepsilon_{\rm ion}S$ in the electron energy equation;
- the electron enthalpy flux $\tfrac52k_BT_e\boldsymbol\Gamma_e$.

To compare the model with an LTE total conductivity, I worked out the flux those terms give
in LTE and wrote it as $-\kappa_R\nabla T$. It is computed in `model()` in `lte_check.py`.
The paper figure (`kappa_LTE.pdf`, Appendix A) uses the Mutation++ $\kappa_R$ instead, so this
estimate is no longer plotted; `lte_check.py` prints it next to the Mutation++ value.

## Assumptions

1. LTE: $T_e=T_h=T$, and $n_e(T)$ follows the Saha composition at $p=1$ atm. I used
   PlasmeX's `argon_lte` fit, through `bohm.n_neutral_lte`.
2. Steady, one-dimensional, quasi-neutral ($n_e=n_i$), at constant pressure, and with the
   heavy gas at rest ($\mathbf u=0$). With equal ion and atom masses, the heavy-gas frame
   condition is
   $$n_i\mathbf u_i+n_n\mathbf u_n=0 .$$
3. No current flows, so the flux is ambipolar: $\boldsymbol\Gamma_e=\boldsymbol\Gamma_i$, i.e.
   $\mathbf u_e=\mathbf u_i$. The electron–ion friction therefore vanishes.
4. Neglected: electron–neutral friction (order $m_e/M$), the electron thermal-diffusion force
   ($C_e$, $C_i$), and the mixture force $Y_i\mathbf F$. That force is zero at steady state,
   since $\nabla p=0$ and the plasma is quasi-neutral.

## Ion flux

Under these assumptions, the electron balance (eq. `mom-e`) reduces to

$$
0=-en_e\mathbf E-\nabla p_e ,
$$

so the ambipolar field passes the electron pressure gradient on to the ions. The ion balance
(eq. `mom-i`) becomes

$$
0=en_i\mathbf E-\nabla p_i-Mn_i\bar\nu^M_{in}(\mathbf u_i-\mathbf u_n).
$$

Add the two, set $n_e=n_i$, $p_e+p_i=2n_ek_BT$, and use $\bar\nu^M_{in}=e/(M\mu_{in})$
(eq. `nu-in`). The ion–neutral slip velocity $\mathbf w=\mathbf u_i-\mathbf u_n$ is then

$$
\mathbf w=-\frac{2\mu_{in}}{n_i}\,\frac{k_B}{e}\,\nabla(n_eT).
$$

Combine this with the frame condition $n_i\mathbf u_i+n_n\mathbf u_n=0$, which gives
$\mathbf u_i=(n_n/n_h)\,\mathbf w$ with $n_h=n_i+n_n$. In LTE, $n_e=n_e(T)$, so

$$
\boldsymbol\Gamma_i=n_i\mathbf u_i
=-\frac{n_n}{n_h}\,2\mu_{in}\,\frac{k_B}{e}\,\frac{d(n_eT)}{dT}\,\nabla T .
$$

The factor $n_n/n_h$ appears because the atoms flow the other way. It makes the flux vanish
at full ionization, where no atoms remain for the ions to diffuse through.

## Energy carried

Three energies travel with each ion pair:

- **Ionization energy.** At steady state $\nabla\cdot\boldsymbol\Gamma_i=S$, so the sink
  $-\varepsilon_{\rm ion}S$ equals $-\nabla\cdot(\varepsilon_{\rm ion}\boldsymbol\Gamma_i)$.
  It therefore acts as a heat flux $\varepsilon_{\rm ion}\boldsymbol\Gamma_i$.
- **Electron enthalpy.** The electron energy flux (eq. `heatflux-e`) contains
  $\tfrac52k_BT\,\boldsymbol\Gamma_e=\tfrac52k_BT\,\boldsymbol\Gamma_i$.
- **Heavy-particle enthalpy: no net contribution.** Ions and atoms carry equal and opposite
  fluxes, so $\tfrac52k_BT(\boldsymbol\Gamma_i+\boldsymbol\Gamma_n)=0$. In any case, the
  heavy energy equation (eq. `henergy`) has no species-diffusion enthalpy term.

So the reactive heat flux is $\mathbf q_R=(\varepsilon_{\rm ion}+\tfrac52k_BT)\,\boldsymbol\Gamma_i\equiv-\kappa_R\nabla T$, which gives

$$
\boxed{\;
\kappa_R=\left(\varepsilon_{\rm ion}+\tfrac52k_BT\right)\frac{n_n}{n_h}\,2\mu_{in}\,\frac{k_B}{e}\,\frac{d(n_eT)}{dT}
\;}
$$

with these inputs:

- $\varepsilon_{\rm ion}=15.7596$ eV;
- $\mu_{in}$, the low-field charge-exchange mobility of eq. `mu_i0`. The Wannier field
  correction vanishes because the LTE field is small;
- $d(n_eT)/dT$, a central difference ($\pm1$ K) of the Saha composition at constant $p$.

## Check: weak-ionization limit

For $n_e\ll n_n\approx p/(k_BT)$, Saha gives $n_e^2/n_n\propto T^{3/2}e^{-\varepsilon_{\rm ion}/k_BT}$,
so $n_e\propto T^{1/4}e^{-\varepsilon_{\rm ion}/2k_BT}$ and

$$
\frac{d(n_eT)}{dT}=n_e\,\frac{\varepsilon_{\rm ion}+\tfrac52k_BT}{2k_BT}.
$$

Substituting, with the ambipolar diffusivity $D_a=2\mu_{in}k_BT/e$, gives

$$
\kappa_R\simeq\frac{n_eD_a}{2}\,\frac{(\Delta H)^2}{k_BT^2},\qquad \Delta H=\varepsilon_{\rm ion}+\tfrac52k_BT .
$$

$\Delta H$ is the enthalpy of the reaction ${\rm Ar}\rightleftharpoons{\rm Ar^+}+e$
(translational part). This has the Butler–Brokaw structure: (reaction enthalpy)$^2$ ×
(diffusion coefficient) × (composition derivative). So the estimate is the standard
reactive conductivity, written with the model's own ion–neutral mobility.

## Comparison with Mutation++

The estimate has the same shape as Mutation++'s reactive conductivity and peaks at the same
temperature, about 14.5 kK. It is 11–22% higher from 8 to 18 kK, and the ratio rises to
1.3–1.4 above 20 kK, where $\kappa_R$ itself is very small. Likely contributors, none of
which I have tested separately:

- **Ion–neutral cross section.** $\mu_{in}$ comes from the Benilov charge-exchange fit
  (eq. `Qia`); Mutation++ uses its own collision integrals. $\kappa_R\propto\mu_{in}$, so a
  10–20% difference in $\overline Q_{in}$ shifts $\kappa_R$ by the same amount. This is
  probably the largest contribution.
- **Neglected terms.** The electron thermal-diffusion heat flux and force, and the
  electron–neutral friction, were dropped.
- **First-order theory.** The binary ion–neutral friction is first-order Chapman–Enskog.
  Mutation++ solves the full multicomponent problem to a higher order.
- **Reaction enthalpy.** $\Delta H$ contains only the ionization energy and the translational
  enthalpy, not the internal (excited-state) enthalpies of Ar and Ar$^+$. Ar$^{2+}$ is
  ignored.

None of this affects the simulation. $\kappa_R$ is only a diagnostic: it shows what the
model's diffusion and source terms carry in LTE, so the model can be compared with LTE totals.
