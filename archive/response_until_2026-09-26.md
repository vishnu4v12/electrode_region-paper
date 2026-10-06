# Eq. (7) vs Eq. (20), and the electron pressure in the heavy-gas momentum equation

> **Update:** Option 1 below has since been adopted in the paper. The $-Y_s\mathbf F$ terms are now in eqs. (3)–(5), so their sum is the identity $0=0$. Eq. (7) now only defines the net force $\mathbf F$. The frequencies in (3)–(5) are written $\bar\nu^M$, and $A_i, A_a$ are renamed $C_i, C_n$. The discussion below describes the equations as they stood before that change.

Equation numbers are those in the compiled `main.pdf`:

- (3)–(5): the species momentum balances (electron, ion, neutral)
- (6): the pairwise momentum-conservation relations
- (7): the bulk force balance, $0 = e(n_i-n_e)\mathbf E - \nabla p$
- (20): the heavy-gas momentum equation (closure 2),
  $\partial_t(\rho\mathbf u) + \nabla\cdot(\rho\mathbf u\mathbf u) = -\nabla p + e(n_i-n_e)\mathbf E$

---

## Question 1

Eq. (7) says that $e(n_i-n_e)\mathbf E - \nabla p = 0$. Eq. (20) says that the same quantity equals the rate of change of heavy-gas momentum, $\partial_t(\rho\mathbf u)+\nabla\cdot(\rho\mathbf u\mathbf u)$. Taken together, (7) would force the left side of (20) to zero, so the heavy gas could never accelerate. Don't the two equations contradict each other?

## Response

They do, if (3)–(5) are read as exact. The two equations come from two different levels of approximation, and the text as written doesn't say so clearly enough.

**Where (7) comes from.** Eq. (7) is not an independent physical law. It is simply the sum of (3)–(5). Each of those balances was written with **all** species inertia dropped, including the inertia of the ions and neutrals. When they are added, the friction terms cancel in pairs by (6) and the thermal-diffusion forces cancel because $C_i + C_n = C_e$. What remains is the statement that the total force on the whole plasma is zero, which is exactly what you get when the whole plasma is assumed to have no inertia. So (7) is the fully inertialess, mechanical-equilibrium limit.

**Where (20) comes from.** Eq. (20) is the same sum, but with the inertia of the heavy mixture put back in. Electron inertia is still negligible, since $m_e/M \sim 10^{-5}$. The total force is no longer zero; it accelerates the heavy gas.

**Why both appear.** The two are used for different jobs, based on a separation of time scales:

- The **relative** motion of each species with respect to the heavy gas, i.e. the diffusion velocities $\mathbf u_s - \mathbf u$, relaxes on the collision time, $1/\nu_{in}$ for ions and $1/\nu_e$ for electrons. That is extremely fast, so for relative motion the inertia really is negligible. This is the only information taken from (3)–(5): they are used to derive the fluxes (10) and (13), which contain only velocity differences plus the frame velocity $\mathbf u$.
- The **common** (frame) motion $\mathbf u$ of the heavy gas relaxes on the acoustic time $L/c_s$, which is much slower. Its inertia is not negligible, so it is governed by (20).

The balances (3)–(5) cannot determine $\mathbf u$ in any case. Their sum (7) contains no velocity, which is exactly why the frame condition (8) and an extra equation for $\mathbf u$ are needed. Closure 2 supplies that equation as (20). Eq. (7) is therefore **not imposed** under closure 2; it is only a diagnostic identity of the inertialess sub-model. Closure 0 is the case where (7) *is* imposed, which is equivalent to dropping the left side of (20). In 1D, (7) can then be integrated using the Maxwell-stress identity $e(n_i-n_e)E = \tfrac{d}{dx}\big(\tfrac12\varepsilon_0E^2\big)$, giving $p - \tfrac12\varepsilon_0E^2 = \text{const}$, the Dalton's-law closure with a field correction.

**The remaining inconsistency, and how to remove it.** Strictly, if the heavy gas accelerates, each heavy species has to take its share of that acceleration. In the exact form, the balances for ions and neutrals carry an inertial term $\rho_s\,D\mathbf u/Dt$ with $\rho_s = Mn_s$. Using (20), this is equivalent to subtracting the species' mass fraction $Y_s = n_s/n_h$ of the net force from each balance:

$$
0 = \mathbf F_s - Y_s\big[e(n_i-n_e)\mathbf E - \nabla p\big] + (\text{friction}) + (\text{thermal diffusion}),
$$

where $\mathbf F_s$ is the species' own electric plus pressure force, and $Y_e \approx 0$. With this correction, the sum of the three balances is $0 = 0$ identically. Eq. (7) then disappears as a constraint, and (3)–(5) and (20) are exactly consistent. This is the standard relative-frame (Stefan–Maxwell) form of multicomponent diffusion.

The equations as written omit the $Y_s$ term. For the ions, its size relative to the ion's own driving force is roughly the ionization fraction $n_i/n_h$, because in the sheath $e(n_i-n_e)E \approx en_iE$. So:

- In the weakly ionized near-electrode layer ($n_i/n_h \sim 1$–$10\%$), the omission is a small, consistent approximation.
- Towards the strongly ionized arc core it is not small.

**Recommendation for the paper.** Either

1. add the $-Y_s[\,e(n_i-n_e)\mathbf E-\nabla p\,]$ term to (4)–(5) and note that (7) is then an identity rather than a constraint (exactly consistent; the ion flux (13) gains one extra term); or
2. keep (3)–(5) as they are, but state explicitly that they govern only the diffusion velocities relative to $\mathbf u$. The paper would then say that (7) is the inertialess limit of (20), not an equation of the model, and that the neglected term is of order $n_i/n_h$.

Option 2 changes nothing in the solver. Option 1 is the rigorous form, but it changes the solver.

---

## Question 2

Eq. (20) is the momentum equation of the heavy gas (ions and neutrals). Why does its pressure $p = n_ek_BT_e + n_hk_BT_h$ include the **electron** pressure, when electrons are not part of the heavy fluid?

## Response

Because electrons are effectively massless, they cannot hold on to momentum. Any force applied to the electrons is passed straight on to the heavy particles.

**From the equations.** The electron balance (3) says

$$
\underbrace{-m_en_e\nu_{en}(\mathbf u_e-\mathbf u_n) - m_en_e\nu_{ei}^{\rm eff}(\mathbf u_e-\mathbf u_i) - C_en_ek_B\nabla T_e}_{\text{force on electrons from collisions with heavy particles}}
= e n_e\mathbf E + \nabla p_e .
$$

By Newton's third law, (6) and $C_i + C_n = C_e$, the heavy particles feel exactly the opposite force:

$$
\mathbf F_{\text{heavy}\leftarrow e} = -\big(en_e\mathbf E + \nabla p_e\big).
$$

Adding this to the forces acting directly on the ions and neutrals, $en_i\mathbf E - \nabla p_i - \nabla p_n$, gives the total force on the heavy gas:

$$
en_i\mathbf E - \nabla(p_i + p_n) - en_e\mathbf E - \nabla p_e
= e(n_i-n_e)\mathbf E - \nabla p ,
$$

which is the right side of (20). The electron pressure gradient and the electric force on the electrons both enter the heavy-gas equation through friction, even though the electrons contribute nothing to $\rho$ on the left side.

**Physical picture.** In a quasi-neutral region, electron pressure pushes the hot, light electrons outward. They cannot run ahead of the ions without creating charge separation, so an ambipolar field builds up that holds the electrons back and pulls the ions forward. The electron pressure is therefore communicated to the heavy gas partly through the field acting on the ions, and partly through direct electron–heavy friction. Either way it ends up as a force on the heavy mixture. This is the same mechanism that makes $T_e$, not $T_h$, set the ambipolar diffusion rate and the ion sound speed $c_s \approx \sqrt{k_B(T_e+T_h)/M}$.

**Why it matters here.** Near the electrode $T_e \gg T_h$, and $n_e$ changes by orders of magnitude over the layer, so $\nabla p_e$ is steep even though $p_e$ itself is small compared with $p_n$. In the quasi-neutral part of the layer, $e(n_i-n_e)\mathbf E \approx 0$, and $\nabla p_e$ is one of the forces that can drive heavy-gas motion there. Inside the sheath, the space-charge force $e(n_i-n_e)\mathbf E$ takes over.

**If $p_e$ were left out of (20),** the model would violate total momentum conservation: the force the electron pressure exerts on the electrons would be removed from the system and never passed to any other species. Leaving it out would only be consistent if the electron momentum balance (3) were also dropped, which is not the case in this model.

---

## Question 3

The momentum-transfer cross section $Q^{(1)}(\varepsilon)$ and the collision frequency $\nu_{eH}(C)=n_nCQ^{(1)}$ are easy to understand. When we need an average over energy, $\bar\nu^M$ is defined so that the momentum exchange can be written as $m_en_e\bar\nu^M(\mathbf u_e-\mathbf u_n)$, and eq. (31) gives its value. In the same way, $\bar\nu^E$ is defined so that the energy exchange takes the form $\tfrac32k_B(T_e-T_h)\,\tfrac{2m_e}{M}\,\bar\nu^E_{en}$, and eq. (32) gives its value. Is that correct? M&K eq. (VIII 2.14a) makes $\bar\nu^M$ obvious. Which equation in M&K supports $\bar\nu^E$?

(Equation numbers are those in the compiled `main.pdf`: (31) is `eq:nu-M` and (32) is `eq:nu-ea`. M&K means Mitchner & Kruger, *Partially Ionized Gases*, Chapter VIII (`Chapt_8_Transport_Properties.pdf`).)

## Response

Yes, that is correct. Each averaged frequency is **defined** by the macroscopic term it has to reproduce. The two averages differ because momentum transfer and energy transfer involve different parts of the electron distribution.

### Where $\bar\nu^E$ is in M&K

M&K define it directly, as **eq. (VIII 2.2a)** on p. 371:

$$
\bar\nu_{eh}=\int_0^\infty\frac{\tfrac12m_eC^2}{\tfrac32kT_e}\,\nu^{(1)}_{eh}(C)\,f^0\,4\pi C^2\,dC .
$$

They call it "the energy-weighted mean collision frequency" and refer to "the equation following (II 6.32)" for where it comes from. Chapter II is not in the PDF we have, so I could not check (II 6.32) itself. The derivation below shows where the definition comes from using only Chapter VIII. With $\langle g\rangle_0=\int g f^0 4\pi C^2dC$ and $\tfrac32kT_e=\tfrac12m_e\langle C^2\rangle_0$, (2.2a) becomes

$$
\bar\nu_{eh}=\frac{\langle C^2\nu_{eh}\rangle_0}{\langle C^2\rangle_0},
$$

which is exactly our eq. (32).

If you want to derive it inside Chapter VIII, use **eq. (VIII 2.15)** on p. 374. Its right-hand side is the elastic electron–heavy collision operator acting on $f^0$:

$$
\sum_h\frac{m_e}{m_h}\frac{1}{C^2}\frac{d}{dC}\Big[C^3\nu^{(1)}_{eh}\Big(f^0+\frac{kT_h}{m_eC}\frac{df^0}{dC}\Big)\Big].
$$

(M&K write (2.15) multiplied through by $C^2$.) Take its energy moment for a Maxwellian $f^0$ at $T_e$. For a Maxwellian, $df^0/dC=-(m_eC/kT_e)f^0$, so the bracket reduces to $C^3\nu f^0(1-T_h/T_e)$. Multiply by $\tfrac12m_eC^2\,4\pi C^2$ and integrate by parts once:

$$
\Big(\frac{d\varepsilon}{dt}\Big)_{\rm el}
=-4\pi m_e\frac{m_e}{M}\Big(1-\frac{T_h}{T_e}\Big)\int_0^\infty C^4\nu f^0\,dC
=-\frac{2m_e}{M}\Big(1-\frac{T_h}{T_e}\Big)\big\langle\tfrac12m_eC^2\nu\big\rangle_0 .
$$

By (2.2a), $\langle\tfrac12m_eC^2\nu\rangle_0=\tfrac32kT_e\,\bar\nu_{eh}$. So

$$
\Big(\frac{d\varepsilon}{dt}\Big)_{\rm el}=-\frac32k_B(T_e-T_h)\,\frac{2m_e}{M}\,\bar\nu^E_{eh}
$$

per electron. Multiplied by $n_e$, this is the loss term in eq. (14) (`eq:eenergy`) and the gain term in eq. (16). This is the energy-exchange analogue of how (2.14a) gives $\bar\nu^M$.

### Where $\bar\nu^M$ comes from, for comparison

In (2.14a), put $-\partial f^0/\partial C=(m_eC/kT_e)f^0$:

$$
\sigma=\frac{4\pi n_ee^2}{3kT_e}\int_0^\infty\frac{C^4f^0}{\nu_{eH}}dC
=\frac{n_ee^2}{m_e}\,\frac{\langle C^2/\nu_{eH}\rangle_0}{\langle C^2\rangle_0}
\equiv\frac{n_ee^2}{m_e\bar\nu^M_{en}} ,
$$

which is eq. (31).

### Why one is a harmonic mean and the other an arithmetic mean

The difference comes from where $\nu$ sits in each calculation.

- **Momentum:** the drift comes from the anisotropic part $f^1$. From (2.11a), each electron's response to the field is proportional to $1/\nu(C)$. Collisions limit the drift, so the current adds up $1/\nu$ over electrons, which gives a harmonic mean. Slow electrons dominate it, and in argon those are the ones trapped in the Ramsauer minimum.
- **Energy:** elastic energy transfer comes from the isotropic part $f^0$. Each electron loses a fraction $\approx2m_e/M$ of its energy per collision, and it collides at rate $\nu(C)$. The loss therefore adds up $\varepsilon\,\nu$ over electrons, which gives an arithmetic mean weighted by energy. Fast electrons dominate it.

Both averages use the same weight, $C^2f^0$. So $\bar\nu^M\le\bar\nu^E$, with equality only when $\nu_{eH}$ does not depend on energy. For argon the ratio is 4–4.5.

### A caveat on "defined so that the friction is $m_en_e\bar\nu^M(\mathbf u_e-\mathbf u_n)$"

In short: $\bar\nu^M$ is fixed by how fast the electrons **drift** for a given force, not by how much momentum they **exchange**. The friction $m_en_e\bar\nu^M(\mathbf u_e-\mathbf u_n)$ gets the momentum exchange right when the electrons are pushed by $\mathbf E$ or $\nabla p_e$. It gets it wrong when they are pushed by $\nabla T_e$, and the term $C_en_ek_B\nabla T_e$ in eq. (3) is the correction for that case. The true momentum exchange is the **sum** of the two terms.

The argument has four steps.

**Step 1: the friction is fixed by the force, not by $\nu$.** Electron inertia is negligible, so the electron momentum balance is just "collisional momentum transfer + applied force = 0". Call the momentum the neutrals take from the electrons per unit volume and time $\mathbf R$. Then

$$
\mathbf R=-en_e\mathbf E-\nabla p_e
$$

exactly, whatever the cross section. The collision frequency does not decide how large the friction is; the force does. What the collision frequency decides is **how fast the electrons have to drift** to produce that friction. So the only meaningful way to define an averaged frequency is through the ratio

$$
\nu_{\rm eff}=\frac{|\mathbf R|}{m_en_e|\mathbf u_e-\mathbf u_n|}
=\frac{\text{momentum exchanged}}{m_en_e\times\text{drift}} .
$$

A single number $\bar\nu$ can only work if this ratio is the same for every way of driving the electrons.

**Step 2: a two-group picture shows why the average is harmonic.** Take two groups of electrons with equal numbers, a slow group with frequency $\nu_s$ and a fast group with $\nu_f$. Push every electron with the same force $\mathbf F$. Each group reaches its own terminal drift, $\mathbf u_k=\mathbf F/(m_e\nu_k)$. The mean drift is

$$
\bar{\mathbf u}=\frac{\mathbf F}{m_e}\cdot\tfrac12\Big(\frac1{\nu_s}+\frac1{\nu_f}\Big)
\quad\Rightarrow\quad
\nu_{\rm eff}=\frac{1}{\tfrac12(1/\nu_s+1/\nu_f)} ,
$$

the harmonic mean. The group that collides less drifts faster, so it carries most of the flux. That is the origin of eq. (31), $\bar\nu^M=\langle C^2\rangle_0/\langle C^2/\nu\rangle_0$. The kinetic theory adds the weight $C^2f^0$, but the structure is the same.

**Step 3: the kinetic result, written out.** With $\partial f_M/\partial C=-(m_eC/k_BT_e)f_M$, all three driving terms in eq. (`eq:lorentz-f1`) carry the same factor $(C/\nu)f_M$:

$$
f^1=-\frac{C}{\nu_{eH}}f_M\Big[\frac{e\mathbf E}{k_BT_e}+\nabla\ln p_e+\Big(x-\tfrac52\Big)\nabla\ln T_e\Big].
$$

The drift is $\mathbf u_e-\mathbf u_n=\tfrac{4\pi}{3}\int_0^\infty C^3f^1\,dC$ (M&K eq. VIII 2.12). The integral that appears is $\tfrac{4\pi}{3}\int C^4f_M/\nu\,dC=\tfrac13\langle C^2/\nu\rangle_0=k_BT_e/(m_e\bar\nu^M)$, so

$$
\mathbf u_e-\mathbf u_n=
-\frac{1}{m_e\bar\nu^M}\Big(e\mathbf E+\frac{\nabla p_e}{n_e}\Big)
-\frac{k_B}{m_e\bar\nu^M}\Big(\langle x\rangle_\nu-\tfrac52\Big)\nabla T_e .
$$

Multiplying by $m_en_e\bar\nu^M$ and rearranging gives

$$
0=-en_e\mathbf E-\nabla p_e-m_en_e\bar\nu^M(\mathbf u_e-\mathbf u_n)-C_en_ek_B\nabla T_e,
\qquad C_e=\langle x\rangle_\nu-\tfrac52 ,
$$

which is exactly eq. (3). So eq. (3) is the kinetic solution rewritten to look like a force balance, not a friction model that was assumed.

For $\mathbf E$ and $\nabla p_e$, the shape of $f^1$ is $(C/\nu)f_M$, and the ratio of Step 1 is exactly $\bar\nu^M$. That is also why the Einstein relation $D_e/\mu_e=k_BT_e/e$ holds with the same $\bar\nu^M$: mobility and diffusion see the same shape of $f^1$.

**Step 4: $\nabla T_e$ drives a flux with no net friction.** Consider $\mathbf E=0$ and uniform $p_e$, so that only $\nabla T_e$ acts. By Step 1, the net force on the electrons is zero, so the momentum they exchange with the neutrals is exactly zero, $\mathbf R=0$. Yet Step 3 says they still drift:

$$
\mathbf u_e-\mathbf u_n=-\frac{k_BC_e}{m_e\bar\nu^M}\nabla T_e\neq0 .
$$

The two-group picture explains this. The factor $(x-\tfrac52)$ in $f^1$ has opposite signs for slow electrons ($x<\tfrac52$) and fast electrons ($x>\tfrac52$). $\nabla T_e$ therefore pushes the two groups in **opposite** directions, with forces that cancel, $\mathbf F_s=-\mathbf F_f$, so there is no net friction. Their drifts are $\mathbf F_s/(m_e\nu_s)$ and $-\mathbf F_s/(m_e\nu_f)$, which cancel only if $\nu_s=\nu_f$. When $\nu$ depends on energy, the group that collides less wins, and a net flux appears with zero net momentum exchange. In argon the slow electrons sit in the Ramsauer minimum, collide less, and dominate.

A pure friction law $m_en_e\nu(\mathbf u_e-\mathbf u_n)$ cannot represent this for **any** choice of $\nu$: zero friction would force zero drift. This is why eq. (3) needs the extra term $C_en_ek_B\nabla T_e$. In the $\nabla T_e$-only case, eq. (3) reads $0=-m_en_e\bar\nu^M(\mathbf u_e-\mathbf u_n)-C_en_ek_B\nabla T_e$. The two terms cancel, so their sum, which is the real momentum exchange, is zero as it must be. Neither term alone is the physical friction.

**Summary.**

| | your statement | more precise statement |
|---|---|---|
| what $\bar\nu^M$ reproduces | the friction $m_en_e\bar\nu^M(\mathbf u_e-\mathbf u_n)$ | the drift per unit force (mobility), for $\mathbf E$- and $\nabla p_e$-driven motion |
| true momentum exchange | $m_en_e\bar\nu^M(\mathbf u_e-\mathbf u_n)$ | $m_en_e\bar\nu^M(\mathbf u_e-\mathbf u_n)+C_en_ek_B\nabla T_e$ |

**Check with a constant $\nu$.** If $\nu$ does not depend on energy, the weight $C^4f_M/\nu$ is proportional to $C^4f_M$, which gives $\langle x\rangle_\nu=\tfrac52$ and $C_e=0$. At the same time $\bar\nu^M=\bar\nu^E=\nu$. The caveat disappears, and your statement is then exactly right. Both the thermal force and the $M$/$E$ distinction come from the energy dependence of $Q^{(1)}$.

### Suggestion for the paper

The appendix currently asserts that eq. (32) is the energy-exchange average without a source. You could cite M&K eq. (VIII 2.2a) at eq. (32), and eq. (VIII 2.14a) at eq. (31). In the paragraph that states $\tfrac{2m_e}{M}\langle\varepsilon\nu_{eH}\rangle_0(1-T_h/T_e)=\tfrac32k_B(T_e-T_h)\tfrac{2m_e}{M}\bar\nu^E_{en}$, you could add "obtained as the energy moment of the elastic collision operator (M&K eq. VIII 2.15)".

---

## Question 4

Where has each collision frequency been calculated from?

- $\bar\nu^M_{en}$: from eq. (31), using the cross-section data as a function of energy? If so, write eq. (31) out with $\nu_{eH}$ expanded, in the paper too.
- $\bar\nu^M_{ei}$: from M&K (5.1a) and (5.1d). Has it been checked against the eq. (31) form of averaging?
- $\bar\nu^M_{ie}$: from the momentum balance.
- $\bar\nu^M_{in}$: from the mobility relation.
- $\bar\nu^E_{ei}$: from M&K (5.1d), which I believe is the averaging of eq. (32).
- $\bar\nu_{ee}$: from M&K (2.2d), which I believe is the averaging of eq. (32).

(Equation numbers are those in the compiled `main.pdf`. M&K numbers refer to Chapter VIII. The code references are to the PlasmeX solver, `../PlasmeX`.)

## Response

Your list is right for all six, with one important correction for $\bar\nu^M_{ei}$ and one clarification for $\bar\nu_{ee}$. The summary comes first, then the details.

| frequency | where it comes from | type of average |
|---|---|---|
| $\bar\nu^M_{en}$ | eq. (31) on the LXCat argon $Q^{(1)}(\varepsilon)$ | harmonic ($M$), Lorentz gas |
| $\bar\nu^E_{en}$ | eq. (32) on the same $Q^{(1)}(\varepsilon)$ | arithmetic ($E$) |
| $\bar\nu^E_{ei}$ | M&K (5.1d) = eq. (35) | arithmetic ($E$): eq. (32) applied to the Coulomb $Q^{(1)}$ |
| $\bar\nu^M_{ei}$ | $0.5\,\bar\nu^E_{ei}$, from M&K (5.1a) | Spitzer–Härm, **not** eq. (31) |
| $\bar\nu^M_{ie}$ | momentum conservation, eq. (6) | inherited from $\bar\nu^M_{ei}$ |
| $\bar\nu^M_{in}$ | $e/(M\mu_{in})$, eq. (47), with the Wannier mobility | defined by the mobility ($M$) |
| $\bar\nu^M_{ni},\ \bar\nu^M_{ne}$ | momentum conservation, eq. (6) | inherited |
| $\bar\nu^E_{eh}$ | $\bar\nu^E_{en}+\bar\nu^E_{ei}$, eq. (42) | arithmetic ($E$) |
| $\bar\nu_{ee}$ | not used in the model equations; appears only in the scaling estimates, from the NRL formulary | M&K (2.2b) would be the $E$-type average |

### 1. $\bar\nu^M_{en}$: yes, eq. (31) on the cross-section data

The input is the argon momentum-transfer cross section $Q^{(1)}(\varepsilon)$ from LXCat/BOLSIG+ (`PlasmeX/exe/ar_Q_en.dat`). With $\nu_{eH}(C)=n_nCQ^{(1)}(\varepsilon)$ and $\langle C^2\rangle_0=3k_BT_e/m_e$, eq. (31) becomes

$$
\bar\nu^M_{en}=\frac{\langle C^2\rangle_0}{\langle C^2/\nu_{eH}\rangle_0}
=n_n\,\frac{3k_BT_e/m_e}{4\pi\int_0^\infty C^3f_M\,[Q^{(1)}]^{-1}\,dC}.
$$

Changing variable to $\varepsilon=\tfrac12m_eC^2$ (so $C^3dC=2\varepsilon\,d\varepsilon/m_e^2$) and using $f_M=(m_e/2\pi k_BT_e)^{3/2}e^{-\varepsilon/k_BT_e}$ gives

$$
\bar\nu^M_{en}=n_n\sqrt{\frac{2}{m_e}}\;\frac{3\sqrt\pi}{4}\;
\frac{(k_BT_e)^{5/2}}{\displaystyle\int_0^\infty\varepsilon\,e^{-\varepsilon/k_BT_e}\,\frac{d\varepsilon}{Q^{(1)}(\varepsilon)}} .
$$

The same steps on eq. (32) give, for comparison,

$$
\bar\nu^E_{en}=n_n\sqrt{\frac{2}{m_e}}\;\frac{4}{3\sqrt\pi}\;
\frac{\displaystyle\int_0^\infty\varepsilon^2e^{-\varepsilon/k_BT_e}\,Q^{(1)}(\varepsilon)\,d\varepsilon}{(k_BT_e)^{5/2}} .
$$

The contrast is clear in these forms. $\bar\nu^M$ integrates $1/Q$ under the weight $\varepsilon e^{-\varepsilon/kT}$, so it is dominated by low energies, where argon's $Q^{(1)}$ has its Ramsauer minimum. $\bar\nu^E$ integrates $Q$ under the weight $\varepsilon^2e^{-\varepsilon/kT}$, so it is dominated by high energies.

**Check with a constant cross section.** With $Q^{(1)}=Q$, the integrals give $\bar\nu^E=\tfrac43n_n\bar c_eQ$ and $\bar\nu^M=\tfrac{9\pi}{32}\cdot\tfrac43n_n\bar c_eQ\approx0.88\,\bar\nu^E$. They are still not equal, because a constant $Q$ means $\nu\propto C$, which is not a constant $\nu$.

**How the code does it.** `gen_lorentz_tables.py` does not evaluate eq. (31) literally. It evaluates M&K (2.14a) for $\sigma_{ea}$ on the same $Q^{(1)}$ and tabulates it (`ar_e_lorentz.dat`, column 2). The solver then forms $\bar\nu^M_{en}=e/(m_e\mu_{en})=n_ee^2/(m_e\sigma_{ea})$. As shown in Question 3, putting a Maxwellian into (2.14a) gives exactly eq. (31), so the two are the same number. The script also checks $\sigma_{ea}/(n_ee)$ against the BOLSIG+ mobility `emob.dat`. They agree to within 1% for $T_e\le2$ eV (6% at 0.1 eV), and differ by 7–23% above 3 eV, where BOLSIG+ includes inelastic channels.

**Change to the paper.** Both expanded forms of eq. (31) above are now in `sections/appendix.tex`, followed by a sentence explaining the steps, the constant-$Q$ limit, and the equivalence with M&K (2.14a). The paper compiles and the equation numbers are unchanged. I did not add the energy form of eq. (32); it is given above if you want it for symmetry.

### 2. $\bar\nu^E_{ei}$: yes, M&K (5.1d) is the eq. (32) average

M&K (5.1d) gives

$$
\bar\nu_{ei}=n_i\frac{4\sqrt{2\pi}}{3}\Big(\frac{m_e}{k_BT_e}\Big)^{3/2}\Big(\frac{e^2}{4\pi\varepsilon_0m_e}\Big)^2\ln\Lambda ,
$$

which is the paper's eq. (35). It is exactly the eq. (32) / M&K (2.2a) energy-weighted average of the Coulomb momentum-transfer frequency. To show this, take the Rutherford cross section $Q^{(1)}_{ei}=4\pi b_0^2\ln\Lambda$ with $b_0=e^2/(4\pi\varepsilon_0m_eC^2)$. This gives $\nu_{ei}(C)=n_iCQ^{(1)}_{ei}=4\pi n_iK\ln\Lambda/C^3$, with $K=(e^2/4\pi\varepsilon_0m_e)^2$. With $a=\sqrt{2k_BT_e/m_e}$, the Maxwellian moments are $\langle C^{-1}\rangle_0=2/(\sqrt\pi a)$ and $\langle C^2\rangle_0=\tfrac32a^2$. Then

$$
\frac{\langle C^2\nu_{ei}\rangle_0}{\langle C^2\rangle_0}
=4\pi n_iK\ln\Lambda\,\frac{\langle C^{-1}\rangle_0}{\langle C^2\rangle_0}
=4\pi n_iK\ln\Lambda\,\frac{4}{3\sqrt\pi a^3}
=n_i\frac{4\sqrt{2\pi}}{3}\Big(\frac{m_e}{k_BT_e}\Big)^{3/2}K\ln\Lambda ,
$$

which is (5.1d). The same check is done numerically in `gen_lorentz_tables.py`: the column-1 integral on the Coulomb $Q^{(1)}$ reproduces (5.1d) to within $10^{-3}$. The numerical prefactor in eq. (35), $2.906\times10^{-12}$ (with $n_i$ in m$^{-3}$ and $T_e$ in eV), is also correct.

### 3. $\bar\nu^M_{ei}$: from (5.1a) and (5.1d). Eq. (31) would give the wrong value

You are right about the source. M&K (5.1a) gives the Spitzer–Härm conductivity $\sigma=1.975\,n_ee^2/(m_e\bar\nu_{ei})$. Defining $\bar\nu^M_{ei}$ through $\sigma=n_ee^2/(m_e\bar\nu^M_{ei})$, the same definition as eq. (`eq:onsager-map`), gives

$$
\bar\nu^M_{ei}=\frac{\bar\nu^E_{ei}}{1.975}=0.506\,\bar\nu^E_{ei}.
$$

The paper and the code (`src/coulomb.H`, `TransportFactor`) round this to 0.5, a difference of 1.2%.

**Checked against eq. (31)?** It has been checked in the sense that the eq. (31) machinery, applied to the Coulomb cross section, reproduces M&K's own Lorentz-gas number. But the value used is deliberately **not** the eq. (31) value, and it must not be. Applying eq. (31) to $\nu_{ei}\propto C^{-3}$ gives

$$
\bar\nu^M_{ei}\big|_{\rm eq.\,31}
=\frac{\langle C^2\rangle_0^2}{\langle C^5\rangle_0\langle C^{-1}\rangle_0}\,\bar\nu^E_{ei}
=\frac{(\tfrac32a^2)^2}{(12a^5/\sqrt\pi)(2/\sqrt\pi a)}\,\bar\nu^E_{ei}
=\frac{3\pi}{32}\,\bar\nu^E_{ei}=0.295\,\bar\nu^E_{ei},
$$

that is, a conductivity $\sigma=3.40\,n_ee^2/(m_e\bar\nu_{ei})$. This is the Lorentzian value in M&K Table 1 (p. 413): $\sigma_{\rm Lorentz}/\sigma_{\rm exact}=1.72$, and $1.72\times1.975=3.40$. `gen_lorentz_tables.py` checks this as its first acceptance test ($\sigma_L=1.975/0.582=3.39$).

**Why eq. (31) fails here.** Eq. (31) is the Lorentz-gas result, and a Lorentz gas neglects electron–electron collisions. That is justified for electron–neutral collisions in a weakly ionized gas. It is not justified for electron–ion collisions, because for singly charged ions $\bar\nu_{ee}\approx\sqrt2\,\bar\nu_{ei}$ (see item 6). For Coulomb collisions $\nu_{ei}\propto C^{-3}$, so the harmonic mean is dominated by the fast electrons, which barely collide with ions. In a Lorentz gas they would carry most of the current. In a real plasma, electron–electron collisions keep scattering those fast electrons and remove that advantage. The conductivity drops from $3.40$ to $1.975$ in units of $n_ee^2/(m_e\bar\nu_{ei})$; the ratio $0.582$ in (5.1a) is the Spitzer correction factor. So the correct $\bar\nu^M_{ei}$ is $0.506\,\bar\nu^E_{ei}$ (Spitzer–Härm), not $0.295\,\bar\nu^E_{ei}$ (Lorentz, eq. 31).

A side note: the first Sonine approximation in Table 1 gives $\sigma^{(1)}/\sigma_{\rm exact}=0.506$, i.e. $\sigma^{(1)}=n_ee^2/(m_e\bar\nu_{ei})$. The lowest-order (displaced-Maxwellian) friction therefore uses $\bar\nu^E_{ei}$ itself. The factor $1/1.975$ is the correction from the higher Sonine terms.

**Suggestion for the paper.** At eq. (38) (`eq:coulomb`) the paper says only that $0.5\approx1/1.975$ is the Spitzer conductivity weight. You could add one sentence: "This is not the Lorentz average (31) applied to the Coulomb cross section, which would give $3\pi/32\approx0.29$, because electron–electron collisions, negligible for electron–neutral transport, are as frequent as electron–ion collisions in the Coulomb limit (M&K Table 1)."

### 4. $\bar\nu^M_{ie}$: yes, from momentum conservation

Eq. (48) comes from the second relation in eq. (6), $m_en_e\bar\nu^M_{ei}=Mn_i\bar\nu^M_{ie}$:

$$
\bar\nu^M_{ie}=\frac{m_en_e}{Mn_i}\,\bar\nu^M_{ei}.
$$

This is not an independent collision frequency. It is defined so that the electron–ion friction in eqs. (3) and (4) is equal and opposite, and it carries the same Spitzer factor 0.5. Its size is $\sim(m_e/M)\,\bar\nu^M_{ei}\approx1.4\times10^{-5}\,\bar\nu^M_{ei}$ (for $n_e\approx n_i$). The code builds the equivalent ratio $r=\bar\nu^M_{ie}/\bar\nu^M_{in}=(1-\text{coul\_fac})\,\mu_{in}n_e/(|\mu_e|n_i)$ in `src/aider.H`, in which the masses cancel. The same applies to $\bar\nu^M_{ne}=m_en_e\bar\nu^M_{en}/(Mn_n)$ and $\bar\nu^M_{ni}=n_i\bar\nu^M_{in}/n_n$, the first and third relations of eq. (6).

### 5. $\bar\nu^M_{in}$: yes, from the mobility relation

The chain is eq. (43) → (44) → (45)–(46) → (47):

1. The charge-exchange cross section $\overline Q_{ia}(T_h)$ comes from the Benilov fit, eq. (43).
2. The low-field mobility is $\mu_i^0=1.104\times10^5/(\overline Q_{ia}\sqrt{T_h}\,n_n)$, eq. (44). This agrees to within 1% with the first-approximation Chapman–Enskog mobility for a constant cross section, $3e\sqrt\pi/(8n_n\overline Q\sqrt{Mk_BT_h})=1.113\times10^5/(n_n\overline Q\sqrt{T_h})$.
3. The Wannier correction for field heating, $\mu_{in}=\mu_i^0\sqrt{T_h/T_{\rm eff}}$, is eqs. (45)–(46).
4. Finally $\bar\nu^M_{in}=e/(M\mu_{in})$, eq. (47).

It is an $M$-type frequency by construction, since it is defined to reproduce the ion mobility. That is the ion analogue of $\bar\nu^M_{en}=e/(m_e\mu_{en})$. No separate velocity average is taken, because the mobility data already contain it.

### 6. $\bar\nu_{ee}$: M&K (2.2b), not "(2.2d)", and it is not used in the model equations

There is no (2.2d) in M&K Chapter VIII; p. 371 has only (2.2a) and (2.2b). You presumably mean **(2.2b)**:

$$
\bar\nu_{ee}=n_e\frac{8\sqrt\pi}{3}\Big(\frac{m_e}{k_BT_e}\Big)^{3/2}\Big(\frac{e^2}{4\pi\varepsilon_0m_e}\Big)^2\ln\Lambda
=4.110\times10^{-12}\,\frac{n_e[\mathrm{m^{-3}}]\ln\Lambda}{T_e[\mathrm{eV}]^{3/2}}\ \mathrm{s^{-1}} .
$$

You are right that it is the same kind of average as eq. (32): M&K introduce (2.2a) and (2.2b) together as "energy-weighted mean collision frequencies". The difference is that for an electron–electron pair the average is over the **relative** motion, with reduced mass $m_e/2$ and relative speed $g$. The Rutherford cross section then has $b_0=e^2/(4\pi\varepsilon_0\,\tfrac{m_e}{2}g^2)$, which is twice as large at the same speed, so the cross section is 4 times larger, and the relative-speed Maxwellian is $\sqrt2$ wider. Repeating the item-2 calculation with these changes gives

$$
\bar\nu_{ee}=\sqrt2\,\bar\nu_{ei}\big|_{n_i\to n_e},
$$

and $8\sqrt\pi/(4\sqrt{2\pi})=\sqrt2$ confirms it.

**Where it enters the paper.** $\bar\nu_{ee}$ does **not** appear in any model equation. It appears only in the scaling section: the Maxwellian-validity condition, eq. (52) (`eq:maxw1`), and the mean free path, eq. (`eq:mfp`). There, the value $\nu_{ee}\sim\nu_{ei}\sim5\times10^{11}$ s$^{-1}$ is quoted from the NRL formulary with $\ln\Lambda\approx3$, not from M&K (2.2b). The two are consistent for an order-of-magnitude estimate: at $T_e=1$ eV, $n_e=5\times10^{22}$ m$^{-3}$ and $\ln\Lambda=3$, (5.1d) gives $4.4\times10^{11}$ and (2.2b) gives $6.2\times10^{11}$ s$^{-1}$. In the model itself, electron–electron collisions enter only implicitly, through the Spitzer–Härm factors 1.975, 3.203 and 1.389 in eq. (34), which is exactly why $\bar\nu^M_{ei}\neq$ eq. (31) in item 3.

**Minor consistency note.** In the same scaling paragraph, the electron–neutral estimate $\nu_{en}=\bar c_eN\bar Q_{ea}\approx8\times10^9$ s$^{-1}$ is neither $\bar\nu^M$ nor $\bar\nu^E$. With a constant $Q$, $\bar\nu^E=\tfrac43\bar c_eNQ$ and $\bar\nu^M\approx0.88\,\bar\nu^E$. This does not affect the conclusion, since the condition is met by a factor of $10^6$. If you want the paper to be consistent throughout, you could cite M&K (2.2b) for $\bar\nu_{ee}$ and note that $\bar\nu_{eh}$ in eq. (52) is the $E$-type average of M&K (2.2a).
