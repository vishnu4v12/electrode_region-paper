# Appendix 3.1 and 3.2 explained from the basics

This note explains Sections 3.1 ("Electron transport: general framework") and 3.2 ("Partially ionized limit: electron–neutral transport") of the paper's Appendix A, `sections/appendix.tex`. It assumes only basic mechanics, the ideal gas law, and calculus. Equation numbers such as (27) refer to the compiled `main.pdf`.

Contents

0. Why this appendix exists
1. Background you need
2. Section 3.1: two forces, two fluxes (the Onsager framework)
3. Section 3.2: computing the coefficients when electrons hit neutral atoms
4. Numbers for argon
5. Summary: what goes where in the model
6. Notes for the authors (possible corrections)
7. Symbol table

---

## 0. Why this appendix exists

The model in Section 2 solves conservation equations for electrons, ions, and atoms. Two of them concern the electrons:

* the **electron flux** $\boldsymbol\Gamma_e$ (how many electrons per m² per s cross a surface), eq. (9);
* the **electron energy equation**, eq. (13), with its energy flux $\mathbf h_e$, eq. (14).

These equations contain coefficients the model cannot produce by itself:

| coefficient | appears in | physical meaning |
|---|---|---|
| $\bar\nu^M_e$ | electron flux (9) | how often an electron's momentum is randomized ("friction") |
| $C_e$ | electron flux (9), heat flux (14) | extra push on electrons from a temperature gradient (thermal diffusion) |
| $\kappa_e$ | heat flux (14) | electron thermal conductivity |
| $\bar\nu^E_{eh}$ | energy equation (13) | how fast electrons hand their energy to the heavy particles |

The appendix supplies them. Section 3.1 explains **which** combinations of microscopic quantities these coefficients must be. Section 3.2 **computes** them in the case where electrons mostly collide with neutral argon atoms.

---

## 1. Background you need

### 1.1 The electron distribution function

There are about $10^{22}$ electrons per m³, far too many to follow one by one. Instead we use a **distribution function** $f(\mathbf r,\mathbf C)$:

$$
f(\mathbf r,\mathbf C)\,d^3r\,d^3C = \text{number of electrons in the small volume } d^3r \text{ whose velocity is within } d^3C \text{ of } \mathbf C .
$$

Here $\mathbf C$ is the electron velocity and $C=|\mathbf C|$ its speed. Everything macroscopic is an average, or "moment", of $f$:

$$
n_e=\int f\,d^3C,\qquad
\boldsymbol\Gamma_e=\int \mathbf C\,f\,d^3C,\qquad
\mathbf J_e=-e\,\boldsymbol\Gamma_e,\qquad
\mathbf q_e=\int \tfrac12 m_eC^2\,\mathbf C\,f\,d^3C .
$$

These are the density, the particle flux, the electric current (the minus sign is because the electron charge is $-e$), and the energy flux (kinetic energy carried across a surface per unit area per unit time).

### 1.2 The Maxwellian and the variable $x$

If electrons collide often enough with each other and with everything else, $f$ settles to the **Maxwellian**:

$$
f_M(C)=n_e\left(\frac{m_e}{2\pi k_BT_e}\right)^{3/2}\exp\!\left(-\frac{m_eC^2}{2k_BT_e}\right).
$$

It depends only on the speed $C$, not on direction. This is **isotropic**: as many electrons move left as right, so the Maxwellian carries **no current and no heat flux**. All transport comes from a small *deviation* from it (Section 3).

Two numbers recur throughout:

* the dimensionless energy $x=\dfrac{m_eC^2}{2k_BT_e}=\dfrac{\varepsilon}{k_BT_e}$, which says how energetic an electron is compared with a typical one;
* the averages $\langle C^2\rangle_0=3k_BT_e/m_e$ (the mean energy is $\tfrac32k_BT_e$) and the mean speed $\bar c_e=\sqrt{8k_BT_e/\pi m_e}$.

At $T_e=1$ eV (about 11 600 K), $\bar c_e\approx 6.7\times10^5$ m/s.

### 1.3 Cross sections and collision frequency

Picture a neutral atom as a target of area $Q$ (m²). An electron moving at speed $C$ through atoms of density $n_n$ sweeps out a volume $CQ$ per second, so it hits atoms at the rate

$$
\nu_{eH}(C)=n_n\,C\,Q^{(1)}(\varepsilon)\quad[\text{s}^{-1}].
$$

The subscript "eH" means "electron–heavy particle". The superscript $(1)$ marks the **momentum-transfer cross section**. A glancing collision barely deflects the electron, so it should count less than a head-on bounce. $Q^{(1)}$ weights each scattering angle $\chi$ by $(1-\cos\chi)$, the fraction of forward momentum lost. With this weighting, $\nu_{eH}$ is "the rate at which the electron forgets its direction".

**Key point:** $Q^{(1)}$ depends on the electron energy $\varepsilon$, so $\nu_{eH}$ is different for slow and fast electrons. That single fact drives everything in Section 3.2.

### 1.4 Why the electron mass ratio matters

An electron is about 73 000 times lighter than an argon atom. When it bounces off an atom it is like a ping-pong ball hitting a bowling ball:

* its **direction** (momentum) is randomized in roughly **one** collision;
* its **energy** hardly changes. It loses only a fraction of about $2m_e/M\approx2.7\times10^{-5}$ per collision.

So the electron velocity distribution quickly becomes nearly **isotropic** (direction is forgotten fast), while its energy distribution is set by slower processes. This separation is what makes the "two-term" approximation of Section 3.2 work.

### 1.5 Argon is unusual: the Ramsauer minimum

For argon, $Q^{(1)}$ is not constant. It drops to a deep minimum near 0.23 eV, where electrons pass through atoms almost freely (a quantum interference effect), and rises steeply above it. From the cross-section table the code uses (`~/work/PlasmeX/exe/ar_Q_en.dat`, LXCat):

| $\varepsilon$ [eV] | 0.1 | **0.23** | 0.5 | 1 | 2 | 5 | 10 |
|---|---|---|---|---|---|---|---|
| $Q^{(1)}$ [$10^{-20}$ m²] | 0.61 | **0.098** | 0.48 | 1.39 | 2.70 | 7.46 | 15.3 |

Between 0.23 eV and 10 eV, $Q^{(1)}$ rises about 150-fold. Fast electrons therefore collide far more often than slow ones, and different averages of $\nu_{eH}$ give very different numbers (Section 3.6).

### 1.6 Why "averaged" collision frequencies are needed

The fluid equations have one friction coefficient per species. The microscopic picture has a whole function $\nu_{eH}(C)$. To turn a function into one number we must average over the electrons, but the right average depends on what is being computed. For example:

* the **current** is carried mostly by electrons that travel far between collisions;
* the **energy loss** is dominated by electrons that collide often.

Each question weights the electrons differently, so each needs its own average. That is the origin of the superscripts $M$ (momentum or mobility) and $E$ (energy exchange).

---

## 2. Section 3.1: two forces, two fluxes

### 2.1 What pushes electrons: the generalized field

Two things push a fluid element of electrons:

1. the electric field: force $-e\mathbf E$ per electron;
2. the pressure gradient: force $-\nabla p_e$ per unit volume, that is $-\nabla p_e/n_e$ per electron, with $p_e=n_ek_BT_e$.

Together the force per electron is

$$
-e\mathbf E-\frac{\nabla p_e}{n_e}=-e\left(\mathbf E+\frac{\nabla p_e}{en_e}\right)\equiv-e\,\boldsymbol{\mathcal E}.
$$

So the electrons respond to the **generalized field** $\boldsymbol{\mathcal E}=\mathbf E+\nabla p_e/(en_e)$ exactly as they would to a pure electric field. A pressure gradient is equivalent to an electric field of size $\nabla p_e/(en_e)$. This is why drift (from $\mathbf E$) and diffusion (from $\nabla n_e$) always come with the same coefficient: the Einstein relation.

The **temperature gradient** $\nabla T_e$ is a second, genuinely different driving force, as Section 2.3 shows.

### 2.2 Linear response

When the driving forces are small (the usual situation in a collisional plasma), each flux is a **linear** combination of the forces. With two forces, $\boldsymbol{\mathcal E}$ and $\nabla T_e$, and two fluxes, the current $\mathbf J_e$ and the energy flux $\mathbf q_e$, we need a 2×2 table of coefficients. Eq. (27)–(28):

$$
\boxed{\mathbf J_e=\sigma\,\boldsymbol{\mathcal E}+\varphi\,\nabla T_e}
$$

$$
\boxed{\mathbf q_e=-\frac52\frac{k_BT_e}{e}\mathbf J_e-\lambda'\nabla T_e-T_e\varphi\,\boldsymbol{\mathcal E}}
$$

Reading each term:

* $\sigma\boldsymbol{\mathcal E}$ is ordinary **Ohm's law**. $\sigma$ is the electrical conductivity [S/m].
* $\varphi\nabla T_e$ is a **current driven by a temperature gradient**, even with no field (the thermoelectric effect). Electrons on the hot side are faster and diffuse differently from those on the cold side, so a net flow appears. $\varphi$ is the thermal-diffusion coefficient.
* $-\frac52\frac{k_BT_e}{e}\mathbf J_e$ is energy carried **along with** the moving electrons, at $\tfrac52k_BT_e$ per electron (see 2.4). Since $\mathbf J_e=-e\boldsymbol\Gamma_e$, this term equals $+\tfrac52k_BT_e\,\boldsymbol\Gamma_e$.
* $-\lambda'\nabla T_e$ is **Fourier heat conduction**: heat flows from hot to cold.
* $-T_e\varphi\boldsymbol{\mathcal E}$ is **heat driven by the field** (the Peltier effect).

### 2.3 Why the same $\varphi$ appears twice (Onsager reciprocity)

The coefficient coupling "temperature gradient → current" ($\varphi$) and the one coupling "field → heat" ($T_e\varphi$) are not independent. Onsager's reciprocity theorem, a consequence of the time-reversibility of collisions, says these cross-coefficients are equal once the fluxes and forces are written in matching form. That is why the equations need only **three** independent coefficients ($\sigma,\varphi,\lambda'$), not four. Section 3.2 confirms this directly: computing both cross-terms from the kinetic theory gives the same integral.

### 2.4 Why $\tfrac52k_BT_e$ and not $\tfrac32k_BT_e$

Each electron has mean kinetic energy $\tfrac32k_BT_e$. But when a parcel of gas flows into a region, it must also push the gas ahead of it out of the way. That pushing is work, $p_e$ per unit volume, or $k_BT_e$ per electron. The total energy carried per particle is the **enthalpy** per particle, $\tfrac32k_BT_e+k_BT_e=\tfrac52k_BT_e$. You can also see it directly in the kinetic calculation (Section 3.5): the factor $(x-\tfrac52)$ appears there.

### 2.5 Eliminating the field: eq. (28) second form and eq. (29)

The field is often not what you control or know; the current is. So we rewrite $\mathbf q_e$ in terms of $\mathbf J_e$ and $\nabla T_e$. Solve the current equation for the field:

$$
\boldsymbol{\mathcal E}=\frac{\mathbf J_e-\varphi\nabla T_e}{\sigma}.
$$

Substitute into $\mathbf q_e$:

$$
\mathbf q_e=-\frac52\frac{k_BT_e}{e}\mathbf J_e-\lambda'\nabla T_e-\frac{T_e\varphi}{\sigma}\mathbf J_e+\frac{T_e\varphi^2}{\sigma}\nabla T_e .
$$

Collect the $\nabla T_e$ terms and the $\mathbf J_e$ terms:

$$
\mathbf q_e=-\underbrace{\left(\lambda'-T_e\frac{\varphi^2}{\sigma}\right)}_{\lambda}\nabla T_e+\underbrace{\left(-\frac52\frac{k_BT_e}{e}-T_e\frac{\varphi}{\sigma}\right)}_{\psi}\mathbf J_e .
$$

This is eq. (28) second form, with $\lambda$ and $\psi$ as in eq. (29).

### 2.6 Two thermal conductivities: $\lambda'$ versus $\lambda$

* $\lambda'$ is the heat conductivity measured while holding the field $\boldsymbol{\mathcal E}=0$.
* $\lambda$ is the heat conductivity measured while holding the current $\mathbf J_e=0$.

They differ for the following reason. Suppose you heat one end of an isolated piece of plasma. The temperature gradient tries to drive a current $\varphi\nabla T_e$. With nowhere for the current to go, charge piles up and builds a field $\boldsymbol{\mathcal E}=-(\varphi/\sigma)\nabla T_e$ that exactly stops the current. That field also drives a heat flux, $-T_e\varphi\boldsymbol{\mathcal E}=+T_e(\varphi^2/\sigma)\nabla T_e$, which flows *against* the conduction. So the net conduction at zero current is reduced:

$$
\lambda=\lambda'-T_e\varphi^2/\sigma\ \le\ \lambda'.
$$

**Which one does the model need?** The electron energy equation (13)–(14) already contains, as separate terms:

* the energy carried by the electron flux ($\tfrac52k_BT_e\boldsymbol\Gamma_e$ in $\mathbf h_e$),
* the work done by the field (the Joule term $-e\boldsymbol\Gamma_e\cdot\mathbf E$).

Everything to do with current and field is therefore already accounted for. What is left over is the conduction that happens *with the current treated separately*, which is $\lambda$, not $\lambda'$. Using $\lambda'$ would count the field-driven heat twice. Hence the appendix's statement "$\kappa_e=\lambda$, not $\lambda'$".

### 2.7 Matching the Onsager form to the model's variables: eq. (30)

The model does not use $\sigma,\varphi,\lambda$ directly. It uses a friction frequency $\bar\nu^M_e$ and a thermal-diffusion coefficient $C_e$ in the electron momentum balance (3). Take that balance with the heavy particles at rest (they move much more slowly than electrons) and drop the tiny inertia term $Y_e\mathbf F$:

$$
0=-en_e\mathbf E-\nabla p_e-m_en_e\bar\nu^M_e\,\mathbf u_e-C_en_ek_B\nabla T_e .
$$

In words: **electric force + pressure force = friction + thermal-diffusion force**. Solve for $n_e\mathbf u_e=\boldsymbol\Gamma_e$ and multiply by $-e$:

$$
\mathbf J_e=-en_e\mathbf u_e=\frac{n_ee^2}{m_e\bar\nu^M_e}\left(\mathbf E+\frac{\nabla p_e}{en_e}\right)+\frac{en_ek_BC_e}{m_e\bar\nu^M_e}\nabla T_e .
$$

Compare with $\mathbf J_e=\sigma\boldsymbol{\mathcal E}+\varphi\nabla T_e$, term by term:

$$
\sigma=\frac{n_ee^2}{m_e\bar\nu^M_e}\ \Rightarrow\ \bar\nu^M_e=\frac{n_ee^2}{m_e\sigma},\qquad
\varphi=\sigma\frac{k_B}{e}C_e\ \Rightarrow\ C_e=\frac{e}{k_B}\frac{\varphi}{\sigma}.
$$

So the **friction frequency is just the conductivity in disguise**, and $C_e$ is the thermal-diffusion coefficient made dimensionless.

Next, the energy flux (14) with the heavy particles at rest, and $C_i+C_n=C_e$:

$$
\mathbf h_e=\tfrac52k_BT_e\boldsymbol\Gamma_e-\kappa_e\nabla T_e+k_BT_eC_e\boldsymbol\Gamma_e
=-\left(\tfrac52+C_e\right)\frac{k_BT_e}{e}\mathbf J_e-\kappa_e\nabla T_e .
$$

Compare with $\mathbf q_e=-\lambda\nabla T_e+\psi\mathbf J_e$:

$$
\kappa_e=\lambda,\qquad \tfrac52+C_e=-\frac{e\psi}{k_BT_e}.
$$

The last relation is a **consistency check**, not new information. Put $\psi$ from eq. (29) into it: $-e\psi/(k_BT_e)=\tfrac52+(e/k_B)(\varphi/\sigma)=\tfrac52+C_e$. It holds automatically. It confirms that the same $C_e$ in the force balance (thermal-diffusion force) and in the heat flux (thermal-diffusion heat flux) is the Onsager pair of Section 2.3.

**Result of 3.1.** If we can compute the three numbers $\sigma$, $\varphi$, and $\lambda$, we know all of the model's coefficients:

$$
\bar\nu^M_e=\frac{n_ee^2}{m_e\sigma},\qquad C_e=\frac{e}{k_B}\frac{\varphi}{\sigma},\qquad \kappa_e=\lambda .
$$

Sections 3.2 and 3.3 compute them in two extreme cases (mostly neutrals; fully ionized), and 3.4 bridges between the two.

---

## 3. Section 3.2: when electrons hit mostly neutral atoms

### 3.1 The "Lorentz gas" picture

In a weakly ionized gas, an electron almost always collides with a neutral atom. Lorentz's simplification is:

* the atoms are heavy, immobile scatterers (by 1.4, they barely recoil);
* electron–electron collisions are ignored;
* each electron interacts only with the atoms and the fields.

Each electron then moves independently through a forest of fixed obstacles, and the kinetic equation becomes linear and solvable.

*Caveat on applicability:* at LTE 10 000 K argon is only about 2 % ionized, but Coulomb (electron–ion) collisions are so much stronger than electron–atom ones that they already dominate ($\bar\nu^E_{ei}/\bar\nu^E_{en}\approx 12$, see `ask_claude.md`). The pure Lorentz result is the cold, weakly ionized end of the problem. Section 3.4 of the appendix blends it with the Coulomb result of Section 3.3.

### 3.2 The two-term expansion

By 1.4, collisions scramble directions quickly, so $f$ is nearly isotropic. We write it as a large isotropic part plus a small part that depends on direction:

$$
f(\mathbf r,\mathbf C)=f_M(C)+\frac{\mathbf C}{C}\cdot\mathbf f^1(C),\qquad |\mathbf f^1|\ll f_M .
$$

$\mathbf C/C$ is the unit vector along the velocity. So the second term is positive for electrons moving one way and negative for those moving the other way. It is the "tilt" of the distribution that makes the current. Only $\mathbf f^1$ carries flux. The Maxwellian part integrates to zero flux, as noted in 1.2.

### 3.3 Solving for the tilt $\mathbf f^1$

The **Boltzmann equation** in steady state says: drift in space plus acceleration by forces equals the change due to collisions,

$$
\mathbf C\cdot\nabla f-\frac{e\mathbf E}{m_e}\cdot\frac{\partial f}{\partial\mathbf C}=\left(\frac{\partial f}{\partial t}\right)_{\rm coll}.
$$

For the tilted part, collisions with fixed scatterers simply destroy the tilt at the rate $\nu_{eH}(C)$:

$$
\left(\partial f/\partial t\right)_{\rm coll}\big|_{\rm aniso}=-\nu_{eH}\,\frac{\mathbf C}{C}\cdot\mathbf f^1 .
$$

Insert $f\approx f_M$ on the left, which is the dominant part, and use $\partial f_M/\partial\mathbf C=(\mathbf C/C)\,\partial f_M/\partial C$. Matching the terms proportional to $\mathbf C/C$ gives

$$
C\,\nabla f_M-\frac{e\mathbf E}{m_e}\frac{\partial f_M}{\partial C}=-\nu_{eH}\,\mathbf f^1
\quad\Longrightarrow\quad
\mathbf f^1=\frac{e}{m_e\nu_{eH}}\frac{\partial f_M}{\partial C}\,\mathbf E-\frac{C}{\nu_{eH}}\nabla f_M .
$$

This is a balance: the field and the gradients tilt the distribution, and collisions un-tilt it. A strong tilt needs strong forcing or rare collisions (small $\nu_{eH}$).

Now evaluate $\nabla f_M$. Since $\ln f_M=\ln n_e-\tfrac32\ln T_e-x+\text{const}$ with $x=m_eC^2/2k_BT_e$, and $\nabla x=-x\nabla\ln T_e$,

$$
\nabla\ln f_M=\nabla\ln n_e+\left(x-\tfrac32\right)\nabla\ln T_e
=\nabla\ln p_e+\left(x-\tfrac52\right)\nabla\ln T_e ,
$$

using $\ln p_e=\ln n_e+\ln T_e+\text{const}$. Also $\partial f_M/\partial C=-(m_eC/k_BT_e)f_M$. Then

$$
\mathbf f^1=\underbrace{\frac{e}{m_e\nu_{eH}}\frac{\partial f_M}{\partial C}\,\mathbf E-\frac{C}{\nu_{eH}}f_M\nabla\ln p_e}_{\text{both}\ \propto\ C f_M/\nu_{eH}}-\frac{C}{\nu_{eH}}f_M\left(x-\tfrac52\right)\nabla\ln T_e .
$$

The first two terms have the same shape. Writing the pressure term as $\frac{e}{m_e\nu_{eH}}\frac{\partial f_M}{\partial C}\cdot\frac{\nabla p_e}{en_e}$ combines them into the generalized field of 2.1:

$$
\boxed{\mathbf f^1=\frac{e}{m_e\nu_{eH}}\frac{\partial f_M}{\partial C}\,\boldsymbol{\mathcal E}-\frac{C}{\nu_{eH}}f_M\left(x-\tfrac52\right)\nabla\ln T_e}
$$

This is the microscopic origin of the generalized field: the pressure gradient tilts the distribution in exactly the same way an electric field does.

> **Note:** the paper's eq. (31) writes $\boldsymbol{\mathcal E}$ **and** keeps $\nabla\ln p_e$ in the bracket, which counts the pressure force twice. Either $\boldsymbol{\mathcal E}\to\mathbf E$ or the $\nabla\ln p_e$ should be dropped (Section 6). The coefficients in eq. (36), and the PlasmeX tables (`exe/gen_lorentz_tables.py`), use the correct form above.

The factor $(x-\tfrac52)$ deserves a comment. A temperature gradient at fixed pressure does not push "all electrons" one way. It pushes **fast** electrons ($x>\tfrac52$) down the gradient and **slow** ones ($x<\tfrac52$) up it. If all electrons had the same collision frequency, the two would cancel exactly (shown in 3.7). Any net current from $\nabla T_e$ is therefore a measure of how unequal fast-electron and slow-electron collision rates are.

### 3.4 From the tilt to the fluxes

The current is $\mathbf J_e=-e\int\mathbf C f\,d^3C$. Only the tilted part contributes, and the angular integral of $\mathbf C(\mathbf C\cdot\mathbf f^1)/C$ gives a factor $\tfrac13$:

$$
\mathbf J_e=-\frac{4\pi e}{3}\int_0^\infty C^3\,\mathbf f^1\,dC,\qquad
\mathbf q_e=\frac{4\pi}{3}\int_0^\infty \tfrac12m_eC^2\,C^3\,\mathbf f^1\,dC .
$$

Every $\mathbf f^1$ carries $1/\nu_{eH}$, and the moment adds $C^3$. So every flux is an integral of the form

$$
\int_0^\infty(\text{something})\;\frac{C^4f_M}{\nu_{eH}}\,dC .
$$

This shape gives rise to the two kinds of average defined next.

### 3.5 The two kinds of average: eqs. (32) and (33)

**(a) The plain Maxwellian average**, eq. (32): the ordinary mean of a quantity $g(C)$ over all electrons (here $f_M$ is normalized to 1):

$$
\langle g\rangle_0=\int g\,f_M\,d^3C=4\pi\int_0^\infty g\,C^2f_M\,dC .
$$

The $4\pi C^2$ is the area of a sphere of radius $C$ in velocity space: there are more ways to have speed $C$ when $C$ is large.

**(b) The transport average**, eq. (33): the mean of $g$ weighted by how much each electron contributes to transport,

$$
\langle g\rangle_\nu=\frac{\langle g\,C^2/\nu_{eH}\rangle_0}{\langle C^2/\nu_{eH}\rangle_0}=\frac{\int g\,C^4f_M\nu_{eH}^{-1}\,dC}{\int C^4f_M\nu_{eH}^{-1}\,dC}.
$$

Why the weight $C^2/\nu_{eH}$? It is (speed) × (speed / collision rate) = speed × distance per collision, which is the electron's **diffusivity**. An electron that moves fast and rarely collides carries much more current than one that is constantly being knocked around. In energy variables the weight is $\varepsilon\,e^{-\varepsilon/k_BT_e}/Q^{(1)}(\varepsilon)$. For argon this strongly favours electrons near the Ramsauer minimum, where $Q^{(1)}$ is tiny.

The appendix stresses that the transport average is used **only** for the two dimensionless numbers $\langle x\rangle_\nu$ and $\langle x^2\rangle_\nu$, never to define a collision frequency.

### 3.6 The two collision frequencies: eqs. (34) and (35)

Both are averages of $\nu_{eH}$ over the same weight $C^2f_M$ (electrons weighted by their kinetic energy), but taken in two different ways.

#### (i) Momentum-transfer (mobility) average, $\bar\nu^M_{en}$: a harmonic mean

Take only the field term of $\mathbf f^1$ and compute the current. Using $-\partial f_M/\partial C=(m_eC/k_BT_e)f_M$:

$$
\mathbf J_e=\frac{4\pi e^2}{3m_e}\int_0^\infty\frac{C^3}{\nu_{eH}}\left(-\frac{\partial f_M}{\partial C}\right)dC\;\boldsymbol{\mathcal E}
=\frac{n_ee^2}{3k_BT_e}\left\langle\frac{C^2}{\nu_{eH}}\right\rangle_0\boldsymbol{\mathcal E}\equiv\sigma_{en}\boldsymbol{\mathcal E}.
$$

Now **define** $\bar\nu^M_{en}$ as the single constant frequency that gives this same conductivity in the fluid formula $\sigma=n_ee^2/(m_e\bar\nu)$. With $\langle C^2\rangle_0=3k_BT_e/m_e$:

$$
\bar\nu^M_{en}=\frac{\langle C^2\rangle_0}{\langle C^2/\nu_{eH}\rangle_0}.
$$

This is eq. (34), first form. It averages $1/\nu$ and then inverts, which is a **harmonic mean**, the same rule used for resistors in parallel. Each electron speed is a parallel "channel" for current, so electrons with *small* $\nu$ (the easy channels) dominate.

*Second form.* Put $\nu_{eH}=n_nCQ^{(1)}$, so $C^2/\nu_{eH}=C/(n_nQ^{(1)})$, and $\langle\cdot\rangle_0=4\pi\int C^2(\cdot)f_M\,dC$:

$$
\bar\nu^M_{en}=n_n\frac{3k_BT_e/m_e}{4\pi\int_0^\infty C^3f_M\,[Q^{(1)}]^{-1}dC}.
$$

*Third form.* Change variable to energy, $\varepsilon=\tfrac12m_eC^2$, so $C^3dC=(2\varepsilon/m_e^2)\,d\varepsilon$. Insert $f_M=(m_e/2\pi k_BT_e)^{3/2}e^{-\varepsilon/k_BT_e}$. The constants collect into

$$
\bar\nu^M_{en}=n_n\sqrt{\frac{2}{m_e}}\;\frac{3\sqrt\pi}{4}\;\frac{(k_BT_e)^{5/2}}{\int_0^\infty\varepsilon\,e^{-\varepsilon/k_BT_e}\,[Q^{(1)}(\varepsilon)]^{-1}d\varepsilon}.
$$

This is what is evaluated numerically: one integral over the tabulated cross section.

*Check with a constant cross section.* If $Q^{(1)}=Q$ is constant, the integral is $(k_BT_e)^2/Q$, so $\bar\nu^M_{en}=n_nQ\sqrt{2k_BT_e/m_e}\cdot\tfrac{3\sqrt\pi}{4}$. Divide by the "naive" value $\tfrac43n_n\bar c_eQ=\tfrac43n_nQ\sqrt{8k_BT_e/\pi m_e}$:

$$
\frac{\tfrac{3\sqrt\pi}{4}\sqrt2}{\tfrac43\sqrt{8/\pi}}=\frac{9\pi}{32}\approx0.88 ,
$$

which is the number quoted in the appendix. Even with a constant cross section the two averages are not identical, because $\nu=n_nCQ$ still varies with speed.

#### (ii) Energy-exchange average, $\bar\nu^E_{en}$: an arithmetic mean

In each collision an electron of energy $\varepsilon$ loses on average $(2m_e/M)\,\varepsilon\,(1-T_h/T_e)$ to the atom (1.4; the $T_h/T_e$ part accounts for the atoms' own thermal motion, which gives some energy back). The total loss rate per electron adds this up over all electrons: $\tfrac{2m_e}{M}\langle\varepsilon\,\nu_{eH}\rangle_0(1-T_h/T_e)$. Writing $\varepsilon=\tfrac12m_eC^2$ and $\langle\varepsilon\rangle_0=\tfrac32k_BT_e$,

$$
\frac{2m_e}{M}\langle\varepsilon\,\nu_{eH}\rangle_0\left(1-\frac{T_h}{T_e}\right)=\frac32k_B(T_e-T_h)\,\frac{2m_e}{M}\,\bar\nu^E_{en},
\qquad
\bar\nu^E_{en}=\frac{\langle C^2\nu_{eH}\rangle_0}{\langle C^2\rangle_0}.
$$

This is eq. (35), and it is exactly the elastic term of the energy equation (13). Here $\nu$ itself is averaged: an **arithmetic mean**. Electrons that collide *often* dominate, which is the opposite of the mobility average.

In energy variables this becomes $\bar\nu^E_{en}=\tfrac43n_n\bar c_e\,\overline Q^{(1,1)}_{en}$ with $\overline Q^{(1,1)}_{en}=\tfrac12\int_0^\infty x^2Q^{(1)}(xk_BT_e)\,e^{-x}dx$. That is the standard "collision integral" form. For a constant $Q$ it gives $\overline Q^{(1,1)}=Q$.

#### (iii) Why they differ, and why it matters

For any positive function, the **harmonic mean ≤ arithmetic mean** over the same weights, with equality only if the function is constant. So always $\bar\nu^M_{en}\le\bar\nu^E_{en}$.

For argon the gap is large because of the Ramsauer minimum:

* the **mobility** average is dominated by slow electrons sitting in the minimum, where they barely collide, so $\bar\nu^M$ is small;
* the **energy** average is dominated by fast electrons above it, where they collide a lot, so $\bar\nu^E$ is large.

The model therefore uses $\bar\nu^M$ in the force balances (3)–(5) and $\bar\nu^E$ in the energy equation (13). Swapping them would make the friction about 4× too large, or the elastic energy loss about 4× too small (numbers in Section 4).

### 3.7 The three Lorentz coefficients: eq. (36)

With the averages defined, the three coefficients needed by 3.1 follow from the integrals of 3.4.

**Conductivity.** From 3.6(i):

$$
\sigma_{en}=\frac{n_ee^2}{m_e\bar\nu^M_{en}} .
$$

**Thermal diffusion.** Take the $\nabla T_e$ term of $\mathbf f^1$ and compute its current:

$$
\varphi_{en}\nabla T_e=\frac{en_e}{3T_e}\left\langle\frac{C^2}{\nu_{eH}}\left(x-\tfrac52\right)\right\rangle_0\nabla T_e
\quad\Rightarrow\quad
\left(\frac{\varphi}{\sigma}\right)_{en}=\frac{k_B}{e}\left(\langle x\rangle_\nu-\tfrac52\right).
$$

So $C_e=\langle x\rangle_\nu-\tfrac52$: **the average energy of the current-carrying electrons, minus 5/2**.

*Check:* if $\nu_{eH}$ is constant, the transport weight in energy is $x^{3/2}e^{-x}$, whose mean is $\Gamma(7/2)/\Gamma(5/2)=5/2$. So $C_e=0$, as stated in the appendix. A nonzero $C_e$ exists **only because** $Q^{(1)}$ depends on energy.

**Zero-current thermal conductivity.** Compute the heat flux from the $\nabla T_e$ term at $\boldsymbol{\mathcal E}=0$, which gives $\lambda'$. Then subtract $T_e\varphi^2/\sigma$ (2.6). The algebra collapses neatly:

$$
\lambda'_{en}=\frac{n_ek_B^2T_e}{m_e\bar\nu^M_{en}}\left\langle\left(x-\tfrac52\right)^2\right\rangle_\nu,
\qquad
\lambda_{en}=\lambda'_{en}-T_e\frac{\varphi^2}{\sigma}=\frac{n_ek_B^2T_e}{m_e\bar\nu^M_{en}}\left[\langle x^2\rangle_\nu-\langle x\rangle_\nu^2\right].
$$

The zero-current conductivity is proportional to the **variance** of the electron energy under the transport weight. The interpretation: at zero current, as many electrons go each way, so heat can be conducted only by fast electrons going one way while slow ones go the other. The more spread in energy among the current carriers, the more heat that exchange carries. For constant $\nu$, the variance of $x^{3/2}e^{-x}$ is $5/2$, giving the classical Lorentz value $\lambda=\tfrac52\,n_ek_B^2T_e/(m_e\nu)$.

### 3.8 How it is evaluated in practice

All of eqs. (33)–(36) are one-dimensional integrals over the measured argon cross section (LXCat). They depend only on $T_e$, and every collision frequency is proportional to $n_n$. So they are computed once, **per unit neutral density**, and tabulated against $T_e$ (`~/work/PlasmeX/exe/gen_lorentz_tables.py` → `ar_e_lorentz.dat`). During a simulation the code just looks up the table and multiplies by the local $n_n$.

---

## 4. Numbers for argon (from `ar_Q_en.dat`)

| $T_e$ [eV] | $\bar\nu^M_{en}/n_n$ [m³/s] | $\bar\nu^E_{en}/n_n$ [m³/s] | $\bar\nu^E/\bar\nu^M$ | $\langle x\rangle_\nu$ | $C_e=\langle x\rangle_\nu-\tfrac52$ | variance $\langle x^2\rangle_\nu-\langle x\rangle_\nu^2$ |
|---|---|---|---|---|---|---|
| 0.5 | $3.0\times10^{-15}$ | $1.28\times10^{-14}$ | 4.3 | 1.00 | −1.50 | 0.72 |
| 1.0 | $8.8\times10^{-15}$ | $3.85\times10^{-14}$ | 4.4 | 0.89 | −1.61 | 0.82 |
| 1.5 | $1.7\times10^{-14}$ | $7.26\times10^{-14}$ | 4.3 | 0.88 | −1.62 | 0.87 |
| 2.0 | $2.7\times10^{-14}$ | $1.09\times10^{-13}$ | 4.0 | 0.88 | −1.62 | 0.93 |
| *constant $\nu$* | | | *1* | *2.50* | *0* | *2.50* |

Reading the table:

* **$\bar\nu^E/\bar\nu^M\approx4$–4.4.** This confirms the appendix's "factor of 4–4.5". The two averages truly are different numbers.
* **$\langle x\rangle_\nu\approx0.9$, well below 2.5.** The current is carried by electrons with only about a third of the typical thermal energy: the slow ones trapped near the Ramsauer minimum.
* **$C_e\approx-1.6$, negative.** In the electron flux (9), the temperature-gradient term is $-(1+C_e)\,\frac{n_ek_B}{m_e\bar\nu^M_e}\nabla T_e$, and $1+C_e\approx-0.6<0$. So in a purely electron–neutral gas, a temperature gradient at fixed density pushes electrons **toward the hotter region**. The physical reason: on the hot side, electrons are faster but collide much more (rising $Q^{(1)}$), so they diffuse *less* than the cooler electrons. The net flow is from cold to hot.
* **Variance ≈ 0.8 instead of 2.5.** Electron heat conduction in argon is about 3× weaker than the textbook constant-$\nu$ formula would give with the same $\bar\nu^M$.

---

## 5. Summary: what goes where in the model

```
measured Q^(1)(ε) for argon (LXCat)
        │
        ▼  Lorentz-gas kinetic theory (3.2)
   ┌───────────────────────────────┬─────────────────────────────┐
   │ harmonic mean  ν̄^M_en (34)    │ arithmetic mean ν̄^E_en (35) │
   │ ⟨x⟩_ν, ⟨x²⟩_ν (33)            │                             │
   └──────────────┬────────────────┴──────────────┬──────────────┘
                  ▼ (36)                          │
        σ_en, φ_en, λ_en                          │
                  │ bridged with Coulomb (3.3–3.4)│
                  ▼ mapping (30), Section 3.1     ▼
   ν̄^M_e, C_e  → electron flux (9)      ν̄^E_eh → elastic loss in energy eq. (13)
   κ_e = λ, C_e → electron energy flux (14)
```

The ideas in one paragraph: electrons are pushed by an effective field (electric plus pressure) and by a temperature gradient. Their response (current and heat) needs three numbers, and these map one-to-one onto the friction frequency, thermal-diffusion coefficient, and thermal conductivity of the fluid model. When electrons mostly hit neutral atoms, the three numbers come from a short kinetic calculation in which collisions un-tilt a nearly isotropic distribution. Because argon's cross section varies 150-fold with energy, "the" collision frequency splits into a small mobility value and a large energy-exchange value, about 4× apart. Each must be used in its own equation.

---

## 6. Notes for the authors (possible corrections; the paper was not changed)

1. **Eq. (31) double-counts the pressure gradient.** It has $\boldsymbol{\mathcal E}=\mathbf E+\nabla p_e/(en_e)$ in the first term **and** $\nabla\ln p_e$ in the bracket. Since $\frac{e}{m_e\nu_{eH}}\frac{\partial f_M}{\partial C}\frac{\nabla p_e}{en_e}=-\frac{C}{\nu_{eH}}f_M\nabla\ln p_e$, the pressure force appears twice. The correct form is either

   $$\mathbf f^1=\frac{e}{m_e\nu_{eH}}\frac{\partial f_M}{\partial C}\boldsymbol{\mathcal E}-\frac{C}{\nu_{eH}}f_M\left(x-\tfrac52\right)\nabla\ln T_e$$

   or the same with $\boldsymbol{\mathcal E}\to\mathbf E$ and $\nabla\ln p_e$ kept. Eq. (36) and the code (`gen_lorentz_tables.py`, whose $\varphi$ integral has only $(x-\tfrac52)$) already correspond to the correct form. Only the printed equation is affected.
2. **"Heat flux" $\mathbf q_e$ in 3.1** is the total kinetic-energy flux, including the $\tfrac52k_BT_e\boldsymbol\Gamma_e$ enthalpy part. That is the same object as $\mathbf h_e$ in eq. (14). Saying so explicitly ("energy flux $\mathbf q_e\equiv\mathbf h_e$") would make the comparison in eq. (30) easier to follow.
3. **"eH" in $\nu_{eH}$** is never spelled out. A short "(e = electron, H = heavy particle, here the neutral atom)" would help.
4. **3.1 applies the mapping (30) with the heavy species at rest.** That is fine (the drag terms are added separately in eq. (9)), but it is not stated. One clause would remove a likely reader question.

---

## 7. Symbol table

| symbol | meaning | units |
|---|---|---|
| $f,\ f_M,\ \mathbf f^1$ | distribution function; its Maxwellian part; its first-order (tilt) part | s³ m⁻⁶ |
| $\mathbf C,\ C$ | electron velocity, speed | m/s |
| $x=\varepsilon/k_BT_e$ | electron energy relative to thermal energy | – |
| $Q^{(1)}(\varepsilon)$ | electron–atom momentum-transfer cross section | m² |
| $\nu_{eH}(C)=n_nCQ^{(1)}$ | speed-dependent electron–atom collision frequency | s⁻¹ |
| $\bar\nu^M_{en}$ | mobility (harmonic) average of $\nu_{eH}$, used in force balances | s⁻¹ |
| $\bar\nu^E_{en}$ | energy-exchange (arithmetic) average of $\nu_{eH}$, used in the energy equation | s⁻¹ |
| $\overline Q^{(1,1)}_{en}$ | collision integral, $\bar\nu^E_{en}=\tfrac43n_n\bar c_e\overline Q^{(1,1)}_{en}$ | m² |
| $\langle\cdot\rangle_0,\ \langle\cdot\rangle_\nu$ | Maxwellian average; transport (diffusivity-weighted) average | – |
| $\boldsymbol{\mathcal E}=\mathbf E+\nabla p_e/(en_e)$ | generalized field | V/m |
| $\sigma$ | electrical conductivity | S/m |
| $\varphi$ | thermal-diffusion (thermoelectric) coefficient | A m⁻¹ K⁻¹ |
| $\lambda',\ \lambda$ | thermal conductivity at zero field; at zero current ($=\kappa_e$) | W m⁻¹ K⁻¹ |
| $\psi$ | heat carried per unit current | V |
| $C_e=(e/k_B)\varphi/\sigma$ | dimensionless electron thermal-diffusion coefficient | – |
| $\bar c_e=\sqrt{8k_BT_e/\pi m_e}$ | mean electron speed | m/s |
