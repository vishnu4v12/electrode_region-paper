# Governing Equations and Transport Coefficients of the Drift–Diffusion Plasma Model

*Reference notes for the PlasmeX solver — three-species (electron, Ar$^+$, Ar) drift–diffusion
plasma model with an electrostatic field solve and a semi-implicit, operator-split species
treatment.*

## Abstract

This note collects, in closed form, the governing equations solved by the code and the
transport-coefficient formulas that close them. The species momentum balances are inertialess
(Stefan–Maxwell) force balances that are solved algebraically for the drift velocities; the
electron transport coefficients are obtained from the electron Boltzmann equation in two
limits — a Lorentz gas dominated by electron–neutral collisions, and a fully ionized plasma
dominated by electron–ion (Coulomb) collisions — bridged by a resistance-type interpolation
formula. Ion transport uses a field-dependent (Wannier) charge-exchange mobility with an
electron-friction correction. The neutral density is not an independent unknown; it is closed
algebraically or by its own continuity/momentum equation, selectable at run time.

---

## 1. Notation and conventions

Three species are carried: electrons, singly charged argon ions Ar$^+$, and neutral argon
atoms Ar. Subscripts $\mathrm e$, $\mathrm i$, $\mathrm n$ denote electron, ion, neutral;
$\mathrm h$ denotes the heavy-particle (ion + neutral) mixture.

| symbol | meaning | unit |
|---|---|---|
| $n_{\mathrm e},\,n_{\mathrm i},\,n_{\mathrm n}$ | electron, ion, neutral number density | m$^{-3}$ |
| $n_{\mathrm h}=n_{\mathrm i}+n_{\mathrm n}$ | heavy-particle number density | m$^{-3}$ |
| $\rho = M n_{\mathrm h}$ | heavy-particle mass density | kg m$^{-3}$ |
| $T_{\mathrm e}$ | electron temperature | K |
| $T_{\mathrm g}=T_{\mathrm i}=T_{\mathrm n}$ | heavy-particle (gas) temperature, shared by ions and neutrals | K |
| $\mathcal E = \tfrac32 n_{\mathrm e} k_{\mathrm B} T_{\mathrm e}$ | electron thermal energy density | J m$^{-3}$ |
| $\phi$, $\boldsymbol E = -\nabla\phi$ | electrostatic potential, field | V, V m$^{-1}$ |
| $\boldsymbol u_{\mathrm e},\,\boldsymbol u_{\mathrm i},\,\boldsymbol u_{\mathrm n}$ | species mean (fluid) velocities | m s$^{-1}$ |
| $\boldsymbol u \equiv \boldsymbol u_{\mathrm h}$ | mass-averaged heavy-particle velocity | m s$^{-1}$ |
| $\boldsymbol\Gamma_{\mathrm e}=n_{\mathrm e}\boldsymbol u_{\mathrm e}$, etc. | number flux | m$^{-2}$ s$^{-1}$ |
| $p_{\mathrm e}=n_{\mathrm e}k_{\mathrm B}T_{\mathrm e}$, $p_{\mathrm i}=n_{\mathrm i}k_{\mathrm B}T_{\mathrm g}$, $p_{\mathrm n}=n_{\mathrm n}k_{\mathrm B}T_{\mathrm g}$ | partial pressures | Pa |
| $p = p_{\mathrm e}+p_{\mathrm i}+p_{\mathrm n}$ | total pressure | Pa |
| $e$, $m_{\mathrm e}$, $M$, $k_{\mathrm B}$, $\varepsilon_0$ | elementary charge, electron mass, heavy-particle (argon) mass, Boltzmann constant, vacuum permittivity | SI |

The stored electron mobility $\mu_{\mathrm e}$ is *negative* (it carries the electron's charge
sign); the ion mobility $\mu_{\mathrm i}$ is positive. All transport coefficients below are
functions of the local state $(n_{\mathrm e},n_{\mathrm i},n_{\mathrm n},T_{\mathrm e},T_{\mathrm g},|\boldsymbol E|)$
evaluated pointwise (in the numerical scheme, on cell faces from the two adjacent cell values).

Two geometries are supported: 1-D planar ($x$) and 2-D axisymmetric $r$–$z$. In the latter every
divergence and gradient carries the cylindrical metric,
$\nabla\cdot\boldsymbol F=\tfrac1r\partial_r(rF_r)+\partial_zF_z$.

---

## 2. Governing equations

### 2.1 Charged-species continuity

$$
\partial_t n_{\mathrm e}+\nabla\cdot\boldsymbol\Gamma_{\mathrm e} = S,
\qquad
\partial_t n_{\mathrm i}+\nabla\cdot\boldsymbol\Gamma_{\mathrm i} = S,
\tag{1}
$$

$$
S = k_{\mathrm i}(T_{\mathrm e})\,n_{\mathrm e}n_{\mathrm n} - k_{\mathrm r}(T_{\mathrm e})\,n_{\mathrm e}^2 n_{\mathrm i},
\tag{2}
$$

with $k_{\mathrm i}$ and $k_{\mathrm r}$ the (tabulated) ionization and three-body recombination
rate coefficients, §3.10. Electron-impact ionization and its inverse, three-body recombination,
are the only reactions; there is no radiative loss channel.

### 2.2 Momentum balances and the Stefan–Maxwell flux closure

The three species momentum equations are taken in the inertialess (drift–diffusion) limit — the
momentum flux and time-derivative terms are dropped, leaving a local force balance between the
electric force, the pressure gradient, interspecies friction, and (for the electron) a
thermal-diffusion force:

$$
0 = -e n_{\mathrm e}\boldsymbol E - \nabla p_{\mathrm e}
    - m_{\mathrm e} n_{\mathrm e}\nu_{\mathrm{en}}(\boldsymbol u_{\mathrm e}-\boldsymbol u_{\mathrm n})
    - m_{\mathrm e} n_{\mathrm e}\nu_{\mathrm{ei}}^{\mathrm{eff}}(\boldsymbol u_{\mathrm e}-\boldsymbol u_{\mathrm i})
    - C_{\mathrm e} n_{\mathrm e} k_{\mathrm B}\nabla T_{\mathrm e},
\tag{3}
$$

$$
0 = +e n_{\mathrm i}\boldsymbol E - \nabla p_{\mathrm i}
    - M n_{\mathrm i}\nu_{\mathrm{in}}(\boldsymbol u_{\mathrm i}-\boldsymbol u_{\mathrm n})
    - M n_{\mathrm i}\nu_{\mathrm{ie}}(\boldsymbol u_{\mathrm i}-\boldsymbol u_{\mathrm e})
    + A_{\mathrm i} n_{\mathrm e} k_{\mathrm B}\nabla T_{\mathrm e},
\tag{4}
$$

$$
0 = -\nabla p_{\mathrm n}
    - M n_{\mathrm n}\nu_{\mathrm{ni}}(\boldsymbol u_{\mathrm n}-\boldsymbol u_{\mathrm i})
    - M n_{\mathrm n}\nu_{\mathrm{ne}}(\boldsymbol u_{\mathrm n}-\boldsymbol u_{\mathrm e})
    + A_{\mathrm a} n_{\mathrm e} k_{\mathrm B}\nabla T_{\mathrm e}.
\tag{5}
$$

Here $\nu_{\mathrm{en}}$, $\nu_{\mathrm{ei}}^{\mathrm{eff}}$, $\nu_{\mathrm{in}}$, $\nu_{\mathrm{ie}}$,
$\nu_{\mathrm{ni}}$, $\nu_{\mathrm{ne}}$ are momentum-transfer collision frequencies (§3),
$C_{\mathrm e}$, $A_{\mathrm i}$, $A_{\mathrm a}$ are the electron thermal-diffusion coefficient
and its momentum-conserving reactions on the ion and neutral, and
$\nu_{\mathrm{ei}}^{\mathrm{eff}}=0.5\,\bar\nu_{\mathrm{ei}}$ is the conductivity-weighted
electron–ion friction frequency (distinct from the full equilibration rate $\bar\nu_{\mathrm{ei}}$
used in the energy exchange, §3.8). The frequencies are tied pairwise by momentum conservation
(Newton's third law on each collision pair):

$$
m_{\mathrm e} n_{\mathrm e}\nu_{\mathrm{en}} = M n_{\mathrm n}\nu_{\mathrm{ne}},\qquad
m_{\mathrm e} n_{\mathrm e}\nu_{\mathrm{ei}}^{\mathrm{eff}} = M n_{\mathrm i}\nu_{\mathrm{ie}},\qquad
n_{\mathrm i}\nu_{\mathrm{in}} = n_{\mathrm n}\nu_{\mathrm{ni}},
\tag{6}
$$

and the thermal-diffusion forces sum to zero, $A_{\mathrm i}+A_{\mathrm a}=C_{\mathrm e}$ (the total
force the temperature gradient exerts on the electrons is exactly balanced by its reaction on the
two heavy species). Adding (3)–(5) and using (6) removes every collisional and
thermal-diffusion term:

$$
0 = e(n_{\mathrm i}-n_{\mathrm e})\boldsymbol E - \nabla p,
\tag{7}
$$

the bulk force balance between the space-charge electric force and the total pressure gradient.
Because it contains no velocity, (3)–(5) determine only two independent relative
velocities; the common translational velocity of all three species (the "frame") is left
undetermined by the momentum balances alone and must be fixed separately.

**The heavy-particle frame.** The frame is fixed by heavy-particle (ion + neutral) mass
continuity. Summing the ion and neutral continuity equations makes the reaction source $S$
cancel, $\partial_t n_{\mathrm h}+\nabla\cdot(n_{\mathrm i}\boldsymbol u_{\mathrm i}+n_{\mathrm n}\boldsymbol u_{\mathrm n})=0$;
writing the heavy mass flux as $\rho\boldsymbol u$ (mass-averaged velocity, electron mass
neglected) gives

$$
n_{\mathrm i}(\boldsymbol u_{\mathrm i}-\boldsymbol u) + n_{\mathrm n}(\boldsymbol u_{\mathrm n}-\boldsymbol u) = 0
\quad\Longleftrightarrow\quad
\boldsymbol u_{\mathrm n} = \boldsymbol u - \frac{n_{\mathrm i}}{n_{\mathrm n}}(\boldsymbol u_{\mathrm i}-\boldsymbol u).
\tag{8}
$$

This closes the system: two of (3)–(5) determine the electron and ion velocities relative
to the heavy gas, and (8) determines the neutral velocity once the heavy-gas frame
velocity $\boldsymbol u$ itself is known (from the heavy momentum equation, §2.6, or
simply zero if the heavy gas is taken to be at rest).

**Diffusion velocities.** Writing $\boldsymbol V_s=\boldsymbol u_s-\boldsymbol u$ for the velocity
of species $s$ relative to the heavy-gas frame, equations (3)–(5) are unchanged in form
with $\boldsymbol u_s\to\boldsymbol V_s$ (the frame velocity cancels out of every friction term,
which depends only on velocity differences), and (8) becomes
$n_{\mathrm i}\boldsymbol V_{\mathrm i}+n_{\mathrm n}\boldsymbol V_{\mathrm n}=0$. The frame velocity
re-enters only through $\boldsymbol\Gamma_s=n_s\boldsymbol u=n_s\boldsymbol u+n_s\boldsymbol V_s$.

**Electron flux.** Let $\nu_{\mathrm e}=\nu_{\mathrm{en}}+\nu_{\mathrm{ei}}^{\mathrm{eff}}$, so that

$$
\mu_{\mathrm e} = -\frac{e}{m_{\mathrm e}\nu_{\mathrm e}} < 0,\qquad
D_{\mathrm e} = \frac{k_{\mathrm B}T_{\mathrm e}}{m_{\mathrm e}\nu_{\mathrm e}} = \frac{|\mu_{\mathrm e}|k_{\mathrm B}T_{\mathrm e}}{e},
\tag{9}
$$

and let $w_{\mathrm{en}}=\nu_{\mathrm{en}}/\nu_{\mathrm e}$, $w_{\mathrm{ei}}=\nu_{\mathrm{ei}}^{\mathrm{eff}}/\nu_{\mathrm e}=1-w_{\mathrm{en}}$
be the electron's friction weights on the neutrals and the ions. Solving (3) for
$\boldsymbol V_{\mathrm e}$, expanding $\nabla(n_{\mathrm e}T_{\mathrm e})=T_{\mathrm e}\nabla n_{\mathrm e}+n_{\mathrm e}\nabla T_{\mathrm e}$,
multiplying by $n_{\mathrm e}$ and adding back $n_{\mathrm e}\boldsymbol u$ gives the electron flux
in the form the code evaluates:

$$
\boxed{\;
\boldsymbol\Gamma_{\mathrm e} = \mu_{\mathrm e}n_{\mathrm e}\boldsymbol E
 - D_{\mathrm e}\nabla n_{\mathrm e}
 - (1+C_{\mathrm e})\,D_{\mathrm e}\,\frac{n_{\mathrm e}}{T_{\mathrm e}}\nabla T_{\mathrm e}
 + n_{\mathrm e}\Big[\boldsymbol u + w_{\mathrm{en}}(\boldsymbol u_{\mathrm n}-\boldsymbol u) + w_{\mathrm{ei}}(\boldsymbol u_{\mathrm i}-\boldsymbol u)\Big].
\;}
\tag{10}
$$

The first two terms are ordinary drift and diffusion, with $D_{\mathrm e}$ the Einstein
diffusivity. The third is the "unit" thermal-diffusion term, $-(|\mu_{\mathrm e}|/e)\nabla p_{\mathrm e}$
split into a density and a temperature part, plus the additional $C_{\mathrm e}$ contribution from
the Soret (thermal-diffusion) force in (3). The bracket is the collisional drag exerted by
the other two species, weighted by the electron's own friction fractions.

**Ion flux.** Let $\nu_{\mathrm i}=\nu_{\mathrm{in}}+\nu_{\mathrm{ie}}$, $r=\nu_{\mathrm{ie}}/\nu_{\mathrm{in}}$,
so that the ion's *total* mobility (both collision partners) and the friction weights on the
neutral and electron are

$$
\mu_{\mathrm{tot}} = \frac{e}{M\nu_{\mathrm i}} = \frac{\mu_{\mathrm{in}}}{1+r},\qquad
w_{\mathrm{in}} = \frac1{1+r},\qquad
w_{\mathrm{ie}} = \frac{r}{1+r},
\tag{11}
$$

with $\mu_{\mathrm{in}}=e/(M\nu_{\mathrm{in}})$ the ion–neutral (charge-exchange) mobility, §3.7. The
same steps applied to (4) give an ion flux in which the drag on the neutral,
$w_{\mathrm{in}}\boldsymbol V_{\mathrm n}$, is *algebraic* in $\boldsymbol V_{\mathrm i}$ through the frame
condition $\boldsymbol V_{\mathrm n}=-(n_{\mathrm i}/n_{\mathrm n})\boldsymbol V_{\mathrm i}$. Solving for
$\boldsymbol V_{\mathrm i}$ and defining

$$
\alpha = \frac1{1+w_{\mathrm{in}}\,n_{\mathrm i}/n_{\mathrm n}},\qquad
\mu_{\mathrm i} = \alpha\,\mu_{\mathrm{tot}},\qquad
D_{\mathrm i} = \mu_{\mathrm i}\frac{k_{\mathrm B}T_{\mathrm g}}{e},
\tag{12}
$$

the ion flux is

$$
\boxed{\;
\boldsymbol\Gamma_{\mathrm i} = \mu_{\mathrm i}n_{\mathrm i}\boldsymbol E
 - D_{\mathrm i}\nabla n_{\mathrm i}
 - \mu_{\mathrm i}n_{\mathrm i}\frac{k_{\mathrm B}}{e}\nabla T_{\mathrm g}
 + A_{\mathrm i}\,\mu_{\mathrm i}n_{\mathrm e}\frac{k_{\mathrm B}}{e}\nabla T_{\mathrm e}
 + n_{\mathrm i}\Big[\boldsymbol u + \alpha\,w_{\mathrm{ie}}(\boldsymbol u_{\mathrm e}-\boldsymbol u)\Big].
\;}
\tag{13}
$$

Every implicit ion coefficient carries the back-reaction factor $\alpha\le1$ (the ion loses
mobility because it must also drag the neutral gas along with it through the frame
condition); only the cross term with the electron velocity remains explicit.

**Neutral flux.** Multiplying (8) by $n_{\mathrm n}$,

$$
\boldsymbol\Gamma_{\mathrm n} = n_{\mathrm n}\boldsymbol u_{\mathrm n} = n_{\mathrm h}\boldsymbol u - \boldsymbol\Gamma_{\mathrm i}.
\tag{14}
$$

At $\boldsymbol u=0$ this is simply $\boldsymbol\Gamma_{\mathrm n}=-\boldsymbol\Gamma_{\mathrm i}$: every ion
that reaches the wall returns as a neutral atom, and no independent neutral flux need ever be
evaluated.

**Limits.** With the collisional-drag terms switched off, (10) and (13) collapse to
the ordinary drift–diffusion fluxes $\boldsymbol\Gamma_{\mathrm e}=\mu_{\mathrm e}n_{\mathrm e}\boldsymbol E-D_{\mathrm e}\nabla n_{\mathrm e}-(1+C_{\mathrm e})D_{\mathrm e}(n_{\mathrm e}/T_{\mathrm e})\nabla T_{\mathrm e}$,
$\boldsymbol\Gamma_{\mathrm i}=\mu_{\mathrm i}n_{\mathrm i}\boldsymbol E-D_{\mathrm i}\nabla n_{\mathrm i}-\mu_{\mathrm i}n_{\mathrm i}(k_{\mathrm B}/e)\nabla T_{\mathrm g}+A_{\mathrm i}\mu_{\mathrm i}n_{\mathrm e}(k_{\mathrm B}/e)\nabla T_{\mathrm e}$,
with $\alpha=1$. With the electron thermal-diffusion switch also off, $C_{\mathrm e}=A_{\mathrm i}=0$
and the classical drift–diffusion pair is recovered.

### 2.3 Electrostatic field

$$
\nabla\cdot(\varepsilon_0\nabla\phi) = -e(n_{\mathrm i}-n_{\mathrm e}),\qquad
\boldsymbol E = -\nabla\phi.
\tag{15}
$$

At an electrode boundary held at fixed current, the applied current density $J_{\mathrm{app}}$
is enforced through a time-integrated boundary condition on the field itself: the total current
there, conduction plus displacement, $\boldsymbol j+\varepsilon_0\partial_t\boldsymbol E$, is held
equal to $J_{\mathrm{app}}$, and the field value used to seed that ODE is only an initial
condition.

### 2.4 Electron energy

$$
\partial_t\mathcal E + \nabla\cdot\boldsymbol h_{\mathrm e}
= -e\,\boldsymbol\Gamma_{\mathrm e}\cdot\boldsymbol E
 - \frac32 n_{\mathrm e}k_{\mathrm B}(T_{\mathrm e}-T_{\mathrm g})\,\frac{2m_{\mathrm e}}{M}\,\nu_{\mathrm{ex}}
 - \varepsilon_{\mathrm i}\,S,
\qquad \varepsilon_{\mathrm i}=15.76\ \mathrm{eV},
\tag{16}
$$

$$
\boldsymbol h_{\mathrm e} = \frac52 k_{\mathrm B}T_{\mathrm e}\,\boldsymbol\Gamma_{\mathrm e}
 - \kappa_{\mathrm e}\nabla T_{\mathrm e}
 + k_{\mathrm B}T_{\mathrm e}\,n_{\mathrm e}\Big[A_{\mathrm i}(\boldsymbol u_{\mathrm e}-\boldsymbol u_{\mathrm i}) + A_{\mathrm a}(\boldsymbol u_{\mathrm e}-\boldsymbol u_{\mathrm n})\Big].
\tag{17}
$$

The right-hand side of (16) is Joule heating by the electric field acting on the total
(drift + diffusion) electron flux, elastic energy loss to both heavy species at the
Maxwellian-averaged exchange frequency $\nu_{\mathrm{ex}}$ (§3.8), and the ionization-energy
sink of the *net* ionization rate (three-body recombination returns the energy it removed).
There is no radiative loss term. The heat flux (17) is an enthalpy (advective) flux at
$5/2\,k_{\mathrm B}T_{\mathrm e}$ per electron, Fourier conduction with the electron thermal
conductivity $\kappa_{\mathrm e}$ (§3.6), and a thermal-diffusion (Onsager reciprocal) heat flux
carrying the same coefficients $A_{\mathrm i},A_{\mathrm a}$ that appear as forces in (4)–(5),
with $A_{\mathrm i}+A_{\mathrm a}=C_{\mathrm e}$.

Equations (16)–(17) can be rewritten entirely in the evolved energy variable
$\mathcal E=\tfrac32 n_{\mathrm e}k_{\mathrm B}T_{\mathrm e}$. With
$\kappa_{\mathrm e}=K\,n_{\mathrm e}k_{\mathrm B}D_{\mathrm e}$ (§3.6), the advective (enthalpy)
coefficient $h_{\mathrm{adv}}=\tfrac53+\tfrac23C_{\mathrm e}$, and the drift-plus-drag part of the
electron velocity $\boldsymbol v_{\mathrm e}=\mu_{\mathrm e}\boldsymbol E+w_{\mathrm{en}}(\boldsymbol u_{\mathrm n}-\boldsymbol u)+w_{\mathrm{ei}}(\boldsymbol u_{\mathrm i}-\boldsymbol u)$,

$$
\boldsymbol h_{\mathrm e} =
 h_{\mathrm{adv}}\,\mathcal E\,\boldsymbol v_{\mathrm e} + \tfrac53\mathcal E\boldsymbol u
 -\Big[h_{\mathrm{adv}}(1+C_{\mathrm e})+\tfrac23K\Big]D_{\mathrm e}\nabla\mathcal E
 +\Big[h_{\mathrm{adv}}C_{\mathrm e}+\tfrac23K\Big]D_{\mathrm e}\frac{\mathcal E}{n_{\mathrm e}}\nabla n_{\mathrm e}
 -\tfrac23\mathcal E\Big[A_{\mathrm i}(\boldsymbol u_{\mathrm i}-\boldsymbol u)+A_{\mathrm a}(\boldsymbol u_{\mathrm n}-\boldsymbol u)\Big].
\tag{18}
$$

At $C_{\mathrm e}=0$, $K=5/2$ (thermal diffusion and the harmonic conductivity switched off), this
reduces to $\boldsymbol h_{\mathrm e}=\tfrac53\mathcal E\boldsymbol v_{\mathrm e}-\tfrac{10}{3}D_{\mathrm e}\nabla\mathcal E+\tfrac53D_{\mathrm e}(\mathcal E/n_{\mathrm e})\nabla n_{\mathrm e}$,
the ordinary $5/2\,k_{\mathrm B}T_{\mathrm e}$-enthalpy, Fourier-conduction electron energy
equation.

### 2.5 Heavy-particle (gas) energy

$$
n_{\mathrm h}c_v\big(\partial_t T_{\mathrm g}+\boldsymbol u\cdot\nabla T_{\mathrm g}\big)
= \nabla\cdot(\kappa_{\mathrm{hp}}\nabla T_{\mathrm g})
 + \frac32 n_{\mathrm e}k_{\mathrm B}(T_{\mathrm e}-T_{\mathrm g})\frac{2m_{\mathrm e}}{M}\nu_{\mathrm{ex}}
 + Q_{\mathrm{in}}
 - n_{\mathrm h}k_{\mathrm B}T_{\mathrm g}\,\nabla\cdot\boldsymbol u,
\tag{19}
$$

with the elastic-exchange term the exact negative of the electron loss in (16) (same
$\nu_{\mathrm{ex}}$, same local values — the gas gains exactly what the electrons lose), and the
ion–neutral frictional (viscous) heating

$$
Q_{\mathrm{in}} = \frac{e\,n_{\mathrm i}\,|\boldsymbol u_{\mathrm i}-\boldsymbol u_{\mathrm n}|^2}{\mu_{\mathrm{in}}}
\qquad(\to 0\text{ if the ion–neutral drift is not resolved}).
\tag{20}
$$

$Q_{\mathrm{in}}$ is positive definite and reduces to the drift-only Joule dissipation on the ion
channel when the neutral gas is at rest; it is *not* the same as $\boldsymbol E\cdot\boldsymbol j_{\mathrm i}$,
which in a strong sheath field can be negative (the field does work against the ions moving up
their own density gradient). The advective term $\boldsymbol u\cdot\nabla T_{\mathrm g}$ and the
compressive-work term $-n_{\mathrm h}k_{\mathrm B}T_{\mathrm g}\nabla\cdot\boldsymbol u$ are present
only when the heavy gas has its own momentum equation (closure 2 below); they vanish
identically when $\boldsymbol u=0$. There is no ionization-energy term here — it is entirely
charged to the electrons in (16) — and $c_v$, $\kappa_{\mathrm{hp}}$ are given in §3.9.

### 2.6 The neutral density: three closures

The atom (neutral) density is never advanced by its own transport operator by default; it is a
*background* species recovered from the other unknowns. Three closures are available:

**Closure 0 — mechanical equilibrium (Dalton's law).** The neutral density is read off the
total-pressure balance,

$$
n_{\mathrm n} = \frac{p - n_{\mathrm e}k_{\mathrm B}T_{\mathrm e}}{k_{\mathrm B}T_{\mathrm g}} - n_{\mathrm i},
\qquad
p = p_0 + \frac{\varepsilon_0}{2}\big(E^2 - E_{\mathrm{arc}}^2\big),
\tag{21}
$$

with $p_0$ a reference pressure and $E_{\mathrm{arc}}$ the prescribed field at the domain's inflow
boundary. This is the spatial integral, from that boundary, of the bulk momentum balance (7)
including the electrostatic (Maxwell) stress correction; it imposes $\nabla p=e(n_{\mathrm i}-n_{\mathrm e})\boldsymbol E$
pointwise and treats the total pressure as slaved to the local field, i.e. the neutral gas is
assumed to respond to a pressure imbalance instantaneously (infinite sound speed).

**Closure 1 — neutral continuity, heavy-particle-frame flux.** The heavy-particle density is
frozen (conserved cell by cell together with the ion), and the pressure is allowed to float:

$$
\partial_t n_{\mathrm n} + \nabla\cdot\boldsymbol\Gamma_{\mathrm n} = -S,
\qquad
\boldsymbol\Gamma_{\mathrm n} = -\boldsymbol\Gamma_{\mathrm i}.
\tag{22}
$$

Adding this to the ion continuity equation gives $\partial_t n_{\mathrm h}=0$ exactly: this is the
$\boldsymbol u=0$ limit of the general frame condition (8) treated as a genuine flux equation
rather than an algebraic constraint.

**Closure 2 — the heavy-gas momentum equation.** The heavy gas is given its own
mass-averaged velocity $\boldsymbol u$, advanced by an inviscid single-fluid momentum equation:

$$
\partial_t n_{\mathrm n} + \nabla\cdot\boldsymbol\Gamma_{\mathrm n} = -S,\qquad
\boldsymbol\Gamma_{\mathrm n} = n_{\mathrm h}\boldsymbol u - \boldsymbol\Gamma_{\mathrm i},
\tag{23}
$$

$$
\partial_t(\rho\boldsymbol u) + \nabla\cdot(\rho\boldsymbol u\boldsymbol u)
 = -\nabla p + e(n_{\mathrm i}-n_{\mathrm e})\boldsymbol E,
\qquad p = n_{\mathrm e}k_{\mathrm B}T_{\mathrm e} + n_{\mathrm h}k_{\mathrm B}T_{\mathrm g}.
\tag{24}
$$

Under this closure $\partial_t n_{\mathrm h}+\nabla\cdot(n_{\mathrm h}\boldsymbol u)=0$ and the
pressure relaxes acoustically rather than instantaneously; closures 0 and 1 are the two
limits of (24) with the inertial (and hence time-derivative and advective) terms dropped —
the algebraic constraint $\nabla p=e(n_{\mathrm i}-n_{\mathrm e})\boldsymbol E$ must then be either
imposed pointwise (closure 0) or left unenforced (closure 1). No viscosity is included in (24).

### 2.7 Boundary conditions

At the inflow (arc) boundary: Dirichlet values of $n_{\mathrm e}$, $n_{\mathrm i}$, $T_{\mathrm e}$,
$T_{\mathrm g}$ (prescribed profiles or scalars); $n_{\mathrm n}$ from the pressure closure at
$p_0$; the field set by the current-driven condition of §2.3.

At the electrode/wall boundary: electrons leave at the (one-sided Maxwellian) thermal effusion
flux $\tfrac14 n_{\mathrm e}\bar c_{\mathrm e}$ (with $\bar c_{\mathrm e}=\sqrt{8k_{\mathrm B}T_{\mathrm e}/(\pi m_{\mathrm e})}$
the electron mean speed), reduced by any thermionic emission current when enabled, applied as an
inhomogeneous-Neumann (flux) condition; ions leave at their outward drift flux when that flux
already points into the wall, and otherwise at their own thermal effusion flux
$\tfrac14n_{\mathrm i}\bar c_{\mathrm i}$; the electron energy carries $2k_{\mathrm B}T_{\mathrm e}$
per absorbed electron; $T_{\mathrm g}$ is held at the wall temperature; $n_{\mathrm n}$ is
extrapolated (no wall condition, since it is not solved); $\phi=0$; the normal heavy-gas
velocity (closure 2) is zero (impermeable wall).

At a symmetry axis (2-D): reflective (zero normal flux); at an outer radial boundary:
homogeneous Neumann.

---

## 3. Transport coefficients

### 3.1 General framework for electron transport

The electron current density $\boldsymbol J_{\mathrm e}=-e\boldsymbol\Gamma_{\mathrm e}$ and the
electron heat flux $\boldsymbol q_{\mathrm e}$ obey a linear-response (Onsager) relation with the
generalized field $\boldsymbol{\mathcal E}=\boldsymbol E+\nabla p_{\mathrm e}/(en_{\mathrm e})$ and the
temperature gradient as the two driving forces:

$$
\boldsymbol J_{\mathrm e} = \sigma\,\boldsymbol{\mathcal E} + \phi\,\nabla T_{\mathrm e},
\tag{25}
$$

$$
\boldsymbol q_{\mathrm e} = -\frac52\frac{k_{\mathrm B}T_{\mathrm e}}{e}\boldsymbol J_{\mathrm e}
 - \lambda'\nabla T_{\mathrm e} - T_{\mathrm e}\phi\,\boldsymbol{\mathcal E}
 = -\lambda\nabla T_{\mathrm e} + \psi\,\boldsymbol J_{\mathrm e},
\tag{26}
$$

$$
\lambda = \lambda' - T_{\mathrm e}\frac{\phi^2}{\sigma},\qquad
\psi = -\frac52\frac{k_{\mathrm B}T_{\mathrm e}}{e} - T_{\mathrm e}\frac{\phi}{\sigma},
\tag{27}
$$

where $\sigma$ is the electrical conductivity, $\phi$ the thermal-diffusion (Nernst/Seebeck-type)
coefficient, $\lambda'$ the "frozen-field" thermal conductivity (at fixed $\boldsymbol{\mathcal E}$),
and $\lambda$ the ordinary (zero-current) thermal conductivity. Reading (25)–(27) as
the electron momentum balance and energy flux,

$$
\mu_{\mathrm e} = -\frac{\sigma}{e n_{\mathrm e}},\qquad
C_{\mathrm e} = \frac{e}{k_{\mathrm B}}\frac{\phi}{\sigma},\qquad
\kappa_{\mathrm e} = \lambda,\qquad
\Big(\text{enthalpy coefficient} = \frac52+C_{\mathrm e} = -\frac{e\psi}{k_{\mathrm B}T_{\mathrm e}}\Big).
\tag{28}
$$

The conductivity that belongs in the *energy* equation is the zero-current $\lambda$, not
$\lambda'$: a temperature gradient alone drives the current $\phi\nabla T_{\mathrm e}$, and
holding $\boldsymbol J_{\mathrm e}=0$ requires $\boldsymbol{\mathcal E}=-(\phi/\sigma)\nabla T_{\mathrm e}$,
whose own heat flux $-T_{\mathrm e}\phi\boldsymbol{\mathcal E}=+T_{\mathrm e}\phi^2/\sigma\,\nabla T_{\mathrm e}$
reduces the apparent conduction from $\lambda'$ down to $\lambda$; since the energy equation (16)
already separately accounts for $\boldsymbol E\cdot\boldsymbol J$ and the enthalpy flux, its
conductivity must be the zero-current one.

$\sigma$, $\phi$, $\lambda'$, $\lambda$ are evaluated in two limits and bridged, §3.2–3.4.

### 3.2 Partially ionized limit: electron–neutral (Lorentz gas) transport

When electron–neutral collisions dominate over electron–electron and electron–ion collisions,
the electron distribution is a Lorentz gas: the isotropic part of the distribution is a
Maxwellian $f_{\mathrm M}(T_{\mathrm e})$ and only the small anisotropic (first-order) correction
responds to the field and temperature gradient. With $C$ the electron speed,
$\varepsilon=\tfrac12m_{\mathrm e}C^2$ the electron energy, and
$\nu_{\mathrm{eH}}(C)=n_{\mathrm n}\,C\,Q^{(1)}(\varepsilon)$ the (energy-dependent)
electron–neutral momentum-transfer collision frequency built from the momentum-transfer cross
section $Q^{(1)}(\varepsilon)$, a first-order (Chapman–Enskog-type) solution of the electron
Boltzmann equation gives the anisotropic correction

$$
f^1 = \frac{e}{m_{\mathrm e}\nu_{\mathrm{eH}}}\frac{\partial f_{\mathrm M}}{\partial C}\boldsymbol{\mathcal E}
 - \frac{C}{\nu_{\mathrm{eH}}}f_{\mathrm M}\nabla\ln(n_{\mathrm e}f_{\mathrm M}),
\tag{29}
$$

using $\nabla\ln(n_{\mathrm e}f_{\mathrm M})=\nabla\ln p_{\mathrm e}+(x-\tfrac52)\nabla\ln T_{\mathrm e}$,
$x=m_{\mathrm e}C^2/(2k_{\mathrm B}T_{\mathrm e})$, for the Maxwellian case. Substituting into the
current $\boldsymbol J_{\mathrm e}=-e n_{\mathrm e}\tfrac{4\pi}{3}\int C^3 f^1\,dC$ and the conduction
heat flux $\tfrac{2\pi n_{\mathrm e}m_{\mathrm e}}{3}\int(C^2-\tfrac{5k_{\mathrm B}T_{\mathrm e}}{m_{\mathrm e}})C^3 f^1\,dC$
and reading off the coefficients of $\boldsymbol{\mathcal E}$ and $\nabla T_{\mathrm e}$ in
(25)–(26) gives, with all integrals sharing the common weight $C^4f_{\mathrm M}/\nu_{\mathrm{eH}}$,

$$
\sigma = \frac{4\pi n_{\mathrm e}e^2}{3k_{\mathrm B}T_{\mathrm e}}\int_0^\infty\frac{C^4f_{\mathrm M}}{\nu_{\mathrm{eH}}}\,dC,
\qquad
\phi = \frac{4\pi e n_{\mathrm e}}{3T_{\mathrm e}}\int_0^\infty\frac{C^4f_{\mathrm M}}{\nu_{\mathrm{eH}}}\Big(x-\frac52\Big)dC,
\tag{30}
$$

$$
\lambda' = \frac{2\pi n_{\mathrm e}m_{\mathrm e}}{3T_{\mathrm e}}\int_0^\infty\Big(C^2-\frac{5k_{\mathrm B}T_{\mathrm e}}{m_{\mathrm e}}\Big)\frac{C^4f_{\mathrm M}}{\nu_{\mathrm{eH}}}\Big(x-\frac52\Big)dC
 = n_{\mathrm e}k_{\mathrm B}D_{\mathrm{en}}\Big\langle\Big(x-\frac52\Big)^2\Big\rangle,
\tag{31}
$$

$$
\lambda = \lambda' - T_{\mathrm e}\frac{\phi^2}{\sigma}
 = n_{\mathrm e}k_{\mathrm B}D_{\mathrm{en}}\big[\langle x^2\rangle - \langle x\rangle^2\big],
\tag{32}
$$

where $\langle\cdot\rangle$ denotes the average under the weight $C^4f_{\mathrm M}/\nu_{\mathrm{eH}}$
(equivalently, in energy, weight $\varepsilon\,e^{-\varepsilon/k_{\mathrm B}T_{\mathrm e}}/Q^{(1)}(\varepsilon)$),
and

$$
D_{\mathrm{en}} = \frac{4\pi}3\int_0^\infty\frac{C^4f_{\mathrm M}}{\nu_{\mathrm{eH}}}\,dC = \frac{\sigma\,k_{\mathrm B}T_{\mathrm e}}{n_{\mathrm e}e^2}
\tag{33}
$$

is the Lorentz electron–neutral diffusivity (so $D_{\mathrm{en}}/|\mu_{\mathrm{en}}|=k_{\mathrm B}T_{\mathrm e}/e$
exactly, the Einstein relation). The Maxwellian-averaged electron–neutral momentum-transfer
(friction) frequency is

$$
\bar\nu_{\mathrm{ea}} = \frac{\langle C^2\nu_{\mathrm{eH}}\rangle_{f_{\mathrm M}}}{\langle C^2\rangle_{f_{\mathrm M}}}
 = \frac43\,n_{\mathrm n}\bar c_{\mathrm e}\,\bar Q^{(1,1)}_{\mathrm{ea}}(T_{\mathrm e}),
\qquad
\bar c_{\mathrm e}=\sqrt{\frac{8k_{\mathrm B}T_{\mathrm e}}{\pi m_{\mathrm e}}},
\tag{34}
$$

evaluated under the *same* $C^2$ weight (not the weight of (30)–(32)): this is the rate at
which a slowly drifting Maxwellian is damped by friction, $-m_{\mathrm e}n_{\mathrm e}\bar\nu_{\mathrm{ea}}\boldsymbol u$,
and — because the same weight also governs energy transfer — it is simultaneously the elastic
energy-loss rate, $\tfrac{2m_{\mathrm e}}{M}\langle\varepsilon\nu_{\mathrm{eH}}\rangle(1-T_{\mathrm g}/T_{\mathrm e})=\tfrac32k_{\mathrm B}(T_{\mathrm e}-T_{\mathrm g})\tfrac{2m_{\mathrm e}}{M}\bar\nu_{\mathrm{ea}}$
used in §2.4–2.5. It is emphatically **not** the mobility (harmonic-mean) collision frequency
$e/(m_{\mathrm e}|\mu_{\mathrm{en}}|)=\langle C^2\rangle/\langle C^2/\nu_{\mathrm{eH}}\rangle$: for
argon between 0.5 and 2 eV the two differ by a factor of 4–4.5, because the energy-exchange rate
is dominated by fast electrons above the Ramsauer minimum in $Q^{(1)}(\varepsilon)$, while the
mobility is dominated by the slow electrons trapped inside it.

Equations (30)–(34) are evaluated once, off-line, on the tabulated argon
electron–neutral momentum-transfer cross section and stored as functions of $T_{\mathrm e}$ per
unit neutral density (for $\sigma,\phi,\lambda',\lambda$) or per neutral density directly (for
$\bar\nu_{\mathrm{ea}}$), for interpolation at run time.

### 3.3 Strongly ionized limit and the Coulomb logarithm

When electron–ion (Coulomb) collisions dominate, the classical strongly-ionized
(Spitzer–Härm) solution for a singly charged ion species gives

$$
\sigma = 1.975\,\frac{n_{\mathrm e}e^2}{m_{\mathrm e}\bar\nu_{\mathrm{ei}}},\qquad
\lambda = 3.203\,\frac{n_{\mathrm e}k_{\mathrm B}^2T_{\mathrm e}}{m_{\mathrm e}\bar\nu_{\mathrm{ei}}},\qquad
\phi = 1.389\,\frac{n_{\mathrm e}ek_{\mathrm B}}{m_{\mathrm e}\bar\nu_{\mathrm{ei}}},\qquad
\psi = -3.203\,\frac{k_{\mathrm B}T_{\mathrm e}}{e},
\tag{35}
$$

so that $(e/k_{\mathrm B})\phi/\sigma=1.389/1.975=0.7033$ (the collision frequency cancels), and the
consistency check of (28), $-e\psi/(k_{\mathrm B}T_{\mathrm e})=3.203=\tfrac52+0.703$, holds
exactly. Here $\bar\nu_{\mathrm{ei}}$ is the Maxwellian-averaged electron–ion momentum-transfer
frequency,

$$
\bar\nu_{\mathrm{ei}} = n_{\mathrm i}\,\frac{4\sqrt{2\pi}}{3}\Big(\frac{m_{\mathrm e}}{k_{\mathrm B}T_{\mathrm e}}\Big)^{3/2}
 \Big(\frac{e^2}{4\pi\varepsilon_0 m_{\mathrm e}}\Big)^2\ln\Lambda
 = 2.9063\times10^{-12}\,\frac{n_{\mathrm i}\,\ln\Lambda}{T_{\mathrm e}[\mathrm{eV}]^{3/2}}\ \ \mathrm{s^{-1}},
\tag{36}
$$

$$
\ln\Lambda = \max\Big(2,\ \ln\frac{\lambda_{\mathrm D}}{b_{90}}\Big),\qquad
\lambda_{\mathrm D}=\sqrt{\frac{\varepsilon_0 T_{\mathrm e}[\mathrm{eV}]}{n_{\mathrm e}e}},\qquad
b_{90}=\frac{e}{12\pi\varepsilon_0 T_{\mathrm e}[\mathrm{eV}]},
\tag{37}
$$

with $\lambda_{\mathrm D}$ the electron Debye length and $b_{90}$ the classical 90°-deflection
impact parameter; the logarithm is floored at 2 to keep $\bar\nu_{\mathrm{ei}}$ finite when the
plasma is not weakly coupled. Note that the rate (36) is set by the *ion* density (the
collision partner), while the Debye length in (37) uses the *electron* density; the two
coincide only under quasi-neutrality, which need not hold pointwise (e.g. inside a sheath).
$\bar\nu_{\mathrm{ei}}=0$ when $T_{\mathrm e}$ is below a numerical floor or when either species is
absent.

### 3.4 Bridging the partially and strongly ionized limits

Both limits are expressed in the same three quantities $\sigma,\phi,\lambda$ (or, equivalently,
$C_{\mathrm e}=(e/k_{\mathrm B})\phi/\sigma$ and $\lambda$ directly), so the physical transport
coefficient at intermediate ionization is obtained by an additive (resistance-type) combination
of the electron–neutral and electron–ion contributions, in the ratio of their momentum-transfer
frequencies,

$$
P = \frac{\bar\nu_{\mathrm{ea}}}{\bar\nu_{\mathrm{ei}}}.
\tag{38}
$$

**Thermal-diffusion coefficient.** The electron–neutral and electron–ion values of $\phi/\sigma$
are combined as a weighted mean in $1/P$ (equivalently, as two resistances in a mixture rule
that is exact in both limits, $P\to\infty$ recovering the electron–neutral value and $P\to0$ the
electron–ion value):

$$
\Big(\frac{\phi}{\sigma}\Big)_{\mathrm e} = \frac{(\phi/\sigma)_{\mathrm{ei}}/P + (\phi/\sigma)_{\mathrm{ea}}}{1+1/P},
\qquad
\Big(\frac{\phi}{\sigma}\Big)_{\mathrm{ei}} = \frac{1.389}{1.975}\frac{k_{\mathrm B}}{e},
\tag{39}
$$

$$
C_{\mathrm e} = \frac{e}{k_{\mathrm B}}\Big(\frac{\phi}{\sigma}\Big)_{\mathrm e},
\qquad
A_{\mathrm i} = \frac{e}{k_{\mathrm B}}\frac{(\phi/\sigma)_{\mathrm{ei}}/P}{1+1/P},
\qquad
A_{\mathrm a} = C_{\mathrm e} - A_{\mathrm i}.
\tag{40}
$$

$A_{\mathrm i}$ is the ion's momentum-conserving reaction to the electron–ion (Coulomb) part of the
thermal-diffusion force alone, since only the electron–ion collision directly couples momentum
between electrons and ions; the electron–neutral part reacts on the neutral instead, giving
$A_{\mathrm a}=C_{\mathrm e}-A_{\mathrm i}$, which appears in no separate equation because the
neutral has no momentum operator of its own.

**Thermal conductivity.** The zero-current conductivities combine as parallel thermal
resistances (an additive-inverse mixture rule, again exact in both limits):

$$
\frac1{\lambda_{\mathrm e}} = \frac1{\lambda_{\mathrm{ea}}} + \frac1{\lambda_{\mathrm{ei}}},
\qquad
\kappa_{\mathrm e} = \lambda_{\mathrm e} = K\,n_{\mathrm e}k_{\mathrm B}D_{\mathrm e},
\tag{41}
$$

with $D_{\mathrm e}$ the total (bridged) diffusivity of (42) below, so that the dimensionless
factor $K$ is independent of $n_{\mathrm e}$ (both $\lambda_{\mathrm{ea}}$ and $\lambda_{\mathrm{ei}}$
scale linearly with it). In the electron–neutral-only limit $K=\lambda_{\mathrm{ea}}/(n_{\mathrm e}k_{\mathrm B}D_{\mathrm e})$;
in the fully ionized limit $K=3.203/2=1.60$.

### 3.5 What the code forms: putting the pieces together

The mobility and diffusivity actually used in the fluxes (10), (13) are

$$
\mu_{\mathrm{en}} = -\frac{\sigma_{\mathrm{ea}}(T_{\mathrm e})}{e\,n_{\mathrm e}}\Big|_{\text{electron--neutral only}},
\qquad
D_{\mathrm{en}} = |\mu_{\mathrm{en}}|\frac{k_{\mathrm B}T_{\mathrm e}}{e},
\qquad
\nu_{\mathrm{en}} = \frac{e}{m_{\mathrm e}|\mu_{\mathrm{en}}|},
\tag{42}
$$

$$
c_f = \frac{\nu_{\mathrm{en}}}{\nu_{\mathrm{en}}+0.5\,\bar\nu_{\mathrm{ei}}}
\qquad\big(=w_{\mathrm{en}}\text{ of }(10)\big),
\qquad
\mu_{\mathrm e} = c_f\,\mu_{\mathrm{en}},\qquad
D_{\mathrm e} = c_f\,D_{\mathrm{en}},
\tag{43}
$$

the same additive-resistivity form applied directly to $\sigma$, with the $0.5$ on
$\bar\nu_{\mathrm{ei}}$ the strongly-ionized conductivity weight of (35). $C_{\mathrm e}$,
$A_{\mathrm i}$, $A_{\mathrm a}$, and $\kappa_{\mathrm e}$ follow from §3.4.

### 3.6 Ion mobility, diffusivity and friction weights

The ion–neutral (charge-exchange) mobility uses the classical resonant charge-exchange cross
section for Ar$^+$ in Ar and its low-field limit corrected for field heating through the
Wannier/Viehland–Mason two-temperature theory:

$$
\bar Q_{\mathrm{ia}} = 2\Big(7\times10^{-10} - 0.6\times10^{-10}\ln\frac{2k_{\mathrm B}T_{\mathrm g}}{e}\Big)^{\!2}\ \mathrm{m^2},
\qquad
\mu_{\mathrm i}^0 = \frac{1.104\times10^{5}}{\bar Q_{\mathrm{ia}}\sqrt{T_{\mathrm g}}\;n_{\mathrm n}},
\tag{44}
$$

$$
v_{\mathrm d}^2 = \frac{\sqrt{1+4AB}-1}{2B},\qquad
A = (\mu_{\mathrm i}^0|\boldsymbol E|)^2,\qquad
B = \frac{M}{3k_{\mathrm B}T_{\mathrm g}},
\tag{45}
$$

$$
T_{\mathrm{eff}} = T_{\mathrm g} + \frac{M v_{\mathrm d}^2}{3k_{\mathrm B}},
\qquad
\mu_{\mathrm{in}} = \mu_{\mathrm i}^0\sqrt{\frac{T_{\mathrm g}}{T_{\mathrm{eff}}}}.
\tag{46}
$$

$v_{\mathrm d}$ is the (self-consistent) mean ion drift speed under field heating, and
$T_{\mathrm{eff}}$ is the ion's effective kinetic temperature including that drift energy; at low
field $v_{\mathrm d}\to\mu_{\mathrm i}^0|\boldsymbol E|$, $T_{\mathrm{eff}}\to T_{\mathrm g}$,
$\mu_{\mathrm{in}}\to\mu_{\mathrm i}^0$.

The ion's electron-friction fraction, needed for $r=\nu_{\mathrm{ie}}/\nu_{\mathrm{in}}$ in
(11)–(12), follows from momentum conservation between the electron and ion momentum
balances, $n_{\mathrm i}M\nu_{\mathrm{ie}}=n_{\mathrm e}m_{\mathrm e}\nu_{\mathrm{ei}}^{\mathrm{eff}}$
with $\nu_{\mathrm{ei}}^{\mathrm{eff}}=(1-c_f)\,e/(m_{\mathrm e}|\mu_{\mathrm e}|)$, and
$\nu_{\mathrm{in}}=e/(M\mu_{\mathrm{in}})$, giving (masses cancelling)

$$
r = (1-c_f)\,\frac{\mu_{\mathrm{in}}\,n_{\mathrm e}}{|\mu_{\mathrm e}|\,n_{\mathrm i}},
\qquad
w_{\mathrm{ie}} = \frac{r}{1+r},\qquad
\mu_{\mathrm{tot}} = \frac{\mu_{\mathrm{in}}}{1+r},
\tag{47}
$$

$$
\alpha = \frac1{1+w_{\mathrm{in}}\,n_{\mathrm i}/n_{\mathrm n}},\qquad
\mu_{\mathrm i} = \alpha\,\mu_{\mathrm{tot}},\qquad
D_{\mathrm i} = \mu_{\mathrm i}\frac{k_{\mathrm B}T_{\mathrm g}}{e}.
\tag{48}
$$

As $n_{\mathrm n}\to0$ (near-total ionization) the ion mobility saturates at the
Coulomb-collision-limited value set by the electron mobility, $\mu_{\mathrm i}\to|\mu_{\mathrm e}|n_{\mathrm i}/((1-c_f)n_{\mathrm e})$,
rather than diverging.

### 3.7 Elastic electron–heavy energy exchange and ion–neutral friction heating

The total elastic momentum/energy-exchange frequency that appears in (16) and (19) is
the sum of the electron–neutral and electron–ion Maxwellian-averaged momentum-transfer
frequencies (the *same* quantities as (34) and (36), used here in full — not the
conductivity-weighted $0.5\bar\nu_{\mathrm{ei}}$ of (35)):

$$
\nu_{\mathrm{ex}} = \bar\nu_{\mathrm{ea}} + \bar\nu_{\mathrm{ei}}.
\tag{49}
$$

Both terms are, by construction, the Maxwellian-averaged elastic energy-loss rate
$\langle\varepsilon\nu\rangle/\langle\varepsilon\rangle$ for a near-Maxwellian electron
population on a cold heavy-particle background, so that the energy the electrons lose in
(16) is exactly the energy the gas gains in (19).

The ion–neutral frictional heating of the gas energy equation, (20), uses the
charge-exchange mobility $\mu_{\mathrm{in}}$ of (46) (not the reduced ion mobility
$\mu_{\mathrm i}$, which already contains the ion's electron-friction correction).

### 3.8 Heavy-particle (Ar, Ar$^+$) heat capacity and thermal conductivity

The heat capacity per heavy particle (translational degrees of freedom only — the electron
energy is carried separately, and the ionization energy is a source term, not a heat-capacity
contribution) depends on whether the neutral closure treats the gas as isobaric or isochoric:

$$
c_v = \frac52 k_{\mathrm B}\ \ (\text{closure 0, isobaric}),
\qquad
c_v = \frac32 k_{\mathrm B}\ \ (\text{closures 1, 2, isochoric}).
\tag{50}
$$

The heavy-particle thermal conductivity is the sum of an atom-dominated and an ion-dominated
channel, each following a Chapman–Enskog-type first-approximation formula built from the
appropriate collision cross sections:

$$
\kappa_{\mathrm{hp}} = \kappa_{\mathrm a} + \kappa_{\mathrm i},
\tag{51}
$$

$$
\kappa_{\mathrm a} = \frac{75k_{\mathrm B}}{64}\sqrt{\frac{\pi k_{\mathrm B}T_{\mathrm g}}{M}}\,
 \frac{n_{\mathrm n}}{n_{\mathrm n}\bar Q_{\mathrm{aa}}+n_{\mathrm i}\bar Q_{\mathrm{ia}}},
\qquad
\kappa_{\mathrm i} = \frac{75k_{\mathrm B}}{64}\sqrt{\frac{\pi k_{\mathrm B}T_{\mathrm g}}{M}}\,
 \frac{n_{\mathrm i}}{n_{\mathrm n}\bar Q_{\mathrm{ia}}+n_{\mathrm i}\bar Q_{\mathrm{ii}}},
\tag{52}
$$

$$
\bar Q^{(2,2)}_{\mathrm{aa}} = 1.12\times10^{-18}\,T_{\mathrm g}^{-0.2},\qquad
\bar Q^{(2,2)}_{\mathrm{ia}} = 3.6\times10^{-18}\,T_{\mathrm g}^{-0.3},\qquad
\bar Q^{(2,2)}_{\mathrm{ii}} = \frac{e^4\ln\Lambda}{36\pi\varepsilon_0^2(k_{\mathrm B}T_{\mathrm g})^2}
\quad[\mathrm{m^2},\ T_{\mathrm g}\text{ in K}],
\tag{53}
$$

with $\ln\Lambda$ the same Coulomb logarithm as (37) (evaluated on $n_{\mathrm e}$, $T_{\mathrm e}$).
This form is regular in both limits: $\kappa_{\mathrm i}\to$ a pure ion–ion Coulomb value as
$n_{\mathrm n}\to0$, and $\kappa_{\mathrm a}\to$ the pure-gas atom–atom value as $n_{\mathrm i}\to0$.

### 3.9 Reaction rate coefficients

The electron-impact ionization rate $k_{\mathrm i}(T_{\mathrm e})$ and the three-body
recombination rate $k_{\mathrm r}(T_{\mathrm e})$ of (2) are tabulated effective rate
coefficients, log-linearly interpolated in $T_{\mathrm e}$ (eV). The ionization energy charged
to the electron energy equation per net ionization event is $\varepsilon_{\mathrm i}=15.76$ eV,
the first ionization potential of argon.

---

## Summary of the closed system

The state advanced in time is $\{n_{\mathrm e},n_{\mathrm i},\mathcal E,T_{\mathrm g},\phi,(\boldsymbol u)\}$,
governed by (1)–(2) with fluxes (10), (13); Poisson's equation
(15); the electron energy equation (16)–(18); the heavy-particle energy
equation (19)–(20); and, depending on the selected closure, either the algebraic
relation (21) or the additional transport equations (22)–(24) for the
neutral density (and, under closure 2, the heavy-gas velocity). Every transport coefficient
entering these equations — $\mu_{\mathrm e},D_{\mathrm e},C_{\mathrm e},\kappa_{\mathrm e}$ for the
electron; $\mu_{\mathrm i},D_{\mathrm i},A_{\mathrm i}$ for the ion; $\nu_{\mathrm{ex}},Q_{\mathrm{in}}$
for the interspecies exchange; $c_v,\kappa_{\mathrm{hp}}$ for the heavy gas — is a closed-form or
tabulated function of the local state, given in §3.
