%For claude
%enrgy scale between ions and neutral = 1/(n_i Q_in + n_n Q_nn) Qin ~ Qni 
%energy scale between electrons and heavy = 1/(n_i Q_ei + n_n Q_en) sqrt(2 m_e/M) 
%time scale for recombination = 1/krn^2, so length scale for recombination = J/(e*n) * 1/krn^2 
%
Give me those scales for 10000K, J=1e6 J=1e7 and 20000K J=1e7, J=1e8
---

## Response (Claude, 2026-10-05)

**Inputs.** Argon, $p = 1$ atm, LTE ($T_e = T_h = T$). Densities from PlasmeX `exe/ar_eq.dat` (Hoffert–Lien column; Saha differs by < 2 %), with $n_n = p/k_BT - 2n_e$. Rates $k_i, k_r$ from `hoff_f.dat` / `hoff_r.dat` (the code default). $\overline Q_{in}$ = Benilov–Naidis fit (the one in the paper), $\overline Q_{nn} = 1.12\times10^{-18}\,T^{-0.2}$ (the $Q_{aa}$ of `src/heavy_transport.H`). $\overline Q_{en}$ = Maxwellian energy-exchange average of `ar_Q_en.dat`, $\overline Q_{ei}$ = energy-exchange Coulomb average of appendix eq. (nu-ei). Script: scratchpad `scales.py` (not saved in the repo).

| | $T=10\,000$ K | $T=20\,000$ K |
|---|---|---|
| $n_e = n_i$ [m$^{-3}$] | $1.52\times10^{22}$ | $1.82\times10^{23}$ |
| $n_n$ [m$^{-3}$] | $7.04\times10^{23}$ | $3.77\times10^{21}$ |
| ionization degree | 0.021 | 0.98 |
| $\overline Q_{in}$ [m$^2$] | $8.9\times10^{-19}$ | $7.8\times10^{-19}$ |
| $\overline Q_{nn}$ [m$^2$] | $1.8\times10^{-19}$ | $1.5\times10^{-19}$ |
| $\overline Q_{ei}$ [m$^2$] ($\ln\Lambda$) | $2.0\times10^{-17}$ (4.6) | $4.8\times10^{-18}$ (4.4) |
| $\overline Q_{en}$ [m$^2$] | $3.7\times10^{-20}$ | $7.6\times10^{-20}$ |

### 1. Ion–neutral energy scale — formula is **not quite right**

$1/(n_i\overline Q_{in} + n_n\overline Q_{nn})$ is the mean free path of a **neutral** (against ions and other neutrals). Neutral–neutral collisions only redistribute energy among the neutrals; they do not move any energy between ions and neutrals, so they cannot relax $T_i - T_n$. Writing the exchange $\dot q_{n\leftarrow i} \propto n_in_n\overline Q_{in}\bar c\,(T_i-T_n)$ into both energy equations, the difference relaxes at

$$
\frac{d(T_i-T_n)}{dt} = -(\nu_{in}+\nu_{ni})(T_i-T_n),\qquad \nu_{in}+\nu_{ni}=(n_n+n_i)\,\overline Q_{in}\,\bar c ,
$$

and since $M_i = M_n$ a collision transfers an $O(1)$ fraction of the energy. So the correct scale is

$$
\lambda_{in}^{E} \approx \frac{1}{(n_i+n_n)\,\overline Q_{in}} ,
$$

which is exactly the $\lambda_i$ already in the paper. ($Q_{in}\sim Q_{ni}$ is fine: it is the same collision.)

| | $10\,000$ K | $20\,000$ K |
|---|---|---|
| your form, $1/(n_i\overline Q_{in}+n_n\overline Q_{nn})$ | $7.2\times10^{-6}$ m | $7.0\times10^{-6}$ m |
| **correct**, $1/((n_i+n_n)\overline Q_{in})$ | $\mathbf{1.6\times10^{-6}}$ **m** | $\mathbf{6.9\times10^{-6}}$ **m** |

At 10 000 K your form is 4.6× too long, because it replaces the dominant $n_n\overline Q_{in}$ by the smaller $n_n\overline Q_{nn}$. At 20 000 K the two agree because $n_i$ dominates.

⚠ **Paper check (not changed):** `model.tex` line 13 says this gives $\lambda_i\sim10^{-7}$ m. The formula gives $1.6\times10^{-6}$ m at 10 000 K and $7\times10^{-6}$ m at 20 000 K, so it is $10^{-6}$–$10^{-5}$ m. The paper's own $\lambda_e$ formula (line 9) gives $8.5\times10^{-7}$ m and $3.0\times10^{-7}$ m. That means **$\lambda_i > \lambda_e$**, not $\lambda_i \ll \lambda_e$ as written in the hierarchy on line 15. The conclusion ($\lambda_i \ll \lambda_T \sim 10^{-4}$ m, so $T_i = T_n$ holds) is unaffected.

### 2. Electron–heavy energy scale — formula is right if $\sqrt{2m_e/M}$ **divides**

As you wrote it, it reads like multiplication, which would give $\sim10^{-8}$ m, a meaningless number. The intended (and paper's) form is

$$
\lambda_T \approx \frac{\lambda_{eh}}{\sqrt{2m_e/M}},\qquad \lambda_{eh}=\frac{1}{n_i\overline Q_{ei}+n_n\overline Q_{en}},\qquad \sqrt{2m_e/M}=5.24\times10^{-3}.
$$

The square root comes from a diffusive estimate: an electron loses energy after $\sim M/2m_e$ collisions, so it random-walks a distance $\lambda_{eh}\sqrt{M/2m_e}$. A sharper version of the same estimate uses separate momentum and energy frequencies, $\lambda_T=\sqrt{D_e\,\tau_E}=\sqrt{k_BT_e/(m_e\bar\nu^M_e)\,/\,(\tfrac{2m_e}{M}\bar\nu^E_{eh})}$. It gives the same order and is shown as a check:

| | $10\,000$ K | $20\,000$ K |
|---|---|---|
| $\lambda_{eh}$ (code averages) | $3.0\times10^{-6}$ m | $1.1\times10^{-6}$ m |
| **$\lambda_T=\lambda_{eh}/\sqrt{2m_e/M}$** | $\mathbf{5.7\times10^{-4}}$ **m** | $\mathbf{2.2\times10^{-4}}$ **m** |
| same, with the paper's $\overline Q_{en}=2\times10^{-20}$ and MK $\overline Q_{ei}$ | $4.5\times10^{-4}$ m | $1.6\times10^{-4}$ m |
| check: $\sqrt{D_e\tau_E}$ | $3.9\times10^{-4}$ m | $1.4\times10^{-4}$ m |

Electron–ion collisions dominate at both temperatures ($\bar\nu^E_{ei}/\bar\nu^E_{en}\approx 12$ at 10 000 K and $3000$ at 20 000 K), so $\lambda_T$ is mostly a Coulomb number. This is consistent with the paper's "$\sim10^{-4}$ m".

### 3. Recombination scale

**Time.** $\tau_{rec}=1/(k_rn_e^2)$ is correct as the **recombination time** (three-body recombination with an electron as the third body: loss per ion $=k_rn_en_e$). However, the time to relax a departure from ionization equilibrium is set by the net rate. Linearizing $\dot n_e=k_in_en_n-k_rn_e^3$ about equilibrium at fixed $p$ and $T$ (so $\delta n_n=-2\,\delta n_e$) gives

$$
\tau_{chem}=\frac{1}{2k_rn_e^2+2k_in_e} .
$$

At 10 000 K, $\tau_{chem}\approx\tau_{rec}/2$. At 20 000 K the ionization term is 48× larger than the recombination term, because only 2 % neutrals are left and $k_i n_e$ is fast. So $\tau_{rec}$ overestimates the relaxation time there by a factor of ~100.

**Length.** $J/(en_e)\cdot\tau$ is the distance an *individual electron* drifts before it recombines. It is **not** the thickness of the non-equilibrium layer. In 1-D steady state, $\nabla\cdot\mathbf J=0$ and quasi-neutrality give $d\Gamma_e/dx=d\Gamma_i/dx$. The current-carrying part of the flux is divergence-free, so it cannot create or carry away a density departure. With $\Gamma_i=\mu_iJ/[e(\mu_e+\mu_i)]-D_a\nabla n$, the drift part is constant, and only ambipolar diffusion transports the deviation. The layer scale is therefore the Damköhler (ionization-layer) length

$$
\lambda_{rec}=\sqrt{D_a\,\tau_{chem}},\qquad D_a\simeq\frac{k_B(T_e+T_h)}{M\bar\nu_{in}},\quad \bar\nu_{in}\ \text{evaluated with } n_h=n_i+n_n ,
$$

The total heavy density $n_h$ enters because the ion diffuses relative to the heavy gas as a whole. This length does **not** depend on $J$.

| | $10\,000$ K | $20\,000$ K |
|---|---|---|
| $k_r$ [m$^6$/s], $k_i$ [m$^3$/s] | $2.3\times10^{-41}$, $7.6\times10^{-21}$ | $1.2\times10^{-42}$, $1.1\times10^{-17}$ |
| $\tau_{rec}=1/(k_rn_e^2)$ | $1.9\times10^{-4}$ s | $2.5\times10^{-5}$ s |
| $\tau_{chem}$ | $9.1\times10^{-5}$ s | $2.5\times10^{-7}$ s |
| $D_a$ | $3.0\times10^{-3}$ m$^2$/s | $1.9\times10^{-2}$ m$^2$/s |
| **$\sqrt{D_a\tau_{chem}}$** | $\mathbf{5.2\times10^{-4}}$ **m** | $\mathbf{6.9\times10^{-5}}$ **m** |

Your requested $J$-dependent values, $L=\dfrac{J}{en_e}\,\tau$:

| case | $u_e=J/(en_e)$ | $L$ with $\tau_{rec}$ (your formula) | $L$ with $\tau_{chem}$ |
|---|---|---|---|
| 10 000 K, $J=10^6$ A/m$^2$ | 411 m/s | $7.7\times10^{-2}$ m | $3.8\times10^{-2}$ m |
| 10 000 K, $J=10^7$ A/m$^2$ | 4110 m/s | $7.7\times10^{-1}$ m | $3.8\times10^{-1}$ m |
| 20 000 K, $J=10^7$ A/m$^2$ | 344 m/s | $8.6\times10^{-3}$ m | $8.7\times10^{-5}$ m |
| 20 000 K, $J=10^8$ A/m$^2$ | 3440 m/s | $8.6\times10^{-2}$ m | $8.7\times10^{-4}$ m |

At 10 000 K these come out at centimetres to tens of centimetres, which shows they are not the layer thickness. The diffusive $\sqrt{D_a\tau_{chem}}$ gives $\sim5\times10^{-4}$ m, comparable to $\lambda_T$ and to the domain length ($2\times10^{-4}$ m).

**Caveat on the rates.** These numbers use Hoffert–Lien, the code default. Hinnov–Hirschberg ($k_r=5.6\times10^{-39}T_e[\mathrm{eV}]^{-4.5}$, the `saha` option) is 470× faster at 10 000 K and 400× faster at 20 000 K. It would shrink $\tau$ by about that factor and $\sqrt{D_a\tau}$ by about 20×, to roughly $2.5\times10^{-5}$ m at 10 000 K. At 10 000 K, chemical non-equilibrium over the $2\times10^{-4}$ m domain therefore depends on which rate set is used.

### Summary

| scale | formula OK? | $10\,000$ K | $20\,000$ K |
|---|---|---|---|
| ion–neutral energy | use $(n_i+n_n)\overline Q_{in}$, drop $Q_{nn}$ | $1.6\times10^{-6}$ m | $6.9\times10^{-6}$ m |
| electron–heavy energy | yes, **divide** by $\sqrt{2m_e/M}$ | $5.7\times10^{-4}$ m | $2.2\times10^{-4}$ m |
| recombination time | $1/(k_rn_e^2)$ OK as recombination time; relaxation time is $1/(2k_rn_e^2+2k_in_e)$ | $9.1\times10^{-5}$ s | $2.5\times10^{-7}$ s |
| recombination length | $J/(en)\cdot\tau$ is not the layer scale; use $\sqrt{D_a\tau_{chem}}$ | $5.2\times10^{-4}$ m | $6.9\times10^{-5}$ m |

---

## Follow-up: why $D_a$ in Section 3 uses $n_h=n_i+n_n$ and not $n_n$ alone

**Short answer.** An ion cannot diffuse through a gas that stays where it is. In 1-D, with the heavy gas at rest, every ion that moves toward the wall must be balanced by a neutral moving the other way. That is eq. (11) with $\mathbf u=0$: $\boldsymbol\Gamma_n=-\boldsymbol\Gamma_i$. So the ion rubs against neutrals that are moving *towards* it, and the relative speed $u_i-u_n$ is larger than $u_i$ by the factor $n_h/n_n$. The ion–neutral friction per unit ion speed therefore scales with $n_n\times(n_h/n_n)=n_h$, not $n_n$. This is the standard result that a **binary diffusion coefficient depends on the total density of the mixture**, not on the density of one partner.

### Step by step

The ion–neutral friction force on the ions, from eq. (4), is

$$
\mathbf R_{in}=-Mn_i\,\bar\nu_{in}\,(\mathbf u_i-\mathbf u_n),\qquad \bar\nu_{in}=n_n\overline Q_{in}\bar c\ \ (\propto n_n).
$$

**Step 1: neutrals move too.** In steady state the heavy gas is at rest (mechanical equilibrium), so by eq. (8)

$$
n_i\mathbf u_i+n_n\mathbf u_n=n_h\mathbf u=0\quad\Rightarrow\quad \mathbf u_n=-\frac{n_i}{n_n}\mathbf u_i .
$$

**Step 2: relative velocity.**

$$
\mathbf u_i-\mathbf u_n=\mathbf u_i\left(1+\frac{n_i}{n_n}\right)=\mathbf u_i\,\frac{n_h}{n_n}.
$$

**Step 3: friction.**

$$
\mathbf R_{in}=-Mn_i\,(n_n\overline Q_{in}\bar c)\,\frac{n_h}{n_n}\,\mathbf u_i=-Mn_i\,\underbrace{(n_h\overline Q_{in}\bar c)}_{\bar\nu_{in}\,n_h/n_n}\,\mathbf u_i .
$$

The $n_n$ cancels. The effective friction frequency for ion motion *relative to the gas at rest* is $\bar\nu_{in}\,n_h/n_n$.

**Step 4: ambipolar diffusion.** Electrons and ions move together (quasi-neutrality), so the electric field cancels when their force balances are added. Electron–ion friction is internal to the pair, and electron–neutral friction is smaller by $\sim m_e/M$. That leaves

$$
-\nabla(p_e+p_i)=Mn_h\overline Q_{in}\bar c\;\boldsymbol\Gamma_i
\quad\Rightarrow\quad
\boldsymbol\Gamma_i=-D_a\nabla n,\qquad D_a=\frac{k_B(T_e+T_h)}{M\,n_h\overline Q_{in}\bar c}=\frac{k_B(T_e+T_h)}{M\bar\nu_{in}}\cdot\frac{n_n}{n_h}
$$

(temperature gradients ignored). This is what I used: $\bar\nu_{in}$ evaluated with $n_h$ in place of $n_n$.

### The paper's own equations give the same result

This is not an extra assumption. It is already inside eqs. (10)–(11). Take eq. (10) with the gas at rest. Then $\mathbf F=\rho\,D\mathbf u/Dt=0$, so the $Y_i\mathbf F$ term drops, and neglect $\bar\nu_{ie}$. The drag term $+\frac{n_i}{\bar\nu_i}\bar\nu_{in}\mathbf u_n$ then becomes $-\frac{n_i}{n_n}\boldsymbol\Gamma_i$ because of eq. (11). Moving it to the left-hand side:

$$
\boldsymbol\Gamma_i\left(1+\frac{n_i}{n_n}\right)=\frac{en_i\mathbf E-\nabla p_i}{M\bar\nu_{in}}
\quad\Rightarrow\quad
\boldsymbol\Gamma_i=\frac{en_i\mathbf E-\nabla p_i}{M\bar\nu_{in}\,n_h/n_n}.
$$

The appendix mobility (eq. 50, $\bar\nu_{in}\propto n_n$) is correct as written, because the code keeps $\mathbf u_n$ explicitly in the drag term. The $n_h$ appears only once the neutral counter-flow is substituted in. For a back-of-the-envelope $D_a$, which has no separate $\mathbf u_n$, you must put it in by hand, and that is what $n_h$ does.

### Why it matters only at 20 000 K

| | $n_n$ [m$^{-3}$] | $n_h$ [m$^{-3}$] | $n_h/n_n$ | $D_a$ with $n_n$ | $D_a$ with $n_h$ |
|---|---|---|---|---|---|
| 10 000 K | $7.04\times10^{23}$ | $7.19\times10^{23}$ | 1.02 | $3.04\times10^{-3}$ m$^2$/s | $2.97\times10^{-3}$ m$^2$/s |
| 20 000 K | $3.77\times10^{21}$ | $1.85\times10^{23}$ | 49 | $0.91$ m$^2$/s | $1.85\times10^{-2}$ m$^2$/s |

* **Weakly ionized (10 000 K, 2 %).** The neutrals barely need to move to make room for the ions ($u_n\approx-0.02\,u_i$), so $n_h\approx n_n$ and the choice makes no difference.
* **Nearly fully ionized (20 000 K, 98 %).** With $n_n$ alone, the ions would feel only 2 % of the gas as obstacles, and $D_a\to\infty$ as $n_n\to0$. But for the plasma to shift, the few remaining neutrals must stream through the dense ions at $u_n=-49\,u_i$, and they collide with ions at the rate $n_i\overline Q_{in}\bar c$. The bottleneck is the neutrals pushing through the ions, and that is controlled by $n_i$, which dominates $n_h$. Using $n_n$ would overestimate $D_a$ by 49× and $\sqrt{D_a\tau_{chem}}$ by 7×: $4.8\times10^{-4}$ m instead of $6.9\times10^{-5}$ m.

**In the limit** $n_n\to0$ there is nothing left to diffuse relative to. The "plasma" *is* the whole gas, and a gas at rest has no internal diffusion. The $n_h$ form handles this correctly: it describes neutrals diffusing through ions, so it stays finite. The $n_n$ form wrongly diverges.
