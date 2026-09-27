# Response to questions

## Question 1

> **Update:** the suggested fix below has been applied in `sections/model.tex` (Section 2.1): $\overline Q_{ee}$ now uses $12\sqrt2\,\pi$ and is defined as $\bar\nu_{ee}/(n_e\bar c_e)$, and the Coulomb logarithm is written as $\ln(1.24\times10^7T_e^{3/2}n_e^{-1/2})$. With LTE argon at 1 atm and $\overline Q_{ea}=0.2\times10^{-19}$ m$^2$, $\lambda_e$ becomes $1.7\times10^{-6}$ m at 10 000 K (was $3.2\times10^{-6}$) and $0.6\times10^{-6}$ m at 20 000 K (was $1.1\times10^{-6}$), so the text's "$\lambda_e\sim10^{-6}$ m" is unchanged.

Section 2.1 of the paper gives $\overline Q_{ee}=e^4\ln\Lambda/[32\pi\varepsilon_0^2(k_BT_e)^2]$, citing M&K. This seems to contradict the value obtained from M&K eq. (VIII 2.2b). Why?

## Response

You are right: the two do not agree. The expression in Section 2.1 is **about half** the value that follows from M&K (2.2b). No averaging convention consistent with (2.2b) gives the $32\pi$.

### What (2.2b) implies for $\overline Q_{ee}$

M&K (2.2b) is

$$
\bar\nu_{ee}=n_e\frac{8\sqrt\pi}{3}\Big(\frac{m_e}{k_BT_e}\Big)^{3/2}\Big(\frac{e^2}{4\pi\varepsilon_0m_e}\Big)^2\ln\Lambda .
$$

Section 2.1 writes the mean free path as $\lambda_e=1/(n_e\overline Q_{ee}+n_a\overline Q_{ea})$, which is the same as $\lambda_e=\bar c_e/(\bar\nu_{ee}+\bar\nu_{ea})$ (eq. `eq:mfp`), provided each averaged cross section is defined by

$$
\overline Q_{ee}\equiv\frac{\bar\nu_{ee}}{n_e\,\bar c_e},\qquad \bar c_e=\sqrt{\frac{8k_BT_e}{\pi m_e}} .
$$

Substituting (2.2b):

$$
\overline Q_{ee}
=\frac{8\sqrt\pi}{3}\,\frac{e^4\ln\Lambda}{16\pi^2\varepsilon_0^2}\,
\frac{m_e^{3/2}}{m_e^2(k_BT_e)^{3/2}}\,\sqrt{\frac{\pi m_e}{8k_BT_e}} .
$$

The powers of $m_e$ cancel, $(k_BT_e)$ appears as $(k_BT_e)^{-2}$, and the number is $\dfrac{8\sqrt\pi\cdot\sqrt\pi}{3\cdot16\pi^2\cdot2\sqrt2}=\dfrac{1}{12\sqrt2\,\pi}$. So

$$
\boxed{\;\overline Q_{ee}=\frac{e^4\ln\Lambda}{12\sqrt2\,\pi\,\varepsilon_0^2(k_BT_e)^2}\approx\frac{e^4\ln\Lambda}{16.97\,\pi\,\varepsilon_0^2(k_BT_e)^2}\;}
$$

compared with the paper's $e^4\ln\Lambda/(32\pi\varepsilon_0^2(k_BT_e)^2)$. The ratio is

$$
\frac{\overline Q_{ee}^{\rm paper}}{\overline Q_{ee}^{(2.2b)}}=\frac{12\sqrt2}{32}=0.530 ,
$$

independent of $T_e$, $n_e$ and $\ln\Lambda$. I checked this numerically in two ways. Evaluating (2.2b) directly and dividing by $n_e\bar c_e$ gives the boxed value to machine precision. Separately, averaging the e–e Rutherford cross section $Q^{(1)}=4\pi b_0^2\ln\Lambda$, with $b_0=e^2/(4\pi\varepsilon_0\mu g^2)$ and $\mu=m_e/2$, over the relative-speed Maxwellian with the energy weight of (2.2a) reproduces $\bar\nu_{ee}$ from (2.2b) to within $10^{-4}$.

Per unit $\ln\Lambda$:

| $T_e$ | paper, $32\pi$ | from (2.2b), $12\sqrt2\,\pi$ |
|---|---|---|
| 10 000 K | $4.39\times10^{-18}$ m² | $8.27\times10^{-18}$ m² |
| 20 000 K | $1.10\times10^{-18}$ m² | $2.07\times10^{-18}$ m² |

### Other conventions do not rescue the $32\pi$

- Using the **relative** mean speed of an e–e pair, $\sqrt2\,\bar c_e$, instead of $\bar c_e$ gives $e^4\ln\Lambda/(24\pi\varepsilon_0^2(k_BT_e)^2)$. That is still not $32\pi$, and it would be inconsistent with $\lambda_e=\bar c_e/\nu$ anyway.
- Using $\bar\nu_{ei}$ (M&K 5.1d) with $n_i\to n_e$ instead of $\bar\nu_{ee}$ also gives $24\pi$, because $\bar\nu_{ee}=\sqrt2\,\bar\nu_{ei}$ (for an e–e pair the reduced mass is $m_e/2$ and the relative-speed Maxwellian is $\sqrt2$ wider; $8\sqrt\pi/(4\sqrt{2\pi})=\sqrt2$).

The paper's expression can be written as $\tfrac{\pi}{2}\big(e^2/4\pi\varepsilon_0k_BT_e\big)^2\ln\Lambda$. That is the Coulomb cross section evaluated at one representative thermal energy, not an average over the Maxwellian. It probably comes from a "cross section at thermal energy" estimate (possibly from another chapter of M&K or another textbook), with a numerical prefactor that is not the one implied by (2.2b). Chapter II/VII are not in the PDF we have, so I could not check which formula the citation was meant to point to. What is certain is that it is not consistent with (VIII 2.2b), and (2.2b) is M&K's own formula for $\bar\nu_{ee}$.

### Does it matter for the conclusions?

Only at the factor-of-2 level. With the (2.2b) value, $\overline Q_{ee}$ is 1.89 times larger, so the e–e contribution to $1/\lambda_e$ is 1.89 times larger. Where e–e collisions dominate (the hotter, more ionized side), $\lambda_e$ is up to about 1.9 times **shorter** than the Section 2.1 estimate. $\lambda_e\sim10^{-6}$ m remains correct as an order of magnitude, and the hierarchy $\lambda_{De}\ll\lambda_i\ll\lambda_e\ll\lambda_d\sim\lambda_T$ is unchanged.

### A second issue in the same sentence

The paper says "$\ln\Lambda=1.24\times10^7T_e^{3/2}n_e^{-1/2}$". That expression is $\Lambda$ itself ($\Lambda=12\pi n_e\lambda_{De}^3$), not $\ln\Lambda$. It should read

$$
\ln\Lambda=\ln\!\big(1.24\times10^7\,T_e^{3/2}n_e^{-1/2}\big),\qquad T_e\ \text{in K},\ n_e\ \text{in m}^{-3}.
$$

If the non-logarithmic value were used, $\overline Q_{ee}$ would be off by orders of magnitude. The quoted $\lambda_e\sim10^{-6}$ m suggests the logarithm was actually taken, so this looks like a typo only.

### Suggestion for the paper (not applied)

Replace the Section 2.1 expression with

$$
\overline Q_{ee}=\frac{\bar\nu_{ee}}{n_e\bar c_e}=\frac{e^4\ln\Lambda}{12\sqrt2\,\pi\,\varepsilon_0^2(k_BT_e)^2}
\quad\text{(from M\&K eq. VIII 2.2b)},
$$

fix $\ln\Lambda=\ln(1.24\times10^7T_e^{3/2}n_e^{-1/2})$, and recompute the quoted $\lambda_e$. It will change by less than a factor of 2. This also makes Section 2.1 consistent with eq. (`eq:mfp`) in the length-scales subsection, which uses $\lambda_e=\bar c_e/(\nu_{en}+\nu_{ee})$.

---

## Question 2

Do not use eq. (10), the effective frequency. Keep the ion flux, eq. (11), in the same form as the electron flux, eq. (9), to avoid confusing readers. Also correct Section 2.1 for grammar and flow, and do not use a subsection for each governing equation.

## Response

All three changes are in `sections/model.tex` (and one paragraph of `sections/appendix.tex`). The paper compiles with no errors or undefined references. Equation numbers below are the **new** ones in the compiled `main.pdf`. Everything after the electron flux has moved down by one, because the old eq. (10) is gone.

### 1. The ion flux now has the same form as the electron flux

The effective frequency $\tilde\nu^M_i=(n_h/n_n)\bar\nu^M_{in}+\bar\nu^M_{ie}$ (old eq. 10) and every use of $\tilde\nu^M_i$ have been removed. As for the electrons, the ion flux is obtained by solving the ion balance (4) directly for $\mathbf u_i$. With the total ion momentum-transfer frequency

$$
\bar\nu^M_i=\bar\nu^M_{in}+\bar\nu^M_{ie},
$$

the flux (now eq. 10) is

$$
\boldsymbol\Gamma_i=\frac{en_i}{M\bar\nu^M_i}\mathbf E
-\frac{k_BT_h}{M\bar\nu^M_i}\nabla n_i
-\frac{n_ik_B}{M\bar\nu^M_i}\nabla T_h
+C_i\frac{n_ek_B}{M\bar\nu^M_i}\nabla T_e
-\frac{n_i}{n_h}\frac{\mathbf F}{M\bar\nu^M_i}
+\frac{n_i}{\bar\nu^M_i}\big(\bar\nu^M_{in}\mathbf u_n+\bar\nu^M_{ie}\mathbf u_e\big),
$$

which matches the electron flux (9) term by term:

$$
\boldsymbol\Gamma_e=-\frac{en_e}{m_e\bar\nu^M_e}\mathbf E
-\frac{k_BT_e}{m_e\bar\nu^M_e}\nabla n_e
-(1+C_e)\frac{n_ek_B}{m_e\bar\nu^M_e}\nabla T_e
+\frac{n_e}{\bar\nu^M_e}\big(\bar\nu^M_{en}\mathbf u_n+\bar\nu^M_{ei}\mathbf u_i\big).
$$

**Is this consistent with the solver?** Yes. PlasmeX also writes the ion balance with the plain sum $\bar\nu^M_{in}+\bar\nu^M_{ie}$ (`src/aider.H`: `mu_tot = mu_in/(1+r)`, with $r=\bar\nu^M_{ie}/\bar\nu^M_{in}$). It then substitutes $n_n\mathbf u_n=n_h\mathbf u-n_i\mathbf u_i$ into the drag term, which gives the factor `alpha`. The old $\tilde\nu^M_i$ was exactly that substitution done in the paper. The new form leaves $\mathbf u_n$ in the drag term instead. A new sentence after the fluxes says so: the drag terms couple the two fluxes, and $\mathbf u_n$ is fixed by the frame condition, $n_n\mathbf u_n=n_h\mathbf u-\boldsymbol\Gamma_i$.

**Other places updated to match:**

- **Discussion after the ion flux.** It now says that $\bar\nu^M_i\approx\bar\nu^M_{in}$ in the weakly ionized bulk, and that $\bar\nu^M_i\to\bar\nu^M_{ie}$ as $n_n\to0$. The counter-flow ($n_h/n_n$) explanation was dropped.
- **Classical limit.** The ion mobility is now $e/(M\bar\nu^M_i)$. The extra condition $n_i\ll n_n$ is no longer needed.
- **Ion boundary condition (19).** The drift term is now $\dfrac{en_i}{M\bar\nu^M_i}\mathbf E\cdot\hat{\mathbf n}$.
- **Appendix, ion frequencies.** The paragraph after eq. (`eq:nu-ie`) now defines $\bar\nu^M_i=\bar\nu^M_{in}+\bar\nu^M_{ie}$ instead of referring to the deleted equation.

### 2. No subsections for individual governing equations

The seven `\subsubsection` headings (species continuity, momentum balances, fluxes, electrostatic field, electron energy, heavy-particle energy, heavy-gas momentum) have been removed. Section 2.1 is now running text. Each equation is introduced by a lead-in sentence, for example "The species fluxes follow from the species momentum equations", "The electric field is obtained from Poisson's equation", or "Finally, the heavy-gas velocity is advanced by…".

The notation paragraph (partial pressures, $n_h$, $\rho$, $\mathbf u_s$, $\boldsymbol\Gamma_s$, $\mathbf u$) now comes **before** the continuity equations, so every symbol is defined before it is used. It previously ended with "These frequencies and the remaining transport coefficients…" before any frequency had been introduced. That now reads "All transport coefficients introduced below…".

### 3. Grammar and flow in Section 2.1

Typos fixed: "considerin th", "transiet", "therfore", "inplies", "apperars", "a oppposite", "are therefore are", "shows the same", "basis as well a necessity". Also fixed: missing spaces ("temperature.Use", "\(E_z\)from") and "2e-4 m" (now $2\times10^{-4}$ m). The broken $\overline Q_{ee}$ sentence ("for which the averaged … cross section, defined as …") now reads "for which the averaged momentum-transfer cross section is $\overline Q_{ee}=\dots$". $k_B$ is now defined where it is first used.

Beyond grammar, a few corrections were needed for consistency. Please check these:

| where | before | now | why |
|---|---|---|---|
| mass fractions | $Y_e+Y_i+Y_n=1$ | $Y_i+Y_n=1$, $Y_e\sim m_e/M$ | $\rho=Mn_h$ excludes electrons, so $Y_i+Y_n=1$ exactly. The sum of (3)–(5) is $0=0$ to within $m_e/M$. |
| electron flux | $Y_e\mathbf F$ silently absent | "the negligible term $Y_e\mathbf F$ has been dropped" | (3) contains it and (9) does not |
| neutral flux | "multiplying the frame condition by $n_n$" | "by $n_h$" | $n_h\mathbf u=n_i\mathbf u_i+n_n\mathbf u_n$ |
| heavy mass conservation | "equation (heavy-momentum) implies $\partial_tn_h+\nabla\cdot(n_h\mathbf u)=0$" | follows from adding the ion and neutral continuity equations with (11) | the momentum equation does not imply mass conservation |
| $\overline Q_{ia}$ fit | $2[7\times10^{-10}-0.6\times10^{-10}\ln(2k_BT_h)]$ | $2[\dots\ln(2k_BT_h/e)]^2~\mathrm m^2$ | now matches the appendix |
| $\lambda_d$ | $E_z$ undefined | "$E_z$ is the axial electric field in the arc column" | **please confirm this is what you meant** |
| $\lambda_T$ | stated without reason | one clause: each elastic collision transfers only a fraction $\sim2m_e/M$ of the electron energy | explains why $\lambda_T\gg\lambda_e$ |
| momentum eq. | "$\mathbf F_s$ are all the forces", $\mathbf R_s$ undefined | $\mathbf F_s$ = electric + pressure force, $\mathbf R_s$ = momentum from collisions | both symbols are now defined |

### Two things I did **not** change, which you should look at

**(a) $\lambda_i\sim10^{-7}$ m does not follow from LTE densities.** With the fit above:

| $T$ | $n_a$ (LTE, 1 atm) | $\overline Q_{ia}$ | $\lambda_i=1/(n_a\overline Q_{ia})$ |
|---|---|---|---|
| 10 000 K | $7.0\times10^{23}$ m$^{-3}$ | $8.9\times10^{-19}$ m² | $1.6\times10^{-6}$ m |
| 20 000 K | $3.9\times10^{21}$ m$^{-3}$ | $7.8\times10^{-19}$ m² | $3.3\times10^{-4}$ m |
| 3 000 K (cold gas at wall) | $2.4\times10^{24}$ m$^{-3}$ | $1.1\times10^{-18}$ m² | $3.7\times10^{-7}$ m |

At LTE conditions, $\lambda_i\gtrsim\lambda_e$, which contradicts the stated hierarchy $\lambda_{De}\ll\lambda_i\ll\lambda_e$. The value $\sim10^{-7}$ m holds only for the cold, dense heavy gas next to the electrode. Either the text should say that $\lambda_i$ is evaluated at the near-wall heavy-gas temperature, or the hierarchy needs revising.

**(b) The ion boundary condition and the code now differ slightly.** *(Superseded by Question 3, part 3: the boundary condition differs from the code in more than the mobility.)* Eq. (19) now uses $\bar\nu^M_i=\bar\nu^M_{in}+\bar\nu^M_{ie}$. The solver's wall mobility (`mob(ARP_ID)`) still includes the frame substitution, which amounts to $e/[M((n_h/n_n)\bar\nu^M_{in}+\bar\nu^M_{ie})]$. The difference is the factor $n_h/n_n=1+n_i/n_n$, which is close to 1 at the wall, where $n_i\ll n_n$. If you want the paper and the code to be identical, either note this in the text or change the code's boundary drift to use `mu_tot`.

Backups of the previous `model.tex` and `appendix.tex` are in the session scratchpad.

---

## Question 3

Even though a different mobility is used, the governing equations solved in the paper and in the code are the same, right? Explain.

## Response

**In the interior of the domain: yes, exactly the same.** The paper's new ion flux, eq. (10), and the code's ion flux are the same equation. They are just arranged differently. **At the electrode: no.** The paper's ion boundary condition, eq. (19), does not match the code. The main reason is not the mobility. The mismatch is older than this edit (details in part 3, which replaces item (b) of Question 2).

### 1. Why the interior equations are identical

Both start from the same ion momentum balance, eq. (4). Collect every force that is not friction into

$$
\mathbf G_i=en_i\mathbf E-\nabla p_i+C_in_ek_B\nabla T_e-Y_i\mathbf F .
$$

Then (4) reads

$$
0=\mathbf G_i-Mn_i\bar\nu^M_{in}(\mathbf u_i-\mathbf u_n)-Mn_i\bar\nu^M_{ie}(\mathbf u_i-\mathbf u_e).
$$

**Paper form (eq. 10).** Move the $\mathbf u_i$ terms to the left, with $\bar\nu^M_i=\bar\nu^M_{in}+\bar\nu^M_{ie}$:

$$
\mathbf u_i=\frac{\mathbf G_i}{Mn_i\bar\nu^M_i}+\frac{\bar\nu^M_{in}\mathbf u_n+\bar\nu^M_{ie}\mathbf u_e}{\bar\nu^M_i}. \tag{A}
$$

Here $\mathbf u_n$ is still on the right. It is not independent, though. The frame condition (8) fixes it in terms of $\mathbf u_i$:

$$
\mathbf u_n=\frac{n_h\mathbf u-n_i\mathbf u_i}{n_n}. \tag{B}
$$

**Code form.** The code substitutes (B) into (A) before solving. The $\mathbf u_i$ in (B) moves to the left:

$$
\mathbf u_i\Big(\bar\nu^M_i+\frac{n_i}{n_n}\bar\nu^M_{in}\Big)
=\frac{\mathbf G_i}{Mn_i}+\frac{n_h}{n_n}\bar\nu^M_{in}\mathbf u+\bar\nu^M_{ie}\mathbf u_e .
$$

The bracket on the left is

$$
\bar\nu^M_{in}+\bar\nu^M_{ie}+\frac{n_i}{n_n}\bar\nu^M_{in}
=\frac{n_h}{n_n}\bar\nu^M_{in}+\bar\nu^M_{ie}\equiv\tilde\nu^M_i ,
$$

which is exactly the old effective frequency. Since $\tilde\nu^M_i\mathbf u=\big(\tfrac{n_h}{n_n}\bar\nu^M_{in}+\bar\nu^M_{ie}\big)\mathbf u$, this becomes

$$
\mathbf u_i=\mathbf u+\frac{\mathbf G_i}{Mn_i\tilde\nu^M_i}+\frac{\bar\nu^M_{ie}}{\tilde\nu^M_i}(\mathbf u_e-\mathbf u). \tag{C}
$$

(C) is the old eq. (11), and it is what PlasmeX solves. In `src/aider.H`:

- `mu_tot = mu_in/(1+r)` is $e/(M\bar\nu^M_i)$, with $r=\bar\nu^M_{ie}/\bar\nu^M_{in}$;
- `alpha = 1/(1 + w_in n_i/n_n)`, with $w_{in}=\bar\nu^M_{in}/\bar\nu^M_i$, is the factor from substituting (B);
- `mob(ARP_ID) = alpha*mu_tot` $=\dfrac{e}{M(\bar\nu^M_i+\tfrac{n_i}{n_n}\bar\nu^M_{in})}=\dfrac{e}{M\tilde\nu^M_i}$.

The header comment in `aider.H` spells out the same step: "substituting the frame condition into the ION line makes its $u_n$ term algebraic".

So (A)+(B) and (C) are the same pair of equations. (C) is (A) with (B) already substituted. For a given state, they give the same $\boldsymbol\Gamma_i$ and the same $\boldsymbol\Gamma_n=n_h\mathbf u-\boldsymbol\Gamma_i$. The electron flux (9), the continuity, Poisson, energy and heavy-gas momentum equations are unchanged, so the whole interior system in the paper is the one the code solves.

### 2. So what does "the mobility" mean in each form?

The two coefficients on $\mathbf E$ answer different questions. Both are correct.

| | coefficient on $\mathbf E$ | what it is |
|---|---|---|
| paper, eq. (10) | $\dfrac{e}{M\bar\nu^M_i}$ | ion drift **relative to the other species**, with $\mathbf u_n$ and $\mathbf u_e$ held fixed. The neutral part, $e/(M\bar\nu^M_{in})$, is the measured charge-exchange mobility $\mu_{in}$. |
| code, (C) | $\dfrac{e}{M\tilde\nu^M_i}$ | ion drift **relative to the heavy-gas frame** $\mathbf u$, after the neutrals' counter-flow (B) has been accounted for |

When ions drift through the neutrals, conservation of heavy-particle mass forces the neutrals to move the opposite way. This raises the ion–neutral relative velocity, so the net ion drift for a given field is smaller than $e/(M\bar\nu^M_i)$ by the factor $\bar\nu^M_i/\tilde\nu^M_i$. In the paper form, that reduction is not in the coefficient; it comes back through the $\bar\nu^M_{in}\mathbf u_n$ drag term. In the code form, it is built into the coefficient. The size of the effect is

$$
\frac{\tilde\nu^M_i}{\bar\nu^M_i}=1+\frac{\bar\nu^M_{in}}{\bar\nu^M_i}\,\frac{n_i}{n_n}\approx1+\frac{n_i}{n_n},
$$

which is close to 1 wherever the gas is weakly ionized. It only matters towards the strongly ionized arc side.

The two forms are also treated slightly differently numerically. The code handles the $\mathbf u_n$ back-reaction implicitly, through `alpha`, and lags only $\mathbf u_e$ (and $\mathbf u_i$ in the electron drag) to the previous step. That is a time-discretization choice, not a different equation.

### 3. Where paper and code do differ: the ion wall flux

The paper's eq. (19) says

$$
\boldsymbol\Gamma_i\cdot\hat{\mathbf n}=\tfrac14n_i\bar c_i+\max\!\Big(\frac{en_i}{M\bar\nu^M_i}\mathbf E\cdot\hat{\mathbf n},\,0\Big).
$$

The code (`exe/bc.H`, `w.i_wall`) does

$$
\boldsymbol\Gamma_i\cdot\hat{\mathbf n}=
\begin{cases}
\boldsymbol\Gamma_i^{\rm expl}\cdot\hat{\mathbf n}, & \text{if }\boldsymbol\Gamma_i^{\rm expl}\cdot\hat{\mathbf n}>0\ \text{(towards the wall)},\\[4pt]
\tfrac14n_i\bar c_i, & \text{otherwise,}
\end{cases}
\qquad \bar c_i=\sqrt{\frac{8k_BT_{\rm wall}}{\pi M}},
$$

where $\boldsymbol\Gamma_i^{\rm expl}$ is the **whole** explicit ion flux on the wall face. That is $n_i\big[\tfrac{e}{M\tilde\nu^M_i}\mathbf E+\text{thermal, thermal-diffusion, drag and mixture-force velocities}\big]$.

There are three differences, in decreasing order of importance:

1. **Either/or versus sum.** The code uses the drift flux *or* the thermal flux. The paper adds them. When the field drives ions into the wall, the paper's value is larger by $\tfrac14n_i\bar c_i$. At a cold wall ($T_{\rm wall}\sim10^3$ K, $\bar c_i\approx700$ m/s) this is small compared with the sheath drift, but it is a different formula.
2. **Field drift only versus total explicit flux.** The paper's drift term contains only $\mathbf E$. The code's includes the thermal-diffusion, drag and mixture-force contributions, and the "towards the wall" test is applied to that total.
3. **Mobility.** The paper now has $e/(M\bar\nu^M_i)$ and the code has $e/(M\tilde\nu^M_i)$. At the wall this differs by the factor $1+n_i/n_n\approx1$, so it is the least important of the three. It is also the only one my edit introduced; points 1 and 2 were already there.

### 4. Recommendation

The interior equations need no change. To make eq. (19) match the code, write it as

$$
\boldsymbol\Gamma_i\cdot\hat{\mathbf n}=
\begin{cases}
\boldsymbol\Gamma_i\cdot\hat{\mathbf n}\ \text{from eq. (10)}, & \text{if directed into the wall},\\
\tfrac14n_i\bar c_i, & \text{otherwise,}
\end{cases}
$$

that is, "ions reach the wall at the flux given by eq. (10) when that flux is directed into the wall, and at the thermal-effusion flux otherwise". Because it points to eq. (10) itself, this also removes the question of which mobility to quote. Evaluating eq. (10) with $\mathbf u_n$ from (B) gives exactly the code's $\tilde\nu^M_i$ coefficient. I have **not** made this change; tell me if you want it.

---

## Question 4

Fix the ion boundary condition in the paper.

## Response

Done. Eq. (19) in `sections/model.tex` and the electrode-BC paragraph in `sections/methods.tex` now describe what `exe/bc.H` (`w.i_wall`) does. The paper compiles with no errors or undefined references.

### The new boundary condition

Let $\boldsymbol\Gamma_i^{d}$ be the ion flux (10) **without its density-diffusion term**:

$$
\boldsymbol\Gamma_i^{d}=\boldsymbol\Gamma_i+\frac{k_BT_h}{M\bar\nu^M_i}\nabla n_i .
$$

It contains the field drift and the thermal, thermal-diffusion, mixture-force and drag contributions. Then

$$
\boldsymbol\Gamma_i\cdot\hat{\mathbf n}=
\begin{cases}
\boldsymbol\Gamma_i^{d}\cdot\hat{\mathbf n}, & \boldsymbol\Gamma_i^{d}\cdot\hat{\mathbf n}>0,\\[4pt]
\tfrac14n_i\bar c_i, & \text{otherwise},
\end{cases}
\qquad
\bar c_i=\sqrt{\frac{8k_BT_{\rm wall}}{\pi M}} .
$$

The text after the equation explains the two cases. In the first (sheath field driving ions into the electrode), the ion density gradient at the wall is zero. In the second, the thermal-effusion flux still absorbs ions where transport alone does not carry them into the wall.

### Why "without the density-diffusion term"

This is a detail I missed in Question 3. The code's "drift" is `expl_flx(ARP_ID)`, the **explicit** part of the ion flux: $n_i(\mu_i\mathbf E+\text{extra velocities})$. The density-diffusion term $-D_i\nabla n_i$ is treated implicitly and is not part of it. `species_bc` then sets the implicit diffusive flux at the wall to `i_wall - i_drift`:

- **drift into the wall:** `i_wall = i_drift`, so the diffusive flux is zero. That is the zero-gradient statement.
- **otherwise:** the diffusive flux is whatever makes the total equal $\tfrac14n_i\bar c_i$.

So the test and the value both use $\boldsymbol\Gamma_i^{d}$, not the full $\boldsymbol\Gamma_i$. The recommendation at the end of Question 3 said "the flux given by eq. (10)"; the paper now says the precise version.

The mobility question from Question 3 also disappears. $\boldsymbol\Gamma_i^{d}$ is defined from eq. (10), and evaluating it with $\mathbf u_n$ from the frame condition gives exactly the code's $e/(M\tilde\nu^M_i)$ coefficient.

### Other changes in the same paragraphs

| file | removed | why |
|---|---|---|
| `model.tex` | "This treatment is consistent with the field-aware treatment of the counter-diffusing sheath-edge current used elsewhere [Benilov2016, LiangTrelles]." | That comparison was for the sum form $\tfrac14n_i\bar c_i+\max(\mu_in_iE,0)$. The either/or form is a different closure, and I could not verify that it matches those papers. |
| `methods.tex` | "Hagelaar-type drift-aware closure, $\Gamma_i=\tfrac14n_i\bar c_i+\max(\mu_in_iE,0)$ … structurally the same closure used for … [Benilov2016] (their $J_{CD}$, eq. A.7) and … [LiangTrelles] (their $J_{bd}$, eq. 11) …" | Same reason. It is replaced by a description of eq. (19) and of how the wall value is imposed on the implicit diffusion flux. |

Both citations are still used elsewhere in the paper, so the bibliography is unchanged. If you want a literature comparison back, it should be rechecked against the either/or form.

The rest of the Methods paragraph is unchanged. That includes the two open items: the field-independent electron wall flux and the missing emission term.

Backups of the previous `model.tex` and `methods.tex` are in the session scratchpad.

---

## Question 5

Using the PlasmeX code now in `./PlasmeX`, update the boundary-condition section for all governing equations. For the cathode region, electrode emission is on; for the anode region, it is off.

## Response

Done. Section 2.2 (`sections/model.tex`) is rewritten, and the electrode-BC paragraph in `sections/methods.tex` is updated to match. The paper compiles with no errors, undefined references or overfull boxes. The pasted `./PlasmeX` is identical to `../PlasmeX`. The sources used are `exe/bc.H` (`extern_bc`, `potential_bc`, `neutral_bc`, `electrode_wall`, `electrode_wall_temp`, `species_bc`), `exe/prob.H` (`arc_state`), `src/solve_neutrals.cpp` (heavy-gas boundary faces), `src/aider.H` (`UpdateCurrentDrivenBC`, `UpdateWallTemperature`), and the two decks `exe/anode_inputs` and `exe/cathode_inputs`.

### What the section now contains, equation by equation

Equation numbers are the new ones in `main.pdf`. $\hat{\mathbf n}$ points into the wall.

| governing equation | electrode, $x=0$ | bulk edge, $x=L$ |
|---|---|---|
| electron continuity | $\boldsymbol\Gamma_e\cdot\hat{\mathbf n}=\tfrac14n_e\bar c_e-\Gamma_{\rm em}$ (18) | $n_e=n_e^{\rm bulk}$ |
| ion continuity | drift-or-thermal switch (21) | $n_i=n_i^{\rm bulk}$ |
| neutral continuity | $\boldsymbol\Gamma_n\cdot\hat{\mathbf n}=-\boldsymbol\Gamma_i\cdot\hat{\mathbf n}$ (22) | Dalton at $p_0$ (25) |
| Poisson | $\phi=0$ | total current $=J_{\rm app}$ (26) |
| electron energy | $\mathbf h_e\cdot\hat{\mathbf n}=2k_BT_e\,\tfrac14n_e\bar c_e-2k_BT_{\rm wall}\Gamma_{\rm em}$ (20) | $T_e=T_{\rm bulk}$ |
| heavy-particle energy | $T_h=T_{\rm wall}$; fixed (anode) or heat balance (23)–(24) (cathode) | $T_h=T_{\rm bulk}$ |
| heavy-gas momentum | $\mathbf u\cdot\hat{\mathbf n}=0$, wall pressure from the interior | open: bulk state (pressure $p_0$) outside, $\partial\mathbf u/\partial n=0$ |

**Thermionic emission (19)**, Richardson–Dushman with Schottky lowering (`electrode_wall`, `Jem`):

$$
\Gamma_{\rm em}=\frac{A_RT_{\rm wall}^2}{e}\exp\!\Big(-\frac{W-\Delta W}{k_BT_{\rm wall}}\Big),\qquad
\Delta W=\sqrt{\frac{e^3E_w}{4\pi\varepsilon_0}},\qquad E_w=\max(\mathbf E\cdot\hat{\mathbf n},0).
$$

It is applied only when `prob.electrode_emission = 1`: on for the cathode, off for the anode ($\Gamma_{\rm em}=0$). The paper says so explicitly.

**Cathode heat balance (23)–(24)** (`electrode_wall_temp`, applied when `prob.wall_temp_balance = 1`, cathode only):

$$
q_w=\frac{k_c}{L_c}(T_{\rm wall}-T_0),
$$

$$
q_w=(2k_BT_e+W)\tfrac14n_e\bar c_e-(2k_BT_{\rm wall}+W)\Gamma_{\rm em}+(\varepsilon_{\rm ion}-W)\,\boldsymbol\Gamma_i\cdot\hat{\mathbf n}-\kappa_h\nabla T_h\cdot\hat{\mathbf n},
$$

which is term by term `q_wall = eeng_maxw + e_maxw*W + jem_e_app + jem_app*W + q_ion + cond`. The two emission terms are negative (cooling). The code solves it for $T_{\rm wall}$ by bisection at the start of every step. The text explains why the root is unique: emission cooling rises steeply with $T_{\rm wall}$, while the conduction term rises linearly.

*(Table 1 was later removed from the paper; see Question 6. The values below are kept here for reference only.)*

**New table 1** gathers what differs between the two cases, taken from the decks:

| | anode | cathode |
|---|---|---|
| $T_{\rm bulk}$ | 10 000 K | 20 000 K |
| $n_e^{\rm bulk}=n_i^{\rm bulk}$ | $1.4\times10^{22}$ m$^{-3}$ | $1.8\times10^{23}$ m$^{-3}$ |
| $\lvert J_{\rm app}\rvert$ | $1\times10^7$ A/m² | $7.8\times10^7$ A/m² |
| $p_0$ | 101 325 Pa | 101 325 Pa |
| emission | off | on, $W=4.55$ eV (W), $A_R=6\times10^5$ A m$^{-2}$ K$^{-2}$ |
| $T_{\rm wall}$ | 1000 K fixed | heat balance: $k_c=100$ W/m/K, $L_c=1$ cm, $T_0=1000$ K |

### What changed relative to the previous text

- **Emission.** Removed "the surface is taken to be absorbing and non-emitting" and "Thermionic and secondary electron emission are not included … under-predicted". Thermionic emission is now modelled for the cathode. Only secondary emission by ion impact is still missing; there is no such term in the code.
- **Electron energy.** Now a numbered equation (20) that includes the energy brought in by emitted electrons, $2k_BT_{\rm wall}$ each (`Jem_E = Jem*2 k_B T_w`).
- **Wall temperature.** Previously a single "$T_{\rm wall}$". It is now fixed for the anode and from the heat balance for the cathode.
- **Bulk edge.** It is now explicit that $T_e=T_h=T_{\rm bulk}$ (`arc_state`: `gastemp = arc_Te`). The heavy-gas edge is open at pressure $p_0$ (`solve_neutrals.cpp`: "the ARC face is open … its pressure is prob.pres").
- **Methods paragraph.** Now mentions emission for the cathode. It explains that each wall flux, including the energy flux, is imposed as prescribed-minus-explicit on the implicit diffusion flux, and that $T_{\rm wall}$ is re-solved every step. Open items: the field-independent electron wall flux, and no secondary emission.

### Things to check

**(a) The anode deck does not solve the heavy-gas momentum equation.** `anode_inputs` has `plasmex.neutral_closure = 0` and `cathode_inputs` has `= 2`. Under closure 0 there is no $\mathbf u$ ($\mathbf u=0$), and $n_n$ comes from Dalton's law with an electrostatic-pressure correction at every point:

$$
p=p_0+\tfrac12\varepsilon_0\big(E^2-E_{\rm arc}^2\big),\qquad
n_n=\frac{p-n_ek_BT_e}{k_BT_h}-n_i ,
$$

(`solve_neutrals.cpp`, "CLOSURE 0"). Section 2.1 presents only the heavy-gas momentum equation (17), and the new BC section gives its conditions without saying they apply only to the cathode run. So the paper does not describe what the anode run actually solves. Two options:

1. Add one sentence plus the equation above to Section 2.1: the anode case uses the inertialess limit of (17), with $\mathbf u=0$ and the Dalton closure, because at 0.6% ionization the acoustic relaxation is fast compared with anything of interest. Then mark the heavy-gas rows of the BC section "cathode only".
2. Rerun the anode with `neutral_closure = 2`, so both cases use the equations as written.

I have not changed this. It is a decision about the model, not the BCs.

**(b) Anode current sign.** The anode deck has `prob.arc_Japp = -1e7` and the cathode has `+7.8e7`. The code's convention is a coordinate-oriented current with the arc at $x$-lo. The table lists magnitudes, $|J_{\rm app}|$, because the paper puts the electrode at $x=0$ (the opposite end from the code). If you want signed values in the paper, the sign convention needs a sentence.

**(c) Ion kinetic energy in $q_w$.** The code's ion heat term is $(\varepsilon_{\rm ion}-W)\Gamma_i$. It does not include the kinetic energy ions gain crossing the sheath, which in many cathode models is the dominant heating term. In this model, that energy is dissipated into the heavy gas through $Q_{in}$ and reaches the wall through the conduction term. That is consistent as long as the drift-diffusion ions are collisional right up to the wall. You may want one sentence saying so, since reviewers of cathode papers look for this term.

Backups of the previous `model.tex` and `methods.tex` are in the session scratchpad.

---

## Question 6

Do not give values for the boundary conditions in Table 1, because I may run at different values.

## Response

Done. **Table 1 has been removed** from `sections/model.tex`, together with every numerical boundary value in Section 2.2:

- the table itself (bulk temperature and density, $|J_{\rm app}|$, $p_0$, $W$, $A_R$, anode $T_{\rm wall}$);
- the heat-balance constants after eq. (23): "$k_c=100$ W m$^{-1}$ K$^{-1}$ … $L_c=1$ cm … $T_0=1000$ K" now reads "where $k_c$ is the thermal conductivity of the cathode, $L_c$ its length, and $T_0$ the temperature of its cooled base";
- the two references to the table.

With the values gone, the table only had two qualitative rows left (emission on/off, wall temperature fixed/heat balance). The opening paragraph of Section 2.2 already states both, so the table was deleted rather than kept with symbols only. That sentence now reads: "The anode and cathode cases share the same set of conditions and differ in two respects: electron emission from the electrode is switched on for the cathode and off for the anode, and the cathode temperature is obtained from a surface heat balance, whereas the anode temperature is prescribed."

The paper compiles with no errors, undefined references or overfull boxes. When you fix the run parameters, they would fit naturally in the results section, next to each case.

---

## Question 7

> **Superseded** by `literature_review.md`. This answer relied on the old, now-deleted `todo.txt`, and wrongly treated the heavy-gas momentum equation as a physics novelty. It is only a transient-relaxation device: the steady state has $\mathbf u=0$.

Almeida also solves a similar problem. What novelty should I propose in this paper?

## Response

**Short answer:** do not claim "a unified model that resolves the sheath without dividing the layer into sub-layers". That is exactly Almeida, Benilov and Naidis (2008), and the paper already credits them for it. The novelty has to come from what PlasmeX does that Almeida 2008 does **not**. Based on a side-by-side reading of Almeida 2008 and the model as now written in Section 2 and implemented in PlasmeX, there are five real differences. Three are strong enough to be headline claims, provided the results back them up.

### 1. What Almeida 2008 actually does (the baseline to beat)

From the paper itself (`validation_fig/Almeida_2008_...pdf`):

| aspect | Almeida 2008 |
|---|---|
| electrode | **cathode only** (Hg at 100 bar, Ar at 1 atm); the anode is in a separate paper, Almeida 2009, for very-high-pressure Hg |
| time | **steady** boundary-value problem, solved by Newton + Petukhov's 4th-order method on a variable grid |
| control parameter | current density $j_c$, prescribed |
| heavy-gas motion | **none**: "convective effects … are neglected"; mechanical equilibrium $-\nabla p+e(n_i-n_e)\mathbf E=0$ (their eq. 9). They note that ion/atom inertia in the ionization layer "is not described by the diffusion equations used in this work". |
| cathode temperature | **given parameter** $T_w$. Their conclusion says the model must still be "combined with the model of nonlinear surface heating". |
| ion wall BC | $n_i=0$ (Knudsen-order argument) |
| outer BC | local LTE with a radiation-dominated energy balance, solved algebraically at $x=L$ |
| transport | Stefan–Maxwell with Chapman–Enskog coefficients, interpolated between the weakly and strongly ionized limits |
| geometry | 1D spherically symmetric, $r_c=1$ mm |
| outputs | distributions across the layer, current–voltage characteristics, $q_c(T_w)$, and a comparison with spectroscopy |

### 2. Where PlasmeX genuinely differs

Ranked by how defensible and how interesting each difference is:

**(A) Heavy-gas momentum inside the electrode layer: relaxing Almeida's mechanical-equilibrium assumption. Strongest, most "physics" claim.**
Almeida imposes $\nabla p=e(n_i-n_e)\mathbf E$ everywhere. Your model advances $\rho D\mathbf u/Dt=e(n_i-n_e)\mathbf E-\nabla p$ (eq. 17). Your species balances carry the consistent mixture-force shares $Y_s\mathbf F$, so they stay exact in an accelerating frame. That is the piece that makes this more than bolting on a Navier–Stokes solver. The physical hook: in the cathode sheath the electrostatic pressure $\tfrac12\varepsilon_0E^2$ reaches $\sim10^4$ Pa (code comment in `solve_neutrals.cpp`), about 10% of atmospheric. The question "does the heavy gas actually reach mechanical equilibrium, and on what time scale; what flow and density depletion does the sheath force drive?" is one Almeida explicitly left open.
*Needed to claim it:* a cathode result that shows $\mathbf u$, the $n_h$ perturbation and the pressure field. Ideally also a direct comparison: the same case with closure 0 (Almeida's assumption) against closure 2, with the change in fall voltage and heat flux. **The anode deck currently runs closure 0** (Question 5, item a). For this claim to cover both electrodes, run the anode with closure 2 too, or state clearly that it is a cathode result.

**(B) Time-dependent, current-driven (galvanostatic) formulation. Strong, and unique in this class of models.**
Almeida solves a steady problem at fixed $j_c$. PlasmeX integrates in time, with the total current (conduction + displacement) held at $J_{\rm app}$ (eq. 26). This gives you things a steady solver cannot:
- the transient formation of the sheath and the relaxation hierarchy (dielectric relaxation, ionization, electron energy, heavy-gas acoustic, cathode thermal);
- the response to current changes (ramps, steps, pulsed or AC arcs, arc ignition), which matter for welding and for pulsed and switching arcs;
- field reversals and the sign of the anode fall *emerging* rather than being assumed. The intro already makes this argument; this formulation is what makes it true.

*Needed:* at least one genuinely transient result, for example a current step, with the time scales identified. A steady state reached by time-marching alone does not demonstrate this; it only shows the method.

**(C) One formulation for both electrodes of an atmospheric-pressure argon arc. Strong, if the literature check holds.**
Your cathode and anode cases use identical equations and differ only in two boundary conditions (emission on/off, wall temperature balanced/fixed). Almeida 2008 is cathode-only. Almeida 2009 treats the anode, but for ~100 bar Hg, not atmospheric argon. This gives a clean narrative: *the same physics produces a positive cathode fall with strong ion heating at one electrode and a negative anode fall at the other.* The anode result can be validated against three independent sources already in your bibliography or todo list: Dinulescu & Pfender 1980 (−1.5 to −2.7 V), Sanders & Pfender 1984 (measured −1.4 to −2.1 V), and Liang & Trelles. That is the most citable validation you can have (`todo.txt`, item 5).
**Check before claiming "first":** I believe Baeva, Benilov, Almeida and Uhrlandt published 2D unified non-equilibrium models of atmospheric argon arcs that include both electrode layers, around 2016 in J. Phys. D. They are not in `references.bib`. Search Baeva/Benilov 2016–2020 and cite them. Frame (C) as "a single time-dependent 1D formulation applied to both electrodes, with heavy-gas dynamics", not "first model of both electrodes".

**(D) Self-consistent cathode temperature from a surface heat balance. Moderate.**
Almeida takes $T_w$ as an input and says coupling to surface heating is future work. You solve $q_w(T_w)=k_c(T_w-T_0)/L_c$ every step (eqs. 23–24), with the heat-flux partition explicit: electron heating, emission cooling, ion recombination, gas conduction. The 1D rod model is simple, but it closes the loop Almeida left open, and the partition of $q_w$ by channel is a useful output in itself.
*Caveat:* the ion kinetic energy from the sheath enters $q_w$ only indirectly, through $Q_{in}$ and gas conduction (Question 5, item c). Reviewers will ask about this, so show the partition and argue it.

**(E) An electrode-region model designed as a boundary module for LTE arc codes. Practical, and it matches your intro.**
Your outer boundary takes its state from an LTE bulk solution. The natural deliverable is maps of the fall voltage $U$, the heat flux $q_w$ (and its partition), and the ion flux and impact energy, as functions of $(J_{\rm app},T_{\rm bulk},T_w)$ for each electrode. An LTE CFD code can then use these as electrode boundary conditions. Benilov 2016 and Liang & Trelles do "sheath coupling" with simplified sheath models; you would supply the same interface from a unified, collisional-sheath model. This is the claim that justifies the paper's framing ("LTE bulk + resolved electrode region").
*Needed:* a parametric sweep, even a coarse one, plus a table or correlation.

**Differences that are real but not headline material:**
- *Ion wall condition:* flux-type (drift or thermal) instead of $n_i=0$. A sensitivity comparison of the two is a nice short subsection, not a main claim.
- *Transport closure:* M&K-averaged LXCat cross sections joined to Spitzer–Härm, with Onsager-consistent thermal diffusion and field-heated (Wannier) ion mobility. Almeida's Stefan–Maxwell/Chapman–Enskog closure is similar in spirit, so present yours as careful, not new.
- *Numerics:* AMR to sub-Debye resolution with semi-implicit time integration, and an r–z capable code. This is an enabler; claim it only as "extensible to multi-D" unless you show a 2D result.

### 3. Claims to avoid, or to fix first

| claim | problem |
|---|---|
| "first unified model of the electrode region" | Almeida 2008 (and 2009 for the anode) |
| "first to resolve the sheath without a separate sheath model" | same |
| **ion impact energy** (it is in the **title**) | `todo.txt` P0 item 3: the Wannier estimate $\tfrac32kT_{\rm eff}=4.92$ eV currently **exceeds** the whole 1.70 V sheath fall. Until that is fixed, either drop "ion impact energy" from the title or report it only with the collisionless bound $e\Delta\phi$ alongside. |
| any quantitative fall voltage | `todo.txt` P0 items 1–2: the dielectric-relaxation time-step fix and a proper $\Delta x$/$\Delta t$ convergence study come first |

### 4. Suggested positioning

**Recommended headline:** (A) + (B), applied through (C), and delivered as (E). In one sentence:

> *A time-dependent, current-driven unified model of the near-electrode region that includes heavy-gas dynamics, applied with a single set of equations to both the cathode and the anode of an atmospheric-pressure argon arc, and used to provide the electrode boundary data (fall voltage, heat-flux partition, species fluxes) needed by LTE arc models.*

**Draft contribution paragraph** for the end of the introduction (edit once the results exist):

> The present work extends the unified approach of Almeida *et al.* in four respects. First, the heavy gas is not assumed to be in mechanical equilibrium: a heavy-gas momentum equation is solved together with species momentum balances that remain consistent in the accelerating frame, which allows the flow and density perturbation driven by the space-charge force in the sheath to be quantified. Second, the equations are integrated in time with the arc current imposed through the total (conduction plus displacement) current, so that the formation of the sheath, the hierarchy of relaxation times, and the response of the electrode layer to changes in current are obtained directly, and the polarity of the electrode fall emerges from the solution. Third, the same formulation is applied to both the cathode and the anode of an atmospheric-pressure argon arc, the two cases differing only in the electron emission and thermal conditions at the electrode surface; at the cathode, the surface temperature is obtained self-consistently from the surface heat balance. Fourth, the results are cast as electrode boundary data—fall voltage, heat flux and its partition, and species fluxes as functions of current density and bulk-plasma state—for use in LTE models of the arc column.

### 5. Minimum results to support this

1. Fix the time step (`todo.txt` P0 #1) and do the $\Delta x$/$\Delta t$ convergence study (#2). Nothing quantitative is publishable before this.
2. Cathode, closure 0 vs closure 2 at the same conditions → claim (A).
3. A current step or ramp at one electrode, with the relaxation times identified → claim (B).
4. Anode with closure 2 (or a stated reason for 0); the fall voltage and its sign compared with Dinulescu & Pfender, Sanders & Pfender, and Liang & Trelles → claim (C).
5. Cathode $T_w$ and the $q_w$ partition at 2–3 current densities, and a comparison with Almeida's $q_c(T_w)$ for argon → claim (D) and a direct benchmark against the baseline paper.
6. A coarse sweep in $J_{\rm app}$ (and $T_{\rm bulk}$) → claim (E).
7. The ion impact energy either fixed or removed from the title.

If you can only do part of this, (A) with item 2 and (C) with item 4 give the clearest distinction from Almeida at the lowest cost, because the code already supports both.

## Question 8

See the plotfiles in `PlasmeX/exe`, write a fast script that evaluates the time scales of study T1 (`literature_review.md` §7.4), and judge whether the result looks good and publishable.

*(Rewritten on request: the first version used terms such as "mode", "coupled mode" and "diffusive mode" without explaining them. This version builds from the basics.)*

## Response

### 1. What we are trying to measure, and why

For the paper's positioning (`literature_review.md` §7) we want one number: **how long the near-electrode layer takes to adjust when something changes**, for example when the arc current changes. Call it the **response time** $\tau$ of the layer.

Why it matters: the LTE arc model changes the conditions at the edge of the layer over time (current, bulk temperature).

- If the layer adjusts **much faster** than those changes, it is always in its steady state. A table of steady PlasmeX results is then enough as a boundary condition for the LTE model.
- If the layer adjusts **about as fast as, or slower than**, those changes, a steady table is wrong, and a transient electrode model is needed.

So $\tau$ is the number that decides whether your "transient electrode patch" argument is needed or not.

### 2. How the script measures $\tau$ (the basic idea)

Think of a cup of hot coffee cooling in a room. Its temperature does not drop linearly. It drops fast at first and then slows down, approaching room temperature more and more gently. Mathematically:
$$
T(t)=T_{\rm room}+\big(T_{\rm start}-T_{\rm room}\big)\,e^{-t/\tau}.
$$
$\tau$ is the "response time": after each interval $\tau$, the *remaining* difference shrinks by a factor $e\approx2.7$. After $3\tau$ only 5% of the difference is left.

The simulation behaves the same way. It starts from a guessed initial profile and moves towards the steady solution. For any quantity $q$ we track (the voltage drop, the wall electric field, $n_e$ at the wall, …), the approach at late times looks like
$$
q(t)=q_{\rm steady}+A\,e^{-t/\tau}.
$$
The script (`PlasmeX/exe/timescales.py`):
1. reads every plotfile and extracts about ten such quantities (voltage drop, wall field, wall $n_e$, wall $T_e$, ion fraction of the wall current, gas temperature near the wall, total gas content, gas velocity, …);
2. fits the formula above to the **last part** of the run, which gives $q_{\rm steady}$ and $\tau$ for each quantity;
3. draws each quantity on a log scale as distance-to-steady-value vs time. A pure exponential shows as a **straight line**, so you can see at a glance whether the formula fits.

**Why fit only the end of the run?** A real system has several processes adjusting at the same time, for example charged particles rearranging, gas cooling and gas flowing in. Each process has its own response time. Early on they all contribute and the curve is messy. The fast processes die out quickly, so near the end **only the slowest process is left**. The end of the run therefore shows the slowest response time of the system, whatever the initial guess was. The slowest time is also the one that matters for the question in §1.

(In the earlier version I called each of these separate processes a **"mode"**. That is all the word means: one independent way in which the system relaxes, with its own $\tau$.)

Usage:
```bash
cd PlasmeX/exe
python3 timescales.py plt*                  # table + timescales.png + timescales.csv
python3 timescales.py --tmin 8e-6 step/plt* # step test: ignore frames before the step
```
It reads the AMReX files directly, without yt, and takes about 2 s for 41 files. I checked that it reads exactly the same numbers as yt.

### 3. What your current run gives

The run is `anode_inputs`: current density $-10^7$ A m$^{-2}$, bulk $T_e=12\,000$ K, wall at 1000 K, domain length $L=200$ µm, simulated from 0 to 8 µs.

| quantity | response time $\tau$ | how much of the change is still left at 8 µs |
|---|---|---|
| voltage drop across the domain (2.95 → 2.96 V) | 2.70 µs | 0.4% |
| electric field at the wall | 2.55 µs | 3% |
| $n_e$ at the wall | 2.79 µs | ~0 |
| $T_e$ at the wall | 3.05 µs | 0.2% |
| ion share of the wall current | 2.48 µs | 3% |
| total amount of gas in the domain | 2.67 µs | 5% |
| gas velocity at the arc end (→ 0) | 2.65 µs | 2% |
| gas temperature 20 µm from the wall | 2.57 µs | 23% (not finished) |

**What is good:**
- All quantities give **the same** $\tau\approx2.5$–$3$ µs. So they are all being held back by one and the same slow process.
- The straight lines on the log plot are very clean (fit quality $R^2>0.9999$). $\tau$ changes by less than 3% if I fit a shorter piece of the run.
- Checks: the current at the wall equals the imposed current, and the gas that enters the domain equals the increase of gas inside it. So nothing is being created or lost numerically.
- 2.7 µs is inside the 1–10 µs range I estimated in `literature_review.md` §7.2.

**What the slow process is, physically:** the charged particles (electrons, ions) rearrange very quickly, within the first 0.2 µs, which is the first saved frame. What takes microseconds is the **gas**:
- near the cold wall the gas cools;
- cooler gas at the same pressure is denser, so it needs more atoms;
- those atoms flow in slowly from the arc side;
- as the gas changes, the ion and electron profiles follow it, so the voltage and wall field keep changing slightly until the gas has settled.

(The neutral-atom velocity `UNX` rises late in the run. That is not a problem: in steady state the neutral atoms drift slightly against the ions, so `UNX` settles at a small nonzero value. The script uses the overall gas velocity `HVELX`, which correctly goes to zero.)

### 4. The main problem: the result may depend on the size of our box

**Basic fact about diffusion and heat conduction.** When something spreads by diffusion (particles) or conduction (heat), the time it takes to cross a distance $L$ is roughly
$$
t\sim\frac{L^2}{D},
$$
where $D$ is the diffusion coefficient (or, for heat, the thermal diffusivity). The important part is the **square**: double the distance and the time is **four times** longer. This is why a thin steak cooks much faster than a thick one.

**Why this matters here.** PlasmeX simulates a box of length $L=200$ µm. At the far end (the arc side) we *impose* the arc conditions: $n_e$, $T_e$, gas temperature. That is only correct if the layer has **fully ended** before 200 µm, i.e. the plasma there has really reached arc conditions. In the final profile it has not:

| distance from the wall | $n_e$ | gas temperature |
|---|---|---|
| 200 µm (imposed arc value) | $6.46\times10^{22}$ | 11 994 K |
| 150 µm | $5.76\times10^{22}$ (89%) | 10 510 K |
| 100 µm | $4.89\times10^{22}$ (76%) | 8 657 K |
| 50 µm | $3.73\times10^{22}$ (58%) | 6 107 K |

The profiles are still changing when they reach the box edge. So the layer is **wider than the box**, and the edge of the box is cutting it off.

The consequence: the slow process we measured is the gas settling **over the whole box**. If that is a diffusion or conduction process, its time grows like $L^2$. Then the 2.7 µs would be **a property of our choice of box size**, not of the physics:
- with $L=100$ µm we might get ~0.7 µs;
- with $L=400$ µm, ~11 µs.

A referee would immediately ask "how does $\tau$ change with the domain length?", and at present we cannot answer.

**Can we tell which process it is?** I calculated the response times of the two obvious single processes on the final profile, using the code's own coefficients:
- heat conduction in the gas alone, across the box: **0.52 µs**;
- diffusion of ions and electrons together alone (called *ambipolar* diffusion), across the box: **0.63 µs** if the arc end is held fixed, **5.2 µs** if it is left free.

Neither equals 2.7 µs. The likely reason is that these processes do not act separately. Gas cooling changes the gas density, which changes ion diffusion, which changes the field, which changes the heating, and so on. When several processes drive each other like this, the combined response time differs from each one alone. (That is what I meant by a **"coupled mode"**: one slow relaxation produced by several processes acting together.) Each of them, however, gets slower with $L^2$, so the combined one very probably does too.

**What to do.** Run the same case with $L=100$ µm and $L=400$ µm, keeping the same cell size (i.e. change `geometry.prob_hi` and `amr.n_cell` together), and compare $\tau$.
- **$\tau$ stays about the same:** the response time belongs to the layer, and 2.7 µs is a real, publishable number.
- **$\tau$ grows roughly as $L^2$:** the slow response is set by how thick the gas layer is. That thickness is not decided by PlasmeX but by the arc (the LTE solution). This is still a useful result for your positioning: it means the electrode region responds on the time scale of the arc's own thermal layer, which is exactly where quasi-stationary coupling fails. But the box length must then be taken from an LTE solution (where the LTE temperature profile near the anode ends), not chosen by hand.

### 5. The other problems, in plain terms

**(a) We measured a start-up, not a response to a change.**
The run starts from a guessed profile. The *rate* at which it settles ($\tau$) is meaningful, but the *size* of each change is not. For the paper we want a real test: the arc is steady, the current is suddenly raised by 10%, and we watch how the voltage and the heat flux follow. This needs **no code change**, because a restart re-reads the `prob.*` values:
```bash
cd PlasmeX/exe
mkdir step_up && cd step_up && cp ../*.dat .
mpirun -np 1 ../main1d.gnu.MPI.ex ../anode_inputs amr.restart=../chk400000000 \
    prob.arc_Japp=-1.1e7 max_step=900000000 amr.plot_int=2500000 \
    plasmex.steady_check_int=0 plasmex.check_solution_int=1
python3 ../timescales.py --tmin 8e-6 'plt*'
```
- `prob.arc_Japp=-1.1e7` is the +10% step.
- `amr.plot_int=2500000` saves a frame every 50 ns.
- `plasmex.steady_check_int=0` stops the steady-state detector from ending the run early.

Also run a −10% step (to check that the response is symmetric) and a step in `prob.arc_Te` (a change of the bulk plasma instead of the current).

**(b) The fast part is invisible.**
Everything the electrons and ions do happens before the first saved frame (0.2 µs). We see only the slow gas part. Both times matter for the paper: the fast one says how quickly the sheath reacts, the slow one how quickly the whole layer settles. For the step runs, save frames much more often at the start (every 1–10 ns for the first ~0.5 µs, i.e. `amr.plot_int` of 50 000–500 000), then less often.

**(c) The time step has not been shown to be small enough.**
PlasmeX's `CLAUDE.md` warns that one time-step limit (the dielectric relaxation limit) is printed but not enforced, and that too large a step gave a wrong sheath voltage in one of their tests (0.715 V instead of 1.702 V). For a *time-dependent* result we must show that the answer does not change when the step is made smaller: repeat one step test with `plasmex.time_step=1e-14` (half the current value), and one with one more refinement level (`amr.max_level=4`). $\tau$ and the voltage change should stay the same.

**(d) The electrode boundary condition is simplified.**
At the wall, the electron flux is fixed at its thermal value, and `CLAUDE.md` notes that the electrode "cannot regulate the current it accepts". For a steady result this can be defended, but a paper about how the voltage *responds in time* will be questioned on it. Either show that $\tau$ hardly changes when the ion wall condition is switched (`plasmex.ion_bc_hi` 1 vs 3), or state the limitation clearly. The cathode also needs its own run. There the wall temperature is recomputed instantly each step (the electrode has no heat capacity in the model), so the slow heating of the electrode itself is not included yet.

### 6. Verdict

- **The script works and the first result is clean.** One well-defined response time in every quantity, a good fit, and current and mass conserved.
- **The number 2.7 µs cannot be published yet,** mainly because:
  1. the layer does not fit in the 200 µm box, so $\tau$ may just reflect the box size (§4);
  2. it comes from a start-up rather than a controlled change (§5a).
- **The fix is cheap:** about 6–8 restarts from the existing checkpoint, and no code changes:
  - two box lengths;
  - ±10% current steps and one bulk-temperature step, with frequent early frames;
  - one smaller-time-step check and one finer-mesh check.
- **With those runs, T1 becomes a clear one- or two-figure section of the paper:** the fast (sheath) and slow (gas) response times against current, how the voltage and heat flux follow a current step compared with the steady prediction, and from that the rule for when a steady-table boundary condition is good enough for a transient LTE model.

## Question 9

For a diffusion process the relaxation time depends on the domain? If so, since I solve a diffusion equation for the gas temperature, how can the result not depend on the box? Shouldn't every problem depend on the domain length?

## Response

**Short answer: no.** A diffusion problem has an $L$-dependent relaxation time only when **nothing but diffusion** acts, so that the profile has to stretch across the whole box. If the equation also contains a **local term that pulls the solution back to a fixed value at every point** (a reaction, a sink, an exchange), that term sets its own length scale $\delta$. Once the box is much longer than $\delta$, the box size no longer matters. Your gas-temperature equation has such a term, but only for part of the problem. That is the whole story, and below it is shown with a formula and with numbers from your run.

### 1. Pure diffusion: the time *does* depend on $L$

Take
$$
\frac{\partial n}{\partial t}=D\,\frac{\partial^2 n}{\partial x^2},\qquad n \text{ fixed at } x=0 \text{ and } x=L .
$$
The steady profile is a straight line from one end to the other, so it fills the whole box. A disturbance can only disappear by diffusing out through the two ends, and that takes the time to diffuse across the box. The slowest disturbance is a half sine wave, $\sin(\pi x/L)$, and it decays as $e^{-t/\tau}$ with
$$
\tau=\frac{L^2}{\pi^2 D}.
$$
Double $L$ and $\tau$ is four times longer. Here the answer depends on the box, as you say.

### 2. Diffusion plus a local restoring term: it does *not* (once the box is large enough)

Now add a term that pulls $n$ back towards a local value $n_0$ at a rate $1/\tau_r$ at every point (for example recombination, or energy exchange with another species):
$$
\frac{\partial n}{\partial t}=D\,\frac{\partial^2 n}{\partial x^2}-\frac{n-n_0}{\tau_r}.
$$
Put the same half sine wave in. Each point now loses its disturbance in two ways: by diffusion (as before) and directly through the local term. The two rates simply add:
$$
\frac{1}{\tau}=\frac{1}{\tau_r}+\frac{\pi^2 D}{L^2}.
$$
The two terms balance at the length
$$
\delta=\sqrt{D\,\tau_r},
$$
which is a property of the physics (of $D$ and $\tau_r$), not of the box. So:
- **Box short** ($L\ll\delta$): diffusion wins, $\tau\approx L^2/(\pi^2D)$, and the answer depends on the box.
- **Box long** ($L\gg\delta$): the local term wins, $\tau\approx\tau_r$, and **the box drops out**. The steady profile also stops filling the box: it is a boundary layer of thickness $\sim\delta$ next to the wall, flat beyond it.

A concrete example: for $L=10\,\delta$, $\tau$ differs from $\tau_r$ by only 10%; for $L=3\delta$, by about 100%.

This is why a sheath does not depend on the box: its thickness is set by the Debye length, a physical length. It is also why a boundary layer in fluid mechanics has a thickness that does not depend on the size of the room.

**So the correct statement is:** a relaxation time is independent of the domain **when the domain is much larger than the physical length scales of the problem**. The whole point of the $L$ test is to check whether we are in that regime.

### 3. Your gas-temperature equation

In PlasmeX the gas temperature obeys (simplified):
$$
n_h c_v\frac{\partial T_h}{\partial t}=\frac{\partial}{\partial x}\Big(k_h\frac{\partial T_h}{\partial x}\Big)+C_{el}\,(T_e-T_h)+Q_{in}.
$$
The middle term, the collisional energy transfer from electrons to the gas, is a local restoring term: it pulls $T_h$ towards $T_e$ at every point, with a rate
$$
\frac{1}{\tau_{ex}}=\frac{C_{el}}{n_h c_v}.
$$
That gives the gas temperature its own physical length, $\delta_T=\sqrt{\alpha\,\tau_{ex}}$, with $\alpha=k_h/(n_h c_v)$ the thermal diffusivity of the gas.

Rough numbers from your final profile. These use the textbook (NRL) electron–ion collision frequency and neglect the electron–neutral part, so they are order of magnitude, not the code's exact value:

| distance from the wall | $T_e$ | $T_h$ | $\tau_{ex}$ | $\alpha$ [m²/s] | $\delta_T=\sqrt{\alpha\tau_{ex}}$ |
|---|---|---|---|---|---|
| 200 µm (arc end) | 12 001 K | 11 994 K | 0.4 µs | 0.020 | ~90 µm |
| 100 µm | 10 046 K | 8 657 K | 0.8 µs | 0.012 | ~100 µm |
| 50 µm | 8 521 K | 6 107 K | 1.7 µs | 0.0067 | ~105 µm |
| 20 µm | 7 124 K | 3 781 K | 4.6 µs | 0.0030 | ~115 µm |

Two things follow:
1. **The exchange times, 0.4–4.6 µs, bracket the measured $\tau\approx2.7$ µs.** So the slow response we measured may well be this electron–gas energy exchange, a physical time that does not depend on the box. That is good news.
2. **But the physical length $\delta_T\approx100$ µm is not much smaller than the box, $L=200$ µm.** The box is only about $2\delta_T$ long, which is the "box short or comparable" regime of §2. The box can therefore still shift $\tau$ noticeably. That matches what we see: $T_h$ is still 10 500 K at 150 µm, i.e. the layer has not ended inside the box.

### 4. The part that has no restoring term

The exchange term only pulls $T_h$ towards $T_e$. It controls the **difference** $T_e-T_h$ (how far the gas is from equilibrium with the electrons). It does not decide the **overall temperature level** of the plasma. That level is set by the balance of Joule heating, conduction and radiation in the arc. In a real arc it varies over millimetres, the scale of the arc column that the LTE model computes. In PlasmeX that level is imposed at the arc end of the box.

So the near-electrode problem really has two parts:
- **The non-equilibrium part** ($T_e\neq T_h$, $n_e$ away from its Saha value, the sheath). It has physical lengths (Debye length, $\delta_T\sim100$ µm, the ionization-layer thickness) and physical times (0.4–5 µs, and much less for the sheath). If the box is several times longer than these lengths, the times do **not** depend on the box. This part belongs to PlasmeX.
- **The equilibrium part** (the overall temperature profile of the arc near the electrode). It has no local restoring term, so its response time depends on how far it extends, like pure diffusion. That extent is set by the arc, not by the electrode region. This part belongs to the LTE model.

That split is exactly the division of labour between an LTE arc model and an electrode-layer model, and it supports your positioning: the layer model should supply the non-equilibrium response times, and the LTE model the slow equilibrium one.

### 5. What this changes in the plan (Question 8, §4)

- The $L$ test is still the right test, but the box lengths should be chosen against $\delta_T\approx100$ µm: run $L=200$, $400$ and $800$ µm (about $2$, $4$ and $8\,\delta_T$).
- **Expected outcome if the picture above is right:** $\tau$ changes noticeably from 200 to 400 µm and hardly at all from 400 to 800 µm, settling at a value of a few µs set by the energy exchange. That settled value is the publishable number.
- **If instead $\tau$ keeps growing like $L^2$** up to 800 µm, then the slow response is the equilibrium part of §4. In that case the result to report is that the slow response of the electrode region is set by the arc's thermal layer, and the box must be sized from an LTE solution.
- A cheap extra check that needs no new runs: compute $\tau_{ex}$ with the code's own exchange frequency (`ElasticExchangeFreq` in `src/coulomb.H`) along the profile and compare it with the fitted $\tau$. I used a textbook formula above, which can be off by a factor of about 2.

## Question 10

I have already shown with length scales that 200 µm is a typical electrode region. Why not just accept $\tau$ as the time scale for this type of domain? Or is it better to verify that it is really domain independent?

## Response

**Verify it.** The length-scale argument is the right reason for *choosing* 200 µm, but it does not show that 200 µm is *large enough* for $\tau$ to be independent of it. There are four reasons, and the check is cheap.

### 1. Your own length scales say the box is comparable to the layer, not much larger

The paper (`sections/model.tex`, §"Length scales and the thickness of the near-electrode layer") gives:

| scale | value in the paper |
|---|---|
| electron–heavy energy-relaxation length $\lambda_T$ | $\approx2\times10^{-4}$ m |
| diffusion length $\Lambda=(kT_e/e)/E$ | $\approx5\times10^{-4}$ m |
| boundary-layer thickness from Dinulescu & Pfender | 0.2–0.6 mm |
| chosen domain | $2\times10^{-4}$ m, "comparable to both $\Lambda$ and $\lambda_T$" |

"Comparable" is the key word. Question 9 showed that a result stops depending on the box only when the box is **several times** longer than the physical length. With the box about equal to that length, the box still matters. Your numbers even put $\Lambda$ at 2.5 times the box.

The simple model of Question 9 shows how much it can matter. With the arc end held fixed,
$$
\frac{1}{\tau}=\frac{1}{\tau_r}\Big(1+\frac{\pi^2\delta^2}{L^2}\Big),
$$
so the measured $\tau$, as a fraction of the true physical $\tau_r$, is:

| box length $L$ | $\tau/\tau_r$ |
|---|---|
| $1\,\delta$ | 0.09 |
| $2\,\delta$ | 0.29 |
| $4\,\delta$ | 0.62 |
| $8\,\delta$ | 0.87 |

With $\delta_T\approx100$ µm (Question 9) and $L=200$ µm, we are near the $2\delta$ row. The box could therefore be making $\tau$ look several times shorter than the physical value. This model is simplified, so the table shows the size of the risk, not a correction to apply. A factor of 3 matters for the question the paper wants to answer ("is $\omega\tau\ll1$?").

### 2. The simulation itself shows the layer is cut off

At the arc end the boundary imposes $T_h=12\,000$ K, but at 150 µm $T_h$ is still 10 500 K, and $n_e$ is 89% of its arc value. The profiles have not flattened out before the boundary. So the boundary is not simply "outside the layer": it is holding up the end of the layer. Anything that depends on the whole layer (the slow $\tau$, the heavy-particle heat conducted towards the wall, the total voltage drop) can feel where that boundary is placed.

### 3. The length-scale estimate is only an order of magnitude

The estimates give $10^{-4}$ m to within a factor of a few (0.2–0.6 mm in the literature you cite). If $\tau$ scales anywhere near $L^2$, choosing 200 µm instead of 500 µm changes it by up to a factor of 6. A time that changes by a factor of 6 with an arbitrary modelling choice cannot be quoted as "the response time of the electrode region". It would have to be quoted as "the response time of our 200 µm model", which is much weaker.

### 4. A referee will ask, and the same check also protects your steady results

Domain-length sensitivity is the standard question for any result computed on a truncated domain, just as mesh convergence is. The same runs also answer a question about the **steady** results the paper already reports: do the fall voltage and the wall heat flux change when the box grows? If they do not, the runs **numerically confirm** the 200 µm choice made in `model.tex`, and you can say so in one sentence. That strengthens the length-scale section instead of leaving it as an estimate.

### What to do

- **Runs:** the same `anode_inputs` case with $L=400$ µm and $L=800$ µm. Keep the refinement near the wall identical: the same `tag_electrode_depth`, `max_level` and finest cell size, and `amr.n_cell` scaled with $L$ so the base cell size stays the same, or double the base cell and add one refinement level. In 1D the cost grows roughly with the number of cells, so 2–4 times per run. Run each until the slow quantities have settled (about 5τ, 15–20 µs), with `plasmex.steady_check_int=0`, then apply `timescales.py`.
- **Compare across $L$ = 200, 400, 800 µm:**
  1. $\tau$ from `timescales.py`;
  2. the steady fall voltage, $E_{\rm wall}$, $n_e$ and $T_e$ at the wall, and the wall heat flux;
  3. where $n_e$ and $T_h$ reach, say, 99% of their arc values.
- **How to read it:**
  - $\tau$ and the wall quantities change little from 400 to 800 µm → that value is physical and publishable. If 200 µm agrees with it too, the paper's domain is confirmed as it stands.
  - They change from 200 to 400 µm but not from 400 to 800 µm → use 400 µm (or more) for the paper and state the domain study.
  - They keep changing up to 800 µm → the slow response belongs to the equilibrium arc profile (Question 9, §4). Then the right statement is that the electrode region's slow response is set by the arc's thermal layer, and the box has to be sized from an LTE solution.

**If you cannot afford the runs,** the fallback is to quote $\tau$ explicitly as a property of the 200 µm model and to state the $L$ sensitivity as a limitation. You could also support it with the box-independent estimate $\tau_{ex}\approx0.4$–$5$ µs, computed from the electron–gas exchange rate (Question 9). That is defensible but clearly weaker, and it leaves the T1 claim of the paper resting on an estimate.

## Question 11

Results of `timescales.py` for the run at $J=6\times10^6$ A m$^{-2}$ (pasted output; same 200 µm box assumed).

## Response

### 1. The run itself: very clean

- It was run long enough: every diagnostic has less than 0.5% of its change left (only $p_{\rm wall}$ is at 2%). The earlier $10^7$ run stopped at 8 µs with the gas temperature still 23% away.
- All the independent quantities agree on one response time:

| quantity | $\tau$ [µs] |
|---|---|
| voltage drop | 3.591 |
| $E_{\rm wall}$ | 3.588 |
| $n_e$ wall | 3.595 |
| $T_e$ wall | 3.602 |
| $T_h$ 20 µm from the wall | 3.506 |
| gas content | 3.517 |
| ion share of the wall current | 3.574 |
| gas velocity (→ 0) | 3.520 |

  **Median $\tau=3.59$ µs, spread ±1.3%.** This is as clean a single relaxation as one can get from a simulation.

- **Script fix.** The last line of the output ("tau = 3.881 us") picked the *largest* value, which here was $p_{\rm wall}$, the least reliable fit. The script now prints the median and range instead. It excludes $p_{\rm wall}$ (a tiny change riding on sound-wave ringing) and $J_{\rm wall}$ (imposed by the boundary condition, so its fit is to noise). For this output it would print about 3.59 µs, range 3.51–3.60 µs.

### 2. Comparison with $J=10^7$ A m$^{-2}$

| | $J=10^7$ | $J=6\times10^6$ |
|---|---|---|
| $\tau$ (median) | 2.68 µs (±10%, run too short) | 3.59 µs (±1.3%) |
| voltage drop across the box | 2.96 V | 2.09 V |
| $E_{\rm wall}$ | $2.5\times10^6$ V/m | $0.89\times10^6$ V/m |
| $n_e$ at the wall | $5.4\times10^{20}$ | $3.2\times10^{20}$ m$^{-3}$ |
| $T_e$ at the wall | 5 590 K | 5 600 K |
| $T_h$ 20 µm from the wall (steady) | ~3 740 K | ~2 900 K |
| gas content $\int n_h\,dx$ | $2.19\times10^{20}$ | $2.85\times10^{20}$ m$^{-2}$ (+30%) |

**Lower current → slower response:** $\tau$ grows by a factor 1.34 when $J$ drops by a factor 0.6. With two points only, that is roughly $\tau\propto J^{-0.6}$. Do not quote this exponent until there are at least three or four currents.

**Plain-language reason:** less current means less Joule heating, so the gas near the wall is cooler (2 900 K instead of 3 740 K at 20 µm) and denser (30% more gas in the box). In cooler, denser gas:
- the electrons pass their energy to the gas more slowly, per gas particle (the exchange time $\tau_{ex}$ of Question 9 grows towards the cold wall);
- heat spreads more slowly (lower thermal diffusivity $\alpha$).

Both make the layer respond more slowly. So the trend makes physical sense.

### 3. What this does *not* settle: the box question

Both explanations above fit the trend:
- the **physical** one (the exchange time grows in cooler gas), for which $\tau$ is independent of the box;
- the **box-limited** one ($\tau\sim L^2/\alpha$, and $\alpha$ is lower in cooler gas), for which $\tau$ depends on the box.

Changing $J$ cannot tell them apart; only changing $L$ can (Question 10).

One more reason to do the $L$ test: **at lower current the layer gets thicker**. The paper's diffusion length is $\Lambda=(kT_e/e)/E$, and the arc field $E=J/\sigma$ is proportional to $J$. So with the paper's numbers,
$$
\Lambda\approx5\times10^{-4}\ \text{m at }10^7\ \text{A m}^{-2}\quad\longrightarrow\quad\Lambda\approx8\times10^{-4}\ \text{m at }6\times10^6\ \text{A m}^{-2}.
$$
That is four times the 200 µm box. The box question is therefore **more** serious at low current, not less.

### 4. Next steps, in order

1. **Box test at $J=6\times10^6$** (this setup already runs cleanly to steady state): $L=400$ and $800$ µm, the same near-wall refinement, and the same run length and frame spacing. Compare $\tau$, the steady wall quantities, and where $n_e$ and $T_h$ reach their arc values (Question 10).
2. **Rerun $J=10^7$ for longer** (≥ 20 µs, like this run), so that both currents are at the same quality.
3. **Current scan at the box length validated in step 1:** at least 4 currents, e.g. $3\times10^6$, $6\times10^6$, $10^7$ and $2\times10^7$ A m$^{-2}$. This gives $\tau(J)$, which is the T1 figure for the paper.
4. **Fast part:** your frames are 1 µs apart, and $n_e$ at the wall reaches 90% by the first frame, so the fast (sheath) response is still invisible. For the step tests (Question 8, §5a), save frames every 1–10 ns for the first ~0.5 µs.

If you point me to the directory of the $6\times10^6$ plotfiles, I can check at once how far the layer extends in that run (the values of $n_e$ and $T_h$ at 150 µm compared with the arc end). That shows directly how much the 200 µm box is cutting it off.
