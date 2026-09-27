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
