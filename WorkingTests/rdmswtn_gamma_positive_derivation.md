# Gamma drift: positive analytic CDF/survivor representations

2026-10-19. Follow-up to `rdmswtn_gamma_feasibility.md`. Research only;
production kernels and parameter mappings are unchanged.

## Finding

There is a stronger result than the original alternating Taylor series:
**two complementary positive series cover every positive real gamma shape
and rate**. They give exact point-start CDF/survivor representations without
quadrature. One uses ordinary incomplete gamma functions; the other uses a
short recurrence of scaled gamma/Gaussian moments. In particular, the
survivor need not be calculated by cancellation of enormous signed terms.

This is an exact analytic-series solution, not a finite expression in a fixed
number of ordinary special functions. I have not identified such a finite
expression for unrestricted real shape. Nor do these results yet establish
a universally inexpensive likelihood: term counts can grow with shape, and
tiny CDFs need a direct branch rather than subtraction from the survivor.
The direct early-time CDF expansion below has an analytic remainder bound
and succeeds in several of the previously failing cases.

The practical recommendation is to pursue these representations rather than
the original entire-function expansion. A compiled, adaptive implementation
and a performance gate remain necessary before adopting gamma variability.

## Definitions

Let `V ~ Gamma(k, rate tau)`, with `k,tau,d,s,t > 0`. Set

\[
 \lambda=\frac d{s^2},\quad q=\frac\tau\lambda=\frac{\tau s^2}{d},
 \quad A=\frac{d^2}{2s^2t},\quad p=\frac{\tau s}{\sqrt t},
 \quad z=\frac{(\tau-\lambda)s}{\sqrt t}.
\]

Write `F_k,S_k,f_k` for the point-distance gamma-mixture law at fixed `tau`.
For any positive real `nu`, define

\[
 j_\nu=E_{\Gamma(\nu,\tau)}
 \left[\frac1{s\sqrt t}\phi\left(\frac{d-Vt}{s\sqrt t}\right)\right],
 \qquad f_\nu=\frac dt j_\nu.
\]

It is useful to retain the two conditional Wald CDF components separately:

\[
 T_\nu=E\,\Phi\left(\frac{Vt-d}{s\sqrt t}\right),\qquad
 H_\nu=E\left[e^{2\lambda V}
        \Phi\left(-\frac{d+Vt}{s\sqrt t}\right)\right],
 \quad F_\nu=T_\nu+H_\nu.
\]

`H` is the reflection component here, not the entire function called `H` in
the earlier feasibility note. All these probabilities/components are
nonnegative. The formulas below allow any real `k>0`, including shapes below
one. Mean/SD coordinates are still `k=(v/sv)^2`, `tau=v/sv^2`.

## 1. A surprising exact inverse-gamma case

At `tau=lambda`, `z=0`. The duplication formula for the gamma function
simplifies the existing cylinder density to

\[
 \boxed{T_{\rm decision}\sim
       \operatorname{InvGamma}\left(\frac{k+1}{2},\frac{d^2}{2s^2}\right).}
\]

Here the inverse-gamma density is
`a^alpha/Gamma(alpha) * t^(-alpha-1) exp(-a/t)`.
If `P,Q` denote lower/upper regularized incomplete gamma functions,

\[
 F_k(t)=Q((k+1)/2,A),\qquad S_k(t)=P((k+1)/2,A).
\]

This is a finite closed form for **arbitrary real shape**, on the parameter
surface `tau=d/s^2`. It supplies a useful validation identity as well as the
base family for the first general representation.

## 2. Positive inverse-gamma mixture for q <= 1

Expand a lower-rate gamma density into higher-rate gamma densities:

\[
 \Gamma(v;k,\tau)=\sum_{n=0}^\infty
 w_n\Gamma(v;k+n,\lambda),\qquad
 w_n=\frac{(k)_n}{n!}q^k(1-q)^n.
\]

For `0<q<1` these are negative-binomial probabilities and sum to one.
Every component is the inverse-gamma law just derived. Consequently

\[
 \boxed{F_k(t)=\sum_{n=0}^\infty w_n
 Q\left(\frac{k+n+1}{2},A\right),\quad
 S_k(t)=\sum_{n=0}^\infty w_n
 P\left(\frac{k+n+1}{2},A\right).}
\]

This gives **both tails directly**, with no subtraction or numerical
integration. Tonelli's theorem justifies the interchange. Generate weights
by `w_(n+1)/w_n=(1-q)(k+n)/(n+1)` and use log weights or start near their mode
when needed. Do not separately form `Gamma(k+n)` or `q^k` in ordinary scale.

The omitted negative-binomial mass bounds the absolute error of either tail.
For a relative tail guarantee, the omitted bound must be small compared with
the accumulated tail, or a sharper bound using the monotone `P,Q` factors
must be derived. The prototype uses an absolute bound; it does not certify
relative accuracy for arbitrarily tiny probabilities.

The incomplete-gamma shape advances by one on each parity subsequence:
`P(alpha+1,A)=P(alpha,A)-exp(-A) A^alpha/Gamma(alpha+1)` and the corresponding
`Q` update has a plus sign. Thus two parity streams can replace repeated
special-function calls. Update the smaller tail in its stable direction;
an unstable subtraction would defeat the purpose. This optimization is not
implemented in the research script.

## 3. Positive moment series for q >= 1

This is the main new result for the complementary domain, including the
shape-100 overflow cases.

For a differentiable bounded `h`, gamma integration by parts gives

\[
 E_{\nu+1}h-E_\nu h=\frac1\tau E_{\nu+1}h'.
\]

The derivative of the **combined conditional Wald CDF** with respect to
drift is particularly simple:

\[
 \partial_v F(t\mid v)=2\lambda e^{2\lambda v}
           \Phi\left(-\frac{d+vt}{s\sqrt t}\right).
\]

The Gaussian derivatives of the two conditional terms cancel. Therefore

\[
 F_{\nu+1}-F_\nu=\frac{2\lambda}{\tau}H_{\nu+1},\quad
 S_\nu-S_{\nu+1}=\frac{2\lambda}{\tau}H_{\nu+1},\quad
 T_{\nu+1}-T_\nu=\frac t\tau j_{\nu+1}.
\]

The derivative of the reflection term is `2 lambda g - t j(v)`.
Another application of the same identity yields

\[
 (\tau-2\lambda)H_{\nu+1}=\tau H_\nu-tj_{\nu+1}.
\]

Set `r=1-2/q`. Solve this recurrence **backwards in shape**:

\[
 H_\nu=\frac t\tau\sum_{n=1}^\infty r^{n-1}j_{\nu+n}.
\]

At fixed rate, increasing shape sends drift to infinity, so `S_(k+N)` and
`H_(k+N)` tend to zero. For `q>1`, `|r|<1`, so the backward boundary term
vanishes. Telescoping the survivor recurrence and combining the two sums
gives

\[
 \boxed{S_k(t)=\frac t\tau\sum_{m=2}^\infty
             \bigl[1-r^{m-1}\bigr]j_{k+m}.}
\]

For **every `q>=1`, all coefficients are nonnegative**. This remains true
when `1<q<2` and `r` is negative: the survivor series is positive even though
the separate reflection series alternates. At `q=1`, `r=-1`, only even `m`
contribute, with coefficient two; the result agrees with the inverse-gamma
identity. The sums still converge there because the moments decay faster
than a geometric sequence as shape tends to infinity.

At `q=2`, `r=0`: `H_nu=(t/tau)j_(nu+1)` and the survivor formula simplifies
further. There is no singular quotient at `tau=2d/s^2`.

The exact CDF is `1-S`, evaluated as `-expm1(logS)` when appropriate. This is
well conditioned when the CDF is not tiny. Positivity of the survivor series
does **not** make complementing it accurate for tiny CDFs; Section 5 addresses
some of that domain directly.

### A truncation bound

Let `J_nu(z)=integral y^(nu-1) exp(-y^2/2-zy) dy` and
`R_nu=J_(nu+1)/J_nu`. Log convexity of positive moments shows that `R_nu`
increases with `nu`. Integration by parts gives

\[
 R_\nu=\frac\nu{z+R_{\nu+1}},\qquad
 \frac{j_{\nu+1}}{j_\nu}=\frac p{z+R_{\nu+1}}.
\]

For `z>=0` the latter ratios decrease. Once `rho_N=j_(k+N+2)/j_(k+N+1)<1`,
the remainder after `m=N` is bounded by

\[
 0\le\mathcal R_N\le
 \frac t\tau\frac{B j_{k+N+1}}{1-\rho_N},\qquad
 B=\begin{cases}1&q\ge2,\\2&1\le q<2.\end{cases}
\]

This is an **exact-arithmetic** truncation bound, separate from roundoff
and continued-fraction initialization errors. The prototype records it.
It gives a principled stopping test, unlike the old series' small-term rule.

## 4. Prevent cylinder underflow and prefactor overflow

Do not evaluate a large-order raw `D_-k(z)` and multiply it by powers and
exponentials. Instead write `k=rho+M`, with integer `M>=0`, `0<rho<=1`.
Use the positive continued fraction above to obtain moment ratios backwards
from a sufficiently large order. Initialize only the **fractional-order**
scaled function

\[
 \frac{J_\rho(z)}{\Gamma(\rho)}=e^{z^2/4}D_{-\rho}(z).
\]

For `0<=z<=2`, its raw cylinder factor is representable. A direct SciPy
fractional-order call is not reliable throughout the moderate range: the
expanded checks found about `5.4e-10` relative error near `z=sqrt(40)`.
For large positive `z`, use

\[
 e^{z^2/4}D_{-\rho}(z)
 =z^{-\rho}\left[\sum_{n=0}^{N-1}
 \frac{(-1)^n(\rho)_{2n}}{n!(2z^2)^n}+\mathcal R_N\right].
\]

The next term bounds this remainder: expand `exp(-y^2/2)` under the positive
Laplace integral and integrate its Taylor remainder. Since `rho<=1`, a
large-`z` seed needs few terms even when the target `k` is very large.
For `2<z<=20`, the prototype instead continues from `z*=32` using the
**positive convergent identity**

\[
 J_\rho(z)=\sum_{n=0}^\infty
       \frac{(z_*-z)^n}{n!}J_{\rho+n}(z_*).
\]

Its moment ratios come from the same positive continued fraction. Because
`R_(rho+n)(z*) <= (rho+n)/z*` and `rho<=1`, successive summand ratios are at
most `(z*-z)/z*<1`; a geometric bound controls the omitted terms. This fixes
the moderate-argument seed error without quadrature or arbitrary precision.
At `z=0` use the gamma-ratio identity; near zero use its small-`z` series.

Generate the **normalized** `j` values directly:

\[
 \log j_\rho=-\log(s\sqrt{2\pi t})-A+\rho\log p
       +\log[J_\rho(z)/\Gamma(\rho)],
\]
\[
 \log j_{\nu+1}=\log j_\nu+
                   \log\left(\frac{pR_\nu}{\nu}\right).
\]

This prevents cancellation between separately accumulated `nu log p` and
`log(J_nu/Gamma(nu))`. In the prototype that normalization improved the
worst survivor relative discrepancy from `2.7e-11` to `5.7e-14`.
Log-sum-exp handles the final positive series. The prototype uses extended
precision for cumulative log additions, ordinary doubles for the backward
ratios, and depth doubling to check the continued fraction. That depth test
is empirical; a fully certified floating-point implementation remains work.

For `q<1`, `z<0`, the first positive incomplete-gamma representation avoids
needing this positive-`z` continued fraction altogether.

The cylinder integral and its standard large-variable expansions can be
checked against [DLMF 12.5](https://dlmf.nist.gov/12.5) and
[DLMF 12.9](https://dlmf.nist.gov/12.9). The positive mixture and survivor
identities here are derived above, not attributed to those references.

## 5. A direct CDF expansion with a remainder bound

For `delta=tau-lambda>0`, expand only `exp(-u V^2/(2s^2))` in the density.
Integrate each term in time **analytically** from zero to `t`. Define
`a=d^2/(2s^2)` and

\[
 K=\frac{d\tau^k}{s\sqrt{2\pi}\delta^k},\qquad
 B_n=K\frac{(k)_{2n}}{n!(2s^2\delta^2)^n}
       a^{n-1/2}\Gamma(1/2-n,a/t)>0.
\]

Then for each finite `N`,

\[
 \boxed{F_k(t)=\sum_{n=0}^{N-1}(-1)^n B_n+R_N,
      \quad |R_N|\le B_N,\quad (-1)^N R_N\ge0.}
\]

The remainder sign and bound follow directly from Taylor's theorem for
`exp(-w)`, `w>=0`, before either positive integration. Thus this is not an
uncontrolled asymptotic assertion. If the partial sum is `C_N>B_N`, its
relative error is at most `B_N/(C_N-B_N)` in exact arithmetic.

**The infinite series generally diverges.** Stop when the bound meets the
required accuracy and reject the branch if it does not. It cannot replace
the exact positive survivor representation everywhere. This is a finite
analytic special-function calculation, with no time or rate quadrature.

Negative-order upper incomplete gamma is well defined at positive argument.
To avoid large factors, use
`a^(n-1/2) Gamma(1/2-n,A) = t^(n-1/2) E_(n+1/2)(A)` and evaluate the scaled
generalized exponential integral with a continued fraction. The current
research prototype instead uses 60-digit `mpmath.gammainc`; its observed
term counts are not timings for a compiled double-precision implementation.
An alternating-sum condition check and roundoff allowance must accompany
the analytic truncation bound.

## 6. Validation and actual remaining cost

Run `python WorkingTests/rdmswtn_gamma_positive.py`. The JSON beside it
records 25 cases on the original grid, six additional inverse-gamma
identity checks varying `d,s` and real shape, and 36 survivor stress cases
crossing real shapes `.3,3.7,100`, `q=1.0001,1.5,2,20`, and `A=.1,1,20`
with `d=1.3,s=.7`. Independent adaptive log-drift
integration is a **reference only**; the candidate functions do not call it.

Observed results on the original grid:

| Check | Maximum discrepancy |
|---|---:|
| Positive NB mixture, either probability | `1.8e-14` absolute |
| Positive moment survivor series | `5.7e-14` relative |
| Scaled density recurrence | `1.4e-14` relative |

The 36 added survivor stress cases had maximum relative discrepancy
`5.7e-14` after replacing the inaccurate moderate-argument cylinder seed.

For `d=s=1`, mean drift 2, shape 100 (`q=50`):

| t | Positive S terms | S relative discrepancy | Direct CDF bounded? |
|---:|---:|---:|---|
| .02 | 7,477 | `6.7e-15` | 8 terms; remainder bound `5.7e-16` relative |
| .2 | 1,809 | `1.7e-14` | 13 terms; bound `9.7e-15` |
| 1 | 785 | `7.6e-15` | 25 terms; bound `2.0e-14` |
| 5 | 395 | `1.1e-14` | bound not reached; use survivor |
| 100 | 149 | `5.7e-14` | bound not reached; use survivor |

In the three bounded shape-100 CDF cases, relative discrepancies against
independent integration were at most `5.4e-15`. At `t=100`, `logS` is about
`-69.7153`; the earlier series lost approximately 42 decimal digits there.
The new positive representation avoids that cancellation entirely.

These grids are feasibility checks, not a universal domain certificate.
The survivor truncation bound covers truncation conditional on correct
moments; it does not certify the special functions, continued fractions,
or the independent numerical reference.

There is still a real performance issue. The survivor expansion's terms
are concentrated roughly near total shape `nu=tau*d/t` with a width of
order `p`, so summing from zero can waste many terms. Centering and geometric
tail bounds can remove that work. But both moment seeding by an integer
shape recurrence and the necessary width can grow with shape. The
`sv -> 0` limit needs a controlled fixed-drift expansion or another uniform
asymptotic treatment. Substituting the fixed-drift law at an arbitrary
shape threshold would be an unquantified approximation.

## 7. What is solved and what remains

The mathematical point-start survivor now has an exact, cancellation-free,
quadrature-free representation for the entire positive parameter domain.
The CDF is directly positive for `q<=1`, is the exact complement for `q>=1`,
and has a separately bounded direct expansion in part of the small-CDF
domain. This removes the specific overflow/cancellation obstruction in the
existing derivation tests without changing the model.

It does **not** yet provide a complete, uniformly fast double-precision
`log f/log F/log S` kernel. Remaining gates are a direct CDF branch wherever
the bounded expansion is insufficient, relative-tail stopping for the NB
branch, certified/scaled special-function recurrences, and compiled speed
measurements including the constant-drift limit. For example, `k=2.5`,
`q=1.25`, `t=.02` is outside the successful bounded-CDF cases: a tiny CDF
there should not simply be obtained from the otherwise accurate survivor.

Uniform distance still requires a separate derivation: `q` and `A` change
with distance, so these formulas cannot be averaged by reusing fixed
negative-binomial weights. I have not established that uniform start either
requires quadrature mathematically or is weakly identified statistically.
Those are separate questions; neither is needed for this point-start result.
