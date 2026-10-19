# Gamma drift variability for RDMSWTN: feasibility investigation

2026-10-19. Research only; no production model, dispatch, or parameter mapping
has been changed. The implementation requirement is **no quadrature in the
likelihood**, because the existing RSUB quadrature is already too costly in
practice. Quadrature below is an independent validation reference only.

**Follow-up:** [Positive analytic representations](rdmswtn_gamma_positive_derivation.md)
derive complementary positive series covering all positive shapes/rates,
a scaled moment evaluator that avoids the failures below, and a separately
bounded direct CDF expansion. These supersede the direct alternating survivor
series as the main candidate; a universally fast kernel is still unproven.

## Assessment

The fixed-distance density has a compact exact expression in a parabolic
cylinder function. This is considerably simpler than implementing the
Laguerre-function expression in Steingroever et al. directly. The practical
obstacle is the CDF and especially log survivor needed by every losing racer.

There is an exact, convergent, quadrature-free survivor series with a short
coefficient recurrence (derived below). Direct summation is fast and accurate
in part of the domain, but catastrophically ill-conditioned elsewhere. This
establishes an analytic avenue, **not** a production-ready universal kernel.
No impossibility of a fast analytic implementation is claimed. Nor is there
evidence yet that this option can match the current Gaussian kernel's speed.

The appropriate next investigation is a stable, regime-dependent analytic
CDF/survivor evaluator, with accuracy and speed gates before wiring a new
model option. Simply adding a gamma quadrature would not meet the requirement.

## What the paper establishes

Steingroever, Wabersich & Wagenmakers (2021; published online 2020),
*Modeling across-trial variability in the Wald drift rate parameter*, derive
both mixtures. Their gamma analytic expression uses generalized Laguerre
functions. They implement the gamma model with trial-specific latent drift
draws in JAGS rather than that analytic expression. This is an implementation
choice, not evidence that numerical quadrature is mathematically necessary.

Source: https://link.springer.com/article/10.3758/s13428-020-01448-7
(SW-GAM mixture section and appendix).

The formulas and recurrence below are this investigation's derivations from
the mixture integral, not formulas attributed to that paper.

## Relation to the actual repository

* `src/model_FRQ.h`: FRQ uses regularized incomplete beta/gamma functions and
  beta-prime rate mixing. Its current kernel explicitly needs no Bessel
  function or quadrature.
* `src/model_OPR.h`: fixed-rate OPR uses modified Bessel I, positive Poisson/
  incomplete-gamma series, and asymptotic expansions. Its gamma-rate branch
  uses positive beta-prime series. This is a good standard for the numerical
  properties we want, but those recurrences do not directly evaluate this law.
* `src/model_RSUB.h`: packet models have Bessel/Poisson structure; the IG
  gamma-mixture branch (`rsub::igg`) uses localized Gamma-times-Gaussian
  quadrature. Its envelope is mathematically very close to this problem,
  but reusing its quadrature would reproduce the performance problem.
* `src/model_RDM.cpp`, `drdmswtn`/`prdmswtn`: the current no-clock Gaussian
  density/CDF/survivor use `rdmswtn_spv.h` when appropriate.
* `src/rdmswtn_spv.h`: the Gaussian distance recurrence is already in
  production. The older `WorkingTests/rdmswtn_analytic_derivation.md` is
  historical and should not be read as current dispatch documentation.
* `src/wald_functions.h`: fixed-drift uniform-distance Wald functions are
  available, but averaging them over gamma drift by quadrature is not the
  proposed production solution.

The key structural difference is that RDMSWTN keeps diffusion SD fixed when
drift varies. It is not a pure time-scale mixture: changing drift changes
both conditional mean and conditional relative variability. For example,
`Var(T | V)/E(T | V)^2 = s^2/(d V)`. In RSUB's gamma mixtures the trial rate
is a pure time scale. Thus beta-prime conjugacy does not transfer directly.

## Exact point-distance density

Let decision time be `t = RT - t0 > 0`, distance `d > 0`, diffusion SD `s > 0`,
and trial drift `V ~ Gamma(k, rate tau)`. Preserve the mean/SD coordinates by

\[
k=(v/s_v)^2,\qquad \tau=v/s_v^2.
\]

Here `v` and `sv` are the actual gamma mean and SD. In the existing truncated
normal model they are the mean and SD of the underlying normal, before
truncation; equal coordinate values therefore do not imply matched actual
drift moments across the two families.

The conditional density is

\[
f(t\mid V,d)=\frac{d}{s\sqrt{2\pi t^3}}
 \exp\!\left[-\frac{(d-Vt)^2}{2s^2t}\right].
\]

Set `z = (tau - d/s^2) s / sqrt(t)`. Integrating drift gives

\[
\boxed{f(t)=\frac{d}{s\sqrt{2\pi t^3}}
 \left(\frac{\tau s}{\sqrt t}\right)^k
 \exp\!\left[-\frac{d^2}{2s^2t}+\frac{z^2}{4}\right]D_{-k}(z).}
\]

This follows from

\[
J_k(z)=\int_0^\infty y^{k-1}e^{-y^2/2-zy}\,dy
       =\Gamma(k)e^{z^2/4}D_{-k}(z).
\]

NIST's integral representation uses `U(k-1/2,z) = D_{-k}(z)`:
https://dlmf.nist.gov/12.5.E1 .

The exponential integrand is quadratic plus linear in drift. The usual
modified-Bessel-K integral instead has a linear and reciprocal term. Bessel
relations at special orders or alternative transforms do not provide a
general-purpose incomplete-Bessel evaluator for this density and its CDF.

Numerically, the object to evaluate is the **log/scaled** cylinder expression,
including its prefactor. A raw `D_{-k}(z)` call can underflow while the final
density is ordinary-sized. In the 20-case check, scipy's raw cylinder function
returned zero in three cases at shape 100. Arbitrary-precision evaluation of
the combined expression agreed with independent quadrature to relative
error `8.0e-15`. This verifies the identity, not a fast double-precision kernel.

## An exact survivor series without quadrature

Define

\[
x=t^{-1/2},\quad a=\frac{d^2}{2s^2},\quad
b=\tau s-\frac d s,\quad
C=\frac{d\tau^k s^{k-1}}{\sqrt{2\pi}\Gamma(k)},
\]

and write the entire function

\[
H(x)=e^{-ax^2}J_k(bx)=\sum_{n\ge0}h_n x^n.
\]

Then `f(t) = C t^{-(k+3)/2} H(t^-1/2)`. Since the positive-drift model is
proper, integrating this density from `t` to infinity gives the exact series

\[
\boxed{S(t)=2C x^{k+1}\sum_{n\ge0}\frac{h_n x^n}{k+n+1}.}
\]

It converges at every finite positive `t`: H is entire and its series can be
integrated termwise on the compact interval `[0,x]`. It is more than a formal
large-time asymptotic expansion. Convergence alone does not ensure stable
floating-point summation.

Integration by parts gives `J_k'' - z J_k' - k J_k = 0`. Consequently,

\[
H''+(4a-b^2)xH'+[2a-kb^2+(4a^2-2ab^2)x^2]H=0,
\]

which yields the short recurrence

\[
h_{n+2}=-\frac{[(4a-b^2)n+2a-kb^2]h_n
 +(4a^2-2ab^2)h_{n-2}}{(n+2)(n+1)},
\]

with negative-index coefficients zero and

\[
h_0=2^{k/2-1}\Gamma(k/2),\qquad
h_1=-b\,2^{(k-1)/2}\Gamma((k+1)/2).
\]

Production would normalize by h0 and recur directly on `h_n x^n` to avoid
unnecessary powers and gamma overflow. The prototype does this. Apart from
initial gamma/log-gamma evaluations, each term costs constant arithmetic.
This has the sort of inexpensive analytic structure requested.

### Numerical limitation of direct summation

`rdmswtn_gamma_series.py` compares the recurrence against independent
log-rate integration at 25 combinations of shape and time, with `d=s=1`,
mean drift 2. Illustrative double-precision results:

| Shape k | Time | Terms including conservative stopping margin | Relative survivor error |
|---:|---:|---:|---:|
| 1 | 1 | 48 | about machine precision |
| 2.5 | 1 | 50 | about machine precision |
| 10 | 5 | 75 | 9.1e-12 |
| 10 | 1 | 134 | 5.7e-6 |
| 10 | 0.2 | 291 | 7.9e12; invalid probability |
| 100 | 1 | overflow | failure |

At `k=10, t=.2`, the correctly evaluated sum loses about 29 decimal digits
through cancellation; 100-digit arithmetic restores agreement with the
reference. At `k=100,t=100`, about 42 digits cancel. At early times even a
nearly accurate S would still be unsuitable for obtaining tiny F by
subtraction. A separate direct log-CDF treatment is necessary.

The prototype's small-term stopping rule is deliberately **not** a certified
error bound. It reports false convergence in cancellation regimes. The term
counts and failures must not be presented as a production algorithm.

Promising analytic work remains: local series continuation with stability
analysis, separate leading-edge and survivor representations, scaled cylinder
recurrences/continued fractions, and a controlled near-constant-drift branch.
A local Taylor integrator for the density still needs an accurate CDF seed;
it does not solve the integration-constant problem by itself. Neither scaling
alone nor compensated summation can recover dozens of cancelled digits.
These are development directions, not validated fixes.

## Closed-form special cases

For exponential drift (`k=1`) the point-distance CDF itself reduces to normal
probabilities and an exponential/Gaussian integral. Let

\[
c=s\sqrt t,\quad P=\Phi(-d/c),\quad
J=E[\phi((d-Vt)/c)/c],\quad \delta=\tau-2d/s^2.
\]

J is just `t f(t)/d`, hence available from the cylinder expression (at k=1,
it reduces to an erfc/scaled-erfc expression). Then

\[
T=P+tJ/\tau,\qquad
R=(\tau P-tJ)/\delta,\qquad F=T+R.
\]

At delta=0 use its removable limit, not the displayed quotient. There are
additional cancellation issues near that point and in tiny survivors.
Four arbitrary-precision evaluations agree with independent integration to
relative CDF error at most `3.6e-15`.

Positive integer shapes can in principle be generated by differentiating the
unnormalized exponential-drift integral with respect to tau:

\[
F_n=\frac{\tau^n}{(n-1)!}
 (-\partial_\tau)^{n-1}(F_1/\tau).
\]

This suggests an optional fixed-shape Erlang family. It does **not** give a
continuously free real shape without a further method, and its numerical
recurrences have not been implemented or benchmarked here. These results
are for point distance, not a solved analytic uniform-start family.

## Tail and model interpretation

For fixed positive s and shape k,

\[
f(t)\sim C h_0\,t^{-(k+3)/2},\qquad
S(t)\sim\frac{2C h_0}{k+1}\,t^{-(k+1)/2}.
\]

This differs from the pure time-scale gamma mixtures, whose survivor power
is k. It also supplies an inexpensive far-tail branch once remainder control
is developed. Numerical checks at times 1e4 and 1e6 approach both constants.
The limit here holds at fixed k and s>0; it is not uniform in k or s->0.

The accumulator is proper for every k>0. Its mean decision time exists for
k>1 and equals `d tau/(k-1)`. Its second moment exists for k>3, because
`E[T^2 | V] = d^2/V^2 + d s^2/V^3`. Thus finite variance requires more than
the k>2 condition for pure time-scale mixtures. A gamma law need not be
restricted to finite-moment shapes, but its interpretation should be clear.

The gamma option therefore changes more than positive support: it controls
near-zero drift mass and the RT tail exponent. A positive-truncated normal
has a nonzero density at zero, corresponding to the k=1 tail exponent here.

## Integration scope and decision gates

An eventual interface could be `RDMSWTN(drift_dist=c("normal","gamma"))`,
with normal as the default, existing names and Gaussian behavior preserved,
and the existing fixed-drift branch retained at sv=0. Gamma with sv>0 needs
v>0; `posdrift=FALSE` does not acquire a meaningful gamma analogue and should
be rejected explicitly for that combination.

Do not automatically route correlated Gaussian drift draws to gamma.
Finishing-time copulas can conceptually use new marginals; gamma drift
correlation requires a separately specified joint law. Neither should be in
the first performance feasibility target.

Uniform start distance A>0 remains additional work. The Gaussian distance
recurrence cannot just be reused unchanged. A gamma-specific recurrence or
antiderivative would need derivation and validation; otherwise an outer
distance quadrature would defeat the stated requirement. Kill clocks also
need separate cumulative treatment. Independent guess clocks combine
algebraically once f and S are available, but do not solve their evaluation.

The next concrete gate should be a compiled, **quadrature-free point-distance**
kernel for log f, log F and log S, with controlled accuracy across real k and
the sv->0 limit. Benchmark it against today's `rdmswtn_spv` kernel and a
representative race likelihood, including parameter-varying workloads. Only
then extend A>0 and clocks. No speed claim is made from the Python prototypes.

## Reproducible checks and limits

Run from the repository root:

```
python WorkingTests/rdmswtn_gamma_feasibility.py
python WorkingTests/rdmswtn_gamma_series.py
```

Dependencies: Python, numpy, scipy, mpmath. The JSON files hold actual results.
The first script verifies the density identity, checks F+S (maximum absolute
error `1.1e-14`), checks the time derivative using the smaller probability
tail (maximum relative discrepancy `6.9e-7`, finite-difference step `1e-4 t`),
and checks convergence toward constant drift at shapes up to 1e6. It also
demonstrates why fixed prior Gauss-Laguerre nodes are not sufficient: worst
relative errors on its selected grid were 32%, 15%, and 1.6% with 32, 64,
and 128 nodes. This quadrature is a diagnostic, not the recommendation.

The second script records the convergent survivor recurrence's successes
and failures and checks the exponential-drift CDF identity. It intentionally
does not assert that the failing candidate is acceptable. All grids are
small feasibility experiments, not a production domain or error certificate.
