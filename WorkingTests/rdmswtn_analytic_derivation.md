# Analytic threshold integration for RDMSWTN

Research result, 2026-09-21. Production kernels have not been changed.

The useful result is an exact convergent analytic series for the CDF with both
`sv > 0` and `A > 0`. Its coefficients obey a small linear recurrence. This
eliminates the outer threshold quadrature. A backwards evaluation is substantially
more stable than a single midpoint Taylor expansion. This is **not** a claimed
finite expression in elementary functions or ordinary normal CDFs.

## What is currently expensive

In `src/model_RDM.h`, `drdmswtn_joint_A_sv_density_postrunc` already analytically
integrates the density over both drift and distance. The no-clock CDF helper
`prdmswtn_joint_A_sv_cdf_postrunc` instead averages 20 fixed-distance CDFs. Each
calls two bivariate normal CDFs. The untruncated CDF in `src/model_RDM.cpp`
instead averages fixed-drift, uniform-distance CDFs over drift quantiles.

The important target is therefore the **CDF**, which enters every losing
accumulator's survivor, rather than another derivation of the density.

## Definitions

Write decision time as t = RT - t0 > 0, drift mean as mu, drift SD as sigma = sv,
diffusion SD as s, and distance r in [b-A,b], with b=B+A. For positive drifts,
let p(v) be the *untruncated* N(mu,sigma^2) density and Z=Phi(mu/sigma).
Define c=s sqrt(t), Q=s^2+sigma^2 t, and S=sqrt(t Q).

The conditional Wald CDF is

    F(t | v,r) = Phi((vt-r)/c) + exp(2rv/s^2) Phi(-(r+vt)/c).

Define unnormalised terms

    T(r) = integral_0^infinity p(v) Phi((vt-r)/c) dv,
    H(r) = integral_0^infinity p(v) exp(2rv/s^2) Phi(-(r+vt)/c) dv.

Then the desired answer is

    F_A(t) = 1/(A Z) integral_(b-A)^b [T(r)+H(r)] dr.

## Several approaches examined

1. **Integrate drift first / Gaussian geometry.** This gives the existing
   fixed-distance bivariate-normal expression. T has a tractable antiderivative
   through Gaussian rectangle moments. H contains an exponential of a positive
   quadratic times a moving bivariate-normal probability. A Gaussian-CDF
   substitution does not automatically integrate that product.

2. **Integrate distance first.** The fixed-drift primitive introduces 1/v times
   differences of normal probabilities. The apparent singularity at v=0 is
   removable in the combined expression, but averaging the separated pieces over
   normal drift does not produce the desired simple Gaussian moment. Dropping
   the neighbourhood of zero would change the model.

3. **Special functions / integrating factor.** The ODE below has a quadratic
   exponential integrating factor. Its subsequent antiderivative naturally
   involves erfi/Dawson-type terms and products with normal probabilities. Merely
   naming a new integral or permitting complex Gaussian parameters does not
   establish a fast, stable implementation. No impossibility theorem for a
   finite special-function representation is asserted here.

4. **Differentiate in distance, then integrate the resulting series.** This
   closes exactly in six functions with polynomial coefficients. It yields an
   inexpensive analytic recurrence, including a useful stability direction.

## The key identity

Let g(v,r)=exp(2rv/s^2) Phi(-(r+vt)/c) and
j(v,r)=phi((r-vt)/c)/c. Direct differentiation gives

    g_r = (2v/s^2) g - j,
    g_v = (2r/s^2) g - t j.

Both identities were checked symbolically with SymPy. Since
p'(v)=-(v-mu)p(v)/sigma^2, integration by parts on [0,infinity) gives

    integral v p(v) g(v,r) dv
      = mu H + sigma^2 integral p(v) g_v(v,r) dv
        + sigma^2 p(0) Phi(-r/c).

The boundary term is essential: omitting it changes positive-truncated drift
into a different model. Consequently

    T' = -J,
    H' = (a + beta r) H - k J + d L,

where

    a = 2mu/s^2,                beta = 4sigma^2/s^4,
    k = 1 + 2sigma^2 t/s^2,    d = 2sigma phi(mu/sigma)/s^2,
    L(r) = Phi(-r/c),
    J(r) = phi((r-mu t)/S)/S * Phi(w(r)),
    w(r) = (sigma^2 r + mu s^2)/(s sigma sqrt(Q)).

Thus H has a first-order linear differential equation driven by simple
Gaussian functions, rather than needing independent bivariate-normal evaluations
at each distance.

## Closing the recurrence without convolutions

Put G=phi((r-mu t)/S)/S, K=G phi(w), E=phi(r/c), and
u=w'=sigma/(s sqrt(Q)). Then

    J' = -(r-mu t)/S^2 J + u K,
    K' = [-(r-mu t)/S^2 - u w(r)] K,
    L' = -E/c,
    E' = -r E/c^2.

Together with T and H these form y'=(M0+r M1)y for a six-component vector.
All coefficients are at most affine in r. No symbolic differentiation of
high-order bivariate normal CDFs is needed; no coefficient convolution is needed.

For r=r0+h x, write y(r0+h x)=sum y_n x^n. If

    dy/dx = (C0+x C1)y,

then

    y_(n+1) = [C0 y_n + C1 y_(n-1)]/(n+1),   y_(-1)=0.

The exact average on this interval is

    1/Z sum_(n=0)^infinity (T_n+H_n)/(n+1).

For a midpoint interval [-1,1], only even powers contribute, giving
sum (T_(2n)+H_(2n))/(2n+1)/Z. The midpoint version is fast but less stable.

The linear system has entire coefficients and entire solutions for finite
initial values. The power series converges uniformly on every compact interval,
so termwise integration is justified. An infinite series is an exact analytic
representation; a fixed 24-term implementation is a numerical approximation to
that representation, with a truncation error still requiring control.

## Evaluate backwards

Start at r=b and take h<0, propagating towards b-A. The homogeneous error in H
between two distances is multiplied by

    exp[a (r-r0) + beta (r^2-r0^2)/2].

For mu>=0 and r>=0 this is at most one when moving backwards. A forward or
midpoint propagation can instead amplify the initial special-function error.
For negative mu the coefficient can change sign; the same stability assertion
must not be applied without checking it.

The prototype limits the dimensionless step coefficients, takes 24 terms per
step, analytically sums the integral, and propagates T and H to the next step.
J,K,L,E are recomputed cheaply at each step. There is only **one** expensive
initialisation for the whole distance interval, even when it needs many steps.

## Initial values and their numerical trap

Let rho=sigma sqrt(t/Q), mu_p=mu+2r sigma^2/s^2, h1=(mu t-r)/S,
h2=(-mu_p t-r)/S, and ell=2r mu/s^2+2r^2 sigma^2/s^4. Then

    T(r) = Phi2(mu/sigma,h1;rho),
    H(r) = exp(ell) Phi2(mu_p/sigma,h2;-rho).

These are exact. Directly calculating H can lose accuracy through an extremely
small bivariate probability multiplied by a large exponential. For mu>=0,
the alternative

    H = exp(ell + log Phi(h2))
        - exp(ell) Phi2(-mu_p/sigma,h2;rho)

avoids the worst anticorrelation cancellation. The prototype additionally uses
a tail fallback when ell>20 or rho>0.9. This fallback computes **one** smooth
posterior-normal expectation:

    m* = (sigma^2 r+mu s^2)/Q,    tau* = s sigma/sqrt(Q),
    H = c G E[ 1_(V*>0) R((r+V*t)/c) ],
    V* ~ N(m*,tau*^2),           R(x)=Phi(-x)/phi(x).

R is evaluated using erfc or a large-x Mills expansion. The expectation uses
64-point GL over standardised posterior coordinates bounded by +/-12. It is
numerical quadrature, explicitly **not** a new closed form. Its omitted normal
tail is negligible in the tested positive-mu domain; the finite-node error is
empirically validated, not rigorously certified. This fallback removes the
overflow/cancellation seen in the initial prototype, without reintroducing
threshold quadrature. Bivariate normal functions themselves also use numerical
special-function algorithms internally.

## Untruncated drift and optional clocks

For `posdrift=FALSE`, replace the v integration range by the real line, set
Z=1 and d=0, and use J=G, K=0. The same recurrence applies. Initial values
require only ordinary normal CDFs:

    T=Phi(h1), H=exp(ell+log Phi(h2)).

The resulting CDF remains defective where appropriate. It must not be
renormalised by eventual hit probability.

This derivation addresses the no-kill decision CDF. Independent guess clocks
can reuse it algebraically: S_response=S_decision S_guess. A kill clock
multiplies the response density by S_kill at raw time, for every Erlang shape;
that multiplication commutes with averaging over drift and distance. A killed
CDF requires a time integral of this product and is **not** solved by the
distance recurrence alone. The operational-time model can reuse the no-clock
CDF after its existing clock transformation. Correlated-drift outer integration
also remains; making each marginal cheaper does not eliminate it.

## Validation and performance

Run from the repository root:

    Rscript WorkingTests/rdmswtn_analytic_series.R

This machine required a temporary R Makevars override for its missing
unversioned gfortran linker name:

    FLIBS = /usr/lib/x86_64-linux-gnu/libgfortran.so.5 -lm -lquadmath

The driver uses seed 20260921. It checks 2,000 ordinary positive-drift cases,
2,000 untruncated cases, independent integration over drift, and a 648-case
positive-drift stress grid. The stress grid spans t=.005 to 5, mu=.05 to 5,
B=.05 to 3, A=.05 to 2, s=.5 or 1, and sv=.1 to 3. Independent references
average the existing analytic fixed-drift uniform-distance Wald CDF using R's
adaptive integrator, changing the integration order and avoiding the new
reflection seed. The driver also checks the time derivative against the
existing analytic density and independently integrated density. In 200 ordinary
cases, the backwards CDF agreed with independent drift integration to 2.70e-13
absolute error. Its central time difference (step 1e-5) agreed with independent
densities to 8.88e-10. The installed density routine differed from those
references by up to 9.75e-7, so matching that routine alone would be an
inappropriate accuracy test. That discrepancy was recorded, not patched here.

Observed stress-grid maximum absolute CDF discrepancy against independent
drift integration: 4.46e-11, with no nonfinite answers. Backwards recurrence
versus 80-node threshold quadrature using the same robust initializer:
1.19e-14. These are empirical results, not all-parameter relative-error bounds.

A representative timing run, 40,000 evaluations per repetition and median of
five repetitions, gave:

| Kernel | Seconds |
|---|---:|
| Original fixed-threshold algebra | 0.019 |
| Original 20-node threshold algebra | 0.313 |
| One midpoint expansion, 32 terms | 0.045 |
| Backwards stepping, 24 terms | 0.061 |

The safer backwards method was about **5.1x faster** than the original threshold
quadrature and **3.2x** the original single-threshold cost. The midpoint method
was about 7x faster but fails badly on broad stress intervals and is not the
recommended universal dispatch. The original algebra benchmark copies the
ordinary-range CDF calculation and its hybrid BVN dispatch; it excludes its
zero-result fallback and some wrapper/log bookkeeping. This is a kernel
benchmark, not an end-to-end likelihood or fit benchmark. Timing includes the
new tail-seed fallback where triggered. Earlier comparisons against a slower
common robust initializer gave larger speedups; those are not the appropriate
headline comparison with current production algebra.

## Production implications

Use the backwards analytic recurrence as the main candidate for replacing
outer threshold quadrature. Keep the existing A=0 and sv=0 paths. Before
production dispatch, add adaptive remainder/error checks, bounds on step count,
direct/log survival handling, tail-relative-error validation, and robust
handling of t<=0, t=infinity, and extremely small Z. Validate smoothness across
step-count and initializer branch changes. The current prototype deliberately
does not provide that complete interface. Very large dimensionless
A*sv^2*b/s^4 can require many recurrence steps, so the measured speedup is not
uniform over all parameters. A cost-based fallback or an asymptotic branch is
appropriate there.

## Background source

Steingroever, Wabersich & Wagenmakers (2021), *Modeling across-trial variability
in the Wald drift rate parameter*, establishes the truncated-normal drift
mixture density through completing the square:
https://pmc.ncbi.nlm.nih.gov/articles/PMC8219596/ . The distance-ODE recurrence
and the numerical experiments above were derived for this repository during
this investigation; that paper is background, not a citation for this result.
