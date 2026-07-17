# Plan: fast fallback and exact two-racer correlated LBA/BAwL kernel

Date: 2026-07-17  
Target branch: `bawl-correlated-race`  
Current head when revised: `ad94c480`

## Outcome

Replace the common correlated-BAwL likelihood workload with an exact
two-correlated-racer component, while retaining a substantially cheaper
shared-factor Gauss-Hermite (GH) implementation as the general fallback.

The primary target is the n-back + prospective-memory (PM) design:

- correct and error form one correlated pair;
- PM has `rho=0` and remains independent;
- `RACE=2` disables PM and `RACE=3` enables it;
- either pair member or PM may be the observed winner;
- both `calc_ll_oo()` and `calc_ll_oo_pw()` use the same evaluator;
- the all-finite hot path reads `ParamTable` columns directly.

The central abstraction is a race made of components: one correlated pair and
zero or more independent singleton racers. At time `t`:

```text
pair member w wins:
    g_pair,w(t) * product S_independent(t)

independent racer q wins (including PM):
    f_q(t) * S_pair(t) * product S_other_independent(t)

total race survivor:
    S_pair(t) * product S_independent(t)
```

This factorisation is exact. Pair rows are discovered from active nonzero
loadings, not hard-coded accumulator positions.

## Scope

In scope:

1. Correctness-preserving fallback fixes:
   - zero-loading racers do not increase the positivity dimension;
   - an independent PM does not activate the 64-node denominator path;
   - correlated hot-path log-CDF calls use `pnorm_log_direct()`;
   - positive-drift node normalisers are cancelled algebraically.
2. A prepared, fused no-clock fallback kernel:
   - LBA/leak geometry is prepared once per row and particle;
   - only conditional means and normal endpoints vary over GH nodes;
   - no full `fast_vq` staging or repeated mask scans.
3. An exact pair kernel for:
   - standard LBA (`k=0`);
   - leak-only BAwL (`k>=0`, no active kill/guess clocks);
   - positive and unrestricted Gaussian drifts;
   - point starts (`A=0`) through dedicated exact branches;
   - unequal pair loading magnitudes;
   - arbitrary additional independent racers;
   - `RACE` masks, finite responses, race survivors, truncation,
     censoring, omissions, and leaky defective tails.
4. A direct one-dimensional pair integration (`numeric_pair`) used as the
   fallback for numerically unstable rectangles and as the authoritative
   test oracle for the exact kernel.
5. A cheaper analytic GH centre for loaded-winner integrands with selective
   scan fallback.

Out of scope for the first exact implementation:

- more than two mutually correlated active racers;
- multiple non-independent pairs in one trial;
- exact active local/global kill or guess clocks;
- changes to the public `rho` parameterisation or simulator semantics.

Those cases use the improved fallback. Routing is per trial: one unsupported
trial must not force all other trials in a particle off the exact path.

---

## Design rule: build shared final pieces once

The implementation should be organised around the following final units. A
work package must extend these units rather than add a temporary alternative
that a later package deletes.

### A. `BAwLTimeGeometry`

One authoritative fixed-time LBA/leak geometry object, shared by:

- the prepared GH node evaluator;
- exact pair survivor and cause-density formulas;
- exact component GSL integrands for known-winner censored intervals;
- existing scalar LBA/BAwL PDF/CDF kernels, but only in the final cleanup
  step. Until then the ordinary scalar kernels stay untouched: they are the
  bit-for-bit `rho=0` reference the new routes are validated against, and
  refactoring them before the new likelihood is proven adds a large
  regression surface for no immediate gain. `BAwLTimeGeometry` must still be
  built on the same authoritative helpers and parity-tested against the
  scalar kernels (T1) from the start.

It owns quantities independent of the drift mean and SD:

```text
tau = t-t0
b = B+A
exp(-k*tau), 1-exp(-k*tau)
C1, C2
drift bounds L and U
piecewise survivor coefficients alpha and beta
cause-density coefficients gamma0 and gamma1
valid/not-started/point-start status.
```

Use `bawl_leak_factors()` and the exact `k -> 0` limit. There must be no second
implementation of `exp(-k*tau)` or the small-`k*tau` branch in correlated code.

### B. `BAwLConditionalGeometry`

A prepared view for a row conditional on a particle's rho, used only by the GH
fallback. It composes `BAwLTimeGeometry` with:

```text
factor slope
residual sv and inverse residual sv
standardised threshold offset
m, span, Jacobian and stable scale terms
natural/log endpoint eligibility.
```

At node `z`, only

```text
v_q = v + slope*z
c = c_threshold + v_q/residual_sv
c_hi = c+span
```

and the corresponding normal terms are evaluated.

### C. `BAwLCorrTrialLayout` and `BAwLCorrSharedState`

One classifier and one shared state are used by ordinary, exact, and fallback
routes. They resolve once:

- `lM` role mapping;
- active rows and `RACE_mask`;
- winner rows and unique-trial/expansion maps;
- parameter column indices and direct `ParamTable` pointers;
- LT/UT/LC/UC and finite/nonfinite partitions;
- contaminant column;
- reusable per-particle masks, routes, pair indices, and output buffers.

The per-particle layout classifies active rows as loaded or independent and
produces one route:

```text
ordinary_zero_or_one_loaded
exact_pair
numeric_pair
gh_no_clock
gh_generic_clock
invalid
```

`numeric_pair` is the direct one-dimensional pair integration, selected only
when the exact rectangle route reports `unstable`; it should be rare in
production but keeps a pair trial off the general GH path.

This classifier is also the sole authority for positivity dimension. Do not
write a separate loaded-row count inside the denominator code.

Numerator and denominator routes are independent invariants of the layout:

```text
numerator route:   ordinary | exact_pair | numeric_pair | gh_prepared | generic_clock
positivity route:  none | univariate | bivariate | gh
```

A quadrature positivity denominator must not alter the numerator node
schedule or disable numerator reuse. This was the root of the earlier
64-node coupling and must hold as an explicit invariant, not an accident of
the current code shape.

### D. Probability and moment primitives

One numerical layer provides:

- `pnorm_std()` for guarded natural central calculations;
- `pnorm_log_direct()` for log/tail calculations;
- incomplete univariate normal moments through order three;
- bivariate rectangle probability, first moments, and cross moment,
  assembled from a shared per-time-point boundary grid of BVN corner values
  rather than independent per-rectangle evaluations;
- compensated/signed combinations and explicit stability status.

The analytic GH centre, exact pair kernel, and `numeric_pair` conditional
expectations all consume this layer.

### E. Component evaluator

One per-unique-trial evaluator returns pair/singleton component likelihoods.
Both total and pointwise APIs consume its output. Contaminants, flooring, and
expansion happen once after the complete per-trial likelihood is assembled.

The normalisation contract is explicit: internal component methods return the
unnormalised log value together with the component's own log drift
normaliser,

```cpp
struct LogComponentValue {
  double log_unnormalised;
  double log_drift_normaliser;  // log D_pair, log q, or 0 (unrestricted)
  ComponentStatus status;
};
```

and per-trial assembly subtracts each normaliser exactly once. This prevents
double division by `D_pair` (or by a singleton `q`) in truncation, endpoint
survivors, and PM-winner assembly, where the same component value appears in
several expressions.

---

## Current issues to capture in the baseline

The current correlated likelihood in `src/particle_ll.cpp`:

- materialises a `NumericMatrix` per particle;
- recomputes time/leak geometry inside PDF and CDF callbacks at every node;
- performs a 10+10 node scan/refinement through `|rho|<=.6`, 12+12 above;
- calls `R::pnorm(..., log=TRUE)` directly at four positivity sites;
- counts every active racer in the closed-form positivity dimension.

The last point is especially important. The denominator loop increments
`n_active` for an exact `rho=0` PM row. Correct/error + independent PM is thus
misclassified as three-dimensional and activates the 64-node denominator
sweep, even though the numerator already factors PM out of the latent node
loop.

Record timings and route/node counts before changing this behaviour.

---

## Mathematical specification

### 1. Fixed-time piecewise geometry

Let racer `i` have absolute evaluation time `t` and decision time

```text
tau_i = t-t0_i.
```

Different racers may have different `t0`, `A`, `B`, `sv`, and `k`. If
`tau_i<=0`, its survivor is one and cause density zero.

For `tau_i>0`, define the drift interval associated with starts in `[0,A_i]`.

For LBA (`k_i=0`, `b_i=B_i+A_i`):

```text
L_i = B_i/tau_i
U_i = (B_i+A_i)/tau_i.
```

For leaky BAwL (`k_i>0`):

```text
E_i  = exp(-k_i*tau_i)
G_i  = 1-E_i
C1_i = k_i*b_i/G_i
C2_i = k_i*E_i/G_i
L_i  = C1_i-C2_i*A_i
U_i  = C1_i.
```

Conditional on drift `V_i`, the no-hit probability is piecewise affine:

```text
s_i(V_i,t) = 1                         for V_i<=L_i
             alpha_i+beta_i*V_i       for L_i<V_i<U_i
             0                         for V_i>=U_i.
```

Coefficients:

```text
LBA:   alpha_i = b_i/A_i
       beta_i  = -tau_i/A_i

BAwL:  alpha_i = C1_i/(C2_i*A_i)
       beta_i  = -1/(C2_i*A_i).
```

The winner cause-density weight over `[L_i,U_i]` is affine:

```text
h_i(V_i,t) = gamma0_i+gamma1_i*V_i

LBA:   gamma0_i = 0
       gamma1_i = 1/A_i

BAwL:  gamma0_i = -exp(k_i*tau_i)*k_i*b_i/A_i
       gamma1_i =  exp(k_i*tau_i)/A_i.
```

The BAwL expression follows from
`|da(V,t)/dt|=exp(k*t)*(V-k*b)` and has the LBA expression as its limit.

### 2. Bivariate rectangle moments

For correlated underlying drifts `(V_1,V_2)` and rectangle
`R=[a_1,b_1] x [a_2,b_2]`, define:

```text
P(R)   = E[1(V in R)]
M1(R)  = E[V_1 1(V in R)]
M2(R)  = E[V_2 1(V in R)]
M12(R) = E[V_1 V_2 1(V in R)].
```

These suffice because every fixed-time integrand contains at most one affine
factor in each drift.

For unrestricted drifts the domain lower bound is `-Inf`; for jointly positive
drifts it is zero. Under positive drifts divide pair quantities by:

```text
D_pair = P(V_1>0,V_2>0).
```

Independent positive-drift racers continue through ordinary individually
conditioned BAwL kernels. Since they are independent of the pair, the full
normaliser factorises exactly.

### 3. Exact pair survivor

Let `d_i=0` for positive drifts and `-Inf` otherwise. Split each axis into:

```text
C_i = [d_i,L_i]     # survivor factor 1
Q_i = [L_i,U_i]     # survivor factor alpha_i+beta_i V_i.
```

The unnormalised survivor is:

```text
S*_pair(t) =
    P(C_1 x C_2)
  + alpha_1 P(Q_1 x C_2) + beta_1 M1(Q_1 x C_2)
  + alpha_2 P(C_1 x Q_2) + beta_2 M2(C_1 x Q_2)
  + alpha_1 alpha_2 P(Q_1 x Q_2)
  + beta_1 alpha_2 M1(Q_1 x Q_2)
  + alpha_1 beta_2 M2(Q_1 x Q_2)
  + beta_1 beta_2 M12(Q_1 x Q_2).
```

Then `S_pair=S*_pair/D_pair` for positive drifts. This value is required when
an independent PM racer wins.

### 4. Exact pair cause density

If racer 1 wins, define:

```text
W_1 = [L_1,U_1]
C_2 = [d_2,L_2]
Q_2 = [L_2,U_2].
```

The unnormalised cause density is:

```text
g*_pair,1(t) =
    gamma0_1 P(W_1 x C_2)
  + gamma1_1 M1(W_1 x C_2)
  + gamma0_1 alpha_2 P(W_1 x Q_2)
  + gamma1_1 alpha_2 M1(W_1 x Q_2)
  + gamma0_1 beta_2 M2(W_1 x Q_2)
  + gamma1_1 beta_2 M12(W_1 x Q_2).
```

Swap indices for racer 2. Divide by `D_pair` under positive drifts.

For LBA this reduces to:

```text
g*_pair,1(t) =
    M1(W_1 x C_2)/A_1
  + [b_2 M1(W_1 x Q_2)-tau_2 M12(W_1 x Q_2)]/(A_1 A_2).
```

### 5. Rectangle moment computation

Implement one allocation-free primitive:

```cpp
struct BvnRectMoments {
  double p, m1, m2, m12;
  enum Status : uint8_t { ok, zero_mass, unstable, invalid } status;
};
```

Use derivatives of rectangle probability with respect to the mean. If
`P`, `g=dP/dmu`, and `H=d2P/dmu dmu'`, then:

```text
E[V 1_R] = mu P + Sigma g
E[(V-mu)(V-mu)' 1_R] = P Sigma + Sigma H Sigma.
```

For a BVN CDF corner:

```text
d/dx Phi2 = phi(x) Phi((y-rho*x)/sqrt(1-rho^2))
d/dy Phi2 = phi(y) Phi((x-rho*y)/sqrt(1-rho^2))
d2/dxdy   = phi2(x,y;rho).
```

Use corresponding analytic same-axis derivatives, the standardisation chain
rule, and signed four-corner sums. Infinite-bound derivatives are zero.

The dominant cost of the exact kernel is BVN CDF evaluation, and the
rectangles within one pair evaluation share boundaries. For pair survival the
boundaries form the `3x3` grid `{d_1,L_1,U_1} x {d_2,L_2,U_2}`; for a pair
cause a `2x3` grid `{L_w,U_w} x {d_l,L_l,U_l}`. Compute one boundary grid of
corner CDF values and boundary derivatives per time point,

```cpp
struct BvnBoundaryGrid {
  double x[3], y[3];      // standardised axis boundaries
  double cdf[3][3];       // Phi2 corner values
  double dx[3][3];        // d/dx corner derivatives
  double dy[3][3];        // d/dy corner derivatives
  double dxy[3][3];       // phi2 corner values
  uint8_t nx, ny;
  BvnStatus status;
};
```

and assemble every rectangle probability and moment as signed combinations of
cached corners. Never evaluate the same corner twice within one pair
evaluation.

Also cache per-pair quantities that do not depend on time — `rho`,
`sqrt(1-rho^2)`, `D_pair`, `log D_pair`, and standardised positivity
boundaries — once per particle. This matters most inside outer-time
integration for censored known-winner observations, where the same pair is
evaluated at many time points.

Numerical requirements:

- reuse `norm_cdf_2d()`/`norm_ucdf_2d()` from `src/gaussian.h`;
- choose lower, upper, or complemented forms to avoid subtracting values near
  one;
- use compensated/signed accumulation;
- include stable log/scaled handling for tiny rectangles, especially
  `rt` near `t0`;
- add explicit `rho` limits near `0`, `+1`, and `-1`;
- return `unstable` rather than silently changing sign or clamping a material
  error;
- route an unstable trial to `numeric_pair` first, and to the fused GH
  fallback only if that also fails.

### 6. Degenerate starts

Point starts are handled exactly in the first implementation. In practice
`A=0` is essentially unreachable (the parameter transform bounds `A` above
zero), so this is a robustness and routing-continuity measure, not a
performance one — but the point-start formulas are simpler than the
interval-start ones, so exact support is cheaper than a GH detour. Use a
dedicated `A=0` branch; never divide by a small artificial `A`.

For a point-start winner the drift is pinned at the unique threshold-crossing
value (note `v*_w = U_w`, the interval-start upper drift bound):

```text
LBA:   v*_w(t) = b_w/tau_w          |dv*_w/dt| = b_w/tau_w^2
BAwL:  v*_w(t) = k_w*b_w/G_w        |dv*_w/dt| = k_w^2*b_w*E_w/G_w^2
```

and the unnormalised pair cause is a boundary density:

```text
g*_pair,w(t) = f_Vw(v*_w) |dv*_w/dt| E[s_l(V_l,t) 1(V_l>d_l) | V_w=v*_w].
```

The conditional loser expectation is an incomplete univariate-normal moment
calculation (piecewise-affine survivor under a conditional normal); no BVN
rectangle is required.

A point-start loser has an indicator survivor `1(V_l<=v*_l)` in drift space
(`L_l=U_l=v*_l`), so pair survivor and cause expressions degenerate to
univariate or conditional-normal terms. Support point-start winner,
point-start loser, and both-point-start in the first exact version, and test
`A -> 0` continuity into these branches (T9).

### 7. Direct one-dimensional pair integration (`numeric_pair`)

The rectangle formulas are mathematically exact, but tiny rectangles at
`rt` near `t0` with `|rho|` near one require cancellation between very small
CDF and moment terms — exactly where the latent-factor GH also struggles.
Between the exact route and the general GH fallback sits a pair-specific
numerical route that integrates over one drift only, using the fact that the
conditional expectation over the other drift is analytic.

Conditional on `V_w=v`, `V_l` is normal with known mean and SD, and
`s_l(V_l,t) 1(V_l>d_l)` is piecewise affine, so
`e_l(v,t) = E[s_l(V_l,t) 1(V_l>d_l) | V_w=v]` uses only the univariate
incomplete-moment primitives. The unnormalised pair cause is then one
finite-interval integral over the winner drift:

```text
g*_pair,w(t) = integral_{max(d_w,L_w)}^{U_w}
    (gamma0_w+gamma1_w*v) phi((v-mu_w)/sv_w)/sv_w e_l(v,t) dv,
```

and the unnormalised pair survivor splits over the constant and affine
regions of `s_w`:

```text
S*_pair(t) = integral_{d_w}^{L_w} phi_w(v) e_l(v,t) dv
           + integral_{L_w}^{U_w} phi_w(v) (alpha_w+beta_w*v) e_l(v,t) dv.
```

Divide by `D_pair` as usual. The integrands are smooth within each affine
piece, so fixed high-order Gauss-Legendre per piece (or the existing GSL
adaptive integrator for the semi-infinite piece) suffices.

This route has three uses:

1. robust fallback for rectangles reporting `unstable`;
2. authoritative test oracle for the exact kernel (T6) — dense latent-factor
   GH is not accurate enough near `|rho|=.95`, `rt` near `t0` to certify the
   closed form, since it can share the convergence problem being replaced;
3. diagnostic separating a rectangle-moment bug from a GH convergence
   problem.

---

## Final architecture

### 1. Shared state and direct column pointers

Add `BAwLCorrSharedState` outside the particle loop, parallel to
`RaceSharedState`. Suggested contents:

```cpp
struct BAwLCorrSharedState {
  bool valid;

  // R holders keeping backing memory alive
  Rcpp::NumericVector rt, LT, UT, LC, UC;
  Rcpp::IntegerVector race_nacc;
  Rcpp::LogicalVector race_mask;

  // Data-fixed arrays
  std::vector<int> active;
  std::vector<int> winner;
  std::vector<int> role_is_correct;
  std::vector<int> finite_unique_idx;
  std::vector<int> other_unique_idx;

  // Reused per-particle layout and scratch
  std::vector<BAwLCorrTrialLayout> layout;
  std::vector<int> singleton_winner_mask;
  std::vector<int> singleton_loser_mask;
  std::vector<int> isok_int;
  std::vector<BAwLConditionalGeometry> conditional_geometry;
  std::vector<double> singleton_result;
  std::vector<double> pair_result;
  std::vector<double> ll_unique;
};
```

Resolve `keep_names` into pointers to `ParamTable::base` once, as the ordinary
raw path already does. The correlated entry point consumes:

```cpp
const double* const* cols
```

and must not materialise, clone, repack, or perform string lookups on the
all-finite exact or fused no-clock paths.

The legacy generic-clock fallback may materialise only when it is actually
selected. Do not maintain separate total and pointwise correlated evaluators.

### 2. Canonical trial classification

After applying the existing `lM` role mapping, scan active rows once per
particle and trial:

```text
0 nonzero effective loadings:
    ordinary independent race

1 nonzero effective loading:
    marginally ordinary; treat all racers as independent

2 nonzero loadings, valid pair, no clocks:
    exact pair when numerically supported
    numeric_pair when the exact route reports unstable
    otherwise fused no-clock GH

>2 nonzero loadings, valid no-clock factor model:
    fused no-clock GH

active clocks:
    generic GH/race fallback
```

Rules:

- count only `isok && RACE_mask` rows;
- exact-zero PM is an active independent singleton, not a loaded dimension;
- inactive PM is absent entirely;
- pair rows may have any indices;
- two row variance shares imply direct correlation
  `sign(r1*r2)*sqrt(abs(r1*r2))`;
- unequal loading magnitudes (`|r1| != |r2|`) are fully supported by the
  exact kernel — it needs only the marginal means, marginal SDs, and the
  induced direct correlation. Do not route a numerically valid unequal pair
  off the exact path; if a model family expects equal magnitudes, record a
  construction-time diagnostic instead;
- exact `rho=0` uses the ordinary path for bit-for-bit nesting;
- use the current `1e-14` loading threshold initially and test continuity.

This layout is reused by denominator, numerator, exact, singleton, and fallback
code. No downstream routine recounts active/loaded rows.

### 3. Fused no-clock GH fallback

Prepare one `BAwLConditionalGeometry` per loaded row and particle. Independent
rows are evaluated once with ordinary singleton kernels.

Use the unnormalised joint-positive identities for loaded rows:

```text
q(z)  = Phi(mu_z/sd_z)
F0(t) = unrestricted no-clock CDF

weighted winner = f0(t)
weighted loser  = q(z)-F0(t).
```

For `posdrift=FALSE`, the same prepared endpoint evaluator instead returns:

```text
winner = f0(t)
loser  = 1-F0(t),
```

and there is no positivity denominator. Select this mode once per likelihood
route, not with a repeated model branch inside every endpoint calculation.

Truncation windows go through the total race survivor. Independent survivor
terms depend on the endpoint time, so they cannot be factored outside the
endpoint difference: with

```text
S*_all(t) = E_z[product_loaded (q_i(z)-F0_i(t|z))] * product_indep (q_j-F0_j(t))
```

the window probability is

```text
P(LT<T<=UT) = [S*_all(LT) - S*_all(UT)] / (D_loaded * product_indep q_j).
```

Every survivor factor — loaded and independent — is evaluated at both `LT`
and `UT` inside the subtraction. Only the time-independent positivity
constants `q` cancel between numerator and denominator. `log S_all(t)`
(architecture section 7) is the sole truncation interface for both the exact
and GH routes; do not write a separate loaded-group-only window formula.

The final no-clock loop should evaluate a complete trial over nodes:

```cpp
for (int j : gh_no_clock_trials) {
  for (int q = 0; q < n_nodes; ++q) {
    double node_log = node_log_weight[j,q];
    for (int row : layout[j].loaded_rows) {
      const double vq = geom[row].v + geom[row].slope*z[j,q];
      node_log += row_is_winner
        ? prepared_log_pdf_unrestricted(geom[row], vq)
        : prepared_log_survivor(geom[row], vq, drift_mode);
    }
    node_value[q] = node_log;
  }
  output[j] = stable_log_sum_exp(node_value);
}
```

This removes:

- repeated `rt-t0`, `b`, leak exponential, `m`, span, and Jacobian work;
- separate PDF and CDF callback entry;
- full winner/loser mask scans per node;
- full-length `fast_vq` and node result arrays.

Factor existing scalar endpoint algebra into shared prepared evaluators; do not
create temporary unnormalised raw callbacks that this fused loop later removes.
A scalar one-row wrapper may exist for tests and generic reuse, but it must call
the same prepared evaluator.

Benchmark prepared trial-major/node-inner against prepared node-major/row-inner
before choosing the final ordering. Both must share the same geometry and
endpoint code, so this comparison does not create two maintained kernels.

Use a fixed-capacity stack node buffer with a validated maximum-256 override
branch. Retain stable max-then-sum log reduction.

### 4. Positivity denominator and normal CDF

The positivity dimension equals `layout.loaded_rows.size()`.

- Two loaded rows retain the exact BVN denominator.
- One loaded row uses `pnorm_log_direct()`.
- More than two loaded rows use GH.
- Zero-loading rows never set `den_by_quadrature` or activate 64 nodes.
- With truncation, conditioned independent survivors remain outside the
  loaded-group denominator.

Audit all correlated positivity calls together:

- replace `R::pnorm(x,0,1,1,1)` with `pnorm_log_direct(x,true)`;
- use `pnorm_std()` only for guarded natural central values;
- test generic-clock and fused paths against authoritative R values.

### 5. Analytic GH centre

Use the final fused evaluator with one recentered pass. Do not implement this
on the legacy buffer/callback layout.

Ignoring loser tilt, the winner contribution is an affine-weighted marginal
normal restricted to `[L_w,U_w]`:

```text
(gamma0_w+gamma1_w V_w)
  Normal(V_w|mu_w,sv_w) 1(L_w<V_w<U_w).
```

Use the shared incomplete-normal moments to obtain its mean and variance. If

```text
a_w = sign(row_rho_w)*sqrt(abs(row_rho_w)),
```

then:

```text
E[z]   = a_w*(E[V_w]-mu_w)/sv_w
Var[z] = 1-abs(row_rho_w)
         + abs(row_rho_w)*Var[V_w]/sv_w^2.
```

This derivation assumes the integrand contains a loaded winner. Some GH
integrands do not: an independent/PM winner over a loaded survivor group,
unknown-winner censoring, truncation normalisers, and exact-pair failures on
PM-winner trials are survivor-only. There is no winner tilt from which to
derive the centre, so the routing is:

```text
loaded winner in integrand:
    analytic centre -> one fine pass -> scan only on failed diagnostics

survivor-only integrand:
    existing scan/recentre (or a separately derived survivor proposal later)
```

Never apply the winner-centre formula to a survivor-only integral.

Run one fine pass. Invoke the existing scan+recenter only for trials whose
already-computed nodes show an edge peak, invalid integral, invalid proposal,
or failed lower-order/error check. Masks are per trial. Compact fallback rows
only if profiling shows masked traversal remains significant.

### 6. Exact pair plus singleton assembly

Core allocation-free calls:

```cpp
bawl_corr_pair_log_survival(..., absolute_t, row1, row2, ...)
bawl_corr_pair_log_cause(..., absolute_t, winner_row, loser_row, ...)
```

They return value plus status (`ok`, `zero`, `unsupported`, `unstable`,
`invalid`). `unstable` selects `numeric_pair` for that trial;
`unsupported` (or a `numeric_pair` failure) selects fused GH.

For all-finite known responses, call `dbawl_raw()`/`pbawl_raw()` once over
singleton masks, excluding pair and GH-generic rows. Assemble:

```text
winner in pair:
    log g_pair,w(rt) + sum log S_singleton_loser(rt)

winner singleton/PM:
    log f_singleton_w(rt)
    + log S_pair(rt)
    + sum log S_other_singleton(rt).
```

Canonical coverage:

```text
RACE=2: pair active, PM absent
RACE=3, n-back winner: pair cause * PM survivor
RACE=3, PM winner: PM density * pair survivor.
```

Any additional true zero-loading racers use the same singleton assembly.

### 7. Nonfinite and truncated data

One component survivor supplies:

```text
log S_all(t) = log S_pair(t)+sum log S_singleton(t).
```

Use it for:

- `log[S_all(LT)-S_all(UT)]` truncation normalisers;
- unknown-winner censoring intervals;
- missing-RT interval unions;
- pair survival when an independent racer wins;
- omission and leaky never-finish mass.

Known-winner intervals still require:

```text
integral_low^high cause_w(t) dt.
```

Use the exact component cause as the one-dimensional GSL integrand. It reads
raw columns and row indices and builds `BAwLTimeGeometry` once per outer time
evaluation. There is no nested latent quadrature on exact trials.

For finite RT with unknown winner, log-sum all component cause densities.
Evaluate `S_pair(+Inf)` with explicit leaky limits rather than passing infinity
through finite-time expressions.

Active clocks use the generic fallback and existing clock semantics. If the
adapter resolves all rates off and sets `kill_active=false`, the no-clock route
may be used.

### 8. Finalisation

After a complete unique-trial value is assembled:

- apply `pContaminant` once using existing finite/`+Inf` rules;
- apply truncation correction once;
- floor only the final trial likelihood;
- expand/reduce for total output or write pointwise output.

No intermediate pair, singleton, rectangle, or node term is floored at
`min_ll`.

---

## File organisation

Suggested final files:

```text
src/model_LBA.h                 scalar LBA/BAwL formulas consuming shared geometry
src/bawl_geometry.h             BAwLTimeGeometry and BAwLConditionalGeometry
src/bvn_rect_moments.h          univariate/BVN probability and moment primitives
src/model_bawl_corr_exact.h     pair survivor/cause and component helpers
src/particle_ll.cpp             shared state, classification, dispatch, finalisation
src/utils.h                     retained generic scalar/raw callback contracts
tests/testthat/test-bawl-correlated.R
tests/testthat/test-bvn-rect-moments.R
tools/benchmark-bawl-corr-exact.R
```

If include constraints make a separate header undesirable, keep the same
logical separation in existing files. Do not duplicate functions to satisfy a
file layout preference.

---

## Test matrix

Tests are defined once here. Work packages refer to these IDs rather than
restating overlapping grids.

### T1. Shared geometry and endpoint parity

Compare shared/prepared evaluators with existing scalar LBA/BAwL PDF/CDF
values over:

- LBA, small positive `k`, and moderate `k`;
- central and tail endpoints;
- `rt` near `t0` and ordinary RT;
- point and interval starts;
- positive and unrestricted drifts;
- natural and forced log fallbacks.

### T2. Normal CDF audit

Compare `pnorm_std()` and `pnorm_log_direct()` with R's `pnorm` over the
production mean/SD/rho range and extreme tails. Verify all correlated generic
and fused call sites use the shared helpers.

### T3. Layout, RACE, and positivity dimension

Test:

- zero, one, two, and three genuinely loaded active rows;
- active versus inactive independent PM;
- pair rows at non-default indices;
- malformed pair magnitudes;
- `lM` sign roles and both cell-correlation signs;
- mixed routes in one particle.

Correct/error + active independent PM must report loaded dimension two and
must not activate the 64-node denominator. Giving PM nonzero rho must restore
the genuine three-loaded fallback.

### T4. Fused GH fallback

Compare:

- old positive-weighted versus new unnormalised node expressions;
- prepared versus scalar endpoints;
- prepared trial-major versus prepared node-major values;
- fused final integrals versus dense 200/256-node reference;
- analytic-centre one-pass versus scan/refine;
- forced difficult cases that must scan again;
- more than two loaded racers;
- active clocks bypassing fused no-clock code.

Cover both rho signs, positive/unrestricted drift, LBA/BAwL, RACE masks, PM,
and truncation.

### T5. Rectangle moments

Compare probability against a deterministic bivariate reference — direct
conditional-normal integration, or `mvtnorm::pmvnorm(..., algorithm =
TVPACK())`; do not use the default randomized GenzBretz algorithm — and
moments against independent conditional-normal integration for:

- rho in `{-.95,-.8,-.5,0,.5,.8,.95}`;
- central, one-tail, two-tail, narrow, and semi-infinite rectangles;
- asymmetric means/SDs;
- continuity near rho 0 and limits near +/-1;
- upper/complemented and scaled-log paths.

### T6. Exact pair scalar values

Compare pair survivor and each pair-member cause density against the direct
one-dimensional `numeric_pair` integration as the primary reference, with
dense 200/256-node GH as a secondary comparison only. Dense GH at
`rho=.95`, `rt` near `t0` is not accurate enough to certify the closed form
on its own. Grid:

```text
rho       = {-0.95,-0.8,-0.5,-0.1,0,0.1,0.5,0.8,0.95}
posdrift  = {FALSE,TRUE}
k         = {0,small positive,moderate positive}
rt-t0     = {very small,ordinary,upper tail}
v/sv      = {negative,near zero,ordinary positive}
A         = {0,ordinary} for winner, loser, and both.
```

Include asymmetric racer parameters, unequal loading magnitudes, and
degenerate-route cases. Report three error measures rather than one blanket
criterion: absolute error `|p-p_ref|`, relative error where the reference is
not tiny, and log-domain error `|log p - log p_ref|` for tail values. A
single `1e-4` log-likelihood criterion hides large relative errors on
ordinary probabilities and is unassessable after underflow. Targets:
tight relative error in the ordinary region; `<=1e-4` log error everywhere
intended for exact evaluation.

Additionally test internal probability identities, which catch sign,
Jacobian, normaliser, and assembly errors that two agreeing implementations
can share:

```text
-dS_pair/dt = g_pair,1(t)+g_pair,2(t)        (numerical derivative,
                                              away from t0 boundaries)
0 <= S_pair(t) <= 1
S_pair(t+eps) <= S_pair(t)
g_pair,i(t) >= 0
integral[g_1+g_2] = 1-S_pair(Inf)
-dS_all/dt = sum_i g_i(t)                     (full race)
```

### T7. Pair plus independent racers

Explicit canonical tests:

1. `RACE=2`: pair only.
2. `RACE=3`, pair winner: pair cause times PM survivor.
3. `RACE=3`, PM winner: PM density times pair survivor.
4. Changing PM parameters affects only enabled trials.
5. Two correlated plus two independent racers, every winner identity.
6. Independent PM does not change pair normalisation.

Compare with dense factor integration and ordinary singleton values.

### T8. Data-path coverage

For exact pair + independent PM:

- all-finite known responses;
- uniform and varying LT/UT;
- `rt=-Inf`, `rt=+Inf`, and `rt=NA`;
- known and unknown response labels;
- finite RT with missing winner;
- leaky defective upper-tail mass;
- compressed and uncompressed data;
- nonzero contaminant;
- `calc_ll_oo()==sum(calc_ll_oo_pw())`.

### T9. Nesting and routing continuity

- exact rho zero remains bit-for-bit ordinary LBA/BAwL;
- sweep rho through zero and exact/fallback thresholds;
- sweep `A -> 0` into the exact point-start branches;
- force exact->numeric_pair->GH route switches and check value continuity;
- sweep `k` through the LBA limit;
- sweep RT toward `t0` at high rho;
- inspect first and second likelihood differences for routing kinks.

### T10. Full regression

Run focused tests from a fresh temporary installation, then the full package
suite. Preserve active-clock and more-than-two-loaded behaviour.

---

## Benchmark and observability

Add `tools/benchmark-bawl-corr-exact.R` with fixed seeds and at least 20,000
trials per scenario:

```text
plain LBA, two racers
correlated LBA, two racers
correlated LBA + independent PM, mixed RACE 2/3
leak-only BAwL, two racers
leak-only BAwL + independent PM
positive and unrestricted drifts
rho .5, .8, .95
ordinary RT and RT near t0
forced no-clock GH (>2 loaded)
forced generic GH (active clocks)
prepared trial-major versus prepared node-major.
```

Report absolute time, ratio to ordinary LBA/BAwL, speedup from baseline, route
counts, geometry preparations, and node evaluations.

Environment-gated test counters:

```text
ordinary_zero_rho_trials
ordinary_single_loaded_trials
exact_pair_trials
exact_pair_pair_winner_trials
exact_pair_independent_winner_trials
exact_pair_point_start_trials
numeric_pair_trials
gh_no_clock_trials
gh_generic_clock_trials
loaded_dimension_0/1/2/3plus
prepared_rows
fused_node_evaluations
bvn_corner_evaluations
analytic_center_eligible_trials
analytic_center_success_trials
survivor_scan_trials
scan_refinement_trials.
```

Counters are compile-time disabled or non-atomic outside test/benchmark mode.

Initial performance gates:

- active independent PM never causes a 64-node denominator sweep;
- enabling PM costs approximately one ordinary singleton evaluation;
- exact canonical pair is at least 3x faster than the current correlated path;
- positive-drift fallback improves materially from normaliser cancellation;
- no-clock fallback prepares leak geometry once and shows a clear end-to-end
  win over legacy node callbacks;
- one-pass centred trials avoid approximately half the full node work;
- ordinary noncorrelated likelihoods do not regress materially.

Record actual results in this file; treat these as architecture gates, not
promises of a particular final multiplier.

---

## Progress log

### Step 1 — baseline and observability (DONE, 2026-07-17)

Added `src/bawl_corr_counters.h` (env-gated via `EMC2_BAWLCORR_COUNTERS`,
non-atomic, exported accessors `bawl_corr_counter_values()` /
`bawl_corr_counters_reset()`) and `tools/benchmark-bawl-corr-exact.R`.

Baseline (this machine, 1 core, 20,000 unique trials, 1 particle, median of
3 calls, `devtools::load_all()` build with src/Makevars flags):

```text
scenario                          seconds   ratio vs plain LBA raw path
plain_lba_2                        0.007      1.0
corr_lba_2_rho50                   0.139     19.9
corr_lba_2_rho80                   0.196     28.0
corr_lba_2_rho95                   0.226     32.3
corr_lba_2_rho08_unrestricted      0.120     17.1
corr_lba_2_rho08_near_t0           0.215     30.7
leak_bawl_2_rho08                  0.214     30.6
plain_bawl_pm_race23               0.010      1.4
corr_lba_pm_race23_rho08           0.382     54.6
leak_bawl_pm_race23_rho08          0.398     56.9
forced_gh_3loaded_rho08            0.500     71.4
forced_generic_clock_rho08         4.107    586.7
```

Route/node counters confirm the baseline issues named above:

- `corr_lba_pm_race23_rho08`: 9,826 of 20,000 trials (every RACE=3 trial)
  are misclassified as three-dimensional and run the 64-node quadrature
  denominator (628,864 node evaluations) even though PM has exact rho=0.
- Numerator work: 10-node scan at |rho|<=.6 (200,000 node evals), 12-node
  above; refinement re-runs nearly all trials at high rho (19,762/20,000 at
  rho=.95), i.e. the scan-reuse shortcut rarely helps exactly where cost is
  highest.
- Generic clock fallback is ~14x the no-clock GH route at the same rho.

### Step 2 — shared geometry, state, pointers, layout (DONE, 2026-07-17)

- `src/bawl_geometry.h`: `BAwLTimeGeometry` (statuses valid / point_start /
  not_started / infinite / invalid) built on `bawl_leak_factors()` and the
  exact `k -> 0` limit; `BAwLPreparedRow` conditional view; prepared
  unnormalised endpoint evaluators (`f0`, `F0`, `q`, `1-F0`, `q-F0`) with the
  raw kernels' RAW-acceptance natural branch and authoritative log fallback.
  Parity-tested against `dleakyba`/`pleakyba` over the T1 grid (k x tau x
  A x v x sv x rho x z, 648 cases) plus survivor/cause coefficient
  identities and degenerate statuses.
- `BAwLCorrSharedState` built once per likelihood call: lM role mapping,
  RACE masks, truncation windows, winner rows, expansion map, and direct
  ParamTable column pointers (keep_names order).  The correlated entry point
  no longer receives a materialised matrix; the generic-clock fallback and
  the rho=0 ordinary route materialise lazily via a callback, and the fast
  node path reads `ParamTable` columns directly.
- `bawl_corr_classify_particle()`: canonical per-particle
  `BAwLCorrTrialLayout` (active/loaded rows, pair candidates, winner
  loading, route enum).  Routing behaviour unchanged in this step; layout
  feeds the loaded-dimension counters (T3 test asserts dimension 2 for
  correct/error + independent PM, 3+ when PM is loaded, ordinary at rho=0).
- Verified: benchmark lls bit-identical pre/post refactor; full
  `test-bawl-correlated.R` green from a fresh temp-library install.

### Step 3 — denominator and normal-CDF fixes (DONE, 2026-07-17)

- Positivity dimension now comes from the layout's loaded-row count; the
  closed-form denominator builds its 1-D/BVN arguments from loaded rows
  only.  Zero-loading rows' constant q factors are cancelled algebraically
  from both the numerator node integrands (fast and generic evaluators) and
  the denominator.
- Correct/error + independent PM: `den_quadrature_trials` is now 0 (was
  every RACE=3 trial); PM with nonzero rho still routes to the quadrature
  denominator (asserted in the T3 test).
- All correlated positivity sites use `pnorm_log_direct()` instead of
  `R::pnorm(..., log=TRUE)` (closed-form 1-D denominator, generic
  `log_positive_trial`, fast node q loop).
- Benchmark (N=2000): PM scenarios 0.039/0.040 s -> 0.018/0.020 s (~2.1x);
  all correlated scenarios gained ~15-25% from the pnorm switch;
  `forced_gh_3loaded` 0.050 -> 0.035 s.  Verified against a dense 800-node
  GH reference on a PM trial: |diff| = 1.5e-6 (the removed 64-node sweep
  carried ~1.5e-3 of quadrature error per trial).  Full correlated test
  file green from a fresh temp-library install.

### Step 4 — final fused, unnormalised no-clock fallback (DONE, 2026-07-17)

- The no-clock GH route now prepares one `BAwLPreparedRow` per active
  nonzero-loading row from direct `ParamTable` columns.  GH nodes update only
  the conditional mean; fixed-time/leak geometry, endpoint spans, and scale
  terms are not rebuilt inside the node loop, and the old full conditional
  parameter-table staging has been removed.
- Prepared endpoint modes now provide the unnormalised winner density `f0`
  and loser survivor `q-F0`; the product of loaded-row `q` terms is retained
  only for the positive-drift denominator.  Exact-zero-loading rows are
  evaluated once as independent singleton factors and remain out of the
  shared-factor denominator dimension.
- Prepared natural-CDF and survivor paths preserve the scalar kernel's
  central, underflow, near-one, and saturated-tail acceptance behaviour.
  Active local/global clock variants remain on the generic GH fallback.
- Route observability reports prepared rows and fused node evaluations.  At
  N=2000 and one timing repetition, prepared no-clock scenarios measured
  0.014--0.040 s per call; the forced generic-clock scenario measured 0.870 s.
  The three-loaded scenario reported 6,000 prepared rows and 176,000 fused
  node evaluations per call.
- Full `tests/testthat/test-bawl-correlated.R`: 5,033 passing, 2 expected
  skips, 0 failures.  The new T4 test checks prepared-vs-generic parity and
  confirms that the generic-clock route does not consume prepared geometry.

### Step 5 — shared probability and moment primitives (DONE, 2026-07-17)

- Added `src/bawl_corr_exact.h` with guarded univariate normal interval
  moments through order three and deterministic BVN boundary grids.
- Rectangle probability, first moments, and cross moment use signed Kahan
  corner combinations and the mean-derivative identities.  The mixed
  derivative uses the matching marginal/conditional density pair, including
  asymmetric rectangles.
- Tiny/cancelled rectangles return `unstable` rather than being clamped; true
  zero-mass rectangles retain a separate status.  The pair-positive
  normalizer has an explicit nearly singular Gaussian limit.
- Added `tests/testthat/test-bvn-rect-moments.R`; conditional-normal reference
  integration passes central, tail, narrow, and semi-infinite rectangle cases
  at rho `-.8`, `0`, and `.8`.

### Step 6 — exact pair scalar component (DONE, 2026-07-17)

- Implemented exact pair survivor and member-cause formulas from the shared
  `BAwLTimeGeometry`, including LBA, leak-only BAwL, unrestricted and jointly
  positive drifts, unequal loadings, and independent numeric-pair fallback.
- Point starts use the dedicated boundary-density/Jacobian branch; point
  starts for either racer and both racers are covered by regression tests.
- Exact rectangle values agree with the independent 64-node one-dimensional
  pair integration throughout the stable test grid.  At narrow or strongly
  negative tail rectangles the status route selects `numeric_pair`.

### Step 7 — exact finite component wiring (DONE, 2026-07-17)

- Two loaded no-clock trials now dispatch through one exact per-trial
  component evaluator.  Pair winners, independent/PM winners, arbitrary
  singleton survivors, RACE masks, and expansion use the same assembly for
  total and pointwise APIs.
- The exact-only path returns before GH work; mixed particles retain masked
  per-trial fallback routing.  Exact paths use direct ParamTable column
  pointers and do not materialize a parameter matrix.

### Step 8 — nonfinite, truncation, and censoring wiring (DONE, 2026-07-17)

- Added component-survivor assembly for truncation normalization, unknown
  winner censoring, missing-RT interval unions, and known-winner outer-time
  integration.  `+Inf` uses the explicit leaky-tail limit.
- Contamination, truncation correction, flooring, and expansion occur once at
  final per-trial assembly.  Focused tests cover finite, `+Inf`, `-Inf`,
  missing-RT, known/unknown, and pointwise cases.

### Step 9 — analytic GH centre (DONE, 2026-07-17)

- Added the shared order-three moment centre for loaded-winner GH fallback
  integrands, with sign-aware latent-factor mean and variance and a selective
  fine pass.  Survivor-only and unknown-winner integrands retain scan-based
  centres.
- The centre is applied only on the fused fallback; exact and singleton
  evaluators are unchanged.  Benchmark counters report successful centres on
  the three-loaded fallback route.

### Step 10 — numerical routing and fallback cleanup (DONE, 2026-07-17)

- Exact rectangle instability routes first to `numeric_pair`; numeric failure
  falls through to the fused or generic GH route.  Near `|rho|=1` uses the
  continuous positive-normalizer limit and numeric pair integration.
- Exact-only early return removes fallback-node work for the canonical pair
  path, while route counters now count only actual fallback trials.
- N=2,000 one-particle benchmark (one repetition, this machine): plain LBA
  0.001 s; exact correlated pair 0.008--0.009 s; correlated pair plus PM
  0.009--0.010 s; three-loaded GH 0.041 s; generic-clock GH 0.874 s.  The
  exact path used 1,998--2,000 exact trials in the pair scenarios, with only
  2 numeric-pair trials at rho `.5` and 250 in the near-`t0` stress case.

### Step 11 — regression, documentation, and cleanup (DONE, 2026-07-17)

- Focused correlated regression: 5,032 passing, 2 expected skips.  The new
  BVN/pair/data-path file also passes.
- The workspace-wide sequential test run reaches the correlated and new exact
  tests successfully.  Remaining failures are outside this correlated change:
  an existing map helper lookup, stochastic trend snapshots, and missing Wald
  helper symbols in the current test environment.
- `R CMD INSTALL --no-multiarch --no-test-load` succeeds in a fresh temporary
  library, `Rcpp::compileAttributes()` is up to date, and `git diff --check`
  is clean.  A full `R CMD check` is environment-blocked by the unavailable
  suggested package `DiagrammeR`; the source objects generated during the
  verification build were removed afterward.

### Step 12 — production hardening from the Annika preburn profile (DONE, 2026-07-17)

Profiling a real hierarchical preburn (Annika control data, 114 subjects,
LT=.1/UC=2.75 truncation-censoring, 11,718 unique trials) showed the fitted
run paying ~14x over independent BAwL rather than the benchmarked ~8x, with
15.7% of pair trials bouncing to the 64-node numeric route at sampled
parameters and the truncation normalizer recomputed per unique trial.

- `bawl_corr_rect_from_grid`: the p<=0 zero-mass test now accepts rectangles
  whose boundary-derivative scale is below 1e-14 (was 1e-300).  The mass of
  such a rectangle is bounded far below the 1e-10 likelihood floor, so these
  are genuine tail rectangles at small tau, not cancellations; previously
  each one forced its whole trial numeric.
- Per-particle cell memoization in `BAwLCorrSharedState`: the positive-drift
  orthant q_AB (keyed on mu/sd pair + rho) and the truncation normalizer
  log Z (keyed on pair row parameters + LT/UT; pure pair trials only) are
  now computed once per design cell instead of once per unique trial.
  Caches are cleared per particle because the ParamTable base is refilled in
  place.
- Verified from a fresh temp-library install: benchmark lls identical across
  all scenarios; near-t0 stress 13x -> 8x with 0 numeric-pair trials (was
  250/2,000); focused correlated tests 5,035 + 251 passing, 0 failures.
  Real-data parity vs the pre-fix build: max |dll| 6.4e-5 (old numeric-route
  quadrature error), numeric share 9-15% -> 0%.  End-to-end 3-iteration
  preburn benchmark: 527 s -> 201 s (independent BAwL: 37 s).

## Implementation sequence

This order deliberately creates final shared infrastructure before either the
fused or exact kernel. No step should add a buffer, callback, classifier, or
formula that a later step replaces.

1. **Baseline and observability**
   - Add benchmark and gated counters.
   - Record current nodes, 64-node PM incidence, accuracy, and timings.

2. **Shared geometry, state, pointers, and layout**
   - Add `BAwLTimeGeometry` built on the same authoritative helpers
     (`bawl_leak_factors()`, the `k -> 0` limit) and verify parity against
     the existing scalar kernels (T1). Do not rewrite the ordinary scalar
     kernels themselves; they remain the bit-for-bit `rho=0` reference until
     step 11.
   - Add `BAwLConditionalGeometry` and shared endpoint helpers.
   - Add final `BAwLCorrSharedState`, direct column pointers, and canonical
     `BAwLCorrTrialLayout` (T3).
   - Route all trials to existing behaviour initially, but remove duplicated
     role/RACE/loaded-row discovery.

3. **Easy denominator and normal-CDF fixes on shared layout**
   - Count only layout-loaded rows.
   - Cancel zero-loading positivity constants consistently.
   - Keep correct/error + PM out of the 64-node path.
   - Replace direct correlated `R::pnorm` calls (T2, T3).

4. **Final fused, unnormalised no-clock fallback**
   - Implement `f0` and `q-F0` modes inside the shared prepared endpoint
     evaluator; do not add temporary raw callbacks.
   - Implement fused trial/node reduction and benchmark both loop orders.
   - Route no-clock GH through it; retain generic clocks (T4).

5. **Shared probability/moment primitives**
   - Add incomplete univariate moments used by the future centre and the
     `numeric_pair` conditional expectations.
   - Add deterministic BVN boundary grids, rectangle derivatives, and
     moments (T5).

6. **Exact pair scalar component**
   - Implement pair survivor/cause with shared time geometry and cached
     boundary grids.
   - Implement the point-start (`A=0`) branches.
   - Implement the direct one-dimensional `numeric_pair` fallback.
   - Add positive denominator, rho limits, and status routing.
   - Validate against `numeric_pair` and the probability identities (T6).

7. **Exact finite component wiring**
   - Activate exact-pair layout route.
   - Combine exact pair with batched singleton masks for pair and PM winners.
   - Feed one unique-trial evaluator to total and pointwise APIs (T7).

8. **Nonfinite/truncation/censoring wiring**
   - Add total component survivor and leaky `+Inf` limits.
   - Add exact component outer-time integrand for known-winner intervals.
   - Finalise contaminant/truncation/floor/expansion once (T8).

9. **Analytic centre on the final fused fallback**
   - Use shared incomplete moments and time geometry.
   - Apply only to loaded-winner integrands; survivor-only integrands keep
     scan/recentre.
   - Add one-pass evaluation and selective scan/refine masks.
   - Do not modify exact or singleton evaluators (T4, T9).

10. **Numerical routing and optional compaction**
    - Tune exact/numeric_pair/GH stability thresholds using T5/T6/T9 and
      observed route counters.
    - Compact fallback trials only if profiling demonstrates a benefit.

11. **Regression, documentation, and cleanup**
    - Run T10 and final benchmarks.
    - Only now consider refactoring the ordinary scalar kernels onto
      `BAwLTimeGeometry`, and only if profiling or maintainability justifies
      it; rerun T1/T10 if done.
    - Remove temporary probes, not shared primitives.
    - Update comments that describe GH as the only correlated method.
    - Record accuracy, route counts, and performance here.

---

## Acceptance criteria

1. Correct/error uses the exact path with both `RACE=2` and `RACE=3`.
2. A PM winner uses `f_PM*S_pair`.
3. Any number of zero-loading racers factor through singleton kernels.
4. Independent PM never increases positivity dimension or invokes 64 nodes.
5. Numerator and positivity-denominator routes are independent: a quadrature
   denominator never alters the numerator node schedule or disables
   numerator reuse.
6. Inactive RACE rows affect neither likelihood nor normalisation.
7. Rho zero remains bit-for-bit ordinary LBA/BAwL.
8. Positive drift retains jointly truncated correlated-Gaussian semantics.
9. Exact LBA/leak values meet the accuracy grid or route smoothly through
   numeric_pair to GH, including point starts and unequal loading
   magnitudes.
10. Truncation and censoring go through `log S_all(t)` with every survivor
    factor evaluated at each endpoint; no loaded-group-only window formula
    exists.
11. The no-clock GH path prepares time/leak geometry once per row/particle.
12. Active clocks and >2 loaded racers retain accurate fallback behaviour.
13. Total and pointwise likelihoods share one per-trial evaluator.
14. Exact/fused all-finite paths use direct column pointers without Rcpp
    matrix clone/materialisation.
15. No duplicated LBA/leak geometry, role mapping, loaded-row classification,
    exact formulas, or total/pointwise assembly remains within the
    correlated routes; the ordinary scalar kernels may retain their own
    geometry until the optional step-11 refactor.
16. Benchmarks show a practical improvement large enough to justify fitting.

## Verification notes

- Use a fresh `R CMD INSTALL` into a temporary library; the default installed
  EMC2 package may be stale.
- Use compiler flags matching `src/Makevars`.
- Preserve GH node-count environment overrides for convergence tests.
- Do not floor intermediate terms.
- Keep public model, simulator, and transform behaviour unchanged.
