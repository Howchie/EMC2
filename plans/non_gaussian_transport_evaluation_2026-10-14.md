The appropriate replacement is a full-dimensional, non-Gaussian conditional
move with an explicit reversible construction. However, the constructions
screened here do not yet justify integration. A fixed, nonlinear likelihood
correction was inexpensive but inaccurate. Exact conditional bridges improved
acceptance beyond the reported Gaussian-refit ceiling, but their extra work
did not produce a convincing efficiency advantage. Following the stopping
condition in `shear_lsd_transport_2026-10-14.md`, this evaluation stops before
production changes, continuations, or stop-rule benchmarks.

This is a negative screening result and a concrete architecture specification,
not a claim that an efficient replacement has been established. The sampler in
`R/`, `src/`, and `tests/` is unchanged.

The evaluation uses HEAD `fb1f3f45`, the restored C++ likelihood path, and
`/data/work/EMC2_workstreams/fixtures/lbaE_replay_rebuild_s1.rds`. Timed work was
serialized under the shared benchmark lock. Existing files and local changes
were preserved. The graph provided an initial architecture index; the recent
shear implementation and experiments were inspected directly because they
postdate much of that index.

The evidence changes the design requirements in several ways. The old curved
mechanism's occasional successful seeds are a reason to investigate geometry,
but do not establish which mechanism is needed or that baseline failure within
a time budget implies lack of ergodicity. More importantly, the additional
`Claude Result 2026-10-14.md` and line-collapse experiment already establish a
failure mode absent from the main plan's candidate list: removing a subject's
cliff along one line can still leave its other coordinates pinning the group.
That rules out repeating that construction. It does not prove that every
one-dimensional nonlinear correction to a full-dimensional map must fail.

There are three distinct obstacles: nonlinear conditional shape, accumulation
of subject errors, and genuine posterior dependence between the proposed group
coordinate and the remaining group coordinates. A proposal can address the
first while failing on either of the others. Even an exact subject transport
cannot make the true marginalized group density flat. The reported refit
ceiling and Gaussian full-collapse widths are useful empirical diagnostics,
not mathematical upper bounds for non-Gaussian transports.

All numerical material below is in `shear_lsd_transport_2026-10-14/`.

The first construction represents a subject's likelihood by

\[
 \widehat L_s(a)=\exp\{h_s^Ta-\tfrac12a^TH_sa+f_s(v_s^Ta)\},\qquad H_s\succeq0.
\]

The nonlinear function is a continuous piecewise-linear log-likelihood
correction with linear tails. Multiplying by any finite positive-definite
Gaussian group prior gives a proper conditional surrogate. Unlike the removed
atlas, this approximates the likelihood as a function of subject parameters;
it is not a conditional mode fitted at one group state. Group dependence enters
through the exact prior, so no likelihood optimization is required at candidate
states.

Let the Gaussian part of that conditional have mean m and covariance S. Set
x=v'a, m_x=v'm, s_x²=v'Sv, and b=Sv/s_x². The marginal of x is proportional to
`N(x | m_x, s_x²) exp(f(x))`. Each segment integrates analytically by completing
the square. Conditional on x, all remaining coordinates are Gaussian with mean
`m + b(x-m_x)` and covariance `S - s_x² bb'` on the complement of v.

The map first preserves x's nonlinear marginal CDF, then preserves the
standardized Gaussian complement residual. Every subject coordinate can move.
This ordering is the substantive distinction from line collapse and from a
map that moves a Gaussian complement first and corrects a nonlinear line last.
The log Jacobian is the old normalized scalar log density minus the new one,
plus the difference of the complement covariance factors' log determinants.
Equivalently it is `log q_old(a) - log q_new(a_new)`. The old shear selection
rule remains applicable when these are fixed, state-indexed charts that
compose and invert.

Two fitting methods were screened. A joint quadratic and hinge regression on
saved window likelihoods selected t0 for 73 of 80 subjects. On draw 1001 it
increased summed subject orbit spread from 1,797 to 159,977 and nearly
eliminated acceptance. This was stopped after its smoke replay. The second
method retained the Laplace quadratic and fitted nonlinear corrections to
*likelihood residual differences along transports*. It used draws 600, 700,
800, and 900 as probes, excluded the four evaluation draws, and selected among
parameter axes and threshold contrasts. Its fit cost was 0.469 elapsed seconds
for all 80 subjects, including fresh C++ likelihood probes.

On chain 1, evaluation draws 1001, 1000, 980, and 950, the work-fitted method
gave the following arithmetic means. The spread is the sum over subjects of
their orbit log-density variances on the grid -2 to +2 in increments of 0.5
window standard deviations. It excludes the separate group prior term.

| Map | Subject orbit spread | Acceptance -1 | -0.5 | +0.5 | +1 |
|---|---:|---:|---:|---:|---:|
| Existing fixed frames and walls | 21,824 | .015 | .173 | .257 | .0010 |
| Same Gaussian likelihoods, nonlinear-first ordering with correction zero | 95,461 | .047 | .265 | .134 | .00033 |
| Work-fitted nonlinear correction | 117,403 | .014 | .138 | .324 | .00037 |

The zero-correction control matters: changing the ordering of Gaussian
whitening already changes the orbit on a non-Gaussian target. Improvements
against one ordering alone would not establish an improvement against the
current sampler. Evaluation draw 950 dominates the spread, and the correction
does not repair it. The nonlinear correction itself is therefore rejected,
not merely its fitting cost. The prototype took roughly 0.38–0.46 serial CPU
seconds for nine states, about 0.043–0.051 seconds per state including exact
likelihood validation. These are R prototype timings, not production claims.

The second construction uses exact likelihoods throughout. It is a transported
version of [Neal's conditional dragging move](https://arxiv.org/abs/math/0502099).
It gives a useful reference for how much nonlinear conditional relaxation can
buy without fitting curves, and supplies a precise fallback design if a future
proposal makes such relaxation cheap enough.

For group endpoints g0 and g1, propose a symmetric increment in a frozen group
coordinate scale. Choose their group midpoint and construct its subject frames.
Use the same midpoint under endpoint reversal. In that chart let z_s denote a
subject residual and B_{s,g}(z_s) its transport to endpoint g, using the existing
wall-aware maps. Define the *unnormalized* endpoint density with respect to z:

\[
 F_{s,g}(z)=L_s(B_{s,g}(z))\,
           N(B_{s,g}(z)\mid\mu_g,V_g)\,|D_zB_{s,g}(z)|.
\]

The midpoint whitening determinant is included. The group density P(g)
includes the group priors and the covariance-coordinate Jacobian already used
by shear. Huang–Wand auxiliary variables remain integrated out during the move
and are redrawn afterward as in the existing implementation.

For subject s choose a bridge length K_s. Starting from z_{s,0}, at intermediate
index i=1,...,K_s-1 apply a reversible transition for

\[
 \rho_{s,i}(z)\ \propto\
 F_{s,g0}(z)^{1-i/K_s}F_{s,g1}(z)^{i/K_s}.
\]

The prototype uses Metropolis-corrected pCN proposals, with a fixed mixture of
three scales. Thus its intermediate kernels use the actual non-Gaussian
likelihood, not an approximation to it. Return the final subject state and

\[
 W_s={1\over K_s}\sum_{i=0}^{K_s-1}
       [\log F_{s,g1}(z_{s,i})-\log F_{s,g0}(z_{s,i})].
\]

Accept the joint group and subject proposal with

\[
 \min\{1,\exp[\log P(g1)-\log P(g0)+\sum_sW_s]
                 q(g0\mid g1)/q(g1\mid g0)\}.
\]

This follows by taking the ratio of forward and reverse path probabilities.
Reversibility of each intermediate transition cancels the intermediate density
terms, leaving the stated average endpoint log ratios. No conditional
normalizing constants are estimated. A finite bridge is exact despite poor
mixing; poor mixing reduces acceptance. K_s=1 is the ordinary transported
Metropolis proposal. Different subjects can have different fixed K_s because
their path ratios factorize. An identity transition for a subject is also
valid, although it does not remove that subject's rejection contribution.

This is a pairwise pathwise Metropolis kernel. Its stochastic endpoints cannot
be inserted into the existing eight-candidate orbit weighting formula without
another extended-target derivation. Treating them as deterministic orbit
candidates would invalidate the move.

Uniform bridges were evaluated at the same eight chain-1 draws used to reproduce
the original fixed-frame acceptance figures: 1001, 1000, 980, 950, 900, 800, 700,
and 600. Each nontrivial bridge has 12 independent randomized paths per draw
and jump. These are replay states from one chain and seed, not independent
convergence experiments. Monte Carlo error and dependence between draws remain.

| Method | Acceptance -1 | -0.5 | +0.5 | +1 |
|---|---:|---:|---:|---:|
| Fixed transport, K=1 | .031 | .200 | .393 | .128 |
| Full Gaussian refit, historical plan result | .10 | .37 | .45 | .16 |
| Exact bridge, K=4 | .084 | .282 | .427 | .179 |
| Exact bridge, K=8 | .110 | .322 | .502 | .209 |
| Exact bridge, K=16 | .170 | .392 | .538 | .241 |

The K=1 results reproduce the plan's fixed-frame row. The refit row is copied
from the plan rather than rerun here. The K=16 bridge improves all four mean
acceptance figures over that historical ceiling. This confirms that the
Gaussian-refit result is not a fundamental limit. It does not establish a
competitive sampler.

For a symmetric proposal choosing equally among these four jumps, the expected
squared accepted displacement in window variance units is .0768, .1102, .1312,
and .1607 at K=1,4,8,16 respectively. This includes Metropolis rejection. A
separate replay of the actual eight-candidate fixed-frame round, using 12
candidate ensembles per draw and a scale grid from .05 to 1.6 window standard
deviations, attained a mean selection-weighted displacement of .1018 at scale
.8. This is a scale-grid replay comparison, not a tuned production continuation.
The bridge therefore buys only about 1.08x, 1.29x, and 1.58x displacement per
proposal at K=4,8,16 against that existing round.

The likelihood cost explains the problem. Fresh timing on five representative
subjects, 100 calls per subject and batch size, gave mean costs of .230 ms for
one row, .298 ms for two, .638 ms for eight, and 1.098 ms for sixteen. A cached
current likelihood lets a bridge use one fresh endpoint row plus K-1 calls with
two endpoint rows each. That is K serial calls per subject, versus one
eight-row call for the existing round. Applying these measured means to all
80 subjects estimates likelihood-only CPU costs of .090, .185, and .376 core
seconds for K=4,8,16, compared with .051 core seconds for the existing round.
These extrapolations exclude transport, dispatch, and queue overhead and do not
pretend that five subjects measure the whole workload exactly. They are already
unfavorable: displacement grows much more slowly than likelihood work.

The actual fixed-frame replay used about .005 seconds of master candidate work
and .05–.06 seconds of serial subject evaluation per round. At nine reserved
cores, charging the master work as required by the plan makes its total price
roughly .10 core seconds. Even granting the bridge *zero* additional setup or
dispatch cost does not demonstrate the required clear gain: K=4 is at best a
modest improvement under this generous comparison, while longer bridges lose
ground. A C++ implementation could remove R overhead; it would not remove the
extra likelihood calls or their dependence within a subject path.

A budgeted bridge was also screened. From separate training draws 600,700,800,
and 900, it measured each subject's squared non-Gaussian likelihood-error change
on transported paths and greedily allocated bridge steps to reduce
`sum_s v_s / K_s`. Lengths were frozen before evaluation. This is the discrete
analogue of K_s proportional to sqrt(v_s/c_s), under an independent-work
variance approximation and fixed cost budget. With an average K of 2,4,8, it
gave extra steps to 20,45,68 subjects respectively; the maximum length was
27,64,64. On evaluation draws 1001,1000,980,950, its squared displacements were
.0378,.0453,.0697, compared with .0308 for K=1. At average K=8, its four jump
acceptances were .102,.287,.311,.028; uniform K=8 on those same draws achieved
.092,.333,.407,.086. Concentrating work from this pilot did not generalize well.
These experiments do not prove that every cost allocation or bridge proposal
will fail, but reject this concrete inexpensive allocation rule.

The bridge also does not eliminate the large-n problem by construction. It
replaces deterministic transport error by path work and dissipation. In a
well-mixed approximation, independent subject work variances still add, often
like sum_s v_s/K_s. If an extensive fraction of subjects needs relaxation, the
required total work grows with the population. Selecting only a few currently
bad subjects and repairing them after rejection is not an exact shortcut: that
selection must appear in a reversible kernel or an auxiliary-variable density.
It cannot be treated as free screening.

Correctness checks passed at the scope tested. The ridge map's inverse and
identity errors were below 1e-15 on the selected subjects; its analytic log
Jacobian matched finite-difference determinants within 4.8e-10. A separate
finite-state enumeration sums all bridge paths for two subjects and two group
states. Both uniform lengths and unequal lengths, including (1,4), satisfy
detailed balance within 2.8e-17, and transition rows sum to one within 2.3e-16.
These validate the formula and representative numerical operations, not every
tail case of an eventual C++ implementation. The prototype CDF helper has
finite numerical search limits; production would require robust log-tail
inversion and explicit failure handling rather than silent clipping.

If later work produces a bridge or nonlinear chart that passes the efficiency
screen, integration should follow this concrete boundary:

1. Keep one sampler and a frozen per-coordinate transition plan. Preserve the
   existing centered ensemble and subject transitions. Replace a coordinate's
   shear round only when a directly measured competing kernel wins its funded
   pilot. Every slow coordinate receives pilot evidence independently of the
   rounds it was previously allocated; low allocation must not prevent testing.
2. Put each subject's entire bridge path in one C++ claim-queue task. Share the
   endpoint and midpoint frames for the round. Return the final parameters,
   final authoritative likelihood, path work, and actual cost. Intermediate
   bridge steps require no chain-master barriers or cross-subject communication.
3. Make one joint acceptance decision after all subjects return. On rejection,
   discard every path endpoint. On acceptance, update alpha, its likelihood
   cache, group state, and covariance inverse together. Use the existing
   integrated-prior auxiliary refresh. Keep subject RNG streams reproducible
   and independent of which worker claims a task; the offline prototype fixes
   its worker count and is not a scheduling-independent RNG implementation.
4. Freeze bridge lengths and proposal scales during each move. A plan may be
   fixed or selected symmetrically from the two group endpoints; it must not
   be shortened according to the realized forward path. Fit likelihood charts
   outside retained sampling, or provide a separate valid adaptive-MCMC
   argument. A fixed kernel preserving the posterior does not by itself prove
   exactness of a history-dependent adaptation scheme. Existing production
   rebuild adaptation is a separate architectural obligation, not resolved by
   these pathwise balance checks.
5. Price master work at the number of reserved cores and subject work at its
   measured cost. Select by accepted displacement per total cost initially,
   then by lag-k improvement against checkpoint variance in the prescribed
   400-iteration continuations. Acceptance, small inverse error, or a good
   likelihood fit alone are insufficient. Only after that run locked 3-chain,
   9-core benchmarks on at least three seeds and on the targeted other fixtures.

For deterministic charts, the same requirements imply fixed group-state
functions, genuine bijections with explicit tails and Jacobians, and the exact
likelihood in the outer correction. Historical subject draws may train a frozen
chart. In a conditional move, a map may also depend on coordinates held
invariant along that move. Thus the plan's statement that maps must not depend
on subjects should be read as a prohibition on unaccounted state-dependent
fitting, not on all conditional transport. The mathematical basis for approximate
triangular maps with exact outer correction is
[Parno and Marzouk](https://arxiv.org/abs/1412.5492).

I would not spend the next budget on another state-indexed Gaussian atlas or
on a production bridge implementation. A narrowly scoped next experiment, if
undertaken, should test whether a fixed low-dimensional *joint* nonlinear
likelihood correction captures the changing threshold/rate/scale trade-offs,
while transporting the full Gaussian complement. Such a model must be trained
and validated on transport-error changes, include interaction terms, and be
compared with its zero-correction ordering control. This remains an unvalidated
hypothesis; neither the present rank-one failure nor the older line-collapse
failure establishes that a rank-two version will succeed. No implementation or
benchmark success is claimed for it.

Reproduction, from the repository root, uses the commands below. Set
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1` for timing. The RDS outputs include
results and the fitted models where applicable. The CSV files provide tables
without restoring the fixture. The bridge logs contain randomized path results;
the balance check independently enumerates a finite target.

```bash
python3 plans/shear_lsd_transport_2026-10-14/check_bridge_balance.py

Rscript plans/shear_lsd_transport_2026-10-14/check_ridge_map.R \
  plans/shear_lsd_transport_2026-10-14/ridge_work_replay.rds \
  /data/work/EMC2_workstreams/fixtures/lbaE_replay_rebuild_s1.rds

/data/work/EMC2_workstreams/bench_lock.sh env \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 RIDGE_FIT=work \
  Rscript plans/shear_lsd_transport_2026-10-14/ridge_first.R \
  /data/work/EMC2_workstreams/fixtures/lbaE_replay_rebuild_s1.rds \
  /tmp/ridge-work.rds 1 1001,1000,980,950

/data/work/EMC2_workstreams/bench_lock.sh env \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  Rscript plans/shear_lsd_transport_2026-10-14/dragging.R \
  /data/work/EMC2_workstreams/fixtures/lbaE_replay_rebuild_s1.rds \
  /tmp/bridges.rds 1001,1000,980,950,900,800,700,600 1,4,8,16 12 3

/data/work/EMC2_workstreams/bench_lock.sh env \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 BRIDGE_ALLOCATION=work \
  Rscript plans/shear_lsd_transport_2026-10-14/dragging.R \
  /data/work/EMC2_workstreams/fixtures/lbaE_replay_rebuild_s1.rds \
  /tmp/allocated-bridges.rds 1001,1000,980,950 1,2,4,8 12 3

/data/work/EMC2_workstreams/bench_lock.sh env \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  Rscript plans/shear_lsd_transport_2026-10-14/shear_esjd.R \
  /data/work/EMC2_workstreams/fixtures/lbaE_replay_rebuild_s1.rds \
  /tmp/shear-displacement.rds
```
