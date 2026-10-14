# WS4: Test problems whose geometry we know

Read `README.md` in this directory first. It covers the shared rules (worktree, private library, machine lock, comment style, exactness) and the current state of the code.

## Goal

Build a small, cheap set of test problems for the hierarchical sampler. For each problem:
- we *know*, by construction or by measurement, which kind of slow mixing it contains;
- we know how slowly a plain sampler should mix on it;
- we have a trustworthy reference posterior.

Then a mechanism can be judged on the cases it was designed for, and on cases where it should cost almost nothing.

The set has two parts:
1. **Mock models:** Gaussian-type hierarchical models where every geometric feature is a dial.
2. **A principled way to place real EAM datasets (LBA, RDM, LNR, DDM) on the same map:** measured geometry descriptors, and generators that hit chosen points on that map.

## Why

Every bottleneck found so far came from one fixture, and fixing it exposed the next:
- correlated group means (T30);
- curved subject ridges (lbaT);
- log-sd funnels (T30, lbaE);
- floors at bounds (N-back t0);
- unidentified parameters. lbaT leaves A unidentified, so it tested identification, not mixing.

Fixtures mix several geometries, so no result could be attributed to one mechanism, and comparisons drowned in seed noise (continuation ESS varied 7x across seeds on N-back). The user will not fund an ablation matrix until the test set makes its results interpretable. This workstream is that prerequisite. WS1 (allocator) and WS3 (particle step) will consume its output.

## The geometries that matter, and the mechanism meant for each

| Geometry | What creates it | Mechanism designed for it |
|---|---|---|
| Strongly correlated group parameters | Correlated subject effects, few subjects | Group ensemble update |
| Group means coupled to subjects (centred-parameterisation slowness) | Subjects weakly informed relative to the group sd: shrinkage near 1 | Shear rounds (interweaving) along mean directions |
| Log-sd funnel | Small group sd relative to subject posterior sd; few subjects | Shear log-sd rounds; collapsed group move |
| Curved subject posteriors (ridges, bananas) | Non-linear trade-offs, e.g. LBA A/B/t0, threshold against rate | Subject HMC; curved (state-indexed) shear frames |
| Mass against declared bounds | Parameters near a bound, e.g. t0 near its floor or a minimum rt | Bound walls in the shear/collapsed transport |
| Non-identified or multi-modal subjects | Too little data per parameter | None. The test set must *screen these out* or label them, because no sampler change fixes identification |
| Easy | Informative subjects, moderate group sd | None should engage; measures pure overhead |

## Work

### 1. Mock model framework

A hierarchical model where subject s has a likelihood in alpha_s we specify exactly. It must run through the *full* sampler: worker pool, particle step, shear, collapsed move, HMC, allocator. Otherwise it tests nothing.

**Build the mocks as proper EMC2 models, not as a side hook.**
- Each mock is a new model structure that goes through the package exactly as an EAM does: `design()`, `make_emc()`, `fit()`, the C++ likelihood path, the worker pool, every move. It just has a trivially simple likelihood.
- The template is REXG: a thin model that wraps an ex-Gaussian density.
  - Model definition: `REXG()` in `R/model_LNR.R` (around line 431).
  - C++ dispatch: `src/race_dispatch.cpp`.
  - Column registration: `src/col_registry.h`.
  - Read `project_emc2_race_kernel_contract.md` in memory first: the `p_types` prefix is the kernel's column order, and `Ttransform` is R-only.
- Write, say, a Gaussian model whose observation column has a normal density with mean and sd from mapped parameters. The design formulas then give the linear structure: which parameters vary over which within-subject cells, which sets X and the information matrix.
- Work out how an observation that is not a response time travels through the data columns and checks. Negative values, `R` with one level, and rt validation are the things to check. Choose the cleanest encoding the existing contract allows; do not special-case the sampler.
- Add a multivariate or non-linear variant (product ridge, banana) as further small models of the same kind, if the geometry needs a likelihood that a mapped univariate normal cannot express.
- Do not use the custom-likelihood hook (`design_model_custom_ll`). It bypasses the C++ path that the sampler's moves use.
- Give each new model the standard model tests (density against a numerical reference, `make_data` simulation, a recovery smoke test).

**The dials:**
- **Linear-Gaussian:** `y_sj ~ N(x_j' alpha_s, sigma^2)`. The subject's information matrix `I_s = X'X / sigma^2` is set exactly by the design X and the trial count. This gives correlations, informativeness and shrinkage. With a fixed group covariance the posterior is conjugate, so mu given V is exact. The full posterior with the IW / diagonal-gamma priors can be computed to high accuracy by numerical integration or a long independent reference.
- **Non-linear links** for curvature: a product ridge `y ~ N(a * b, sigma)`, a banana `y ~ N((a, b + k a^2), Sigma)`, a soft bound via `exp`.
- **Group structure:** number of subjects, group sd (funnel strength), group correlations.
- **Bounds:** declared parameter bounds with mass near them.
- **Multimodality:** for the screening check, as a negative control.

**Reference posteriors:**
- exact where conjugate;
- otherwise independent long runs, plus agreement between deliberately different samplers.

Include a small simulation-based calibration (SBC) harness on two or three mock models, checking that the *full adaptive* sampler is calibrated. The package has `run_sbc()` in `R/SBC.R`; check whether it fits, and keep each SBC run cheap.

### 2. Geometry descriptors (the principled part)

Define a few computable descriptors that predict how a plain centred sampler mixes on a given hierarchical posterior. Compute them from the model and data at a reference state (the truth or a posterior mode), via subject Laplace fits and the group covariance. Candidates; derive and check each:

- **Shrinkage**, per subject: `B_s = (I_s + V^-1)^-1 V^-1`, the prior's share of the subject's posterior precision.
  - In the linear-Gaussian model, the lag-one autocorrelation of the group mean under centred Gibbs along direction u is u' (mean_s B_s) u (Papaspiliopoulos, Roberts and Sköld, 2007, and related work).
  - It is near 1 when subjects are weakly informed: centred is slow there, and the shear/interweaving was built for it.
- **Funnel strength** for each log sd: how much of the log sd's marginal variance survives conditioning on the subjects. In Gaussian approximation that depends on n_subjects and B; derive it.
- **Group correlation:** the condition number and leading correlations of the group posterior covariance. This is what the ensemble update targets.
- **Non-Gaussianity** of each subject posterior:
  - the gap between the Laplace quadratic and the true log density at ±1 and ±2 posterior sds along principal axes (as in the lbaE orbit scan, `project_emc2_lbaE_orbit_scan.md`);
  - and/or an importance-sampling estimate of KL(Laplace || posterior).
- **Bound proximity:** posterior mass within a small margin of each declared bound.
- **Identification:** multi-start optimum spread and Hessian condition per subject. `benchmarks/subject-hmc/screen_identification.R` exists.

For each descriptor, validate its prediction on the mock models, where the truth is known. For example, predicted versus measured lag-one autocorrelation of the group means under a group-moves-off sampler, which is available through `emc2.*` switches; see `harness/fit_arm2.R` and the arms in the README.

A descriptor that does not predict is dropped or fixed, not kept for decoration.

### 3. Placing real EAM data on the map

Write a generator for LBA, RDM, LNR and DDM hierarchical data whose dials are:
- design (which parameters vary by condition);
- trials per cell;
- number of subjects;
- group means, sds and correlations.

Compute the descriptors for each generated dataset at the truth: subject information by finite differences of the C++ likelihood, or expected information by simulation. Then choose dial settings that land each EAM fixture in one target regime from the table above, with the others mild.

Also compute the descriptors for the existing real fixtures (`fixture_nback31.rds`, `lba7r_*`, `lbaE_fit.rds`, Forstmann, T30). That says which regimes they actually occupy; we have been guessing.

### 4. The catalogue and the runner

**A catalogue of about 8-12 cheap fixtures** (a full-sampler fit in minutes, not hours). Each entry records:
- the generator settings and seed;
- its descriptor values;
- its intended regime and the mechanism designed for it;
- the predicted plain-sampler mixing (from the descriptors);
- its reference posterior;
- whether it is mock or EAM.

Include at least one easy control in each family.

**A runner** that, for a given build and option set, runs a fixed short protocol on each fixture and reports:
- correctness: agreement with the reference posterior in MCSE units;
- worst-coordinate ESS per CPU-second and per wall-second;
- which mechanisms engaged, with counts (rounds, acceptances, HMC transitions, collapsed moves);
- wall and CPU.

The protocol is fixed-length continuations from a shared stationary checkpoint; see `project_emc2_benchmark_protocol.md` and `benchmarks/subject-hmc/continue_arm.R`. Measure each fixture's null spread (the same arm over several seeds) so the number of seeds needed to separate a given effect is known before anyone runs an ablation.

## Interim checkpoints (other workstreams wait on these)

Publish progress in `/data/work/EMC2_workstreams/geometry/STATUS.md`: what is ready, where it is, how to call it.
1. The mock framework running through the full sampler, plus the linear-Gaussian family with exact references. **Publish this as soon as it works.** WS1 and WS3 need it for calibration.
2. Descriptors validated on the mocks.
3. EAM generators and the catalogue.
4. The runner, with each fixture's null spread.

## Deliverables

- Code under `benchmarks/geometry/`: the framework, descriptors, generators, catalogue, runner and SBC harness, with a README.
- The mock models as proper EMC2 model structures (R definition, C++ kernel, column registration), each with its own tests.
- `ws4_report.md`, containing:
  - the descriptor derivations and validation plots or tables;
  - the catalogue table;
  - where the existing real fixtures sit on the map;
  - null spreads;
  - what each regime predicts for each mechanism;
  - what the runner cannot yet tell apart.

## Boundaries

- Do not change sampler behaviour. If a mock exposes a sampler bug, document it with a reproducer in your report and in `STATUS.md`, and leave the fix to the owning workstream (WS1 allocator, WS3 particle step) or to the user.
- Keep fixtures cheap. Every fit and continuation that uses more than 4 cores, or whose time you report, goes through `bench_lock.sh`.
- No stop-rule fits to explore.
