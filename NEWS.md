# EMC2 3.4.0
## Development

-   **New model `OPR()`:** the opponent Poisson race. Each accumulator
    counts positive minus negative Poisson events to a continuous
    criterion `K`, with the exact Bessel first-passage density and the
    noncentral chi-squared CDF evaluated analytically in the compiled
    likelihood. `parameterization = "ratio"` (default) samples the net
    rate `v` and opposing-activity ratio `rho`; `"background"` samples `v`
    and a background rate `lambda0` in both streams; `"rates"` samples
    `lambda_p` and `lambda_m` directly and allows defective accumulators.

-   **Bounds are part of the transform (model change):** a declared bound
    side that the transform can represent now sets the transform's limit,
    so the sampler maps the whole real line inside the bound instead of
    hitting a likelihood cliff at it. For example, t0 with the default
    bound of 0.05 is now `t0 = 0.05 + exp(x)` (it was `exp(x)`, floored
    below 0.05), and `pContaminant` is `0.001 + 0.998 * pnorm(x)`. Lower
    sides of `exp` and both sides of `pnorm` are linked; `identity` and the
    upper side of `exp` keep the strict bound check. Sides declared in a
    model's or user's `transform` are kept as given, and a
    `design(bound = ...)` override moves the link with it. Parameters with
    a constant keep the plain link, so constants and exception values
    (`t0 = log(0)`, `A = log(0)`, ...) mean what they did, and a sampled
    `log(0)` or `qnorm(0)` still maps to 0 on a linked parameter.
    This changes the model: the group distribution of t0 is now normal on
    `log(t0 - 0.05)`, so its implied density goes to zero at the bound
    rather than being cut off there. Fits whose t0 posteriors sat on 0.05
    will differ. Sampled-scale numbers (priors, start values, `make_data`
    parameter vectors, `get_pars(map = FALSE)`) for linked parameters are on
    the shifted scale; `summary()` of a design prints each link. Designs,
    priors and fits made before this change keep the transform stored in
    them and map exactly as before.

-   **One sampler:** ensemble updates support standard, diagonal-gamma,
    factor, SEM and infinite-factor hierarchies without nuisance parameters.
    Shear, collapsed moves and subject HMC additionally require standard or
    diagonal-gamma without group designs, nuisance parameters or a joint
    likelihood. Standard models with covariates on `mu` use the ensemble and
    particle updates. The `emc2.*` options that
    switched them off or changed their schedules (`emc2.ensemble*`,
    `emc2.group_move*`, `emc2.shear_*`, `emc2.collapsed_move`,
    `emc2.subject_hmc*`, `emc2.alloc_*`, `emc2.independent_chains`,
    `emc2.sample_rebuilds`, `emc2.burn_mode_start`, `emc2.laplace_refresh`,
    `emc2.laplace_conditional`, `emc2.pooled_candidates`,
    `emc2.iwish_rcond_tol`, `emc2.epsilon_bounds`,
    `emc2.epsilon_scale_bounds`) are removed and now have no effect. Worker, profiling and failure-policy
    options are unchanged.

-   **Collapsed transport surrogate:** empirical quadratic fits now remove
    linear slopes in directions with no positive fitted curvature. Mode-based
    and older saved fits receive the same safeguard before use. This prevents
    the collapsed Gaussian integral from growing without bound along a flat
    direction. Sanitized base fits are cached across iterations, with pooled
    coefficients re-applied at the current value.
    A fractional allocation that selects zero rounds skips collapsed-context
    setup. Each chain keeps its own collapsed-round allocation instead of
    inheriting the largest rate selected by another chain.
    During burn and adaptation, a proposal rebuild with too few collapsed
    moves to measure retains the last measured movement score within the same
    stage. This avoids treating an empty block as fresh evidence that an
    unmeasured move is promising; production-stage movement continues to pool
    over the sample stage.

-   **Shear transport at declared bounds:** subjects sitting near a
    declared lower bound (for example `t0` in the LBA) are transported with a
    map that keeps the bound's half-space, so candidates near the bound are
    no longer lost to the likelihood floor; the collapsed group move uses the
    same walls. The map keeps the move's exact Jacobian, composition and
    inversion. The benchmark harness gained a fixed-length continuation
    protocol (`benchmarks/subject-hmc/continue_arm.R`) that scores arms by
    effective sample size per CPU second on the slowest group direction from
    a shared checkpoint over several seeds.

-   **Checkpoint continuation:** production options absent from older saved
    checkpoints inherit the package defaults when a fit is resumed.

-   **`make_random_effects()` defaults:** the default spread
    (`variance_proportion`) now applies to each design cell rather than to
    each regression coefficient, so a cell reached through a contrast
    (intercept plus contrast) no longer gets the sum of both spreads. Draws
    are checked against the likelihood's exclusive bounds; a parameter
    exactly on a bound was previously accepted and then gave every trial the
    likelihood floor.

-   **Likelihood of binned response times:** with `rt_resolution` each
    response time now stands for its bin instead of the bin's lower edge.
    Race models use one midpoint evaluation away from support boundaries and
    integrate within the fractional bin and its onset transition (using a
    $C^1$ smoothstep blend and floor-aware low-order quadrature), so `t0` can
    enter the earliest bin without discarding its remaining probability.
    DDM/BOU, stop-signal, and LogicalRules likelihoods now apply the same
    clipped-bin integration in their dedicated kernels. Bins are clipped to
    censoring and truncation bounds; timed mixtures, missing responses and
    correlated races use the same bin convention.
    Local Erlang guess clocks retain their support from true zero, including
    the guess mass before the decision onset in bins that straddle `t0`.
    An RT exactly at a grid-aligned upper bound belongs to the preceding bin.
    This removes the first-order bias of evaluating at the lower edge. The
    bin width is stored on the data, and
    `score_held_out()` checks against it.

-   **DDM start-point and non-decision-time variability:** the Wiener density
    with `SZ` or `st0` is now a fixed-node integral instead of an adaptive
    cubature. The cubature's tolerance was absolute, so small densities (fast
    responses, wide start-point ranges, small boundaries) were off by up to a
    nat or more, and its decision to subdivide made the log likelihood jump in
    the parameters. The start point is integrated in closed form term by term
    at short decision times and by a fixed rule over the large-time series
    otherwise; `st0` uses a fixed rule whose nodes follow the process time
    scale. The density is now smooth in every parameter and close to a
    high-precision reference. With `rt_resolution` and `st0 > 0` a binned
    response takes the bin-midpoint density plus a curvature correction for
    each side of the bin, and integrates a side exactly where it lies on the
    density's leading edge (near onset and near the end of the onset window),
    where the plain midpoint biased `st0` and `t0`. The DDM CDF used for
    censoring and truncation uses the same fixed rules in place of its own
    adaptive cubature, which had the same absolute-tolerance error.

-   **R-hat stopping rules per selection:** `mean_gd` and `max_gd` now apply
    to each element of `selection` separately, so the default burn rule needs
    a mean R-hat below 1.1 for the subject parameters and, separately, for the
    group means. Previously the R-hats were pooled, so the far more numerous
    subject parameters could pass burn with badly mixed group means. A vector
    named by selection sets thresholds per element and leaves unnamed elements
    without that rule, e.g. `mean_gd = c(alpha = 1.1), max_gd = c(mu = 1.2)`.
    Progress messages report each rule.

-   **Laplace-conditional efficient proposal:** each chain keeps a Gaussian
    approximation of every subject's likelihood (the Hessian at its burn or
    refresh mode, which the Laplace fits already computed and discarded). Every
    iteration it is combined with the current group mean and covariance, so the
    efficient proposal is the subject's conditional at the current group state
    instead of a regression fitted to earlier draws, which lagged behind the
    chain.

-   **Default sampler:** every hierarchical fit now runs one architecture:
    the ensemble group update, shear rounds, subject Hamiltonian Monte Carlo
    and the collapsed group move, with one shared allocation that sets each
    move's rate (including each subject's particle-step rate) from its
    measured mixing per unit of the chain's wall time. Hierarchies a move
    does not support keep the others.

-   **Shear group update:** the ensemble group update is followed by joint moves of the group level and every subject, in which
    each subject keeps its standardised position in its own conditional
    posterior (a well-informed subject barely moves). The
    ensemble draws the group given the subjects; the shear moves them
    together; alternating the two is interweaving (Yu & Meng 2011). Mean steps
    follow the directions in which the ensemble's group mean is predicted to
    mix slowly, and scale steps rescale one parameter's group sd with the
    subjects' deviations (the hierarchical funnel). Each step is taken with
    the probability that the ensemble alone mixes slowly in it, so a
    well-informed hierarchy gets few steps; each direction and scale adapts
    its own step size, and every step is accepted on the posterior ratio
    times the map's exact Jacobian, screened first on the subjects' Gaussian
    likelihood approximations (delayed acceptance). Under the standard
    hierarchy the Huang-Wand auxiliary `a_half` is integrated out of each
    step and redrawn after it, so a weakly informed group sd is no longer
    held in place by it. Standard (no group designs) and diagonal-gamma
    hierarchies; others keep the ensemble update.

-   **Bug fixes (hierarchy Gibbs steps):** factor and SEM loading rows that
    mix free and fixed entries now regress on the response net of the fixed
    entries. The default factor constraint (diagonal fixed at 1) has such rows
    whenever there are two or more factors, so those loadings were biased.
    SEM K, G and B rows get the same fix, and each B/G row's residual is
    refreshed from the complete row. Per-column SEM prior variances
    (`lambda_var`, `K_var`, `B_var`, `G_var`) are now aligned by column in
    the Gibbs step, prior sampling and bridge sampling. The infinite-factor
    shrinkage update refreshes its cumulative products after every `delta`
    draw (wrong from the third factor on).
-   **SEM:** cyclic structural models (`B_mat`) are rejected; the conjugate
    B update is exact only for recursive models. Bridge sampling and the
    group information criteria now score each subject under its
    covariate-conditional prior, the density the sampler targets. Population
    summaries (`Sigma`, `mu_implied`) include the K x A G cross terms, and
    `mu_implied` no longer drops the K contribution. Efficient proposals
    condition on each subject's covariate-conditional mean on every backend.
-   **Sampler robustness:** pooled-parameter candidate likelihoods follow the
    failure policy (classified, counted under a "pooled" source, retried row
    by row, aborted under `strict`); NA/NaN/+Inf are no longer silently
    zero density. Laplace diagnostics always report a status, treat infinite
    or missing Pareto k as poor, and drop stale curvature when a refresh fails.
    A non-positive-definite group covariance is a rejection, never evaluated
    through a repaired factor. Covariance regularisation is equivariant to
    parameter units. Window covariances are no longer shrunk towards a
    diagonal (Schafer-Strimmer): the shrinkage removed most correlations at
    p of about 70, so every production rebuild of the chain proposal cut
    subject move rates by more than half and never recovered.
    Production proposal rebuilds count
    transitions rather than stored draws, and the schedule is documented as
    adaptation at increasingly rare times (Chimisov, Latuszynski & Roberts
    2018), not diminishing adaptation.

-   **Sampler proposal fixes:** the chain proposal is centred on the plain
    posterior mean of its window (it was a likelihood-weighted mean, which
    degenerates to about one draw in high dimension), with the plain window
    covariance. The efficient proposal's mean now follows the current group mean.
    Burn's Laplace proposals use the widened prior from the first block, with
    a mean that matches their covariance, and the Laplace fit diagnostic is
    rechecked at every refresh. Particle mixture weights adapt to squared jump
    per likelihood evaluation, and the particle count to the ESS of the
    independence step. The ensemble pool's defensive tail scale is
    1 + 3 / sqrt(p). Failed candidate likelihoods in the particle step count
    as zero density. With a correlated prior on mu, `type = "diagonal-gamma"`
    now draws mu jointly.

-   **New (pooled parameters):** `design(pooled = c("t0"))` declares
    coefficients that every subject shares and that are still sampled: one
    value per iteration, drawn from the summed likelihood of all subjects
    under a normal prior taken from the named `pmean`/`psd` entries of
    `prior()`. They sit between `constants` and hierarchical parameters.
    `sampled_pars()` leaves them out unless `include_pooled = TRUE`, the
    group-level samplers and the ensemble update never see them, and
    `get_pars()` returns them with `mu` and with every subject's `alpha`.
    Each iteration adds one move that shifts the pooled value, every subject
    and the group mean together, at about 5-10% extra time per iteration.
    Not available with `type = "single"`, joint models, custom likelihoods,
    `group_design`, nuisance samplers or `marginalise`; bridge sampling and
    Bayes factors are not yet supported. The move's slope and step come from
    the summed curvature of the subject likelihoods (alpha-alpha block
    included), refreshed once per warm-up block, and the collapsed group move
    conditions its subject approximations on the current pooled value. On a
    simulated 40-subject RDM with a truly shared `t0`, the pooled fit reached
    the stopping rule as fast as the hierarchical fit.
-   **New default (ensemble group update):** every hierarchical fit now
    updates the group level with an ensemble MCMC step (Neal 2011;
    Shestopaloff & Neal 2013). Each sweep draws a pool of candidate states per
    subject from a fixed proposal built from that subject's posterior, then
    alternates the existing group Gibbs step with reselecting each subject's
    state from its pool. No likelihood is evaluated after the pool is built.
    The update is exact for the standard, diagonal-gamma, factor,
    infinite-factor and SEM hierarchies (checked against conjugate and Stan
    references) and is skipped for fits with nuisance parameters and for
    single-subject fits. On simulated RDM fits it gave
    5-28x the effective samples per second on the worst group quantity and
    about 10x on subject parameters, at 1.25-1.65x the time per production
    iteration; warm-up is shorter, so total fit time rose 5-46%. On the real
    forstmann data and the N-back witness it converged where the previous
    sampler did not.
-   **Changed (production proposals):** particle-step proposals are rebuilt
    during the `sample` stage from all production draws so far, once the
    stage holds 100, 200, 400, ... iterations (diminishing adaptation, Roberts
    & Rosenthal 2007). A subject whose warm-up window had not converged no
    longer keeps an off-centre proposal for the whole fit. Particle mixture
    weights and scales stay fixed during `sample` (they previously kept
    adapting at a constant rate).
-   **Fix (factor, infinite-factor and SEM hierarchies):** the subject step
    integrates the factor scores out, but the group step then conditioned on
    scores drawn for the previous subject states. Each latent-factor Gibbs
    step now redraws the scores first (partially collapsed Gibbs; van Dyk &
    Park 2008).
-   **Fix (SEM subject prior):** every subject was given the population
    moments at the average covariate, which biased group SDs and correlations
    whenever covariates varied (z up to 27 against a Stan reference). The
    subject prior now conditions on each subject's own covariates.
-   **Fix (adaptation stopping):** `adapt` no longer loops or stops with a
    misleading message when conditional proposals cannot be built, and the
    default `min_unique` is 150 whether or not `stop_criteria` is supplied.
-   **Removed (group moves):** the interweaving group move (legacy, regularized
    adaptive Metropolis and robust adaptive Metropolis variants), its options
    and `group_move_diagnostics()`. The ensemble update
    replaces it: at realistic dimension a rigid group translation gained
    nothing over the Gibbs step.

-   **Fix (failed AM candidates):** group-move candidate likelihood failures
    are retried on the master after worker errors, recorded through the shared
    failure policy, and treated as zero-density rejections unless strict mode
    requests an abort.

-   **Fix (group-move warmup resume):** group-move schema checks no longer
    compare the prior's auxiliary `design` attribute, whose model closures are
    reconstructed when custom-kernel pointers are restored between stages.
    Numeric prior fields and other prior attributes remain checked. Existing
    version-1 warmup checkpoints with the old signature are normalized when
    resumed.
-   **Truncating defective omissions:** `design(TC = list(filter_defective =
    TRUE))` declares that data were filtered to observed responses. A finite
    `UT` then also truncates a defective model's never-finish outcome, so the
    likelihood normalises by `P(LT <= RT <= UT)` and conditions on a response
    being observed; omissions left in the data can only come from
    `pContaminant`, and `make_emc()` rejects them otherwise. The default
    (`FALSE`) keeps the outcome in the retained sample space,
    `P(LT <= RT <= UT) + P(RT = Inf)`. `make_data()` and `predict()` follow the
    design's setting (`filter_defective = NULL`), or an explicit argument.
-   **Fix (truncated defective races):** an upper-censored trial with unknown
    response and `UC < UT < Inf` counted the finite mass above `UT`, giving
    probabilities above one. The correlated BAwL pair, drift-factor and
    finishing-time copula routes left the never-finish outcome out of their
    finite-`UT` normaliser while the `rho = 0` route kept it, scored retained
    omissions as `min_ll`, and (drift factor) mixed both conventions; every
    route now applies the declared one. For correlated fits to data with
    omissions filtered out, `filter_defective = TRUE` reproduces the previous
    likelihood of every trial on a correlated route exactly (trials in
    `rho = 0` cells now match them). Defective LogicalRules races refuse the
    option.
-   **Fix (correlated BAwL likelihood spikes):** the exact two-racer route
    builds each pair value by subtracting bivariate normal CDF corners, whose
    errors are absolute. Deep in a racer's tail it returned cancellation
    residuals, and a truncated trial divided one by another: with a small
    winner `sv` and strong negative `rho`, ordinary trials scored up to +344
    log-likelihood where the true value floors, and the region absorbed the
    sampler. Every exact pair value now carries an error bound. A trial whose
    numerator or truncation normaliser does not clear its bound is retried
    with tvpack corners, then on a numeric route that works in log space
    wherever the mass is (it previously scanned a fixed +/-12 SD window).
    Typical posterior draws are unchanged apart from far-tail trials, which
    move by up to 2e-3.
-   **Fix (normal upper tails):** under the fast normal CDF,
    `pnorm_std(x, lower = FALSE)` was computed as `1 - Phi(x)`, which has no
    relative precision beyond `x ~ 5` and is zero beyond `x ~ 8.3`. It now
    uses `Phi(-x)`. This affects the natural-space helpers in
    `wald_functions.h` (normal intervals, stop-loss terms, the point-start
    Wald CDF, e.g. -3% at an early RT), the lognormal survivor and the
    correlated BAwL numeric route. Results for `x <= 0` are unchanged.

-   **BOU performance:** The default bounded-OU solver grid is now 256 cells
    with a 1 ms target step and five drift/start quadrature nodes. BOU now
    anchors decay at the interval midpoint by default, integrating start-point
    variability with an exact uniform seed and requiring only drift-node
    marches. The previous start-point anchor remains available with
    `BOU(anchor = "start")`; fits using the two anchors are different models.

-   **Breaking:** BTAwL no longer accepts `chart = "endpoint"`; sample
    `tau`/`tau_t` directly.
-   **New:** Ballistic models now expose a trailing `eta` parameter on the
    unbounded identity scale (default `0`), excluded from `p_types_canonical`.
    For the supported LBA, clock-free uncorrelated BAwL, BAwD/BAwDp, BAwF,
    BAwR, and BTAwL constructors, `eta` warps operational accumulation time;
    `eta = 0` is the exact parent model. The CDF/survivor are composed in
    operational time, densities include the time-warp Jacobian, `t0` remains
    additive, and simulation applies the inverse warp. BAwL clock/correlated
    variants, logical-rules, and non-ballistic models reject `eta`.

## Performance

-   Particle workers now receive only the current random effects and
    subject-specific population means for the subjects assigned to them.
    Population covariance factors are still broadcast once per worker, but the
    former full-subject-state broadcast no longer multiplies traffic by the
    worker count. The covariance retained as the factor cache's validation
    reference is reused rather than serialised a second time. Dynamic
    longest-processing-time scheduling and per-subject reproducibility are
    unchanged.

-   Preallocated sample histories are now filled through native slice writers.
    This preserves the public matrix/array representation while avoiding R's
    copy-on-modify duplication of every accumulated iteration on each write.

-   Clean-template workers now recycle every 50 iterations by default; legacy
    workers forked from the chain retain the conservative 10-iteration default.
    `options(emc2.worker_recycle = )` continues to override either value.

-   Set `options(emc2.sampler_profile = TRUE)` to attach per-iteration sampler
    phase timings and worker-message sizes as the `sampler_profile` attribute of
    a chain's sample store. Profiling message sizes performs extra serialisation
    work and is intended for architecture benchmarks, not production fits.

-   The particle step no longer re-factorises its proposal covariances on every iteration. `chains_var`/`eff_var` are constant for a whole block and the group covariance is shared by all subjects within an iteration, so each is now decomposed once instead of `n_subjects * n_proposals` times per iteration. The saving scales with the cube of the number of parameters: negligible for small models, roughly 5 ms per iteration at 50 parameters with 30 subjects, and ~21 ms at 100 parameters.

-   Likelihood calls now hand the C++ mapper compressed design matrices instead of materialising a full-length copy of every design on every call. Log-likelihoods are bit-identical.

-   Each chain now forks its workers once per *block* and keeps them alive across its iterations, instead of forking a fresh set on every iteration. Subjects are partitioned longest-first from their measured times and rebalanced each iteration, rather than being split into the contiguous equal-*count* chunks `mc.preschedule` uses, and a chain that finishes its block early releases its cores to the chains still running, which grow their pools onto them. On a real 112-subject, 17-parameter RDMSWTN fit (3 chains x 8 cores, 40-iteration blocks, 3 interleaved repetitions) this cuts the median block from 37.3 s to 32.4 s, a 13% saving, and puts worker utilisation at 94%; the remainder is master-side dispatch (3.9% of an iteration) and residual subject imbalance (2.0%). Where a pool cannot be started -- Windows, no `mkfifo`, or a single worker -- the previous `mcmapply` path runs unchanged.

    Workers talk to their chain over named pipes. A persistent *socket* cluster is the obvious alternative and is slower than the fork it replaces: on the same probe, `mclapply` fork/join costs 33.5 ms per round, `makeForkCluster` 44.1 ms, and named pipes 1.0 ms. Forked workers inherit the block-constant state, so only the group-level draw and the per-subject bookkeeping cross the pipe each iteration -- about 5 kB on a small model, growing with the square of the number of parameters because the group covariance travels with it (about 90 kB at 40 parameters).

    Each *subject* carries its own L'Ecuyer stream, in the start points as well as in the particle step, so a fit is now reproducible from a fixed seed regardless of how many cores it runs on, which the previous path was not -- it derived its children's streams from `mc.cores`. This is also what makes it safe for a chain to take over a finished sibling's cores mid-block. Fitting no longer leaves the calling session switched to L'Ecuyer-CMRG. Draws from an existing fit are not reproduced exactly by this version; the two are statistically equivalent, not identical.

-   Each chain's workers are now re-forked every 10 iterations rather than
    living for a whole block. A freshly forked worker shares all its pages with
    its chain, but R's garbage collector writes to the header of every object it
    marks, and each of those writes un-shares a page permanently -- so over a
    long block a worker converges on a private copy of the chain's heap. On a
    250 MB heap with 8 workers, total physical memory rose from 0.28 GB (the
    previous per-iteration forks, which died before collecting anything) to
    0.75 GB after 10 iterations and 3.8 GB after 200. Re-forking on a period
    keeps most of the sharing: on that heap it holds the peak to 1.0 GB, a 72%
    saving. A re-fork costs a roughly fixed ~0.19 s, so what it adds to a block
    depends on what an iteration costs -- about 2% on a fit at ~800 ms an
    iteration, more on a small model whose whole block is only seconds long.
    Tune with `options(emc2.worker_recycle = )`; 0 disables it. Results are
    unaffected -- each subject's RNG stream lives in the master, so replacing
    the workers cannot move a draw.

-   The pool's start handshake backs off from 0.2 ms instead of sleeping a flat
    10 ms per worker, which it paid serially on every start, grow and recycle.

## Bug fixes

-   Fixed some maths that was wrong about the ECDF plot for SBC

-   `c_log_likelihood_race` could overflow the stack for models with more than
    64 parameter columns. Several of its branches stage one accumulator row in a
    64-entry buffer, and the guard that kept those in bounds had been widened to
    1024 without widening them. All of the fixed buffers in the race and
    logical-rules likelihoods are now `thread_local` vectors sized to the model,
    so there is no column limit at all -- which matters because a trend model
    extends the parameter-type count one name at a time. This also takes
    `c_log_likelihood_race`'s stack frame from 280 kB back to 3 kB.

-   A `TC$guess_window` (or `TC$w_outlier`) narrower than the observed RTs is now
    an error rather than being applied anyway. The uniform guess density is zero
    outside its window, so trials past the edge were receiving a guess component
    they should not have had -- for `w_outlier = 0.1` on two responses, every RT
    beyond 5 s. The derived default cannot hit this: it is built from the same
    bounds the RTs were already checked against.

-   The ROU equilibrium parameterization now expresses `theta = v / k` and
    `chi = s * sqrt(tk)` in physical units. Threshold effects in `B` therefore
    remain identifiable when `theta` and `chi` are shared across conditions.
-   Joint and trend samplers now align component metadata with the sampled
    parameter names before building proposal covariance blocks. This prevents
    preburn failures caused by logical component indices whose length differed
    from the group covariance dimension.

## New features
-   REXG now conditions each ex-Gaussian process time on being strictly
    positive before adding its accumulator-specific `t0` shift. R and C++
    densities, CDFs, survivors, and simulation now agree on the support
    `rt > t0` and include the zero-truncation normalizer.

-   `RDMSWTNcorr()` and `RDMSWTN_TTcorr()` gain `correlate`, which selects what
    `rho` actually correlates. The default `"times"` is the existing Gaussian
    copula on the two participating finishing times. The new `"drifts"` makes
    `rho` the correlation of the *between-trial drift draws*, exactly as in
    `BAwLcorr()`: rates come from an equicorrelated normal (truncated to the
    positive orthant under `posdrift = TRUE`) and the accumulators then race
    independently on the drawn rates. It requires `sv > 0` on the correlated
    rows, since at `sv = 0` there are no draws to correlate.

    The two are not nested and neither is a limit of the other, so the setting
    has to be chosen deliberately. `"times"` remains the default, and every
    existing correlated RDMSWTN fit is unaffected.

    Both models now share one implementation. `src/drift_factor.h` owns the
    one-factor construction
    `v_q = v + sign(rho) sv sqrt(|rho|) z`, `sv_q = sv sqrt(1 - |rho|)`, the
    joint positive-orthant weight, and the correlated-draw simulator; the
    adaptive Gauss-Hermite driver that integrates the factor out (including
    the truncation normaliser evaluated inside the factor integral, the
    contaminant split, and the orthant denominator) is now parameterised by
    the model's drift mean/SD column positions rather than hard-coded to
    BAwL's. What is *not* shared is BAwL's exact two-racer rectangle kernel
    and its fused node evaluator: both rely on the fixed-time survivor being
    affine in the drift, which is true of a ballistic accumulator and false of
    a Wald one, so a Wald kernel always takes the generic shared-factor route.
    BAwLcorr log-likelihoods are bit-identical to the previous version across
    all of its routes.

    Because each shared-factor node is a numerically integrated kernel on that
    route, its Gauss-Hermite node count is raised above BAwL's for `|rho| >
    0.6` (28 rather than 12); measured against an independent 160-node
    reference on a 24-trial design this takes the worst case from 1.3e-2 to
    6.2e-4 in total log-likelihood.

    One behaviour change falls out of the sharing: `BAwLcorr()` now also
    rejects `sv = 0` on a correlated row, where it previously accepted a
    specification in which `rho` could have no effect.

-   `RDMSWTNcorr()` and `RDMSWTN_TTcorr()` now accept `posdrift = FALSE` with
    `sv = 0`, so a Gaussian copula can couple the finishing times of
    unrestricted-drift accumulators. A negative mean rate makes the marginal
    defective (`F(Inf) = p < 1`); the copula construction is unchanged and
    stays exact, because uniforms at or above the marginal plateau represent
    the atom at infinity. `rho` therefore also couples the intrinsic
    omissions: the probability that neither of a coupled pair ever finishes is
    `Phi2(qnorm(p1), qnorm(p2); rho)` instead of `(1-p1)(1-p2)`. `sv > 0` on a
    correlated row still errors under `posdrift = FALSE`; that combination is
    reserved for a correlated-drift model.

-   **Breaking:** `RDMSWTN()`, `RDMSWTN_TT()` and `LogicalRulesRDMSWTN()` with
    `posdrift = FALSE` now sample `v` on the natural scale
    (`transform = "identity"`, bounds `(-Inf, Inf)`, default `1`) rather than
    the log scale, so the negative mean rates the unrestricted-drift model is
    defined over are actually reachable. `posdrift = TRUE` is unchanged.
    Sampled `v` values from earlier `IO` fits are on the old log scale and are
    not comparable; refit, or map them with `exp()`.

-   New `pGuess` parameter: a uniform "guess" (outlier) contaminant, in the
    spirit of Ratcliff & Tuerlinckx (2002) and HDDM's `w_outlier`. Unlike
    `pContaminant`, which is a Bernoulli **omission** mixture contributing mass
    only at `rt = Inf`, `pGuess` mixes a flat density directly into observed RT
    densities, giving fast and slow outliers a likelihood floor. It is available
    on the race families and, for the first time, on the DDM family (`DDM`,
    `DDMGNG`, `BOU`), which also gain `pContaminant`.

    The two contaminants are nested rather than competing, so neither can push
    the other out of `[0, 1]`:

        P(omission) = pContaminant
        P(guess)    = (1 - pContaminant) * pGuess
        P(process)  = (1 - pContaminant) * (1 - pGuess)

    At `pGuess = 0` every likelihood is bit-for-bit identical to before, so
    `pContaminant`'s behaviour is unchanged.

    Two things worth knowing. First, both contaminants are proportions among
    **retained** trials, not generated ones: they are applied after truncation
    renormalisation, matching how `make_missing()` contaminates after the
    truncation cut. Second, the guess window defaults to the effective
    truncation/censoring window `[max(LT, LC), min(UC, UT)]`, widening to
    `max(5, floor(max(rt)) + 1)` when that has no finite upper edge -- so, unlike
    HDDM's fixed 5 s, a 7 s outlier still gets mixture protection. Override it
    with `TC$guess_window` (or `TC$w_outlier` for the HDDM spelling), and see
    `?resolve_guess_window`. `make_data()`/`make_missing()` simulate guesses via
    the `pGuess` and `guess_window` arguments.

    `pGuess` is applied only in the compiled likelihoods; `design()` raises an
    error rather than silently ignoring a free `pGuess` on a model with no
    compiled likelihood.

-   Finished general trends implementation (stay tuned; tutorial still on the way)


# EMC2 3.3.0

## Bug fixes

-   Addressed a major issue in `map = TRUE`, so that now mapping of population level parameters is now done correctly. See issue #119

## New features

-   More general support for `group_design` (tutorial on the way), also using `map = TRUE`

-   Start of a more general trends implementation (stay tuned; tutorial still on the way)

-   Some more SEM/FA functionality added i.e. `rotate_loadings`

# EMC2 3.2.1

## Bug fix

-   Added a warning that in hierarchical models, using `map = TRUE` in functions like `map = TRUE` or `map = TRUE` does not return the population-level marginal mean and variances on the original scale for group-level parameters. See issue #119

# EMC2 3.2.0

## New features

-   group_design specification (tutorial coming up)

-   trends specification (tutorial also coming up)

-   broader continuous covariates support with `map = TRUE`

-   fMRI joint modelling (tutorial on <https://osf.io/preprints/psyarxiv/rhfk3_v1>)

-   Made changes to how lR was used in compressed likelihood for race models. Run `update2version` to reuse older samples

## Bug fixes

-   Made legend including or excluding more flexible for `plot_cdf`, `plot_density` and `plot_stat`

# EMC2 3.1.1

## New features

-   added thin to fit/run_emc which can either be set to TRUE to automatically thin based on ESS, or on a numeric to only keep 1/x samples

-   added probit/SDT model for bimanual choices

## Bug fixes

-   Rare bug in sampling removed

-   Small bug fixes in plot_data to make it more flexible

-   cleared up argumentation of run_emc/fit

# EMC2 3.1.0

## New features

-   model_averaging function, which allows you to compare evidence for an effect across a set of models

## Bug fixes

-   Small hotfix in which creating proposals and the start of burn would sometimes fail for large number of subjects

-   Patched up old error in which model bounds weren't considered in data generation

-   Fixed error in which compare_subject would return IC for whole dataset for every subject.

# EMC2 3.0.0

## New features

-   IMPORTANT: to keep your old samples compatible with current EMC2, run update2version(<name of old samples>)

-   IMPORTANT: Design and prior are now also S3 methods with their own S3 classes, see EMC2 paper

-   plot_fit is deprecated and has branched of into plot_density, plot_cdf and plot_stat

-   sampled_p_vector is deprecated and is now named sampled_pars

-   Added a design_plot function, which makes a plot of the proposed accumulation process

-   Sampling is completely reworked. Adaptive tuning of the number of particles and more stable convergence

## Bug Fixes

-   Fixed rare case where conditional MVN would break

-   Fixed bug in predict on joint models

-   Suppressed unwanted print statements in DDM estimation/prediction

-   Added more checks to a wide array of functions to ensure proper input format

# EMC2 2.1.0

## New features

-   Added a website with vignettes, changelog and a reference

-   Added `run_sbc()` function to perform simulation-based calibration for a design

-   Added `prior_help()` to get more information on the prior for a certain `type`

-   Changed DDM implementation, which is faster and more accurate

-   Bridge sampling now also works for `type = "blocked"`

## Bug Fixes

-   Made bridge sampling for inverse-gamma and inverse-wishart more robust

-   Made `prior()` function work more generally
