# EMC2 3.4.0

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
