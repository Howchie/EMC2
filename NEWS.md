# EMC2 3.4.0

## Performance

-   The particle step no longer re-factorises its proposal covariances on every iteration. `chains_var`/`eff_var` are constant for a whole block and the group covariance is shared by all subjects within an iteration, so each is now decomposed once instead of `n_subjects * n_proposals` times per iteration. The saving scales with the cube of the number of parameters: negligible for small models, roughly 5 ms per iteration at 50 parameters with 30 subjects, and ~21 ms at 100 parameters.

-   Likelihood calls now hand the C++ mapper compressed design matrices instead of materialising a full-length copy of every design on every call. Log-likelihoods are bit-identical.

-   Each chain now forks its workers once per *block* and keeps them alive across its iterations, instead of forking a fresh set on every iteration. Subjects are partitioned longest-first from their measured times and rebalanced each iteration, rather than being split into the contiguous equal-*count* chunks `mc.preschedule` uses, and a chain that finishes its block early releases its cores to the chains still running, which grow their pools onto them. On a real 112-subject, 17-parameter RDMSWTN fit (3 chains x 8 cores, 40-iteration blocks, 3 interleaved repetitions) this cuts the median block from 37.3 s to 32.4 s, a 13% saving, and puts worker utilisation at 94%; the remainder is master-side dispatch (3.9% of an iteration) and residual subject imbalance (2.0%). Where a pool cannot be started -- Windows, no `mkfifo`, or a single worker -- the previous `mcmapply` path runs unchanged.

    Workers talk to their chain over named pipes. A persistent *socket* cluster is the obvious alternative and is slower than the fork it replaces: on the same probe, `mclapply` fork/join costs 33.5 ms per round, `makeForkCluster` 44.1 ms, and named pipes 1.0 ms. Forked workers inherit the block-constant state, so only the group-level draw and the per-subject bookkeeping cross the pipe each iteration -- about 5 kB on a small model, growing with the square of the number of parameters because the group covariance travels with it (about 90 kB at 40 parameters).

    Each *subject* carries its own L'Ecuyer stream, in the start points as well as in the particle step, so a fit is now reproducible from a fixed seed regardless of how many cores it runs on, which the previous path was not -- it derived its children's streams from `mc.cores`. This is also what makes it safe for a chain to take over a finished sibling's cores mid-block. Fitting no longer leaves the calling session switched to L'Ecuyer-CMRG. Draws from an existing fit are not reproduced exactly by this version; the two are statistically equivalent, not identical.

## Bug fixes

-   Fixed some maths that was wrong about the ECDF plot for SBC

-   The ROU equilibrium parameterization now expresses `theta = v / k` and
    `chi = s * sqrt(tk)` in physical units. Threshold effects in `B` therefore
    remain identifiable when `theta` and `chi` are shared across conditions.

## New features

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
