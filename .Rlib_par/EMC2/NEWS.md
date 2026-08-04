# EMC2 3.4.0

## Performance

-   The particle step no longer re-factorises its proposal covariances on every iteration. `chains_var`/`eff_var` are constant for a whole block and the group covariance is shared by all subjects within an iteration, so each is now decomposed once instead of `n_subjects * n_proposals` times per iteration. The saving scales with the cube of the number of parameters: negligible for small models, roughly 5 ms per iteration at 50 parameters with 30 subjects, and ~21 ms at 100 parameters.

-   Likelihood calls now hand the C++ mapper compressed design matrices instead of materialising a full-length copy of every design on every call. Log-likelihoods are bit-identical.

-   Experimental and **off by default**: `options(emc2.worker_pool = TRUE)` forks each chain's workers once per *block* and keeps them alive across its iterations, instead of forking a fresh set on every iteration. Subjects are then partitioned longest-first from their measured times and rebalanced each iteration, rather than being split into the contiguous equal-*count* chunks `mc.preschedule` uses. On a real 112-subject, 17-parameter RDMSWTN fit (3 chains x 8 cores, 40-iteration blocks, 3 interleaved repetitions) it cuts the median block from 37.3 s to 32.9 s, a 12% saving, and does the same work for 14% less CPU (481 down to 415 core-seconds) -- the removed fork and the removed chunk imbalance respectively.

    Workers talk to their chain over named pipes. A persistent *socket* cluster is the obvious alternative and is slower than the fork it replaces: on the same probe, `mclapply` fork/join costs 33.5 ms per round, `makeForkCluster` 44.1 ms, and named pipes 1.0 ms. Forked workers inherit the block-constant state, so only the group-level draw and the per-subject bookkeeping (about 5 kB) cross the pipe each iteration.

    Each *subject* carries its own L'Ecuyer stream, so a pooled fit is reproducible from a fixed seed regardless of how many cores it runs on -- which the default path is not, since it derives its children's streams from `mc.cores`. A pooled fit does not reproduce the default path's draws; the two are statistically equivalent, not identical. With `emc2.dynamic_cores` also set, cores released by a finished chain grow the pool instead of feeding `mc.cores`, which is safe here for the same reason.

-   Experimental and **off by default**: `options(emc2.dynamic_cores = TRUE)` lets a chain that finishes its block early hand its cores to the chains still running, instead of leaving them idle until the block ends. Chains keep their independent per-block batching, so unlike `emc2.flat_parallel` this adds no synchronisation: a chain publishes a marker when it exits and the others re-read, once per iteration, how many are still alive. On a real 112-subject RDMSWTN fit it halves the straggle (idle core-seconds 20-36 down to 16-20) but is only ~1% faster overall (35.3s to 35.0s per 40-iteration block), because a straggler absorbs donated cores at roughly 20% efficiency. It is off by default because `mclapply` derives its children's L'Ecuyer streams from `mc.cores`, so a timing-dependent core count makes a fit irreproducible from a fixed seed.

-   Experimental and **off by default**: `options(emc2.flat_parallel = TRUE)` schedules all (chain, subject) particle updates of an iteration as one flat, longest-first, load-balanced task list instead of nesting a fork over subjects inside a fork over chains. It was written to stop early-finishing chains from leaving cores idle, but on a standard 24-subject, 3-chain RDM fit it is 1.3-2.1x *slower* than the existing path at every core count tested, so it is not recommended. The reason is structural and is documented in `R/parallel_pool.R`: the Gibbs step forces a barrier every iteration, and an iteration is too little work to amortise per-iteration dispatch. The same measurements show the existing path does not benefit from more cores than there are chains either.

    It is kept because it may still pay off where a single subject's likelihood is expensive enough to dwarf the barrier (the PDE-backed models) or where subject costs are very uneven. It also gives the *same* draws regardless of core count, because each (chain, subject) carries its own L'Ecuyer RNG stream; the default path does not have that property.

## Bug fixes

-   Fixed some maths that was wrong about the ECDF plot for SBC

## New features

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
