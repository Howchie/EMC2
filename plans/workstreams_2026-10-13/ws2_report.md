# WS2 rebuild overhead

Baseline: `0a5773dd`. Implementation: `/data/work/EMC2_ws2`, branch
`ws2/rebuild-overhead`. Private libraries and raw results live in
`/data/work/EMC2_ws2_lib`. The original working tree, allocator, particle step
and in-iteration ensemble code are untouched.

## Measurement scope

Further runs use **nine CPUs, 3 chains × 3 cores** for sampler continuations,
with affinity `23-31`. The user requested relative comparisons alongside other
work and superseded the shared-lock rule. Wall results therefore include
shared-machine contention; they are not an idle-machine or 27-core comparison.

The component screen runs baseline and candidate back to back twice with seed
713 and checks complete results and caller RNG using `identical()`. The
installed-package profiler measures exclusive wall, master CPU, worker CPU,
actual dispatch width and idle core-seconds. Worker CPU is returned with each
fork's output: R's child CPU counters can attribute a completed fork to the
following component when it is reaped late. Sampling CPU from persistent
workers remains partial.

The first source-overlay screen defined its own unwrapped `auto_mclapply`;
its candidate worker CPU was not captured. Only its paired **wall times** and
identity checks are used here. CPU and idle figures come from the separately
installed-package profiles. The committed screen routes both implementations
through the same telemetry wrapper.

The fixture named `lba7r_after_preburn.rds` has `init = FALSE` and no sampler
history. It was skipped by the screen. `prepare.R` creates a separate 20-draw
checkpoint for an early-burn profile. Completed Forstmann, T30 and N-back fits
retain sample draws only, so their rebuild profiles use those saved histories.
LBA adapt and sample-entry profiles use the genuine post-adapt checkpoint.
The sample-next profile advances the saved transition counter to its next
rebuild gate; no intervening draws are simulated. Other fixtures' sample-entry
and sample-next routes are workload profiles on saved sample histories, not
reconstructed first and second production blocks. Earlier adapt histories on
those completed fits cannot be recovered.

## Implemented changes

| Change | Classification | Exactness and lifetime |
|---|---|---|
| Omit unused pooled alpha history for independent chains; reuse covariance checks | Redundant work | Same means, covariance fallback and proposal output; no retained cache |
| Claim jobs in descending subject cost for Laplace refresh and fit checks | Parallel across subjects | Each job uses its recorded RNG stream; original result order restored |
| Claim atlas node fits by subject cost | Parallel across subject/node pairs | Same deterministic fits and assembly order |
| Claim orbit jobs; pack one context per participating chain | Parallel across coordinates; redundant packing | Same likelihood evaluations, grids and widths |
| Fork sufficiently wide collapsed regressions | Parallel across subjects | Same quadratic fits; narrow/small workloads stay serial |
| Slice diagnostic arrays directly | Redundant indexing and copying | Same dimensions, names, values and subject grouping |
| Distribute large R-hat/ESS groups | Parallel across diagnostic groups | Existing multivariate grouping retained; custom functions and small histories stay serial |

Rebuild forks inherit their inputs and exit after the dispatch. They do not
send new model definitions to the persistent sampling workers. The queue cache
contains only one shared eight-byte claim map and a busy flag; it does not
retain sampler objects or callbacks. Nested calls use the original apply path.

`block_proposals` now includes planning, allocation and proposal bookkeeping,
excluding the separately recorded efficient proposal preparation.
`block_eff_proposals` covers efficient preparation in burn, adapt and sample.
The two fields are disjoint and together cover `add_proposals`; a skipped
sample rebuild records zero in both.

## Timing evidence

One installed-package LBA burn rebuild on the post-adapt state fell from
**79.154 to 65.147 seconds** at nine cores, a 14.007-second difference (17.7%).
This is one profile pair, not a whole-fit claim. Repeated component screens
separately support the Laplace, orbit, covariance and regression improvements.
The original under-ten-second target is not demonstrated. The measured
remaining work includes roughly 188 worker CPU seconds in Laplace and 96 in
orbit measurements alone, plus the allocator's roughly five seconds of wall.
These phases already perform substantial parallel computation.

### Installed rebuild profiles

| History / route | Rebuild wall before / after (s) | Idle before / after (core s) |
|---|---:|---:|
| forstmann_burn | 2.331 / 2.255 | 11.57 / 10.38 |
| forstmann_sample_entry_route | 2.883 / 2.780 | 15.69 / 14.77 |
| forstmann_sample_next_route | 3.067 / 3.096 | 17.11 / 16.88 |
| lba_adapt_state_adapt | 70.259 / 62.346 | 263.46 / 190.57 |
| lba_adapt_state_burn | 79.154 / 65.147 | 339.43 / 213.25 |
| lba_preburn20_burn | 64.719 / 57.058 | 230.98 / 161.91 |
| lba_sample_first | 141.726 / 108.140 | 564.66 / 266.97 |
| lba_sample_next | 66.088 / 59.232 | 233.21 / 167.36 |
| nback31_burn | 50.987 / 42.991 | 193.19 / 119.38 |
| nback31_sample_entry_route | 306.335 / 252.538 | 1103.50 / 641.81 |
| nback31_sample_next_route | 298.484 / 240.151 | 1030.81 / 513.04 |
| T30_burn | 17.168 / 16.349 | 119.52 / 113.27 |
| T30_sample_entry_route | 13.547 / 12.201 | 83.81 / 73.89 |
| T30_sample_next_route | 13.390 / 12.646 | 83.30 / 73.69 |

Two additional Forstmann sample-next pairs measured 3.202 / 3.259 and
3.347 / 3.355 seconds (before / after). Their means are 3.275 / 3.307 seconds,
a roughly 1% slowdown on this small workload. The queue avoids small job
dispatches, but the combined changes do not improve every fixture and route.
These extra pairs use the final `deliver8` library; their timings are exported
in `small-model-repeats.csv`.

Timing and memory exports are in [results](../../benchmarks/rebuild-overhead/results/).
`components.csv` contains exclusive wall, master/worker CPU, idle core-seconds,
busy cores and dispatch widths for each fixture/stage route. Component idle
figures use `max(0, 9 × wall - master_cpu - worker_cpu)` and include time lost
to contention or unmeasured fork overhead. Exclusive rows can be summed without
counting nested work twice. Profile pairs are one pass each; the following
component table contains means of two paired repeats.

| Fixture | Component | Stage | Wall before / after (s) |
|---|---|---|---:|
| forstmann_full_fit.rds | `.emc_check_likelihood_fits` | refresh | 0.946 / 0.865 |
| forstmann_full_fit.rds | `.emc_collapsed_fits` | refresh | 0.039 / 0.059 |
| forstmann_full_fit.rds | `.emc_laplace_refresh` | refresh | 1.120 / 1.087 |
| forstmann_full_fit.rds | `.emc_shear_atlas_add` | refresh | 0.000 / 0.000 |
| forstmann_full_fit.rds | `.emc_shear_orbit_steps` | refresh | 0.580 / 0.609 |
| forstmann_full_fit.rds | `create_chain_proposals` | adapt | 0.093 / 0.049 |
| forstmann_full_fit.rds | `create_chain_proposals` | burn | 0.121 / 0.048 |
| forstmann_full_fit.rds | `create_chain_proposals` | sample | 0.088 / 0.048 |
| forstmann_full_fit.rds | `split_along_dim` | refresh | 0.005 / 0.001 |
| lba7r_after_adapt.rds | `.emc_check_likelihood_fits` | refresh | 4.289 / 3.937 |
| lba7r_after_adapt.rds | `.emc_collapsed_fits` | refresh | 0.951 / 0.386 |
| lba7r_after_adapt.rds | `.emc_laplace_refresh` | refresh | 35.283 / 31.560 |
| lba7r_after_adapt.rds | `.emc_shear_atlas_add` | refresh | 10.917 / 10.968 |
| lba7r_after_adapt.rds | `.emc_shear_orbit_steps` | refresh | 19.829 / 17.300 |
| lba7r_after_adapt.rds | `create_chain_proposals` | adapt | 0.787 / 0.388 |
| lba7r_after_adapt.rds | `create_chain_proposals` | burn | 0.856 / 0.439 |
| lba7r_after_adapt.rds | `create_chain_proposals` | sample | 0.772 / 0.373 |
| lba7r_after_adapt.rds | `split_along_dim` | refresh | 0.075 / 0.019 |
| lba7r_sample_stage.rds | `.emc_check_likelihood_fits` | refresh | 4.637 / 4.151 |
| lba7r_sample_stage.rds | `.emc_collapsed_fits` | refresh | 0.824 / 0.442 |
| lba7r_sample_stage.rds | `.emc_laplace_refresh` | refresh | 35.160 / 32.996 |
| lba7r_sample_stage.rds | `.emc_shear_atlas_add` | refresh | 0.003 / 0.002 |
| lba7r_sample_stage.rds | `.emc_shear_orbit_steps` | refresh | 18.833 / 16.588 |
| lba7r_sample_stage.rds | `create_chain_proposals` | adapt | 0.665 / 0.364 |
| lba7r_sample_stage.rds | `create_chain_proposals` | burn | 0.657 / 0.366 |
| lba7r_sample_stage.rds | `create_chain_proposals` | sample | 0.676 / 0.370 |
| lba7r_sample_stage.rds | `split_along_dim` | refresh | 0.063 / 0.012 |
| nback31_fullc_1_emc.rds | `.emc_check_likelihood_fits` | refresh | 2.284 / 2.167 |
| nback31_fullc_1_emc.rds | `.emc_collapsed_fits` | refresh | 0.455 / 0.237 |
| nback31_fullc_1_emc.rds | `.emc_laplace_refresh` | refresh | 13.614 / 12.149 |
| nback31_fullc_1_emc.rds | `.emc_shear_atlas_add` | refresh | 50.117 / 46.852 |
| nback31_fullc_1_emc.rds | `.emc_shear_orbit_steps` | refresh | 6.425 / 5.542 |
| nback31_fullc_1_emc.rds | `create_chain_proposals` | adapt | 0.226 / 0.114 |
| nback31_fullc_1_emc.rds | `create_chain_proposals` | burn | 0.283 / 0.322 |
| nback31_fullc_1_emc.rds | `create_chain_proposals` | sample | 0.208 / 0.117 |
| nback31_fullc_1_emc.rds | `split_along_dim` | refresh | 0.066 / 0.009 |
| T30_seed1_emc.rds | `.emc_check_likelihood_fits` | refresh | 1.038 / 1.065 |
| T30_seed1_emc.rds | `.emc_collapsed_fits` | refresh | 0.211 / 0.150 |
| T30_seed1_emc.rds | `.emc_laplace_refresh` | refresh | 1.984 / 1.851 |
| T30_seed1_emc.rds | `.emc_shear_atlas_add` | refresh | 0.002 / 0.002 |
| T30_seed1_emc.rds | `.emc_shear_orbit_steps` | refresh | 1.588 / 1.567 |
| T30_seed1_emc.rds | `create_chain_proposals` | adapt | 0.235 / 0.124 |
| T30_seed1_emc.rds | `create_chain_proposals` | burn | 0.257 / 0.120 |
| T30_seed1_emc.rds | `create_chain_proposals` | sample | 0.227 / 0.122 |
| T30_seed1_emc.rds | `split_along_dim` | refresh | 0.040 / 0.007 |

### Stopping diagnostics

Means of two runs in separately installed libraries. LBA sample-stage checks
fall from **8.068 to 2.779 seconds** (65.6%); T30 and N-back also benefit.
The post-adapt check shows no gain because most of its work stays below the
parallel workload threshold. Small Forstmann summaries remain serial.

| Fit | Stopping check wall before / after (s) | Idle before / after (core s) |
|---|---:|---:|
| forstmann_full_fit.rds | 0.375 / 0.367 | 3.00 / 2.94 |
| lba7r_after_adapt.rds | 7.243 / 7.316 | 57.76 / 58.25 |
| lba7r_sample_stage.rds | 8.067 / 2.779 | 64.54 / 13.00 |
| nback31_fullc_1_emc.rds | 2.877 / 2.155 | 23.02 / 13.01 |
| T30_seed1_emc.rds | 2.156 / 1.333 | 17.25 / 8.25 |

## Identity and memory

- Initial full validation: **121 cases passed**, no failures, errors or skips.
- After adding fit-check/atlas queues and the regression workload guard:
  **48 relevant cases passed**, no failures, errors or skips.
- All **90 paired component comparisons** across five initialized saved fits
  returned bit-identical complete results and caller RNG state.
- Two burn blocks on both Forstmann and LBA returned bit-identical draws,
  chain RNG streams, proposal matrices and recorded plan fields. These use the
  existing tests' fixed allocator so wall-time probes cannot select different
  move counts. Ordinary timing-based allocation can vary between runs; default
  whole-fit bit identity is not claimed.
- The array-slice diagnostic comparisons returned identical complete stopping
  results and RNG states in two repeats on each of five fits, checking `mu`,
  `sigma2` and `alpha` with R-hat and ESS criteria.
- Eleven Forstmann burn blocks, with 3 chains × 3 cores: the three persistent
  templates each kept exactly constant RSS (about **229 MiB**, zero range).
  Sampling workers were recreated each block and remained around 226–233 MiB.
  Master RSS rose from about 268 to 380 MiB during profiling. This establishes
  stable worker memory over eleven blocks, not bounded master memory.

The final diagnostic change passes **18 targeted cases**, including old/new
R-hat and ESS, multivariate R-hat to fifteen digits, argument forwarding to
custom functions, R-hat rules, adaptation stopping, proposal schedules and
sample storage. Ten comparisons on shared inputs return strictly identical
complete stopping results, sampler state and RNG. Separately deserialized
results recreate model/prior function environments, so raw `identical()` across
those files compares environment identity as well as values; their numerical
fields match exactly. `decisions.R` provides the shared-input strict check.

Twelve repeated sample-stage diagnostic checks keep the master's file descriptor
count at **six** throughout. Their fork workers exit after each dispatch.
After warm-up their mean RSS fluctuates around 560–579 MiB, with no monotonic
increase; the final mean is 573 MiB. Collection occurs between checks. This
addresses repeated diagnostic dispatch lifetime separately from the eleven
sampler-block template test. Complete burn-plus-adapt convergence runs
have not been timed; saved-state differences cannot be multiplied into a
whole-fit saving without the actual rebuild count and trajectory. No claim is
made that the original under-ten-second target has been met.

## Overlapping proposal rebuilding with sampling (not implemented)

A proposal built from the previous completed block could be computed while the
next block samples. The saving is at most the rebuild wall that actually fits
inside that next sampling block. The plan's roughly 30 seconds of preparation
in a roughly 120-second sampling block suggests a ceiling near 30 seconds per
eligible block, or roughly 20% of combined time, before resource contention.
Convergence checks, concatenation and checkpoint writes still need their
current-block inputs. This is a ceiling, not a measured prediction: the
current sampler already uses nearly all available cores during likelihood
work, and giving background rebuilding cores can slow sampling.

Each applied proposal must depend only on completed past draws and independent
adaptation randomness, remain fixed throughout its sampling block, and be
applied on a predetermined block schedule. Every resulting kernel must retain
its target invariant distribution. Those conditions alone do not prove
convergence of adaptive MCMC. Production would need diminishing changes to the
kernel and containment, or a finite adaptation period followed by a frozen, irreducible and aperiodic
kernel. See [Roberts and Rosenthal (2007), Theorem 2](https://www2.stat.duke.edu/homeweb/scs/Courses/Stat376/Papers/AdaptiveMC/RobRosCouplingErgod2007JAP.pdf)
for the diminishing-adaptation and containment argument. Increasing the interval between rebuilds alone does not prove
diminishing adaptation: non-vanishing changes could persist at the sparse
rebuild times. Discrete frame and curve gates need particular attention.

The current Laplace checks advance subject RNG streams. An asynchronous
rebuild cannot overwrite streams that sampling has already advanced; it would
need separate adaptation streams or a carefully defined stream allocation.
Checkpointing must preserve the pending proposal's history cutoff and application
boundary. Broadcasting mutable sampler objects while a block updates them is
unsafe; inputs need an immutable snapshot. Extra snapshot memory, worker
contention, queue cancellation and reproducible restart behaviour would all
need measurement. This option changes draws and adaptation timing and has
therefore been left for the user to decide.

## HMC pilot failures (added from WS1)

The HMC calibration pilots in `.emc_hmc_calibrate` used to turn every
likelihood error into a NULL. A broken kernel or worker therefore looked like a
rejected step, and the subject's step size and HMC price were set from it,
with no count and no warning. WS1's replays hit this when saved states were
replayed without `restore_custom_kernel_pointers()`, because forked pilots
raced on Rcpp's build cache. Pilot errors are now classified in the same way as
a production HMC transition's. Under `strict`, a non-numerical failure stops
the fit. Each pilot returns its tally, and the master counts it under the
`pilot` source and announces it, because a forked pilot's own counters and
warnings end with the process. The change does not alter sampling when no
likelihood fails. `test-subject-hmc.R` checks it with a mocked failing
likelihood.

## Commits and reproducibility

| Commit | Scope |
|---|---|
| `bcdec264` | Complete proposal timing and attributed component profiler |
| `302eca82` | Unused pooled histories and duplicate covariance checks |
| `c656be6f` | Shared rebuild claims and balanced Laplace refresh |
| `91ad2797` | Balanced orbit jobs and per-chain context reuse |
| `eea9d0dc` | Parallel collapsed regressions |
| `5b17e9eb` | Direct diagnostic array slices |
| `79019246` | Balanced likelihood fit checks |
| `6a3f6d7f` | Balanced atlas node fits |
| `46e56829` | Regression workload/width guard |
| `3372d5df` | Large stopping diagnostics over unchanged MCMC groups |
| `983e4120` | HMC pilot likelihood failures counted and under the failure policy (from WS1) |

The benchmark [README](../../benchmarks/rebuild-overhead/README.md) describes
commands and interpretation. Raw inputs, sampler states and full result
comparisons remain outside the repository. Source snapshots and loaded-library
paths are recorded in the private logs. Rebuild profiles use `deliver6`; that
contains every rebuild change listed above. The final stopping-diagnostic
comparison uses `deliver8`, which additionally contains `3372d5df`.
