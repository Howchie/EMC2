# Subject proposal investigation

Starting revision: `abe2c850` (`playground`). Acceptance criteria: maximum
R-hat below 1.1, subject ESS in the hundreds for 1,500 production iterations
per chain on `wide`, and no ESS/second regression (including the worst subject)
on `diag`, `corr`, `strong`, and `tight`. Fits run sequentially, 3 chains by
9 cores. ESS is a diagnostic, not a stopping rule.

## First candidate: lower the local warm-up scale floor

`update_pm_settings()` previously imposed epsilon floors of 0.1 (preburn)
and 0.4 (burn/adapt). These multiply a covariance factor that may still be
on the group scale or reflect separated chains. They therefore prevent the
proposal from reaching a sufficiently small scale for a tight posterior.
`update_epsilon_scale()` also imposed a 0.05 floor at proposal refreshes.
This candidate replaced the local proposal floors with a numerical bound of
0.0001. Independent components retain their existing adaptation floors:
shrinking an off-centre independent proposal can reduce overlap further.
Proposal refreshes preserve valid small scales down to 0.0001.
Production still freezes epsilon and mixture weights; proposal moments retain
their existing expanding-window doubling schedule.

This changes the tuning range of the existing conditional importance sampler,
not its weights or state transition construction. The relevant particle Gibbs
construction is Andrieu, Doucet and Holenstein (2010),
[Particle Markov chain Monte Carlo methods](https://www.stats.ox.ac.uk/~doucet/andrieu_doucet_holenstein_PMCMC.pdf).

Inspection also established that `weighted_moments()` likelihood-weights the
proposal **mean**, but computes an ordinary unweighted covariance. Covariance
double weighting is therefore not a supported explanation of this failure.

## Controlled scale probe

Run `Rscript benchmarks/subject-move/scale_probe.R`. This probes one stationary
update on an exactly known Gaussian posterior; it does not fit the RDM data.
The auxiliary density and local proposal density cancel, so reselection weights
are the target density. The particle count matches 40% of `50 * sqrt(p)`.
The proposal base SD is 0.5. Each row uses 2,000 independent repetitions.

For 24 parameters and posterior SD 0.05, epsilon 0.4 yielded zero moves;
epsilon 0.04 yielded movement in 86.3% of updates, with mean squared jump
0.240 posterior variances per coordinate. Full results are in
`results/scale-probe.csv`. These measurements demonstrate a possible failure
mechanism, not successful repair of the RDM fixture.

## Reproduction

Install into a separate library and check the `EMC2 from` line in each log:

```sh
mkdir -p /tmp/emc-subject-local
R CMD INSTALL --library=/tmp/emc-subject-local .
FIXTURE_ROOT=benchmarks/subject-move/results/scale-local SAVE_EMC=1 SAVE_ALPHA=1 \
  Rscript benchmarks/group-move/fixture_fits.R fit /tmp/emc-subject-local wide baseline 1
```

Copy the existing `fixture-wide.rds` into the result directory first to retain
the original data. `SAVE_EMC=1` now saves each completed stage and logs epsilon
quantiles. Results also record the loaded library and subject R-hat.

## Rejected subject random-walk candidates

Lowering the local scale floor alone completed the wide run with median subject
ESS 17.9, minimum subject ESS 4.1, maximum subject R-hat 2.01, and maximum
group R-hat 1.48. Coordinatewise warm-up Metropolis moves gave median/minimum
subject ESS 13.2/3.4 and maximum subject R-hat 3.12. Continuing those moves
through production gave median/minimum ESS 8.5/3.5 and maximum subject R-hat
2.99. Correlated component block Metropolis moves gave median/minimum ESS
14.6/3.4 and maximum subject R-hat 2.96. These runs still had chains in
different subject-level locations.

## Rejected conditional elliptical slice trial

The published elliptical slice kernel was tested under the Gaussian group prior.
On `wide` it completed 1,500 production draws with subject ESS min/median/max
3.4/7.5/74.9, maximum subject R-hat 2.99, and maximum group R-hat 1.46. It did
not resolve chain separation.

## Incomplete preconditioned MALA trial

A temporary installed build added a MALA move using the group-prior covariance
as its preconditioner, finite differences for the likelihood gradient, an
exact forward/reverse Metropolis-Hastings correction, and warm-up-only scale
adaptation. Its correlated-Gaussian unit test and the required particle,
ensemble, and proposal-schedule tests passed against that installed snapshot.
The wide run exited with code 137 after saving its adapt checkpoint and before
writing production draws, so it has no final ESS or R-hat result. The MALA
implementation is not in the current worktree source.

## Current worktree candidate: regularized empirical proposals

The current source regularizes poorly conditioned empirical proposal
covariances and rescales proposal multipliers per subject when proposal
moments are refreshed. It also adapts a state-centred proposal covariance
during preburn/burn, then freezes it. The code labels this latter update RAM,
but it uses a multi-particle conditional-IS pool statistic rather than the
single-proposal Metropolis acceptance probability in Vihola's published RAM
algorithm; I am treating it as experimental pending review and the wide-fit
diagnostics.

The regularized-covariance build is installed at
`/tmp/emc-subject-ram1`. Its wide baseline completed with subject ESS
min/median/max 3.9/6.0/160.3, maximum subject R-hat 2.19, and maximum group
R-hat 1.27. The high adapt-stage movement did not translate into mixing. The
same snapshot passed `test-particle-invariance`, `test-ensemble-move`,
`test-sample-proposal-schedule`, `test-sampler-guards`, and
`test-subject-proposal-scale`.

## Metric-preconditioned MALA trial

A separate installed snapshot at `/tmp/emc-subject-mala-metric1` uses the
regularized subject covariance as a fixed MALA proposal metric, with the group
Gaussian retained in the target density. Metric and step-scale adaptation end
before production. Its correlated-Gaussian test and the full particle,
ensemble, proposal-schedule, guard, and scale suites have passed. The wide
baseline, wide ensemble, and 13-parameter baseline fixture fits are running
from `results/mala-metric*`. The MALA integration is not currently present in
shared `R/sampling.R`; that file changed during the isolated build, so results
are tied to the installed snapshot until the source edit can be reconciled.
Concurrent-run timings are excluded from comparisons.
