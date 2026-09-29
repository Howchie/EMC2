# Handoff: burn stage mixes poorly on an LBA fit; adapt and sample are fine (EMC2 build 2026-09-29 00:11)

## Symptom
The first fit of a new LBA on the new build ran 2,000 burn iterations and still left burn
with mean R-hat 1.115. An RDM fit to the same data with the same sampler needed 500.
At the end of burn the chains were separated along a population-wide speed-caution
direction. Over the last 500 burn iterations the chain means differed by 2-3
within-chain SDs:

| parameter | chain means |
|---|---|
| v_LoadLow | 2.74 / 2.64 / 2.54 |
| B_LoadLow | 0.11 / 0.00 / -0.01 |

t0_LoadLow, v.fam and v.e_low were also separated. Subject-level chain gaps reached 5
within-chain SDs (e.g. 52176 t0_LoadHigh -2.46 / -2.23 / -2.48, SD 0.05). The mode start
was clean: 135/135 BFGS runs converged and 135/135 Hessians were usable.

## Key finding: the problem is confined to burn

| window | subject acceptance | standardised sq jump / coord / iter | ensemble moved / sweep | max subject chain gap (within-SD) |
|---|---:|---:|---:|---:|
| burn, last 500 (read live at 02:05) | 0.17-0.78 | 0.08-0.26 | ~0.09 (cumulative to that point) | 5.0 |
| adapt + sample (snapshot) | 0.69-0.85 | 0.6-1.5 | **0.68** over sample sweeps 2301-2400 | 1.0 |
| RDM comparator (converged, sample) | | | 0.57 | 0.3 |

The ensemble figure for the sample sweeps is chain 1's cumulative `moved` counter,
216.0 at 2,300 sweeps and 284.4 at 2,400. The jump sizes are standardised by the pooled
SD of the window.

Once adapt builds the chain proposals and `ensemble_q`, both the subject kernel and the
ensemble move work very well on this model. Sample-stage max R-hat is falling:
1.90, 1.33, 1.21 at 100/200/300 draws, as the burn-era chain offsets wash out. The cost is
the burn stage: 2,000 iterations at low ensemble efficiency, with the subject moves an
order of magnitude smaller than later.

On the frozen snapshot, `ensemble_q` covers the posterior about as well as the RDM's does.
Mahalanobis m2 median 18.1 vs 16.1 (ideal 17); max variance ratio median 2.0 vs 1.7. So
nothing points at the sample-stage q_s.

## Questions for the package
1. What proposals do the subject kernel and the ensemble pools use during burn, and when
   are they refreshed? The burn-era ensemble acceptance (~9%) and jump sizes are far below
   post-adapt. Is burn relying on the mode-start Laplace covariance, and is it refreshed
   often enough ("Laplace refresh: 135/135" appears 4 times in the log)?
2. Why is the burn exit criterion so slow here when adapt fixes the mixing within ~200
   iterations? Consider building the empirical/ensemble proposals earlier in burn, or
   leaving burn on a weaker criterion.
3. Record the ensemble move rate and weight ESS per stage (or per window), not only
   cumulatively; that gap is why this took a while to see. The live checkpoint also
   drops burn draws once sampling starts, so burn-stage diagnostics are unrecoverable
   after the fact.

## Files
Snapshot and scripts are in this folder (`plans/convergence/`). Other paths are relative
to the project root `/data/work/PM/NirvanaHons_Nback`; run the model scripts from there,
since they need `CleanData.RData` and the custom `.cpp` kernels.
- Snapshot at 02:18: `samples_lba_rates_d45_snapshot_0218.RData`. It holds adapt 200 +
  sample 200 only (burn draws were already discarded by the checkpoint) and includes
  `ensemble_q`, cumulative `ensemble_stats` and `burn_mode_stats`. Log at the same time:
  `lba_rates_d45_snapshot_0218.log`.
- Run script (exact model): `Fits/run_scripts/lba_rates_t0R_vs2_mt_d45.R`. BAwL
  uncorrelated, posdrift TRUE; custom v kernel `match_joint_trend_rates.cpp`; A ~ 1;
  sv ~ lM with the sv intercept fixed at 1; t0 ~ 0 + Load + repeatTrial;
  diagonal-gamma shape 2; 3 chains × 4 cores; 45 subjects. Rerunning it reproduces the
  burn (default MT RNG), and a burn-stage checkpoint copy would recover the burn draws.
- Live fit: `Samples/samples_control_Exp1_lba_rates_t0R_vs2_mt_d45.RData` and
  `Fits/lba_rates_t0R_vs2_mt_d45.log`.
- Comparator (same sampler, converged):
  `Fits/Match_Control_Exp1_rdm7s_mass_massR_noA_rates_t0R_UT0_vs2_mt_d45_ens.RData`.
- Diagnostics (absolute paths): `subj_moves.R` (per-subject acceptance and jump over the
  stored window), `q_fit.R` (ensemble_q coverage), `subj_stuck.R <samples> <stage>`
  (per-subject chain gaps).
- Unreliable: my offline replay of pool log-likelihoods
  (`.emc_ensemble_ll_candidate` on the saved object) returned -Inf for every draw,
  including for the RDM. That is a calling error, not evidence.

## Other package bugs found this session
- `make_emc(type = "diagonal")` without `par_groups` gives a single full-covariance
  block: `par_group` is all 1, whereas `R/fitting.R:1269` intends `1:n`. The type seems
  to be mapped to "standard" before that branch.
- `predict()` errors ("p matrix columns must include: <par>") when a parameter is
  constant across all posterior draws.

## Outcome (2026-09-29)
- **Burn efficient proposals: kept.** From 200 burn iterations, every burn block builds
  `create_eff_proposals(from_burn = TRUE)`, and the subject kernel runs the adapt-style
  four-component mixture (R/fitting.R `add_proposals`, R/sampling.R `run_stage` /
  `new_particle`).
  - Evidence: `repro_burn_wide_eff_from_200.R`, from the 200-iteration checkpoint
    `burn_pool_200.RData`, with the 3x ensemble widening. Burn exited at mean R-hat 1.033
    after 4 tries; without the change it hit max tries at 1.117.
  - That is one seed on one model. More seeds and fixtures are still to come.
- **3x widening of the burn `ensemble_q`:** tried narrower, reverted, still 3x.
- **Pathfinder and the other mode optimisers** (BFGS, nlminb, L-BFGS-B, Pathfinder):
  evaluated and removed. BFGS stays; see `benchmarks/mode-start/README.md`.
- **Bugs in the list above:**
  - The `type = "diagonal"` claim was wrong as stated. `par_group` is `1:p` and the Gibbs
    step is diagonal. The actual bugs were a full inverse-Wishart start `theta_var` and
    prior sampling that ignored `par_group`; both are fixed.
  - `predict()` with a free parameter constant across draws is fixed.
  - Also fixed: an infinite recursion in `get_data()` / `sampled_pars()` on joint models.
- Data files here (`*.RData`) are local only and git-ignored.
