# Allocator replay and calibration harness

Tools for auditing the move allocator (`R/move_allocation.R`,
`src/alloc_need.cpp`). They check whether its predictions hold, whether its
plans are stable, how close its solver gets to its own model's optimum, and
what a rebuild costs. The specification is in
`plans/workstreams_2026-10-13/ws1_allocator_spec.md` and the results are in
`ws1_report.md` beside it.

Every script loads the package from `EMC_LIB`. Set `R_LIBS` to the same path
so the forked workers load it too, and check the `EMC2 from` banner. Any run
over 4 cores, or any timed run, goes under
`/data/work/EMC2_workstreams/bench_lock.sh` with `taskset -c 0-26`. Don't
reinstall a library, or edit a script, that a queued run reads: queue frozen
copies instead.

## Pipeline

1. **Calibration run**:
   `calibrate.R <ckpt.rds> <out_dir> <seed>`.
   - Continues a sample-stage checkpoint for `ITERS` production iterations
     (default 500), applying the checkpoint's recorded `production_options`.
   - Forces a rebuild at the start, then rebuilds at the offsets in
     `REBUILDS` (default `100,300`), giving blocks of 100 and 200.
   - Sources `trace.R`, which records each rebuild's inputs and outputs
     without changing the sampler:
     - `pre_<k>.rds`: the state `add_proposals()` received, including the
       block's per-iteration trace;
     - `alloc_<k>.rds`: per-chain measurements, units, records and the solve.
   - Saves the draws in `run.rds`.
   - `ARM_R` names an R file defining `arm_options()` and/or `arm_hook(emc)`,
     for counterfactual arms. `arm_force.R` is one: with `FORCE=up` it puts
     rounds on unbought shear coordinates and HMC on unbought subjects; with
     `FORCE=down` it divides the bought rates. The recorded predictions are
     rescaled to the forced plan.
2. **Calibration tables**: `analyse.R <run_dir> [label]`.
   - Compares every rebuild's predicted lag-one rho per coordinate with what
     the next block realised. The realised value is measured three ways:
     - the allocator's ESJD estimator, against the rebuild's pooled variance;
     - the standard lag-one ACF;
     - IACT, from `posterior::ess_basic`.
   - Prediction columns:
     - `pred` is the prediction as solved, using the confidence bound on the
       rest of the iteration;
     - `pred_pt` replaces that bound with its point estimate;
     - `pred0` is the prediction with nothing bought.
   - `rest_in` and `rest_real` give the rest of the iteration's factor as
     measured and as realised, which separates move-factor error from drift
     in the rest.
   - Writes `calib.rds`.
3. **Plans and the binding coordinate**: `plans.R <run_dir>`. For each
   rebuild and chain it prints tau, regime, rounds, collapsed rate, HMC
   count and particle rate, plus the coordinate that binds the minimax and
   how that coordinate realised.
4. **Decision stability**: `bootstrap.R <run_dir> <k> [n_rep] [variants]`.
   - Needs `k >= 2`, so that `pre_<k>` carries a full block trace.
   - First replays the rebuild without resampling and checks that it
     reproduces the recorded plans. It also checks that the reconstruction
     of the allocator's inputs from the trace is the identity.
   - Then moving-block resamples (block length 10) one of four input sets:
     - the rest-of-iteration jumps;
     - the moves' evidence;
     - the costs;
     - all three together.
   - Reruns the rebuild after each resample and writes `boot_<k>.rds`. The
     spread of tau, regime, rounds and HMC counts across replicates is the
     plan's sampling noise, attributed to each input set.
5. **Exact inputs and cost** (any script that replays a saved state must
   pass it through `restore_custom_kernel_pointers()` first; otherwise every
   likelihood call reloads a custom kernel, which inflates times and, when
   pilots fork, races on Rcpp's build cache so that calls fail at random):
   - `capture_post.R <pre_k.rds> <post_k.rds> [n_cores]` captures the exact
     state `.emc_alloc_rebuild()` receives, including seed and options.
   - `prof_rebuild.R <post.rds> [reps]` times and profiles one rebuild on that
     state and saves its plans, so two libraries can be compared plan for
     plan.
6. **Solver quality**:
   - `oracle.R` evaluates the allocator's own objective for any plan. Group
     coordinates are credited at the subjects' rates the plan chooses, with
     the particle floor applied.
   - `ws1_oracle()` grid-searches small unit sets.
   - `oracle_gap.R [n] [seed]` compares the solver's plan with the oracle's
     on random problems: one shear coordinate and two subjects. It reports
     the gap, and whether the solver's claimed objective is the objective of
     the plan it returns.
7. **Pricing equivalence**: `verify_batch.R <out.rds> <alloc_*.rds...>`. It
   wraps `.emc_alloc_units_need` and, on every call made while solving the
   captured rebuilds, compares it with the per-unit `.emc_alloc_unit_need`.
   With `REF=1` it saves the plans, so a run on another library can check
   them.

## Fixtures

The calibration fixtures are sample-stage checkpoints:

- `forstmann_full_fit`, `lbaE_fit` and `lba7r_sample_stage` in
  `/data/work/EMC2_workstreams/fixtures/`;
- T30 and nback31 fullc checkpoints. Their paths are in `ws1_report.md`.

A 500-iteration run at 3 chains × 9 cores takes 3–25 minutes, depending on
the fixture.
