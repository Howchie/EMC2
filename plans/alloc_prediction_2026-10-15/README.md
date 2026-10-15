# Allocator prediction: evaluate options for predicting a coordinate's tau

You are picking up an evaluation task in the EMC2 sampler (`/data/work/EMC2_dev_oo`,
branch `playground`). Do not commit unless the user asks. Read the header comment of
`R/move_allocation.R` first: it defines every quantity below.

## What the allocator does

Every 100-iteration block (`step_size`, one "try" in the console output), each chain
records per-iteration jumps of its group coordinates (mu along the shear's mean
directions, then the log group sds) and of each subject. At the rebuild it turns them
into a lag-one rho per coordinate, `rho = 1 - ESJD / (2 var)`, where var is pooled over
every chain's window draws. "Lag" always means sampler iterations within the current
block, not tries and not particle proposals. Moves multiply rho by measured factors:
shear rounds `g`, collapsed moves, subject HMC and particle rates. The plan then
minimises `max_k tau_k * cost` with `tau = (1 + rho) / (1 - rho)`, which treats each
coordinate as a single AR(1). Chains choose one shared HMC regime: none, solver rates,
0.5 or 1 for every subject.

## What was found (2026-10-14/15)

1. **The lower-bound bias is fixed in the working tree (uncommitted).** rho was a lower
   bound: ESJD was inflated by 2 sqrt(2/n), about 28% for a block. Every coordinate,
   and every plan's prediction, looked faster than it was. The point estimate is now
   used (`.emc_alloc_rho`).
2. **The solver fix is in the working tree too.** A plan re-solved with subject HMC
   rates held fixed (`r_fixed`, from `.emc_alloc_settle`) could fall back to the
   nothing-chosen plan. That plan runs no HMC, but it kept the group credit for the
   held rates. `.emc_alloc_solve` now treats held HMC rates like held particle rates.
3. **The remaining error is the AR(1) mapping, not other chains.** Pooled and
   within-chain variance agree within 1-14% (`calib2.R`). Lag-1 predicted from one
   block matches the next block without bias once the lower bound is gone. At lags 5-10
   the next block's autocorrelation is 0.09-0.22 higher than AR(1) predicts, so tau is
   understated about 2-3x, most on the slowest coordinates (`predtest.R`):

   | next-block tau, predicted / realized (median) | A: lower bound | B: point (now) | C: two-rate fit, lags 1-3 |
   |---|---|---|---|
   | PM burn | 1/3.1 | 1/2.3 | 1/1.5 |
   | PM screening arms | 1/2.3 to 1/3.8 | 1/1.7 to 1/2.8 | 1/1.4 to 1/1.8 |
   | slowest 10% of coordinates | up to 1/7.7 | up to 1/5.7 | up to 1/2.9 |
   | Forstmann (stationary, fast) | ~1 | ~1.1 | ~1.1 |

   C is a Prony fit from lags 1-3: two decay rates and a weight. It is noisy: the median
   error is about 2x either way. The realized tau of slow coordinates from 100-draw
   blocks is itself an underestimate.
4. **A tried and rejected fix is in `slow_share_experiment.patch` (not applied).** It
   split each group coordinate into a fast part, which shear and collapsed moves act
   on, and a share riding a subject statistic at that statistic's lag-1 s. The share
   was identified from lag-2: `w = (r2 - r1^2) / ((s - r1)^2 + r2 - r1^2)`, and the
   subject plan moved `1 - s` by the subjects' relative movement.
   - Version 1 took s from the subjects' mean or log-sd. That statistic moves with the
     shear, so the ridge was attributed to the subjects, its shear rounds were cut, and
     R-hat reached 2.1 (`arms/slow_1.rds`).
   - Version 2 took s from positions in each subject's conditional,
     `y_s = Q_s (a_s - m_s)`, which the shear holds fixed. The ridge got w = 0 as
     intended, but the PM arm was still worse than without it (`arms/slow2_1.rds`:
     min ESS 5, R-hat 1.64).
   - Either way, tau remains mispredicted. Read the patch for its tests and the shared
     cache pattern; do not resurrect it as is.

## Your task

Evaluate options for predicting each group coordinate's tau under a plan, and recommend
one, with evidence. Candidates include:
- the two-rate fit (C);
- shrinking C's noise, for example by pooling the long-lag excess across chains or
  coordinates of a kind;
- a per-coordinate correction factor taken from the block's own measured lag-k
  autocorrelation;
- anything better you find.

Hard constraints from the user:
- **Within the current block only.** No multi-block history and no forecasting of
  convergence. Burn is where this matters most, and older blocks are stale because the
  allocation itself changes the chain.
- **Moves change only what they act on.** A prediction must say how a plan's moves
  (shear rounds, collapsed moves, subject HMC and particle rates) change the coordinate.
  Shear and collapsed moves hold subjects' positions in their conditionals fixed (see
  `R/shear_move.R`); interweaving with the ensemble update means the subject updates
  and the group moves both matter.
- **No serial work in chain processes, and no warning cop-outs.** If the extra
  computation is more than trivial, it goes in C++.
- **Comments state what and why only.** No narrative, fixture names or numbers in code.

Evaluate first on the saved draws, in the style of `predtest.R`: one block predicts the
next. Score lag-k autocorrelation and tau, and report the slowest coordinates
separately. Only then change the allocator. Any check that needs fits must run one fit
at a time, 3 chains x 9 cores, and only when the user has no fit running (check
`uptime` and ask if unsure).

## Data and scripts here

- `predtest.R` / `predtest.rds`: the block-to-next-block test and its per-coordinate
  results.
- `calib2.R`: pooled versus within-chain variance and predicted versus realized lag-1.
- `cont.R`: a 400-iteration screening continuation of the PM checkpoint. Run it as
  `EMC_LIB=<lib> Rscript cont.R <label> <seed> 400 <out.rds>`; it saves group draws
  and plans.
- `replay2.R`, `replay3.R`: replay the allocator on the checkpoint offline. Run them
  from `/data/work/PM/NirvanaHons_Nback` with `LIB=<lib>`.
- `arms/*.rds`: the screening arms' draws and plans.
  - `old_1`: allocator before bf677e14.
  - `new_1`: bf677e14.
  - `slow_1` / `slow2_1`: the two slow-share versions.

Fixtures:
- PM checkpoint:
  `/data/work/PM/NirvanaHons_Nback/Samples/samples_control_Exp1_frqk7r_blk_gL_sepA_sepL_sepV_sv_t0T_mt.RData`.
  FRQ model with the old `delta` column, 900 burn + 100 adapt iterations, unconverged,
  with a K-v-t0 ridge.
  - Load it with a library built from bf677e14 plus the allocator files, not the
    current tree: 4d629b37 replaced `delta` with `plo`.
  - Saved fits keep shared fields on chain 1 only, so copy them to the other chains.
- Stationary easy fixture: `/data/work/EMC2_workstreams/fixtures/forstmann_full_fit.rds`
  (LNR). Score it with `benchmarks/subject-hmc/continue_arm.R`.

Install into a private library with
`MAKEFLAGS=-j24 nice -n 19 R CMD INSTALL --no-test-load --library=<lib> .`, and delete
`src/alloc_need.o` and `src/RcppExports.o` after changing exports. Tests:
`tests/testthat/test-move-allocation.R`, `test-collapsed-move.R`, `test-shear-move.R` and
`test-subject-hmc.R`, run with `NOT_CRAN=true EMC2_TEST_LEVEL=full`.
