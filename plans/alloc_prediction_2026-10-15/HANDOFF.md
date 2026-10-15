# Allocator handoff (2026-10-15, end of day)

Worktree `/data/work/EMC2_dev_oo`, branch `playground`. Nothing is committed. Read `README.md` (the original
brief) and the header of `R/move_allocation.R` first. Earlier notes: `TRANSPORT_FINDINGS.md` and `eval_kappa/`.

## Where things stand

**HEAD (bf677e14 allocator, "new") is known to be wrong.**
- It rates every coordinate with a lower-bound lag-one rho.
- It looks better only because two errors partly cancel (below).
- Do not revert to it as a fix. The user's view: a correctly scored allocator should never be beaten by
  "turn the allocator off and run HMC every iteration". On the broken builds, that is roughly what the bad
  scoring produces anyway.

**The working tree is "point + κ floor"** (`handoff_2026-10-15/working_tree_point_kappa.patch`):
1. `.emc_alloc_rho` is the point estimate, which replaces `.emc_alloc_rho_lo`.
2. The held-HMC solve fix: the nothing-chosen plan is not the fallback while HMC rates are held above zero.
   This is documented in `README.md`.
3. **The κ floor: written but not built, tested or screened.**
   - It multiplies a group coordinate's tau by κ = (Geyer tau over the stage so far, every chain) /
     (AR(1) tau of the window's whole-iteration lag one), with κ ≥ 1.
   - Code: `.emc_alloc_stage_tau`, `.emc_alloc_whole_rho`, `.emc_alloc_kappa`, `.emc_alloc_eff`.
   - Wiring: unit `kappa` fields, `range_parts`, `scaled_parts`, `unit_need`, `units_need` and `predict`.
   - Tests are appended to `tests/testthat/test-move-allocation.R`.
   - **Expect it to fail as a fix:** the stage tau needs more than about 300 draws before it carries signal
     (see below).

**Other material:**
- The transport-scored allocator is shelved: `transport_allocator.patch`, and
  `eval_transport/move_allocation_transport.R`.
- `R/sampling.R` is back at HEAD; per-move path recording lives only in the transport patch.

## The user's goal and constraints (keep these in view)

- **Get into the posterior as fast as possible.** Burn is about 3 blocks (300 draws), and the default
  sampling target is 10 tries. Long production runs are not the goal.
- **Multi-block information is allowed**, but it is useless unless it carries signal within roughly the
  first 300 draws. Anything that needs more history arrives too late.
- **The allocator should buy whatever stops chains drifting apart.** That may be a particle step, subject
  HMC, a shear round or a collapsed move.
- From the README: moves change only what they act on; no serial work in chain processes (non-trivial work
  goes in C++); no warning cop-outs; comments state what and why only; fits run one at a time on
  3 chains x 9 cores, never while the user is fitting; do not commit unless asked.

## What is established

**1. Lag-one is the wrong score for group coordinates.**
- Lag-one prediction is unbiased, but lag-one does not track tau.
- On PM, the slowest group coordinates have run tau 70–200, while plans predict about 10.
- With pooled variance, rho1 = B/V + (W/V) ρ_w, so the measured lag-one is mostly bounce plus the
  between-chain offset.

**2. Who actually moves the slow directions** (instrumented path recording, run-level attribution;
`attrib.R`, `attrib2.R`):
- **PM ridge** (v_LoadLow, v_LoadHigh, K, sv against t0_LoadHigh): shear carries 56–86% of the persistent
  transport. Subject HMC and particle steps carry about zero, yet take about 70% of compute.
- **T30:** the group update, shear and collapsed moves carry it; particle steps carry ≤ 0.
- **Forstmann (null):** clean.
- Defect identified but **never tested in isolation:** `.emc_alloc_dev_ratio` credits the group for
  subject HMC and particle movement measured at lag one, which is bounce.

**3. Why "point" is a systematic regression.** Three paired seeds, 400-iteration continuations of the PM
checkpoint (`handoff_2026-10-15/arms/`):

| | CPU (s) | ridge ESS per thousand CPU-s | median ESS per thousand CPU-s |
|---|---|---|---|
| new | 10.9–11.1k | 0.72–1.13 | 3.0–4.1 |
| point | 16.5–17.8k | 0.30–0.54 | 1.1–2.5 |

- **Mechanism, from replays of saved states** (`replay_state.R`):
  1. The group's tau is predicted at about 10 when it is about 150, so the minimax target is about 10.
  2. The point estimate (correctly) rates subjects slower, so each chain buys more subject moves to reach
     that target. In one replay, one chain chose every particle step plus HMC 0.28.
  3. The shared particle rate is the maximum over chains, so every chain runs particle rate 0.9–1.0.
  4. The re-solve at that held rate drops HMC to about 0.01–0.13.
- The lower bound had hidden this by making subjects look fast. **More accurate scoring of one input
  exposing a wrong objective is the user's central point: the scoring, not the estimate, is broken.**
- The max rule itself is not the bug. Under the current scoring, alternatives score 3–4x worse
  (`replay_q.R`), because holding a subject's particle rate low where HMC moves it poorly strands it.

**4. Every estimate of group tau from one block fails on the slow tail.** Ratio of predicted to run tau
for the slowest 8 coordinates, on the saved arms:

| estimator | ratio, slowest 8 |
|---|---|
| lag-one AR(1) | 0.12–0.27 |
| κ, Geyer within the block | about 0.3–0.5 |
| transport at L = 10–25 | 0.26–0.48 (best rank, Spearman 0.6–0.9) |
| transport at L ≥ 50 | noise |
| block means across chains | about 0.3, very noisy |

Scripts: `eval_transport/horizon.R`, `scripts/bm.R`.

**5. Stage-so-far Geyer tau** (`scripts/stagecheck.R`, predicting draws 101–400):

| stage length | predicted ÷ realized, slowest 8 |
|---|---|
| 100 draws | 0.20–0.43 |
| 200 draws | 0.14–0.79 |
| 300 draws | 0.48–1.22 (partly in-sample) |

That is too late and too noisy for a 300-draw burn.

**6. "Agreement" scoring** (moves towards the other chains' mean, `scripts/agree.R`):
- Per iteration it equals the existing one-step credit (correlation 0.8–0.99 across coordinates for every
  family).
- Net over a block it is noise: halves of a block don't agree and consecutive blocks don't agree. This holds
  on PM, T30 burn and Forstmann.

**7. Screening.**
- Single-seed 400-iteration arms cannot rank allocators on ridge R-hat: the same shear spending gives
  R-hat 1.2–1.7.
- CPU and ESS per thousand CPU-seconds over 3 paired seeds do separate a 50% effect.
- Use paired seeds (`scripts/cont.R` / `cont_save.R`). Use replays from saved states (`replay_state.R`) to
  compare decisions from an identical state.

## Directions not yet tried (hypotheses, not recommendations)

1. **Remove the group's credit for subject moves** (`.emc_alloc_dev_ratio` / `.emc_alloc_credit`).
   - Point 2 says subject HMC and particle steps do not move the slow group directions. Subject moves would
     then be bought only for the subjects' own coordinates.
   - Cheap to test: a replay from saved states, then 3 paired seeds against new.
2. **Stop the subjects' target being set by an unachievable group target.**
   - Under minimax tau x cost, an underestimated group tau makes the solver push subjects to match it.
   - Possible fixes: subjects bought to their own marginal value only, or a group target that cannot pull
     subject spending.
3. **Baselines to bound the problem** (3 paired seeds each): allocator off with HMC regime 1 every
   iteration, and allocator off with no HMC. The user expects a correct allocator to beat the first.
4. **Multi-block evidence that arrives within about 3 blocks.**
   - For example, pool each family's measured persistent transport across blocks with decay. The machinery
     is in `transport_allocator.patch`: per-family path recording, and `.emc_alloc_transport_shares` for a
     Bartlett lag-sum against the total change.
   - Within one block it ranks families but underpredicts tau; check whether pooling 2–3 blocks stabilises
     the per-move shares on PM and T30 burn.

## Data and tools

**Fixtures:**
- PM checkpoint:
  `/data/work/PM/NirvanaHons_Nback/Samples/samples_control_Exp1_frqk7r_blk_gL_sepA_sepL_sepV_sv_t0T_mt.RData`.
  It has the old `delta` column, so it needs a library built from bf677e14 plus the allocator files. Run
  from `/data/work/PM/NirvanaHons_Nback`.
- T30 and Forstmann: `benchmarks/subject-hmc/fit_arm.R` with `CORES_PER_CHAIN=9 SAVE_EMC=1`. For path
  attribution, use `trim = FALSE` (`scripts/fit_arm_full.R`).

**Scratch** (session scratchpad
`/data/work/EMC2_dev_oo/.Rtmp/claude-1000/-data-work-EMC2-dev-oo/19c4c6df-25fb-4858-8ef0-fdc4d9b6acea/scratchpad`;
it may be cleaned):
- Saved chain states for replays (80–90 MB each): `pair_{new,point}_{2,3}_emc.rds`,
  `arm_pointsave_1_emc.rds`, and the instrumented fits `full_T30_emc.rds` and `full_forst_emc.rds`.
- Libraries:
  - `libbf`: bf677e14 plus the transport-version allocator. Transport is switched off with
    `options(emc2.alloc_transport = FALSE)`; an inert `emc2.alloc_q_rule` scaffold is also compiled in.
  - "new" and "old": in the other session's scratch,
    `/data/work/EMC2_dev_oo/.Rtmp/claude-1000/-data-work-EMC2-dev-oo/d3f36cc0-f8ff-45dc-a75e-971049367f47/scratchpad/scr/lib{new,old}`.
  - A git worktree at bf677e14 lives in `scratchpad/wt_bf`; remove it with `git worktree remove` when done.

**Commands:**
- Install: `MAKEFLAGS=-j24 nice -n 19 R CMD INSTALL --no-test-load --library=<lib> .` After export changes,
  delete `src/alloc_need.o` and `src/RcppExports.o`.
- Tests: `NOT_CRAN=true EMC2_TEST_LEVEL=full` on `test-move-allocation.R`, `test-collapsed-move.R`,
  `test-shear-move.R` and `test-subject-hmc.R`. The current tree's κ tests have not been run.
