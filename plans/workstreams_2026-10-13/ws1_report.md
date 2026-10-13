# WS1 report: is the move allocator making the right decisions?

This covers the branch `ws1/allocator` (from `0a5773dd`). The model as written down from the code is in `ws1_allocator_spec.md`. The harness is in `benchmarks/allocator/` (see its README).

Read the evidence with these limits in mind:
- Every fixture ran with one seed. All numbers are indicative.
- Calibration runs used 3 chains × 9 cores; counterfactual runs used 3 chains × 3 cores.
- The calibration runs used a library built before the `h` fix (commit 1 below). That fix changes plans only for units without `h_n`; production subject units always carry `h_n`, which masked the bug there.

## Verdict (one page)

**Is the model adequate? Yes, as a decision rule, with known noise.**
- **Group coordinates are well predicted.** These are the shear mean directions and log sds, the coordinates that usually bind the minimax.
  - Point predictions are within ±0.06 of the next block's realised lag-one rho on four of five fixtures. The exception is lba7r log sds, at −0.10, where the rest of the iteration itself sped up from 0.55 to 0.49.
  - They correlate 0.54–0.93 with the realisation, 0.85–0.93 on lbaE, nback31 and lba7r means.
  - The predicted move factors match the realised ones within the noise of the yardstick (see 3).
- **Subject coordinates are noisy one at a time.** A single subject coordinate's realised rho over a 100–200-draw block correlates only 0.01–0.23 with its value in the next block. The allocator's prediction correlates 0.34–0.46 with the next block, better than the last measurement does.
  - The subject "worst realised" coordinates are mostly noise: the top 20 regress from about 0.8 to about 0.3 in the next block.
  - The minimax over thousands of subject coordinates is therefore mostly protected by the confidence bound, which the model applies to them.
- **Lag one does not systematically flatter the moves.** IACT/tau runs 1.03–1.51, set by each fixture's own slow modes (lbaE highest). Coordinates with bought moves show no larger ratio than unbought ones.
- **The solver is near-optimal for its own model.** On 40 random problems against a brute-force oracle: median gap 0.999, 90th percentile 1.019, worst 1.117. The worst cases come from the uniform regime floor (see Remaining).

**What it gets wrong**
1. A partial-match bug: the plan without HMC priced HMC from the deviation factors. Fixed in commit 1.
2. No record of what each plan predicted, so it could not be checked. Fixed in commit 2.
3. **Plans vary across resamples of one block's evidence, and the moves' own evidence drives that, not costs or the rest of the iteration.** Under the block bootstrap:
   - the objective's tau typically varies 5–15% (up to 34% on one T30 chain);
   - a family's rounds vary 30–60% (for example, T30 log-sd rounds 0.9–2.9 against a plan of 1.9);
   - the number of HMC subjects varies 2–15 against a plan of 10;
   - the binding coordinate is often a different one (near-ties in a flat minimax);
   - regimes never flipped.
   Smoothing, if wanted, belongs on the move-factor evidence, not on the window of the rest of the iteration.
4. **In a live fit, plans varied run to run only through the pilots' measured CPU time, which fed the HMC cost. `145c30cd` prices pilots by the rows they evaluate instead.** An earlier draft reported a second, unexplained source in multi-core replays. It was an artifact of my replay scripts: they replayed saved states without restoring the custom kernels' pointers. Every pilot likelihood call then rebuilt the kernel through Rcpp's on-disk cache, and concurrent forks raced on that cache. Calls failed at random (`load(index_file)`: zero-byte input; duplicate definitions in the generated source), and the pilots' `tryCatch` turned each failure into an HMC abort. With the pointers restored, a four-core replay is bit-identical across runs, pilots and plans included. Live fits restore pointers on load. The pilots' silent conversion of likelihood errors into aborts still hides failures like this one; that belongs to the HMC code (WS2).
5. **The evidence cap read the evidence against the wrong block (spec S1). Fixed in `a488164e`.** The cap was 2 n_ev / B, with B the next block's length. On the doubling production schedule it allowed about 1.5× the rate the evidence was gathered at, not 2×. I had called this latent. It isn't:
   - a shear cap above the exploration floor binds on at most one coordinate per plan, but that is often the coordinate the minimax binds on;
   - T30's collapsed cap binds too;
   - replayed on captured rebuilds, the fix lowers the plans' modelled tau on 6 of 15 chains (section 7).
6. **S2:** the confidence bound applies to the rest of the iteration only, which makes the plan's predictions as solved optimistic: by up to 0.05 on group coordinates and by 0.07–0.12 on subjects, whose windows are not pooled. This is deliberate protection against max-of-noise. Whether it changes which plan wins was not tested; an earlier draft said it does not, which was unverified.

**What was fixed or made faster:** eight commits, listed in section 1. Six of them change no plan, either bit for bit or to rounding. The S1 fix changes plans by design, and the h fix (commit 1) changes the optimistic pass. On the four-core lba7r input, main-process time per rebuild went from 3.45 s to an estimated 1.6 s:
- solve: 1.9 → 0.75 s;
- measurement: 1.0 → 0.29 s.

The activation fix then takes a further quarter off lba7r's solve (one-core input: 1.25 → 0.94 s).

T30 was the costliest solve, and the first draft understated it by quoting lba7r alone. Activation was recomputed at every collapsed rate, 2,880 times per rebuild, and took 3.9 of the 8.0 s. That is fixed in `89d574e9`: T30's solve is now 4.5 s, and lbaE's fell 0.84 → 0.61 s.

**What remains**
- A burn-stage rebuild profile. The capture script failed (`seq.int` given a vector `from`) and was not pursued.
- Known-answer oracle tests for the collapsed joint grid and for `MIN_GAIN` as separate cases. The new honesty test covers credit, regime, particle floor and activation; the existing tests cover the rest only partly.
- Multi-seed counterfactuals.
- The unbought-coordinate baseline (section 3, last point).
- A coarse regime grid ({0.5, 1} for every subject).

No constant or branch was removed. Calibration on single seeds is not enough to show that one changes no decision.

**What WS4's test set must contain to judge it properly**
- Coordinates whose true lag-one rho under each kernel alone, and in combination, is known, so that factor composition can be tested exactly. My factor-ratio yardstick is biased by about 0.07 even with nothing bought.
- A funnel geometry where IACT is much larger than the lag-one tau, to test the lag-one surrogate where it should fail.
- Correlated group means (the shear's designed-for case), plus a case where a group coordinate binds and one where a subject binds.
- Many near-tied coordinates, to measure flip cost.
- A case where HMC beats the particle step and one where it doesn't.
- A collapsed-move case with a known factor.
- Blocks long enough, or enough replicate seeds, that a single coordinate's realised rho has standard error at most 0.03.

## 1. Commits on `ws1/allocator`

| Commit | What | Test |
|---|---|---|
| `225ce06e` | The plan without HMC kept `h` as a NULL element. `x$h <- NULL` let `x$h` partially match `h_dev`, so the "without HMC" plan bought HMC priced at the deviation factors. This affected the optimistic pilot pass and any unit without `h_n`. | The solver's objective equals the model objective at its own plan (credit, regime floor, particle floor, activation). It fails in 11 places on the old code. |
| `038017cd` | Records `pred_group`, `pred_subject`, `pred0_group` and `pred0_subject` per rebuild (`.emc_alloc_predict`). | Predictions equal the composition the solver assumes and are bounded by its target. |
| `513f5e96` | Movable subject units priced in one C++ call (`emc_alloc_subjects_need_cpp`), with the allocation copied across targets above every level. | Batched pricing equals per-unit pricing, including absent p/h, held q/r, r_min and refresh. Verified on every call of 15 captured rebuilds: 8,826 calls, q identical, x and cost within 2.2e-16. |
| `0df29a3e` | Measurement: row variances via rowSums, log sds by one indexed read, subject draws permuted once. | Row variance and indexed log sd match `stats::var` and the per-draw log sd, including draws = array rank, a bug the test caught before commit. |
| `1295e040` | Pooled draw covariances computed once per rebuild and shared by the chains (`diag(A S A')` per chain). | Shared-cache measurements equal uncached ones and equal the directly transformed draws. |
| `6b7bcfec` | Harness, spec and this report. | n/a |
| `89d574e9` | A family's activation term is computed once per solve unless the collapsed move shifts one of its units (only the HMC subjects' term was hoisted, since the collapsed move shifts the shear units). Column products are in C++ with long-double accumulation, like R's `prod`. | The compiled product is `identical()` to the R expression, including non-finite rates. Plans are `identical()` on the five captured inputs replayed with deterministic pilots. |
| `34aa3eca` | Movable shear coordinates priced in one C++ call, with their fixed fields gathered once per solve; a collapsed rate changes only their rho. Uses R's `R_pow` and no fused multiply-adds, so pricing is bit-identical to R. | Random shear units `identical()` to the per-unit path. All 120 solves of the 15 captured rebuilds identical. T30 solve 4.5 → 3.2 s. |
| `0f667b21` | The solve's per-rate bookkeeping: each unit's level is computed once and only the moved units' are rescaled per rate; shear and fixed results come back from C++ as the lists the solve uses; activation reads those lists; the plan is read out once. | All 120 captured solves identical. T30 solve 3.25 → 1.39 s, rebuild 4.0 → 2.2 s. |
| `145c30cd` | Review Op 7: the window records each subject's HMC likelihood rows. A subject with no HMC in the block is priced as rows per transition × its data rows × the chain's HMC CPU per row, instead of its evidence's CPU, which includes the rebuild's own pilots. | Pricing ignores the evidence's CPU once a time per row exists. Two T30 rebuilds with real pilot CPU: one chain's cost differed before, identical after. A live T30 run records rows, and its next rebuild is identical across two runs on three cores. |
| `a488164e` | Spec S1: shear and collapsed records carry the iterations their evidence spans, discounted with it, and the cap is 2 n_ev / iters. Records without a span keep the old reading. | Steady rate over doubling blocks gives a cap of exactly twice the rate. Without a span, the old reading. The span of a closing block is the window's iteration count. |

Verification of the two measurement commits on six captured production inputs (Forstmann, T30, lbaE, nback31, lba7r ×2):
- measurements agree to about 1e-15;
- with pilot CPU made deterministic, the five single-core inputs give plans all.equal at 1e-10;
- the 4-core input varies under the old code too (item 4 of the verdict).

All five sampler suites pass at the default and full levels after each commit.

## 2. Calibration: predicted against realised

Each fixture's sample-stage checkpoint was continued for 500 iterations, rebuilding at +0, +100 and +300. Rows are chain × coordinate × rebuild, over the two blocks of 100 and 200.

Column definitions:
- **Pred (as solved):** the plan's prediction, including the confidence bound on the rest of the iteration.
- **Pred (point):** the same with the bound removed.
- **Realised (ESJD):** the allocator's own estimator on the block, against the rebuild's pooled variance.
- **Realised (ACF):** the standard lag-one autocorrelation.
- **Bias, SD, Corr:** realised (ESJD) minus, and correlated with, the point prediction.
- **IACT / tau:** the realised IACT over the tau the point prediction implies (median).

| Fixture | Kind | n | Pred (as solved) | Pred (point) | Realised (ESJD) | Realised (ACF) | Bias (point) | SD | Corr (point) | IACT / tau |
|---|---|---|---|---|---|---|---|---|---|---|
| forstmann | group mean dirs | 36 | 0.000 | 0.000 | -0.006 | -0.001 | -0.006 | 0.151 | n/a | 1.03 |
| forstmann | group log sds | 36 | 0.002 | 0.012 | -0.010 | -0.000 | -0.023 | 0.168 | 0.03 | 1.02 |
| forstmann | subjects | 684 | 0.004 | 0.018 | 0.007 | -0.005 | -0.010 | 0.157 | -0.01 | 1.03 |
| T30 | group mean dirs | 117 | 0.214 | 0.242 | 0.256 | 0.241 | +0.014 | 0.158 | 0.70 | 1.13 |
| T30 | group log sds | 117 | 0.421 | 0.450 | 0.465 | 0.457 | +0.015 | 0.179 | 0.54 | 1.10 |
| T30 | subjects | 4680 | 0.110 | 0.197 | 0.214 | 0.206 | +0.017 | 0.203 | 0.34 | 1.19 |
| lbaE | group mean dirs | 99 | 0.367 | 0.403 | 0.407 | 0.360 | +0.005 | 0.142 | 0.91 | 1.37 |
| lbaE | group log sds | 99 | 0.652 | 0.677 | 0.659 | 0.631 | -0.017 | 0.100 | 0.90 | 1.51 |
| lbaE | subjects | 7920 | 0.183 | 0.261 | 0.303 | 0.267 | +0.042 | 0.229 | 0.40 | 1.49 |
| nback31 | group mean dirs | 189 | 0.321 | 0.354 | 0.324 | 0.296 | -0.030 | 0.118 | 0.93 | 1.12 |
| nback31 | group log sds | 189 | 0.313 | 0.363 | 0.313 | 0.301 | -0.050 | 0.110 | 0.72 | 1.23 |
| nback31 | subjects | 5859 | 0.144 | 0.243 | 0.270 | 0.231 | +0.027 | 0.201 | 0.37 | 1.33 |
| lba7r | group mean dirs | 180 | 0.350 | 0.408 | 0.351 | 0.336 | -0.057 | 0.164 | 0.85 | 1.21 |
| lba7r | group log sds | 180 | 0.517 | 0.585 | 0.481 | 0.454 | -0.104 | 0.143 | 0.60 | 1.20 |
| lba7r | subjects | 19080 | 0.235 | 0.312 | 0.301 | 0.277 | -0.011 | 0.176 | 0.50 | 1.16 |

Over- and under-crediting by kind:
- **Group mean directions and log sds:** no consistent sign. Mildly over-credited on T30 (realised slower by +0.015). Under-credited on nback31 and lba7r (realised faster), where the rest of the iteration itself sped up between rebuild and block: lba7r's rest went from 0.55 to 0.49, nback31's from 0.37 to 0.32.
- **Subjects:** over-credited by +0.02 to +0.04 on T30, lbaE and nback31. As solved (with the confidence bound), subject predictions are optimistic by 0.07–0.12 because subject windows are not pooled; that is the spec's S2.

**The binding coordinate.** Per chain and rebuild, the coordinate with the largest point prediction against the largest realised:
- On lbaE and nback31 the binding group coordinate realises near its prediction (for example nback31 mean direction 1: predicted 0.82–0.89, realised 0.83–0.89).
- On T30, lag-one noise on the binding log sd is about ±0.15.
- The realised maximum is often a subject. The persistence check shows those are mostly noise:

| Fixture | Subject realised r (block 2 vs 3) | Point pred vs realised, block 3 | Top 20 by block-2 realised: block 2 → block 3 (predicted) |
|---|---|---|---|
| T30 | 0.23 | 0.38 | 0.76 → 0.32 (0.30) |
| lbaE | 0.13 | 0.34 | 0.80 → 0.47 (0.39) |
| nback31 | 0.15 | 0.44 | 0.98 → 0.32 (0.43) |
| lba7r | 0.01 | 0.46 | 0.81 → 0.29 (0.33) |

For group coordinates, block-to-block r is 0.60–0.92 and prediction-to-realised r is 0.82–0.95.

**Lag one against IACT.** IACT/tau is 1.03–1.51 and is a property of the fixture: lbaE has the highest, about 1.5, on both log sds and subjects. On lbaE, the 20 subjects with the largest block-2 lag-one rho had a mean block-3 IACT of 18. That is the slow-mode failure the brief anticipated. It applies to particular subjects, and it is not caused by the moves bought.

**The ensemble (6b).** "The rest of the iteration" includes the ensemble update. On group coordinates with nothing bought, the realised factor relative to the measured rest is 0.88, where 1 is expected (prediction 0.945). Either part of the iteration after the rest window also moves these coordinates, or the ratio of two noisy estimates is biased. I could not separate the two. If the first, the allocator under-credits whatever runs after the window. WS5 should know about this.

## 3. Counterfactuals

Arm `benchmarks/allocator/arm_force.R`, 300 iterations, rebuild at +0 and +100:
- **up** forces 1 round onto every unbought shear coordinate, and HMC at rate 0.5 onto every unbought subject that has measured HMC factors;
- **down** divides every bought round count and HMC rate by 3.

The predictions are rescaled to the forced plan. Mean rho, point prediction against realised (ESJD), forced coordinates only:

| Arm | Group mean dirs | Group log sds | Subjects (forced HMC) |
|---|---|---|---|
| T30 up | 0.195 → 0.220 (n 9) | 0.317 → 0.336 (n 30) | 0.210 → 0.263 (n 1308) |
| T30 down | 0.296 → 0.326 (n 56) | 0.443 → 0.485 (n 61) | 0.546 → 0.486 (n 110) |
| lbaE up | n/a (none unbought) | 0.463 → 0.400 (n 18) | 0.239 → 0.307 (n 2732) |
| lbaE down | 0.497 → 0.538 (n 51) | 0.695 → 0.699 (n 59) | 0.397 → 0.355 (n 1104) |
| lba7r up | 0.367 → 0.392 (n 6) | 0.319 → 0.294 (n 82) | 0.202 → 0.228 (n 1040) |

- **Off-optimum shear predictions hold within 0.02–0.04 in rho.** The exception is lbaE log sds forced up, which mixed better than predicted (−0.06).
- **HMC forced onto subjects the allocator declined realises 0.03–0.07 slower than predicted.** The allocator abstained even at the optimistic value, so abstention was right.
- **When HMC is cut to a third, subjects mix better than predicted** (−0.04 to −0.06): they lose less than the model expects, so the HMC bought is, if anything, over-credited. One seed; the noise from section 2 applies.

A note on the yardstick: I first compared factor ratios (realised/rest against predicted/rest). They showed gaps of ±0.07 by fixture: T30 +0.07 consistently, lba7r −0.06. The gap is flat in the rest's level, which rules out the composition hypothesis (rest × move overstating the combined effect). But the ratio has a bias of the same size on coordinates with nothing bought (0.88 against 1). Read these tables in rho terms, not as factors.

## 4. Decision stability

Method: moving-block bootstrap (L = 10, 12 replicates) of the last block's trace, followed by a re-solve. Four variants are resampled:
- **rest:** the rest-of-iteration jumps;
- **moves:** shear and collapsed block records, HMC transitions and particle steps;
- **costs:** every recorded cost;
- **all:** all three together.

The bootstrap first checks that replaying the rebuild reproduces the recorded plan, and that rebuilding the inputs from the trace changes nothing. Both checks passed on every run (7/7, 3 chains each). A first set of runs reused the rebuild's RNG seed for resampling, so it held only two distinct replicates per variant. That was fixed in the harness, and the results below are from the rerun.

All-inputs variant (10–90% range across replicates; plan value in brackets):

| Run | tau SD / tau | Mean-dir rounds | Log-sd rounds | Collapsed | HMC subjects | Regime flips | Same binding coord |
|---|---|---|---|---|---|---|---|
| T30 block 3 | 0.21–0.48 / 2.8–4.3 | 0.47–0.73 (0.51–0.72) | 0.9–6.1 (1.9–7.8) | 2–6 (3–6) | 0–15 (3–10) | 0 | 8–67% |
| T30 block 2 | 0.86–1.85 / 5.5–6.8 | 0.11–0.36 (0.12–0.19) | 0.7–2.3 (1.2–1.9) | 2–4 (3–4) | 0–2 (0–1) | 0 | 0–50% |
| lbaE block 3 | 0.96–1.07 / 12.5–13.5 | 0.43–0.65 (0.44–0.62) | 2.2–3.4 (3.15–3.23) | 0.02–0.06 | 0 | 0 | 17–92% |
| lbaE block 2 | 1.55–1.83 / 12.6–14.6 | 0.18–0.60 (0.21–0.48) | 1.9–3.5 (3.1–3.2) | 0–0.04 | 8–16 (17) on one chain | 0 | 8–83% |
| nback31 block 3 | 0.44–0.83 / 12.4–14.7 | 1.0–1.6 (1.1–1.6) | 0.26 (0.26) | 0.03–0.05 | 0–4 (0–1) | 0 | 50–100% |
| forstmann block 3 | 0.06–0.12 / 1.2–1.5 | 0.16–0.20 (0.18) | 0.02–0.06 (0.03) | 0–0.06 | 0–9 (0) | 0 | 8–17% |
| lba7r block 3 | 0.13–1.15 / 6.8–7.2 | 2.0–3.8 (3.2–3.4) | 0.32–0.52 (0.41–0.48) | 0.04–0.05 | 86–89 (88–89) | 0 | 58–75% |

Attribution, comparing the variants:
- **Costs** move little: tau SD is at most 0.24, and rounds vary only on one T30 chain.
- **Rest** moves tau by 1–6% (about 10% on Forstmann, where tau is near 1) and moves rounds a little. On T30 it also moves HMC counts, because subject windows are short.
- **Moves** reproduce most of the all-inputs spread: rounds, collapsed rate and HMC counts. The moves' own evidence (a few hundred rounds and transitions per block) is what makes plans flip.

Material flips:
- T30's log-sd rounds range 0.9–6.1 across resamples.
- On lbaE chain 1, HMC switches between 8 and 16 subjects.
- lba7r is stable except on chain 3, where the moves' evidence alone spreads mean-direction rounds over 2.0–4.2 and gives tau an SD of 1.4. About 89 subjects take HMC in every replicate.
- On T30 block 3 chain 3, the plan's own value (7.8 log-sd rounds) lies outside the all-inputs resample range (2.8–6.1) and the moves-only range (3.8–6.9). The plan sat at the edge of what its own move evidence supports.

**What the flips cost (regret).** Each resampled plan was scored under the model built from the block's full evidence, relative to the plan that evidence gives (`REGRET=1`). The scoring reproduces the solver's own objective on its plan (ratio 1.000 on five chains, 0.999 on one).

| Run | Variant | Regret median / 90% / max |
|---|---|---|
| lbaE block 3 | moves | 6–7% / 12–13% / 12–13% |
| lbaE block 3 | all | 5–6% / 10–13% / 17–20% |
| lba7r block 3 | moves | 4–7% / 7–11% / 11% |
| lba7r block 3 | all | 5–6% / 9–10% / 10–15% |

Both runs used the current library, with the S1 fix.

So a plan made from one block's noisy move evidence typically costs 5–7% of the efficiency the model says the full evidence would buy. That is several times the solver's gap to the oracle (section 5).

On lba7r, a resampled plan also gives 22–50 subjects a particle and HMC rate pair that the full evidence predicts is slower than the particle step alone. The objective hides this, because group coordinates bind.

These figures are the most smoothing of the move evidence (review Op 4) could recover by the model's own reckoning, before smoothing's own lag costs. Evidence is already pooled with decay 0.5. Whether a longer memory, or hysteresis on the plan, recovers a useful share needs a fit-level comparison; it is not done here.

A flat objective near its minimum makes the binding coordinate change often, while tau typically varies 5–15% (3.5–34% across chains). The rounds bought on the non-binding coordinates are what flips, and they cost little at the margin.

## 5. Solver against its own model

`benchmarks/allocator/oracle_gap.R`, 40 random problems (one shear coordinate, two subjects with credit, regime, floors and HMC fixed cost):
- **After the `h` fix:** gap median 0.999, 90th percentile 1.019, maximum 1.117. The claimed objective equals the model objective in every trial. The largest gaps arise where the oracle gives HMC to one subject only, while a regime has to raise every subject.
- **Before the fix:** median gap 1.14; some plans were infeasible; claimed against model objective went as low as 0.74.

## 6. Cost

lba7r sample-stage rebuild input, 3 chains, captured with 4 cores (timing on 9 cores):

| | Rebuild elapsed | Main-process CPU | Solve | Measure |
|---|---|---|---|---|
| `0a5773dd` (per-unit pricing) | 6.6–7.2 s* | 3.45 s | 1.9 s | 1.05 s |
| `513f5e96` (batched pricing) | 5.4–6.0 s | 2.32 s | 0.75–0.79 s | 1.06 s |
| `1295e040` (shared measurement) | not re-timed | about 1.6 s (estimate) | 0.75 s | 0.29 s (same core, 3 chains) |

Solve time per rebuild across fixtures, one core, deterministic pilots, from captured inputs (the lba7r rows above used a four-core input):

| Fixture | Rebuild elapsed | Solve before `89d574e9` | of which activation | Solve after |
|---|---|---|---|---|
| T30 | 8.8 → 4.8 s | 8.0 s | 3.9 s (2,880 calls) | 4.5 s |
| lbaE | 2.3 → 2.0 s | 0.84 s | 0.25 s | 0.61 s |
| lba7r | 6.2 s (after) | 1.25 s | 0.35 s | 0.94 s |
| nback31 | 3.9 s (after) | 0.42 s | 0.09 s | 0.40 s |

Elapsed time is with custom kernels restored. Without restoring, lba7r's rebuild measures 10.4 s and nback31's 5.7 s, because every pilot likelihood call reloads the kernel (item 4).

After `89d574e9`, T30's solve was mostly R pricing of the shear units: those units are re-priced at each of about 40 collapsed rates, in each of 72 solves. `34aa3eca` moved that pricing to C++, taking T30 to 3.2 s and lbaE to 0.54 s, with lba7r and nback31 unchanged within noise. `0f667b21` then removed the per-rate R bookkeeping, taking T30's solve to 1.39 s (from 8.0 s at the start), lbaE's to 0.46 s and lba7r's to 0.74 s. What remains is the C++ shear pricing itself (about 0.56 s on T30) and summing costs.

Measurement per rebuild on the six captured inputs went from 0.27–1.01 s to 0.03–0.29 s. \* These elapsed times, and the roughly 3 s of elapsed time not spent in the main process, were measured on saved states whose custom kernels were not restored. They include a kernel reload on every pilot likelihood call (item 4), so they overstate a live rebuild. The main-process CPU, solve and measurement columns are unaffected. The remaining time is mostly the forked HMC pilot runs, which were not profiled separately. The pilot code lives in `subject_hmc.R`, so speeding it up is a WS2 matter. A burn-stage rebuild profile was attempted, and its capture failed. Burn rebuilds come every 100 iterations, so their split matters more than the sample stage's.

## 7. Other findings

- **Evidence caps (spec S1), fixed in `a488164e`.** Replay on the five captured inputs with the next block twice the last (the production schedule), with deterministic pilots:
  - **Shear caps:** a cap above the exploration floor binds on at most one coordinate per plan. Many more coordinates (up to 10 on nback31) sit at a cap equal to their exploration floor, which is evidence collection, not S1.
  - **Collapsed cap:** binds on T30.
  - **Effect:** the plans' modelled tau falls on 6 of 15 chains:
    - lba7r: 8.21 → 7.47 and 5.96 → 5.39;
    - nback31: 15.1 → 14.2;
    - T30: 4.35 → 4.18, 5.60 → 5.03 and 4.27 → 3.61.
  - **No chain gets worse** beyond 0.04 (lbaE chain 3).
  - **Next-block length:** under the fix, the plans are the same whether the next block is 200 or 400 iterations; before it, they changed with it.
  - **Not measured:** realised mixing.
  - **Cross-chain pooling, unchanged:** a log-sd coordinate's evidence still pools across chains, while its span is the chain's own. The log-sd cap therefore scales with the number of chains, as it did before.
- **Which coordinate sets the target (review Op 3, empirical-Bayes shrinkage of subject rho).** Across the 45 calibration plans, a subject coordinate binds the minimax in 11. Nine of those are Forstmann, where tau is about 1.2 and the binding coordinate barely matters. On the four fixtures with tau above 2, subjects bind 2 of 36 plans:
  - T30 block 3 chain 3: predicted 0.47, realised 0.41;
  - lba7r block 1 chain 2: predicted 0.74, realised 0.64.

  Group coordinates set the target everywhere else. Shrinking subject estimates would therefore rarely move tau. It could only change which subjects buy HMC to meet a target the group sets, and the counterfactuals (section 3) found the declined HMC correctly declined. Not implemented: no plan in these runs shows the noise-chasing the review describes.
- **Exploration floors persist in production.** Steady evidence is about 13.3 against a target of 20 under decay 0.5, so unbought coordinates keep about 6.7 exploration rounds per block. Under the doubling schedule that cost per iteration halves at each rebuild.
- **The particle rate (WS3 interaction).** On lba7r, q is about 0.16 with r about 0.38 for about 90 subjects (q + r < 1). The allocator relies on the particle step's measured factor p, so any WS3 change to the particle step feeds through automatically.
- **Pricing pilots by work (review Op 7), done in `145c30cd`.** In a live fit, the pilots' measured CPU was the only run-to-run variation (verdict item 4). It reaches a plan only through subjects with no HMC transitions in the block, which were priced from `ev$cpu / ev$n`.
  - Before the fix it mattered rarely: in two-run replays of four captured inputs, only one T30 chain's cost and objective differed, and the chosen rates did not.
  - Such subjects are now priced as rows per transition × the subject's data rows × the chain's HMC CPU per row. That time comes from the block's production transitions, or is carried from the last rebuild.
  - Pilot CPU still prices before any production HMC has been timed, at the first rebuild that buys HMC.
- **Regime coarseness.** The regime grid is {0.5, 1} applied to every subject. This causes the oracle's worst gaps. It also produced the one regime switch seen in calibration: at lbaE's first rebuild, one chain put all 80 subjects on HMC while the other chains put 0 and 49.

## 8. Reproducing

See `benchmarks/allocator/README.md`. Fixtures are listed there. The checkpoints for T30 and nback31 are:
- `/data/work/EMC2_dev_oo/.Rtmp/claude-1000/-data-work-EMC2-dev-oo/3baace33-74ca-4333-ba0b-800d175b2295/scratchpad/runs/final/T30_fullc_1_emc.rds`
- `/data/work/EMC2_dev_oo/.Rtmp/claude-1000/-data-work-EMC2-dev-oo/d88c8a86-efae-415e-999b-147668d64e27/scratchpad/bench/A/nback31_fullc_1_emc.rds`
