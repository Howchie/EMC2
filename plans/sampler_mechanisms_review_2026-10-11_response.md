# Response to the sampler mechanisms review of 2026-10-11

Reviewed HEAD was `55b29f4d`; this response covers the commits after it.

## Adopted

1. **Collapsed-only reachability.** A coordinate the joint move lowers
   counts as movable (`.emc_alloc_movable`, f57c955d), and so does one a
   subject-HMC regime lowers (the regimes are built once in
   `.emc_alloc_solve_units` and their lowered levels extend the movable set).
   The movable set is shared by the solves with and without HMC.
2. **Measured zero movement.** The shear and collapsed factors are no longer
   clamped at `1 - 1e-6`: a measured zero jump is a factor of exactly one,
   which `.emc_alloc_factor_inverse` already treats as out of reach, so the
   coordinate is immovable and buys nothing. Test: "measured zero movement
   is a factor of exactly one, and buys nothing" (through
   `.emc_alloc_measure`).
3. **Regimes.** A regime no longer pre-multiplies the subject's rho and leaves
   `h` for a second purchase. It sets `r_min` on each covered subject; the
   unit's need solves one factor `1 - (r_min + x)(1 - h)` at the total rate
   and prices only `x` above the regime's share (`reg$cost`). The group
   credit stays at the regime rate (the supplement only helps). The regime's
   deviation average runs over the same subjects `r_bar` averages, with
   uncovered subjects' deviations unchanged, so a zero-rate regime is the
   identity. The pilot gate uses `regime + x`.
4. **Curve switch.** `emc2.shear_curve = FALSE` now also skips atlas top-ups
   and curve construction for chains that already hold an atlas. The
   continuation arm `static` additionally strips the checkpoint's atlases;
   the first run of that arm (atlas kept, draw-for-draw identical to `def`)
   is saved as `static_atlaskept_1`.
5. **Shared particle rates.** After the cross-chain maximum is taken, a chain
   whose own rates were lower is solved again with those rates held
   (`q_fixed` on its subject units), so its HMC rates, predicted tau, cost
   and objective describe the plan it runs. The record carries
   `particle_fixed`.

## Noted, not changed yet

- **Interleaved HMC accounting (P2).** The interleave is a testing switch
  (default 1) that is a no-op unless a subject HMC rate is bought; it is a
  candidate for removal in the clean-up rather than a fix. If it stays,
  every round's tasks must go through `.emc_hmc_record_tasks` and its time
  out of `shear_el`.
- **Refresh floor is an expected count.** Agreed: under independent Bernoulli
  selection `q + r >= 1` guarantees one refresh in expectation only. The
  arm is being benchmarked as is (lbaT seed 1: slowest-direction ESS per
  kCPU-s 3.4x `def`, R-hat 1.10 vs 1.60); coupling the selection is a
  follow-up if it is adopted.
- **Curve trigger versus allocator variance reference.** Agreed that they
  answer different questions; the allocator's pooled production variance is
  also what makes a drifting chain read every coordinate as immovable
  (lbaT chain 2). Left as documented; a per-gate diagnostic is worth adding.
- **Lag-one versus integrated time.** Known limitation of the AR(1)
  surrogate (it understates the integrated time about 4x on N-back). The
  continuation driver's fixed-direction ESS is the check.

## Checks

`tests/testthat/test-move-allocation.R` (803), `test-collapsed-move.R` (90)
pass on the fixed build; the reviewer's probe script, run with its defect
assertions disabled, shows collapsed-only reachability buying the rounds,
regime and execution agreeing, a zero-rate regime as identity and no curve
built under the switch. Its zero-movement probe supplies the old clamp value
directly and no longer exercises the measurement path.
