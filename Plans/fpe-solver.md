## Post-implementation review (2026-09-02)

The figures below are the original planning estimates, not the landing
benchmarks.  The implementation review changed the conclusions in three
important ways:

- At the shipped ROU resolution (`nx = 512`), the established march remains the
  default.  On the three-subject Forstmann benchmark it took 0.00643 s with
  0.624 summed absolute log-likelihood error, versus 0.00650 s and 0.806 for the
  modal route.  Modal remains available with `options(emc2.rou_modal = TRUE)`;
  it wins on both axes in that benchmark at `nx = 1024`.
- A fixed-boundary singleton must stay scalar: forcing it through the padded
  lane march is a 4–6% regression.  A moving-boundary singleton still uses the
  lane march and is about 1.17× faster.
- Skipping unused RLF diagnostics saves about 8–9% at the shipped `nx = 70` and
  3–5% at `nx = 320`, with bit-identical likelihood outputs.  The larger
  estimates below did not reproduce.  The reviewed implementation also skips
  diagnostics-only lower-censor and conservation work during operator assembly.

The modal propagator is exact in time only in the narrow sense that it has no
time-step discretisation.  It still has Krylov truncation error, and total error
against a spatially converged reference can increase when removing the march's
time error removes a favourable spatial/temporal cancellation.

I did not find a credible package-wide 40–100× improvement still waiting in the current architecture. I did find one substantial OU redesign and two
  smaller, low-risk savings.

   Enhancement                        Applicable scope                                   Measured improvement    Confidence
  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
   Shift-invert spectral OU solver    Fixed-boundary ROU                  ~4–8× for singleton keys; ~1.2–1.7×    Promising prototype; needs adaptive
                                                                                  for full 8-key SIMD batches    mode selection
  ─────────────────────────────────  ──────────────────────────────────  ─────────────────────────────────────  ──────────────────────────────────────
   Skip unused RLF diagnostics        Production graded/Richardson RLF                8–22% across α=1.15–2.0    High; tested outputs were bit-
                                                                                                                 identical
  ─────────────────────────────────  ──────────────────────────────────  ─────────────────────────────────────  ──────────────────────────────────────
   SIMD for singleton OU chunks       One-key solves and SIMD tails        1.41× single fixed ROU; 10–19% for    High; differences only at FMA round-
                                                                                       tails/collapsing cases    off

  ## The significant OU opportunity

  For a fixed boundary, the OU generator built in src/fpe_solver.h:362 is autonomous and tridiagonal. Its adjacent off-diagonal entries are positive,
  so it is diagonally similar to a symmetric tridiagonal matrix:

  [
  S=W^{-1}LW,\qquad
  \frac{W_{i+1}}{W_i}
  =\sqrt{\frac{L_{i+1,i}}{L_{i,i+1}}}.
  ]

  This permits a shifted symmetric Lanczos approximation to

  [
  q(t)=\exp((t-t_0)L)q_0
  ]

  using repeated (O(n_x)) Thomas solves with (I-\gamma S). Survival and boundary flux are projected once, after which every requested time is a short
  modal sum. It uses the existing finite-volume generator and initial seed, not a Wald formula or a different likelihood.

  It removes the (O(n_x n_t)) Crank–Nicolson march currently performed at src/fpe_solver.h:675.

  Accuracy against a 16× finer CN reference (dt=0.00025, nx=512) for a representative case was:

  - 24 modes: maximum density error (8.9\times10^{-6}), CDF error (4.5\times10^{-8}), log-survival error (4.9\times10^{-7}).
  - 28 modes: density error (3.1\times10^{-7}).
  - Removing full Lanczos reorthogonalisation approximately halved runtime without changing these results through 56 tested modes.

  However, 24 modes were not universally sufficient. Low-leak and early-time cases required 32–40 modes. Therefore a fixed dimension would be unsafe.
  A production implementation should:

  1. Extend one Lanczos basis in rungs.
  2. Compare log-density and log-survival only at the actual requested times.
  3. Require positive densities and stable log values.
  4. Fall back to the current CN route if convergence or diagonal scaling is unsafe.

  The SIMD prototype measured approximately 1.7× for eight keys at 32 modes and 1.2–1.3× at 40 modes. Singleton gains were much larger. Thus this
  could materially reduce ROU cost, but it is not honestly a universal order-of-magnitude package speedup.

  It applies only where static_op() is true. Collapsing-boundary ROU and ROUp remain time-dependent. The same mathematics should cover fixed-boundary
  BOU with two flux functionals, but I did not benchmark that extension and have not included it in the estimates.

  ## Immediate RLF saving

  The production Richardson route at src/model_RLF.h:2414 computes several diagnostics that its caller never reads:

  - Eigenbasis condition estimation at src/model_RLF.h:1338.
  - Arnoldi defect construction.
  - Modal state materialisation for positivity probes at src/model_RLF.h:1457.
  - Residual, lower-pressure, flux/mass and density probes at src/model_RLF.h:1691.

  Those quantities are used by the separate adaptive validation path through src/model_RLF.h:1900, but not by rlf_extrapolate() or the production
  graded pair.

  Passing a diagnostics=false flag only from rlf_cache_solve() produced these two-grid improvements:

  - α=1.15: 21.6%
  - α=1.30: 21.7%
  - α=1.70: 15.9%
  - α=1.95: 8.6%
  - α=2.00: 7.8%

  Density, survivor, Krylov dimension, and full output checksums were identical at printed double precision. This is the strongest ready-to-land RLF
  change I found.

  ## Smaller OU correction

  The batch dispatcher currently enters SIMD only when chunk_size > 1 at src/fpe_race.h:2030. On AVX builds, the lane solver already pads unused
  lanes, so singleton chunks should use it too.

  Measured improvements were:

  - One fixed ROU key: about 1.41×.
  - One collapsing key: about 1.19×.
  - 9- and 17-key fixed batches: about 10–17%, by accelerating the final singleton chunk.

  Maximum output changes were (2\times10^{-14})–(6\times10^{-14}), caused by SIMD/FMA evaluation order.

  ## Things I tested and rejected

  - RLF edge-refined meshes worsened representative density errors from about 0.0042 to 0.0046–0.0097 and did not replace Richardson extrapolation.
  - Full dense RLF eigendecomposition took roughly 2.4–2.6 ms versus 0.54–0.80 ms for shift-invert Arnoldi.
  - Coarser RLF Krylov ladders averaged around 1.25× faster but caused log-density errors up to approximately 0.0076 and problematic α=2 tail
    flooring.

  - Direct partial eigenmodes of the OU tridiagonal operator were slower and needed around 64 modes to resolve early density. The initial-state-driven
    shifted Lanczos construction is essential.

  - I found no overlooked repeated solve in the ordinary finite-trial path. Parameter-key grouping, endpoint inclusion, horizon aggregation, and
    deduplication in src/model_ROU.h:161 are already doing the right thing. The entry points in src/fpe_diffusion.cpp:1 are validation-only and are
    not holding back likelihood evaluation.
