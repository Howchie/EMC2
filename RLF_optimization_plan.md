  ### 1. EMC2 Framework Integration & Deduplication (SolveCache)                                                                                   
                                                                                                                                                   
  Currently, RLF (model_RLF.h) is written as a standalone solver that runs complete PDE evaluations independently for each call. Adopting ROU's    
  cache architecture (src/fpe_race.h and src/model_ROU.h) yields immediate order-of-magnitude speedups:                                            
                                                                                                                                                   
  • Parameter Keying & Deduplication (fperace::Key):                                                                                               
  In typical design matrices, hundreds of trial rows share identical parameter combinations. By hashing parameter tuples and solving only once per 
  unique tuple (fperace::cache_get_batch), a dataset with 500 rows across 12 design cells requires 12 PDE solves instead of 500.                   
  • t₀ Decoupling:                                                                                                                                 
  Excluding non-decision time t₀ from the cache key allows trials that differ only in t₀ to share the exact same PDE numerical march (since t₀ only
  shifts the time axis tt = rt - t0).                                                                                                              
  • State Rescaling (σ-nondimensionalization):                                                                                                     
  Scaling state and boundary parameters by noise scale (Y = X/σ) standardizes the diffusion parameter to σ = 1, maximizing key collisions across   
  design cells.                                                                                                                                    
  • Sparse Query Output (query_times):                                                                                                             
  Rather than outputting full dense time grids, record and interpolate log-density (log f) and log-survivor (log S) values only at the exact       
  requested response times (rtᵢ).                                                                                                                  
  ──────                                                                                                                                           
  ### 2. Linear Algebra & Nonlocal Solver Optimizations                                                                                            
                                                                                                                                                   
  In model_RLF.h, the Crank-Nicolson implicit system                                                                                               
                                                                                                                                                   
    ⎛     Δt  ⎞       ⎛     Δt  ⎞                                                                                                                  
    ⎜I - ────L⎟pⁿ⁺¹ = ⎜I + ────L⎟pⁿ                                                                                                                
    ⎝     2   ⎠       ⎝     2   ⎠                                                                                                                  
                                                                                                                                                   
  is solved at every time step using triangular forward/back substitution (RLF_DenseLU::solve), which costs O(N²) and has strict row-by-row data   
  dependencies that prevent SIMD vectorization.                                                                                                    
                                                                                                                                                   
  • **Explicit Inverse Precomputation (                                                                                                            
                                                                                                                                                   
                      -1                                                                                                                           
           ⎛     Δt  ⎞                                                                                                                             
    M    = ⎜I - ────L⎟                                                                                                                             
     inv   ⎝     2   ⎠                                                                                                                             
                                                                                                                                                   
  )**:                                                                                                                                             
  Because the nonlocal operator L and step size Δt are stationary during a march, precompute                                                       
                                                                                                                                                   
    M                                                                                                                                              
     inv                                                                                                                                           
                                                                                                                                                   
  once at setup time (O(N³)).                                                                                                                      
  Each time step then becomes a single dense matrix-vector product (                                                                               
                                                                                                                                                   
    pⁿ⁺¹ = M   bⁿ                                                                                                                                  
            inv                                                                                                                                    
                                                                                                                                                   
  ), which:                                                                                                                                        
                                                                                                                                                   
  1. Replaces recursive forward/back substitution with independent multiply-accumulate operations.                                                 
  2. Enables full SIMD auto-vectorization (AVX2 / AVX-512 / FMA) and cache-line alignment.                                                         
                                                                                                                                                   
  • Exploiting Toeplitz & FFT Structure:                                                                                                           
  The Grünwald-Letnikov spatial discretization on a uniform grid generates a Toeplitz matrix (the transition rate q[|i - j|] depends only on grid  
  index distance |i - j|).                                                                                                                         
      • Matrix-vector multiplications L·p can be computed in O(N log N) using Fast Fourier Transforms (FFT).                                       
      • For fine spatial grids, iterative Krylov solvers (e.g., GMRES or BiCGSTAB with circulant preconditioners) reduce step complexity from O(N²)
      to O(N log N).                                                                                                                               
  • Banded / Truncated Kernel Approximations:                                                                                                      
  The heavy-tailed jump weight                                                                                                                     
                                                                                                                                                   
              -(1+α)                                                                                                                               
    q[d] sim d                                                                                                                                     
                                                                                                                                                   
  drops off with distance d. Truncating long-range interactions beyond a cutoff K ll N converts the dense matrix into a banded matrix of bandwidth 
  K, reducing step complexity to O(K·N) while preserving boundary behavior.                                                                        
  ──────                                                                                                                                           
  ### 3. Time and Space Grid Efficiency                                                                                                            
                                                                                                                                                   
  Currently, rlf_solve defaults to uniform time steps Δt and uniform spatial resolution h, running up to 5 complete PDE solves per evaluation to   
  verify domain expansion and spatial refinement (rlf_fast_path_safe checks).                                                                      
                                                                                                                                                   
  • Non-Uniform Time Step Grading (tgrade / FPE_TimeSchedule):                                                                                     
  Implement ROU's time-grading schedule. Rapid probability flux changes near t ≈ 0 require fine Δt, but later times can use exponentially larger   
  time steps. This reduces total steps                                                                                                             
                                                                                                                                                   
    N                                                                                                                                              
     step                                                                                                                                          
                                                                                                                                                   
  by 5 × to 10 × for a given horizon                                                                                                               
                                                                                                                                                   
    t                                                                                                                                              
     max                                                                                                                                           
                                                                                                                                                   
  .                                                                                                                                                
                                                                                                                                                   
  • Spatial Mesh Grading (grade):                                                                                                                  
  Grading spatial nodes near the starting point z₀ and upper boundary b₀ achieves target accuracy with significantly smaller grid sizes (e.g., N = 
  60 - 80 instead of N = 200 - 400). Since dense matrix ops scale as O(N²), halving N yields a 4 × speedup.                                        
  • Single-Pass Calibrated Solves:                                                                                                                 
  Disable multi-pass adaptive refinement during likelihood calls. Provide user-configurable R options (e.g., emc2.fpe_nx and emc2.fpe_dt) so users 
  can calibrate grid resolution once rather than re-solving multiple grids per particle evaluation.                                                
  ──────                                                                                                                                           
  ### 4. Multi-Lane SIMD Batching                                                                                                                  
                                                                                                                                                   
  ROU interleaves 4 (AVX2) or 8 (AVX-512) independent parameter solves across SIMD vector registers in fpe_solve_batch_ou_lanes.                   
  
  • For RLF, combining Explicit Inverse Precomputation with Multi-Lane Batching turns the inner time loop into a batched Matrix-Matrix             
  multiplication (GEMM of size N × N by N × LANES).
  • Modern CPUs execute batched GEMM at near peak FLOP capacity, enabling 4 to 8 RLF parameter solves to execute concurrently with minimal overhead.
  ──────
  ### Recommended Implementation Roadmap
  
  1. Phase 1 (Framework & Deduplication): Wrap rlf_solve in fperace::SolveCache with fperace::Key deduplication, t₀ decoupling, and sparse output  
  lookup.
  2. Phase 2 (Linear Algebra & Grids): Replace triangular solves with precomputed
  
    M
     inv
  
  matrix-vector products, and introduce time grading (FPE_TimeSchedule).
  3. Phase 3 (SIMD & Nonlocal Fast Methods): Implement 4/8-lane AVX batching for
  
    M   ·p
     inv
  
  , and evaluate Toeplitz FFT / kernel truncation for large N.

---

## Implementation outcome (2026-07-30)

The production implementation follows all three roadmap phases:

- `rlf::SolveCache` hashes the dimensionless tuple
  `(v/s, alpha, (B+A)/s, A/s)`, excludes `t0`, groups repeated rows, and stores
  sparse log-density/log-survivor answers at requested response times. Scalar
  censoring/truncation queries upgrade a sparse entry to a full grid only when
  required.
- The stationary TR--BDF2 system uses one LAPACK inverse per time-step block.
  Its two implicit stages are collapsed into one rational propagator, so each
  full step is one dense BLAS matrix-vector product. The original LU path
  remains available as a numerical reference.
- Piecewise-constant time grading, configurable through
  `emc2.rlf_tgrade`, is implemented with a positivity-aware schedule. Likelihood
  calls default to a calibrated single pass (`emc2.rlf_adaptive = FALSE`);
  adaptive domain and spatial refinement remains available for standalone
  validation.
- An optional structure-of-arrays SIMD path can interleave four AVX2 or eight
  AVX-512 solves. Re-measurement with the extrapolated production path found it
  slower at normal occupancy, so it remains available for diagnostics but is
  off by default.

The optional fast-method ideas were evaluated but not enabled:

- The finite-domain operator is not strictly Toeplitz after the conservative
  lower-boundary closure, and its inverse is dense. FFT multiplication would
  therefore require a different iterative solver and preconditioner, with no
  benefit at the calibrated `nx <= 200` range.
- Truncating the heavy-tailed kernel changes both the stable generator and the
  upper first-passage flux. It is not used without a separate error-controlled
  approximation.
- A nonuniform spatial mesh is incompatible with the shifted Grünwald stencil
  used here. Spatial grading is deferred until a conservative nonuniform
  fractional operator can be validated; the existing uniform-grid adaptive
  refinement remains the accuracy reference.
- Graded time schedules require an inverse for every distinct step block.
  Benchmarks showed that `tgrade = 1` is the fastest calibrated default for this
  dense-inverse backend, while retaining grading as an opt-in resolution knob.

The original `nx = 100` development benchmark measured:

| Optimization | Result |
|---|---:|
| 48 row solves reduced to 8 unique keys | 8.73x |
| SIMD batch versus serial batch (historical, before re-measurement) | 2.46x |
| Explicit inverse versus per-step triangular solve | 2.31x |
| Inverse/reference maximum PDF error | 2.91e-14 |
| Inverse/reference maximum CDF error | 1.99e-14 |

The benchmark source is `WorkingTests/benchmark_rlf.cpp`. Regression coverage is
in `tests/testthat/test-rlf-fht.R` and `tests/testthat/test-rlf-model.R`.

## Production follow-up (2026-07-31)

Profiling the actual single-subject recovery likelihood, rather than isolated
solver calls, exposed two further bottlenecks:

- The time march still used a compiler-generated row dot product. Passing each
  complete propagator multiply to BLAS `dgemv` made the kernel 6-8x faster in
  isolation over 96-192 nodes and reduced a default 123-particle likelihood
  batch from 3.106 to 2.576 seconds without changing the solver result.
- `cores_per_chain` was only spent across subjects. A one-subject fit therefore
  left all but one of those cores idle even though its proposal likelihoods are
  independent. Single-subject initialization and sampling now route that same
  core budget through the existing proposal split; the measured 123-particle
  batch scales from 2.576 seconds on one core to 0.708 seconds on four.

A joint space/time sweep found a better production grid than simply shrinking
the old one. The default is now a 160+200 Richardson pair (`nx = 160`, ratio
1.25) at `dt = 0.016`. It takes 2.116 seconds for the same batch (32% below the
pre-follow-up default, 0.590 seconds on four cores), while improving the
low-alpha profile against a 256+384 reference. At generating alpha 1.1, 1.3,
and 1.7, peak displacements changed from -0.025, -0.020, and +0.016 for the
former 128+192/0.008 grid to -0.005, -0.0004, and +0.006. A separate time-only
profile sweep over alpha 1.1-1.9 bounded the `dt = 0.016` peak displacement
against `dt = 0.004` at 0.0023.

For short time blocks, the solver now skips construction of the collapsed
TR--BDF2 propagator and applies the already-computed inverse twice instead. The
two expressions agree to machine precision; the crossover is selected from the
block length and matrix dimension. Long blocks retain the one-GEMV collapsed
path. A LAPACK triangular-solve alternative was also measured and rejected:
at the production block lengths it was 1.4-2.8x slower.

The same recovery design was exercised through the actual fitting entry points,
not only `calc_ll_manager`. With 123 particles and 2000 expanded race rows,
single-chain initialization measured 2.08 seconds on one core and 0.58 seconds
on four (3.60x). Three preburn updates measured 10.00 and 3.01 seconds (3.32x).
The public manuals intentionally contain only option behavior and usage
guidance; calibration tables and benchmark evidence remain in this maintainer
note and the `WorkingTests` benchmarks.
