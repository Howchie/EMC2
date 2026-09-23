# RNG stream overlap across blocks and chains under L'Ecuyer-CMRG

Date: 2026-09-22

Status: fixed in the build installed 2026-09-22 08:41; verified (see the end of this file). Found while fitting the N-back PM
trend7 RDM/LBA models (`/data/work/PM/NirvanaHons_Nback`). Those fit scripts
call `RNGkind("L'Ecuyer-CMRG")` before `fit()`.

## Summary

When the calling session uses `RNGkind("L'Ecuyer-CMRG")`, the random numbers
driving `fit()` are heavily reused:

1. **Across blocks within a chain.** Each block (`step_size` iterations) restarts
   a chain from an RNG state only a few hundred draws past the previous block's
   start state. In preburn it restarts from the *identical* state. The block
   therefore replays almost all of the previous block's random numbers, shifted
   by a few hundred positions.
2. **Across chains.** Chain *c*'s per-subject particle streams are the same
   streams as chain *c−1*'s, shifted by one subject: chain 2 / subject *s* is
   bit-identical to chain 1 / subject *s+1*. Chain 1's Gibbs stream is chain 2's
   last subject stream.

With the default Mersenne-Twister kind, no overlap was found. Forked children
there reseed from time and PID, so streams are independent but not
reproducible. The bug therefore hits exactly the users who choose L'Ecuyer for
reproducibility.

## Mechanism

The four pieces below combine.

1. **Chains are re-forked every block.** `R/fitting.R:164` runs
   `auto_mclapply(sub_emc, run_stages, ...)` inside the block loop.
   `auto_mclapply()` calls `parallel::mclapply()` with its default
   `mc.set.seed = TRUE`.

2. **Base R seeds children from the parent's current state.** With
   `mc.set.seed = TRUE` and L'Ecuyer, `mclapply()` calls `mc.reset.stream()`
   on *every* call. That snapshots the parent's current `.Random.seed` into
   `LEcuyer.seed`. Child *k* receives `nextRNGStream^(k-1)(LEcuyer.seed)` via
   `mc.set.stream()` / `mc.advance.stream()`. All sampling draws happen in the
   children, so the parent's `.Random.seed` hardly moves between blocks.
   Measured parent advance: 0 draws in preburn, 294–517 draws per block in
   burn/adapt/sample on the test fixture.

3. **Jump-ahead commutes with advancing.** For MRG32k3a, `nextRNGStream` (a
   2^127 jump) and ordinary advancing are powers of the same transition
   matrix, so they commute. If the parent advanced *m* draws between blocks,
   chain *k*'s new entry state is `A^m · J^(k-1) P = A^m · (old entry state)`.
   That is its previous entry state advanced by only *m* draws.

4. **Per-block and per-chain stream derivation reuses that state.**
   - `run_stage()` re-derives the particle streams at the start of every block
     (`R/sampling.R:582`, `wpool_streams <- .emc_subject_streams(...)`). The
     streams are carried forward within a block (`R/sampling.R:715`) but not
     across blocks.
   - `.emc_subject_streams()` (`R/chain_pool.R:1579`) gives subject *s* the state
     `nextRNGStream^s(current seed)`. It then sets the process seed to
     `nextRNGStream^(n+1)`, which `gibbs_step()` (`R/sampling.R:639`) consumes.
   - With chain entry seeds `J^0 P, J^1 P, J^2 P`, subject streams are
     `J^(s+c-1) P`. Chains therefore share streams offset by one subject.
   - `init()` (`R/sampling.R:332`) derives start-point streams the same way, so
     start proposals overlap across chains too. It is called from `run_stages()`
     at `R/fitting.R:282` and from `init_chains()` at `R/sampling.R:396`.
   - The fallback particle path (`R/sampling.R:730`, `parallel::mcmapply` with
     default `mc.set.seed`) seeds subject children from the chain seed the same
     way, so it should show the same cross-chain pattern. Not separately tested.

## Evidence

### Reproduction

Tiny fixture: `forstmann`, 4 subjects, `LNR`, `m ~ lR`, 3 chains,
`cores_for_chains = 3`, `cores_per_chain = 1`, `particle_factor = 5`, run with
`fit(emc, iter = 60, step_size = 20)`.

The wrapper below logs each chain's `.Random.seed` at entry to `run_stages`:

```r
library(EMC2)
RNGkind("L'Ecuyer-CMRG"); set.seed(123)
log_file <- tempfile()
orig <- EMC2:::run_stages
assignInNamespace("run_stages", function(sampler, stage, ...) {
  s <- get(".Random.seed", envir = globalenv())
  cat(sprintf("%s|%d|%s\n", stage, sampler$samples$idx,
              paste(s[2:7], collapse = ",")), file = log_file, append = TRUE)
  orig(sampler, stage = stage, ...)
}, "EMC2")
dat <- forstmann[forstmann$subjects %in% levels(forstmann$subjects)[1:4], ]
dat$subjects <- droplevels(dat$subjects)
des <- design(data = dat, model = LNR, formula = list(m ~ lR, s ~ 1, t0 ~ 1))
emc <- make_emc(dat, des, n_chains = 3)
emc <- fit(emc, iter = 60, step_size = 20, cores_for_chains = 3,
           cores_per_chain = 1, particle_factor = 5)
```

### Results

- **Preburn:** 6 block entries but only 3 distinct seeds. Each chain entered
  every block with the identical state; the first uniform was the same every
  block (e.g. `0.3411063952`).
- **Burn/adapt/sample:** entry seeds are distinct, but each is the same chain's
  previous entry state advanced by *m* draws. Found by matching the first 5
  uniforms of block *b+1* inside the first 10^6 uniforms of block *b*.
  Offsets were 294, 346, 433 (burn); 398, 471, 517 (adapt); 450, 511 (sample).
  Every chain was affected in every block.
- **Across chains:** within every sample block,
  `entry_2 == nextRNGStream(entry_1)` and
  `entry_3 == nextRNGStream^2(entry_1)`. Chain 2 / subject 1 is identical to
  chain 1 / subject 2 (`identical(...)` is `TRUE`).
- **Whole-fit count:** every stream `.emc_subject_streams()` handed out was
  logged, over all blocks, chains and init calls. Ordered pairs where stream
  *j* starts within the first 2×10^5 draws of stream *i*:

  | RNG kind | streams | overlapping pairs |
  |---|---:|---:|
  | `L'Ecuyer-CMRG` | 168 | 2433 / 28056 |
  | Mersenne-Twister (default) | 204 | 0 / 41412 |

### Symptom in a real model

The model is the N-back trend7 RDM (40 subjects, 8137 trials, 3 chains,
1000 sample iterations, `cores_per_chain = 9`). Replicate runs of the
identical model and prior disagree far beyond what the within-run diagnostics
imply:

- lpd agrees within about 2 across runs, but p_WAIC differs by 28–50 (shape-10
  variance prior) and by 50 (shape 2). So WAIC differs by 59–97 between
  identical runs.
- Within a run the three chains agree closely. Per-chain p_WAIC was 485 / 473 /
  481 in one run and 534 / 520 / 535 in the replicate, with almost no
  between-chain excess.
- Every group-level SD was 1–9% larger in the replicate, in the same direction
  for all 16 parameters. The sigma² ESS was 400–1300, so a shift of this
  size is about 8 Monte Carlo SEs.

This is the expected signature. Chains driven by overlapping, recycled
randomness look converged relative to each other, so R-hat and ESS are
overstated, while independent runs land in different places. Because the
input randomness repeats with a short offset, the transition kernel is no
longer driven by fresh draws. Some bias, beyond underestimated MC error,
cannot be ruled out.

## Workaround for users

Do not set `RNGkind("L'Ecuyer-CMRG")` before `fit()`. With the default kind,
children reseed independently, at the cost of reproducibility.

## Proposed fix

The goal is reproducible *and* non-overlapping streams that persist across
blocks.

1. **Own the chain streams.** At chain creation (`make_emc` / first `init`),
   assign chain *c* an independent L'Ecuyer stream, e.g.
   `nextRNGStream^c(base)`, and store it in the sampler (e.g.
   `sampler$rng$chain`). Derive `base` once from the user's seed.
2. **Use substreams inside a chain.** Derive subject streams with
   `parallel::nextRNGSubStream()` from the chain's stream, not
   `nextRNGStream()`. Substreams are 2^76 apart *within* the chain's stream,
   so they can never collide with another chain's stream (2^127 apart).
   Give the Gibbs/group step its own substream too.
3. **Persist state across blocks.** At the end of `run_stage()`, save the
   Gibbs stream state and the per-subject states (`wpool_streams`, and the
   `mcmapply` path's equivalent) into the sampler. At the start of the next
   `run_stage()`, restore them instead of calling `.emc_subject_streams()`
   again. Apply the same to `init()`'s start streams.
4. **Stop base R reseeding.** Call the chain-level `auto_mclapply` (and
   `init_chains()`'s `mclapply`) with `mc.set.seed = FALSE`, so each chain sets
   its own state from `sampler$rng`.
5. **Make it independent of the caller's RNG kind.** Keep the internal streams
   L'Ecuyer regardless of the user's `RNGkind()`, and restore the caller's
   kind and seed afterwards, as `.emc_subject_streams()` already does. That
   way Mersenne-Twister users get reproducibility too.

## Tests to add

- Stream uniqueness: log every stream a short multi-block, multi-chain fit
  uses (as in the reproduction) and assert zero overlapping pairs within a
  large window, under both RNG kinds.
- Block continuity: the RNG state a chain enters block *b+1* with equals the
  state it left block *b* with.
- Reproducibility: two `fit()` calls with the same `set.seed()` give identical
  samples; different seeds give different samples.
- Chain independence: no chain's subject or Gibbs stream equals another
  chain's.

## Verification of the fix (build installed 2026-09-22 08:41)

The fix adds per-chain `rng` state: chains are `nextRNGStream` apart, the
Gibbs step and subjects use `nextRNGSubStream`, `mc.set.seed = FALSE` is set,
and `init()` uses the chain's subject streams. It was checked on the same
fixture (forstmann, 4 subjects, LNR, 3 chains, `iter = 40`,
`step_size = 20`). `run_stages` was wrapped to log `sampler$rng` on entry and
exit of every call:

| check | L'Ecuyer-CMRG | Mersenne-Twister |
|---|---:|---:|
| exit state of block *b* == entry state of block *b+1*, per chain | 57 / 57 | 42 / 42 |
| Gibbs state advances during each block | 100% | 100% |
| distinct stream starts logged | 300 | 225 |
| cross-chain overlaps (within 2×10^5 draws) | 0 | 0 |
| overlaps between different streams (any role) | 0 | – |

All 2850 within-chain matches under L'Ecuyer are the same stream (same subject,
or Gibbs with Gibbs) continuing into the next block, which is the intended
behaviour. Not re-checked: the fallback `mcmapply` particle path.
