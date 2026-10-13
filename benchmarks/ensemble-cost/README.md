# Ensemble cost harness (WS5)

Paths point at the WS5 private libraries (`EMC_LIB`/`R_LIBS`) and fixtures under
`/data/work/EMC2_ws5_libs` and `/data/work/EMC2_workstreams/fixtures`; edit before reuse.

- `prof_ens2.R`: per-function wall time of the ensemble move, pool draws and Gibbs step over a sample-stage run.
- `prof_pair.sh`: `prof_ens2.R` for a base and a head library, alternating, pinned with taskset.
- `sbs_ens.R`: side-by-side identity; every pool draw, move and Gibbs step also computed by the R reference from the same inputs and `.Random.seed`.
- `sweep_mr.R`: pool size x rounds sweep of the move on saved states (in-move group displacement per round, subjects refreshed, fresh-weight odds, cost).
- `cont_m.R`, `queue_m.sh`, `queue_r.sh`: continuations from a saved state with fixed pool size and rounds, exact iteration counts, fresh chain streams per seed, each fit in its own systemd scope so its cgroup CPU is its own.
- `pool_cost.R`: single-core time of the pool likelihood for k fresh members per subject, and the implied pool CPU per iteration at M = 8, 32 and n.
- `diag_lba7r.R`: exact-length continuation printing the per-block record (used to separate per-iteration cost from extra stop-rule tries).
- `mk_t30.sh`: builds the T30 checkpoint the sweeps and continuations start from.
