# Ensemble cost harness (WS5)

Paths point at the WS5 private libraries (`EMC_LIB`/`R_LIBS`) and fixtures under
`/data/work/EMC2_ws5_libs` and `/data/work/EMC2_workstreams/fixtures`; edit before reuse.

- `prof_ens2.R`: per-function wall time of the ensemble move, pool draws and Gibbs step over a sample-stage run.
- `prof_pair.sh`: `prof_ens2.R` for a base and a head library, alternating, pinned with taskset.
- `sbs_ens.R`: side-by-side identity; every pool draw, move and Gibbs step also computed by the R reference from the same inputs and `.Random.seed`.
- `sweep_mr.R`: pool size x rounds sweep of the move on saved states (in-move displacement per round, refreshed share, fresh-weight odds, cost).
- `cont_mr.R`, `cont_queue.sh`: continuations from a saved state with fixed (M, R) or the between-block rule.
- `plan_smoke.R`: short fit from scratch printing the rule's plan history.
- `mk_t30.sh`: builds the T30 checkpoint the sweeps and continuations start from.
