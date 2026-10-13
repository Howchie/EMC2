# Pool size and rounds of the ensemble update (R/ensemble_move.R), set between
# blocks from what the move measured over the last block.
#
# Objective: the group state's expected squared jump per unit of wall time.
# Inside the move, the group state after round r sits J_r from the state the
# move received (squared displacement, standardised by the block's posterior
# variance, averaged over the group means and log sds). J_1 is one more draw
# of the group conditional given the subject states the move received, which
# decorrelates nothing the main Gibbs step had not; the move's gain is
# G(R) = J_R - J_1, the displacement the reselections bring, and its increments
# dJ_r = J_r - J_{r-1} (r >= 2) fall as the rounds mix on the pool-augmented
# target. The chain's own lag-1 displacement of the same coordinates over the
# block, J_tot, against the iteration's wall time C sets the exchange rate:
# a unit of cost must buy at least J_tot / C of displacement.
#   rounds: keep round r while dJ_r / c_round >= J_tot / C, with c_round the
#     wall time of one round; grow while the last round still qualifies.
#   pool: the reselections' displacement depends on M through the share of
#     subjects they refresh, m(x) = mean_s x e_s / (1 + x e_s) for x = M - 1
#     fresh members, e_s being the fresh members' mean weight relative to the
#     carried state's (the odds of displacing it, as for an independence
#     sampler), so G(x) = D m(x) with D = G / m(x_cur). Keep members while
#     D m'(x) / c_member >= J_tot / C, with c_member the pool rows' share of
#     the particle round's wall time per fresh member; grow while the current
#     pool still qualifies.
# Both quantities change only here, between blocks, so the kernel is fixed
# within one. At least two rounds always run, so that one reselection's
# displacement is observed; at least one fresh member is always drawn.

.EMC_ENSEMBLE_ROUNDS_MIN <- 2L
.EMC_ENSEMBLE_ROUNDS_MAX <- 64L
.EMC_ENSEMBLE_POOL_MAX <- 256L

.emc_ensemble_block_new <- function(p, n, M, R) {
  list(n = 0L, jump = matrix(0, 2L * p, R), eff = numeric(n), moved = numeric(n),
       rounds_wall = 0, move_wall = 0, pool_wall = 0, iter_wall = 0, M = M, R = R)
}

# After one move: the per-round displacement from the received group state,
# the fresh-weight odds, and the move's wall time.
.emc_ensemble_block_record <- function(sampler, pars_in, res, M, R, rounds_wall, move_wall) {
  p <- length(pars_in$tmu)
  n <- length(res$eff)
  blk <- sampler$ensemble_stats$block
  if (is.null(blk) || blk$M != M || blk$R != R || nrow(blk$jump) != 2L * p || length(blk$eff) != n)
    blk <- .emc_ensemble_block_new(p, n, M, R)
  if (!is.null(res$traj_mu) && ncol(res$traj_mu) == R) {
    tmu_in <- as.numeric(pars_in$tmu)
    lsd_in <- 0.5 * log(diag(as.matrix(pars_in$tvar)))
    blk$jump <- blk$jump + rbind((res$traj_mu - tmu_in)^2, (res$traj_lsd - lsd_in)^2)
  }
  blk$n <- blk$n + 1L
  blk$eff <- blk$eff + res$eff
  blk$moved <- blk$moved + res$moved
  blk$rounds_wall <- blk$rounds_wall + rounds_wall
  blk$move_wall <- blk$move_wall + move_wall
  sampler$ensemble_stats$block <- blk
  sampler
}

# After the particle round: the pool rows' share of its wall time, taking a
# pool row to cost what a particle row costs (the same likelihood call).
.emc_ensemble_record_rows <- function(sampler, pools, part_wall, n_particles, ran) {
  blk <- sampler$ensemble_stats$block
  if (is.null(pools) || is.null(blk)) return(sampler)
  M1 <- nrow(pools$rows[[1L]])
  part <- ifelse(ran, n_particles, 0)
  part[!is.finite(part)] <- 0
  share <- M1 / (sum(part) + length(part) * M1)
  blk$pool_wall <- blk$pool_wall + part_wall * share * length(part)
  sampler$ensemble_stats$block <- blk
  sampler
}

.emc_ensemble_record_iter <- function(sampler, iter_wall) {
  blk <- sampler$ensemble_stats$block
  if (is.null(blk)) return(sampler)
  blk$iter_wall <- blk$iter_wall + iter_wall
  sampler$ensemble_stats$block <- blk
  sampler
}

# The chain's draws of the group means and log sds over the last `n`
# iterations, as a matrix with one column per iteration.
.emc_ensemble_block_draws <- function(sampler, n) {
  sm <- sampler$samples
  idx <- seq.int(max(1L, sm$idx - n + 1L), sm$idx)
  mu <- sm$theta_mu[, idx, drop = FALSE]
  lsd <- vapply(idx, function(t) 0.5 * log(diag(as.matrix(sm$theta_var[, , t]))),
                numeric(nrow(mu)))
  rbind(mu, matrix(lsd, nrow = nrow(mu)))
}

# The plan for the next block from the last block's record. Returns the
# sampler with `ensemble_plan` set (pool, rounds) and the record appended to
# `ensemble_plan_history`; without a usable record the plan is unchanged.
.emc_ensemble_plan_update <- function(sampler) {
  blk <- sampler$ensemble_stats$block
  if (is.null(blk) || blk$n < 2L || blk$iter_wall <= 0) return(sampler)
  M <- blk$M; R <- blk$R; x <- M - 1L
  draws <- .emc_ensemble_block_draws(sampler, blk$n)
  v <- apply(draws, 1, stats::var)
  use <- is.finite(v) & v > 0
  if (!any(use)) return(sampler)
  J <- colMeans(blk$jump[use, , drop = FALSE] / (blk$n * v[use]))
  J_tot <- mean(colMeans(diff(t(draws[use, , drop = FALSE]))^2) / v[use])
  C <- blk$iter_wall / blk$n
  rate <- J_tot / C
  G <- J[R] - J[1L]
  # Rounds: the increments the reselections bring, against a round's cost.
  c_round <- blk$rounds_wall / (blk$n * R)
  dJ <- if (R >= 2L) diff(J) else numeric(0)
  ok <- is.finite(dJ) & dJ / c_round >= rate
  R_new <- if (length(ok) && all(ok)) min(2L * R, .EMC_ENSEMBLE_ROUNDS_MAX) else
    if (any(ok)) 1L + max(which(ok)) else .EMC_ENSEMBLE_ROUNDS_MIN
  R_new <- max(.EMC_ENSEMBLE_ROUNDS_MIN, R_new)
  # Pool: the refreshed share as a function of the fresh members, against a
  # member's cost.
  e <- blk$eff / blk$n
  e[!is.finite(e) | e < 0] <- 0
  m <- function(x) mean(x * e / (1 + x * e))
  mp <- function(x) mean(e / (1 + x * e)^2)
  c_member <- blk$pool_wall / (blk$n * x)
  D <- if (is.finite(G) && G > 0 && m(x) > 0) G / m(x) else 0
  x_new <- if (D <= 0 || !is.finite(c_member) || c_member <= 0) 1L else
    if (D * mp(x) / c_member >= rate) min(2L * x, .EMC_ENSEMBLE_POOL_MAX - 1L) else {
      grid <- seq_len(x)
      keep <- vapply(grid, function(xx) D * mp(xx) / c_member >= rate, logical(1))
      if (any(keep)) max(grid[keep]) else 1L
    }
  plan <- list(pool = as.integer(x_new + 1L), rounds = as.integer(R_new))
  record <- list(n = blk$n, M = M, R = R, J = J, G = G, J_tot = J_tot, C = C,
                 c_round = c_round, c_member = c_member, eff = mean(e),
                 moved = mean(blk$moved / blk$n), pool = plan$pool, rounds = plan$rounds)
  sampler$ensemble_plan <- plan
  sampler$ensemble_plan_history <- c(sampler$ensemble_plan_history, list(record))
  sampler$ensemble_stats$block <- NULL
  sampler
}
