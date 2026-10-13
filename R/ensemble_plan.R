# Pool size and rounds of the ensemble update (R/ensemble_move.R), set between
# blocks from what the move measured over the last block.
#
# Objective: the slowest coordinate's expected squared jump per unit of wall
# time, over the group means and log sds. Inside the move, coordinate k of
# the group state after round r sits J_kr from the state the move received
# (squared displacement standardised by the block's posterior variance). J_k1
# is one more draw of the group conditional given the subject states the move
# received, which decorrelates nothing the main Gibbs step had not; the move's
# gain on k is G_k(R) = J_kR - J_k1, the displacement the reselections bring.
# The chain's own lag-1 displacement of k over the block, T_k, includes the
# gain of the plan that ran; another plan changes T_k by the change in its
# gain. The gain depends on the pool through the share of subjects the
# reselections refresh, m(x) = mean_s x e_s / (1 + x e_s) for x = M - 1 fresh
# members, e_s being the fresh members' mean weight relative to the carried
# state's (the odds of displacing it, as for an independence sampler):
#   G_k(x', R') = G_k(R') m(x') / m(x).
# The iteration's wall time C changes by c_member per fresh member (the pool
# rows' measured worker CPU spread over the chain's cores) and c_round per
# round (serial in the chain process). The next plan maximises
#   min_k (T_k + G_k(x', R') - G_k(x, R)) / (C + c_member (x' - x) + c_round (R' - R))
# over x' up to twice the current pool and R' up to the current rounds, the
# cheapest plan winning ties; the rounds double when the last round still
# raises the objective.
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

# After the particle round: the pool rows' worker CPU spread over the
# chain's cores, the wall time they add to the round.
.emc_ensemble_record_rows <- function(sampler, pools) {
  blk <- sampler$ensemble_stats$block
  if (is.null(pools) || is.null(blk) || !is.numeric(pools$cpu) || !is.finite(pools$cpu))
    return(sampler)
  cores <- if (is.numeric(pools$cores) && is.finite(pools$cores)) max(1, pools$cores) else 1
  blk$pool_wall <- blk$pool_wall + pools$cpu / cores
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
  out <- rbind(mu, matrix(lsd, nrow = nrow(mu)))
  if (!is.null(rownames(mu))) rownames(out) <- c(rownames(mu), paste0("lsd_", rownames(mu)))
  out
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
  J <- blk$jump[use, , drop = FALSE] / (blk$n * v[use])
  T_k <- colMeans(diff(t(draws[use, , drop = FALSE]))^2) / v[use]
  C <- blk$iter_wall / blk$n
  c_round <- blk$rounds_wall / (blk$n * R)
  c_member <- blk$pool_wall / (blk$n * x)
  # Gain of the reselections by round, negative displacement counting as none.
  gain <- J - J[, 1L]
  gain[!is.finite(gain) | gain < 0] <- 0
  e <- blk$eff / blk$n
  e[!is.finite(e) | e < 0] <- 0
  m <- function(xx) mean(xx * e / (1 + xx * e))
  share <- function(xx) if (m(x) > 0) m(xx) / m(x) else 1
  objective <- function(xx, rr) {
    cost <- C + c_member * (xx - x) + c_round * (rr - R)
    if (!(cost > 0)) return(-Inf)
    min(pmax(0, T_k + gain[, rr] * share(xx) - gain[, R])) / cost
  }
  # Without a measured price for the pool rows the pool stays as it is.
  xs <- if (blk$pool_wall > 0) seq_len(min(2L * x, .EMC_ENSEMBLE_POOL_MAX - 1L)) else x
  rs <- seq.int(min(R, .EMC_ENSEMBLE_ROUNDS_MIN), R)
  f <- outer(xs, rs, Vectorize(objective))
  f[!is.finite(f)] <- -Inf
  # Cheapest first: rows by pool, columns by rounds, so which() in column
  # order takes the fewest rounds, then the smallest pool, among the maxima.
  best <- which(f == max(f), arr.ind = TRUE)
  best <- best[order(best[, 2L], best[, 1L]), , drop = FALSE][1L, ]
  x_new <- xs[best[1L]]
  R_new <- rs[best[2L]]
  if (R_new == R && R >= 2L && objective(x_new, R) > objective(x_new, R - 1L))
    R_new <- min(2L * R, .EMC_ENSEMBLE_ROUNDS_MAX)
  R_new <- max(.EMC_ENSEMBLE_ROUNDS_MIN, R_new)
  plan <- list(pool = as.integer(x_new + 1L), rounds = as.integer(R_new))
  worst <- which.min(T_k)
  record <- list(n = blk$n, M = M, R = R, J = colMeans(J), G = mean(gain[, R]),
                 J_tot = mean(T_k), T_min = T_k[worst], G_min = gain[worst, R],
                 worst = rownames(draws)[use][worst], C = C, c_round = c_round,
                 c_member = c_member, eff = mean(e), moved = mean(blk$moved / blk$n),
                 pool = plan$pool, rounds = plan$rounds)
  sampler$ensemble_plan <- plan
  sampler$ensemble_plan_history <- c(sampler$ensemble_plan_history, list(record))
  sampler$ensemble_stats$block <- NULL
  sampler
}
