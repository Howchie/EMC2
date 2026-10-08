# Allocation of the optional moves: shear rounds on group coordinates
# (R/shear_move.R), collapsed moves of the whole group state
# (R/collapsed_move.R) and subject HMC transitions (R/subject_hmc.R), and of
# how often each subject's particle step runs.
#
# A move is worth its cost only where it raises the effective sample size of
# the fit's slowest coordinate. Every iteration each chain records how far
# the rest of the iteration moved each coordinate:
#   group: (mu, log sd) at the end of the previous iteration -> just before
#     the shear;
#   subjects: alpha_s at the end of the previous iteration -> just before the
#     subject HMC round, the last update of the iteration;
# together with the cost of the rest of the iteration, of each shear round,
# collapsed move and HMC transition, all as wall time on the cores the chain
# holds (.emc_reserved_cost), so work that leaves the chain's workers idle is
# charged for them. At a proposal rebuild each coordinate's lag-one
# correlation under the rest of the iteration is rho = 1 - ESJD / (2 var),
# where var is the coordinate's variance over every chain's window draws: a
# chain held in one region makes jumps small against that pooled spread, and
# shows as slow even though its own draws look settled. An upper confidence
# bound on ESJD makes rho a lower bound, so noise in a short window buys no
# moves. Group coordinates are mu along the shear's mean directions and the
# log group sds; subject coordinates are alpha_s whitened in the reference
# frame of its likelihood fit.
#
# A shear round on group coordinate k multiplies its rho by g_k, and an HMC
# transition taken with probability r multiplies subject coordinate i's by
# 1 - r + r h_i; g and h come from the moves' own jumps, measured the same
# way (a round's selection-weighted squared move, an HMC transition's
# acceptance-weighted squared jump). With tau = (1 + rho) / (1 - rho), the
# allocation minimises
#   max over coordinates of tau * (C_rest + sum_k x_k c_round + sum_s r_s c_s)
# over expected rounds x_k and rates r_s (.emc_alloc_solve). Coordinates the
# moves cannot reach (no shear, no usable fit) keep their rest-of-iteration
# rho and act as a floor: nothing is spent below the slowest of them. A
# collapsed move moves every group coordinate at once: x_c expected moves
# multiply coordinate k's rho by g_c,k^x_c (composed like shear rounds), so
# the solver takes the best of a grid over x_c, each with the per-unit
# allocation on the reduced rhos and x_c c_coll added to C_rest. Its g_c is
# measured at the current rate and extrapolated geometrically, so x_c may
# grow to at most twice that rate per rebuild (or the shared round cap,
# whichever is larger).
#
# Each subject's particle step is a move like the others: it runs with
# probability q_s, multiplies coordinate i's rho by 1 - q_s + q_s p_i and
# costs q_s c_p,s, with p and c_p measured from the steps that ran. The rest
# of the iteration without it has rho t_i = rho_i / (1 - q (1 - p_i)) at the
# rate q in force, so HMC can stand in for the particle step where it moves
# the subject more per unit cost; each subject's (q_s, r_s) is the cheapest
# pair that brings its coordinates to the target. The group's measured
# mixing depends on how the subjects move, which the units do not model, so
# no coordinate may end slower than under its particle step alone. A subject the
# allocation would use without measured HMC moves first gets a bracketed
# step size and pilot transitions on a private copy of its state. The
# allocation changes only at rebuilds, as every adapted quantity does.

.EMC_ALLOC_MAX_ROUNDS <- 4
.EMC_ALLOC_MIN_PAIRS <- 10L
.EMC_ALLOC_SHEAR_PRIOR_G <- 0.5
.EMC_ALLOC_MIN_GAIN <- 0.9
# Rounds per iteration on a shear coordinate whose rounds are not yet
# measured: enough for a measurement in one block, too few to cost much.
.EMC_ALLOC_EXPLORE_ROUNDS <- 0.25
# Particle rates the solver considers, highest first so that ties keep the
# particle step.
.EMC_ALLOC_PARTICLE_RATES <- seq(1, 0, by = -0.125)

.emc_cpu_self <- function() {
  p <- proc.time()
  p[["user.self"]] + p[["sys.self"]]
}

# Worker CPU of every subject likelihood round so far (see
# .emc_wpool_subject_ll).
.emc_subject_ll_cpu <- function() {
  x <- .emc_pool_state$subject_ll_cpu
  if (is.null(x)) 0 else x
}

.emc_cpu_children <- function() {
  p <- proc.time()
  p[["user.child"]] + p[["sys.child"]]
}

# Every cost the allocation compares is wall time on the cores the chain holds:
# a move's computation in the chain's own process and its likelihood rounds'
# round trips keep those cores as busy as they leave them.
.emc_reserved_cost <- function(elapsed, cores) {
  max(1, as.numeric(cores)) * max(0, as.numeric(elapsed))
}

.emc_alloc_wanted <- function(sampler) {
  .emc_shear_active(sampler) || .emc_collapsed_active(sampler) ||
    (isTRUE(.emc_hmc_options()$enabled) && isTRUE(.emc_hmc_support(sampler)$ok))
}

.emc_alloc_lsd <- function(V) 0.5 * log(pmax(diag(as.matrix(V)), .Machine$double.xmin))

.emc_alloc_empty_window <- function(n, p) {
  list(n = 0L, nG = 0L, G = matrix(0, 2L * p, 2L * p), nA = 0L,
       A = replicate(n, matrix(0, p, p), simplify = FALSE),
       cpu_rest = 0, cpu_shear = 0, rounds = 0L, cpu_coll = 0, coll_moves = 0L,
       cpu_hmc = numeric(n), hmc_att = integer(n),
       P = replicate(n, matrix(0, p, p), simplify = FALSE), nP = integer(n),
       cpu_part = numeric(n))
}

# The particle step's own jumps and cost. `cost` is the step's wall time on
# the chain's cores, shared among the subjects that ran by their worker CPU
# `times` (part of the rest-of-iteration cost recorded by
# .emc_alloc_record_cost).
.emc_alloc_record_particle <- function(sampler, before, after, ran, cost, times) {
  w <- .emc_alloc_window(sampler)
  if (is.null(w$nP)) {
    e <- .emc_alloc_empty_window(sampler$n_subjects, sampler$n_pars)
    w$P <- e$P; w$nP <- e$nP; w$cpu_part <- e$cpu_part
  }
  idx <- which(ran)
  if (!length(idx)) return(sampler)
  t_ran <- pmax(0, times[idx])
  share <- if (sum(t_ran) > 0) t_ran / sum(t_ran) else rep(1 / length(idx), length(idx))
  for (k in seq_along(idx)) {
    s <- idx[k]
    d <- after[, s] - before[, s]
    if (!all(is.finite(d))) next
    w$P[[s]] <- w$P[[s]] + tcrossprod(d)
    w$nP[s] <- w$nP[s] + 1L
    w$cpu_part[s] <- w$cpu_part[s] + cost * share[k]
  }
  sampler$alloc$window <- w
  sampler
}

.emc_alloc_window <- function(sampler) {
  w <- sampler$alloc$window
  if (is.null(w)) w <- .emc_alloc_empty_window(sampler$n_subjects, sampler$n_pars)
  w
}

# The group state the rest of the iteration produced, against the previous
# iteration's final (stored) one.
.emc_alloc_record_group <- function(sampler, pars) {
  idx <- sampler$samples$idx
  w <- .emc_alloc_window(sampler)
  d <- c(as.numeric(pars$tmu), .emc_alloc_lsd(pars$tvar)) -
    c(sampler$samples$theta_mu[, idx], .emc_alloc_lsd(sampler$samples$theta_var[, , idx]))
  if (all(is.finite(d))) {
    w$G <- w$G + tcrossprod(d)
    w$nG <- w$nG + 1L
  }
  sampler$alloc$window <- w
  sampler
}

# The subjects just before the HMC round, against the previous iteration's
# final ones.
.emc_alloc_record_subjects <- function(sampler, alpha) {
  idx <- sampler$samples$idx
  w <- .emc_alloc_window(sampler)
  d <- alpha - matrix(sampler$samples$alpha[, , idx], nrow(alpha))
  if (all(is.finite(d))) {
    for (s in seq_len(ncol(d))) w$A[[s]] <- w$A[[s]] + tcrossprod(d[, s])
    w$nA <- w$nA + 1L
  }
  sampler$alloc$window <- w
  sampler
}

.emc_alloc_record_cost <- function(sampler, rest, shear, rounds, tasks, coll = 0,
                                   coll_moves = 0L) {
  w <- .emc_alloc_window(sampler)
  w$n <- w$n + 1L
  w$cpu_rest <- w$cpu_rest + rest
  w$cpu_shear <- w$cpu_shear + shear
  w$rounds <- w$rounds + as.integer(rounds)
  if (is.null(w$cpu_coll)) w$cpu_coll <- w$coll_moves <- 0
  w$cpu_coll <- w$cpu_coll + coll
  w$coll_moves <- w$coll_moves + as.integer(coll_moves)
  for (s in seq_along(tasks)) {
    tk <- tasks[[s]]
    if (is.null(tk)) next
    w$cpu_hmc[s] <- w$cpu_hmc[s] + if (is.null(tk$cost)) tk$cpu else tk$cost
    w$hmc_att[s] <- w$hmc_att[s] + 1L
  }
  sampler$alloc$window <- w
  sampler
}

# Lower bound on rho from a mean squared jump over n pairs and the
# coordinate's variance.
.emc_alloc_rho_lo <- function(esjd, var, n) {
  m <- esjd / (2 * pmax(var, .Machine$double.xmin))
  m_hi <- m * (1 + 2 * sqrt(2 / max(1, n)))
  pmin(1, pmax(0, 1 - m_hi))
}

.emc_alloc_tau <- function(rho) {
  rho <- pmin(rho, 1 - 1e-12)
  (1 + rho) / (1 - rho)
}

# Factor on a shear coordinate's rho from x expected rounds per iteration
# (floor(x) rounds, and one more with probability x - floor(x)), given the
# factor g of one round.
.emc_alloc_shear_factor <- function(g, x) {
  m <- floor(x)
  f <- x - m
  g^m * (1 - f + f * g)
}

# Candidate expected collapsed moves per iteration: eighths up to 4, then
# halves, then the cap itself.
.emc_alloc_joint_grid <- function(cap) {
  x <- seq(0.125, 4, by = 0.125)
  if (cap >= 4.5) x <- c(x, seq(4.5, cap, by = 0.5))
  sort(unique(c(x[x <= cap], cap)))
}

# Choose expected rounds x (group coordinates), and HMC rates r and particle
# rates q (subjects).
#   units: list of list(kind = "fixed" | "shear" | "subject" | "hmc", ...):
#     fixed: rho, a scalar or vector no move changes;
#     shear: rho, g the per-round factor, cost per round, cap (optional);
#     subject: rho (t, the iteration without the particle step), p the
#       particle step's factors or NULL (always run, cost in C_rest), c_p,
#       h the HMC factors or NULL, c_h;
#     hmc: rho, h, cost; a subject unit without a particle factor.
# Each unit's cost for a target rho* is that of the cheapest allocation that
# brings its rho to rho*; the objective tau(rho*) * (C_rest + sum costs) is
# evaluated on one fine grid of targets between the best reachable and the
# unallocated maximum, plus every unit's own level under every candidate
# collapsed rate. The units the collapsed move does not touch are priced once
# on that grid; only the group units are re-priced per candidate rate.
.emc_alloc_solve <- function(units, C_rest, max_rounds = .EMC_ALLOC_MAX_ROUNDS,
                             joint = NULL) {
  units <- lapply(units, .emc_alloc_as_subject)
  n_u <- length(units)
  moved_idx <- if (is.null(joint)) integer(0) else as.integer(joint$idx)
  xcs <- c(0, if (!is.null(joint)) .emc_alloc_joint_grid(joint$cap))
  cand <- lapply(xcs, function(xc) {
    if (xc == 0) return(units)
    fac <- .emc_alloc_shear_factor(joint$g, xc)
    moved <- units
    for (i in seq_along(moved_idx))
      moved[[moved_idx[i]]]$rho <- moved[[moved_idx[i]]]$rho * fac[i]
    moved
  })
  ranges <- lapply(cand, .emc_alloc_range, max_rounds = max_rounds)
  r0 <- ranges[[1L]]
  c_p0 <- sum(vapply(units, function(u)
    if (identical(u$kind, "subject") && !is.null(u$p) && !is.null(u$c_p)) u$c_p else 0,
    numeric(1)))
  zero <- list(x = numeric(n_u), q = rep(1, n_u), rho = r0$hi, tau = .emc_alloc_tau(r0$hi),
               cost = C_rest + c_p0, objective = .emc_alloc_tau(r0$hi) * (C_rest + c_p0),
               joint = 0)
  lo <- min(vapply(ranges, `[[`, numeric(1), "lo"))
  best <- zero
  if (r0$hi > lo + 1e-12 || c_p0 > 0) {
    tau_grid <- exp(seq(log(.emc_alloc_tau(lo)), log(.emc_alloc_tau(r0$hi)),
                        length.out = 400L))
    lv <- unlist(lapply(ranges, function(r) r$levels[r$levels >= r$lo & r$levels <= r$hi]))
    targets <- sort(unique(c((tau_grid - 1) / (tau_grid + 1), lv)))
    tau_t <- .emc_alloc_tau(targets)
    still <- setdiff(seq_len(n_u), moved_idx)
    need <- vector("list", n_u)
    need[still] <- lapply(units[still], .emc_alloc_unit_need, rs = targets,
                          max_rounds = max_rounds)
    base <- C_rest + Reduce(`+`, lapply(need[still], `[[`, "cost"), numeric(length(targets)))
    for (c in seq_along(xcs)) {
      nm <- lapply(cand[[c]][moved_idx], .emc_alloc_unit_need, rs = targets,
                   max_rounds = max_rounds)
      cost <- base + xcs[c] * (if (xcs[c] > 0) joint$cost else 0) +
        Reduce(`+`, lapply(nm, `[[`, "cost"), numeric(length(targets)))
      obj <- tau_t * cost
      k <- which.min(obj)
      if (!length(k) || !is.finite(obj[k]) || obj[k] >= best$objective * (1 - 1e-9)) next
      need[moved_idx] <- nm
      best <- list(x = vapply(need, function(z) z$x[k], numeric(1)),
                   q = vapply(need, function(z) z$q[k], numeric(1)),
                   rho = targets[k], tau = tau_t[k], cost = cost[k],
                   objective = obj[k], joint = xcs[c])
    }
  }
  # The lag-one surrogate is coarse, so a small predicted gain is not worth
  # its moves.
  if (best$objective > .EMC_ALLOC_MIN_GAIN * zero$objective) best <- zero
  best$objective0 <- zero$objective
  best$tau0 <- zero$tau
  best
}

# The unallocated rho of each unit, and the least rho* every unit can reach.
.emc_alloc_range <- function(units, max_rounds) {
  levels <- vapply(units, function(u) max(.emc_alloc_level(u)), numeric(1))
  best_rho <- vapply(units, function(u) switch(u$kind,
    fixed = max(u$rho),
    shear = u$rho * .emc_alloc_shear_factor(u$g, if (is.null(u$cap)) max_rounds else u$cap),
    subject = max(.emc_alloc_level(u) * if (is.null(u$h)) 1 else u$h)), numeric(1))
  list(levels = levels, hi = max(levels, 0), lo = max(best_rho, 0))
}

.emc_alloc_as_subject <- function(u) {
  if (!identical(u$kind, "hmc")) return(u)
  list(kind = "subject", rho = u$rho, p = NULL, c_p = 0, h = u$h, c_h = u$cost)
}

# A unit's rho with nothing allocated: subjects run their particle step.
.emc_alloc_level <- function(u) {
  if (identical(u$kind, "subject") && !is.null(u$p)) u$rho * u$p else u$rho
}

# A unit's allocation for each target rho* in rs: list(x, q, cost), with
# cost Inf where the target is out of its reach.
.emc_alloc_unit_need <- function(u, rs, max_rounds) {
  n_t <- length(rs)
  if (identical(u$kind, "fixed"))
    return(list(x = numeric(n_t), q = rep(1, n_t),
                cost = ifelse(max(u$rho) > rs, Inf, 0)))
  if (identical(u$kind, "shear")) {
    r <- u$rho
    g <- u$g
    cap <- if (is.null(u$cap)) max_rounds else u$cap
    x <- numeric(n_t)
    over <- r > rs
    if (any(over)) {
      if (g >= 1 - 1e-9) {
        x[over] <- Inf
      } else {
        ro <- rs[over]
        # Least x with g^m (1 - f + f g) <= rs / r.
        m <- floor(log(ro / r) / log(g) + 1e-9)
        base <- r * g^m
        xo <- ifelse(base <= ro * (1 + 1e-12), m, m + (1 - ro / base) / (1 - g))
        xo[!is.finite(xo) | xo > cap + 1e-9] <- Inf
        x[over] <- xo
      }
    }
    return(list(x = x, q = rep(1, n_t), cost = ifelse(is.finite(x), x * u$cost, Inf)))
  }
  # Subject: for each particle rate q, the least HMC rate r with
  # rho_i (1 - q (1 - p_i)) (1 - r (1 - h_i)) <= min(rho*, level_i) on every
  # coordinate; then the cheapest q.
  rho <- u$rho
  level <- .emc_alloc_level(u)
  qs <- if (is.null(u$p) || is.null(u$h)) 1 else .EMC_ALLOC_PARTICLE_RATES
  r_need <- matrix(0, length(qs), n_t)
  for (i in seq_along(rho)) {
    K <- pmin(rs, level[i])
    base <- if (is.null(u$p)) rep(rho[i], length(qs)) else rho[i] * (1 - qs * (1 - u$p[i]))
    over <- outer(base, K * (1 + 1e-12), ">")
    if (!any(over)) next
    if (is.null(u$h)) {
      r_need[over] <- Inf
      next
    }
    val <- (1 - outer(1 / base, K)) / (1 - u$h[i])
    r_need <- pmax(r_need, ifelse(over, val, 0))
  }
  r_need[r_need > 1 + 1e-9] <- Inf
  fin <- is.finite(r_need)
  r_need[fin] <- pmin(r_need[fin], 1)
  c_h <- if (is.null(u$c_h)) 0 else u$c_h
  c_p <- if (is.null(u$c_p)) 0 else u$c_p
  cost <- qs * c_p + r_need * c_h
  cost[!fin] <- Inf
  k <- max.col(t(-cost), ties.method = "first")
  pick <- cbind(k, seq_len(n_t))
  list(x = r_need[pick], q = qs[k], cost = cost[pick])
}

# The last m stored draws of a chain (at least two).
.emc_alloc_draws <- function(chain, m) {
  i <- chain$samples$idx
  seq.int(max(1L, i - max(2L, m) + 1L), i)
}

# A chain's measurements at a rebuild, against every chain's window draws.
.emc_alloc_measure <- function(emc, j, hmc_state, Vbar) {
  ch <- emc[[j]]
  p <- ch$n_pars
  n <- ch$n_subjects
  w <- ch$alloc$window
  stored <- is.null(w) || w$n < .EMC_ALLOC_MIN_PAIRS
  m <- if (stored) 200L else w$n
  draws <- lapply(emc, .emc_alloc_draws, m = m)
  if (stored) {
    # Whole-iteration jumps between consecutive stored draws stand in for a
    # window not yet recorded.
    k <- draws[[j]]
    w <- .emc_alloc_empty_window(n, p)
    if (length(k) >= 2L) {
      s <- ch$samples
      for (t in k[-1L]) {
        d <- c(s$theta_mu[, t], .emc_alloc_lsd(s$theta_var[, , t])) -
          c(s$theta_mu[, t - 1L], .emc_alloc_lsd(s$theta_var[, , t - 1L]))
        if (all(is.finite(d))) { w$G <- w$G + tcrossprod(d); w$nG <- w$nG + 1L }
        da <- matrix(s$alpha[, , t] - s$alpha[, , t - 1L], p)
        if (all(is.finite(da))) {
          for (ss in seq_len(n)) w$A[[ss]] <- w$A[[ss]] + tcrossprod(da[, ss])
          w$nA <- w$nA + 1L
        }
      }
    }
  }
  # Group coordinates: mu along the shear's directions (or mu itself), then
  # the log group sds.
  shear <- .emc_shear_active(ch) && is.list(ch$shear$plan) && is.matrix(ch$shear$plan$dirs)
  dirs <- if (shear) ch$shear$plan$dirs else diag(p)
  Dinv <- tryCatch(solve(dirs), error = function(e) NULL)
  if (is.null(Dinv)) {
    shear <- FALSE
    dirs <- diag(p)
    Dinv <- diag(p)
  }
  T_g <- rbind(cbind(Dinv, matrix(0, p, p)), cbind(matrix(0, p, p), diag(p)))
  gx <- do.call(cbind, lapply(seq_along(emc), function(i) {
    s <- emc[[i]]$samples
    vapply(draws[[i]], function(t) c(drop(Dinv %*% s$theta_mu[, t]),
                                     .emc_alloc_lsd(s$theta_var[, , t])), numeric(2L * p))
  }))
  g_var <- apply(gx, 1L, stats::var)
  g_esjd <- if (w$nG > 0) rowSums((T_g %*% w$G) * T_g) / w$nG else rep(0, 2L * p)
  g_rho <- if (w$nG > 0) .emc_alloc_rho_lo(g_esjd, g_var, w$nG) else rep(0, 2L * p)
  g_fac <- rep(NA_real_, 2L * p)
  g_measured <- logical(2L * p)
  if (shear) {
    last <- ch$shear$plan$last
    g_measured <- last$n >= 5L
    g_fac <- ifelse(g_measured,
                    pmin(1 - 1e-6, pmax(0, 1 - (last$E / pmax(1L, last$n)) / (2 * g_var))),
                    .EMC_ALLOC_SHEAR_PRIOR_G)
  }
  # The collapsed move moves every group coordinate at once: its factors come
  # from its acceptance-weighted jumps in the same coordinates.
  coll <- .emc_collapsed_active(ch) && is.list(ch$collapsed$plan) &&
    is.matrix(ch$collapsed$plan$C)
  c_fac <- rep(.EMC_ALLOC_SHEAR_PRIOR_G, 2L * p)
  c_measured <- FALSE
  if (coll) {
    last <- ch$collapsed$plan$last
    c_measured <- is.list(last) && last$n >= 5L
    if (c_measured)
      c_fac <- pmin(1 - 1e-6, pmax(0, 1 - (rowSums((T_g %*% last$E) * T_g) / last$n) /
                                     (2 * g_var)))
  }
  # Subject coordinates: alpha_s whitened in its fit's reference frame.
  lik <- ch$lik_laplace
  W <- lapply(seq_len(n), function(s) {
    if (!is.null(hmc_state$subjects[[s]]$W)) return(hmc_state$subjects[[s]]$W)
    prec <- if (length(lik) >= s && is.list(lik[[s]])) lik[[s]]$prec else NULL
    Ws <- if (is.matrix(prec) && all(dim(prec) == p) && all(is.finite(prec)))
      .emc_hmc_ref_whiten(prec, Vbar)
    if (is.null(Ws)) diag(p) else Ws
  })
  s_rho <- matrix(0, n, p)
  s_var <- matrix(NA_real_, n, p)
  for (s in seq_len(n)) {
    ax <- do.call(cbind, lapply(seq_along(emc), function(i)
      W[[s]] %*% matrix(emc[[i]]$samples$alpha[, s, draws[[i]]], p)))
    s_var[s, ] <- apply(ax, 1L, stats::var)
    if (w$nA > 0) s_rho[s, ] <- .emc_alloc_rho_lo(
      rowSums((W[[s]] %*% w$A[[s]]) * W[[s]]) / w$nA, s_var[s, ], w$nA)
  }
  # Particle steps: their own factors and cost per step, carried over from
  # the last rebuild while the window holds too few of a subject's steps.
  prev <- ch$alloc$particle
  nP <- if (stored || is.null(w$nP)) integer(n) else w$nP
  p_fac <- matrix(NA_real_, n, p)
  c_part <- rep(NA_real_, n)
  for (s in seq_len(n)) {
    if (nP[s] >= .EMC_ALLOC_MIN_PAIRS) {
      esjd <- rowSums((W[[s]] %*% w$P[[s]]) * W[[s]]) / nP[s]
      p_fac[s, ] <- pmin(1, pmax(0, 1 - esjd / (2 * pmax(s_var[s, ], .Machine$double.xmin))))
      c_part[s] <- w$cpu_part[s] / nP[s]
    } else if (is.matrix(prev$fac) && all(dim(prev$fac) == c(n, p)) &&
               all(is.finite(prev$fac[s, ])) && is.finite(prev$cost[s])) {
      p_fac[s, ] <- prev$fac[s, ]
      c_part[s] <- prev$cost[s]
    }
  }
  q_old <- if (is.null(prev$rate)) rep(1, n) else prev$rate
  t_rho <- s_rho
  for (s in which(is.finite(c_part))) {
    den <- 1 - q_old[s] * (1 - p_fac[s, ])
    t_rho[s, ] <- ifelse(den > 1e-6, pmin(1, s_rho[s, ] / pmax(den, 1e-6)), 1)
  }
  # The rest of the iteration, less the particle steps the units now price.
  part_in_rest <- if (stored || is.null(w$cpu_part)) 0 else
    sum(w$cpu_part[is.finite(c_part)])
  C_rest <- if (!stored && w$n > 0 && w$cpu_rest > part_in_rest)
    (w$cpu_rest - part_in_rest) / w$n else NA_real_
  list(stored = stored, pairs = c(group = w$nG, subjects = w$nA), shear = shear,
       g_rho = g_rho, g_var = g_var, g_fac = g_fac, g_measured = g_measured,
       coll = coll, c_fac = c_fac, c_measured = c_measured,
       c_coll = if (!stored && isTRUE(w$coll_moves >= 5L)) w$cpu_coll / w$coll_moves
                else NA_real_,
       c_rate = if (!stored && w$n > 0) w$coll_moves / w$n else 0,
       W = W, s_rho = s_rho, t_rho = t_rho, p_fac = p_fac, c_part = c_part,
       q_old = q_old, s_var = s_var, C_rest = C_rest,
       c_round = if (!stored && w$rounds >= 5L) w$cpu_shear / w$rounds else NA_real_,
       c_hmc = if (!stored) ifelse(w$hmc_att > 0, w$cpu_hmc / pmax(1L, w$hmc_att), NA_real_)
               else rep(NA_real_, n))
}

# HMC factors h for subject s from its evidence; NULL without enough.
.emc_alloc_hmc_factor <- function(ev, W, var, min_n) {
  if (is.null(ev) || ev$n < min_n) return(NULL)
  esjd <- rowSums((W %*% ev$M) * W) / ev$n
  pmin(1, pmax(0, 1 - esjd / (2 * pmax(var, .Machine$double.xmin))))
}

# The units and costs the solver sees for one chain. optimistic: an eligible
# subject without enough HMC evidence counts as h = 0.
.emc_alloc_units <- function(ms, hmc_state, shear_on, optimistic) {
  n <- nrow(ms$s_rho)
  opts <- .emc_hmc_options()
  C_rest <- if (is.finite(ms$C_rest)) ms$C_rest else 1
  # Unmeasured costs are guessed against the whole iteration's.
  C_iter <- C_rest + sum(ms$c_part[is.finite(ms$c_part)])
  c_round <- if (is.finite(ms$c_round)) ms$c_round else 0.1 * C_iter
  units <- lapply(seq_along(ms$g_rho), function(k)
    if (shear_on && is.finite(ms$g_fac[k])) {
      list(kind = "shear", rho = ms$g_rho[k], g = ms$g_fac[k], cost = c_round,
           cap = if (ms$g_measured[k]) NULL else .EMC_ALLOC_EXPLORE_ROUNDS)
    } else list(kind = "fixed", rho = ms$g_rho[k], cost = 0))
  for (s in seq_len(n)) {
    sj <- if (is.null(hmc_state)) NULL else hmc_state$subjects[[s]]
    # A measured particle step is priced in the subject's unit (and left out
    # of C_rest); without an HMC rate to stand in for it, it always runs.
    priced <- !is.null(ms$c_part) && is.finite(ms$c_part[s])
    u <- list(kind = "subject", rho = ms$s_rho[s, ], p = NULL, c_p = 0, h = NULL, c_h = 0)
    if (isTRUE(sj$eligible)) {
      ev <- hmc_state$evidence[[s]]
      h <- .emc_alloc_hmc_factor(ev, ms$W[[s]], ms$s_var[s, ], opts$min_evidence)
      if (is.null(h) && optimistic) h <- numeric(length(ms$s_rho[s, ]))
      if (!is.null(h)) {
        c_s <- if (is.finite(ms$c_hmc[s])) ms$c_hmc[s] else if (!is.null(ev) && ev$n > 0 &&
          ev$cpu > 0) ev$cpu / ev$n else 3 * C_iter / n
        u$h <- h
        u$c_h <- c_s
      }
    }
    if (priced) {
      u$rho <- ms$t_rho[s, ]
      u$p <- ms$p_fac[s, ]
      u$c_p <- ms$c_part[s]
    }
    units[[length(units) + 1L]] <- u
  }
  joint <- if (isTRUE(ms$coll)) {
    list(idx = seq_along(ms$g_rho), g = ms$c_fac,
         cost = if (is.finite(ms$c_coll)) ms$c_coll else 0.1 * C_iter,
         cap = if (ms$c_measured) max(.EMC_ALLOC_MAX_ROUNDS, 2 * ms$c_rate) else
           .EMC_ALLOC_EXPLORE_ROUNDS)
  }
  list(units = units, C_rest = C_rest, c_round = c_round, joint = joint)
}

# Rebuild every chain's shear rounds and HMC rates (in add_proposals(), after
# the likelihood fits are checked and the shear plans built). B is the number
# of transitions until the next rebuild.
.emc_alloc_rebuild <- function(emc, stage, n_cores, B) {
  if (identical(stage, "preburn") || !.emc_alloc_wanted(emc[[1L]])) return(emc)
  emc <- .emc_prepare_chain_rng(emc)
  opts <- .emc_hmc_options()
  hmc_on <- isTRUE(opts$enabled)
  nch <- length(emc)
  Vbar <- lapply(emc, function(ch) {
    k <- .emc_alloc_draws(ch, if (is.null(ch$alloc$window)) 200L else
      max(2L, ch$alloc$window$n))
    apply(ch$samples$theta_var[, , k, drop = FALSE], 1:2, mean)
  })
  hst <- lapply(seq_len(nch), function(j)
    if (hmc_on) .emc_hmc_plan_chain(emc[[j]], stage, Vbar[[j]]) else NULL)
  ms <- lapply(seq_len(nch), function(j) .emc_alloc_measure(emc, j, hst[[j]], Vbar[[j]]))
  shear_on <- vapply(seq_len(nch), function(j) ms[[j]]$shear, logical(1))
  # Calibration for subjects the optimistic allocation would use but that
  # lack measured HMC moves.
  if (hmc_on) {
    jobs <- NULL
    for (j in seq_len(nch)) {
      if (!isTRUE(hst[[j]]$supported)) next
      u <- .emc_alloc_units(ms[[j]], hst[[j]], shear_on[j], optimistic = TRUE)
      x <- .emc_alloc_solve(u$units, u$C_rest, joint = u$joint)$x
      n_g <- length(ms[[j]]$g_rho)
      for (s in seq_len(emc[[j]]$n_subjects)) {
        sj <- hst[[j]]$subjects[[s]]
        if (!isTRUE(sj$eligible) || isTRUE(sj$warmup)) next
        wanted <- if (!is.null(opts$fixed_rate)) as.numeric(opts$fixed_rate) > 0 else
          x[n_g + s] > 0
        if (!wanted) next
        ev <- hst[[j]]$evidence[[s]]
        lacking <- is.null(ev) || ev$n < opts$min_evidence
        bracket <- is.null(opts$fixed_eps) && is.na(sj$eps)
        pilots <- if (lacking && is.null(opts$fixed_rate)) opts$pilots else 0L
        if (bracket || pilots > 0L)
          jobs <- rbind(jobs, data.frame(chain = j, subject = s, bracket = bracket,
                                         pilots = pilots))
      }
    }
    if (!is.null(jobs)) {
      res <- .emc_with_preserved_rng(auto_mclapply(seq_len(nrow(jobs)), function(k) {
        j <- jobs$chain[k]
        .emc_quiet_fork(.emc_hmc_calibrate(emc[[j]], hst[[j]], jobs$subject[k],
                                           bracket = jobs$bracket[k],
                                           pilots = jobs$pilots[k]))
      }, mc.cores = min(nrow(jobs), max(1L, n_cores)), mc.set.seed = FALSE))
      for (k in seq_len(nrow(jobs))) {
        j <- jobs$chain[k]
        s <- jobs$subject[k]
        r <- res[[k]]
        if (!is.list(r) || is.null(r$seed)) next
        emc[[j]]$rng$subjects[[s]] <- r$seed
        hst[[j]] <- .emc_hmc_apply_calibration(hst[[j]], s, r)
      }
    }
  }
  for (j in seq_len(nch)) {
    ch <- emc[[j]]
    m <- ms[[j]]
    st <- hst[[j]]
    if (!is.null(st) && isTRUE(st$supported)) for (s in seq_len(ch$n_subjects)) {
      if (isTRUE(st$subjects[[s]]$eligible) && !is.finite(st$subjects[[s]]$eps)) {
        st$subjects[[s]]$eligible <- FALSE
        st$subjects[[s]]$reason <- "no usable step size"
      }
    }
    u <- .emc_alloc_units(m, st, shear_on[j], optimistic = FALSE)
    sol <- .emc_alloc_solve(u$units, u$C_rest, joint = u$joint)
    n_g <- length(m$g_rho)
    if (shear_on[j]) ch$shear$plan$rounds <- sol$x[seq_len(n_g)]
    if (isTRUE(m$coll)) ch$collapsed$plan$rounds <- sol$joint
    rate <- sol$x[n_g + seq_len(ch$n_subjects)]
    q <- rep(1, ch$n_subjects)
    if (!is.null(st) && isTRUE(st$supported)) {
      for (s in seq_len(ch$n_subjects)) {
        sj <- st$subjects[[s]]
        sj$rate <- if (!isTRUE(sj$eligible)) 0 else if (!is.null(opts$fixed_rate))
          as.numeric(opts$fixed_rate) else if (isTRUE(sj$warmup))
            .EMC_HMC_WARMUP_RATE else rate[s]
        # The particle step stands down only for the rate chosen with it.
        if (isTRUE(sj$eligible) && is.null(opts$fixed_rate) && !isTRUE(sj$warmup) &&
            length(sol$q) >= n_g + s)
          q[s] <- sol$q[n_g + s]
        st$subjects[[s]] <- sj
      }
    }
    rec <- list(
      stage = stage, B = B, iteration = ch$samples$idx,
      transitions = if (is.null(ch$sample_transitions)) NA_integer_ else
        as.integer(ch$sample_transitions),
      stored = m$stored, pairs = m$pairs, C_rest = u$C_rest, c_round = u$c_round,
      c_hmc = m$c_hmc, c_part = m$c_part, tau = sol$tau, tau0 = sol$tau0, cost = sol$cost,
      objective = sol$objective, objective0 = sol$objective0,
      group_rho = m$g_rho, group_g = m$g_fac,
      rounds = if (shear_on[j]) sol$x[seq_len(n_g)] else numeric(n_g),
      collapsed_rounds = sol$joint, collapsed_g = if (isTRUE(m$coll)) m$c_fac,
      collapsed_cost = if (is.null(u$joint)) NA_real_ else u$joint$cost,
      collapsed_eps = if (isTRUE(m$coll)) ch$collapsed$plan$eps else NA_real_,
      subject_rho = apply(m$s_rho, 1L, max),
      hmc_h = vapply(u$units[n_g + seq_len(ch$n_subjects)], function(x)
        if (!is.null(x$h)) max(x$h) else NA_real_, numeric(1)),
      rate = if (is.null(st$subjects)) numeric(ch$n_subjects) else
        vapply(st$subjects, function(x) x$rate, numeric(1)),
      particle_rate = q,
      particle_p = apply(m$p_fac, 1L, function(x) if (all(is.finite(x))) max(x) else NA_real_),
      eps = if (is.null(st$subjects)) rep(NA_real_, ch$n_subjects) else
        vapply(st$subjects, function(x) x$eps, numeric(1)),
      calib_cpu = if (is.null(st$subjects)) 0 else
        sum(vapply(st$subjects, function(x) x$calib_cpu, numeric(1))))
    ch$alloc <- list(window = .emc_alloc_empty_window(ch$n_subjects, ch$n_pars),
                     last = rec, history = c(ch$alloc$history, list(rec)),
                     particle = list(rate = q, fac = m$p_fac, cost = m$c_part))
    if (!is.null(st)) {
      st$alloc <- rec
      ch$subject_hmc <- st
    }
    emc[[j]] <- ch
  }
  # One particle rate per subject for every chain, the largest any chain's
  # allocation chose, so chains carry the same work.
  q <- do.call(pmax, lapply(emc, function(ch) ch$alloc$particle$rate))
  for (j in seq_len(nch)) {
    emc[[j]]$alloc$particle$rate <- q
    emc[[j]]$alloc$last$particle_shared <- q
    h <- emc[[j]]$alloc$history
    h[[length(h)]]$particle_shared <- q
    emc[[j]]$alloc$history <- h
  }
  # One collapsed rate for every chain, the largest any chain's allocation
  # chose: its work runs in each chain's own process, and chains wait for
  # each other at block ends.
  coll <- vapply(seq_len(nch), function(j) isTRUE(ms[[j]]$coll), logical(1))
  if (any(coll)) {
    x_c <- max(vapply(emc[coll], function(ch) ch$collapsed$plan$rounds, numeric(1)))
    for (j in which(coll)) {
      emc[[j]]$collapsed$plan$rounds <- x_c
      emc[[j]]$alloc$last$collapsed_shared <- x_c
      h <- emc[[j]]$alloc$history
      h[[length(h)]]$collapsed_shared <- x_c
      emc[[j]]$alloc$history <- h
    }
  }
  emc
}
