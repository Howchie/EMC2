# Allocation of the optional moves: shear rounds on group coordinates
# (R/shear_move.R) and subject HMC transitions (R/subject_hmc.R).
#
# A move is worth its cost only where it raises the effective sample size of
# the fit's slowest coordinate. Every iteration each chain records how far
# the rest of the iteration moved each coordinate:
#   group: (mu, log sd) at the end of the previous iteration -> just before
#     the shear;
#   subjects: alpha_s at the end of the previous iteration -> just before the
#     subject HMC round, the last update of the iteration;
# together with the CPU of the rest of the iteration, of each shear round
# and of each HMC transition. At a proposal rebuild each coordinate's lag-one
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
# subject the allocation would use without measured HMC moves first gets a
# bracketed step size and pilot transitions on a private copy of its state.
# The allocation changes only at rebuilds, as every adapted quantity does.

.EMC_ALLOC_MAX_ROUNDS <- 4
.EMC_ALLOC_MIN_PAIRS <- 10L
.EMC_ALLOC_SHEAR_PRIOR_G <- 0.5
.EMC_ALLOC_MIN_GAIN <- 0.9
# Rounds per iteration on a shear coordinate whose rounds are not yet
# measured: enough for a measurement in one block, too few to cost much.
.EMC_ALLOC_EXPLORE_ROUNDS <- 0.25

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

.emc_alloc_wanted <- function(sampler) {
  .emc_shear_active(sampler) ||
    (isTRUE(.emc_hmc_options()$enabled) && isTRUE(.emc_hmc_support(sampler)$ok))
}

.emc_alloc_lsd <- function(V) 0.5 * log(pmax(diag(as.matrix(V)), .Machine$double.xmin))

.emc_alloc_empty_window <- function(n, p) {
  list(n = 0L, nG = 0L, G = matrix(0, 2L * p, 2L * p), nA = 0L,
       A = replicate(n, matrix(0, p, p), simplify = FALSE),
       cpu_rest = 0, cpu_shear = 0, rounds = 0L,
       cpu_hmc = numeric(n), hmc_att = integer(n))
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

.emc_alloc_record_cost <- function(sampler, rest, shear, rounds, tasks) {
  w <- .emc_alloc_window(sampler)
  w$n <- w$n + 1L
  w$cpu_rest <- w$cpu_rest + rest
  w$cpu_shear <- w$cpu_shear + shear
  w$rounds <- w$rounds + as.integer(rounds)
  for (s in seq_along(tasks)) {
    tk <- tasks[[s]]
    if (is.null(tk)) next
    w$cpu_hmc[s] <- w$cpu_hmc[s] + tk$cpu
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

# Choose expected rounds x (group coordinates) and rates r (subjects).
#   units: list of list(kind = "fixed" | "shear" | "hmc", rho, g or h, cost),
#     rho a scalar (fixed or shear) or a vector of the subject's coordinates;
#     g the per-round factor of a shear coordinate; h the per-coordinate
#     factor of an HMC transition.
# Each unit's need for a target rho* is the least x (or r) that brings its
# rho to rho*; the objective tau(rho*) * (C_rest + sum cost * x) is evaluated
# on a fine grid of targets between the best reachable and the unallocated
# maximum, plus each unit's own level.
.emc_alloc_solve <- function(units, C_rest, max_rounds = .EMC_ALLOC_MAX_ROUNDS) {
  need <- function(u, rs) {
    r <- u$rho
    if (identical(u$kind, "fixed")) return(if (max(r) > rs) Inf else 0)
    if (identical(u$kind, "shear")) {
      if (r <= rs) return(0)
      g <- u$g
      if (g >= 1 - 1e-9) return(Inf)
      # Least x with g^m (1 - f + f g) <= rs / r.
      m <- floor(log(rs / r) / log(g) + 1e-9)
      base <- r * g^m
      x <- if (base <= rs * (1 + 1e-12)) m else m + (1 - rs / base) / (1 - g)
      cap <- if (is.null(u$cap)) max_rounds else u$cap
      return(if (x > cap + 1e-9) Inf else x)
    }
    act <- r > rs
    if (!any(act)) return(0)
    h <- u$h[act]
    if (any(h >= 1 - 1e-9)) return(Inf)
    x <- max((1 - rs / r[act]) / (1 - h))
    if (x > 1 + 1e-9) Inf else min(1, x)
  }
  best_rho <- vapply(units, function(u) switch(u$kind,
    fixed = max(u$rho),
    shear = u$rho * .emc_alloc_shear_factor(u$g, if (is.null(u$cap)) max_rounds else u$cap),
    hmc = max(u$rho * u$h)), numeric(1))
  hi <- max(vapply(units, function(u) max(u$rho), numeric(1)), 0)
  lo <- max(best_rho, 0)
  x0 <- numeric(length(units))
  zero <- list(x = x0, rho = hi, tau = .emc_alloc_tau(hi), cost = C_rest,
               objective = .emc_alloc_tau(hi) * C_rest)
  zero$objective0 <- zero$objective
  zero$tau0 <- zero$tau
  if (hi <= lo + 1e-12) return(zero)
  tau_grid <- exp(seq(log(.emc_alloc_tau(lo)), log(.emc_alloc_tau(hi)), length.out = 200L))
  levels <- vapply(units, function(u) max(u$rho), numeric(1))
  targets <- sort(unique(c((tau_grid - 1) / (tau_grid + 1), levels[levels >= lo & levels <= hi])))
  best <- zero
  for (rs in targets) {
    x <- vapply(units, need, numeric(1), rs = rs)
    if (any(!is.finite(x))) next
    cost <- C_rest + sum(x * vapply(units, `[[`, numeric(1), "cost"))
    obj <- .emc_alloc_tau(rs) * cost
    if (obj < best$objective * (1 - 1e-9))
      best <- list(x = x, rho = rs, tau = .emc_alloc_tau(rs), cost = cost, objective = obj)
  }
  # The lag-one surrogate is coarse, so a small predicted gain is not worth
  # its moves.
  if (best$objective > .EMC_ALLOC_MIN_GAIN * zero$objective) best <- zero
  best$objective0 <- zero$objective
  best$tau0 <- zero$tau
  best
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
  C_rest <- if (!stored && w$n > 0 && w$cpu_rest > 0) w$cpu_rest / w$n else NA_real_
  list(stored = stored, pairs = c(group = w$nG, subjects = w$nA), shear = shear,
       g_rho = g_rho, g_var = g_var, g_fac = g_fac, g_measured = g_measured,
       W = W, s_rho = s_rho,
       s_var = s_var, C_rest = C_rest,
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
  c_round <- if (is.finite(ms$c_round)) ms$c_round else 0.1 * C_rest
  units <- lapply(seq_along(ms$g_rho), function(k)
    if (shear_on && is.finite(ms$g_fac[k])) {
      list(kind = "shear", rho = ms$g_rho[k], g = ms$g_fac[k], cost = c_round,
           cap = if (ms$g_measured[k]) NULL else .EMC_ALLOC_EXPLORE_ROUNDS)
    } else list(kind = "fixed", rho = ms$g_rho[k], cost = 0))
  for (s in seq_len(n)) {
    sj <- if (is.null(hmc_state)) NULL else hmc_state$subjects[[s]]
    u <- list(kind = "fixed", rho = ms$s_rho[s, ], cost = 0)
    if (isTRUE(sj$eligible)) {
      ev <- hmc_state$evidence[[s]]
      h <- .emc_alloc_hmc_factor(ev, ms$W[[s]], ms$s_var[s, ], opts$min_evidence)
      if (is.null(h) && optimistic) h <- numeric(length(ms$s_rho[s, ]))
      if (!is.null(h)) {
        c_s <- if (is.finite(ms$c_hmc[s])) ms$c_hmc[s] else if (!is.null(ev) && ev$n > 0 &&
          ev$cpu > 0) ev$cpu / ev$n else 3 * C_rest / n
        u <- list(kind = "hmc", rho = ms$s_rho[s, ], h = h, cost = c_s)
      }
    }
    units[[length(units) + 1L]] <- u
  }
  list(units = units, C_rest = C_rest, c_round = c_round)
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
      x <- .emc_alloc_solve(u$units, u$C_rest)$x
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
    sol <- .emc_alloc_solve(u$units, u$C_rest)
    n_g <- length(m$g_rho)
    if (shear_on[j]) ch$shear$plan$rounds <- sol$x[seq_len(n_g)]
    rate <- sol$x[n_g + seq_len(ch$n_subjects)]
    if (!is.null(st) && isTRUE(st$supported)) {
      for (s in seq_len(ch$n_subjects)) {
        sj <- st$subjects[[s]]
        sj$rate <- if (!isTRUE(sj$eligible)) 0 else if (!is.null(opts$fixed_rate))
          as.numeric(opts$fixed_rate) else if (isTRUE(sj$warmup))
            .EMC_HMC_WARMUP_RATE else rate[s]
        st$subjects[[s]] <- sj
      }
    }
    rec <- list(
      stage = stage, B = B, iteration = ch$samples$idx,
      transitions = if (is.null(ch$sample_transitions)) NA_integer_ else
        as.integer(ch$sample_transitions),
      stored = m$stored, pairs = m$pairs, C_rest = u$C_rest, c_round = u$c_round,
      c_hmc = m$c_hmc, tau = sol$tau, tau0 = sol$tau0, cost = sol$cost,
      objective = sol$objective, objective0 = sol$objective0,
      group_rho = m$g_rho, group_g = m$g_fac,
      rounds = if (shear_on[j]) sol$x[seq_len(n_g)] else numeric(n_g),
      subject_rho = apply(m$s_rho, 1L, max),
      hmc_h = vapply(u$units[n_g + seq_len(ch$n_subjects)], function(x)
        if (identical(x$kind, "hmc")) max(x$h) else NA_real_, numeric(1)),
      rate = if (is.null(st$subjects)) numeric(ch$n_subjects) else
        vapply(st$subjects, function(x) x$rate, numeric(1)),
      eps = if (is.null(st$subjects)) rep(NA_real_, ch$n_subjects) else
        vapply(st$subjects, function(x) x$eps, numeric(1)),
      calib_cpu = if (is.null(st$subjects)) 0 else
        sum(vapply(st$subjects, function(x) x$calib_cpu, numeric(1))))
    ch$alloc <- list(window = .emc_alloc_empty_window(ch$n_subjects, ch$n_pars),
                     last = rec, history = c(ch$alloc$history, list(rec)))
    if (!is.null(st)) {
      st$alloc <- rec
      ch$subject_hmc <- st
    }
    emc[[j]] <- ch
  }
  emc
}
