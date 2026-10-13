# Brute-force oracle for the allocator's own model (R/move_allocation.R), for
# small unit sets: the objective a plan really has under the model, with the
# group coordinates credited at the subjects' rates the plan itself chooses.
# Used by the known-answer tests and by oracle_gap.R.

ws1_ns <- function() asNamespace("EMC2")

# The model's objective for a plan: x (per group unit), q and r (per subject
# unit, r the total HMC rate), xc collapsed moves. units as .emc_alloc_units
# returns them (shear/fixed group units first, then subjects); u carries
# g_floor, joint, shear_fixed, hmc_fixed, C_rest. movable: the coordinates the
# minimax ranges over (default: every coordinate some move or credit lowers).
# floor = FALSE scores a plan that breaks the particle-step floor instead of
# rejecting it, with the number of subjects that break it in attr "floor".
ws1_model_objective <- function(u, x, q, r, xc = 0, movable = NULL, credit = TRUE,
                                floor = TRUE) {
  ns <- ws1_ns()
  units <- u$units
  n_u <- length(units)
  n_g <- if (is.null(u$g_floor)) sum(vapply(units, function(z) z$kind %in% c("shear", "fixed"), logical(1))) else length(u$g_floor)
  qq <- rep(1, n_u); rr <- numeric(n_u)
  sidx <- (n_g + 1L):n_u
  if (n_u > n_g) { qq[sidx] <- q; rr[sidx] <- r }
  if (credit && !is.null(u$g_floor)) {
    R <- ns$.emc_alloc_dev_ratio(units, qq, rr)
    if (is.finite(R)) units[seq_len(n_g)] <- ns$.emc_alloc_credit(units[seq_len(n_g)], u$g_floor, R)
  }
  sol <- list(x = c(x, numeric(n_u - n_g)), q = qq, joint = xc)
  # .emc_alloc_predict takes the HMC rate above r_min in x.
  for (i in sidx) sol$x[i] <- rr[i]
  pred <- ns$.emc_alloc_predict(lapply(units, function(z) { z$r_min <- NULL; z }), sol, u$joint)
  rho <- unlist(pred[if (is.null(movable)) seq_len(n_u) else which(movable)])
  # The particle-step floor: no subject coordinate slower than its step alone.
  broken <- 0L
  for (i in sidx) {
    z <- units[[i]]
    if (!is.null(z$p) && any(pred[[i]] > z$rho * z$p + 1e-9)) {
      if (floor) return(Inf)
      broken <- broken + 1L
    }
  }
  cost <- u$C_rest
  for (k in seq_len(n_g)) if (identical(units[[k]]$kind, "shear")) cost <- cost + x[k] * units[[k]]$cost
  for (i in sidx) {
    z <- units[[i]]
    cost <- cost + (if (is.null(z$p)) 0 else qq[i] * z$c_p) + (if (is.null(z$h)) 0 else rr[i] * z$c_h)
  }
  if (xc > 0 && !is.null(u$joint)) cost <- cost + xc * u$joint$cost +
    (if (is.null(u$joint$fixed)) 0 else u$joint$fixed * min(1, xc))
  sf <- if (is.null(u$shear_fixed)) 0 else u$shear_fixed
  hf <- if (is.null(u$hmc_fixed)) 0 else u$hmc_fixed
  if (sf > 0) cost <- cost + sf * (1 - prod(1 - pmin(1, x[vapply(units[seq_len(n_g)], function(z) identical(z$kind, "shear"), logical(1))])))
  if (hf > 0 && n_u > n_g) cost <- cost + hf * (1 - prod(1 - pmin(1, rr[sidx])))
  structure(ns$.emc_alloc_tau(max(rho)) * cost, floor = broken)
}

# Grid oracle over the group units' common rounds x, the collapsed rate and
# each subject's (q, r), vectorised over the subjects' rate combinations.
ws1_oracle <- function(u, xs = seq(0, 4, by = 0.125), xcs = 0, qs = seq(1, 0, by = -0.125),
                       rs = seq(0, 1, by = 0.025)) {
  ns <- ws1_ns()
  units <- u$units
  n_g <- length(u$g_floor)
  sidx <- (n_g + 1L):length(units)
  n_s <- length(sidx)
  qr <- expand.grid(q = qs, r = rs)
  phi <- function(q, r, P, H) {
    b <- pmax(0, q + r - 1)
    (1 - q - r + b) + (q - b) * P + (r - b) * H + b * P * H
  }
  # Per subject and (q, r): its slowest coordinate, feasibility, cost and
  # deviation retention (the credit's numerator).
  per <- lapply(sidx, function(i) {
    z <- units[[i]]
    P <- if (is.null(z$p)) rep(1, length(z$rho)) else z$p
    H <- if (is.null(z$h)) rep(1, length(z$rho)) else z$h
    q <- if (is.null(z$p)) rep(1, nrow(qr)) else qr$q
    r <- if (is.null(z$h)) rep(0, nrow(qr)) else qr$r
    rho <- sapply(seq_along(z$rho), function(k) z$rho[k] * phi(q, r, P[k], H[k]))
    rho <- matrix(rho, nrow(qr))
    ok <- if (is.null(z$p)) rep(TRUE, nrow(qr)) else
      apply(sweep(rho, 2, z$rho * z$p + 1e-9, "<="), 1, all)
    dev <- 0
    if (!is.null(z$d_rho)) {
      Pd <- if (is.null(z$p_dev)) rep(1, length(z$d_rho)) else z$p_dev
      Hd <- if (is.null(z$h_dev)) rep(1, length(z$d_rho)) else z$h_dev
      q0 <- if (is.null(z$q_old)) 1 else z$q_old
      qd <- if (is.null(z$p_dev)) rep(q0, nrow(qr)) else q
      rd <- if (is.null(z$h_dev)) rep(0, nrow(qr)) else r
      base <- z$d_rho / pmax(1 - q0 * (1 - Pd), 1e-6)
      dev <- rowSums(sapply(seq_along(base), function(k) base[k] * phi(qd, rd, Pd[k], Hd[k])))
    }
    list(rho = apply(rho, 1, max), ok = ok, dev = dev, r = r, q = q,
         cost = q * (if (is.null(z$p)) 0 else z$c_p) + r * (if (is.null(z$h)) 0 else z$c_h))
  })
  combos <- as.matrix(do.call(expand.grid, rep(list(seq_len(nrow(qr))), n_s)))
  pick <- function(f) sapply(seq_len(n_s), function(k) per[[k]][[f]][combos[, k]])
  s_rho <- apply(matrix(pick("rho"), nrow(combos)), 1, max)
  ok <- apply(matrix(pick("ok"), nrow(combos)), 1, all)
  s_cost <- rowSums(matrix(pick("cost"), nrow(combos)))
  R_all <- matrix(pick("r"), nrow(combos))
  hf <- if (is.null(u$hmc_fixed)) 0 else u$hmc_fixed
  s_cost <- s_cost + hf * (1 - apply(matrix(1 - pmin(1, R_all), nrow(combos)), 1, prod))
  num <- rowSums(matrix(pick("dev"), nrow(combos)))
  q_all <- matrix(pick("q"), nrow(combos))
  den <- {
    q0 <- rep(1, length(units)); r0 <- numeric(length(units))
    for (i in sidx) { q0[i] <- if (is.null(units[[i]]$q_old)) 1 else units[[i]]$q_old
                      r0[i] <- if (is.null(units[[i]]$r_old)) 0 else units[[i]]$r_old }
    sum(vapply(sidx, function(i) {
      z <- units[[i]]
      if (is.null(z$d_rho)) return(0)
      Pd <- if (is.null(z$p_dev)) rep(1, length(z$d_rho)) else z$p_dev
      Hd <- if (is.null(z$h_dev)) rep(1, length(z$d_rho)) else z$h_dev
      base <- z$d_rho / pmax(1 - q0[i] * (1 - Pd), 1e-6)
      sum(base * phi(q0[i], if (is.null(z$h_dev)) 0 else r0[i], Pd, Hd))
    }, numeric(1)))
  }
  R <- if (den > 0) num / den else rep(1, nrow(combos))
  g_rho <- sapply(seq_len(n_g), function(k) {
    rho <- units[[k]]$rho; f <- u$g_floor[k]; s <- rho - f
    ifelse(R <= 1 | rho >= 1, pmin(1, f + s * R), 1 - (1 - rho) * exp(-s * (R - 1) / (1 - rho)))
  })
  g_rho <- matrix(g_rho, nrow(combos))
  best <- list(obj = Inf)
  for (xc in xcs) {
    fc <- if (xc > 0) ns$.emc_alloc_factor(u$joint$g, xc, if (is.null(u$joint$a)) u$joint$g else u$joint$a) else rep(1, n_g)
    fc <- rep_len(fc, n_g)
    for (x in xs) {
      gr <- sweep(g_rho, 2, fc * vapply(seq_len(n_g), function(k) {
        z <- units[[k]]
        if (identical(z$kind, "shear")) ns$.emc_alloc_factor(z$g, x, if (is.null(z$a)) z$g else z$a) else 1
      }, numeric(1)), "*")
      rho <- pmax(apply(gr, 1, max), s_rho)
      cost <- u$C_rest + s_cost +
        sum(vapply(seq_len(n_g), function(k) if (identical(units[[k]]$kind, "shear")) x * units[[k]]$cost else 0, numeric(1))) +
        (if (xc > 0) xc * u$joint$cost + (if (is.null(u$joint$fixed)) 0 else u$joint$fixed * min(1, xc)) else 0) +
        (if (!is.null(u$shear_fixed) && u$shear_fixed > 0) u$shear_fixed * (1 - (1 - min(1, x))^n_g) else 0)
      obj <- ns$.emc_alloc_tau(rho) * cost
      obj[!ok] <- Inf
      k <- which.min(obj)
      if (obj[k] < best$obj) best <- list(obj = obj[k], x = rep(x, n_g), xc = xc,
                                          q = q_all[k, ], r = R_all[k, ])
    }
  }
  best
}
