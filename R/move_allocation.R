# Allocation of the optional group moves: shear rounds on group coordinates
# (R/shear_move.R) and collapsed moves of the whole group state
# (R/collapsed_move.R). The subject kernel is not allocated: after warm-up an
# eligible subject takes an HMC transition (R/subject_hmc.R) with probability
# hmc_share and its particle step otherwise (fit()'s `hmc_share`), and both
# are part of the iteration's base cost.
#
# A move is worth its cost only where it raises the effective sample size of
# the fit's slowest group coordinate. Every iteration each chain records how
# far the rest of the iteration moved each group coordinate, (mu, log sd) at
# the end of the previous iteration -> just before the shear, together with
# the cost of the rest of the iteration, of each shear round and collapsed
# move, as wall time on the chain's own share of the cores
# (.emc_reserved_cost), so work that leaves the chain's workers idle is
# charged for them; cores taken over from chains that finished their block
# are free. At a proposal rebuild each coordinate's lag-one correlation under
# the rest of the iteration is rho = 1 - ESJD / (2 var), where var is the
# coordinate's variance over every chain's window draws: a chain held in one
# region makes jumps small against that pooled spread, and shows as slow even
# though its own draws look settled. In the sample stage the chain is
# stationary and the kernels settled, so var is taken over every chain's
# production draws, and group jumps (the rest of the iteration's, each shear
# coordinate's, the collapsed move's) pool over the stage's blocks: a short
# window underestimates a slow coordinate's variance, and one slow spell in
# it makes a coordinate look like the slowest. rho is the point estimate: a
# bound below it makes every coordinate, and every plan's prediction, faster
# than measured; a small predicted gain buys nothing (.EMC_ALLOC_MIN_GAIN).
# Group coordinates are mu along the shear's mean directions and the log
# group sds.
#
# A shear round on group coordinate k multiplies its rho by g_k, measured
# from the rounds' selection-weighted squared moves. Shear rounds and
# collapsed moves map the subjects at fixed whitened residuals, so the moves
# within one iteration explore the group state's conditional given those
# residuals, which only the subject updates refresh. As an autoregression
# with correlation a on a share c of the coordinate's variance, k of them
# multiply rho by 1 - c (1 - a^k), with g = 1 - c (1 - a) for one; a comes
# from the net squared jump of each iteration's moves against their per-move
# jumps (.emc_alloc_saturation). c = 1 (a = g) is geometric composition,
# g^k, which is also used until iterations with several moves identify a. A
# measurement belongs to the kernel that made it, so outside production,
# where every rebuild changes the kernels, only the last block measures a
# move. A shear coordinate without a measurement counts at the pooled factor
# and share of the measured rounds of its kind (mean directions or log sds;
# either kind, then a prior factor, without any), and the collapsed move at
# the prior factor until its movement is measured. This gives an unmeasured
# move one allocation window to establish its performance.
# With tau = (1 + rho) / (1 - rho), the allocation minimises
#   max over coordinates of tau * (C_rest + sum_k x_k c_round)
# over expected rounds x_k (.emc_alloc_solve). Coordinates no move can lower
# keep their rest-of-iteration rho and do not bound the target: the
# allocation is over the coordinates the moves can reach. A coordinate a move
# can lower only so far is a floor at that level. A collapsed move moves
# every group coordinate at once: x_c expected moves multiply coordinate k's
# rho by its factor at x_c, so the solver takes the best of a grid over x_c,
# each with the per-unit allocation on the reduced rhos and x_c c_coll added
# to C_rest. Its factors are measured at the current rate, so x_c may grow to
# at most one move per iteration or twice that rate per rebuild (or the
# shared round cap), whichever is larger. The allocation changes only at
# rebuilds, as every adapted quantity does.

.EMC_ALLOC_MAX_ROUNDS <- 4
.EMC_ALLOC_MIN_PAIRS <- 10L
.EMC_ALLOC_SHEAR_PRIOR_G <- 0.5
# Rounds of evidence that weigh as much as a coordinate's reference factor
# (its kind's pooled factor, or the prior): the mean squared jump of n rounds
# has relative standard error about sqrt(2 / n), so twenty rounds pin a
# factor to roughly 30% of the movement it measures.
.EMC_ALLOC_ROUND_EVIDENCE <- 20
# Production evidence of shear rounds and collapsed moves carries across
# rebuilds discounted by this factor, so it is bounded and a move that
# stopped running is measured afresh.
.EMC_ALLOC_EVIDENCE_DECAY <- 0.5
# Move families whose cost splits into the moves' own time and a fixed cost
# per iteration that runs any of them (the frames, the context): the same
# accounting for each.
.EMC_ALLOC_FAMILIES <- c("shear", "collapsed")
.EMC_ALLOC_MIN_GAIN <- 0.9

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

# Every cost the allocation compares is wall time on the cores the chain holds:
# a move's computation in the chain's own process and its likelihood rounds'
# round trips keep those cores as busy as they leave them.
.emc_reserved_cost <- function(elapsed, cores) {
  max(1, as.numeric(cores)) * max(0, as.numeric(elapsed))
}

.emc_alloc_wanted <- function(sampler) {
  .emc_shear_active(sampler) || .emc_collapsed_active(sampler) ||
    isTRUE(.emc_hmc_support(sampler)$ok)
}

.emc_alloc_lsd <- function(V) {
  d <- diag(if (is.matrix(V)) V else as.matrix(V))
  d[d < .Machine$double.xmin] <- .Machine$double.xmin
  0.5 * log(d)
}

# .emc_alloc_lsd of draws k of a p x p x T variance array, as a p x length(k)
# matrix.
.emc_alloc_lsd_draws <- function(V, k) {
  out <- if (getOption("emc2.native_alloc_lsd", TRUE)) emc_alloc_lsd_draws_cpp(V, k)
  if (!is.null(out))
    return(matrix(out, nrow(out), dimnames = list(dimnames(V)[[1L]], NULL)))
  p <- dim(V)[1L]
  # A plain vector index: an index matrix with as many columns as V has
  # dimensions would be read as array subscripts.
  idx <- c(outer(seq_len(p) + (seq_len(p) - 1L) * p, (k - 1) * p * p, "+"))
  matrix(0.5 * log(pmax(V[idx], .Machine$double.xmin)), p,
         dimnames = list(dimnames(V)[[1L]], NULL))
}

# Each row's sample variance (stats::var), NA with fewer than two columns.
.emc_alloc_row_var <- function(x) {
  if (ncol(x) < 2L) return(rep(NA_real_, nrow(x)))
  rowSums((x - rowMeans(x))^2) / (ncol(x) - 1L)
}

.emc_alloc_empty_window <- function(p) {
  list(n = 0L, nG = 0L, G = matrix(0, 2L * p, 2L * p), cpu_rest = 0,
       fam = .emc_alloc_fam_empty())
}

.emc_alloc_window <- function(sampler) {
  w <- sampler$alloc$window
  if (is.null(w)) w <- .emc_alloc_empty_window(sampler$n_pars)
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

# Each move family's whole cost, its moves' own cost, the iterations that
# ran any of them and the moves.
.emc_alloc_fam_empty <- function() {
  stats::setNames(lapply(.EMC_ALLOC_FAMILIES, function(f)
    list(cpu = 0, own = 0, iters = 0L, moves = 0L)), .EMC_ALLOC_FAMILIES)
}

# One family's marginal cost per move and fixed cost per iteration with
# moves, from the window; the fixed cost carries over from the last rebuild
# while the window holds too few such iterations.
.emc_alloc_marginal <- function(w, fam, last = NULL) {
  f <- if (is.null(w$fam)) NULL else w$fam[[fam]]
  if (is.null(f)) f <- list(cpu = 0, own = 0, iters = 0L, moves = 0L)
  carried <- if (is.list(last) && isTRUE(is.finite(last$fixed[[fam]]))) last$fixed[[fam]] else 0
  list(c = if (f$moves >= 5L) f$own / f$moves else NA_real_,
       fixed = if (f$iters >= .EMC_ALLOC_MIN_PAIRS) max(0, f$cpu - f$own) / f$iters else carried)
}

# Costs of one iteration: the rest of the iteration (the subject kernel
# included), and each move family's whole cost, its moves' own cost and
# their number (`moves`, a named list of c(cpu, own, n)).
.emc_alloc_record_cost <- function(sampler, rest, moves = list()) {
  w <- .emc_alloc_window(sampler)
  w$n <- w$n + 1L
  w$cpu_rest <- w$cpu_rest + rest
  if (is.null(w$fam)) w$fam <- .emc_alloc_fam_empty()
  for (fam in intersect(names(moves), .EMC_ALLOC_FAMILIES)) {
    m <- moves[[fam]]
    if (is.null(m) || length(m) < 3L || !is.finite(m[[3]]) || m[[3]] <= 0) next
    f <- w$fam[[fam]]
    f$cpu <- f$cpu + as.numeric(m[[1]])
    f$own <- f$own + as.numeric(m[[2]])
    f$iters <- f$iters + 1L
    f$moves <- f$moves + as.integer(round(m[[3]]))
    w$fam[[fam]] <- f
  }
  sampler$alloc$window <- w
  sampler
}

# rho from a mean squared jump and the coordinate's variance.
.emc_alloc_rho <- function(esjd, var) {
  out <- if (getOption("emc2.native_alloc_rho", TRUE)) emc_alloc_rho_cpp(esjd, var)
  if (!is.null(out)) return(out)
  pmin(1, pmax(0, 1 - esjd / (2 * pmax(var, .Machine$double.xmin))))
}

# The share of each group coordinate's window variance `var` that the group's
# draw given the subjects redraws every iteration: d' S d / N for a mean
# coordinate (row d of Dinv, S the group covariance) and 1 / (2 N) for a log
# group sd. A white part decorrelates at lag one, so it can be no more than
# the measured lag-one decorrelation 1 - rho; it is held to
# .EMC_ALLOC_WHITE_CAP of it, against the window variance's noise.
.EMC_ALLOC_WHITE_CAP <- 0.75
.emc_alloc_white_share <- function(S, Dinv, n, var, rho) {
  p <- nrow(Dinv)
  cv <- c(if (is.matrix(S) && all(dim(S) == p) && all(is.finite(S)))
            rowSums((Dinv %*% S) * Dinv) / n else rep(NA_real_, p),
          rep(1 / (2 * n), p))
  w <- cv / var
  w[!is.finite(w) | w < 0] <- 0
  pmin(w, .EMC_ALLOC_WHITE_CAP * pmax(0, 1 - rho))
}

# The group coordinates as one linear system, for pricing them jointly. A
# coordinate priced alone is an AR(1) series; where the posterior's slow
# direction lies across several coordinates, each of them carries part of it
# and decays slower than its own lag-one correlation says, and a round on
# one of them moves it only within its conditional given the others. With C
# the coordinates' (slow-part) covariance and G the rest of the iteration's
# jump covariance on them, the rest of the iteration acts on linear
# functions as A_rest = I - G C^-1 / 2, x expected rounds on coordinate k as
#   M_k = I - theta_k e_k Q[k, ] / Q[k, k],  Q = C^-1,
# theta_k = 1 - (the rounds' factor within the conditional), and collapsed
# moves as their per-coordinate factors. An iteration is
#   A = M_q ... M_1 diag(f_c) A_rest,
# and coordinate k's integrated autocorrelation time is
#   tau_k = 2 [(I - A)^-1 C]_kk / C_kk - 1.
# With C and G diagonal this is each coordinate's own (1 + rho) / (1 - rho).
#
# cov: the coordinates' within-chain covariance (its diagonal set to `var`);
# G: the rest of the iteration's mean jump outer product; S, Dinv, n: for
# the white part, the group's draw given the subjects (mean coordinates
# Dinv S Dinv' / n, log sds 1 / (2 n)), held to .EMC_ALLOC_WHITE_CAP of the
# jump in every direction. Off-diagonal correlations are shrunk by
# .EMC_ALLOC_COUPLING_SHRINK, and fully where the m jumps behind them are
# fewer than the coordinates. NULL where the system is not usable (the
# coordinates are then priced one by one).
.EMC_ALLOC_COUPLING_SHRINK <- 0.4
.emc_alloc_coupling <- function(cov, var, G, S, Dinv, n, m) {
  q <- length(var)
  p <- nrow(Dinv)
  ok <- function(M) is.matrix(M) && all(dim(M) == q) && all(is.finite(M))
  if (!ok(cov) || !ok(G) || q != 2L * p || !all(is.finite(var)) || any(var <= 0) ||
      any(diag(cov) <= 0) || any(diag(G) <= 0) || !is.finite(m) || m < 2)
    return(NULL)
  lam <- min(1, max(.EMC_ALLOC_COUPLING_SHRINK, q / m))
  off <- function(M) {
    D <- diag(diag(M), q)
    (1 - lam) * M + lam * D
  }
  s <- sqrt(var / diag(cov))
  C <- off(cov * outer(s, s))
  G <- off((G + t(G)) / 2)
  W <- matrix(0, q, q)
  if (is.matrix(S) && all(dim(S) == p) && all(is.finite(S)))
    W[seq_len(p), seq_len(p)] <- Dinv %*% S %*% t(Dinv) / n
  diag(W)[p + seq_len(p)] <- 1 / (2 * n)
  W <- off((W + t(W)) / 2)
  out <- tryCatch({
    R <- chol(G)
    Ri <- backsolve(R, diag(q))
    ev <- eigen(crossprod(Ri, 2 * W %*% Ri), symmetric = TRUE)
    Wc <- crossprod(R, ev$vectors %*% (pmin(pmax(ev$values, 0), .EMC_ALLOC_WHITE_CAP) *
                                          t(ev$vectors)) %*% R) / 2
    Cs <- C - Wc
    Cs <- (Cs + t(Cs)) / 2
    Q <- chol2inv(chol(Cs))
    Ar <- diag(q) - (G - 2 * Wc) %*% Q / 2
    if (!all(is.finite(Ar)) || !all(is.finite(Q)) ||
        max(Mod(eigen(Ar, only.values = TRUE)$values)) >= 1 - 1e-9) NULL else
      list(C = Cs, Q = Q, Ar = Ar, shrink = lam)
  }, error = function(e) NULL)
  out
}

# Each group coordinate's equivalent lag-one rho, (tau - 1) / (tau + 1),
# under the plans in the columns of X (expected rounds per coordinate) with
# the collapsed move's factors fac; g, a: each coordinate's round factor and
# saturation within its conditional (NA: no rounds).
.emc_alloc_group_rho <- function(coupling, g, a, X, fac = NULL) {
  q <- nrow(coupling$C)
  X <- matrix(as.numeric(X), q)
  emc_alloc_group_rho_cpp(coupling$C, coupling$Q, coupling$Ar,
                          if (is.null(fac)) rep(1, q) else rep_len(as.numeric(fac), q),
                          as.numeric(g), as.numeric(a), X)
}

.emc_alloc_tau <- function(rho) {
  rho <- pmin(rho, 1 - 1e-12)
  (1 + rho) / (1 - rho)
}

.emc_alloc_rho_floor <- function(rho, b_share = NULL) {
  if (is.null(b_share)) return(rho)
  b <- rep_len(b_share, length(rho))
  b + (1 - b) * rho
}

.emc_alloc_unit_rho <- function(u, rho = u$rho) {
  .emc_alloc_rho_floor(rho, u$b_share)
}

# Factor on a coordinate's rho from x expected moves per iteration (floor(x)
# moves, and one more with probability x - floor(x)), given the factor g of
# one move and the moves' correlation a within their conditional (see the
# header): k moves have factor 1 - c (1 - a^k) with c = (1 - g) / (1 - a),
# and a = g composes them geometrically, g^k.
.emc_alloc_factor <- function(g, x, a = g) {
  out <- emc_alloc_factor_cpp(g, x, a)
  if (!is.null(out)) return(out)
  m <- floor(x)
  f <- x - m
  (1 - f) * .emc_alloc_factor_k(g, a, m) + f * .emc_alloc_factor_k(g, a, m + 1)
}

.emc_alloc_factor_k <- function(g, a, k) {
  n <- max(length(g), length(a), length(k))
  g <- rep_len(g, n)
  a <- rep_len(a, n)
  k <- rep_len(k, n)
  out <- g^k
  sat <- a < g - 1e-12
  out[sat] <- 1 - (1 - g[sat]) / (1 - a[sat]) * (1 - a[sat]^k[sat])
  out
}

# Least expected moves x with factor(g, x, a) <= t, for targets t < 1 (Inf
# beyond the moves' reach).
.emc_alloc_factor_inverse <- function(g, a, t) {
  out <- rep(Inf, length(t))
  if (g >= 1 - 1e-9) return(out)
  geo <- a >= g - 1e-12
  if (geo) {
    ok <- t > 0
    m <- ceiling(log(t[ok]) / log(g) - 1e-9)
  } else {
    c <- (1 - g) / (1 - a)
    s <- 1 - (1 - t) / c
    ok <- s > 0
    m <- if (a <= 0) rep(1, sum(ok)) else ceiling(log(s[ok]) / log(a) - 1e-9)
  }
  m <- pmax(1, m)
  f0 <- .emc_alloc_factor_k(g, a, m - 1)
  f1 <- .emc_alloc_factor_k(g, a, m)
  out[ok] <- (m - 1) + pmin(1, pmax(0, (f0 - t[ok]) / (f0 - f1)))
  out
}

# Iteration counts by the number of moves k >= 1 they ran (column k), one
# row per coordinate, and the sum of two such tables.
.emc_alloc_khist_add <- function(K, rows, k) {
  if (!length(rows) || k < 1L) return(K)
  if (ncol(K) < k) K <- cbind(K, matrix(0L, nrow(K), k - ncol(K)))
  K[rows, k] <- K[rows, k] + 1L
  K
}

# Iterations in the block that just ended: the allocation window is reset
# only after the move plans have closed their blocks.
.emc_alloc_block_iters <- function(chain) {
  n <- chain$alloc$window$n
  if (is.null(n) || !is.finite(n)) 0 else as.numeric(n)
}

.emc_alloc_khist_sum <- function(A, B) {
  if (is.null(A)) return(B)
  if (is.null(B)) return(A)
  w <- max(ncol(A), ncol(B))
  pad <- function(M) cbind(M, matrix(0L, nrow(M), w - ncol(M)))
  pad(A) + pad(B)
}

# Round factors held to what the rounds can refresh: with `share` of the
# coordinate's variance inside the conditional its rounds sample
# (.emc_shear_orbit_shares), no number of rounds takes the factor below
# 1 - share. One round's factor is raised to that, and the saturation set so
# the composition 1 - (1 - g) (1 - a^k) / (1 - a) levels off there. Factors
# without a measured share are left as they are.
.emc_alloc_cap_rounds <- function(g, a, share) {
  ok <- is.finite(share) & is.finite(g)
  s <- pmin(1, pmax(share, 1e-3))
  g2 <- ifelse(ok, pmax(g, 1 - s), g)
  a0 <- ifelse(is.finite(a), a, g2)
  list(g = g2, a = ifelse(ok, pmin(a0, pmax(0, 1 - (1 - g2) / s)), a))
}

# The moves' correlation a in [0, g] within their conditional, per
# coordinate, from D, the mean over iterations of the net squared jump of
# their moves over twice the coordinate's variance, and K, those iterations'
# move counts (one row per coordinate, or one row for all): D = (1 - g)
# sum_k w_k (1 - a^k) / (1 - a), increasing in a. Without iterations of two
# or more moves a is not identified, and the moves compose geometrically
# (a = g).
.emc_alloc_saturation <- function(g, D, K) {
  a <- g
  if (is.null(K) || !length(D)) return(a)
  K <- as.matrix(K)
  out <- if (getOption("emc2.native_alloc_saturation", TRUE) && is.double(g) && is.double(D))
    emc_alloc_saturation_cpp(g, D, K)
  if (!is.null(out)) return(out)
  for (i in seq_along(g)) {
    w <- K[if (nrow(K) == 1L) 1L else i, ]
    if (sum(w) <= 0 || sum(w[-1L]) <= 0 || !is.finite(D[i]) || g[i] >= 1 - 1e-9) next
    w <- w / sum(w)
    k <- seq_along(w)
    h <- function(x) (1 - g[i]) * sum(w * (1 - x^k)) / (1 - x)
    if (D[i] <= h(0)) { a[i] <- 0; next }
    if (D[i] >= h(g[i])) next
    lo <- 0; hi <- g[i]
    for (it in 1:60) {
      mid <- (lo + hi) / 2
      if (h(mid) < D[i]) lo <- mid else hi <- mid
    }
    a[i] <- (lo + hi) / 2
  }
  a
}

# Shear round factors shrunk towards their kind's (mean directions, log
# sds). The movement m = 1 - g of n rounds has relative standard error about
# sqrt(2 / n), so log m carries sampling variance 2 / n; each coordinate's
# log m is shrunk towards its kind's precision-weighted mean by the spread
# tau^2 between the kind's coordinates beyond that noise. A kind of very
# different coordinates (a slow ridge among quick directions) then barely
# shrinks its well-measured ones; tau^2 is at least 2 / ROUND_EVIDENCE, the
# weight a reference had before. With fewer than three measured coordinates
# in a kind, tau^2 comes from all measured coordinates, else its floor. An
# unmeasured coordinate takes its kind's mean (either kind's, else the prior).
# Movements below .EMC_ALLOC_MIN_MOVE count as that floor.
.EMC_ALLOC_MIN_MOVE <- 1e-3
.emc_alloc_shrink_g <- function(g_obs, n, kind) {
  out <- rep(.EMC_ALLOC_SHEAR_PRIOR_G, length(n))
  meas <- n > 0 & is.finite(g_obs)
  if (!any(meas)) return(out)
  l <- log(pmin(1, pmax(.EMC_ALLOC_MIN_MOVE, 1 - g_obs)))
  v <- 2 / pmax(n, 1e-12)
  t2_min <- 2 / .EMC_ALLOC_ROUND_EVIDENCE
  spread <- function(sel) {
    if (sum(sel) < 3L) return(NA_real_)
    w <- 1 / v[sel]
    lbar <- sum(w * l[sel]) / sum(w)
    sum(w * (l[sel] - lbar)^2) / sum(w) - sum(w * v[sel]) / sum(w)
  }
  t2_all <- spread(meas)
  centre <- function(sel, t2) {
    w <- 1 / (t2 + v[sel])
    sum(w * l[sel]) / sum(w)
  }
  for (kk in unique(kind)) {
    sel <- meas & kind == kk
    t2 <- spread(sel)
    if (!is.finite(t2)) t2 <- t2_all
    t2 <- max(t2_min, if (is.finite(t2)) t2 else t2_min)
    ref <- if (any(sel)) centre(sel, t2) else centre(meas, t2)
    lk <- ifelse(sel, (l / v + ref / t2) / (1 / v + 1 / t2), ref)
    m <- exp(lk[kind == kk])
    # At the movement floor a coordinate does not move: exactly one, which
    # the solver treats as a kernel it cannot use.
    out[kind == kk] <- ifelse(m <= .EMC_ALLOC_MIN_MOVE * (1 + 1e-9), 1, 1 - m)
  }
  out
}

# Candidate expected collapsed moves per iteration: eighths up to 4, then
# halves, then the cap itself.
.emc_alloc_joint_grid <- function(cap) {
  x <- seq(0.125, 4, by = 0.125)
  if (cap >= 4.5) x <- c(x, seq(4.5, cap, by = 0.5))
  sort(unique(c(x[x <= cap], cap)))
}

# Choose expected rounds x for the group coordinates.
#   units: list of list(kind = "fixed" | "shear", ...):
#     fixed: rho, a scalar or vector no move changes;
#     shear: rho, g the per-round factor, cost per round, cap (optional).
# Each unit's cost for a target rho* is that of the cheapest allocation that
# brings its rho to rho*; the objective tau(rho*) * (C_rest + sum costs) is
# evaluated on one fine grid of targets between the best reachable and the
# unallocated maximum, plus every unit's own level under every candidate
# collapsed rate. The units the collapsed move does not touch are priced once
# on that grid; only the group units are re-priced per candidate rate.
# movable: which units count towards the target (default: those a move can
# lower), given so that solves of the same problem with and without a kernel
# are judged on the same coordinates. shear_fixed: the shear's fixed cost
# per iteration that runs any round, charged at the plan's probability of
# such an iteration (.emc_alloc_activation); the collapsed move's
# (joint$fixed) at min(1, x_c).
.emc_alloc_solve <- function(units, C_rest, max_rounds = .EMC_ALLOC_MAX_ROUNDS,
                             joint = NULL, movable = NULL, shear_fixed = 0,
                             coupling = NULL) {
  .emc_alloc_solve_now(units, C_rest, max_rounds, joint, movable, shear_fixed, coupling)
}

# The group units' fields for the joint system (coupling: the leading units,
# one per group coordinate), or NULL where the units are not that block.
.emc_alloc_group_fields <- function(units, coupling, joint, max_rounds) {
  if (is.null(coupling)) return(NULL)
  q <- nrow(coupling$C)
  if (length(units) < q || (!is.null(joint) && !identical(as.integer(joint$idx), seq_len(q))))
    return(NULL)
  us <- units[seq_len(q)]
  kind <- vapply(us, function(u) u$kind, character(1))
  if (!all(kind %in% c("shear", "fixed")) || any(lengths(lapply(us, `[[`, "rho")) != 1L))
    return(NULL)
  sh <- kind == "shear"
  num <- function(f, d, only = TRUE) vapply(seq_len(q), function(k)
    if (only && !sh[k] || is.null(us[[k]][[f]])) d else as.numeric(us[[k]][[f]]), numeric(1))
  g <- num("g", NA_real_)
  a <- num("a", NA_real_)
  a[is.na(a)] <- g[is.na(a)]
  cap <- pmin(num("cap", max_rounds), max_rounds)
  cap[!sh] <- 0
  floor <- pmin(num("floor", 0), cap)
  list(idx = seq_len(q), g = g, a = a, floor = floor, cap = cap, cost = num("cost", 0),
       b = num("b_share", 0, only = FALSE))
}

# The group units' unallocated and least reachable rho under the joint system
# with the collapsed move's factors fac.
.emc_alloc_group_parts <- function(coupling, gf, fac) {
  rho <- .emc_alloc_group_rho(coupling, gf$g, gf$a, cbind(gf$floor, gf$cap), fac)
  list(levels = gf$b + (1 - gf$b) * rho[, 1L], best = gf$b + (1 - gf$b) * rho[, 2L])
}

# The group units' allocation at each target (.emc_alloc_units_need's form)
# under the joint system.
.emc_alloc_group_need <- function(coupling, gf, fac, movable, rs) {
  got <- emc_alloc_group_need_cpp(coupling$C, coupling$Q, coupling$Ar, as.numeric(fac),
                                  gf$g, gf$a, gf$floor, gf$cap, gf$cost, gf$b,
                                  as.logical(movable), as.numeric(rs), 16L)
  n_t <- length(rs)
  lapply(seq_along(gf$idx), function(k) {
    x <- got$x[k, ]
    cost <- x * gf$cost[k]
    if (k == 1L) cost[!got$reached] <- Inf
    list(x = x, cost = cost)
  })
}

.emc_alloc_solve_now <- function(units, C_rest, max_rounds, joint, movable,
                                 shear_fixed, coupling = NULL) {
  n_u <- length(units)
  gf <- .emc_alloc_group_fields(units, coupling, joint, max_rounds)
  if (is.null(gf)) coupling <- NULL
  moved_idx <- if (!is.null(gf)) gf$idx else if (is.null(joint)) integer(0) else
    as.integer(joint$idx)
  # The collapsed move's factor on each group unit at xc expected moves.
  fac_of <- function(xc) if (xc == 0 || is.null(joint)) rep(1, length(moved_idx)) else
    rep_len(.emc_alloc_factor(joint$g, xc, if (is.null(joint$a)) joint$g else joint$a),
            length(moved_idx))
  j_floor <- if (is.null(joint) || is.null(joint$floor)) 0 else min(joint$floor, joint$cap)
  xcs <- if (is.null(joint)) 0 else {
    grid <- .emc_alloc_joint_grid(joint$cap)
    unique(c(j_floor, grid[grid > j_floor]))
  }
  cand <- lapply(xcs, function(xc) {
    if (xc == 0 || !is.null(gf)) return(units)
    fac <- .emc_alloc_factor(joint$g, xc, if (is.null(joint$a)) joint$g else joint$a)
    moved <- units
    for (i in seq_along(moved_idx))
      moved[[moved_idx[i]]]$rho <- moved[[moved_idx[i]]]$rho * fac[i]
    moved
  })
  # Only the moved units' levels change with the collapsed rate.
  parts0 <- .emc_alloc_range_parts(units, max_rounds)
  if (is.null(movable)) movable <- .emc_alloc_movable(units, max_rounds, joint, parts0)
  moved_parts <- if (is.null(gf)) .emc_alloc_scaled_parts(units[moved_idx], max_rounds)
  ranges <- lapply(seq_along(cand), function(c) {
    if (is.null(gf) && (xcs[c] == 0 || !length(moved_idx)))
      return(.emc_alloc_range_of(parts0, movable))
    pm <- if (is.null(gf)) moved_parts(cand[[c]][moved_idx]) else
      .emc_alloc_group_parts(coupling, gf, fac_of(xcs[c]))
    pc <- parts0
    pc$levels[moved_idx] <- pm$levels
    pc$best[moved_idx] <- pm$best
    .emc_alloc_range_of(pc, movable)
  })
  r0 <- ranges[[1L]]
  # The units' fields read once; NULL leaves each to its own pass below.
  us <- if (getOption("emc2.native_alloc_summary", TRUE)) emc_alloc_unit_summary_cpp(units)
  # The plan with nothing chosen still runs every move's exploration.
  x0 <- if (!is.null(us)) us$x0 else vapply(units, function(u)
    if (identical(u$kind, "shear") && !is.null(u$floor)) u$floor else 0, numeric(1))
  shear_rows <- if (!is.null(us)) us$shear else .emc_alloc_act_rows(units)
  activation <- function(X)
    .emc_alloc_act_part(X[shear_rows, , drop = FALSE], 0, shear_fixed, ncol(X))
  j_fixed <- function(xc) if (xc > 0 && !is.null(joint$fixed)) joint$fixed * min(1, xc) else 0
  c_f0 <- sum(x0 * (if (!is.null(us)) us$cost else
    vapply(units, function(u) if (is.null(u$cost)) 0 else u$cost, numeric(1)))) +
    (if (j_floor > 0) j_floor * joint$cost else 0) + j_fixed(j_floor) +
    activation(matrix(x0, ncol = 1L))
  r00 <- if (j_floor > 0) ranges[[1L]] else r0
  zero <- list(x = x0, rho = r00$hi, tau = .emc_alloc_tau(r00$hi),
               cost = C_rest + c_f0,
               objective = .emc_alloc_tau(r00$hi) * (C_rest + c_f0),
               joint = j_floor)
  lo <- min(vapply(ranges, `[[`, numeric(1), "lo"))
  best <- zero
  if (r0$hi > lo + 1e-12) {
    tau_grid <- exp(seq(log(.emc_alloc_tau(lo)), log(.emc_alloc_tau(r0$hi)),
                        length.out = 400L))
    lv <- unlist(lapply(ranges, function(r) r$levels[r$levels >= r$lo & r$levels <= r$hi]))
    targets <- sort(unique(c((tau_grid - 1) / (tau_grid + 1), lv)))
    tau_t <- .emc_alloc_tau(targets)
    still <- setdiff(seq_len(n_u), moved_idx)
    need_of <- function(us, movable, layout = NULL)
      .emc_alloc_units_need(us, movable, targets, max_rounds, layout)
    # Every candidate's moved units differ from the units only in rho.
    moved_layout <- .emc_alloc_need_layout(units[moved_idx], max_rounds)
    need <- vector("list", n_u)
    need[still] <- need_of(units[still], r0$movable[still])
    base <- C_rest + Reduce(`+`, lapply(need[still], `[[`, "cost"), numeric(length(targets)))
    # The shear's activation term is the same at every collapsed rate unless
    # the collapsed move shifts one of its units.
    act_part <- function(nd) {
      if (shear_fixed <= 0 || !length(shear_rows)) return(NULL)
      shear_fixed * emc_alloc_act_term_list_cpp(lapply(nd[shear_rows], `[[`, "x"),
                                                numeric(length(shear_rows)))
    }
    act_held <- if (any(shear_rows %in% moved_idx)) NULL else list(act_part(need))
    best_k <- NULL
    for (c in seq_along(xcs)) {
      nm <- if (is.null(gf))
        need_of(cand[[c]][moved_idx], ranges[[c]]$movable[moved_idx], moved_layout) else
        .emc_alloc_group_need(coupling, gf, fac_of(xcs[c]), ranges[[c]]$movable[moved_idx],
                              targets)
      need_c <- need
      need_c[moved_idx] <- nm
      act <- if (is.null(act_held)) act_part(need_c) else act_held[[1L]]
      c_joint <- xcs[c] * (if (xcs[c] > 0) joint$cost else 0)
      c_fixed <- j_fixed(xcs[c])
      at <- if (getOption("emc2.native_alloc_candidate", TRUE) && is.double(c_joint) &&
                length(c_joint) == 1L && is.double(c_fixed) && length(c_fixed) == 1L)
        emc_alloc_candidate_cpp(base, c_joint, c_fixed, nm, act, tau_t)
      if (is.null(at)) {
        cost <- base + c_joint + c_fixed +
          Reduce(`+`, lapply(nm, `[[`, "cost"), numeric(length(targets))) +
          (if (is.null(act)) 0 else act)
        obj <- tau_t * cost
        k <- which.min(obj)
        at <- if (length(k)) list(k = k, cost = cost[k], objective = obj[k]) else list(k = 0L)
      }
      k <- at$k
      if (k < 1L || !is.finite(at$objective) || at$objective >= best$objective * (1 - 1e-9)) next
      best_nm <- nm
      best_k <- k
      best <- list(x = NULL, rho = targets[k], tau = tau_t[k], cost = at$cost,
                   objective = at$objective, joint = xcs[c])
    }
    # The winning plan's allocation, read once.
    if (!is.null(best_k)) {
      need[moved_idx] <- best_nm
      best$x <- vapply(need, function(z) z$x[best_k], numeric(1))
    }
  }
  # The lag-one surrogate is coarse, so a small predicted gain is not worth
  # its moves.
  if (is.null(best$x) || best$objective > .EMC_ALLOC_MIN_GAIN * zero$objective)
    best <- zero
  best$objective0 <- zero$objective
  best$tau0 <- zero$tau
  # The slowest coordinate of all, moved or not.
  best$tau0_all <- .emc_alloc_tau(r0$hi_all)
  best$pred <- .emc_alloc_predict(units, best, joint, max_rounds)
  if (!is.null(gf)) {
    x <- if (length(best$x) >= length(gf$idx)) best$x[gf$idx] else gf$floor
    x[!is.finite(x)] <- 0
    rho <- .emc_alloc_group_rho(coupling, gf$g, gf$a, x,
                                fac_of(if (is.null(best$joint)) 0 else best$joint))
    best$pred[gf$idx] <- as.list(gf$b + (1 - gf$b) * rho[, 1L])
  }
  best
}

# Each unit's lag-one rho under a solved plan (sol$x, sol$joint): the
# collapsed move's factor at sol$joint on the units it moves, then the
# unit's own rounds.
.emc_alloc_predict <- function(units, sol, joint = NULL, max_rounds = .EMC_ALLOC_MAX_ROUNDS) {
  xc <- if (is.null(sol$joint)) 0 else sol$joint
  idx <- if (is.null(joint)) integer(0) else as.integer(joint$idx)
  fac_c <- if (xc > 0 && length(idx))
    rep_len(.emc_alloc_factor(joint$g, xc, if (is.null(joint$a)) joint$g else joint$a),
            length(idx))
  out <- if (getOption("emc2.native_alloc_predict", TRUE) && is.list(units) &&
             (is.null(sol$x) || is.double(sol$x)) && (is.null(fac_c) || is.double(fac_c)))
    emc_alloc_predict_cpp(units, sol$x, idx, fac_c)
  if (!is.null(out)) return(out)
  lapply(seq_along(units), function(i) {
    u <- units[[i]]
    rho <- u$rho
    b <- rep_len(if (is.null(u$b_share)) 0 else u$b_share, length(rho))
    k <- match(i, idx)
    if (!is.null(fac_c) && !is.na(k)) rho <- rho * fac_c[k]
    x <- if (length(sol$x) >= i && is.finite(sol$x[i])) sol$x[i] else 0
    if (identical(u$kind, "shear"))
      b + (1 - b) * rho * .emc_alloc_factor(u$g, x, if (is.null(u$a)) u$g else u$a) else
      b + (1 - b) * rho
  })
}

# The shear's fixed cost a plan pays per iteration, by target: X holds each
# unit's allocation (rows) at each target (columns). Rounds run whole parts
# plus a Bernoulli fraction independently per coordinate, so an iteration
# builds frames with probability 1 - prod_k (1 - x_k), or 1 once any
# coordinate has a whole round.
.emc_alloc_activation <- function(units, X, shear_fixed) {
  .emc_alloc_act_part(X[.emc_alloc_act_rows(units), , drop = FALSE], 0, shear_fixed, ncol(X))
}

# The shear units, whose rounds activate its fixed cost.
.emc_alloc_act_rows <- function(units)
  which(vapply(units, function(u) identical(u$kind, "shear"), logical(1)))

# A family's fixed cost at each of n_t targets: Xf holds the family's units
# (rows), each running at min(x + add, 1).
.emc_alloc_act_part <- function(Xf, add, fixed, n_t) {
  if (fixed <= 0 || !nrow(Xf)) return(numeric(n_t))
  fixed * emc_alloc_act_term_cpp(Xf, rep_len(as.numeric(add), nrow(Xf)))
}

# The unallocated rho of each unit, the least rho* it can reach, and whether
# any move lowers it at all; the range of targets is over the units that
# can be moved.
.emc_alloc_range <- function(units, max_rounds, movable = NULL)
  .emc_alloc_range_of(.emc_alloc_range_parts(units, max_rounds), movable)

# Each unit's unallocated rho (levels) and least reachable rho (best).
.emc_alloc_range_parts <- function(units, max_rounds) {
  out <- if (getOption("emc2.native_alloc_range", TRUE) && is.list(units))
    emc_alloc_range_parts_cpp(units, max_rounds)
  if (!is.null(out)) return(out)
  .emc_alloc_range_parts_r(units, max_rounds)
}

.emc_alloc_range_parts_r <- function(units, max_rounds) list(
  levels = vapply(units, function(u) max(.emc_alloc_level(u)), numeric(1)),
  best = vapply(units, function(u) switch(u$kind,
    fixed = max(.emc_alloc_unit_rho(u)),
    shear = max(.emc_alloc_rho_floor(
      u$rho * .emc_alloc_factor(u$g, if (is.null(u$cap)) max_rounds else u$cap,
                                if (is.null(u$a)) u$g else u$a), u$b_share))),
    numeric(1)))

# .emc_alloc_range_parts for units that differ from `units` only in a scalar
# rho, as a function of those units: a shear coordinate's factors at its
# floor and cap do not depend on rho. Falls back to the full computation
# for any other units.
.emc_alloc_scaled_parts <- function(units, max_rounds) {
  r_form <- NULL
  function(us) {
    out <- if (getOption("emc2.native_alloc_range", TRUE) && is.list(us))
      emc_alloc_range_parts_cpp(us, max_rounds)
    if (!is.null(out)) return(out)
    if (is.null(r_form)) r_form <<- .emc_alloc_scaled_parts_r(units, max_rounds)
    r_form(us)
  }
}

.emc_alloc_scaled_parts_r <- function(units, max_rounds) {
  kind <- vapply(units, function(u) u$kind, character(1))
  if (!all(kind %in% c("shear", "fixed")) || any(lengths(lapply(units, `[[`, "rho")) != 1L))
    return(function(us) .emc_alloc_range_parts(us, max_rounds))
  sh <- kind == "shear"
  a_of <- function(u) if (is.null(u$a)) u$g else u$a
  f_floor <- vapply(units, function(u) if (identical(u$kind, "shear") && !is.null(u$floor) &&
    u$floor > 0) .emc_alloc_factor(u$g, u$floor, a_of(u)) else NA_real_, numeric(1))
  f_cap <- vapply(units, function(u) if (identical(u$kind, "shear"))
    .emc_alloc_factor(u$g, if (is.null(u$cap)) max_rounds else u$cap, a_of(u)) else NA_real_,
    numeric(1))
  floored <- !is.na(f_floor)
  function(us) {
    rho <- as.numeric(unlist(lapply(us, `[[`, "rho")))
    b <- as.numeric(unlist(lapply(us, function(u)
      if (is.null(u$b_share)) rep(0, length(u$rho)) else u$b_share)))
    levels <- .emc_alloc_rho_floor(rho, b)
    levels[floored] <- .emc_alloc_rho_floor(rho[floored] * f_floor[floored], b[floored])
    best <- .emc_alloc_rho_floor(rho, b)
    best[sh] <- .emc_alloc_rho_floor(rho[sh] * f_cap[sh], b[sh])
    list(levels = levels, best = best)
  }
}

.emc_alloc_range_of <- function(parts, movable = NULL) {
  levels <- parts$levels
  best_rho <- parts$best
  if (is.null(movable) || length(movable) != length(levels))
    movable <- best_rho < levels - 1e-9
  list(levels = levels, best = best_rho, movable = movable,
       hi = if (any(movable)) max(levels[movable], 0) else 0,
       lo = if (any(movable)) max(best_rho[movable], 0) else 0,
       hi_all = max(levels, 0))
}

# Which units a move can lower: by their own moves, or by the joint move.
.emc_alloc_movable <- function(units, max_rounds, joint = NULL, parts = NULL) {
  if (is.null(parts)) parts <- .emc_alloc_range_parts(units, max_rounds)
  movable <- .emc_alloc_range_of(parts)$movable
  if (!is.null(joint)) {
    idx <- as.integer(joint$idx)
    g <- rep_len(joint$g, length(idx))
    movable[idx] <- movable[idx] | (is.finite(g) & g < 1 - 1e-9)
  }
  movable
}

# An unmovable unit's allocation at every target: its exploration floor.
.emc_alloc_unit_still <- function(u, n_t) {
  fl <- if (identical(u$kind, "shear") && !is.null(u$floor)) u$floor else 0
  list(x = rep(fl, n_t), cost = rep(fl * (if (is.null(u$cost)) 0 else u$cost), n_t))
}

# A unit's rho with nothing allocated: a shear coordinate at its exploration
# floor.
.emc_alloc_level <- function(u) {
  if (identical(u$kind, "shear") && !is.null(u$floor) && u$floor > 0)
    return(.emc_alloc_rho_floor(
      u$rho * .emc_alloc_factor(u$g, u$floor, if (is.null(u$a)) u$g else u$a),
      u$b_share))
  .emc_alloc_unit_rho(u)
}

# A unit's allocation for each target rho* in rs: list(x, cost), with cost
# Inf where the target is out of its reach.
.emc_alloc_unit_need <- function(u, rs, max_rounds) {
  n_t <- length(rs)
  if (identical(u$kind, "fixed"))
    return(list(x = numeric(n_t),
                cost = ifelse(max(.emc_alloc_unit_rho(u)) > rs, Inf, 0)))
  if (identical(u$kind, "shear")) {
    r <- u$rho
    b <- rep_len(if (is.null(u$b_share)) 0 else u$b_share, length(r))
    base <- .emc_alloc_rho_floor(r, b)
    g <- u$g
    cap <- if (is.null(u$cap)) max_rounds else u$cap
    floor <- if (is.null(u$floor)) 0 else u$floor
    x <- rep(floor, n_t)
    over <- base > rs
    if (any(over)) {
      den <- (1 - b) * r
      target <- if (den > 0) (rs[over] - b) / den else rep(-Inf, sum(over))
      xo <- .emc_alloc_factor_inverse(g, if (is.null(u$a)) g else u$a, target)
      xo <- pmax(xo, floor)
      xo[!is.finite(xo) | xo > cap + 1e-9] <- Inf
      x[over] <- xo
    }
    return(list(x = x, cost = ifelse(is.finite(x), x * u$cost, Inf)))
  }
  stop("unknown unit kind: ", u$kind)
}

# The units' kinds and their shear coordinates' fixed fields (all but rho),
# for pricing them together; units that differ only in rho share it.
.emc_alloc_need_layout <- function(us, max_rounds) {
  kind <- vapply(us, function(u) u$kind, character(1))
  shear <- kind == "shear" & lengths(lapply(us, `[[`, "rho")) == 1L
  sh <- us[shear]
  field <- function(f, d) vapply(sh, function(u) if (is.null(u[[f]])) d else
    as.numeric(u[[f]]), numeric(1))
  list(kind = kind, shear = shear, g = field("g", NA_real_),
       a = vapply(sh, function(u) as.numeric(if (is.null(u$a)) u$g else u$a), numeric(1)),
       cost = field("cost", NA_real_), cap = field("cap", max_rounds),
       floor = field("floor", 0),
       b_share = field("b_share", 0))
}

# Every unit's allocation at each target (.emc_alloc_unit_need), with the
# movable shear units priced together in one call. layout:
# .emc_alloc_need_layout(us), when the caller prices units of the same layout
# repeatedly.
.emc_alloc_units_need <- function(us, movable, rs, max_rounds, layout = NULL) {
  if (is.null(layout)) layout <- .emc_alloc_need_layout(us, max_rounds)
  out <- vector("list", length(us))
  shear <- movable & layout$shear
  fixed <- movable & layout$kind == "fixed"
  if (any(fixed))
    out[fixed] <- emc_alloc_fixed_need_cpp(vapply(us[fixed], function(u)
      max(.emc_alloc_unit_rho(u)), numeric(1)), as.numeric(rs))
  for (i in which(!shear & !fixed)) out[[i]] <- if (movable[i])
    .emc_alloc_unit_need(us[[i]], rs = rs, max_rounds = max_rounds)
    else .emc_alloc_unit_still(us[[i]], length(rs))
  sh <- which(shear)
  if (length(sh)) {
    k <- match(sh, which(layout$shear))
    out[sh] <- emc_alloc_shears_need_cpp(as.numeric(unlist(lapply(us[sh], `[[`, "rho"))),
                                         layout$b_share[k],
                                         layout$g[k], layout$a[k], layout$cost[k],
                                         layout$cap[k], layout$floor[k], as.numeric(rs))
  }
  out
}

# The last m stored draws of a chain (at least two).
.emc_alloc_draws <- function(chain, m) {
  i <- chain$samples$idx
  seq.int(max(1L, i - max(2L, m) + 1L), i)
}

# Pooled sample covariance of the columns of x (rows are coordinates).
.emc_alloc_cov <- function(x) {
  if (ncol(x) < 2L) return(matrix(NA_real_, nrow(x), nrow(x)))
  xc <- x - rowMeans(x)
  tcrossprod(xc) / (ncol(x) - 1L)
}

.emc_alloc_within_share <- function(chains) {
  out <- if (getOption("emc2.native_alloc_within", TRUE) && is.list(chains))
    emc_alloc_within_share_cpp(chains)
  if (!is.null(out)) {
    rn <- rownames(chains[[1L]])
    if (!is.null(rn)) {
      names(out$total) <- rn
      if (length(rn) > 1L) names(out$var) <- rn
    }
    return(out)
  }
  p <- nrow(chains[[1L]])
  all_draws <- do.call(cbind, chains)
  total <- .emc_alloc_row_var(all_draws)
  within_by_chain <- vapply(chains, .emc_alloc_row_var, numeric(p))
  if (is.null(dim(within_by_chain))) within_by_chain <- matrix(within_by_chain, p)
  within <- rowMeans(within_by_chain, na.rm = TRUE)
  within[!is.finite(within)] <- NA_real_
  between <- pmax(0, total - within)
  b_share <- ifelse(is.finite(total) & is.finite(within) & total == 0, 0, between / total)
  b_share[!is.finite(b_share)] <- NA_real_
  b_share <- pmin(1, pmax(0, b_share))
  list(var = within, total = total, b_share = b_share)
}

# The draws' statistics every chain's measurement shares: mu and the log
# group sds by chain, the mean group variance, and the log sds' within-chain
# variances. Each chain then needs only its own directions. `cache` (an
# environment made for one rebuild, whose draws do not change) holds them per
# set of draw indices, so the rebuild computes them once.
.emc_alloc_draw_stats <- function(emc, vdraws, cache = NULL) {
  st <- if (!is.null(cache) && identical(cache$vdraws, vdraws)) cache$stats else NULL
  if (is.null(st)) {
    mu_by_chain <- lapply(seq_along(emc), function(i)
      emc[[i]]$samples$theta_mu[, vdraws[[i]], drop = FALSE])
    lsd_by_chain <- lapply(seq_along(emc), function(i)
      .emc_alloc_lsd_draws(emc[[i]]$samples$theta_var, vdraws[[i]]))
    S_mean <- tryCatch(Reduce(`+`, lapply(seq_along(emc), function(i)
      apply(emc[[i]]$samples$theta_var[, , vdraws[[i]], drop = FALSE], 1:2, mean))) / length(emc),
      error = function(e) NULL)
    st <- list(mu_by_chain = mu_by_chain, S_mean = S_mean, lsd_by_chain = lsd_by_chain,
               lsd_var = .emc_alloc_within_share(lsd_by_chain)$var)
  }
  if (!is.null(cache)) {
    cache$vdraws <- vdraws
    cache$stats <- st
  }
  st
}

# A chain's measurements at a rebuild, against every chain's window draws.
.emc_alloc_measure <- function(emc, j, stage = NULL, cache = NULL) {
  ch <- emc[[j]]
  p <- ch$n_pars
  w <- ch$alloc$window
  stored <- is.null(w) || w$n < .EMC_ALLOC_MIN_PAIRS
  m <- if (stored) 200L else w$n
  draws <- lapply(emc, .emc_alloc_draws, m = m)
  # Production: variances over every production draw, group jumps pooled
  # over the stage's blocks (see the header).
  production <- identical(stage, "sample")
  pool <- NULL
  if (production && identical(ch$alloc$stage, "sample")) {
    pool <- ch$alloc$pool
    if (!stored) {
      if (is.list(pool)) {
        w$G <- w$G + pool$G
        w$nG <- w$nG + pool$nG
      }
      pool <- list(G = w$G, nG = w$nG)
    }
  }
  vdraws <- draws
  if (production) {
    prod <- lapply(emc, function(x) {
      st <- x$samples$stage[seq_len(x$samples$idx)]
      which(st == "sample")
    })
    if (all(lengths(prod) >= lengths(draws))) vdraws <- prod
  }
  if (stored) {
    # Whole-iteration jumps between consecutive stored draws stand in for a
    # window not yet recorded.
    k <- draws[[j]]
    w <- .emc_alloc_empty_window(p)
    if (length(k) >= 2L) {
      s <- ch$samples
      for (t in k[-1L]) {
        d <- c(s$theta_mu[, t], .emc_alloc_lsd(s$theta_var[, , t])) -
          c(s$theta_mu[, t - 1L], .emc_alloc_lsd(s$theta_var[, , t - 1L]))
        if (all(is.finite(d))) { w$G <- w$G + tcrossprod(d); w$nG <- w$nG + 1L }
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
  ds <- .emc_alloc_draw_stats(emc, vdraws, cache)
  mu_parts <- .emc_alloc_within_share(lapply(ds$mu_by_chain, function(x) Dinv %*% x))
  lsd_parts <- .emc_alloc_within_share(ds$lsd_by_chain)
  g_var <- c(mu_parts$var, lsd_parts$var)
  g_b_share <- c(mu_parts$b_share, lsd_parts$b_share)
  g_esjd <- if (w$nG > 0) rowSums((T_g %*% w$G) * T_g) / w$nG else rep(0, 2L * p)
  g_rho <- if (w$nG > 0) .emc_alloc_rho(g_esjd, g_var) else rep(NA_real_, 2L * p)
  g_valid <- is.finite(g_rho) & is.finite(g_var) & g_var > 0 & is.finite(g_b_share)
  g_rho[!g_valid] <- 1
  g_var[!is.finite(g_var) | g_var <= 0] <- 1
  g_b_share[!is.finite(g_b_share)] <- 0
  # Each iteration redraws the group given the subjects, so a group
  # coordinate is a white part (that draw's variance) on a slow part the
  # moves carry. Its units are priced on the slow part's coefficient, and
  # the moves' jumps are measured against the slow part's variance.
  g_white <- .emc_alloc_white_share(ds$S_mean, Dinv, ch$n_subjects, g_var, g_rho)
  g_white[!g_valid] <- 0
  g_rho <- pmin(1, g_rho / (1 - g_white))
  g_var_slow <- (1 - g_white) * g_var
  g_fac <- g_sat <- rep(NA_real_, 2L * p)
  g_measured <- logical(2L * p)
  g_share <- rep(NA_real_, 2L * p)
  coupling <- NULL
  if (shear) {
    last <- .emc_shear_record(ch$shear$plan$last, p)
    # The log-sd coordinates are the same in every chain, so in production
    # their rounds pool across chains (mean directions are each chain's own).
    if (production && identical(ch$shear$plan$stage, "sample")) {
      ks <- p + seq_len(p)
      for (i in seq_along(emc)) {
        if (i == j) next
        o <- emc[[i]]$shear$plan
        if (!is.list(o) || !identical(o$stage, "sample") || !is.list(o$last)) next
        ol <- .emc_shear_record(o$last, p)
        last$E[ks] <- last$E[ks] + ol$E[ks]
        last$n[ks] <- last$n[ks] + ol$n[ks]
        last$moved[ks] <- last$moved[ks] + ol$moved[ks]
        last$N[ks] <- last$N[ks] + ol$N[ks]
        oK <- ol$K
        oK[seq_len(p), ] <- 0
        last$K <- .emc_alloc_khist_sum(last$K, oK)
      }
    }
    # A coordinate's factor is its own rounds' shrunk towards its kind's by
    # their noise and the kind's spread (.emc_alloc_shrink_g), so a handful
    # of rounds shade the reference rather than replace it.
    nk <- last$n
    last_n <- nk
    last_iters <- last$iters
    g_measured <- nk > 0
    g_obs <- ifelse(g_measured,
                    pmin(1, pmax(0, 1 - (last$E / pmax(nk, 1e-12)) / (2 * g_var_slow))), NA_real_)
    g_fac <- .emc_alloc_shrink_g(g_obs, nk, rep(1:2, each = p))
    g_sat <- .emc_alloc_saturation(
      g_fac, last$N / pmax(1, rowSums(last$K)) / (2 * g_var_slow), last$K)
    # The variance each coordinate's rounds can refresh (plan$cond_var); the
    # log sds' pooled over chains as their rounds are.
    cond <- ch$shear$plan$cond_var
    if (length(cond) != 2L * p) cond <- rep(NA_real_, 2L * p)
    if (production && identical(ch$shear$plan$stage, "sample")) {
      all_cond <- vapply(emc, function(e) {
        cv <- e$shear$plan$cond_var
        if (length(cv) == 2L * p) cv[p + seq_len(p)] else rep(NA_real_, p)
      }, numeric(p))
      pooled <- rowMeans(matrix(all_cond, p), na.rm = TRUE)
      cond[p + seq_len(p)] <- ifelse(is.finite(pooled), pooled, NA_real_)
    }
    g_share <- pmin(1, cond / g_var_slow)
    capped <- .emc_alloc_cap_rounds(g_fac, g_sat, g_share)
    g_fac <- capped$g
    g_sat <- capped$a
    # The group block priced jointly: its coordinates' covariance and the
    # rest of the iteration's operator on them, and each coordinate's rounds
    # as moves within its conditional given the others.
    if (all(g_valid) && w$nG > 0) {
      chains_y <- lapply(seq_along(ds$mu_by_chain), function(i)
        rbind(Dinv %*% ds$mu_by_chain[[i]], ds$lsd_by_chain[[i]]))
      coupling <- .emc_alloc_coupling(
        Reduce(`+`, lapply(chains_y, .emc_alloc_cov)) / length(chains_y), g_var,
        T_g %*% (w$G / w$nG) %*% t(T_g), ds$S_mean, Dinv, ch$n_subjects,
        min(w$nG, sum(lengths(vdraws))))
      if (!is.null(coupling)) {
        cv <- 1 / diag(coupling$Q)
        c_obs <- ifelse(g_measured,
                        pmin(1, pmax(0, 1 - (last$E / pmax(nk, 1e-12)) / (2 * cv))), NA_real_)
        coupling$g <- .emc_alloc_shrink_g(c_obs, nk, rep(1:2, each = p))
        coupling$a <- .emc_alloc_saturation(
          coupling$g, last$N / pmax(1, rowSums(last$K)) / (2 * cv), last$K)
        capped <- .emc_alloc_cap_rounds(coupling$g, coupling$a, pmin(1, cond / cv))
        coupling$g <- capped$g
        coupling$a <- capped$a
      }
    }
  }
  # The collapsed move moves every group coordinate at once: its factors come
  # from its acceptance-weighted jumps in the same coordinates.
  coll <- .emc_collapsed_active(ch) && is.list(ch$collapsed$plan) &&
    is.matrix(ch$collapsed$plan$C)
  c_fac <- c_sat <- rep(.EMC_ALLOC_SHEAR_PRIOR_G, 2L * p)
  c_measured <- FALSE
  if (coll) {
    last <- ch$collapsed$plan$last
    c_measured <- is.list(last) && last$n > 0
    if (c_measured) {
      # A coordinate's moves weigh against the move's pooled factor over
      # all its coordinates, as a shear coordinate's rounds do against its
      # kind's: once measured, a move that does not move earns no credit.
      last <- .emc_collapsed_record(last, p)
      c_obs <- pmin(1, pmax(0, 1 - (rowSums((T_g %*% last$E) * T_g) / last$n) /
                              (2 * g_var_slow)))
      n0 <- .EMC_ALLOC_ROUND_EVIDENCE
      c_fac <- (last$n * c_obs + n0 * mean(c_obs)) / (last$n + n0)
      c_sat <- .emc_alloc_saturation(
        c_fac, rowSums((T_g %*% last$N) * T_g) / max(1, sum(last$K)) / (2 * g_var_slow),
        last$K)
    }
  }
  c_fac[!g_valid] <- 1
  c_sat[!g_valid] <- 1
  # The base cost per iteration: the rest of the iteration, the subject
  # kernel included.
  C_rest <- if (!stored && w$n > 0 && w$cpu_rest > 0) w$cpu_rest / w$n else NA_real_
  if (!is.finite(C_rest)) {
    carried <- ch$alloc$last$C_rest
    if (isTRUE(is.finite(carried))) C_rest <- carried
  }
  # Each move family's marginal and fixed costs, carried over where the
  # window holds too few iterations of it.
  marg <- stats::setNames(lapply(.EMC_ALLOC_FAMILIES, function(f)
    .emc_alloc_marginal(if (stored) NULL else w, f, ch$alloc$last)), .EMC_ALLOC_FAMILIES)
  list(stored = stored, pairs = w$nG, pool = pool,
       shear = shear,
       g_rho = g_rho, g_var = g_var, g_white = g_white, g_b_share = g_b_share, g_valid = g_valid,
       g_fac = g_fac, g_sat = g_sat, g_measured = g_measured, g_share = g_share,
       coupling = coupling,
       g_curved = if (shear && is.list(ch$shear[["curve"]])) ch$shear$curve_rec$k,
       g_n = if (shear) last_n else numeric(2L * p),
       g_iters = if (shear) last_iters else numeric(2L * p),
       c_n = if (coll && c_measured) last$n else 0,
       c_iters = if (coll && c_measured) last$iters else 0,
       coll = coll, c_fac = c_fac, c_sat = c_sat, c_measured = c_measured,
       c_coll = if (stored) NA_real_ else marg$collapsed$c,
       coll_fixed = marg$collapsed$fixed,
       c_rate = if (!stored && w$n > 0 && !is.null(w$fam)) w$fam$collapsed$moves / w$n else 0,
       C_rest = C_rest,
       c_round = if (stored) NA_real_ else marg$shear$c,
       shear_fixed = marg$shear$fixed,
       fixed = list(shear = marg$shear$fixed, collapsed = marg$collapsed$fixed))
}

# The units and costs the solver sees for one chain, one unit per group
# coordinate. B: transitions until the next rebuild. A move with fewer than
# .EMC_ALLOC_ROUND_EVIDENCE moves behind it runs at the exploration rate that
# collects the rest by then, and a measured move runs at up to one move per
# iteration, or twice the rate its evidence was gathered at (its moves over
# the iterations they span, both discounted alike) where that is higher:
# beyond one per iteration the moves share iterations, whose saturation the
# evidence gathered at a lower rate cannot identify. A cap tied to the rate
# alone would ratchet a move down for good after one plan that ran little of
# it. A shear coordinate with no rounds behind it counts at its kind's pooled
# factor and may also run up to once per iteration: held at its exploration
# rate, an unmeasured slow coordinate would set the minimax target, and every
# other coordinate would be bought only down to it.
.emc_alloc_units <- function(ms, shear_on, B = NULL) {
  bounds <- function(n_ev, iters, unmeasured_cap = 0) {
    if (is.null(B) || !is.finite(B) || B <= 0) return(list(floor = 0, cap = Inf))
    n0 <- .EMC_ALLOC_ROUND_EVIDENCE
    floor <- max(0, n0 - n_ev) / B
    span <- if (length(iters) && is.finite(iters) && iters > 0) iters else B
    list(floor = floor, cap = max(floor, 2 * n_ev / span, if (n_ev > 0) 1 else unmeasured_cap))
  }
  has_C_rest <- length(ms$C_rest) == 1L && is.finite(ms$C_rest)
  C_rest <- if (has_C_rest) ms$C_rest else 1
  # Unmeasured costs are guessed against the whole iteration's.
  c_round <- if (is.finite(ms$c_round)) ms$c_round else 0.1 * C_rest
  coupled <- shear_on && is.list(ms$coupling) && length(ms$coupling$g) == length(ms$g_rho) &&
    all(is.finite(ms$coupling$g[is.finite(ms$g_fac)]))
  units <- lapply(seq_along(ms$g_rho), function(k)
    if (shear_on && is.finite(ms$g_fac[k])) {
      b <- bounds(if (is.null(ms$g_n)) 0 else ms$g_n[k],
                  if (is.null(ms$g_iters)) 0 else ms$g_iters[k], unmeasured_cap = 1)
      # A curved coordinate's rounds are priced on lag one, which credits the
      # rest of the iteration with movement the subject updates undo, and
      # cannot see what they do to the coordinates correlated with it; with
      # no rounds bought its factor is never measured. Its curve runs a round
      # in every iteration, and is dropped when that does not pay in realised
      # autocorrelation and cost (.emc_shear_curve_pays).
      if (identical(ms$g_curved, k)) {
        b$floor <- max(b$floor, 1)
        b$cap <- max(b$cap, 1)
      }
      # Under the joint system a round's factor is the one within the
      # coordinate's conditional.
      g_k <- if (coupled) ms$coupling$g[k] else ms$g_fac[k]
      a_k <- if (coupled) ms$coupling$a[k] else if (is.null(ms$g_sat)) g_k else ms$g_sat[k]
      list(kind = "shear", rho = ms$g_rho[k],
           b_share = if (is.null(ms$g_b_share)) NULL else ms$g_b_share[k],
           g = g_k, a = a_k, cost = c_round,
           floor = min(b$floor, .EMC_ALLOC_MAX_ROUNDS), cap = min(b$cap, .EMC_ALLOC_MAX_ROUNDS))
    } else list(kind = "fixed", rho = ms$g_rho[k],
                b_share = if (is.null(ms$g_b_share)) NULL else ms$g_b_share[k], cost = 0))
  joint <- if (isTRUE(ms$coll)) {
    b <- bounds(if (is.null(ms$c_n)) 0 else ms$c_n, if (is.null(ms$c_iters)) 0 else ms$c_iters)
    list(idx = seq_along(ms$g_rho), g = ms$c_fac,
         a = if (is.null(ms$c_sat)) ms$c_fac else ms$c_sat,
         cost = if (is.finite(ms$c_coll)) ms$c_coll else 0.1 * C_rest,
         fixed = if (isTRUE(is.finite(ms$coll_fixed))) ms$coll_fixed else 0,
         floor = b$floor,
         cap = min(b$cap, max(.EMC_ALLOC_MAX_ROUNDS, 2 * ms$c_rate)))
  }
  list(units = units, C_rest = C_rest, c_round = c_round, joint = joint,
       coupling = if (coupled) ms$coupling[c("C", "Q", "Ar")],
       shear_fixed = if (isTRUE(is.finite(ms$shear_fixed))) ms$shear_fixed else 0)
}

# The allocation with the shear's fixed cost (the frames an iteration with
# rounds builds): the best plan with rounds against the best without, which
# it must beat by the minimum gain. The fixed cost is charged at the
# probability that an iteration runs any round under the plan
# (.emc_alloc_solve).
.emc_alloc_solve_units <- function(u) {
  fixed_shear <- if (is.null(u$shear_fixed)) 0 else u$shear_fixed
  if (fixed_shear > 0 && any(vapply(u$units, function(x) identical(x$kind, "shear"), logical(1)))) {
    # Both sides are judged on the coordinates the full kernel can move.
    movable <- .emc_alloc_movable(u$units, .EMC_ALLOC_MAX_ROUNDS, u$joint)
    # Exploration rounds run either way: a floored coordinate keeps its floor
    # (and nothing more) on the side without shear, which then pays the
    # frames too.
    no_shear <- u
    explored <- FALSE
    no_shear$units <- lapply(u$units, function(x) {
      if (!identical(x$kind, "shear")) return(x)
      if (!is.null(x$floor) && x$floor > 0) {
        explored <<- TRUE
        x$cap <- x$floor
        return(x)
      }
      list(kind = "fixed", rho = x$rho, b_share = x$b_share, cost = 0)
    })
    no_shear$shear_fixed <- 0
    no_shear$shear_act <- if (explored) fixed_shear else 0
    no_shear$movable <- movable
    without <- .emc_alloc_solve_units(no_shear)
    on <- u
    on$shear_fixed <- 0
    on$shear_act <- fixed_shear
    on$movable <- movable
    with <- .emc_alloc_solve_units(on)
    if (isTRUE(with$objective < without$objective &&
               with$objective <= .EMC_ALLOC_MIN_GAIN * without$objective0)) {
      with$objective0 <- without$objective0
      with$tau0 <- without$tau0
      return(with)
    }
    return(without)
  }
  movable <- if (!is.null(u$movable)) u$movable else
    .emc_alloc_movable(u$units, .EMC_ALLOC_MAX_ROUNDS, u$joint)
  .emc_alloc_solve(u$units, u$C_rest, joint = u$joint, movable = movable,
                   shear_fixed = if (is.null(u$shear_act)) 0 else u$shear_act,
                   coupling = u$coupling)
}

# Rebuild every chain's shear rounds, collapsed rate and subject kernel (in
# add_proposals(), after the likelihood fits are checked and the shear plans
# built). B is the number of transitions until the next rebuild; hmc_share
# is fit()'s.
.emc_alloc_rebuild <- function(emc, stage, n_cores, B, hmc_share = 0.5) {
  if (identical(stage, "preburn") || !.emc_alloc_wanted(emc[[1L]])) return(emc)
  emc <- .emc_prepare_chain_rng(emc)
  nch <- length(emc)
  Vbar <- lapply(emc, function(ch) {
    k <- .emc_alloc_draws(ch, if (is.null(ch$alloc$window)) 200L else
      max(2L, ch$alloc$window$n))
    apply(ch$samples$theta_var[, , k, drop = FALSE], 1:2, mean)
  })
  hst <- lapply(seq_len(nch), function(j) .emc_hmc_plan_chain(emc[[j]], stage, Vbar[[j]]))
  # A bracketed step size for every eligible subject that has none, on a
  # private copy of its state.
  jobs <- NULL
  if (hmc_share > 0) for (j in seq_len(nch)) {
    if (!isTRUE(hst[[j]]$supported)) next
    for (s in seq_len(emc[[j]]$n_subjects)) {
      sj <- hst[[j]]$subjects[[s]]
      if (isTRUE(sj$eligible) && !isTRUE(sj$warmup) && is.na(sj$eps))
        jobs <- rbind(jobs, data.frame(chain = j, subject = s))
    }
  }
  if (!is.null(jobs)) {
    res <- .emc_with_preserved_rng(auto_mclapply(seq_len(nrow(jobs)), function(k) {
      j <- jobs$chain[k]
      .emc_quiet_fork(.emc_hmc_calibrate(emc[[j]], hst[[j]], jobs$subject[k]))
    }, mc.cores = min(nrow(jobs), max(1L, n_cores)), mc.set.seed = FALSE))
    for (k in seq_len(nrow(jobs))) {
      j <- jobs$chain[k]
      s <- jobs$subject[k]
      r <- res[[k]]
      if (!is.list(r) || is.null(r$seed)) next
      .emc_reject_add(r$rejects, "pilot")
      emc[[j]]$rng$subjects[[s]] <- r$seed
      hst[[j]] <- .emc_hmc_apply_calibration(hst[[j]], s, r)
    }
  }
  # Every chain is measured before any chain's plan changes.
  stats_cache <- new.env()
  ms <- lapply(seq_len(nch), function(j) .emc_alloc_measure(emc, j, stage, cache = stats_cache))
  for (j in seq_len(nch)) {
    ch <- emc[[j]]
    m <- ms[[j]]
    u <- .emc_alloc_units(m, m$shear, B = B)
    sol <- .emc_alloc_solve_units(u)
    n_g <- length(m$g_rho)
    if (m$shear) ch$shear$plan$rounds <- sol$x[seq_len(n_g)]
    if (isTRUE(m$coll)) ch$collapsed$plan$rounds <- sol$joint
    kern <- .emc_hmc_set_rates(hst[[j]], ch$n_subjects, hmc_share)
    st <- kern$st
    rec <- list(
      stage = stage, B = B, iteration = ch$samples$idx,
      transitions = if (is.null(ch$sample_transitions)) NA_integer_ else
        as.integer(ch$sample_transitions),
      stored = m$stored, pairs = m$pairs, C_rest = u$C_rest, c_round = u$c_round,
      shear_fixed = u$shear_fixed, fixed = m$fixed,
      tau = sol$tau, tau0 = sol$tau0, tau0_all = sol$tau0_all, cost = sol$cost,
      objective = sol$objective, objective0 = sol$objective0,
      # Predicted lag-one rho per coordinate under the plan, and with nothing
      # bought, for checking the allocation against the next block.
      pred_group = if (length(sol$pred) >= n_g) vapply(sol$pred[seq_len(n_g)], max, numeric(1)),
      pred0_group = .emc_alloc_rho_floor(m$g_rho, m$g_b_share),
      group_rho = .emc_alloc_rho_floor(m$g_rho, m$g_b_share),
      group_g = m$g_fac, group_a = m$g_sat,
      # Whether the group coordinates were priced jointly, with the
      # conditional round factors and the off-diagonal shrinkage used.
      group_share = m$g_share,
      coupled = !is.null(u$coupling),
      coupled_g = if (!is.null(u$coupling)) m$coupling$g,
      coupled_shrink = if (!is.null(u$coupling)) m$coupling$shrink,
      rounds = if (m$shear) sol$x[seq_len(n_g)] else numeric(n_g),
      curved_rounds = if (m$shear && length(m$g_curved) == 1L) sol$x[m$g_curved] else 0,
      collapsed_rounds = sol$joint, collapsed_g = if (isTRUE(m$coll)) m$c_fac,
      collapsed_a = if (isTRUE(m$coll)) m$c_sat,
      collapsed_cost = if (is.null(u$joint)) NA_real_ else u$joint$cost,
      collapsed_eps = if (isTRUE(m$coll)) ch$collapsed$plan$eps else NA_real_,
      rate = if (is.null(st$subjects)) numeric(ch$n_subjects) else
        vapply(st$subjects, function(x) x$rate, numeric(1)),
      particle_rate = kern$q,
      eps = if (is.null(st$subjects)) rep(NA_real_, ch$n_subjects) else
        vapply(st$subjects, function(x) x$eps, numeric(1)),
      calib_cpu = if (is.null(st$subjects)) 0 else
        sum(vapply(st$subjects, function(x) x$calib_cpu, numeric(1))))
    ch$alloc <- list(window = .emc_alloc_empty_window(ch$n_pars), stage = stage, pool = m$pool,
                     last = rec, history = c(ch$alloc$history, list(rec)),
                     particle = list(rate = kern$q))
    if (!is.null(st)) {
      st$alloc <- rec
      ch$subject_hmc <- st
    }
    emc[[j]] <- ch
  }
  emc
}
