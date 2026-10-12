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
# together with the cost of the rest of the iteration, of each shear round
# and collapsed move, as wall time on the cores the chain holds
# (.emc_reserved_cost), so work that leaves the chain's workers idle is
# charged for them. Subject moves run as rounds over the workers: a particle
# step or HMC transition costs its worker CPU, and the rest of a round's wall
# time is a fixed cost, paid in the rest of the iteration for the particle
# round (which always runs) and once per iteration for the HMC round by any
# plan that uses HMC. At a proposal rebuild each coordinate's lag-one
# correlation under the rest of the iteration is rho = 1 - ESJD / (2 var),
# where var is the coordinate's variance over every chain's window draws: a
# chain held in one region makes jumps small against that pooled spread, and
# shows as slow even though its own draws look settled. In the sample stage
# the chain is stationary and the kernels settled, so var is taken over every
# chain's production draws, and group jumps (the rest of the iteration's,
# each shear coordinate's, the collapsed move's) pool over the stage's
# blocks: a short window underestimates a slow coordinate's variance, and
# one slow spell in it makes a coordinate look like the slowest. An upper
# confidence bound on ESJD makes rho a lower bound, so noise in a short
# window buys no moves. Group coordinates are mu along the shear's mean
# directions and the log group sds; subject coordinates are alpha_s
# whitened in the reference frame of its likelihood fit.
#
# A shear round on group coordinate k multiplies its rho by g_k, and an HMC
# transition taken with probability r multiplies subject coordinate i's by
# 1 - r + r h_i; g and h come from the moves' own jumps, measured the same
# way (a round's selection-weighted squared move, an HMC transition's
# acceptance-weighted squared jump). Shear rounds and collapsed moves map
# the subjects at fixed whitened residuals, so the moves within one
# iteration explore the group state's conditional given those residuals,
# which only the subject updates refresh. As an autoregression with
# correlation a on a share c of the coordinate's variance, k of them
# multiply rho by 1 - c (1 - a^k), with g = 1 - c (1 - a) for one; a comes
# from the net squared jump of each iteration's moves against their
# per-move jumps (.emc_alloc_saturation). c = 1 (a = g) is geometric
# composition, g^k, which is also used until iterations with several moves
# identify a. A measurement belongs to the kernel that made it, so outside
# production, where every rebuild changes the kernels, only the last block
# measures a move. A shear coordinate without a measurement counts at the
# pooled factor and share of the measured rounds of its kind (mean
# directions or log sds; either kind, then a prior factor, without any), and
# the collapsed move at the prior factor until its movement is measured. This
# gives an unmeasured move one allocation window to establish its performance.
# With tau = (1 + rho) / (1 - rho), the
# allocation minimises
#   max over coordinates of tau * (C_rest + sum_k x_k c_round + sum_s r_s c_s)
# over expected rounds x_k and rates r_s (.emc_alloc_solve). Coordinates no
# move can lower (no shear, no usable fit, a subject whose steps do not move
# it) keep their rest-of-iteration rho and do not bound the target: the
# allocation is over the coordinates the moves can reach. A coordinate a
# move can lower only so far is a floor at that level. A
# collapsed move moves every group coordinate at once: x_c expected moves
# multiply coordinate k's rho by its factor at x_c, so the solver takes the
# best of a grid over x_c, each with the per-unit allocation on the reduced
# rhos and x_c c_coll added to C_rest. Its factors are measured at the
# current rate, so x_c may grow to at most twice that rate per rebuild (or
# the shared round cap, whichever is larger).
#
# Each subject's particle step is a move like the others: it runs with
# probability q_s, multiplies coordinate i's rho by 1 - q_s + q_s p_i and
# costs q_s c_p,s, with p and c_p measured from the steps that ran. The rest
# of the iteration without it has rho t_i = rho_i / (1 - q (1 - p_i)) at the
# rate q in force, so HMC can stand in for the particle step where it moves
# the subject more per unit cost; each subject's (q_s, r_s) is the cheapest
# pair that brings its coordinates to the target, and no coordinate may end
# slower than under its particle step alone. The two kernels are selected
# jointly (R/sampling.R). The group mixes on the subjects' deviations from
# the group mean, so a plan's group coordinates are credited with the
# deviation autocorrelation its (q_s, r_s) imply relative to the measured
# one (.emc_alloc_dev_ratio): HMC regimes, the plan without HMC, and a plan
# whose chosen particle rates move the subjects less than it was credited
# for, which is solved again against keeping every particle step
# (.emc_alloc_settle). A subject the
# allocation would use without measured HMC moves first gets a bracketed
# step size and pilot transitions on a private copy of its state. Under
# emc2.alloc_subject_floor = "refresh" the particle step is not a floor: each
# subject is refreshed once per iteration by the particle step or HMC
# (q_s + r_s >= 1) and only has to meet the target. The
# allocation changes only at rebuilds, as every adapted quantity does.

.EMC_ALLOC_MAX_ROUNDS <- 4
.EMC_ALLOC_MIN_PAIRS <- 10L
.EMC_ALLOC_SHEAR_PRIOR_G <- 0.5
# Rounds of evidence that weigh as much as a coordinate's reference factor
# (its kind's pooled factor, or the prior): the mean squared jump of n rounds
# has relative standard error about sqrt(2 / n), so twenty rounds pin a
# factor to roughly 30% of the movement it measures.
.EMC_ALLOC_ROUND_EVIDENCE <- 20
# Production evidence of every move family (shear rounds, collapsed moves,
# subject HMC) carries across rebuilds discounted by this factor, so it is
# bounded and a move that stopped running is measured afresh.
.EMC_ALLOC_EVIDENCE_DECAY <- 0.5
# Move families whose cost splits into the moves' own time and a fixed cost
# per iteration that runs any of them (the frames, the context, the round's
# latency): the same accounting for each.
.EMC_ALLOC_FAMILIES <- c("shear", "collapsed", "hmc")
.EMC_ALLOC_MIN_GAIN <- 0.9
# Subject-HMC rates tried for every subject at once, for the group's credit
# (emc2.alloc_hmc_regimes switches the search off).
.EMC_ALLOC_HMC_REGIMES <- c(0.5, 1)
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
       D = replicate(n, matrix(0, p, p), simplify = FALSE),
       cpu_rest = 0, fam = .emc_alloc_fam_empty(),
       cpu_hmc = numeric(n), hmc_att = integer(n),
       P = replicate(n, matrix(0, p, p), simplify = FALSE), nP = integer(n),
       cpu_part = numeric(n))
}

# The particle step's own jumps and cost. A subject's cost is its marginal
# one, its worker CPU `times`: the round runs every iteration for the
# ensemble candidates too, so its latency and idle cores stay in the
# rest-of-iteration cost. Without worker times, `cost` (the step's wall time
# on the chain's cores) is shared equally.
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
    w$cpu_part[s] <- w$cpu_part[s] + if (sum(t_ran) > 0) t_ran[k] else cost * share[k]
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
# final ones: their jumps, and the jumps of their deviations from the group
# mean (what the subject kernels alone move).
.emc_alloc_record_subjects <- function(sampler, alpha, mu = NULL) {
  idx <- sampler$samples$idx
  w <- .emc_alloc_window(sampler)
  if (is.null(w$D)) w$D <- replicate(ncol(alpha), matrix(0, nrow(alpha), nrow(alpha)),
                                     simplify = FALSE)
  old <- matrix(sampler$samples$alpha[, , idx], nrow(alpha))
  d <- alpha - old
  if (all(is.finite(d))) {
    dd <- if (is.null(mu)) d else
      d - (as.numeric(mu) - as.numeric(sampler$samples$theta_mu[, idx]))
    for (s in seq_len(ncol(d))) {
      w$A[[s]] <- w$A[[s]] + tcrossprod(d[, s])
      w$D[[s]] <- w$D[[s]] + tcrossprod(dd[, s])
    }
    w$nA <- w$nA + 1L
  }
  sampler$alloc$window <- w
  sampler
}

# An HMC transition's cost is its marginal one, the worker CPU of its task;
# the rest of the round's wall time on the chain's cores (`hmc_round`: its
# latency, the chain process's share and idle cores) is a fixed cost of
# running the round at all.
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

# Costs of one iteration: the rest of the iteration, each move family's
# whole cost, its moves' own cost and their number (`moves`, a named list of
# c(cpu, own, n); HMC's own cost is its tasks' worker CPU), and per subject
# the HMC transitions' cost and count.
.emc_alloc_record_cost <- function(sampler, rest, tasks, moves = list()) {
  w <- .emc_alloc_window(sampler)
  w$n <- w$n + 1L
  w$cpu_rest <- w$cpu_rest + rest
  if (is.null(w$fam)) w$fam <- .emc_alloc_fam_empty()
  work <- 0
  attempts <- 0L
  for (s in seq_along(tasks)) {
    tk <- tasks[[s]]
    if (is.null(tk)) next
    c_s <- if (!is.null(tk$work)) tk$work else if (!is.null(tk$cost)) tk$cost else tk$cpu
    w$cpu_hmc[s] <- w$cpu_hmc[s] + c_s
    w$hmc_att[s] <- w$hmc_att[s] + 1L
    work <- work + c_s
    attempts <- attempts + 1L
  }
  if (attempts > 0L) moves$hmc <- c(cpu = if (is.null(moves$hmc)) work else as.numeric(moves$hmc)[1],
                                    own = work, n = attempts)
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

# Factor on a coordinate's rho from x expected moves per iteration (floor(x)
# moves, and one more with probability x - floor(x)), given the factor g of
# one move and the moves' correlation a within their conditional (see the
# header): k moves have factor 1 - c (1 - a^k) with c = (1 - g) / (1 - a),
# and a = g composes them geometrically, g^k.
.emc_alloc_factor <- function(g, x, a = g) {
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

.emc_alloc_khist_sum <- function(A, B) {
  if (is.null(A)) return(B)
  if (is.null(B)) return(A)
  w <- max(ncol(A), ncol(B))
  pad <- function(M) cbind(M, matrix(0L, nrow(M), w - ncol(M)))
  pad(A) + pad(B)
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
# movable: which units count towards the target (default: those a move can
# lower), given so that solves of the same problem with and without a kernel
# are judged on the same coordinates. shear_fixed and hmc_fixed: each
# family's fixed cost per iteration that runs any of its moves, charged at
# the plan's probability of such an iteration (.emc_alloc_activation); the
# collapsed move's (joint$fixed) at min(1, x_c).
.emc_alloc_solve <- function(units, C_rest, max_rounds = .EMC_ALLOC_MAX_ROUNDS,
                             joint = NULL, movable = NULL, shear_fixed = 0,
                             hmc_fixed = 0) {
  units <- lapply(units, .emc_alloc_as_subject)
  n_u <- length(units)
  refresh <- identical(getOption("emc2.alloc_subject_floor", "particle"), "refresh")
  moved_idx <- if (is.null(joint)) integer(0) else as.integer(joint$idx)
  j_floor <- if (is.null(joint) || is.null(joint$floor)) 0 else min(joint$floor, joint$cap)
  xcs <- if (is.null(joint)) 0 else {
    grid <- .emc_alloc_joint_grid(joint$cap)
    unique(c(j_floor, grid[grid > j_floor]))
  }
  cand <- lapply(xcs, function(xc) {
    if (xc == 0) return(units)
    fac <- .emc_alloc_factor(joint$g, xc, if (is.null(joint$a)) joint$g else joint$a)
    moved <- units
    for (i in seq_along(moved_idx))
      moved[[moved_idx[i]]]$rho <- moved[[moved_idx[i]]]$rho * fac[i]
    moved
  })
  if (is.null(movable)) movable <- .emc_alloc_movable(units, max_rounds, joint)
  ranges <- lapply(cand, .emc_alloc_range, max_rounds = max_rounds, movable = movable)
  r0 <- ranges[[1L]]
  c_p0 <- sum(vapply(units, function(u)
    if (identical(u$kind, "subject") && !is.null(u$p) && !is.null(u$c_p)) u$c_p else 0,
    numeric(1)))
  # The plan with nothing chosen still runs every move's exploration.
  x0 <- vapply(units, function(u)
    if (identical(u$kind, "shear") && !is.null(u$floor)) u$floor else 0, numeric(1))
  r_min <- vapply(units, function(u)
    if (identical(u$kind, "subject") && !is.null(u$h) && !is.null(u$r_min)) u$r_min else 0,
    numeric(1))
  activation <- function(X) .emc_alloc_activation(units, X, r_min, shear_fixed, hmc_fixed)
  j_fixed <- function(xc) if (xc > 0 && !is.null(joint$fixed)) joint$fixed * min(1, xc) else 0
  c_f0 <- sum(x0 * vapply(units, function(u) if (is.null(u$cost)) 0 else u$cost, numeric(1))) +
    (if (j_floor > 0) j_floor * joint$cost else 0) + j_fixed(j_floor) +
    activation(matrix(x0, ncol = 1L))
  r00 <- if (j_floor > 0) ranges[[1L]] else r0
  zero <- list(x = x0, q = rep(1, n_u), rho = r00$hi, tau = .emc_alloc_tau(r00$hi),
               cost = C_rest + c_p0 + c_f0,
               objective = .emc_alloc_tau(r00$hi) * (C_rest + c_p0 + c_f0),
               joint = j_floor)
  lo <- min(vapply(ranges, `[[`, numeric(1), "lo"))
  best <- zero
  if (r0$hi > lo + 1e-12 || c_p0 > 0) {
    tau_grid <- exp(seq(log(.emc_alloc_tau(lo)), log(.emc_alloc_tau(r0$hi)),
                        length.out = 400L))
    lv <- unlist(lapply(ranges, function(r) r$levels[r$levels >= r$lo & r$levels <= r$hi]))
    targets <- sort(unique(c((tau_grid - 1) / (tau_grid + 1), lv)))
    tau_t <- .emc_alloc_tau(targets)
    still <- setdiff(seq_len(n_u), moved_idx)
    need_of <- function(us, movable) lapply(seq_along(us), function(i)
      if (movable[i]) .emc_alloc_unit_need(us[[i]], rs = targets, max_rounds = max_rounds,
                                           refresh = refresh)
      else .emc_alloc_unit_still(us[[i]], length(targets)))
    need <- vector("list", n_u)
    need[still] <- need_of(units[still], r0$movable[still])
    base <- C_rest + Reduce(`+`, lapply(need[still], `[[`, "cost"), numeric(length(targets)))
    for (c in seq_along(xcs)) {
      nm <- need_of(cand[[c]][moved_idx], ranges[[c]]$movable[moved_idx])
      need_c <- need
      need_c[moved_idx] <- nm
      cost <- base + xcs[c] * (if (xcs[c] > 0) joint$cost else 0) + j_fixed(xcs[c]) +
        Reduce(`+`, lapply(nm, `[[`, "cost"), numeric(length(targets))) +
        activation(do.call(rbind, lapply(need_c, `[[`, "x")))
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
  # The slowest coordinate of all, moved or not.
  best$tau0_all <- .emc_alloc_tau(r0$hi_all)
  best
}

# The fixed costs a plan pays per iteration, by target: X holds each unit's
# allocation (rows) at each target (columns). Shear rounds run whole parts
# plus a Bernoulli fraction independently per coordinate, so an iteration
# builds frames with probability 1 - prod_k (1 - x_k), or 1 once any
# coordinate has a whole round; subjects take HMC independently at their
# total rates r_min + x, so a round runs with probability 1 - prod_s (1 - r_s).
.emc_alloc_activation <- function(units, X, r_min, shear_fixed, hmc_fixed) {
  out <- numeric(ncol(X))
  if (shear_fixed > 0) {
    sh <- which(vapply(units, function(u) identical(u$kind, "shear"), logical(1)))
    if (length(sh)) {
      Xs <- pmin(X[sh, , drop = FALSE], 1)
      Xs[!is.finite(Xs)] <- 1
      out <- out + shear_fixed * (1 - apply(1 - Xs, 2, prod))
    }
  }
  if (hmc_fixed > 0) {
    hs <- which(vapply(units, function(u) identical(u$kind, "subject") && !is.null(u$h),
                       logical(1)))
    if (length(hs)) {
      R <- pmin(X[hs, , drop = FALSE] + r_min[hs], 1)
      R[!is.finite(R)] <- 1
      out <- out + hmc_fixed * (1 - apply(1 - R, 2, prod))
    }
  }
  out
}

# The unallocated rho of each unit, the least rho* it can reach, and whether
# any move lowers it at all; the range of targets is over the units that
# can be moved.
.emc_alloc_range <- function(units, max_rounds, movable = NULL) {
  levels <- vapply(units, function(u) max(.emc_alloc_level(u)), numeric(1))
  best_rho <- vapply(units, function(u) switch(u$kind,
    fixed = max(u$rho),
    shear = u$rho * .emc_alloc_factor(u$g, if (is.null(u$cap)) max_rounds else u$cap,
                                      if (is.null(u$a)) u$g else u$a),
    subject = max(if (is.null(u$h)) .emc_alloc_level(u) else
                    .emc_alloc_subject_base(u) * u$h)), numeric(1))
  if (is.null(movable) || length(movable) != length(units))
    movable <- best_rho < levels - 1e-9
  list(levels = levels, best = best_rho, movable = movable,
       hi = if (any(movable)) max(levels[movable], 0) else 0,
       lo = if (any(movable)) max(best_rho[movable], 0) else 0,
       hi_all = max(levels, 0))
}

# Which units a move can lower: by their own moves, or by the joint move.
.emc_alloc_movable <- function(units, max_rounds, joint = NULL) {
  movable <- .emc_alloc_range(units, max_rounds)$movable
  if (!is.null(joint)) {
    idx <- as.integer(joint$idx)
    g <- rep_len(joint$g, length(idx))
    movable[idx] <- movable[idx] | (is.finite(g) & g < 1 - 1e-9)
  }
  movable
}

# An unmovable unit's allocation at every target: nothing bought, its
# particle step (if priced) running.
.emc_alloc_unit_still <- function(u, n_t) {
  c_p <- if (identical(u$kind, "subject") && !is.null(u$p) && !is.null(u$c_p)) u$c_p else 0
  fl <- if (identical(u$kind, "shear") && !is.null(u$floor)) u$floor else 0
  list(x = rep(fl, n_t), q = rep(1, n_t), cost = rep(c_p + fl * (if (is.null(u$cost)) 0 else u$cost), n_t))
}

.emc_alloc_as_subject <- function(u) {
  if (!identical(u$kind, "hmc")) return(u)
  list(kind = "subject", rho = u$rho, p = NULL, c_p = 0, h = u$h, c_h = u$cost)
}

# A subject's rho under its particle step alone (before any HMC).
.emc_alloc_subject_base <- function(u) if (is.null(u$p)) u$rho else u$rho * u$p

# A unit's rho with nothing allocated: subjects run their particle step and
# any HMC rate a regime holds them to.
.emc_alloc_level <- function(u) {
  if (identical(u$kind, "shear") && !is.null(u$floor) && u$floor > 0)
    return(u$rho * .emc_alloc_factor(u$g, u$floor, if (is.null(u$a)) u$g else u$a))
  if (!identical(u$kind, "subject")) return(u$rho)
  base <- .emc_alloc_subject_base(u)
  if (!is.null(u$r_min) && !is.null(u$h)) base * (1 - u$r_min * (1 - u$h)) else base
}

# A unit's allocation for each target rho* in rs: list(x, q, cost), with
# cost Inf where the target is out of its reach. refresh: a subject's
# particle step is not a floor, but it is refreshed once per iteration.
.emc_alloc_unit_need <- function(u, rs, max_rounds, refresh = FALSE) {
  n_t <- length(rs)
  if (identical(u$kind, "fixed"))
    return(list(x = numeric(n_t), q = rep(1, n_t),
                cost = ifelse(max(u$rho) > rs, Inf, 0)))
  if (identical(u$kind, "shear")) {
    r <- u$rho
    g <- u$g
    cap <- if (is.null(u$cap)) max_rounds else u$cap
    floor <- if (is.null(u$floor)) 0 else u$floor
    x <- rep(floor, n_t)
    over <- r > rs
    if (any(over)) {
      xo <- .emc_alloc_factor_inverse(g, if (is.null(u$a)) g else u$a, rs[over] / r)
      xo <- pmax(xo, floor)
      xo[!is.finite(xo) | xo > cap + 1e-9] <- Inf
      x[over] <- xo
    }
    return(list(x = x, q = rep(1, n_t), cost = ifelse(is.finite(x), x * u$cost, Inf)))
  }
  # Subject: for each particle rate q, the least HMC rate r with
  # rho_i (1 - q (1 - p_i)) (1 - r (1 - h_i)) <= min(rho*, level_i) on every
  # coordinate (<= rho* under refresh, with r >= 1 - q); then the cheapest q.
  rho <- u$rho
  # The floor is the subject's particle step alone: a regime's HMC rate is
  # bought for the group, and with it a subject may stand its particle step
  # down where the HMC moves it at least as well.
  level <- .emc_alloc_subject_base(u)
  qs <- if (is.null(u$p) || is.null(u$h)) 1 else if (!is.null(u$q_fixed)) u$q_fixed else
    .EMC_ALLOC_PARTICLE_RATES
  r_min <- if (is.null(u$h) || is.null(u$r_min)) 0 else u$r_min
  r_need <- matrix(pmax(r_min, if (refresh) 1 - qs else 0), length(qs), n_t)
  for (i in seq_along(rho)) {
    K <- if (refresh) rs else pmin(rs, level[i])
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
  # The rate a regime holds the subject to is priced with the regime; x is
  # the rate above it.
  cost <- qs * c_p + (r_need - r_min) * c_h
  cost[!fin] <- Inf
  k <- max.col(t(-cost), ties.method = "first")
  pick <- cbind(k, seq_len(n_t))
  list(x = r_need[pick] - r_min, q = qs[k], cost = cost[pick])
}

# The last m stored draws of a chain (at least two).
.emc_alloc_draws <- function(chain, m) {
  i <- chain$samples$idx
  seq.int(max(1L, i - max(2L, m) + 1L), i)
}

# A chain's measurements at a rebuild, against every chain's window draws.
.emc_alloc_measure <- function(emc, j, hmc_state, Vbar, stage = NULL) {
  ch <- emc[[j]]
  p <- ch$n_pars
  n <- ch$n_subjects
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
    vapply(vdraws[[i]], function(t) c(drop(Dinv %*% s$theta_mu[, t]),
                                     .emc_alloc_lsd(s$theta_var[, , t])), numeric(2L * p))
  }))
  g_var <- apply(gx, 1L, stats::var)
  g_esjd <- if (w$nG > 0) rowSums((T_g %*% w$G) * T_g) / w$nG else rep(0, 2L * p)
  g_rho <- if (w$nG > 0) .emc_alloc_rho_lo(g_esjd, g_var, w$nG) else rep(0, 2L * p)
  g_fac <- g_sat <- rep(NA_real_, 2L * p)
  g_measured <- logical(2L * p)
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
    # A coordinate's factor weighs its own rounds against its kind's pooled
    # factor (or the prior) by the rounds behind each, so a handful of
    # rounds shade the reference rather than replace it (see the header).
    nk <- last$n
    last_n <- nk
    g_measured <- nk > 0
    g_obs <- ifelse(g_measured,
                    pmin(1, pmax(0, 1 - (last$E / pmax(nk, 1e-12)) / (2 * g_var))), NA_real_)
    kind <- rep(1:2, each = p)
    g_ref <- rep(.EMC_ALLOC_SHEAR_PRIOR_G, 2L * p)
    for (kk in 1:2) {
      sel <- g_measured & kind == kk
      if (!any(sel)) sel <- g_measured
      if (any(sel)) g_ref[kind == kk] <- sum(nk[sel] * g_obs[sel]) / sum(nk[sel])
    }
    n0 <- .EMC_ALLOC_ROUND_EVIDENCE
    g_fac <- ifelse(g_measured, (nk * g_obs + n0 * g_ref) / (nk + n0), g_ref)
    g_sat <- .emc_alloc_saturation(
      g_fac, last$N / pmax(1, rowSums(last$K)) / (2 * g_var), last$K)
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
                              (2 * g_var)))
      n0 <- .EMC_ALLOC_ROUND_EVIDENCE
      c_fac <- (last$n * c_obs + n0 * mean(c_obs)) / (last$n + n0)
      c_sat <- .emc_alloc_saturation(
        c_fac, rowSums((T_g %*% last$N) * T_g) / max(1, sum(last$K)) / (2 * g_var),
        last$K)
    }
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
  # Deviations from the group mean: their variance and lag-1 autocorrelation
  # under the subject kernels as run (the particle step and any HMC).
  d_var <- matrix(NA_real_, n, p)
  d_rho <- matrix(NA_real_, n, p)
  has_D <- !stored && !is.null(w$D) && w$nA > 0
  for (s in seq_len(n)) {
    ax <- do.call(cbind, lapply(seq_along(emc), function(i)
      W[[s]] %*% matrix(emc[[i]]$samples$alpha[, s, vdraws[[i]]], p)))
    s_var[s, ] <- apply(ax, 1L, stats::var)
    if (w$nA > 0) s_rho[s, ] <- .emc_alloc_rho_lo(
      rowSums((W[[s]] %*% w$A[[s]]) * W[[s]]) / w$nA, s_var[s, ], w$nA)
    if (has_D) {
      dx <- do.call(cbind, lapply(seq_along(emc), function(i)
        W[[s]] %*% (matrix(emc[[i]]$samples$alpha[, s, vdraws[[i]]], p) -
                      emc[[i]]$samples$theta_mu[, vdraws[[i]], drop = FALSE])))
      d_var[s, ] <- apply(dx, 1L, stats::var)
      d_rho[s, ] <- .emc_alloc_rho_lo(
        rowSums((W[[s]] %*% w$D[[s]]) * W[[s]]) / w$nA, d_var[s, ], w$nA)
    }
  }
  # The group's lag-1 autocorrelation under perfectly re-equilibrated
  # subjects (its floor), from the one measured under the subjects' deviation
  # autocorrelation r_bar: rho = f + (1 - f) r_bar.
  r_bar <- if (has_D && any(is.finite(d_rho))) mean(d_rho[is.finite(d_rho)]) else NA_real_
  # Particle steps: their own factors and cost per step, carried over from
  # the last rebuild while the window holds too few of a subject's steps.
  prev <- ch$alloc$particle
  nP <- if (stored || is.null(w$nP)) integer(n) else w$nP
  p_fac <- matrix(NA_real_, n, p)
  # The particle step's factor on the deviations from the group mean (the
  # group mean is fixed during the step, so its jumps are the deviations'),
  # for the group's credit.
  p_dev <- matrix(NA_real_, n, p)
  c_part <- rep(NA_real_, n)
  for (s in seq_len(n)) {
    if (nP[s] >= .EMC_ALLOC_MIN_PAIRS) {
      esjd <- rowSums((W[[s]] %*% w$P[[s]]) * W[[s]]) / nP[s]
      p_fac[s, ] <- pmin(1, pmax(0, 1 - esjd / (2 * pmax(s_var[s, ], .Machine$double.xmin))))
      if (all(is.finite(d_var[s, ])))
        p_dev[s, ] <- pmin(1, pmax(0, 1 - esjd / (2 * pmax(d_var[s, ], .Machine$double.xmin))))
      c_part[s] <- w$cpu_part[s] / nP[s]
    } else if (is.matrix(prev$fac) && all(dim(prev$fac) == c(n, p)) &&
               all(is.finite(prev$fac[s, ])) && is.finite(prev$cost[s])) {
      p_fac[s, ] <- prev$fac[s, ]
      if (is.matrix(prev$fac_dev) && all(dim(prev$fac_dev) == c(n, p)))
        p_dev[s, ] <- prev$fac_dev[s, ]
      c_part[s] <- prev$cost[s]
    }
  }
  q_old <- if (is.null(prev$rate)) rep(1, n) else prev$rate
  r_old <- ch$alloc$last$rate
  if (length(r_old) != n || !all(is.finite(r_old))) r_old <- numeric(n)
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
  # Each move family's marginal and fixed costs, carried over where the
  # window holds too few iterations of it.
  marg <- stats::setNames(lapply(.EMC_ALLOC_FAMILIES, function(f)
    .emc_alloc_marginal(if (stored) NULL else w, f, ch$alloc$last)), .EMC_ALLOC_FAMILIES)
  list(stored = stored, pairs = c(group = w$nG, subjects = w$nA), pool = pool,
       shear = shear,
       g_rho = g_rho, g_var = g_var, g_fac = g_fac, g_sat = g_sat, g_measured = g_measured,
       g_n = if (shear) last_n else numeric(2L * p),
       c_n = if (coll && c_measured) last$n else 0,
       coll = coll, c_fac = c_fac, c_sat = c_sat, c_measured = c_measured,
       c_coll = if (stored) NA_real_ else marg$collapsed$c,
       coll_fixed = marg$collapsed$fixed,
       c_rate = if (!stored && w$n > 0 && !is.null(w$fam)) w$fam$collapsed$moves / w$n else 0,
       W = W, s_rho = s_rho, t_rho = t_rho, p_fac = p_fac, p_dev = p_dev, c_part = c_part,
       q_old = q_old, r_old = r_old, s_var = s_var, C_rest = C_rest,
       d_var = d_var, d_rho = d_rho, r_bar = r_bar,
       g_floor = if (is.finite(r_bar) && r_bar < 1)
         pmin(g_rho, pmax(0, (g_rho - r_bar) / (1 - r_bar))) else g_rho,
       c_round = if (stored) NA_real_ else marg$shear$c,
       shear_fixed = marg$shear$fixed,
       c_hmc = if (!stored) ifelse(w$hmc_att > 0, w$cpu_hmc / pmax(1L, w$hmc_att), NA_real_)
               else rep(NA_real_, n),
       hmc_fixed = marg$hmc$fixed,
       fixed = list(shear = marg$shear$fixed, collapsed = marg$collapsed$fixed,
                    hmc = marg$hmc$fixed))
}

# HMC factors h for subject s from its evidence; NULL without enough.
.emc_alloc_hmc_factor <- function(ev, W, var, min_n) {
  if (is.null(ev) || ev$n < min_n) return(NULL)
  esjd <- rowSums((W %*% ev$M) * W) / ev$n
  pmin(1, pmax(0, 1 - esjd / (2 * pmax(var, .Machine$double.xmin))))
}

# The units and costs the solver sees for one chain. optimistic: an eligible
# subject without enough HMC evidence counts as h = 0.
# B: transitions until the next rebuild. A move with fewer than
# .EMC_ALLOC_ROUND_EVIDENCE moves behind it runs at the exploration rate that
# collects the rest by then, and no move runs faster than twice the rate its
# evidence was gathered at: trust grows with evidence, and a move measured
# as useless fades with it.
.emc_alloc_units <- function(ms, hmc_state, shear_on, optimistic, B = NULL) {
  n <- nrow(ms$s_rho)
  opts <- .emc_hmc_options()
  bounds <- function(n_ev) {
    if (is.null(B) || !is.finite(B) || B <= 0) return(list(floor = 0, cap = Inf))
    n0 <- .EMC_ALLOC_ROUND_EVIDENCE
    list(floor = max(0, n0 - n_ev) / B, cap = max(max(0, n0 - n_ev), 2 * n_ev) / B)
  }
  C_rest <- if (is.finite(ms$C_rest)) ms$C_rest else 1
  # Unmeasured costs are guessed against the whole iteration's.
  C_iter <- C_rest + sum(ms$c_part[is.finite(ms$c_part)])
  c_round <- if (is.finite(ms$c_round)) ms$c_round else 0.1 * C_iter
  units <- lapply(seq_along(ms$g_rho), function(k)
    if (shear_on && is.finite(ms$g_fac[k])) {
      b <- bounds(if (is.null(ms$g_n)) 0 else ms$g_n[k])
      list(kind = "shear", rho = ms$g_rho[k], g = ms$g_fac[k],
           a = if (is.null(ms$g_sat)) ms$g_fac[k] else ms$g_sat[k], cost = c_round,
           floor = min(b$floor, .EMC_ALLOC_MAX_ROUNDS), cap = min(b$cap, .EMC_ALLOC_MAX_ROUNDS))
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
      else if (!is.null(h)) u$h_n <- ev$n
      if (!is.null(h)) {
        c_s <- if (is.finite(ms$c_hmc[s])) ms$c_hmc[s] else if (!is.null(ev) && ev$n > 0 &&
          ev$cpu > 0) ev$cpu / ev$n else 3 * C_iter / n
        u$h <- h
        u$c_h <- c_s
        # Its factor on the deviations, for the group's credit.
        if (!is.null(ms$d_var) && all(is.finite(ms$d_var[s, ]))) {
          u$h_dev <- .emc_alloc_hmc_factor(ev, ms$W[[s]], ms$d_var[s, ], opts$min_evidence)
          if (is.null(u$h_dev)) u$h_dev <- h
        }
      }
    }
    if (!is.null(ms$d_rho) && all(is.finite(ms$d_rho[s, ]))) {
      u$d_rho <- ms$d_rho[s, ]
      if (!is.null(ms$p_dev) && all(is.finite(ms$p_dev[s, ]))) u$p_dev <- ms$p_dev[s, ]
      u$q_old <- if (is.null(ms$q_old)) 1 else ms$q_old[s]
      u$r_old <- if (is.null(ms$r_old)) 0 else ms$r_old[s]
    }
    if (priced) {
      u$rho <- ms$t_rho[s, ]
      u$p <- ms$p_fac[s, ]
      u$c_p <- ms$c_part[s]
    }
    units[[length(units) + 1L]] <- u
  }
  joint <- if (isTRUE(ms$coll)) {
    b <- bounds(if (is.null(ms$c_n)) 0 else ms$c_n)
    list(idx = seq_along(ms$g_rho), g = ms$c_fac,
         a = if (is.null(ms$c_sat)) ms$c_fac else ms$c_sat,
         cost = if (is.finite(ms$c_coll)) ms$c_coll else 0.1 * C_iter,
         fixed = if (isTRUE(is.finite(ms$coll_fixed))) ms$coll_fixed else 0,
         floor = b$floor,
         cap = min(b$cap, max(.EMC_ALLOC_MAX_ROUNDS, 2 * ms$c_rate)))
  }
  # A subject's HMC factor weighs its own transitions against the subjects'
  # pooled factor by the transitions behind each (the shear rule).
  n0 <- .EMC_ALLOC_ROUND_EVIDENCE
  hs <- lapply(units[-seq_along(ms$g_rho)], function(x) x$h)
  ns <- vapply(units[-seq_along(ms$g_rho)], function(x)
    if (is.null(x$h_n)) 0 else x$h_n, numeric(1))
  with_h <- which(ns > 0 & !vapply(hs, is.null, logical(1)))
  if (length(with_h)) {
    h_ref <- Reduce(`+`, lapply(with_h, function(i) ns[i] * hs[[i]])) / sum(ns[with_h])
    # The deviation factors have their own variances, so their own reference.
    hd <- lapply(units[-seq_along(ms$g_rho)], function(x) x$h_dev)
    with_d <- with_h[!vapply(hd[with_h], is.null, logical(1))]
    hd_ref <- if (length(with_d))
      Reduce(`+`, lapply(with_d, function(i) ns[i] * hd[[i]])) / sum(ns[with_d])
    for (i in with_h) {
      k <- length(ms$g_rho) + i
      units[[k]]$h <- (ns[i] * hs[[i]] + n0 * h_ref) / (ns[i] + n0)
      if (!is.null(units[[k]]$h_dev) && !is.null(hd_ref))
        units[[k]]$h_dev <- (ns[i] * units[[k]]$h_dev + n0 * hd_ref) / (ns[i] + n0)
    }
  }
  list(units = units, C_rest = C_rest, c_round = c_round, joint = joint,
       hmc_fixed = if (isTRUE(is.finite(ms$hmc_fixed))) ms$hmc_fixed else 0,
       shear_fixed = if (isTRUE(is.finite(ms$shear_fixed))) ms$shear_fixed else 0,
       r_bar = ms$r_bar, g_floor = ms$g_floor)
}

# A subject's deviation retention per iteration under particle rate q and HMC
# rate r with their joint selection (R/sampling.R): both run with probability
# max(0, q + r - 1), each alone with its rate less that, neither otherwise.
# P and H are the kernels' factors on the deviations.
.emc_alloc_dev_factor <- function(q, r, P, H) {
  both <- max(0, q + r - 1)
  (1 - q - r + both) + (q - both) * P + (r - both) * H + both * P * H
}

# The subjects' deviation autocorrelation under particle rates q and HMC
# rates r (per unit; ignored for other units) relative to the one measured,
# under the rates that ran (q_old, r_old). The measured statistic spans the
# particle step and the group's transports but not the HMC round, so its
# particle factor is divided out and every plan, the measured one included,
# multiplied by its full factor: the ratio is exactly 1 at the measured
# plan. A subject without a measured particle factor keeps its particle
# rate. NA without deviation measurements.
.emc_alloc_dev_ratio <- function(units, q, r) {
  num <- den <- 0
  for (i in seq_along(units)) {
    x <- units[[i]]
    if (!identical(x$kind, "subject") || is.null(x$d_rho) || !all(is.finite(x$d_rho))) next
    k <- length(x$d_rho)
    q0 <- if (is.null(x$q_old)) 1 else x$q_old
    r0 <- if (is.null(x$r_old)) 0 else x$r_old
    H <- if (is.null(x$h_dev)) rep(1, k) else x$h_dev
    P <- x$p_dev
    qi <- q[i]
    if (is.null(P)) {
      P <- rep(1, k)
      qi <- q0
    }
    ri <- if (is.null(x$h_dev)) 0 else r[i]
    base <- x$d_rho / pmax(1 - q0 * (1 - P), 1e-6)
    num <- num + sum(base * .emc_alloc_dev_factor(qi, ri, P, H))
    den <- den + sum(base * .emc_alloc_dev_factor(q0, r0, P, H))
  }
  if (den > 0) num / den else NA_real_
}

# The group coordinates' rho when the subjects' deviation autocorrelation
# moves by ratio from the measured one: rho_k' = f_k + (rho_k - f_k) ratio.
# Where the subjects slow down (ratio > 1) the extrapolation saturates
# instead, 1 - rho_k' = (1 - rho_k) exp(-(rho_k - f_k)(ratio - 1)/(1 - rho_k)):
# the same slope at the measured plan, never above the line and never at one,
# so no plan is priced as a coordinate that never mixes.
.emc_alloc_credit <- function(units, g_floor, ratio) {
  for (k in seq_along(g_floor)) {
    x <- units[[k]]
    s <- x$rho - g_floor[k]
    units[[k]]$rho <- if (ratio <= 1 || x$rho >= 1) min(1, g_floor[k] + s * ratio) else
      1 - (1 - x$rho) * exp(-s * (ratio - 1) / (1 - x$rho))
  }
  units
}

# Subject HMC for the group's sake: every subject with an HMC factor at rate
# at least r (r_min), its cost, and the group coordinates' rho the subjects'
# deviation autocorrelation implies at the particle rates in force
# (.emc_alloc_dev_ratio, credit). NULL without deviation measurements.
.emc_alloc_regime <- function(u, r) {
  units <- u$units
  cost <- 0
  covered <- FALSE
  q <- r_new <- numeric(length(units))
  for (i in seq_along(units)) {
    x <- units[[i]]
    if (!identical(x$kind, "subject")) next
    q[i] <- if (is.null(x$q_old)) 1 else x$q_old
    if (is.null(x$d_rho) || !all(is.finite(x$d_rho)) || is.null(x$h_dev) || is.null(x$h)) next
    covered <- TRUE
    units[[i]]$r_min <- r
    r_new[i] <- r
    cost <- cost + r * x$c_h
  }
  if (!covered) return(NULL)
  ratio <- .emc_alloc_dev_ratio(units, q, r_new)
  if (!is.finite(ratio)) return(NULL)
  list(units = .emc_alloc_credit(units, u$g_floor, ratio), cost = cost, credit = ratio)
}

# A solved plan's particle and total HMC rates per unit.
.emc_alloc_plan_rates <- function(units, sol) {
  n_u <- length(units)
  q <- if (length(sol$q) == n_u) sol$q else rep(1, n_u)
  r <- vapply(seq_len(n_u), function(i) {
    x <- units[[i]]
    if (!identical(x$kind, "subject") || is.null(x$h)) return(0)
    min(1, sol$x[i] + (if (is.null(x$r_min)) 0 else x$r_min))
  }, numeric(1))
  list(q = q, r = r)
}

# The group credit a plan was solved with assumed the particle rates in
# force; the rates it chose may move the subjects less (a particle step
# stood down for HMC that moves the deviations less). Where they do, the
# plan is solved again with the group credited at its own rates, and
# against the same plan with every particle step kept: the cheaper by the
# objective stands.
.emc_alloc_settle <- function(u, sol, units, credit, C_rest, movable, shear_act, hmc_fixed) {
  rates <- .emc_alloc_plan_rates(units, sol)
  ratio <- .emc_alloc_dev_ratio(units, rates$q, rates$r)
  if (!is.finite(ratio) || ratio <= credit * (1 + 1e-6)) return(sol)
  n_g <- length(u$g_floor)
  raw <- u$units
  hold <- function(us, q) {
    for (i in seq_along(us)) if (identical(us[[i]]$kind, "subject") &&
                                 !is.null(us[[i]]$p) && !is.null(us[[i]]$h))
      us[[i]]$q_fixed <- q[i]
    us
  }
  solve <- function(us) .emc_alloc_solve(us, C_rest, joint = u$joint, movable = movable,
                                         shear_fixed = shear_act, hmc_fixed = hmc_fixed)
  own <- units
  own[seq_len(n_g)] <- .emc_alloc_credit(raw[seq_len(n_g)], u$g_floor, ratio)
  own <- solve(hold(own, rates$q))
  # Every particle step kept, the group credited at those rates and the
  # plan's HMC rates (keeping the particle steps need not make up for HMC
  # the plan stood down).
  all_q <- rep(1, length(units))
  kept <- units
  ratio_k <- .emc_alloc_dev_ratio(units, all_q, rates$r)
  if (is.finite(ratio_k))
    kept[seq_len(n_g)] <- .emc_alloc_credit(raw[seq_len(n_g)], u$g_floor, ratio_k)
  kept <- solve(hold(kept, all_q))
  best <- if (isTRUE(own$objective < kept$objective)) own else kept
  best$regime <- sol$regime
  best
}

# The allocation with the HMC round's fixed cost: the best plan that runs the
# round against the best plan without HMC. A plan with HMC must clear the
# minimum gain against the iteration without any moves. The shear's fixed
# cost (the frames an iteration with rounds builds) is weighed the same way:
# the best plan with rounds against the best without. Each fixed cost is
# charged at the probability that an iteration runs any of the family's
# moves under the plan (.emc_alloc_solve).
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
      list(kind = "fixed", rho = x$rho, cost = 0)
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
  off <- lapply(u$units, function(x) {
    if (identical(x$kind, "hmc")) x <- .emc_alloc_as_subject(x)
    if (identical(x$kind, "subject")) x$h <- NULL
    x
  })
  movable <- if (!is.null(u$movable)) u$movable else
    .emc_alloc_movable(u$units, .EMC_ALLOC_MAX_ROUNDS, u$joint)
  has_h <- any(vapply(u$units, function(x) !is.null(x$h) || identical(x$kind, "hmc"),
                      logical(1)))
  # Subject HMC for the group's sake: regimes with every subject at a rate,
  # priced with the group autocorrelations they imply. A coordinate a regime
  # lowers counts as movable.
  regimes <- list()
  if (has_h && isTRUE(getOption("emc2.alloc_hmc_regimes", TRUE))) {
    levels0 <- .emc_alloc_range(u$units, .EMC_ALLOC_MAX_ROUNDS)$levels
    for (r in .EMC_ALLOC_HMC_REGIMES) {
      reg <- .emc_alloc_regime(u, r)
      if (is.null(reg)) break
      reg$r <- r
      regimes[[length(regimes) + 1L]] <- reg
      movable <- movable |
        (.emc_alloc_range(reg$units, .EMC_ALLOC_MAX_ROUNDS)$levels < levels0 - 1e-9)
    }
  }
  shear_act <- if (is.null(u$shear_act)) 0 else u$shear_act
  hmc_fixed <- if (is.null(u$hmc_fixed)) 0 else u$hmc_fixed
  # Without HMC every particle step runs, and the group is credited with the
  # subjects' movement that leaves (the measured group ran under the HMC
  # rates in force).
  if (has_h && !is.null(u$g_floor)) {
    ratio <- .emc_alloc_dev_ratio(u$units, rep(1, length(u$units)), numeric(length(u$units)))
    if (is.finite(ratio)) {
      n_g <- length(u$g_floor)
      off[seq_len(n_g)] <- .emc_alloc_credit(off[seq_len(n_g)], u$g_floor, ratio)
    }
  }
  without <- .emc_alloc_solve(off, u$C_rest, joint = u$joint, movable = movable,
                              shear_fixed = shear_act)
  if (!has_h) return(without)
  with <- .emc_alloc_solve(u$units, u$C_rest, joint = u$joint, movable = movable,
                           shear_fixed = shear_act, hmc_fixed = hmc_fixed)
  with$regime <- 0
  used <- list(units = u$units, credit = 1, C_rest = u$C_rest)
  for (reg in regimes) {
    cand <- .emc_alloc_solve(reg$units, u$C_rest + reg$cost, joint = u$joint,
                             movable = movable, shear_fixed = shear_act,
                             hmc_fixed = hmc_fixed)
    if (isTRUE(cand$objective < with$objective)) {
      cand$regime <- reg$r
      with <- cand
      used <- list(units = reg$units, credit = reg$credit, C_rest = u$C_rest + reg$cost)
    }
  }
  if (!is.null(u$g_floor))
    with <- .emc_alloc_settle(u, with, used$units, used$credit, used$C_rest, movable,
                              shear_act, hmc_fixed)
  if (isTRUE(with$objective < without$objective &&
             with$objective <= .EMC_ALLOC_MIN_GAIN * without$objective0)) {
    with$objective0 <- without$objective0
    with$tau0 <- without$tau0
    return(with)
  }
  without$regime <- 0
  without
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
  ms <- lapply(seq_len(nch), function(j)
    .emc_alloc_measure(emc, j, hst[[j]], Vbar[[j]], stage))
  shear_on <- vapply(seq_len(nch), function(j) ms[[j]]$shear, logical(1))
  # Calibration for subjects the optimistic allocation would use but that
  # lack measured HMC moves.
  if (hmc_on) {
    jobs <- NULL
    for (j in seq_len(nch)) {
      if (!isTRUE(hst[[j]]$supported)) next
      u <- .emc_alloc_units(ms[[j]], hst[[j]], shear_on[j], optimistic = TRUE, B = B)
      sol <- .emc_alloc_solve_units(u)
      x <- sol$x
      regime <- if (is.null(sol$regime)) 0 else sol$regime
      n_g <- length(ms[[j]]$g_rho)
      for (s in seq_len(emc[[j]]$n_subjects)) {
        sj <- hst[[j]]$subjects[[s]]
        if (!isTRUE(sj$eligible) || isTRUE(sj$warmup)) next
        wanted <- if (!is.null(opts$fixed_rate)) as.numeric(opts$fixed_rate) > 0 else
          regime + x[n_g + s] > 0
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
  # One chain's plan; q_fixed holds its particle rates to the shared ones
  # (second pass), so its HMC rates and predictions are solved for the plan
  # it runs.
  plan_chain <- function(j, q_fixed = NULL) {
    ch <- emc[[j]]
    m <- ms[[j]]
    st <- hst[[j]]
    if (!is.null(st) && isTRUE(st$supported)) for (s in seq_len(ch$n_subjects)) {
      if (isTRUE(st$subjects[[s]]$eligible) && !is.finite(st$subjects[[s]]$eps)) {
        st$subjects[[s]]$eligible <- FALSE
        st$subjects[[s]]$reason <- "no usable step size"
      }
    }
    u <- .emc_alloc_units(m, st, shear_on[j], optimistic = FALSE, B = B)
    n_g <- length(m$g_rho)
    if (!is.null(q_fixed)) for (s in seq_len(ch$n_subjects)) {
      k <- n_g + s
      if (!is.null(u$units[[k]]$p) && !is.null(u$units[[k]]$h))
        u$units[[k]]$q_fixed <- q_fixed[s]
    }
    sol <- .emc_alloc_solve_units(u)
    if (shear_on[j]) ch$shear$plan$rounds <- sol$x[seq_len(n_g)]
    if (isTRUE(m$coll)) ch$collapsed$plan$rounds <- sol$joint
    regime <- if (is.null(sol$regime)) 0 else sol$regime
    rate <- pmin(1, regime + sol$x[n_g + seq_len(ch$n_subjects)])
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
      c_hmc = m$c_hmc, c_part = m$c_part, hmc_fixed = u$hmc_fixed, shear_fixed = u$shear_fixed,
      fixed = m$fixed,
      tau = sol$tau, tau0 = sol$tau0,
      tau0_all = sol$tau0_all, cost = sol$cost,
      regime = regime, r_bar = m$r_bar,
      objective = sol$objective, objective0 = sol$objective0,
      group_rho = m$g_rho, group_g = m$g_fac, group_a = m$g_sat,
      rounds = if (shear_on[j]) sol$x[seq_len(n_g)] else numeric(n_g),
      collapsed_rounds = sol$joint, collapsed_g = if (isTRUE(m$coll)) m$c_fac,
      collapsed_a = if (isTRUE(m$coll)) m$c_sat,
      collapsed_cost = if (is.null(u$joint)) NA_real_ else u$joint$cost,
      collapsed_eps = if (isTRUE(m$coll)) ch$collapsed$plan$eps else NA_real_,
      subject_rho = apply(m$s_rho, 1L, max),
      hmc_h = vapply(u$units[n_g + seq_len(ch$n_subjects)], function(x)
        if (!is.null(x$h)) max(x$h) else NA_real_, numeric(1)),
      rate = if (is.null(st$subjects)) numeric(ch$n_subjects) else
        vapply(st$subjects, function(x) x$rate, numeric(1)),
      particle_rate = q,
      particle_fixed = !is.null(q_fixed),
      particle_p = apply(m$p_fac, 1L, function(x) if (all(is.finite(x))) max(x) else NA_real_),
      eps = if (is.null(st$subjects)) rep(NA_real_, ch$n_subjects) else
        vapply(st$subjects, function(x) x$eps, numeric(1)),
      calib_cpu = if (is.null(st$subjects)) 0 else
        sum(vapply(st$subjects, function(x) x$calib_cpu, numeric(1))))
    ch$alloc <- list(window = .emc_alloc_empty_window(ch$n_subjects, ch$n_pars),
                     stage = stage, pool = m$pool, last = rec, history = c(ch$alloc$history, list(rec)),
                     particle = list(rate = q, fac = m$p_fac, fac_dev = m$p_dev, cost = m$c_part))
    if (!is.null(st)) {
      st$alloc <- rec
      ch$subject_hmc <- st
    }
    ch
  }
  for (j in seq_len(nch)) emc[[j]] <- plan_chain(j)
  # One particle rate per subject for every chain, the largest any chain's
  # allocation chose, so chains carry the same work; a chain whose own rates
  # were lower is solved again at the shared ones.
  q <- do.call(pmax, lapply(emc, function(ch) ch$alloc$particle$rate))
  for (j in seq_len(nch))
    if (any(emc[[j]]$alloc$particle$rate < q - 1e-9)) emc[[j]] <- plan_chain(j, q_fixed = q)
  for (j in seq_len(nch)) {
    emc[[j]]$alloc$particle$rate <- q
    emc[[j]]$alloc$last$particle_shared <- q
    h <- emc[[j]]$alloc$history
    h[[length(h)]]$particle_shared <- q
    emc[[j]]$alloc$history <- h
  }
  # Collapsed rounds are allocated from each chain's own measured movement
  # and cost. Keep those rates local: chains have independent transition
  # kernels, and raising every rate to the cross-chain maximum adds work to
  # chains whose own allocation rejected the move.
  emc
}
