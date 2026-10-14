# Shear group update: the ensemble update (R/ensemble_move.R) followed by a
# shear of the group state.
#
# Target: pi(mu, V, a | y) propto p(mu) p(V) prod_s L_s(a_s) N(a_s | mu, V).
# The ensemble update draws the group state given the subjects (the centred
# augmentation). In a direction the data inform weakly it can only move the
# group as far as the subjects already have, and vice versa: the group mean
# along a ridge, and the group sd of a parameter that each subject's data
# barely determine (the hierarchical funnel). The shear moves the group
# state while each subject keeps its standardised position in its own
# conditional (the non-centred augmentation). Alternating the two
# augmentations of one group state is the interweaving strategy of Yu & Meng
# (2011); each step leaves the posterior invariant, so the composition does.
#
# The shear map, for a step u = (u_mu, u_sd):
#   mu -> mu + u_mu,   V -> D V D,   D = diag(exp(u_sd)),
#   a_s -> m_s' + R_s'^-1 R_s (a_s - m_s),
# where N(m_s, (R_s' R_s)^-1) is subject s's conditional posterior under its
# Gaussian likelihood approximation (lik_laplace, precision H_s, see
# .emc_lik_canonical) and prior N(mu, V), and primes denote the moved group
# state. A well-informed subject barely moves; a prior-dominated one is
# shifted and rescaled with the group. Without a usable approximation
# H_s = 0, so a_s -> mu' + D (a_s - mu). The moved subjects depend only on
# the start and end group states, so the maps compose, T_v T_u = T_(u+v),
# and T_-u inverts T_u. The Jacobian of T_u is
#   J_u = exp(sum_i (g_i + 1) u_sd,i) prod_s |R_s| / |R_s'|
# (g_i is the size of the covariance block holding parameter i; the prior
# density of V is taken over the free entries of its blocks), and it composes
# the same way. The approximation shapes the map only; a wrong one makes the
# moves smaller, not invalid.
#
# A round moves along one coordinate: a mean direction d (u = t (d, 0)) or
# the log group sd of one parameter (u = t e_i). Restricted to the orbit
# t -> T_tv(x) the posterior has density f(t) = pi(T_tv x) J_tv(x), and the
# orbit through any point of it gives the same f up to a shift in t. A round
# is the multiple-proposal step of Tjelmeland (2004) and Calderhead (2014) on
# t: an auxiliary point t~ ~ N(0, c^2), N candidates t_j ~ N(t~, c^2), all
# evaluated in one likelihood round, and one of the N + 1 points (the
# current one is t_0 = 0) chosen with probability proportional to f(t_j).
# Given t~ the N + 1 points are exchangeable, so the choice is a Gibbs step
# on their index and leaves the posterior invariant. A likelihood call's
# fixed cost dwarfs its extra rows, so N = 8 candidates per round cost about
# as much as a single Metropolis step and move further.
#
# Which coordinates. Under the approximations the centred Gibbs step moves
# the group mean as an autoregression with matrix
#   A = (n V^-1 + P0)^-1 V^-1 (sum_s Q_s^-1) V^-1,  Q_s = H_s + V^-1,
# whose eigenvectors are the mean directions (columns of dirs, each moving mu
# by one posterior sd of its direction under the approximations) and whose
# eigenvalues lambda_i predict their lag-1 autocorrelations; for a group sd,
# w_i = mean_s (Q_s^-1 V^-1)_ii is the prior share of each subject's
# conditional precision. How many rounds each coordinate gets per iteration
# is set by the move allocation (R/move_allocation.R) from what the rest of
# the iteration already mixes, at proposal rebuilds, and frozen in between.
# Mean steps keep V, and so every R_s: only the centres move.
#
# Under the Huang-Wand prior (standard type), V | a ~ IW(v + d - 1,
# 2 v diag(1 / a)) with a_k ~ IG(1/2, 1 / A_k^2). Holding a fixed would pin
# a scale step to the prior's conditional width in log sd whatever the data
# say. Given V the inverse-Wishart density factorises over the a_k, so a
# integrates out in closed form:
#   log p(V) = -(df + d + 1)/2 log|V| - (v + d)/2 sum_k log(v (V^-1)_kk + 1/A_k^2)
# (for an unblocked parameter, the d = 1 case on its variance). The rounds
# target the posterior with a integrated out, and a is then drawn from its
# conditional given the final V, the draw the Gibbs step makes: together a
# valid update of (V, a).
#
# Curved frames. One Gaussian likelihood approximation per subject shapes
# the map at every group state, so along a curved ridge of the subject
# likelihood (a direction in which the conditional's location and curvature
# change with mu) the transported subjects leave the ridge and the orbit
# collapses to a fraction of the posterior's width. The map is exact with
# any approximation that depends on the group state alone, so each subject
# gets its approximation as a function of mu: the likelihood's canonical
# form (prec, lin) refitted at the conditional mode under group states the
# chain visited (its atlas, accumulated over the rebuilds) and regressed on
# a quadratic in the slow mean directions of mu. The frames then follow each
# subject's ridge wherever the chain is, and the composition and inversion
# of the maps hold as before because the frames at a group state do not
# depend on the path to it. Mean steps with curved frames take the full map
# (the frames' R_s change), so their Jacobian is prod_s |R_s| / |R_s'| like
# a scale step's.
#
# Declared bounds. A model's bounds (bound$minmax under its transforms) are
# half-spaces u'a_s >= b of each subject on the sampled scale, one per
# design cell of a bounded parameter. Beyond one, every trial is at the
# likelihood floor and the subject's conditional is not the Gaussian the
# frame assumes, so a plain transport carries subjects sitting near a bound
# across it and the candidates there are lost. Each subject's wall is the
# half-space its window draws sit closest to; the transport then keeps the
# half-space: the complement of u moves by the affine map between the
# frames' marginals there, and u'a by the quantile map between the frames'
# conditionals given the complement, truncated at the wall. The maps still
# compose and invert, and the Jacobian is exact, so the rounds stay valid;
# far from the wall the map is an affine frame map. The collapsed move
# transports subjects with the same walls.
#
# Supported hierarchies: those in which every subject's prior mean is the
# group mean (standard without group designs, diagonal-gamma); the others
# keep the ensemble update alone.

# Candidates per round, and the share of rounds that should move (the
# optimum of the step scale c for N = 8 on a Gaussian orbit is flat around
# this value).
.EMC_SHEAR_CANDIDATES <- 8L
.EMC_SHEAR_MOVE_TARGET <- 0.68
# Curved frames: atlas nodes fitted per subject at each rebuild and the
# most a sample-stage top-up adds, the drift of the group log sds from the
# atlas's reference covariance that starts it over, the predicted lag-1
# autocorrelation from which a mean direction counts as slow, the measured
# per-round factor from which a slow direction's rounds count as
# ineffective (such directions get quadratic terms, at most this many), the
# nodes per coefficient the quadratic terms need, the ridge penalty on them,
# and the radius (in posterior sds) of nodes used.
.EMC_SHEAR_ATLAS_NODES <- 10L
.EMC_SHEAR_ATLAS_TOPUP <- 80L
.EMC_SHEAR_ATLAS_REANCHOR <- 0.5
.EMC_SHEAR_CURVE_LAM <- 0.5
.EMC_SHEAR_CURVE_TRIGGER <- 0.9
.EMC_SHEAR_CURVE_MAX_DIRS <- 4L
# The pilot that decides on an atlas: subjects and group means it fits, the
# orbit-variance gain it must show, and the rebuilds between repeats.
.EMC_SHEAR_PILOT_SUBJECTS <- 24L
.EMC_SHEAR_PILOT_GRID <- c(-2, -1, 0, 1, 2)
.EMC_SHEAR_CURVE_GAIN <- 2
.EMC_SHEAR_PILOT_EVERY <- 4L
.EMC_SHEAR_ATLAS_PER_COEF <- 3
.EMC_SHEAR_ATLAS_PENALTY <- 10
.EMC_SHEAR_ATLAS_RADIUS <- 12
# The least orbit spread (variance of a subject's orbit contribution over the
# pilot grid, in nats squared) for which a subject gets nodes, and the
# offsets (in window sds) of a log-sd atlas's nodes along each path.
.EMC_SHEAR_ATLAS_MIN_SPREAD <- 0.1
.EMC_SHEAR_ATLAS_PATH <- c(-2, -1, -0.5, 0.5, 1, 2)
# Initial step scale in posterior sds of the coordinate, the rounds whose
# move share weighs as much as a measured orbit width, the grid (in widths)
# a width is measured on, and the bounds every adapted scale stays in.
.EMC_SHEAR_STEP0 <- 2.5
.EMC_SHEAR_STEP_ROUNDS <- 20
.EMC_SHEAR_WIDTH_GRID <- seq(-3, 3, by = 0.5)
.EMC_SHEAR_LOG_STEP_RANGE <- log(c(0.05, 20))

.emc_shear_supported <- function(sampler) {
  isTRUE(sampler$type %in% c("standard", "diagonal-gamma")) &&
    !any(sampler$nuisance) &&
    is.null(sampler$gd) &&
    length(sampler$group_designs) == 0L &&
    isTRUE(sampler$n_subjects >= 2L)
}

# Hierarchies the shear runs on; the collapsed move shares the structural
# check (.emc_collapsed_active).
.emc_shear_active <- function(sampler) .emc_shear_supported(sampler)

# Each subject's likelihood approximation if it is finite, else NULL. The
# particle step's Pareto check does not apply: it judges the Gaussian
# conditional at one group state, while the map is valid with any
# approximation and moves an informative subject far better with a rough one
# than with none (H_s = 0).
.emc_shear_fits <- function(sampler) {
  p <- sampler$n_pars
  lapply(seq_len(sampler$n_subjects), function(s) {
    fit <- if (length(sampler$lik_laplace) >= s) sampler$lik_laplace[[s]] else NULL
    ok <- is.list(fit) && is.matrix(fit$prec) && all(dim(fit$prec) == p) &&
      length(fit$lin) == p && all(is.finite(fit$prec)) && all(is.finite(fit$lin))
    if (ok) fit else NULL
  })
}

# ---- Declared bounds as per-subject half-spaces ---------------------------
# A bounded parameter type (bound$minmax on the natural scale, reached from
# the sampled scale through its monotone transform) bounds each design cell's
# linear predictor, a half-space u'alpha >= b on the sampled scale (constant
# columns shift b). Per subject, the list of (u, b); NULL without a model
# contract. A bound the transform already enforces (an exp lower bound at
# the transform's own floor) is no constraint.
.emc_shear_support <- function(sampler) {
  model <- sampler$model
  if (!is.function(model)) return(NULL)
  ml <- tryCatch(model(), error = function(e) NULL)
  mm <- ml$bound$minmax
  tf <- ml$transform
  if (!is.matrix(mm) || !length(colnames(mm)) || is.null(tf$func)) return(NULL)
  pn <- sampler$par_names
  p <- length(pn)
  limits <- function(pt) {
    func <- tf$func[[pt]]
    lo <- unname(mm[1L, pt]); hi <- unname(mm[2L, pt])
    low <- if (is.null(tf$lower[[pt]]) || !is.finite(tf$lower[[pt]])) 0 else tf$lower[[pt]]
    up <- if (is.null(tf$upper[[pt]]) || !is.finite(tf$upper[[pt]])) 1 else tf$upper[[pt]]
    if (identical(func, "exp")) {
      c(if (lo > low) log(lo - low) else -Inf, if (is.finite(hi)) log(hi - low) else Inf)
    } else if (identical(func, "identity")) {
      c(lo, hi)
    } else if (identical(func, "pnorm")) {
      f <- function(x) stats::qnorm(min(1, max(0, (x - low) / (up - low))))
      c(if (lo > low) f(lo) else -Inf, if (hi < up) f(hi) else Inf)
    } else c(-Inf, Inf)
  }
  lapply(sampler$data, function(dadm) {
    out <- list()
    designs <- attr(dadm, "designs")
    consts <- attr(dadm, "constants")
    if (!is.list(designs)) return(out)
    for (pt in colnames(mm)) {
      X <- designs[[pt]]
      if (!is.matrix(X) || !ncol(X)) next
      lim <- limits(pt)
      if (!any(is.finite(lim))) next
      cols <- match(colnames(X), pn)
      known <- !is.na(cols) | colnames(X) %in% names(consts)
      if (!all(known)) next
      for (r in seq_len(nrow(X))) {
        x <- X[r, ]
        u <- numeric(p)
        u[cols[!is.na(cols)]] <- x[!is.na(cols)]
        if (all(u == 0)) next
        shift <- if (any(is.na(cols))) sum(x[is.na(cols)] * consts[colnames(X)[is.na(cols)]]) else 0
        if (is.finite(lim[1L]))
          out[[length(out) + 1L]] <- list(u = u, b = lim[1L] - shift)
        if (is.finite(lim[2L]))
          out[[length(out) + 1L]] <- list(u = -u, b = -(lim[2L] - shift))
      }
    }
    out
  })
}

# Each subject's transport wall for a block: of its half-spaces, the one its
# window draws sit closest to (smallest mean margin in units of the margin's
# spread). The map is an ordinary frame map where the wall is far, so every
# subject with a half-space gets one. Columns of U are the unit normals,
# b their offsets, Uperp an orthonormal basis of each normal's complement.
.emc_shear_walls <- function(chain, idx) {
  support <- .emc_shear_support(chain)
  if (is.null(support)) return(NULL)
  n <- chain$n_subjects
  p <- chain$n_pars
  U <- matrix(0, p, n)
  b <- rep(NA_real_, n)
  Uperp <- vector("list", n)
  has <- logical(n)
  for (s in seq_len(n)) {
    cs <- support[[s]]
    if (!length(cs)) next
    A <- t(matrix(chain$samples$alpha[, s, idx], p))
    margin <- vapply(cs, function(c) {
      v <- drop(A %*% c$u) - c$b
      mean(v) / max(stats::sd(v), 1e-12)
    }, numeric(1))
    k <- which.min(margin)
    nu <- sqrt(sum(cs[[k]]$u^2))
    U[, s] <- cs[[k]]$u / nu
    b[s] <- cs[[k]]$b / nu
    Uperp[[s]] <- qr.Q(qr(cbind(U[, s], diag(p))))[, -1L, drop = FALSE]
    has[s] <- TRUE
  }
  if (!any(has)) return(NULL)
  list(U = U, b = b, Uperp = Uperp, has = has)
}

# Mean directions and scale weights at group covariance V (see the header).
# Columns of dirs move mu by one posterior sd of their direction under the
# approximations: mu = mu0 + t d with d' (n V^-1 + P0) d = 1 has collapsed
# precision 1 - lambda in t.
.emc_shear_plan <- function(fits, prior, V) {
  p <- nrow(V)
  n <- length(fits)
  Vi <- chol2inv(chol(V))
  S <- matrix(0, p, p)
  w <- numeric(p)
  for (fit in fits) {
    Qi <- chol2inv(chol(if (is.null(fit)) Vi else fit$prec + Vi))
    S <- S + Qi
    w <- w + rowSums(Qi * t(Vi))
  }
  Pc <- n * Vi + as.matrix(prior$theta_mu_invar)
  R <- chol((Pc + t(Pc)) / 2)
  Ri <- backsolve(R, diag(p))
  M <- t(Ri) %*% (Vi %*% S %*% Vi) %*% Ri
  e <- eigen((M + t(M)) / 2, symmetric = TRUE)
  lam <- pmin(pmax(e$values, 0), 1 - 1e-6)
  dirs <- sweep(Ri %*% e$vectors, 2L, sqrt(1 - lam), "/")
  if (!all(is.finite(dirs))) return(NULL)
  list(dirs = dirs, lam = lam, w = pmin(1, pmax(0, w / n)))
}

# The move along coordinate k (1..p mean directions, p+1..2p log group sds)
# as a step direction in u = (u_mu, u_sd).
.emc_shear_coord_dir <- function(plan, k) {
  p <- nrow(plan$dirs)
  if (k <= p) c(plan$dirs[, k], numeric(p)) else replace(numeric(2L * p), k, 1)
}

# Every chain's plan from the proposal window create_chain_proposals() uses.
# Step scales and the last block's round records carry over to the new
# directions they match (|cosine| >= 0.9; a log group sd always matches
# itself); the round counts are set by the move allocation afterwards.
.emc_shear_set_plans <- function(emc, samples_idx, stage, n_cores = 1L) {
  if (!.emc_shear_active(emc[[1L]])) return(emc)
  idx <- .emc_chain_proposal_window(emc, samples_idx, stage)
  if (!length(idx)) return(emc)
  for (j in seq_along(emc)) emc[[j]] <- .emc_shear_set_plan(emc[[j]], idx, stage)
  emc <- .emc_shear_atlas_add(emc, idx, n_cores, stage)
  emc <- .emc_shear_build_curves(emc, idx, stage)
  for (j in seq_along(emc)) emc[[j]] <- .emc_shear_curve_validate(emc[[j]], idx, n_cores)
  .emc_shear_orbit_steps(emc, idx, n_cores)
}

# The orbit of the current state along coordinate k of the plan, at the
# offsets grid (in .emc_shear_coord_path's units): per subject its
# contribution to the orbit's log density (likelihood, hierarchical density
# and the Jacobian of the move's own transport, walls included; NA where its
# transport or likelihood fails), and the group terms (mean and covariance
# priors, and the covariance scaling's Jacobian along a log sd). NULL when
# the coordinate has no path; every contribution NA when a frame fails.
# fits_of gives the fits a path state's frames use.
.emc_shear_orbit_terms <- function(ch, idx, k, grid, fits_of, n_cores = 1L) {
  s <- ch$samples
  i <- s$idx
  p <- ch$n_pars
  n <- ch$n_subjects
  mu <- as.numeric(s$theta_mu[, i])
  V <- s$theta_var[, , i]
  alpha <- matrix(s$alpha[, , i], p)
  spec <- .emc_shear_v_spec(ch)
  walls <- .emc_shear_walls_of(ch)
  g0 <- which(grid == 0)
  par_names <- rownames(s$alpha)
  path <- .emc_shear_coord_path(ch, idx, k, mu, V, grid)
  if (is.null(path) || length(g0) != 1L) return(NULL)
  for (g in seq_along(grid))
    path[[g]]$chol <- tryCatch(chol(path[[g]]$V), error = function(e) NULL)
  if (any(vapply(path, function(st) is.null(st$chol), logical(1)))) return(NULL)
  states <- lapply(path, function(st)
    .emc_shear_state(st$mu, st$V, fits_of(st), ch$prior, spec))
  if (any(vapply(states, is.null, logical(1))))
    return(list(per = matrix(NA_real_, length(grid), n), lp = rep(NA_real_, length(grid))))
  ldn <- function(a, st) {
    z <- backsolve(st$chol, a - st$mu, transpose = TRUE)
    -0.5 * sum(z^2) - sum(log(diag(st$chol)))
  }
  tr <- lapply(states, function(st) .emc_shear_transport(alpha, states[[g0]], st, walls))
  per <- .emc_with_preserved_rng(auto_mclapply(seq_len(n), function(ss) {
    if (!all(vapply(tr, function(x) isTRUE(x$ok[ss]), logical(1))))
      return(rep(NA_real_, length(grid)))
    X <- t(vapply(tr, function(x) x$alpha[, ss], numeric(p)))
    ll <- tryCatch(.emc_burn_ll_batch(X, ch$data[[ss]], ch$model, par_names),
                   error = function(e) rep(NA_real_, length(grid)))
    vapply(seq_along(grid), function(g)
      ll[g] + ldn(X[g, ], path[[g]]) + states[[g0]]$fm$logdet[ss] -
        states[[g]]$fm$logdet[ss] + tr[[g]]$jac[ss], numeric(1))
  }, mc.cores = max(1L, min(n, n_cores)), mc.set.seed = FALSE))
  lp <- vapply(path, function(st) .emc_shear_log_prior_mu(st$mu, ch$prior) +
    (if (k <= p) 0 else tryCatch(.emc_shear_log_prior_v(st$V, spec), error = function(e) NA_real_) +
       sum((spec$size + 1) * st$u_sd)), numeric(1))
  list(per = matrix(unlist(per), length(grid)), lp = lp)
}

# Each subject's orbit spread along the plan's coordinates ks at the current
# state: the variance over the pilot grid of its orbit contribution
# (.emc_shear_orbit_terms), summed over ks, one column per entry of fits_of.
# A frame that fails sets its column to Inf, as does a subject whose
# transport or likelihood fails.
.emc_shear_orbit_spread <- function(ch, idx, ks, fits_of, n_cores = 1L) {
  spread <- matrix(0, ch$n_subjects, length(fits_of))
  for (k in ks) {
    for (kind in seq_along(fits_of)) {
      tm <- .emc_shear_orbit_terms(ch, idx, k, .EMC_SHEAR_PILOT_GRID, fits_of[[kind]], n_cores)
      if (is.null(tm)) break
      v <- apply(tm$per, 2L, function(x) if (all(is.finite(x))) stats::var(x) else Inf)
      spread[, kind] <- spread[, kind] + v
    }
  }
  spread
}

# The sd of the orbit's density along coordinate k (in .emc_shear_coord_path's
# units), on a grid scaled to the width s0 and rescaled until the grid
# resolves the orbit; NA when it cannot be measured.
.emc_shear_orbit_width <- function(ch, idx, k, fits_of, s0, n_cores = 1L) {
  s <- s0
  sd <- NA_real_
  for (it in seq_len(4L)) {
    grid <- s * .EMC_SHEAR_WIDTH_GRID
    tm <- .emc_shear_orbit_terms(ch, idx, k, grid, fits_of, n_cores)
    if (is.null(tm)) return(NA_real_)
    f <- rowSums(tm$per) + tm$lp
    f[!is.finite(f)] <- -Inf
    if (!is.finite(max(f))) return(NA_real_)
    w <- exp(f - max(f))
    w <- w / sum(w)
    m <- sum(w * grid)
    sd <- sqrt(sum(w * (grid - m)^2))
    edge <- w[1L] + w[length(w)]
    if (edge > 0.01) {
      # Wider than the grid: the log density's curvature on it estimates
      # the width.
      ok <- is.finite(f)
      a <- if (sum(ok) >= 3L) stats::coef(stats::lm(f[ok] ~ grid[ok] + I(grid[ok]^2)))[3L] else NA
      s <- if (is.finite(a) && a < 0) max(2 * s, 1 / sqrt(-2 * a)) else 4 * s
    } else if (sd < 0.5 * s) s <- max(sd, s / 8)
    else return(sd)
  }
  sd
}

# Each coordinate's step scale from its orbit under the block's frames: the
# multiple-proposal step is best near .EMC_SHEAR_STEP0 orbit sds, so the
# measured width sets the step, weighed against the step the rounds adapted
# by how many rounds it rests on (.EMC_SHEAR_STEP_ROUNDS count as much as
# the width). A coordinate with few rounds would otherwise keep whatever step it
# last had, and a step far from the orbit's width makes its rounds look
# ineffective, which then buys it fewer rounds still. One job per chain and
# coordinate, each measured within its worker.
.emc_shear_orbit_steps <- function(emc, idx, n_cores = 1L) {
  n0 <- .EMC_SHEAR_STEP_ROUNDS
  jobs <- list()
  for (j in seq_along(emc)) {
    plan <- emc[[j]]$shear$plan
    if (!is.list(plan) || !is.matrix(plan$dirs)) next
    p <- emc[[j]]$n_pars
    n_ev <- if (is.list(plan$last) && length(plan$last$n) == 2L * p) plan$last$n else numeric(2L * p)
    sd_l <- apply(.emc_shear_window_lsd(emc[[j]], idx), 1L, stats::sd)
    for (k in which(n_ev < 4 * n0)) {
      unit <- if (k <= p) 1 else sd_l[k - p]
      if (!is.finite(unit) || unit <= 0) next
      jobs[[length(jobs) + 1L]] <- list(chain = j, k = k, unit = unit, n = n_ev[k],
        s0 = exp(plan$log_step[k]) / .EMC_SHEAR_STEP0 / unit)
    }
  }
  if (!length(jobs)) return(emc)
  contexts <- vector("list", length(emc))
  for (j in unique(vapply(jobs, `[[`, integer(1), "chain")))
    contexts[[j]] <- .emc_shear_context(emc[[j]])
  widths <- .emc_with_preserved_rng(.emc_rebuild_lapply(seq_along(jobs), function(m) {
    jb <- jobs[[m]]
    ch <- emc[[jb$chain]]
    ctx <- contexts[[jb$chain]]
    tryCatch(.emc_shear_orbit_width(ch, idx, jb$k,
                                    function(st) .emc_shear_fits_at(ctx, st$mu, st$V), jb$s0),
             error = function(e) NA_real_)
  }, n_cores = n_cores))
  for (m in seq_along(jobs)) {
    jb <- jobs[[m]]
    w <- widths[[m]]
    if (!is.numeric(w) || length(w) != 1L || !is.finite(w) || w <= 0) next
    plan <- emc[[jb$chain]]$shear$plan
    target <- log(.EMC_SHEAR_STEP0 * w * jb$unit)
    plan$log_step[jb$k] <- min(.EMC_SHEAR_LOG_STEP_RANGE[2], max(.EMC_SHEAR_LOG_STEP_RANGE[1],
      (jb$n * plan$log_step[jb$k] + n0 * target) / (jb$n + n0)))
    emc[[jb$chain]]$shear$plan <- plan
  }
  emc
}

# A curved subject keeps its curve only where its orbit contribution along
# the curve's coordinates (its slow mean directions and log sds), at the
# current state, varies less than under the block's fixed frame: the
# subject's own likelihood at the transported points decides, so a subject
# whose node fits the regression cannot smooth (an unstable likelihood mode)
# stays on the fixed frame. Each frame is a function of the group state, so
# the choice, fixed for the block, keeps the move exact.
.emc_shear_curve_validate <- function(ch, idx, n_cores = 1L) {
  cv <- ch$shear$curve
  if (is.null(cv) || !(length(cv$q) || length(cv$l))) return(ch)
  ctx <- .emc_shear_context(ch)
  if (is.null(ctx$curve)) return(ch)
  spread <- .emc_shear_orbit_spread(
    ch, idx, c(cv$q, ch$n_pars + cv$l),
    list(function(st) .emc_shear_fits_at(ctx, st$mu, st$V), function(st) ctx$packed),
    n_cores)
  drop <- which(!(spread[, 1L] < spread[, 2L]))
  if (!length(drop)) return(ch)
  B <- cv$B
  B[drop] <- list(NULL)
  if (!any(!vapply(B, is.null, logical(1)))) {
    ch$shear$curve <- NULL
    return(ch)
  }
  cv$B <- B
  cv$sets <- .emc_shear_curve_sets(B, cv$q)
  cv$dropped <- drop
  ch$shear$curve <- cv
  ch
}

# One subject's likelihood approximation at the conditional mode under the
# prior N(mu, V), from the start x0; NULL when the mode or its curvature is
# not found.
.emc_shear_node_fit <- function(chain, s, mu, V, x0) {
  par_names <- rownames(chain$samples$alpha)
  ll_batch <- function(x) .emc_burn_ll_batch(x, chain$data[[s]], chain$model, par_names)
  r <- tryCatch(.emc_optimize_subject(x0, list(mu = mu, var = V), ll_batch,
                                      blocks = .emc_joint_ll_blocks(chain, s)),
                error = function(e) NULL)
  if (!is.list(r) || !is.list(r$lik)) return(NULL)
  r$lik[c("prec", "lin")]
}

.emc_shear_set_plan <- function(chain, idx, stage) {
  if (!.emc_shear_active(chain)) return(chain)
  n <- chain$n_subjects
  p <- chain$n_pars
  V <- apply(chain$samples$theta_var[, , idx, drop = FALSE], 1:2, mean)
  plan <- tryCatch(.emc_shear_plan(.emc_shear_fits(chain), chain$prior, V),
                   error = function(e) NULL)
  if (is.null(plan)) return(chain)
  st <- chain$shear
  if (is.null(st)) st <- list()
  old <- st$plan
  base <- c(rep(0, p), rep(-0.5 * log(2 * n), p))
  plan$log_step <- base + log(.EMC_SHEAR_STEP0)
  plan$rounds <- numeric(2L * p)
  plan$stage <- stage
  zero <- .emc_shear_record(NULL, p)
  recent <- plan$pool <- zero
  # Production blocks run settled kernels on a stationary chain, so their
  # rounds pool over the sample stage (from its first block on).
  pooled <- identical(stage, "sample") && identical(old$stage, "sample")
  if (is.list(old) && is.matrix(old$dirs) && all(dim(old$dirs) == p)) {
    scale <- c(rep(NA_real_, p), rep(1, p))
    match <- c(rep(NA_integer_, p), p + seq_len(p))
    for (k in seq_len(p)) {
      d <- plan$dirs[, k]
      cs <- abs(drop(crossprod(old$dirs, d))) /
        (sqrt(colSums(old$dirs^2)) * sqrt(sum(d^2)))
      m <- which.max(cs)
      if (cs[m] >= 0.9) {
        match[k] <- m
        scale[k] <- sqrt(sum(old$dirs[, m]^2) / sum(d^2))
      }
    }
    got <- which(!is.na(match))
    blk <- old$block
    plan$log_step[got] <- old$log_step[match[got]] + log(scale[got])
    carry <- function(rec) {
      out <- zero
      if (is.null(rec)) return(out)
      rec <- .emc_shear_record(rec, p)
      out$E[got] <- rec$E[match[got]] * scale[got]^2
      out$n[got] <- rec$n[match[got]]
      out$moved[got] <- rec$moved[match[got]]
      out$N[got] <- rec$N[match[got]] * scale[got]^2
      out$iters[got] <- rec$iters[match[got]]
      out$K <- matrix(0L, 2L * p, ncol(rec$K))
      out$K[got, ] <- rec$K[match[got], ]
      out
    }
    recent <- carry(blk)
    recent$iters[got] <- .emc_alloc_block_iters(chain)
    if (pooled) {
      prev <- carry(old$pool)
      d <- .EMC_ALLOC_EVIDENCE_DECAY
      plan$pool <- list(E = d * prev$E + recent$E, n = d * prev$n + recent$n,
                        moved = d * prev$moved + recent$moved, N = d * prev$N + recent$N,
                        K = .emc_alloc_khist_sum(d * prev$K, recent$K),
                        iters = d * prev$iters + recent$iters)
    }
  }
  # The allocation measures a coordinate's rounds from the pooled production
  # blocks, or else the last block alone. Step scales follow the last block.
  plan$last <- if (pooled) plan$pool else recent
  # Discarded stages adapt the scales round by round; production only here.
  if (identical(stage, "sample")) plan <- .emc_shear_tune_steps(plan, recent)
  plan$block <- zero
  st$plan <- plan
  st$walls <- tryCatch(.emc_shear_walls(chain, idx), error = function(e) NULL)
  chain$shear <- st
  chain
}

# The slow mean directions whose rounds with the frames the chain has were
# measured to barely move them (the orbit is a fraction of the posterior's
# width): the directions that get quadratic terms, the slowest first.
.emc_shear_curve_dirs <- function(chain, idx) {
  plan <- chain$shear$plan
  rec <- plan$last
  if (!is.list(rec) || is.null(rec$n)) return(integer(0))
  slow <- which(plan$lam >= .EMC_SHEAR_CURVE_LAM)
  if (!length(slow)) return(integer(0))
  Dinv <- tryCatch(solve(plan$dirs), error = function(e) NULL)
  if (is.null(Dinv)) return(integer(0))
  t <- Dinv %*% chain$samples$theta_mu[, idx, drop = FALSE]
  v <- apply(t, 1L, stats::var)
  n <- rec$n[slow]
  g <- 1 - (rec$E[slow] / pmax(1L, n)) / (2 * pmax(v[slow], .Machine$double.xmin))
  q <- slow[n >= 5L & is.finite(g) & g >= .EMC_SHEAR_CURVE_TRIGGER]
  q[seq_len(min(length(q), .EMC_SHEAR_CURVE_MAX_DIRS))]
}

# The group log sds whose rounds were measured to barely move them, the
# slowest first: the log-sd coordinates the curved frames follow. Frames
# fitted at one group state go stale as such a coordinate travels, which no
# pilot at the current state can show, so the measurement alone decides.
# Only the sample stage follows them: the discarded stages refit every
# subject's approximation each block, and their log sds travel too far for
# an atlas to settle.
.emc_shear_curve_lsd <- function(chain, idx, stage) {
  if (!identical(stage, "sample")) return(integer(0))
  rec <- chain$shear$plan$last
  if (!is.list(rec) || is.null(rec$n) || length(idx) < 2L) return(integer(0))
  p <- chain$n_pars
  v <- apply(.emc_shear_window_lsd(chain, idx), 1L, stats::var)
  k <- p + seq_len(p)
  n <- rec$n[k]
  g <- 1 - (rec$E[k] / pmax(1L, n)) / (2 * pmax(v, .Machine$double.xmin))
  l <- which(n >= 5L & is.finite(g) & g >= .EMC_SHEAR_CURVE_TRIGGER)
  l <- l[order(-g[l])]
  l[seq_len(min(length(l), .EMC_SHEAR_CURVE_MAX_DIRS))]
}

# Whether a chain's frames need curving: a slow direction whose rounds are
# ineffective, and a pilot showing that frames following the subjects'
# conditional modes widen that direction's orbit materially. The pilot fits
# the modes of a subsample of subjects at a few group means along the
# direction and compares the orbit's width under those frames with the
# width under the block's fits; it is repeated on a rebuild cadence while
# the rounds stay ineffective. Once an atlas exists it is kept.
.emc_shear_curve_wanted <- function(chain, idx, n_cores) {
  q <- .emc_shear_curve_dirs(chain, idx)
  if (!length(q)) return(list(wanted = FALSE, ran = FALSE))
  pilot <- chain$shear$pilot
  if (is.list(pilot) && pilot$rebuilds < .EMC_SHEAR_PILOT_EVERY)
    return(list(wanted = FALSE, ran = FALSE))
  r <- .emc_shear_curve_pilot(chain, idx, q[1L], n_cores)
  list(wanted = isTRUE(is.finite(r) && r >= .EMC_SHEAR_CURVE_GAIN), ran = TRUE, gain = r)
}

# Rebuilds since the last pilot.
.emc_shear_pilot_tick <- function(chain, ran) {
  pilot <- chain$shear$pilot
  chain$shear$pilot <- if (ran) list(rebuilds = 0L) else
    list(rebuilds = (if (is.list(pilot)) pilot$rebuilds else 0L) + 1L)
  chain
}

# The pilot: the orbit variance along coordinate k (a mean direction or a
# log sd) under exact conditional-mode frames over the variance under the
# block's fits, both over a subsample of subjects at the current state (the
# orbit's log density is a sum over subjects, so the subsample's own orbit is
# compared).
.emc_shear_curve_pilot <- function(chain, idx, k, n_cores) {
  s <- chain$samples
  p <- chain$n_pars
  n <- chain$n_subjects
  i <- s$idx
  mu <- as.numeric(s$theta_mu[, i])
  V <- apply(s$theta_var[, , idx, drop = FALSE], 1:2, mean)
  grid <- .EMC_SHEAR_PILOT_GRID
  path <- .emc_shear_coord_path(chain, idx, k, mu, V, grid)
  if (is.null(path)) return(NA_real_)
  spec <- .emc_shear_v_spec(chain)
  for (g in seq_along(grid)) {
    path[[g]]$chol <- tryCatch(chol(path[[g]]$V), error = function(e) NULL)
    if (is.null(path[[g]]$chol)) return(NA_real_)
    path[[g]]$Vi <- chol2inv(path[[g]]$chol)
    # A log-sd path adds the covariance prior and the scaling's Jacobian;
    # along a mean direction both are constant.
    path[[g]]$lp <- .emc_shear_log_prior_mu(path[[g]]$mu, chain$prior) +
      if (k <= p) 0 else
        tryCatch(.emc_shear_log_prior_v(path[[g]]$V, spec), error = function(e) NA_real_) +
          sum((spec$size + 1) * path[[g]]$u_sd)
  }
  alpha <- matrix(s$alpha[, , i], p)
  subs <- unique(round(seq(1L, n, length.out = min(n, .EMC_SHEAR_PILOT_SUBJECTS))))
  par_names <- rownames(s$alpha)
  jobs <- expand.grid(subject = subs, g = seq_along(grid))
  fits <- .emc_with_preserved_rng(auto_mclapply(seq_len(nrow(jobs)), function(j) {
    ss <- jobs$subject[j]
    st <- path[[jobs$g[j]]]
    .emc_quiet_fork(.emc_shear_node_fit(chain, ss, st$mu, st$V, alpha[, ss]))
  }, mc.cores = max(1L, min(nrow(jobs), n_cores)), mc.set.seed = FALSE))
  frame <- function(fit, st) {
    R <- tryCatch(chol(fit$prec + st$Vi), error = function(e) NULL)
    if (is.null(R)) return(NULL)
    list(m = drop(backsolve(R, forwardsolve(t(R), fit$lin + drop(st$Vi %*% st$mu)))), R = R)
  }
  base <- .emc_shear_fits(chain)
  ok <- subs
  ex <- lapply(seq_along(grid), function(g) lapply(subs, function(ss) {
    f <- fits[[which(jobs$subject == ss & jobs$g == g)]]
    if (is.null(f)) NULL else frame(f, path[[g]])
  }))
  for (g in seq_along(grid)) ok <- ok[!vapply(ex[[g]][match(ok, subs)], is.null, logical(1))]
  if (length(ok) < 2L) return(NA_real_)
  st <- lapply(seq_along(grid), function(g) lapply(subs, function(ss)
    frame(if (is.null(base[[ss]])) list(prec = matrix(0, p, p), lin = numeric(p)) else base[[ss]],
          path[[g]])))
  g0 <- which(grid == 0)
  moved <- function(frames, g, ss) {
    f0 <- frames[[g0]][[match(ss, subs)]]
    f1 <- frames[[g]][[match(ss, subs)]]
    if (is.null(f0) || is.null(f1)) return(NULL)
    list(a = f1$m + drop(backsolve(f1$R, f0$R %*% (alpha[, ss] - f0$m))),
         jac = sum(log(diag(f0$R))) - sum(log(diag(f1$R))))
  }
  ldn <- function(a, st) {
    z <- backsolve(st$chol, a - st$mu, transpose = TRUE)
    -0.5 * sum(z^2) - sum(log(diag(st$chol)))
  }
  rows <- expand.grid(subject = ok, g = seq_along(grid), kind = 1:2)
  ll <- .emc_with_preserved_rng(auto_mclapply(seq_len(nrow(rows)), function(j) {
    ss <- rows$subject[j]
    mv <- moved(if (rows$kind[j] == 1L) ex else st, rows$g[j], ss)
    if (is.null(mv)) return(-Inf)
    v <- tryCatch(.emc_burn_ll_batch(mv$a, chain$data[[ss]], chain$model, par_names),
                  error = function(e) NA_real_)
    if (!is.finite(v)) -Inf else v + ldn(mv$a, path[[rows$g[j]]]) + mv$jac
  }, mc.cores = max(1L, min(nrow(rows), n_cores)), mc.set.seed = FALSE))
  ll <- unlist(ll)
  width <- function(kind) {
    f <- vapply(seq_along(grid), function(g) {
      sum(ll[rows$g == g & rows$kind == kind]) + path[[g]]$lp
    }, numeric(1))
    if (!all(is.finite(f))) return(NA_real_)
    w <- exp(f - max(f))
    w <- w / sum(w)
    m <- sum(w * grid)
    sum(w * (grid - m)^2)
  }
  width(1L) / width(2L)
}

# The group log sds over the draws idx (p x length(idx)).
.emc_shear_window_lsd <- function(chain, idx) {
  0.5 * log(matrix(apply(chain$samples$theta_var[, , idx, drop = FALSE], 3L, diag),
                   chain$n_pars))
}

# The group states at offsets grid along coordinate k of the plan from
# (mu, V), with the offsets in the coordinate's posterior sds over the
# window idx: a mean direction moves mu, a log sd scales V. u_sd is the
# log-sd part of the move.
.emc_shear_coord_path <- function(chain, idx, k, mu, V, grid) {
  p <- chain$n_pars
  if (k <= p) {
    d <- chain$shear$plan$dirs[, k]
    return(lapply(grid, function(t) list(mu = mu + t * d, V = V, u_sd = numeric(p))))
  }
  j <- k - p
  sd_j <- stats::sd(.emc_shear_window_lsd(chain, idx)[j, ])
  if (!is.finite(sd_j) || sd_j <= 0) return(NULL)
  lapply(grid, function(t) {
    u_sd <- replace(numeric(p), j, t * sd_j)
    e <- exp(u_sd)
    list(mu = mu, V = V * outer(e, e), u_sd = u_sd)
  })
}

# Quadratic features of a group state's curve coordinates t, with the
# quadratic terms over the mean directions q.
.emc_shear_atlas_features <- function(t, q) {
  f <- c(1, t)
  if (length(q)) {
    tq <- t[q]
    f <- c(f, outer(tq, tq)[upper.tri(diag(length(q)), diag = TRUE)])
  }
  f
}

# The subjects whose nodes this block fits: those whose orbit contribution
# along the curve coordinates ks, under the block's fixed frames, varies by
# at least .EMC_SHEAR_ATLAS_MIN_SPREAD over the grid, together with any whose
# transport or likelihood fails along the orbit. A subject whose contribution
# is flatter barely shapes the orbit's log density, curved or not. Without
# coordinates to measure every subject is fitted.
.emc_shear_atlas_subjects <- function(ch, idx, ks, n_cores) {
  n <- ch$n_subjects
  if (!length(ks)) return(seq_len(n))
  packed <- .emc_shear_context(ch)$packed
  spread <- .emc_shear_orbit_spread(ch, idx, ks, list(function(st) packed), n_cores)[, 1L]
  which(!is.finite(spread) | spread >= .EMC_SHEAR_ATLAS_MIN_SPREAD)
}

# Add this block's nodes to every chain's atlas: each subject's likelihood
# approximation at its conditional mode under a group state (mean and
# covariance), stored against that state. The approximation is local to the
# mode, whose place on a ridge depends on the group state, so the curves
# regress on the curve coordinates. A mean atlas starts when the pilot shows
# a slow mean direction needs curving (mean = TRUE: the curves are linear in
# every mean direction); its nodes sit at the block's draws, and in the
# discarded stages it starts over when the group log sds drift too far from
# where it started, while the sample stage tops it up to the nodes its terms
# need. A log-sd atlas follows the group log sds whose rounds are measured
# ineffective: its nodes sit on each such coordinate's path from the current
# state, where the orbit travels, at the offsets .EMC_SHEAR_ATLAS_PATH (in
# window sds), so every other group coordinate is held and the regression
# sees only the curve's own coordinates; they are fitted afresh at each
# rebuild. Nodes go only to the subjects whose orbits the curve coordinates
# move (.emc_shear_atlas_subjects).
.emc_shear_atlas_add <- function(emc, idx, n_cores, stage) {
  chains <- integer(0)
  mean_of <- logical(length(emc))
  for (j in seq_along(emc)) {
    ch <- emc[[j]]
    if (!is.list(ch$shear$plan)) next
    atlas <- ch$shear$atlas
    # An atlas without stored log sds only ever followed the mean.
    if (!is.null(atlas) && is.null(atlas$mean)) {
      ch$shear$atlas <- atlas <- NULL
      mean_of[j] <- TRUE
    }
    if (!mean_of[j] && isTRUE(atlas$mean)) mean_of[j] <- TRUE
    if (!mean_of[j]) {
      w <- .emc_shear_curve_wanted(ch, idx, n_cores)
      ch <- .emc_shear_pilot_tick(ch, w$ran)
      mean_of[j] <- isTRUE(w$wanted)
    }
    emc[[j]] <- ch
    # Nodes are fitted only while a curve wants them; a mean atlas no curve
    # uses this block keeps its nodes for when one does.
    if (mean_of[j] || length(.emc_shear_curve_lsd(ch, idx, stage))) chains <- c(chains, j)
  }
  if (!length(chains)) return(emc)
  jobs <- list()
  for (j in chains) {
    ch <- emc[[j]]
    p <- ch$n_pars
    i <- ch$samples$idx
    lsd <- rowMeans(.emc_shear_window_lsd(ch, idx))
    atlas <- ch$shear$atlas
    if (is.null(atlas) || !mean_of[j] ||
        (!identical(stage, "sample") && max(abs(lsd - atlas$lsd_ref)) > .EMC_SHEAR_ATLAS_REANCHOR))
      atlas <- list(lsd_ref = lsd, nodes = vector("list", ch$n_subjects))
    atlas$mean <- mean_of[j]
    q <- if (atlas$mean) .emc_shear_curve_dirs(ch, idx) else integer(0)
    l <- .emc_shear_curve_lsd(ch, idx, stage)
    subs <- .emc_shear_atlas_subjects(ch, idx, c(q, p + l), n_cores)
    ch$shear$atlas <- atlas
    emc[[j]] <- ch
    if (!length(subs)) next
    if (atlas$mean) {
      have <- stats::median(vapply(atlas$nodes[subs], function(a)
        if (is.null(a)) 0L else ncol(a$mu), integer(1)))
      need <- .EMC_SHEAR_ATLAS_NODES
      if (identical(stage, "sample")) {
        nf <- 1L + p + length(l) + length(q) * (length(q) + 1L) / 2
        need <- max(need, min(.EMC_SHEAR_ATLAS_TOPUP, .EMC_SHEAR_ATLAS_PER_COEF * nf - have))
      }
      take <- unique(round(seq(min(idx), max(idx), length.out = need)))
      states <- lapply(take, function(d) list(mu = as.numeric(ch$samples$theta_mu[, d]),
                                              V = ch$samples$theta_var[, , d], draw = d))
    } else {
      mu <- as.numeric(ch$samples$theta_mu[, i])
      V <- ch$samples$theta_var[, , i]
      states <- list(list(mu = mu, V = V, draw = i))
      for (k in p + l) {
        path <- .emc_shear_coord_path(ch, idx, k, mu, V, .EMC_SHEAR_ATLAS_PATH)
        if (!is.null(path)) states <- c(states, lapply(path, function(st)
          list(mu = st$mu, V = st$V, draw = i)))
      }
    }
    for (st in states) for (s in subs)
      jobs[[length(jobs) + 1L]] <- list(chain = j, subject = s, state = st)
  }
  if (!length(jobs)) return(emc)
  emc <- .emc_prepare_chain_rng(emc)
  cost <- .emc_subject_cost(emc[[1L]]$data)
  job_cost <- vapply(jobs, function(jb) cost[jb$subject], numeric(1))
  fits <- .emc_with_preserved_rng(.emc_rebuild_lapply(seq_along(jobs), function(k) {
    jb <- jobs[[k]]
    ch <- emc[[jb$chain]]
    .emc_quiet_fork(.emc_shear_node_fit(
      ch, jb$subject, jb$state$mu, jb$state$V,
      as.numeric(ch$samples$alpha[, jb$subject, jb$state$draw])))
  }, n_cores = n_cores, order = order(job_cost, decreasing = TRUE)))
  for (j in chains) {
    nodes <- emc[[j]]$shear$atlas$nodes
    for (k in which(vapply(jobs, function(jb) jb$chain == j, logical(1)))) {
      f <- fits[[k]]
      if (is.null(f)) next
      s <- jobs[[k]]$subject
      st <- jobs[[k]]$state
      a <- nodes[[s]]
      nodes[[s]] <- list(mu = cbind(a$mu, st$mu),
                         lsd = cbind(a$lsd, 0.5 * log(diag(as.matrix(st$V)))),
                         lin = cbind(a$lin, f$lin), prec = cbind(a$prec, as.vector(f$prec)))
    }
    emc[[j]]$shear$atlas$nodes <- nodes
  }
  emc
}

# Every chain's curved frames from its atlas: per subject the regression of
# (lin, prec) on the features of the group state's curve coordinates about
# the current state (.emc_shear_curve_coords: the plan's mean directions when
# the atlas is a mean one, then the slow log sds standardised by their
# window sds), quadratic in the slow mean directions once the atlas holds
# enough nodes for it, with the quadratic coefficients penalised. Subjects
# without enough nodes keep the block's likelihood fit.
.emc_shear_build_curves <- function(emc, idx, stage) {
  for (j in seq_along(emc)) {
    ch <- emc[[j]]
    ch$shear$curve <- NULL
    plan <- ch$shear$plan
    atlas <- ch$shear$atlas
    if (!is.list(plan) ||
        !is.matrix(plan$dirs) || is.null(atlas$nodes)) {
      emc[[j]] <- ch
      next
    }
    p <- ch$n_pars
    i <- ch$samples$idx
    mu0 <- as.numeric(ch$samples$theta_mu[, i])
    T <- if (isTRUE(atlas$mean)) tryCatch(solve(plan$dirs), error = function(e) NULL) else
      matrix(0, 0L, p)
    l <- .emc_shear_curve_lsd(ch, idx, stage)
    if (is.null(T) || (nrow(T) == 0L && !length(l))) {
      emc[[j]] <- ch
      next
    }
    lsd0 <- 0.5 * log(diag(as.matrix(ch$samples$theta_var[, , i])))[l]
    lsd_scale <- apply(.emc_shear_window_lsd(ch, idx)[l, , drop = FALSE], 1L, stats::sd)
    if (any(!is.finite(lsd_scale) | lsd_scale <= 0)) {
      emc[[j]] <- ch
      next
    }
    q <- if (nrow(T)) .emc_shear_curve_dirs(ch, idx) else integer(0)
    cv0 <- list(mu0 = mu0, T = T, l = l, lsd0 = lsd0, lsd_scale = lsd_scale,
                lsd_range = max(abs(.EMC_SHEAR_ATLAS_PATH)))
    n_lin <- 1L + nrow(T) + length(l)
    n_quad <- length(q) * (length(q) + 1L) / 2
    B <- lapply(atlas$nodes, function(a) {
      if (is.null(a) || is.null(a$lsd)) return(NULL)
      t <- vapply(seq_len(ncol(a$mu)), function(k)
        .emc_shear_curve_coords(cv0, a$mu[, k], a$lsd[, k]), numeric(n_lin - 1L))
      t <- matrix(t, n_lin - 1L)
      keep <- which(sqrt(colSums(t^2)) <= .EMC_SHEAR_ATLAS_RADIUS)
      quad <- length(keep) >= .EMC_SHEAR_ATLAS_PER_COEF * (n_lin + n_quad)
      if (!quad && length(keep) < .EMC_SHEAR_ATLAS_PER_COEF * n_lin) return(NULL)
      qq <- if (quad) q else integer(0)
      X <- t(vapply(keep, function(k) .emc_shear_atlas_features(t[, k], qq),
                    numeric(n_lin + if (quad) n_quad else 0L)))
      Y <- rbind(a$lin[, keep, drop = FALSE], a$prec[, keep, drop = FALSE])
      pen <- diag(c(rep(0, n_lin), rep(.EMC_SHEAR_ATLAS_PENALTY, ncol(X) - n_lin)), ncol(X))
      G <- tryCatch(solve(crossprod(X) + pen), error = function(e) NULL)
      if (is.null(G)) return(NULL)
      coef <- G %*% crossprod(X, t(Y))
      if (!all(is.finite(coef))) return(NULL)
      if (!.emc_shear_curve_predicts(X, t(Y), coef, G)) return(NULL)
      list(coef = coef, quad = quad)
    })
    B <- .emc_shear_curve_proper(ch, idx, cv0, q, B)
    if (any(!vapply(B, is.null, logical(1))))
      ch$shear$curve <- c(cv0, list(q = q, B = B, sets = .emc_shear_curve_sets(B, q)))
    emc[[j]] <- ch
  }
  emc
}

# The curves whose frames exist along their own coordinates: a subject's
# regressed precision must keep its frame precision positive definite at the
# group states the orbit visits along each curve coordinate (.EMC_SHEAR_
# PILOT_GRID), or every candidate there would fail for all subjects; such a
# subject keeps the block's fit.
.emc_shear_curve_proper <- function(ch, idx, cv0, q, B) {
  p <- ch$n_pars
  i <- ch$samples$idx
  mu <- as.numeric(ch$samples$theta_mu[, i])
  V <- ch$samples$theta_var[, , i]
  states <- list()
  for (k in c(q, p + cv0$l)) {
    path <- .emc_shear_coord_path(ch, idx, k, mu, V, .EMC_SHEAR_PILOT_GRID)
    if (!is.null(path)) states <- c(states, path)
  }
  if (!length(states)) return(B)
  for (st in states) {
    Vi <- tryCatch(chol2inv(chol(st$V)), error = function(e) NULL)
    if (is.null(Vi)) next
    t <- .emc_shear_curve_coords(cv0, st$mu, 0.5 * log(diag(as.matrix(st$V))))
    for (s in which(!vapply(B, is.null, logical(1)))) {
      y <- drop(.emc_shear_atlas_features(t, if (B[[s]]$quad) q else integer(0)) %*% B[[s]]$coef)
      prec <- matrix(y[-seq_len(p)], p, p)
      ok <- !inherits(tryCatch(chol((prec + t(prec)) / 2 + Vi), error = function(e) e), "error")
      if (!ok) B[s] <- list(NULL)
    }
  }
  B
}

# Whether a subject's regressed fits predict its nodes better than the
# block's fixed fit would: leave-one-out residuals of the ridge regression
# against those of the nodes' mean, each target standardised, so a subject
# whose node fits are inconsistent (an unstable likelihood mode) keeps its
# fixed frame. Either frame is a function of mu, so the choice is exact.
.emc_shear_curve_predicts <- function(X, Y, coef, G) {
  n <- nrow(X)
  if (n < 3L) return(FALSE)
  sc <- apply(Y, 2, stats::sd)
  use <- which(is.finite(sc) & sc > 0)
  if (!length(use)) return(FALSE)
  Ys <- sweep(Y[, use, drop = FALSE], 2, sc[use], "/")
  h <- rowSums((X %*% G) * X)
  fit <- sweep(X %*% coef[, use, drop = FALSE], 2, sc[use], "/")
  r_curve <- (Ys - fit) / pmax(1 - h, 1e-8)
  r_const <- sweep(Ys, 2, colMeans(Ys)) * n / (n - 1)
  sum(r_curve^2) < sum(r_const^2)
}

# A group state's coordinates in a curve: the mean directions' coordinates
# of mu - mu0, then the curve's log-sd coordinates l standardised about lsd0
# and held within lsd_range, the span of the atlas's path nodes: beyond its
# last nodes a curve keeps the frames it has there rather than extrapolate.
.emc_shear_curve_coords <- function(cv, mu, lsd) {
  t <- drop(cv$T %*% (mu - cv$mu0))
  if (length(cv$l)) {
    r <- if (is.null(cv$lsd_range)) Inf else cv$lsd_range
    t <- c(t, pmin(r, pmax(-r, (lsd[cv$l] - cv$lsd0) / cv$lsd_scale)))
  }
  t
}

# Each subject's likelihood approximation at group state (mu, V), packed for
# the frames: the block's fits, or the curved frames where it has them (one
# product with the stacked coefficients). An interpolated precision that
# leaves the conditional precision indefinite fails that state's frames, so
# the state is invalid (a candidate there is rejected): the frames stay a
# function of the group state.
.emc_shear_fits_at <- function(ctx, mu, V) {
  cv <- ctx$curve
  if (is.null(cv)) return(ctx$packed)
  t <- .emc_shear_curve_coords(cv, mu, 0.5 * log(diag(as.matrix(V))))
  p <- length(mu)
  out <- ctx$packed
  for (set in cv$sets) {
    y <- matrix(drop(.emc_shear_atlas_features(t, set$q) %*% set$coef), p + p * p)
    out$lin[, set$idx] <- y[seq_len(p), , drop = FALSE]
    prec <- array(y[-seq_len(p), , drop = FALSE], c(p, p, length(set$idx)))
    out$cube[, , set$idx] <- (prec + aperm(prec, c(2, 1, 3))) / 2
  }
  out
}

# The curved subjects' coefficients stacked by feature set (quadratic or
# linear) for one product per set.
.emc_shear_curve_sets <- function(B, q) {
  sets <- list()
  for (quad in c(TRUE, FALSE)) {
    idx <- which(vapply(B, function(b) !is.null(b) && identical(b$quad, quad), logical(1)))
    if (!length(idx)) next
    sets[[length(sets) + 1L]] <- list(idx = idx, q = if (quad) q else integer(0),
                                      coef = do.call(cbind, lapply(B[idx], `[[`, "coef")))
  }
  sets
}

# A block's round record: per coordinate the selection-weighted squared
# moves E over n rounds (moved of them moving), the net squared move N of
# its rounds per iteration with K the iterations' round counts, and the
# iterations the record spans (iters, zero where unknown). Fills in the
# fields a record made without them lacks.
.emc_shear_record <- function(rec, p) {
  if (is.null(rec))
    rec <- list(E = numeric(2L * p), n = integer(2L * p), moved = integer(2L * p))
  if (is.null(rec$N)) rec$N <- numeric(2L * p)
  if (is.null(rec$K)) rec$K <- matrix(0L, 2L * p, 1L)
  if (is.null(rec$iters)) rec$iters <- numeric(2L * p)
  rec
}

# A production rebuild moves each step scale towards the target share of
# moving rounds, from the last block's rounds `lst` (rounds within a block
# never change it).
.emc_shear_tune_steps <- function(plan, lst) {
  ok <- lst$n >= 10L
  if (any(ok)) {
    rate <- lst$moved[ok] / lst$n[ok]
    plan$log_step[ok] <- plan$log_step[ok] +
      pmax(-log(2), pmin(log(2), 2 * (rate - .EMC_SHEAR_MOVE_TARGET)))
  }
  plan$log_step <- pmin(.EMC_SHEAR_LOG_STEP_RANGE[2],
                        pmax(.EMC_SHEAR_LOG_STEP_RANGE[1], plan$log_step))
  plan
}

# Covariance prior as each variant's Gibbs step targets it: diagonal-gamma
# has a gamma prior on every precision; standard has an inverse-Wishart per
# blocked parameter group and a gamma prior on the precision of every
# unblocked parameter, given a_half, which is integrated out (see the header).
.emc_shear_v_spec <- function(sampler) {
  p <- sampler$n_pars
  prior <- sampler$prior
  if (identical(sampler$type, "diagonal-gamma")) {
    return(list(iw = list(), gamma_idx = seq_len(p), size = rep(1L, p),
                shape = rep_len(prior$shape, p), rate = rep_len(prior$rate, p),
                uses_a_half = FALSE))
  }
  is_blocked <- sampler$is_blocked[!sampler$nuisance]
  par_group <- sampler$par_group[!sampler$nuisance]
  if (is.null(is_blocked) || length(is_blocked) == 0L) is_blocked <- rep(TRUE, p)
  if (is.null(par_group)) par_group <- rep(1, p)
  blocked <- which(is_blocked)
  iw <- lapply(unique(par_group[blocked]), function(g) blocked[par_group[blocked] == g])
  size <- rep(1L, p)
  for (b in iw) size[b] <- length(b)
  list(iw = iw, gamma_idx = which(!is_blocked), size = size,
       v = prior$v, A = rep_len(prior$A, p), uses_a_half = TRUE)
}

# Log density of V over the free entries of its blocks.
.emc_shear_log_prior_v <- function(V, spec) {
  lp <- 0
  g <- spec$gamma_idx
  if (length(g)) {
    tau <- 1 / V[cbind(g, g)]
    # Density of the variance: the precision density times tau^2.
    lp <- lp + if (spec$uses_a_half) {
      sum((spec$v / 2 + 1) * log(tau) -
            (spec$v + 1) / 2 * log(spec$v * tau + 1 / spec$A[g]^2))
    } else {
      sum((spec$shape[g] + 1) * log(tau) - spec$rate[g] * tau)
    }
  }
  for (b in spec$iw) {
    d <- length(b)
    R <- chol(V[b, b, drop = FALSE])
    df <- spec$v + d - 1
    lp <- lp - (df + d + 1) * sum(log(diag(R))) -
      (spec$v + d) / 2 * sum(log(spec$v * diag(chol2inv(R)) + 1 / spec$A[b]^2))
  }
  lp
}

# a_half given V, as the standard Gibbs step draws it:
# IG((v + d_k) / 2, v (V^-1)_kk + 1 / A_k^2), with d_k the size of k's block.
.emc_shear_draw_a_half <- function(V, spec) {
  prec <- numeric(nrow(V))
  g <- spec$gamma_idx
  if (length(g)) prec[g] <- 1 / V[cbind(g, g)]
  for (b in spec$iw) prec[b] <- diag(chol2inv(chol(V[b, b, drop = FALSE])))
  1 / stats::rgamma(nrow(V), shape = (spec$v + spec$size) / 2,
                    rate = spec$v * prec + 1 / spec$A^2)
}

.emc_shear_log_prior_mu <- function(mu, prior) {
  d <- mu - prior$theta_mu_mean
  -0.5 * drop(crossprod(d, prior$theta_mu_invar %*% d))
}

# Everything the move needs about one group state (mu, V): its prior density
# and each subject's conditional frame (m_s, R_s), with R_s the upper Cholesky
# factor of the conditional precision H_s + V^-1, stacked by subject (fm:
# m by column, R and logdet by slice; .emc_shear_frame takes out one).
.emc_shear_state <- function(mu, V, fits, prior, spec) {
  chol_V <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(chol_V)) return(NULL)
  Vi <- chol2inv(chol_V)
  if (is.null(fits$cube)) fits <- .emc_shear_pack(fits, length(mu))
  out <- emc_shear_frames_cpp(fits$cube, fits$lin, Vi, drop(Vi %*% mu))
  if (!isTRUE(out$ok)) return(NULL)
  .emc_shear_state_of(mu, V, Vi, chol_V, list(m = out$m, R = out$R, logdet = out$logdet),
                      prior, spec)
}

.emc_shear_state_of <- function(mu, V, Vi, chol_V, fm, prior, spec) {
  lp_v <- tryCatch(.emc_shear_log_prior_v(V, spec), error = function(e) NA_real_)
  list(mu = mu, V = V, Vi = Vi, chol_V = chol_V, fm = fm,
       logdet = sum(fm$logdet), lp_v = lp_v,
       lp = .emc_shear_log_prior_mu(mu, prior) + lp_v)
}

# Subject s's frame of a state.
.emc_shear_frame <- function(state, s) {
  list(m = state$fm$m[, s], R = state$fm$R[, , s], logdet = state$fm$logdet[s])
}

# The likelihood approximations as the frames' C++ takes them: a cube of
# precisions (a zero slice for a subject without one) and their linear
# terms by column.
.emc_shear_pack <- function(fits, p) {
  n <- length(fits)
  cube <- array(0, c(p, p, n))
  lin <- matrix(0, p, n)
  for (s in seq_len(n)) {
    fit <- fits[[s]]
    if (is.null(fit)) next
    cube[, , s] <- fit$prec
    lin[, s] <- fit$lin
  }
  list(cube = cube, lin = lin)
}

# The state after a mean step u_mu: V, and so every R_s, is unchanged and
# each centre moves by Q_s^-1 V^-1 u_mu.
.emc_shear_shift_state <- function(cur, u_mu, prior) {
  new <- cur
  new$fm$m <- cur$fm$m + .emc_shear_centre_shift(cur, u_mu)
  new$mu <- cur$mu + u_mu
  new$lp <- .emc_shear_log_prior_mu(new$mu, prior) + cur$lp_v
  new
}

# sum_s log N(a_s | mu, V).
.emc_shear_log_hier <- function(alpha, state) {
  z <- backsolve(state$chol_V, alpha - state$mu, transpose = TRUE)
  -0.5 * sum(z^2) - ncol(alpha) * sum(log(diag(state$chol_V)))
}

# Each frame centre's move under a mean step u_mu with V fixed:
# Q_s^-1 V^-1 u_mu (p x n).
.emc_shear_centre_shift <- function(cur, u_mu) {
  Vi_u <- drop(cur$Vi %*% u_mu)
  R <- cur$fm$R
  matrix(vapply(seq_len(dim(R)[3]), function(s)
    drop(backsolve(R[, , s], forwardsolve(t(R[, , s]), Vi_u))), numeric(length(u_mu))),
    length(u_mu))
}

# The map for a given u, from the group state cur and subjects alpha: the
# moved state, the moved subjects and log|det| of the map's Jacobian. NULL
# when the moved state is invalid.
.emc_shear_propose <- function(cur, u, alpha, fits, prior, spec, curve = NULL,
                               walls = NULL) {
  p <- length(cur$mu)
  u_sd <- u[p + seq_len(p)]
  if (all(u_sd == 0) && is.null(curve) && is.null(walls)) {
    new <- .emc_shear_shift_state(cur, u[seq_len(p)], prior)
    if (!is.finite(new$lp)) return(NULL)
    return(list(state = new, alpha = alpha + new$fm$m - cur$fm$m, log_jac = 0))
  }
  .emc_shear_candidates(cur, u, 1, alpha, fits, prior, spec, curve, walls)[[1L]]
}

# The candidates of one round: for each offset t[k] the map for u = t[k] v
# (as .emc_shear_propose), built in one call (emc_shear_candidates_cpp).
# Per candidate NULL when the moved state is invalid, else the moved state,
# the moved subjects, their summed hierarchical density (hier, as
# .emc_shear_log_hier), log|det| of the map's Jacobian and the walls' part
# of it.
.emc_shear_candidates <- function(cur, v, t, alpha, fits, prior, spec, curve = NULL,
                                  walls = NULL) {
  p <- length(cur$mu)
  if (is.null(fits$cube)) fits <- .emc_shear_pack(fits, p)
  n <- ncol(alpha)
  if (is.null(walls)) {
    has <- logical(0); U <- matrix(0, p, 0); b <- numeric(0)
    Uc <- array(0, c(p, max(p - 1L, 0L), 0))
  } else {
    has <- walls$has; U <- walls$U; b <- ifelse(has, walls$b, 0)
    Uc <- .emc_shear_uperp_cube(walls, p)
  }
  f0 <- .emc_shear_frame_arrays(cur)
  out <- emc_shear_candidates_cpp(t, v, alpha, cur$mu, cur$V, f0$m, f0$R, f0$logdet,
                                  fits$cube, fits$lin, curve, has, U, b, Uc)
  u_sd <- v[p + seq_len(p)]
  lapply(seq_along(t), function(k) {
    o <- out[[k]]
    if (!isTRUE(o$ok)) return(NULL)
    u <- t[k] * v
    d <- exp(u[p + seq_len(p)])
    new <- .emc_shear_state_of(cur$mu + u[seq_len(p)], cur$V * outer(d, d), o$Vi, o$chol_V,
                               list(m = o$m, R = o$R, logdet = o$logdet), prior, spec)
    if (!is.finite(new$lp)) return(NULL)
    list(state = new, alpha = o$alpha, hier = o$hier,
         log_jac = sum((spec$size + 1) * t[k] * u_sd) + o$jac, wall_jac = o$wall_jac)
  })
}

# Every subject transported from group state cur's frames to new's
# (emc_shear_transport_cpp, the C++ form of the frame map and
# .emc_shear_wall_map).
.emc_shear_transport <- function(alpha, cur, new, walls = NULL) {
  f0 <- .emc_shear_frame_arrays(cur)
  f1 <- .emc_shear_frame_arrays(new)
  p <- nrow(alpha)
  if (is.null(walls)) {
    has <- logical(0)
    U <- matrix(0, p, 0)
    b <- numeric(0)
    Uc <- array(0, c(p, max(p - 1L, 0L), 0))
  } else {
    has <- walls$has
    U <- walls$U
    b <- ifelse(has, walls$b, 0)
    Uc <- .emc_shear_uperp_cube(walls, p)
  }
  emc_shear_transport_cpp(alpha, f0$m, f0$R, f0$logdet, f1$m, f1$R, f1$logdet,
                          has, U, b, Uc)
}

# A state's frames stacked as the transport kernel takes them.
.emc_shear_frame_arrays <- function(state) state$fm

# The walls' complement bases as a cube (a zero slice without a wall).
.emc_shear_uperp_cube <- function(walls, p) {
  if (!is.null(walls$Uc)) return(walls$Uc)
  n <- length(walls$has)
  Uc <- array(0, c(p, p - 1L, n))
  for (s in which(walls$has)) Uc[, , s] <- walls$Uperp[[s]]
  Uc
}

# What every round of one move shares: the usable fits and the covariance
# prior.
.emc_shear_context <- function(sampler) {
  cv <- sampler$shear$curve
  ok <- is.list(cv) && is.matrix(cv$T) && ncol(cv$T) == sampler$n_pars &&
    length(cv$B) == sampler$n_subjects && is.list(cv$sets)
  fits <- .emc_shear_fits(sampler)
  list(fits = fits, packed = .emc_shear_pack(fits, sampler$n_pars),
       spec = .emc_shear_v_spec(sampler), curve = if (ok) cv else NULL,
       walls = .emc_shear_walls_of(sampler))
}

# The block's transport walls when they match the sampler's dimensions.
.emc_shear_walls_of <- function(sampler) {
  w <- sampler$shear$walls
  ok <- is.list(w) && is.matrix(w$U) && nrow(w$U) == sampler$n_pars &&
    ncol(w$U) == sampler$n_subjects && length(w$has) == sampler$n_subjects &&
    any(w$has)
  if (!ok) return(NULL)
  w$Uc <- .emc_shear_uperp_cube(w, sampler$n_pars)
  w
}

# A mean step's candidates lie on a line: with Delta_s = Q_s^-1 V^-1 d, the
# subjects move to a_s + t Delta_s and mu to mu + t d, so the hierarchical
# and mean-prior terms are quadratics in t. Returns, per candidate, the moved
# subjects and the change of those terms (V, and so the Jacobian, fixed).
.emc_shear_line <- function(cur, d, alpha, t, prior) {
  Delta <- .emc_shear_centre_shift(cur, d)
  Z0 <- backsolve(cur$chol_V, alpha - cur$mu, transpose = TRUE)
  Zd <- backsolve(cur$chol_V, Delta - d, transpose = TRUE)
  P0 <- as.matrix(prior$theta_mu_invar)
  m0 <- cur$mu - prior$theta_mu_mean
  a1 <- sum(Z0 * Zd) + sum(m0 * (P0 %*% d))
  a2 <- sum(Zd^2) + sum(d * (P0 %*% d))
  list(Delta = Delta, delta = -(a1 * t + 0.5 * a2 * t^2))
}

# ---- Wall-aware transport --------------------------------------------------
# Quantile map between standard normals truncated below at c0 and at c1:
# the point z (z >= c0) keeps its conditional upper-tail probability. Below
# the wall (z < c0) the mirror map between the normals truncated above at c0
# and at c1 applies, so each side maps onto itself and the map is a
# bijection of the line.
.emc_shear_qmap <- function(z, c0, c1) {
  if (is.finite(z) && z < c0) {
    m <- .emc_shear_qmap(-z, -c0, -c1)
    return(list(z = -m$z, log_jac = m$log_jac))
  }
  ls0 <- stats::pnorm(c0, lower.tail = FALSE, log.p = TRUE)
  ls1 <- stats::pnorm(c1, lower.tail = FALSE, log.p = TRUE)
  lf <- stats::pnorm(z, lower.tail = FALSE, log.p = TRUE) - ls0
  z1 <- stats::qnorm(lf + ls1, lower.tail = FALSE, log.p = TRUE)
  list(z = z1, log_jac = stats::dnorm(z, log = TRUE) - stats::dnorm(z1, log = TRUE) - ls0 + ls1)
}

# Transport of one subject from frame f0 to frame f1 (centres m, upper
# Cholesky factors R of the precisions) that keeps the half-space u'a >= b:
# the complement of u moves by the affine map between the frames' marginals
# there, and u'a by the quantile map between the frames' conditionals given
# the complement (Gaussian, precision u'Qu). The maps compose, and where b
# is far below the map is affine: the frame map whose whitening takes the
# wall's normal last. A state outside its wall takes the mirror quantile map
# and stays outside, so the map is a bijection whichever side the state is
# on. Returns the moved subject and the log Jacobian, or NULL when the map
# cannot be computed (the candidate is then invalid).
.emc_shear_wall_map <- function(a, f0, f1, u, b, Uperp) {
  cond <- function(f, y) {
    Qu <- drop(crossprod(f$R, f$R %*% u))
    q <- sum(u * Qu)
    # centre of u'a given the complement coordinates y (relative to m)
    list(mc = sum(u * f$m) - sum(Qu * drop(Uperp %*% y)) / q, sc = 1 / sqrt(q))
  }
  y0 <- drop(crossprod(Uperp, a - f0$m))
  S0 <- crossprod(Uperp, chol2inv(f0$R) %*% Uperp)
  S1 <- crossprod(Uperp, chol2inv(f1$R) %*% Uperp)
  L0 <- tryCatch(t(chol(S0)), error = function(e) NULL)
  L1 <- tryCatch(t(chol(S1)), error = function(e) NULL)
  if (is.null(L0) || is.null(L1)) return(NULL)
  y1 <- drop(L1 %*% forwardsolve(L0, y0))
  c0 <- cond(f0, y0); c1 <- cond(f1, y1)
  t0 <- sum(u * a)
  z0 <- (t0 - c0$mc) / c0$sc
  w0 <- (b - c0$mc) / c0$sc
  if (!is.finite(z0)) return(NULL)
  w1 <- (b - c1$mc) / c1$sc
  qm <- .emc_shear_qmap(z0, w0, w1)
  t1 <- c1$mc + c1$sc * qm$z
  if (!is.finite(t1)) return(NULL)
  list(a = f1$m + drop(Uperp %*% y1) + (t1 - sum(u * f1$m)) * u,
       log_jac = sum(log(diag(L1))) - sum(log(diag(L0))) + qm$log_jac - log(c0$sc) + log(c1$sc))
}

# One round along step direction v with scale step (see the header); x
# carries the current state, its subject likelihoods and the worker pool.
# Adds E, the selection-weighted squared move in t, and whether it moved.
.emc_shear_round <- function(x, v, step, sampler, ctx) {
  # Elapsed time building the candidates and in their likelihood round.
  build_el0 <- proc.time()[["elapsed"]]
  t_aux <- step * stats::rnorm(1L)
  t <- t_aux + step * stats::rnorm(.EMC_SHEAR_CANDIDATES)
  p <- nrow(x$alpha)
  n <- ncol(x$alpha)
  # Mean steps lie on a line when the frames are fixed in mu and no
  # subject has a wall.
  mean_step <- all(v[p + seq_len(p)] == 0) && is.null(ctx$curve) && is.null(ctx$walls)
  if (mean_step) {
    d <- v[seq_len(p)]
    line <- .emc_shear_line(x$cur, d, x$alpha, t, sampler$prior)
    moved_alpha <- function(i) x$alpha + t[i] * line$Delta
    prior_delta <- line$delta
    ok <- which(is.finite(prior_delta))
  } else {
    props <- .emc_shear_candidates(x$cur, v, t, x$alpha, ctx$packed, sampler$prior, ctx$spec,
                                   ctx$curve, ctx$walls)
    ok <- which(!vapply(props, is.null, logical(1)))
    moved_alpha <- function(i) props[[i]]$alpha
    prior_delta <- vapply(seq_along(t), function(i) {
      pr <- props[[i]]
      if (is.null(pr)) return(-Inf)
      pr$hier + pr$state$lp + pr$log_jac - x$cur_hier - x$cur$lp
    }, numeric(1))
  }
  lw <- rep(-Inf, length(t))
  ll_new <- NULL
  if (length(ok)) {
    pn <- sampler$par_names
    cand <- lapply(ok, moved_alpha)
    rows <- lapply(seq_len(n), function(s)
      matrix(vapply(cand, function(a) a[, s], numeric(p)), length(ok), byrow = TRUE,
             dimnames = list(NULL, pn)))
    ll_el0 <- proc.time()[["elapsed"]]
    x$build_el <- x$build_el + ll_el0 - build_el0
    evaluated <- .emc_wpool_subject_ll(x$wpool, rows, sampler, x$wpool_part,
                                      source = "shear")
    x$ll_el <- x$ll_el + proc.time()[["elapsed"]] - ll_el0
    x$wpool <- evaluated$pool
    x$rounds <- x$rounds + 1L
    ll_new <- matrix(evaluated$ll, n, length(ok))
    v_ok <- colSums(ll_new) - sum(x$ll) + prior_delta[ok]
    lw[ok] <- ifelse(is.finite(v_ok), v_ok, -Inf)
  } else x$build_el <- x$build_el + proc.time()[["elapsed"]] - build_el0
  w <- c(0, lw)
  w <- exp(w - max(w))
  w <- w / sum(w)
  j <- sample.int(length(w), 1L, prob = w) - 1L
  x$E <- sum(w[-1L] * t^2)
  x$moved <- j > 0L
  x$t <- if (x$moved) t[j] else 0
  if (x$moved) {
    if (mean_step) {
      x$cur <- .emc_shear_shift_state(x$cur, t[j] * d, sampler$prior)
      x$alpha <- moved_alpha(j)
    } else {
      x$cur <- props[[j]]$state
      x$alpha <- props[[j]]$alpha
    }
    x$cur_hier <- .emc_shear_log_hier(x$alpha, x$cur)
    x$ll <- ll_new[, match(j, ok)]
  }
  x
}

.emc_shear_move <- function(sampler, pars, proposals, stage, wpool = NULL,
                             wpool_part = NULL) {
  keep <- list(sampler = sampler, pars = pars, proposals = proposals, wpool = wpool)
  if (stage == "preburn" || !.emc_shear_active(sampler)) return(keep)
  st <- sampler$shear
  plan <- st$plan
  n <- sampler$n_subjects
  p <- sampler$n_pars
  if (!is.list(plan) || !is.matrix(plan$dirs) || any(dim(plan$dirs) != p) ||
      length(plan$rounds) != 2L * p || !any(plan$rounds > 0)) return(keep)
  if (is.null(st$iterations)) st$iterations <- st$rounds <- 0L
  if (!is.null(sampler$rng$gibbs))
    assign(".Random.seed", sampler$rng$gibbs, envir = globalenv())
  # The allocated expected rounds per coordinate: the whole part, plus one
  # more with the fractional part as probability, in random order.
  whole <- floor(plan$rounds)
  n_rounds <- whole + (stats::runif(2L * p) < plan$rounds - whole)
  take <- rep.int(seq_len(2L * p), n_rounds)
  st$iterations <- st$iterations + 1L
  if (!length(take)) {
    sampler$rng$gibbs <- get(".Random.seed", envir = globalenv())
    sampler$shear <- st
    keep$sampler <- sampler
    return(keep)
  }
  ctx <- .emc_shear_context(sampler)
  alpha <- matrix(proposals[seq_len(p), ], p)
  mu <- as.numeric(pars$tmu)
  cur <- .emc_shear_state(mu, as.matrix(pars$tvar), .emc_shear_fits_at(ctx, mu, pars$tvar),
                          sampler$prior, ctx$spec)
  if (is.null(cur) || !is.finite(cur$lp)) {
    sampler$rng$gibbs <- get(".Random.seed", envir = globalenv())
    sampler$shear <- st
    keep$sampler <- sampler
    return(keep)
  }
  x <- list(cur = cur, cur_hier = .emc_shear_log_hier(alpha, cur),
            alpha = alpha, ll = proposals[p + 1L, ], rounds = 0L, build_el = 0, ll_el = 0,
            wpool = wpool, wpool_part = wpool_part)
  adapt <- stage != "sample"
  any_moved <- FALSE
  blk <- .emc_shear_record(plan$block, p)
  net <- numeric(2L * p)
  round_el0 <- proc.time()[["elapsed"]]
  for (k in take[sample.int(length(take))]) {
    x <- .emc_shear_round(x, .emc_shear_coord_dir(plan, k), exp(plan$log_step[k]),
                          sampler, ctx)
    net[k] <- net[k] + x$t
    blk$E[k] <- blk$E[k] + x$E
    blk$n[k] <- blk$n[k] + 1L
    blk$moved[k] <- blk$moved[k] + as.integer(x$moved)
    if (adapt) {
      plan$log_step[k] <- min(.EMC_SHEAR_LOG_STEP_RANGE[2], max(.EMC_SHEAR_LOG_STEP_RANGE[1],
        plan$log_step[k] + (x$moved - .EMC_SHEAR_MOVE_TARGET) / sqrt(blk$n[k])))
    }
    any_moved <- any_moved || x$moved
  }
  round_el <- proc.time()[["elapsed"]] - round_el0
  # Each coordinate's net move over its rounds in this iteration, by how many
  # it had (the rounds' saturation, R/move_allocation.R).
  ran <- which(n_rounds > 0L)
  blk$N[ran] <- blk$N[ran] + net[ran]^2
  for (kk in unique(n_rounds[ran]))
    blk$K <- .emc_alloc_khist_add(blk$K, ran[n_rounds[ran] == kk], kk)
  plan$block <- blk
  st$plan <- plan
  st$rounds <- st$rounds + x$rounds
  if (ctx$spec$uses_a_half)
    pars$a_half <- stats::setNames(.emc_shear_draw_a_half(x$cur$V, ctx$spec),
                                   names(pars$a_half))
  sampler$rng$gibbs <- get(".Random.seed", envir = globalenv())
  if (any_moved) {
    pars$tmu <- stats::setNames(x$cur$mu, names(pars$tmu))
    pars$tvar <- `dimnames<-`(x$cur$V, dimnames(pars$tvar))
    pars$tvinv <- `dimnames<-`(x$cur$Vi, dimnames(pars$tvar))
    # Every subject's prior mean is the group mean here (see
    # .emc_shear_supported); updates later in the iteration read it.
    if (is.matrix(pars$subj_mu)) pars$subj_mu[] <- x$cur$mu
    proposals[seq_len(p), ] <- x$alpha
    proposals[p + 1L, ] <- x$ll
    pars$alpha <- x$alpha
  }
  sampler$shear <- st
  list(sampler = sampler, pars = pars, proposals = proposals, wpool = x$wpool,
       rounds = x$rounds, round_el = round_el, build_el = x$build_el, ll_el = x$ll_el)
}
