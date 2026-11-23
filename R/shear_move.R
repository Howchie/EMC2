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
# compose and invert, and the Jacobian is exact, so the rounds stay valid.
# Far from the wall the map is an affine frame map with a different whitening
# order from the unbounded transport. Assigning or dropping a wall at a
# rebuild can therefore change a subject's orbit. The collapsed move uses the
# same walls.
#
# Curved frames. The frames may be any function of the group state: the
# map, its composition and its Jacobian hold whatever likelihood
# approximation each state uses. A fit at the likelihood's mode carries a
# subject along a straight line, and where the subject's posterior ridge
# bends (a shift traded against a rate and a threshold) the candidates leave
# it, so the orbit along the group coordinate is a fraction of the
# posterior's width. For one mean direction k that mixes slowly, each
# subject's likelihood is refitted at its conditional mode at a few knots of
# the coordinate c = (D^-1 mu)_k, and the fits (the precision and linear
# term) are interpolated linearly in c and held at the end knots beyond them
# (.emc_shear_curve_subject, emc_shear_curve_fits_cpp). A round whose
# direction leaves c alone runs as before, on the fits at the current c; one
# that moves c rebuilds the frames per candidate.
#
# A round resamples its coordinate within the orbit only, so with an orbit of
# r posterior sds of the coordinate it refreshes at most the share r^2 of
# the coordinate's variance, and a round at the step target
# .EMC_SHEAR_ROUND_REFRESH of that: its factor on the coordinate's
# autocorrelation is 1 - .EMC_SHEAR_ROUND_REFRESH r^2
# (.emc_shear_round_factor). A curve is worth its fits where it widens a
# narrow orbit, not where the fixed fits' orbit already spans the margin. At
# a rebuild a coordinate that is slow over the window is a candidate
# whatever rounds it has had; its curve is kept only if a round under it is
# worth several under the fixed fits (.emc_shear_curve_worth), and reused
# over later rebuilds while it covers the window (.emc_shear_set_curves).
# Where no coordinate is slow the curves cost the window's autocorrelations.
#
# Supported hierarchies: those in which every subject's prior mean is the
# group mean (standard without group designs, diagonal-gamma); the others
# keep the ensemble update alone.

# Candidates per round, and the share of rounds that should move (the
# optimum of the step scale c for N = 8 on a Gaussian orbit is flat around
# this value).
.EMC_SHEAR_CANDIDATES <- 8L
.EMC_SHEAR_MOVE_TARGET <- 0.68
# Initial step scale in posterior sds of the coordinate, the rounds whose
# move share weighs as much as a measured orbit width, the grid (in widths)
# a width is measured on, and the bounds every adapted scale stays in.
.EMC_SHEAR_STEP0 <- 2.5
.EMC_SHEAR_STEP_ROUNDS <- 20
.EMC_SHEAR_WIDTH_GRID <- seq(-3, 3, by = 0.5)
.EMC_SHEAR_LOG_STEP_RANGE <- log(c(0.05, 20))
# The share of a Gaussian orbit's variance one round refreshes at a step of
# .EMC_SHEAR_STEP0 orbit sds (the lag-one autocorrelation of the N = 8 round
# on a standard normal is one minus it).
.EMC_SHEAR_ROUND_REFRESH <- 0.66
# Curved frames: the integrated autocorrelation time over the window from
# which a mean coordinate counts as slow (a thousand iterations then hold
# fewer than a hundred effective draws of it) and the lags it is summed
# over, the knots (in posterior sds of the coordinate around the window's
# mean), and the fixed-frame rounds one round under a curve must be worth
# for the curve to be kept.
.EMC_SHEAR_CURVE_TAU <- 10
.EMC_SHEAR_CURVE_LAGS <- 249L
.EMC_SHEAR_CURVE_KNOTS <- c(-3, -1.5, 0, 1.5, 3)
.EMC_SHEAR_CURVE_WORTH <- 2
# Per chain and rebuild, the coordinates whose fixed-frame orbit may be
# measured and those whose curve may be fitted; and the cosine from which two
# coordinates count as the same one.
.EMC_SHEAR_CURVE_SCREENS <- 4L
.EMC_SHEAR_CURVE_TRIES <- 2L
.EMC_SHEAR_CURVE_SAME <- 0.8

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
      found <- if (getOption("emc2.native_shear_support", TRUE) && is.double(lim) &&
                   (all(!is.na(cols)) || is.numeric(consts)))
        emc_shear_half_spaces_cpp(
          X, cols, if (any(is.na(cols))) as.numeric(consts[colnames(X)[is.na(cols)]]) else numeric(0),
          lim, p)
      if (!is.null(found)) {
        out <- c(out, found)
        next
      }
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
# directions (.emc_shear_carry_weights; a log group sd always matches
# itself); the round counts are set by the move allocation afterwards.
.emc_shear_set_plans <- function(emc, samples_idx, stage, n_cores = 1L) {
  if (!.emc_shear_active(emc[[1L]])) return(emc)
  idx <- .emc_chain_proposal_window(emc, samples_idx, stage)
  if (!length(idx)) return(emc)
  for (j in seq_along(emc)) emc[[j]] <- .emc_shear_set_plan(emc[[j]], idx, stage)
  emc <- .emc_shear_set_curves(emc, idx, n_cores)
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
  if (any(vapply(path, function(st) is.null(st$chol), logical(1))))
    return(list(per = matrix(NA_real_, length(grid), n), lp = rep(NA_real_, length(grid))))
  states <- lapply(path, function(st)
    .emc_shear_state(st$mu, st$V, fits_of(st), ch$prior, spec))
  if (any(vapply(states, is.null, logical(1))))
    return(list(per = matrix(NA_real_, length(grid), n), lp = rep(NA_real_, length(grid))))
  tr <- lapply(states, function(st) .emc_shear_transport(alpha, states[[g0]], st, walls))
  ok <- Reduce(`&`, lapply(tr, function(x) x$ok %in% TRUE))
  # Hierarchical density of every subject at one path state.
  ldn <- vapply(seq_along(grid), function(g) {
    st <- path[[g]]
    z <- backsolve(st$chol, tr[[g]]$alpha - st$mu, transpose = TRUE)
    -0.5 * colSums(z^2) - sum(log(diag(st$chol)))
  }, numeric(n))
  validated <- isTRUE(ch$ll_cache_validated)
  per <- .emc_with_preserved_rng(auto_mclapply(seq_len(n), function(ss) {
    if (!ok[ss]) return(rep(NA_real_, length(grid)))
    X <- t(vapply(tr, function(x) x$alpha[, ss], numeric(p)))
    ll <- tryCatch(.emc_burn_ll_batch(X, ch$data[[ss]], ch$model, par_names,
                                      ll_cache_validated = validated),
                   error = function(e) rep(NA_real_, length(grid)))
    vapply(seq_along(grid), function(g)
      ll[g] + ldn[ss, g] + states[[g0]]$fm$logdet[ss] -
        states[[g]]$fm$logdet[ss] + tr[[g]]$jac[ss], numeric(1))
  }, mc.cores = max(1L, min(n, n_cores)), mc.set.seed = FALSE))
  lp <- vapply(path, function(st) .emc_shear_log_prior_mu(st$mu, ch$prior) +
    (if (k <= p) 0 else tryCatch(.emc_shear_log_prior_v(st$V, spec), error = function(e) NA_real_) +
       sum((spec$size + 1) * st$u_sd)), numeric(1))
  list(per = matrix(unlist(per), length(grid)), lp = lp)
}

# The sd of the orbit's density along coordinate k (in .emc_shear_coord_path's
# units), on a grid scaled to the width s0 and rescaled until the grid
# resolves the orbit; NA when it cannot be measured.
.emc_shear_orbit_width <- function(ch, idx, k, fits_of, s0, n_cores = 1L) {
  if (!is.finite(s0) || s0 <= 0) return(NA_real_)
  s <- s0
  for (it in seq_len(8L)) {
    grid <- s * .EMC_SHEAR_WIDTH_GRID
    tm <- .emc_shear_orbit_terms(ch, idx, k, grid, fits_of, n_cores)
    if (is.null(tm)) return(NA_real_)
    f <- rowSums(tm$per) + tm$lp
    # A failed frame/transport is missing numerical information, not zero
    # posterior density. Retry closer in rather than measuring an artificial
    # truncation (or discarding an otherwise useful coordinate altogether).
    if (anyNA(f) || any(f == Inf)) {
      s <- s / 4
      next
    }
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
      # Nearly flat fitted curvature must not launch a probe tens of log-sd
      # units away. Expand gradually, checking the intervening density.
      s <- if (is.finite(a) && a < 0) min(4 * s, max(2 * s, 1 / sqrt(-2 * a))) else 4 * s
    } else if (sd < 0.5 * s) s <- max(sd, s / 8)
    else return(sd)
  }
  # No resolved width: leave the existing step alone.
  NA_real_
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
  packed <- vector("list", length(emc))
  for (j in unique(vapply(jobs, `[[`, integer(1), "chain"))) {
    packed[j] <- list(.emc_shear_context(emc[[j]])$packed)
    # Data a sampling block has prepared needs no further validation.
    emc[[j]]$ll_cache_validated <- .emc_ll_data_prepared(emc[[j]]$data)
  }
  widths <- .emc_with_preserved_rng(.emc_rebuild_lapply(seq_along(jobs), function(m) {
    jb <- jobs[[m]]
    ch <- emc[[jb$chain]]
    fits <- packed[[jb$chain]]
    tryCatch(.emc_shear_orbit_width(ch, idx, jb$k, function(st) fits, jb$s0),
             error = function(e) NA_real_)
  }, n_cores = n_cores))
  for (j in seq_along(emc)) emc[[j]]$ll_cache_validated <- NULL
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

# The integrated autocorrelation time of one coordinate from its draws in
# every chain (`xs`, a series per chain): the autocovariances about each
# chain's own mean, averaged over the chains, summed by Geyer's initial
# monotone sequence over at most .EMC_SHEAR_CURVE_LAGS lags (a lower bound
# for a coordinate slower than they resolve). NA when it cannot be measured.
.emc_shear_window_iact <- function(xs) {
  lag <- min(min(lengths(xs)) - 2L, .EMC_SHEAR_CURVE_LAGS)
  lag <- lag - (lag + 1L) %% 2L
  if (lag < 1L) return(NA_real_)
  acov <- rowMeans(matrix(vapply(xs, function(x)
    drop(stats::acf(x, lag.max = lag, type = "covariance", plot = FALSE)$acf),
    numeric(lag + 1L)), lag + 1L))
  if (!all(is.finite(acov)) || acov[1L] <= 0) return(NA_real_)
  odd <- seq.int(1L, lag, by = 2L)
  pairs <- (acov[odd] + acov[odd + 1L]) / acov[1L]
  stop <- which(pairs <= 0)
  if (length(stop)) pairs <- pairs[seq_len(stop[1L] - 1L)]
  max(0, 2 * sum(cummin(pairs)) - 1)
}

# The mean coordinates that are slow by the window's own draws: an
# integrated autocorrelation time (`tau`, returned beside the flags) of
# .EMC_SHEAR_CURVE_TAU iterations or more over every chain. `cs_all` holds
# each chain's window draws of the coordinates (a row each); fewer than 50
# draws measure nothing. The rounds a coordinate has had play no part, so
# one the allocation gives no rounds can still ask for a curve.
.emc_shear_curve_slow <- function(cs_all) {
  p <- nrow(cs_all[[1L]])
  if (min(vapply(cs_all, ncol, integer(1))) < 50L)
    return(list(slow = logical(p), tau = rep(NA_real_, p)))
  tau <- vapply(seq_len(p), function(k)
    .emc_shear_window_iact(lapply(cs_all, function(x) x[k, ])), numeric(1))
  list(slow = is.finite(tau) & tau >= .EMC_SHEAR_CURVE_TAU, tau = tau)
}

# The factor one round at the step target takes off a coordinate's
# autocorrelation when its orbit is `r` of the coordinate's posterior sd:
# the round refreshes .EMC_SHEAR_ROUND_REFRESH of the orbit's variance, the
# share r^2 of the coordinate's. An orbit cannot be wider than the margin on
# average, so a measured r above one counts as one.
.emc_shear_round_factor <- function(r) {
  1 - .EMC_SHEAR_ROUND_REFRESH * pmin(1, pmax(0, r))^2
}

# The fixed-frame rounds one round under a curve is worth: with orbits
# r_fixed and r_curve of the coordinate's posterior sd, the n at which n
# fixed-frame rounds take off what the curved round does,
# factor(r_fixed)^n = factor(r_curve). A round's factor is bounded below by
# that of an orbit spanning the margin, so where the fixed fits' orbit
# already covers most of the margin the worth stays near one however wide
# the curve's orbit is measured; the width measure's noise cannot make a
# curve look worth having there. NA when either is not measured.
.emc_shear_curve_worth <- function(r_fixed, r_curve) {
  out <- log(.emc_shear_round_factor(r_curve)) / log(.emc_shear_round_factor(r_fixed))
  out[!(is.finite(r_fixed) & is.finite(r_curve) & r_fixed > 0)] <- NA_real_
  out
}

# Whether a coordinate whose curve failed is still to be left alone. `rec`
# records the failure: the stage, the window's last draw (`at`) and length
# (`len`), and the failures in a row (`fails`). It is tried again once the
# evidence has changed: in another stage, or when the draws since outnumber
# those it was judged on, twice over for each earlier failure.
.emc_shear_curve_waits <- function(rec, now, stage) {
  identical(rec$stage, stage) && now >= rec$at &&
    now - rec$at < 2^(rec$fails - 1L) * rec$len
}

# The failure records with the coordinate `row` added (`failed` TRUE: its
# record counts one more failure, from this window) or taken out. Two rows
# closer than .EMC_SHEAR_CURVE_SAME are one coordinate.
.emc_shear_curve_note <- function(recs, row, failed, now, len, stage) {
  same <- vapply(recs, function(r)
    abs(sum(r$row * row)) / sqrt(sum(r$row^2) * sum(row^2)) >= .EMC_SHEAR_CURVE_SAME,
    logical(1))
  fails <- max(0L, vapply(recs[same], function(r) as.integer(r$fails), integer(1)))
  recs <- recs[!same]
  if (failed)
    recs <- c(recs, list(list(row = row, at = now, len = len, stage = stage,
                              fails = fails + 1L)))
  utils::tail(recs, 2L * length(row))
}

# Whether a chain's curve still serves the window: the window's draws of the
# curve's own coordinate (`cc`, every chain's) centre within its inner knots
# and spread on the knots' scale. A curve that does is not refitted.
.emc_shear_curve_covers <- function(cv, cc) {
  K <- length(cv$knots)
  m <- mean(cc)
  scale <- (cv$knots[K] - cv$knots[1L]) / diff(range(.EMC_SHEAR_CURVE_KNOTS))
  ratio <- stats::sd(cc) / scale
  is.finite(m) && is.finite(ratio) && m >= cv$knots[2L] && m <= cv$knots[K - 1L] &&
    ratio >= 0.5 && ratio <= 2
}

# One subject's likelihood fits along a curve: at each knot's group mean the
# fit at the subject's conditional mode there (prior N(mu, V)), each started
# from its inner neighbour's mode so that the knots follow one mode. NULL
# when any knot has no finite fit.
.emc_shear_curve_subject <- function(chain, s, mus, V, x0) {
  K <- ncol(mus)
  mid <- (K + 1L) %/% 2L
  out <- vector("list", K)
  validated <- isTRUE(chain$ll_cache_validated)
  fit <- function(q, x) {
    r <- tryCatch(.emc_burn_subject_mode(chain, s, x0 = x,
                                         prior = list(mu = mus[, q], var = V),
                                         ll_cache_validated = validated),
                  error = function(e) NULL)
    ok <- is.list(r) && is.list(r$lik) && is.matrix(r$lik$prec) &&
      all(is.finite(r$lik$prec)) && all(is.finite(r$lik$lin)) &&
      length(r$alpha) == length(x) && all(is.finite(r$alpha))
    if (ok) list(prec = r$lik$prec, lin = r$lik$lin, mode = as.numeric(r$alpha))
  }
  out[[mid]] <- fit(mid, x0)
  if (is.null(out[[mid]])) return(NULL)
  for (side in list(seq_len(K - mid) + mid, rev(seq_len(mid - 1L)))) {
    x <- out[[mid]]$mode
    for (q in side) {
      out[[q]] <- fit(q, x)
      if (is.null(out[[q]])) return(NULL)
      x <- out[[q]]$mode
    }
  }
  out
}

# The curve of mean coordinate k from every subject's knot fits (`nodes`, NULL
# for a subject without them): a subject listed in `keep` follows its knots,
# any other has its fixed fit at every knot.
.emc_shear_curve_pack <- function(chain, k, knots, nodes, keep) {
  p <- chain$n_pars
  n <- chain$n_subjects
  K <- length(knots)
  base <- .emc_shear_pack(.emc_shear_fits(chain), p)
  cube <- array(base$cube, c(p, p, n, K))
  lin <- array(base$lin, c(p, n, K))
  for (s in keep) for (q in seq_len(K)) {
    cube[, , s, q] <- nodes[[s]][[q]]$prec
    lin[, s, q] <- nodes[[s]][[q]]$lin
  }
  list(row = unname(solve(chain$shear$plan$dirs)[k, ]), knots = knots, cube = cube, lin = lin)
}

# The mean directions D with direction k replaced by `dir` and its coordinate
# by `row`' mu: `dir` scaled to move that coordinate by one, and the other
# directions cleared of it, so that row k of the inverse is `row`. NULL when
# the result is not a basis.
.emc_shear_pin_dir <- function(D, k, dir, row) {
  if (length(dir) != nrow(D) || length(row) != nrow(D)) return(NULL)
  scale <- sum(row * dir)
  if (!is.finite(scale) || abs(scale) < 1e-12) return(NULL)
  dir <- dir / scale
  out <- D - outer(dir, drop(row %*% D))
  out[, k] <- dir
  if (!all(is.finite(out)) || !is.finite(rcond(out)) || rcond(out) < 1e-12) return(NULL)
  out
}

# The direction along which the group mean covaries with the coordinate
# `row`' mu over the window's draws (every chain's, within chains): moving
# along it leaves the other coordinates at their conditional means, so the
# coordinate's conditional along it is as wide as its margin. `fallback`
# when the window does not give one.
.emc_shear_curve_dir <- function(emc, idx, row, fallback) {
  W <- tryCatch(Reduce(`+`, lapply(emc, function(ch)
    stats::cov(t(ch$samples$theta_mu[, idx, drop = FALSE])))) / length(emc),
    error = function(e) NULL)
  if (is.null(W) || !all(is.finite(W))) return(fallback)
  u <- drop(W %*% row)
  v <- sum(row * u)
  if (!is.finite(v) || v <= 0) return(fallback)
  u / v
}

# The mean coordinate of `dirs` closest to the coordinate `row`' mu (the
# cosine between `row` and the rows of the inverse), and that cosine.
.emc_shear_curve_match <- function(dirs, row) {
  Dinv <- tryCatch(solve(dirs), error = function(e) NULL)
  if (is.null(Dinv) || length(row) != ncol(Dinv)) return(list(k = NA_integer_, cos = 0))
  cs <- abs(drop(Dinv %*% row)) / sqrt(rowSums(Dinv^2) * sum(row^2))
  k <- which.max(cs)
  list(k = k, cos = cs[k])
}

# The coordinates a chain without a curve may try, in order: those another
# chain's curve runs along (the same ridge, already shown to bend), then
# `own`, its slow coordinates. `failed` holds the coordinates (as rows) that
# failed here in earlier rebuilds and are still waited on; candidates along
# them stay out.
.emc_shear_curve_candidates <- function(emc, j, own, failed = NULL) {
  dirs <- emc[[j]]$shear$plan$dirs
  rows <- lapply(emc[-j], function(ch) ch$shear[["curve"]]$row)
  same <- function(r) {
    m <- .emc_shear_curve_match(dirs, r)
    if (m$cos >= .EMC_SHEAR_CURVE_SAME) m$k
  }
  along <- unlist(lapply(rows[!vapply(rows, is.null, logical(1))], same))
  setdiff(unique(c(along, own)), unlist(lapply(failed, same)))
}

# Each chain's curve for the block (see the header). A chain with a curve
# keeps it, at no cost, while the curve covers the window
# (.emc_shear_curve_covers); otherwise it refits the coordinate's curve. A
# chain without one tries its candidates (.emc_shear_curve_candidates) one at
# a time: those another chain's curve runs along, then its own coordinates
# that are slow over the window's draws since burn-in
# (.emc_shear_curve_slow), those that have not failed before the others and
# the slowest first. The curved coordinate's direction
# becomes the one the window measures (.emc_shear_curve_dir), and stays the
# plan's direction k while the curve lasts (.emc_shear_set_plan). The knot
# fits follow the window: its mean and covariance, and its posterior sd of
# the coordinate.
#
# A candidate is judged by what a round along it can deliver. The orbit's
# width is measured at five of the window's states (their geometric mean,
# over the states measured under both sets of fits), as a share r of the
# coordinate's posterior sd over the window, first under the fixed fits. If
# a curve whose orbit spanned the whole margin would still not make a round
# worth .EMC_SHEAR_CURVE_WORTH fixed-frame rounds (.emc_shear_curve_worth),
# no curve is fitted. Otherwise every subject's knots are fitted, the orbit
# is measured again under the curve, and the curve stays only if its round
# is worth that many: a curved coordinate runs a round in every iteration,
# so a round there has to do clearly more than the rounds it replaces.
#
# The work of a rebuild is bounded per chain: the fixed-frame orbits of at
# most .EMC_SHEAR_CURVE_SCREENS coordinates and the knot fits of at most
# .EMC_SHEAR_CURVE_TRIES, the fits as one job per subject and the widths as
# one per state over the rebuild's cores. A coordinate that failed is left
# alone until the evidence has changed (.emc_shear_curve_waits), so the
# failures of one coordinate thin out geometrically. A refit that fails
# leaves the chain its curve once. A coordinate that gains or loses its
# curve starts its round record afresh, at a step from the measured orbit:
# its rounds are another kernel's.
.emc_shear_set_curves <- function(emc, idx, n_cores = 1L) {
  if (!isTRUE(getOption("emc2.shear_curve", TRUE)) ||
      isTRUE(getOption("emc2.shear_delayed", FALSE)) || length(idx) < 2L) return(emc)
  nch <- length(emc)
  p <- emc[[1L]]$n_pars
  now <- max(idx)
  settled <- idx[!emc[[1L]]$samples$stage[idx] %in% c("preburn", "burn")]
  gm <- function(w) if (length(w)) exp(mean(log(w))) else NA_real_
  # Per chain: the window's draws of its coordinates in every chain, its slow
  # coordinates, and the coordinates still waited on.
  view <- lapply(emc, function(ch) {
    plan <- ch$shear$plan
    if (!is.list(plan) || !is.matrix(plan$dirs)) return(NULL)
    Dinv <- tryCatch(solve(plan$dirs), error = function(e) NULL)
    if (is.null(Dinv)) return(NULL)
    draws <- function(at) lapply(emc, function(c2)
      Dinv %*% c2$samples$theta_mu[, at, drop = FALSE])
    slow <- .emc_shear_curve_slow(draws(settled))
    recs <- Filter(function(r) is.list(r) && length(r$row) == p && length(r$fails) == 1L,
                   ch$shear$curve_rec$failed)
    waits <- vapply(recs, .emc_shear_curve_waits, logical(1), now = now, stage = plan$stage)
    # Coordinates that have not failed come before those tried again.
    fails <- integer(p)
    for (r in recs) {
      m <- .emc_shear_curve_match(plan$dirs, r$row)
      if (m$cos >= .EMC_SHEAR_CURVE_SAME) fails[m$k] <- as.integer(r$fails)
    }
    own <- which(slow$slow)
    list(Dinv = Dinv, cs = draws(idx), tau = slow$tau, recs = recs,
         own = own[order(fails[own], -slow$tau[own])],
         blocked = lapply(recs[waits], `[[`, "row"))
  })
  tried <- attempts <- vector("list", nch)
  builds <- integer(nch)
  done <- logical(nch)
  for (attempt in seq_len(.EMC_SHEAR_CURVE_SCREENS)) {
    todo <- vector("list", nch)
    for (j in which(!done)) {
      done[j] <- TRUE
      vw <- view[[j]]
      if (is.null(vw)) next
      st <- emc[[j]]$shear
      plan <- st$plan
      k <- st$curve_rec$k
      had <- length(k) == 1L && is.list(st[["curve"]])
      if (had) {
        if (attempt > 1L ||
            .emc_shear_curve_covers(st[["curve"]], unlist(lapply(vw$cs, function(x) x[k, ]))))
          next
      } else {
        if (builds[j] >= .EMC_SHEAR_CURVE_TRIES) next
        cands <- setdiff(.emc_shear_curve_candidates(emc, j, vw$own, vw$blocked), tried[[j]])
        if (!length(cands)) next
        k <- cands[1L]
        tried[[j]] <- c(tried[[j]], k)
        done[j] <- FALSE
      }
      sd_c <- stats::sd(unlist(lapply(vw$cs, function(x) x[k, ])))
      if (!is.finite(sd_c) || sd_c <= 0) next
      row <- vw$Dinv[k, ]
      dirs <- .emc_shear_pin_dir(plan$dirs, k,
                                 .emc_shear_curve_dir(emc, idx, row, plan$dirs[, k]), row)
      if (is.null(dirs)) dirs <- plan$dirs
      ch <- emc[[j]]
      ch$shear$plan$dirs <- dirs
      ch$shear["curve"] <- list(NULL)
      ch$ll_cache_validated <- .emc_ll_data_prepared(ch$data)
      mu_ref <- rowMeans(ch$samples$theta_mu[, idx, drop = FALSE])
      todo[[j]] <- list(
        k = k, had = had, chain = ch, row = row, sd_c = sd_c, tau = vw$tau[k],
        fixed = .emc_shear_pack(.emc_shear_fits(ch), p),
        s0 = exp(plan$log_step[k]) / .EMC_SHEAR_STEP0,
        at = unique(c(ch$samples$idx, idx[round(length(idx) * (1:4) / 5)])),
        knots = sum(row * mu_ref) + sd_c * .EMC_SHEAR_CURVE_KNOTS,
        mus = mu_ref + outer(dirs[, k], sd_c * .EMC_SHEAR_CURVE_KNOTS),
        V = apply(ch$samples$theta_var[, , idx, drop = FALSE], 1:2, mean),
        x0 = apply(ch$samples$alpha[, , idx, drop = FALSE], 1:2, mean))
    }
    chains <- which(!vapply(todo, is.null, logical(1)))
    if (!length(chains)) break
    started <- proc.time()[["elapsed"]]
    # The orbit widths of the chains js at their states, under the fits
    # fits_of gives each.
    widths <- function(js, fits_of) {
      jobs <- do.call(rbind, lapply(js, function(j)
        data.frame(chain = j, at = todo[[j]]$at)))
      w <- .emc_with_preserved_rng(.emc_rebuild_lapply(seq_len(nrow(jobs)), function(m) {
        td <- todo[[jobs$chain[m]]]
        ch <- td$chain
        ch$samples$idx <- jobs$at[m]
        fits <- fits_of(td)
        tryCatch(.emc_shear_orbit_width(ch, idx, td$k, function(st) fits, td$s0),
                 error = function(e) NA_real_)
      }, n_cores = n_cores))
      w <- vapply(w, function(x) if (is.numeric(x) && length(x) == 1L) x else NA_real_,
                  numeric(1))
      out <- vector("list", nch)
      for (j in js) out[[j]] <- w[jobs$chain == j]
      out
    }
    w_fixed <- widths(chains, function(td) td$fixed)
    open <- chains[vapply(chains, function(j) {
      w <- w_fixed[[j]]
      worth <- .emc_shear_curve_worth(gm(w[is.finite(w)]) / todo[[j]]$sd_c, 1)
      is.finite(worth) && worth >= .EMC_SHEAR_CURVE_WORTH
    }, logical(1))]
    w_curve <- vector("list", nch)
    screen_el <- fit_el <- proc.time()[["elapsed"]] - started
    if (length(open)) {
      jobs <- do.call(rbind, lapply(open, function(j)
        data.frame(chain = j, subject = seq_len(emc[[j]]$n_subjects))))
      nodes <- .emc_with_preserved_rng(.emc_rebuild_lapply(seq_len(nrow(jobs)), function(m) {
        td <- todo[[jobs$chain[m]]]
        s <- jobs$subject[m]
        .emc_shear_curve_subject(td$chain, s, td$mus, td$V, td$x0[, s])
      }, n_cores = n_cores))
      fit_el <- proc.time()[["elapsed"]] - started
      for (j in open) {
        td <- todo[[j]]
        nd <- nodes[jobs$chain == j]
        have <- which(!vapply(nd, is.null, logical(1)))
        builds[j] <- builds[j] + 1L
        todo[[j]]$subjects <- length(have)
        if (length(have))
          todo[[j]]$curve <- c(.emc_shear_curve_pack(td$chain, td$k, td$knots, nd, have),
                               list(dir = unname(td$chain$shear$plan$dirs[, td$k])))
      }
      fitted <- open[vapply(open, function(j) is.list(todo[[j]]$curve), logical(1))]
      if (length(fitted))
        w_curve <- widths(fitted, function(td) c(td$fixed, list(curve = td$curve)))
    }
    el <- proc.time()[["elapsed"]] - started
    for (j in chains) {
      td <- todo[[j]]
      st <- emc[[j]]$shear
      rec <- st$curve_rec
      wf <- w_fixed[[j]]
      wc <- w_curve[[j]]
      built <- is.list(td$curve) && length(wc) == length(wf)
      both <- if (built) is.finite(wf) & is.finite(wc) else is.finite(wf)
      w <- c(fixed = gm(wf[both]), curve = if (built) gm(wc[both]) else NA_real_)
      r <- w / td$sd_c
      worth <- .emc_shear_curve_worth(r[["fixed"]], r[["curve"]])
      passed <- built && is.finite(worth) && worth >= .EMC_SHEAR_CURVE_WORTH
      strikes <- if (passed) 0L else (if (length(rec$strikes)) rec$strikes else 0L) + 1L
      kept <- passed || (td$had && strikes < 2L)
      info <- list(k = td$k, passed = passed, kept = kept, tau = td$tau,
                   r_fixed = r[["fixed"]], r_curve = r[["curve"]], worth = worth,
                   w_fixed = w[["fixed"]], w_curve = w[["curve"]],
                   fitted = j %in% open, subjects = if (j %in% open) td$subjects else 0L,
                   screen_seconds = screen_el, fit_seconds = fit_el - screen_el,
                   seconds = el)
      if (passed) {
        st$curve <- td$curve
        st$plan$dirs <- td$chain$shear$plan$dirs
      } else if (!kept) st["curve"] <- list(NULL)
      # A chain that keeps its curve is done; one without tries its next
      # candidate.
      if (kept) done[j] <- TRUE
      failed <- if (kept && !passed) view[[j]]$recs else
        .emc_shear_curve_note(view[[j]]$recs, td$row, !kept, now, length(idx), st$plan$stage)
      view[[j]]$recs <- failed
      attempts[[j]] <- c(attempts[[j]], list(info))
      count <- function(x) if (length(x)) x else 0
      st$curve_rec <- c(
        if (kept) list(k = td$k),
        list(failed = failed, strikes = strikes, last = info, tried = attempts[[j]],
             screens = count(rec$screens) + 1L,
             builds = count(rec$builds) + (j %in% open),
             seconds = count(rec$seconds) + el / length(chains)))
      if (kept != td$had)
        st$plan <- .emc_shear_restart_coord(
          st$plan, td$k, if (kept) w[["curve"]] else w[["fixed"]])
      emc[[j]]$shear <- st
    }
  }
  emc
}

# Coordinate k's round record cleared, and its step set from an orbit width
# when one was measured.
.emc_shear_restart_coord <- function(plan, k, width = NA_real_) {
  for (nm in c("pool", "last")) {
    rec <- plan[[nm]]
    if (!is.list(rec)) next
    for (f in c("E", "n", "moved", "N", "iters")) if (length(rec[[f]]) >= k) rec[[f]][k] <- 0
    if (is.matrix(rec$K) && nrow(rec$K) >= k) rec$K[k, ] <- 0
    plan[[nm]] <- rec
  }
  if (is.finite(width) && width > 0)
    plan$log_step[k] <- min(.EMC_SHEAR_LOG_STEP_RANGE[2], max(.EMC_SHEAR_LOG_STEP_RANGE[1],
                                                               log(.EMC_SHEAR_STEP0 * width)))
  plan
}

# How a rebuilt plan's coordinates inherit the old plan's step scales and
# round records. A new mean direction is a combination sum_m c_m d_m of the
# old ones, each of which moves mu by one posterior sd of its direction; it
# takes old direction m's per-round statistics with weight c_m^2 / sum c^2,
# rescaled by the length ratio |d_m| / |d| (round moves are in units of
# the direction). A direction the old plan already had inherits its own
# record alone. Directions with nearly equal lambda are fixed only up to a
# rotation among themselves, so a rebuild can turn them into new
# combinations of the old ones; matching one old direction each would leave
# them unmeasured. w and scale are 2p x 2p (new by old coordinate); a log
# group sd inherits its own record unscaled, and an old basis that cannot
# be inverted passes nothing to the mean directions.
.emc_shear_carry_weights <- function(old_dirs, new_dirs) {
  p <- ncol(new_dirs)
  w <- scale <- matrix(0, 2L * p, 2L * p)
  sd_k <- p + seq_len(p)
  w[cbind(sd_k, sd_k)] <- 1
  scale[] <- 1
  cf <- tryCatch(solve(old_dirs, new_dirs), error = function(e) NULL)
  if (!is.null(cf) && all(is.finite(cf))) {
    wm <- t(cf^2)
    tot <- rowSums(wm)
    wm[tot > 0, ] <- wm[tot > 0, , drop = FALSE] / tot[tot > 0]
    w[seq_len(p), seq_len(p)] <- wm
    scale[seq_len(p), seq_len(p)] <- outer(1 / sqrt(colSums(new_dirs^2)),
                                           sqrt(colSums(old_dirs^2)))
  }
  list(w = w, scale = scale)
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
  same_stage <- identical(stage, old$stage)
  if (is.list(old) && is.matrix(old$dirs) && all(dim(old$dirs) == p)) {
    cw <- .emc_shear_carry_weights(old$dirs, plan$dirs)
    # A curve keeps its direction and coordinate: they take the place of the
    # new direction that is mostly the old one (.emc_shear_pin_dir).
    k_old <- st$curve_rec$k
    st$curve_rec$k <- NULL
    if (length(k_old) == 1L && is.list(st[["curve"]])) {
      k_new <- which.max(cw$w[seq_len(p), k_old])
      pinned <- .emc_shear_pin_dir(plan$dirs, k_new, st$curve$dir, st$curve$row)
      if (!is.null(pinned)) {
        plan$dirs <- pinned
        st$curve_rec$k <- k_new
        cw <- .emc_shear_carry_weights(old$dirs, plan$dirs)
      }
    }
    if (is.null(st$curve_rec$k)) st["curve"] <- list(NULL)
    got <- which(rowSums(cw$w) > 0)
    blk <- old$block
    plan$log_step[got] <- rowSums(cw$w * (rep(old$log_step, each = 2L * p) +
                                            log(cw$scale)))[got]
    carry <- function(rec) {
      out <- zero
      if (is.null(rec)) return(out)
      rec <- .emc_shear_record(rec, p)
      sq <- cw$w * cw$scale^2
      out$E <- drop(sq %*% rec$E)
      out$n <- drop(cw$w %*% rec$n)
      out$moved <- drop(cw$w %*% rec$moved)
      out$N <- drop(sq %*% rec$N)
      out$iters <- drop(cw$w %*% rec$iters)
      out$K <- cw$w %*% rec$K
      out
    }
    recent <- carry(blk)
    recent$iters[got] <- .emc_alloc_block_iters(chain)
    production_pool <- identical(stage, "sample") && same_stage
    if (production_pool || !identical(stage, "sample")) {
      prev <- carry(if (same_stage) old$pool else old$last)
      if (!same_stage && !identical(stage, "sample")) {
        low <- prev$n < 5L
        prev$E[low] <- 0
        prev$n[low] <- 0L
        prev$moved[low] <- 0L
        prev$N[low] <- 0
        prev$K[low, ] <- 0L
        prev$iters[low] <- 0
      }
      d <- .EMC_ALLOC_EVIDENCE_DECAY
      plan$pool <- list(E = d * prev$E + recent$E, n = d * prev$n + recent$n,
                        moved = d * prev$moved + recent$moved, N = d * prev$N + recent$N,
                        K = .emc_alloc_khist_sum(d * prev$K, recent$K),
                        iters = d * prev$iters + recent$iters)
    } else {
      # Entering production the pool starts empty and the last block alone
      # measures a coordinate; one it barely ran keeps its discounted
      # warm-up evidence instead of reverting to the unmeasured prior.
      prev <- carry(old$last)
      low <- recent$n < 5L
      if (any(low)) {
        d <- .EMC_ALLOC_EVIDENCE_DECAY
        for (f in c("E", "n", "moved", "N", "iters"))
          recent[[f]][low] <- d * prev[[f]][low] + recent[[f]][low]
        prev$K[!low, ] <- 0L
        recent$K <- .emc_alloc_khist_sum(d * prev$K, recent$K)
      }
    }
  }
  plan$last <- if (identical(stage, "sample") && !same_stage) recent else plan$pool
  # Discarded stages adapt the scales round by round; production only here.
  if (identical(stage, "sample")) plan <- .emc_shear_tune_steps(plan, recent)
  plan$block <- zero
  st$plan <- plan
  st$walls <- tryCatch(.emc_shear_walls(chain, idx), error = function(e) NULL)
  chain$shear <- st
  chain
}

# The group log sds over the draws idx (p x length(idx)).
.emc_shear_window_lsd <- function(chain, idx) {
  out <- if (getOption("emc2.native_window_lsd", TRUE))
    emc_window_lsd_cpp(chain$samples$theta_var, idx, as.integer(chain$n_pars))
  if (!is.null(out)) return(out)
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
  if (is.null(fits$cube)) fits <- .emc_shear_pack(fits, length(mu))
  fits <- .emc_shear_fits_at(fits, mu)
  out <- emc_shear_state_cpp(mu, V, fits$cube, fits$lin, spec)
  if (!isTRUE(out$ok)) return(NULL)
  fm <- list(m = out$m, R = out$R, logdet = out$logdet)
  list(mu = mu, V = V, Vi = out$Vi, chol_V = out$chol_V, fm = fm,
       logdet = sum(fm$logdet), lp_v = out$lp_v,
       lp = .emc_shear_log_prior_mu(mu, prior) + out$lp_v)
}

# The likelihood approximations at group mean mu: the packed fits, or their
# curve's interpolation there.
.emc_shear_fits_at <- function(fits, mu) {
  if (is.null(fits$curve)) return(fits)
  emc_shear_curve_fits_cpp(as.numeric(mu), fits$curve, ncol(fits$lin))
}

# Whether a mean step u_mu changes the fits (moves the curve's coordinate).
.emc_shear_curve_moves <- function(fits, u_mu) {
  cv <- fits$curve
  if (is.null(cv)) return(FALSE)
  abs(sum(cv$row * u_mu)) > 1e-9 * sqrt(sum(cv$row^2) * sum(u_mu^2))
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
.emc_shear_propose <- function(cur, u, alpha, fits, prior, spec, walls = NULL) {
  p <- length(cur$mu)
  u_sd <- u[p + seq_len(p)]
  if (all(u_sd == 0) && is.null(walls) && !.emc_shear_curve_moves(fits, u[seq_len(p)])) {
    new <- .emc_shear_shift_state(cur, u[seq_len(p)], prior)
    if (!is.finite(new$lp)) return(NULL)
    return(list(state = new, alpha = alpha + new$fm$m - cur$fm$m, log_jac = 0))
  }
  .emc_shear_candidates(cur, u, 1, alpha, fits, prior, spec, walls)[[1L]]
}

# The candidates of one round: for each offset t[k] the map for u = t[k] v
# (as .emc_shear_propose), built in one call (emc_shear_candidates_cpp).
# Per candidate NULL when the moved state is invalid, else the moved state,
# the moved subjects, their summed hierarchical density (hier, as
# .emc_shear_log_hier), log|det| of the map's Jacobian and the walls' part
# of it.
.emc_shear_candidates <- function(cur, v, t, alpha, fits, prior, spec, walls = NULL) {
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
                                  fits$cube, fits$lin, has, U, b, Uc, spec, fits$curve)
  u_sd <- v[p + seq_len(p)]
  u_mu <- v[seq_len(p)]
  lapply(seq_along(t), function(k) {
    o <- out[[k]]
    if (!isTRUE(o$ok)) return(NULL)
    fm <- list(m = o$m, R = o$R, logdet = o$logdet)
    V <- o$V
    dimnames(V) <- dimnames(cur$V)
    new <- list(mu = cur$mu + t[k] * u_mu, V = V, Vi = o$Vi, chol_V = o$chol_V,
                fm = fm, logdet = sum(fm$logdet), lp_v = o$lp_v)
    new$lp <- .emc_shear_log_prior_mu(new$mu, prior) + new$lp_v
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
  key <- sampler[c("lik_laplace", "n_pars", "n_subjects", "type", "prior",
                   "is_blocked", "nuisance", "par_group")]
  key$walls <- sampler$shear$walls
  key$curve <- .emc_shear_curve_of(sampler)
  .emc_context_memo("shear", key, function() {
    fits <- .emc_shear_fits(sampler)
    packed <- .emc_shear_pack(fits, sampler$n_pars)
    packed$curve <- key$curve
    walls <- .emc_shear_walls_of(sampler)
    list(fits = fits, packed = packed, spec = .emc_shear_v_spec(sampler), walls = walls,
         args = .emc_shear_round_args(packed, walls, sampler$n_pars))
  })
}

# The block's curve when it matches the sampler's dimensions. The two-stage
# round selects on one fixed quadratic fit per subject, so it runs without.
.emc_shear_curve_of <- function(sampler) {
  cv <- sampler$shear[["curve"]]
  p <- sampler$n_pars
  n <- sampler$n_subjects
  ok <- is.list(cv) && length(cv$row) == p && length(cv$knots) >= 2L &&
    length(cv$cube) == p * p * n * length(cv$knots) &&
    length(cv$lin) == p * n * length(cv$knots) &&
    !isTRUE(getOption("emc2.shear_delayed", FALSE))
  if (ok) cv else NULL
}

# The fits and walls as the native rounds take them.
.emc_shear_round_args <- function(packed, walls, p) {
  if (is.null(walls)) {
    has <- logical(0); U <- matrix(0, p, 0); b <- numeric(0)
    Uc <- array(0, c(p, max(p - 1L, 0L), 0))
  } else {
    has <- walls$has; U <- walls$U; b <- ifelse(has, walls$b, 0)
    Uc <- .emc_shear_uperp_cube(walls, p)
  }
  list(cube = packed$cube, lin = packed$lin, has = has, U = U, b = b, Uc = Uc)
}

# The last value built for a slot, reused while its inputs are unchanged: a
# move's context depends only on proposal state that is fixed between
# rebuilds.
.emc_context_store <- new.env(parent = emptyenv())
.emc_context_memo <- function(slot, key, build) {
  hit <- .emc_context_store[[slot]]
  if (!is.null(hit) && identical(hit$key, key)) return(hit$value)
  value <- build()
  .emc_context_store[[slot]] <- list(key = key, value = value)
  value
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
.emc_shear_round <- function(x, v, step, sampler, ctx, t = NULL) {
  # Elapsed time building the candidates and in their likelihood round.
  build_el0 <- proc.time()[["elapsed"]]
  # `t` holds the candidate steps when they have been drawn already.
  if (is.null(t)) {
    t_aux <- step * stats::rnorm(1L)
    t <- t_aux + step * stats::rnorm(.EMC_SHEAR_CANDIDATES)
  }
  p <- nrow(x$alpha)
  n <- ncol(x$alpha)
  # Mean steps lie on a line when no subject has a wall.
  mean_step <- all(v[p + seq_len(p)] == 0) && is.null(ctx$walls) &&
    !.emc_shear_curve_moves(ctx$packed, v[seq_len(p)])
  if (mean_step) {
    d <- v[seq_len(p)]
    line <- .emc_shear_line(x$cur, d, x$alpha, t, sampler$prior)
    moved_alpha <- function(i) x$alpha + t[i] * line$Delta
    prior_delta <- line$delta
    ok <- which(is.finite(prior_delta))
  } else {
    if (.emc_shear_round_native_ok(x, v, sampler)) {
      done <- .emc_shear_round_native(x, v, t, sampler, ctx, build_el0)
      if (!is.null(done)) return(done)
    }
    props <- .emc_shear_candidates(x$cur, v, t, x$alpha, ctx$packed, sampler$prior, ctx$spec,
                                   ctx$walls)
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
  if (isTRUE(getOption("emc2.shear_delayed", FALSE)))
    return(.emc_shear_round_delayed(x, t, ok, moved_alpha, prior_delta, mean_step,
                                    if (mean_step) d else NULL,
                                    if (mean_step) NULL else props, sampler, ctx,
                                    build_el0))
  if (length(ok)) {
    pn <- sampler$par_names
    cand <- lapply(ok, moved_alpha)
    rows <- emc_candidate_rows_cpp(cand, pn)
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

# Whether the round can run in one native call (emc_shear_round_cpp): a
# direction without a mean part, subject likelihoods computed in this process
# on the stored-argument path, and the one-stage selection.
.emc_shear_round_native_ok <- function(x, v, sampler) {
  p <- nrow(x$alpha)
  isTRUE(getOption("emc2.native_shear_round", TRUE)) &&
    !isTRUE(getOption("emc2.shear_delayed", FALSE)) &&
    all(v[seq_len(p)] == 0) &&
    (is.null(x$wpool) || !isTRUE(x$wpool$alive)) &&
    isTRUE(sampler$ll_cache_validated) && is.null(sampler$marginalise) &&
    is.list(sampler$data) && is.data.frame(sampler$data[[1L]]) &&
    is.null(attr(sampler$data[[1L]], "pooled_values"))
}

# The round of .emc_shear_round() in one native call, on the same random
# stream. NULL when a subject's likelihood needs the full path; no random
# number has been drawn then.
.emc_shear_round_native <- function(x, v, t, sampler, ctx, build_el0) {
  p <- nrow(x$alpha)
  cur <- x$cur
  fits <- ctx$packed
  if (is.null(fits$cube)) fits <- .emc_shear_pack(fits, p)
  # A log-sd round keeps mu, and with it the curve's fits.
  fits <- .emc_shear_fits_at(fits, cur$mu)
  walls <- ctx$walls
  if (is.null(walls)) {
    has <- logical(0); U <- matrix(0, p, 0); b <- numeric(0)
    Uc <- array(0, c(p, max(p - 1L, 0L), 0))
  } else {
    has <- walls$has; U <- walls$U; b <- ifelse(has, walls$b, 0)
    Uc <- .emc_shear_uperp_cube(walls, p)
  }
  f0 <- .emc_shear_frame_arrays(cur)
  o <- emc_shear_round_cpp(t, v, x$alpha, cur$mu, cur$V, f0$m, f0$R, f0$logdet,
                           fits$cube, fits$lin, has, U, b, Uc, ctx$spec,
                           .emc_shear_log_prior_mu(cur$mu, sampler$prior),
                           x$cur_hier, cur$lp, x$ll, sampler$data, sampler$model,
                           sampler$par_names)
  if (isTRUE(o$fallback)) return(NULL)
  el <- proc.time()[["elapsed"]] - build_el0
  x$build_el <- x$build_el + max(0, el - o$ll_seconds)
  x$ll_el <- x$ll_el + o$ll_seconds
  if (o$evaluated) x$rounds <- x$rounds + 1L
  x$E <- o$E
  x$moved <- o$j > 0L
  x$t <- if (x$moved) t[o$j] else 0
  if (x$moved) {
    V <- o$V
    dimnames(V) <- dimnames(cur$V)
    fm <- list(m = o$m, R = o$R, logdet = o$logdet)
    x$cur <- list(mu = o$mu, V = V, Vi = o$Vi, chol_V = o$chol_V, fm = fm,
                  logdet = sum(fm$logdet), lp_v = o$lp_v, lp = o$lp)
    x$alpha <- o$alpha
    x$cur_hier <- .emc_shear_log_hier(x$alpha, x$cur)
    x$ll <- o$ll
  }
  x
}

# A run of rounds along the log-sd coordinates `ks` in one native call
# (emc_shear_rounds_cpp), each from the state the last one left and on the
# same random stream as .emc_shear_round() round by round. Returns x after
# the rounds that ran, their count `done`, each one's net step, expected
# squared jump and whether it moved, the adapted step sizes, and `t_next`,
# the steps already drawn for a round the native path declined.
.emc_shear_rounds_native <- function(x, ks, log_step, blk_n, adapt, sampler, ctx) {
  el0 <- proc.time()[["elapsed"]]
  p <- nrow(x$alpha)
  cur <- x$cur
  fits <- ctx$packed
  if (is.null(fits$cube)) fits <- .emc_shear_pack(fits, p)
  # A log-sd round keeps mu, and with it the curve's fits.
  fits <- .emc_shear_fits_at(fits, cur$mu)
  walls <- ctx$walls
  if (is.null(walls)) {
    has <- logical(0); U <- matrix(0, p, 0); b <- numeric(0)
    Uc <- array(0, c(p, max(p - 1L, 0L), 0))
  } else {
    has <- walls$has; U <- walls$U; b <- ifelse(has, walls$b, 0)
    Uc <- .emc_shear_uperp_cube(walls, p)
  }
  f0 <- .emc_shear_frame_arrays(cur)
  o <- emc_shear_rounds_cpp(as.integer(ks), log_step, as.integer(blk_n), adapt,
                            .EMC_SHEAR_LOG_STEP_RANGE, .EMC_SHEAR_MOVE_TARGET,
                            .EMC_SHEAR_CANDIDATES, x$alpha, cur$mu, cur$V, cur$chol_V,
                            f0$m, f0$R, f0$logdet, x$cur_hier, cur$lp,
                            .emc_shear_log_prior_mu(cur$mu, sampler$prior), x$ll,
                            fits$cube, fits$lin, has, U, b, Uc, ctx$spec,
                            sampler$data, sampler$model, sampler$par_names)
  el <- proc.time()[["elapsed"]] - el0
  x$build_el <- x$build_el + max(0, el - o$ll_seconds)
  x$ll_el <- x$ll_el + o$ll_seconds
  x$rounds <- x$rounds + o$evaluated
  if (o$any_moved) {
    V <- o$V
    dimnames(V) <- dimnames(cur$V)
    fm <- list(m = o$m, R = o$R, logdet = o$logdet)
    x$cur <- list(mu = o$mu, V = V, Vi = o$Vi, chol_V = o$chol_V, fm = fm,
                  logdet = sum(fm$logdet), lp_v = o$lp_v, lp = o$lp)
    x$alpha <- o$alpha
    x$cur_hier <- o$cur_hier
    x$ll <- o$ll
  }
  list(x = x, done = o$done, t = o$t, E = o$E, moved = o$moved, log_step = o$log_step,
       t_next = if (length(o$t_next)) o$t_next)
}

# The round in two stages. The candidate (or the current state) is selected
# with the subjects' likelihoods replaced by their quadratic fits, a step that
# is reversible for the posterior under those fits; the selected candidate is
# then accepted with the ratio of the likelihood to its fit at the candidate
# against the same ratio at the current state, which restores the posterior
# itself. Only the selected candidate's likelihoods are computed. x$E is the
# expected squared jump given the selection.
.emc_shear_round_delayed <- function(x, t, ok, moved_alpha, prior_delta, mean_step,
                                     d, props, sampler, ctx, build_el0) {
  p <- nrow(x$alpha)
  n <- ncol(x$alpha)
  x$E <- 0
  x$moved <- FALSE
  x$t <- 0
  if (!length(ok)) {
    x$build_el <- x$build_el + proc.time()[["elapsed"]] - build_el0
    return(x)
  }
  fits <- ctx$packed
  if (is.null(fits$cube)) fits <- .emc_shear_pack(fits, p)
  cand <- lapply(ok, moved_alpha)
  fit_ll <- emc_shear_surrogate_cpp(c(list(x$alpha), cand), fits$cube, fits$lin)
  fit_delta <- fit_ll[-1L] - fit_ll[1L]
  lw <- fit_delta + prior_delta[ok]
  lw[!is.finite(lw)] <- -Inf
  w <- c(0, lw)
  w <- exp(w - max(w))
  w <- w / sum(w)
  pick <- sample.int(length(w), 1L, prob = w) - 1L
  x$rounds <- x$rounds + 1L
  if (pick == 0L) {
    x$build_el <- x$build_el + proc.time()[["elapsed"]] - build_el0
    return(x)
  }
  j <- ok[pick]
  rows <- emc_candidate_rows_cpp(cand[pick], sampler$par_names)
  ll_el0 <- proc.time()[["elapsed"]]
  x$build_el <- x$build_el + ll_el0 - build_el0
  evaluated <- .emc_wpool_subject_ll(x$wpool, rows, sampler, x$wpool_part,
                                     source = "shear")
  x$ll_el <- x$ll_el + proc.time()[["elapsed"]] - ll_el0
  x$wpool <- evaluated$pool
  ll_new <- as.numeric(evaluated$ll)
  log_accept <- sum(ll_new) - sum(x$ll) - fit_delta[pick]
  u <- stats::runif(1L)
  if (!is.finite(log_accept)) return(x)
  x$E <- min(1, exp(log_accept)) * t[j]^2
  if (log(u) >= log_accept) return(x)
  x$moved <- TRUE
  x$t <- t[j]
  if (mean_step) {
    x$cur <- .emc_shear_shift_state(x$cur, t[j] * d, sampler$prior)
    x$alpha <- cand[[pick]]
  } else {
    x$cur <- props[[j]]$state
    x$alpha <- props[[j]]$alpha
  }
  x$cur_hier <- .emc_shear_log_hier(x$alpha, x$cur)
  x$ll <- ll_new
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
  cur <- .emc_shear_state(mu, as.matrix(pars$tvar), ctx$packed,
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
  order <- take[sample.int(length(take))]
  one_round <- function(k, t = NULL) {
    x <<- .emc_shear_round(x, .emc_shear_coord_dir(plan, k), exp(plan$log_step[k]),
                           sampler, ctx, t)
    net[k] <<- net[k] + x$t
    blk$E[k] <<- blk$E[k] + x$E
    blk$n[k] <<- blk$n[k] + 1L
    blk$moved[k] <<- blk$moved[k] + as.integer(x$moved)
    if (adapt) {
      plan$log_step[k] <<- min(.EMC_SHEAR_LOG_STEP_RANGE[2], max(.EMC_SHEAR_LOG_STEP_RANGE[1],
        plan$log_step[k] + (x$moved - .EMC_SHEAR_MOVE_TARGET) / sqrt(blk$n[k])))
    }
    any_moved <<- any_moved || x$moved
  }
  # Runs of log-sd coordinates go to the native path in one call each.
  sd_run <- order > p
  i <- 1L
  while (i <= length(order)) {
    k <- order[i]
    if (!sd_run[i] || !.emc_shear_round_native_ok(x, .emc_shear_coord_dir(plan, k), sampler)) {
      one_round(k)
      i <- i + 1L
      next
    }
    last <- i
    while (last < length(order) && sd_run[last + 1L]) last <- last + 1L
    ks <- order[i:last]
    o <- .emc_shear_rounds_native(x, ks, plan$log_step, blk$n, adapt, sampler, ctx)
    x <- o$x
    for (r in seq_len(o$done)) {
      kr <- ks[r]
      net[kr] <- net[kr] + o$t[r]
      blk$E[kr] <- blk$E[kr] + o$E[r]
      blk$n[kr] <- blk$n[kr] + 1L
      blk$moved[kr] <- blk$moved[kr] + as.integer(o$moved[r])
    }
    plan$log_step <- o$log_step
    any_moved <- any_moved || any(o$moved[seq_len(o$done)])
    i <- i + o$done
    if (!is.null(o$t_next)) {
      one_round(order[i], o$t_next)
      i <- i + 1L
    }
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
