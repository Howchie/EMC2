# ---------------------------------------------------------------------------
# Stage 0 (option 3): model-agnostic marginalization of the shared subject-level
# non-decision time t0 out of the per-subject race likelihood.
#
#   Lbar_s(theta) = int p(t0 | eta) * L_s(t0, theta) dt0
#
# where L_s is the COMPLETE subject likelihood (response + omission trials) and
# the SAME t0 is inserted into every trial at each quadrature node (coherent
# marginalization, not the incoherent trialwise prod-of-integrals).
#
# t0 is sampled on the log scale in EMC2, so we model p(t0|eta) as lognormal:
#   log t0 ~ N(mu, sigma^2)   (eta = (mu, sigma))
# and integrate on the log-t0 axis with Gauss-Hermite quadrature against that
# normal. This matches EMC2's hierarchical prior (normal on the sampled scale),
# so Stage 1/2 slot in without changing the marginalization semantics.
#
# The wrapper is model-agnostic: it overwrites the named t0 column(s) of the
# proposal matrix and re-calls the already-registered race likelihood via
# EMC2:::calc_ll_manager. Nothing here is specific to LBA/RDM/LNR or GNG.
# ---------------------------------------------------------------------------

# Gauss-Legendre rule on [-1, 1], Golub-Welsch on the Legendre Jacobi matrix
# (companion to src/gl_quad.h). We integrate on the log-t0 axis over the
# feasible interval (-Inf, log m_s]; there the integrand p(t0|eta)*L(t0) is
# smooth, so GL converges fast. Above m_s = min_i r_i - eps a response would be
# infeasible (density 0, min_ll-floored), which the proposal excludes by capping
# t0 at m_s -- capping both removes that mass AND removes the integrand kinks
# that make a plain Gauss-Hermite rule over the whole normal non-convergent.
gl_rule <- function(n) {
  stopifnot(n >= 2L)
  i <- seq_len(n - 1L)
  offdiag <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- offdiag
  J[cbind(i + 1L, i)] <- offdiag
  e <- eigen(J, symmetric = TRUE)
  x <- e$values
  w <- 2 * (e$vectors[1, ]^2)          # sum(w) = 2 on [-1, 1]
  ord <- order(x)
  list(x = x[ord], w = w[ord])
}

# Feasibility ceiling m_s = min finite rt - eps (t0 must be below every response
# RT). Returned on the log scale to match the sampled t0 axis.
log_ms_from_dadm <- function(dadm, eps = 1e-3) {
  rt <- dadm$rt
  mn <- suppressWarnings(min(rt[is.finite(rt)]))
  if (!is.finite(mn)) return(Inf)      # no response trials -> no feasibility cap
  log(max(mn - eps, .Machine$double.eps))
}

# Natural-scale LOWER bound on t0, read from the model's bound object (default
# 0.05 for race models; user-overridable via design()'s `bound` arg, which
# map.R merges over the default). Below this bound the compiled likelihood's
# c_do_bound (particle_ll.cpp:2988) marks EVERY trial not-ok and floors the whole
# subject to n_trials*min_ll -- a dead, flat region. Returned on the log scale to
# match the sampled t0 axis. Assumes t0 ~ 1 (identity map, transform = exp), so
# the coefficient axis is log(natural t0); for structured t0 designs the mapping
# is affine and this bound is a conservative floor.
log_lo_from_model <- function(model, t0_ptype = "t0") {
  mm <- model$bound$minmax
  if (is.null(mm) || !(t0_ptype %in% colnames(mm))) return(-Inf)
  lo <- mm[1, t0_ptype]
  if (!is.finite(lo) || lo <= 0) return(-Inf)  # no informative lower bound
  log(lo)
}

# Hoisted per-subject likelihood closure: builds the design/constant/type payload
# ONCE and calls the compiled kernel directly, so quadrature/MCMC loops don't pay
# calc_ll_manager's per-call model()-instantiation + design-rebuild overhead.
# Returns function(proposals matrix) -> vector of summed subject log-liks.
make_ll_fun <- function(dadm, model_gen, min_ll = log(1e-10)) {
  model <- model_gen()
  if (is.null(model$c_name))
    stop("make_ll_fun: compiled path only (model$c_name is NULL).")
  dadm <- EMC2:::.cache_ll_data_attrs(dadm)
  p_types   <- names(model$p_types)
  designs   <- EMC2:::.oo_expanded_designs(dadm)
  constants <- attr(dadm, "constants"); if (is.null(constants)) constants <- NA
  ctype <- model$c_name; bnd <- model$bound; tr <- model$transform
  ptr <- model$pre_transform; trend <- model$trend
  fun <- function(proposals) {
    if (is.null(dim(proposals)))
      proposals <- matrix(proposals, nrow = 1L, dimnames = list(NULL, names(proposals)))
    EMC2:::calc_ll_oo(proposals, dadm, constants = constants, designs = designs,
                      type = ctype, bnd, tr, ptr, p_types = p_types,
                      min_ll = min_ll, trend = trend)
  }
  attr(fun, "log_ms") <- log_ms_from_dadm(dadm)  # upper: feasibility (min rt)
  attr(fun, "log_lo") <- log_lo_from_model(model) # lower: t0 param bound
  fun
}

# log-sum-exp over rows-fixed columns: given a matrix M (n_particles x n_nodes)
# of log-terms, return length-n_particles log(sum_k exp(M[,k])).
.lse_rows <- function(M) {
  m <- apply(M, 1L, max)
  m[!is.finite(m)] <- 0
  m + log(rowSums(exp(M - m)))
}

# Marginal log-likelihood, integrating the shared log-t0 over N(mu, sigma^2)
# truncated at the feasibility ceiling log m_s (unnormalized, matching the
# proposal's int_l^{m_s} p(t0|eta) L(t0) dt0).
#   proposals : n_particles x n_pars matrix (must contain the t0 column;
#               its values are IGNORED -- t0 is integrated out).
#   ll_fun    : closure from make_ll_fun(dadm, model): proposals -> vector of lls.
#   mu, sigma : hyperparameters eta on the log-t0 (sampled) scale. Length 1, or
#               length n_particles for a per-particle eta.
#   n_nodes   : number of Gauss-Legendre nodes on the feasible interval.
#   span      : low-tail half-width in sigma units (upper end is capped at m_s).
#   log_ms    : feasibility ceiling on the log scale (default: attr(ll_fun,"log_ms")).
#   t0_name   : name of the t0 coefficient column (default "t0"; must be t0~1).
# Returns a length-n_particles vector of log Lbar_s.
marginal_ll_t0 <- function(proposals, ll_fun, mu, sigma,
                           n_nodes = 64L, span = 6, log_ms = NULL, log_lo = NULL,
                           t0_name = "t0") {
  if (is.null(dim(proposals)))
    proposals <- matrix(proposals, nrow = 1L, dimnames = list(NULL, names(proposals)))
  stopifnot(t0_name %in% colnames(proposals))
  np <- nrow(proposals)
  mu    <- rep_len(mu, np)
  sigma <- rep_len(sigma, np)
  if (is.null(log_ms)) log_ms <- attr(ll_fun, "log_ms")
  if (is.null(log_ms)) stop("marginal_ll_t0: log_ms unknown (pass it or use make_ll_fun).")
  if (is.null(log_lo)) log_lo <- attr(ll_fun, "log_lo")
  if (is.null(log_lo)) log_lo <- -Inf                   # no param lower bound

  lo <- pmax(mu - span * sigma, log_lo)                 # t0 param-bound floor
  hi <- pmin(mu + span * sigma, log_ms)                 # feasibility cap, per particle
  gl <- gl_rule(n_nodes)
  half <- (hi - lo) / 2; mid <- (hi + lo) / 2           # affine map [-1,1] -> [lo,hi]

  terms <- matrix(-Inf, nrow = np, ncol = n_nodes)      # log( |J| w_k N(x_k) L(x_k) )
  ok <- half > 0
  for (k in seq_len(n_nodes)) {
    xk <- mid + half * gl$x[k]                          # log-t0 node, per particle
    pk <- proposals; pk[ok, t0_name] <- xk[ok]
    llk <- rep(-Inf, np)
    if (any(ok)) llk[ok] <- ll_fun(pk[ok, , drop = FALSE])
    terms[, k] <- log(gl$w[k]) + log(half) + dnorm(xk, mu, sigma, log = TRUE) + llk
  }
  out <- .lse_rows(terms)
  out[!ok] <- -Inf
  out
}

# Brute-force reference: fine trapezoid over the same truncated interval.
# Used only to validate marginal_ll_t0(); not for production.
grid_marginal_ll_t0 <- function(proposals, ll_fun, mu, sigma,
                                 n_grid = 800L, span = 6, log_ms = NULL, log_lo = NULL,
                                 t0_name = "t0") {
  if (is.null(dim(proposals)))
    proposals <- matrix(proposals, nrow = 1L, dimnames = list(NULL, names(proposals)))
  stopifnot(nrow(proposals) == 1L)                      # reference: single particle
  if (is.null(log_ms)) log_ms <- attr(ll_fun, "log_ms")
  if (is.null(log_lo)) log_lo <- attr(ll_fun, "log_lo")
  if (is.null(log_lo)) log_lo <- -Inf
  hi <- min(mu + span * sigma, log_ms)
  xs <- seq(max(mu - span * sigma, log_lo), hi, length.out = n_grid)
  ld <- dnorm(xs, mu, sigma, log = TRUE)
  ll <- numeric(n_grid)
  for (g in seq_len(n_grid)) {
    pk <- proposals; pk[, t0_name] <- xs[g]
    ll[g] <- ll_fun(pk)
  }
  logterm <- ld + ll
  m <- max(logterm)
  dx <- xs[2] - xs[1]
  wt <- rep(dx, n_grid); wt[1] <- wt[n_grid] <- dx / 2   # trapezoid
  m + log(sum(wt * exp(logterm - m)))
}
