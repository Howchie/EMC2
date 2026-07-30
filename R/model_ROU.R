# ============================================================================
# ROU: the racing Ornstein-Uhlenbeck process (leaky accumulator)
#
# Each accumulator is  dX = (v - k X) dt + s dW,  X(0) ~ U(0, A),  absorbed at
# b = B + A.  There is no closed form for the first-passage density, so the
# likelihood is a Fokker-Planck solve (src/fpe_solver.h) cached per distinct
# parameter tuple (src/fpe_race.h).
#
# The R-side dROU/pROU below call the SAME C++ cache as the sampled likelihood,
# so make_data()/predict() and the fit cannot disagree.
# ============================================================================

# Resolution knobs, exposed as options so accuracy can be swept against a real
# fit without a recompile.  Defaults come from the step-10 sweep
# (race_ou_integration_plan.md section 8) and MUST match FPE_Grid in
# src/fpe_race.h, which is what the sampled likelihood uses.
#
# nx is the binding constraint: once the time grid is graded (tgrade = 32),
# halving dt buys 7-16% while raising nx buys a factor of two.  tgrade beyond 32
# buys nothing measurable.
.rou_grid <- function() {
  list(
    nx = getOption("emc2.fpe_nx", 512L),
    dt_target = getOption("emc2.fpe_dt", 4e-3),
    grade = getOption("emc2.fpe_grade", 8),
    tgrade = getOption("emc2.fpe_tgrade", 32)
  )
}

# Boundary form codes, matching FPE_BoundaryKind in src/fpe_models.h.
.ROU_BND <- c(fixed = 0L, weibull = 1L, exponential = 2L, linear_additive = 3L,
              linear_multiplicative = 4L)

# The suffix on c_name selects the form in resolve_race_model_adapter(); the
# shape parameters themselves travel as ordinary optional columns.
.ROU_SUFFIX <- c(fixed = "", weibull = "_BWEIB", exponential = "_BEXP",
                 linear_additive = "_BLIN_ADD", linear_multiplicative = "_BLIN_MULT")

.rou_cols <- function(pars) {
  n <- nrow(pars)
  cn <- dimnames(pars)[[2]]
  # The boundary form is inferred from which shape columns are present, so the
  # R-side dfun/pfun/rfun agree with the C++ likelihood without needing to be
  # told the model variant separately.
  bkind <- if (!("Binf" %in% cn) || !("tau" %in% cn)) .ROU_BND[["fixed"]]
           else if ("pw" %in% cn) .ROU_BND[["weibull"]]
           else .ROU_BND[["exponential"]]
  list(
    v = pars[, "v"],
    k = if ("k" %in% cn) pars[, "k"] else rep(0, n),
    B = pars[, "B"],
    A = if ("A" %in% cn) pars[, "A"] else rep(0, n),
    t0 = pars[, "t0"],
    s = if ("s" %in% cn) pars[, "s"] else rep(1, n),
    bkind = bkind,
    Binf = if (bkind == 0L) numeric(0) else pars[, "Binf"],
    tau = if (bkind == 0L) numeric(0) else pars[, "tau"],
    pw = if (bkind == 1L) pars[, "pw"] else numeric(0)
  )
}

# .rou_cols() cannot distinguish "exponential" from "linear_additive" /
# "linear_multiplicative" -- they take the same columns -- so a linear-collapse
# model has to say so explicitly. Set by the model's dfun/pfun/rfun closures,
# which do know.
.rou_with_kind <- function(p, kind) {
  if (!is.null(kind) && p$bkind != 0L) p$bkind <- .ROU_BND[[kind]]
  p
}

.rou_pdf_cdf <- function(rt, pars, kind = NULL) {
  p <- .rou_with_kind(.rou_cols(pars), kind)
  g <- .rou_grid()
  # NA rates mark accumulators that are not in the race on this trial; the
  # solver has nothing to say about them, so keep them out and zero them after.
  bad <- is.na(p$v) | is.na(rt) | !is.finite(rt)
  if (any(bad)) {
    # Substitute a valid dummy tuple so the grouping pass never sees NA; these
    # rows share one key, so they cost at most one extra solve, and their
    # results are discarded below.
    p$v[bad] <- 1; p$k[bad] <- 0; p$B[bad] <- 1
    p$A[bad] <- 0; p$t0[bad] <- 0; p$s[bad] <- 1
    rt[bad] <- 1
    if (p$bkind != 0L) {
      p$Binf[bad] <- 0.5; p$tau[bad] <- 1
      if (length(p$pw)) p$pw[bad] <- 1
    }
  }
  out <- rou_pdf_cdf_vec(rt, p$v, p$k, p$B, p$A, p$t0, p$s,
                         as.integer(g$nx), g$dt_target, g$grade, g$tgrade,
                         p$bkind, p$Binf, p$tau, p$pw)
  if (any(bad)) {
    out$pdf[bad] <- 0
    out$cdf[bad] <- 0
  }
  out
}

dROU <- function(rt, pars, kind = NULL) .rou_pdf_cdf(rt, pars, kind)$pdf

pROU <- function(rt, pars, kind = NULL) .rou_pdf_cdf(rt, pars, kind)$cdf

#### random

rROU <- function(lR, pars, ok = rep(TRUE, nrow(pars)), kind = NULL,
                 dt = getOption("emc2.rou_sim_dt", 1e-3),
                 t_max = getOption("emc2.rou_sim_tmax", 30)) {
  p <- .rou_with_kind(.rou_cols(pars), kind)
  nr <- length(levels(lR))
  bad <- rep(NA, length(lR) / nr)
  out <- data.frame(R = bad, rt = bad)

  dt_mat <- matrix(Inf, nrow = nr, ncol = nrow(pars) / nr)
  ok2 <- ok & !is.na(p$v)
  if (any(ok2)) {
    dt_mat[ok2] <- rou_hit_times_vec(p$v[ok2], p$k[ok2], p$B[ok2], p$A[ok2],
                                     p$s[ok2], dt, t_max, p$bkind,
                                     if (p$bkind == 0L) numeric(0) else p$Binf[ok2],
                                     if (p$bkind == 0L) numeric(0) else p$tau[ok2],
                                     if (length(p$pw)) p$pw[ok2] else numeric(0))
  }
  bad_col <- apply(dt_mat, 2, function(x) all(is.infinite(x)))
  R <- max.col(-t(dt_mat), ties.method = "first")
  pick <- cbind(R, seq_len(ncol(dt_mat)))
  rt <- matrix(p$t0, nrow = nr)[pick] + dt_mat[pick]
  R <- factor(levels(lR)[R], levels = levels(lR))
  R[bad_col] <- NA
  rt[bad_col] <- Inf
  out$R <- factor(levels(lR)[R], levels = levels(lR))
  out$rt <- rt
  out <- .apply_timed_guess_winner(out, levels(lR))
  out
}

# Pure R reference simulator fallback for .rfun_ROU when emc2.cpp_rfun is FALSE
.rfun_ROU_R <- function(lR, pars, ok = rep(TRUE, nrow(pars)), kind = NULL,
                        dt = getOption("emc2.rou_sim_dt", 1e-3),
                        t_max = getOption("emc2.rou_sim_tmax", 30)) {
  p <- .rou_with_kind(.rou_cols(pars), kind)
  nr <- length(levels(lR))
  bad <- rep(NA, length(lR) / nr)
  out <- data.frame(R = bad, rt = bad)

  n_trials <- nrow(pars) / nr
  dt_mat <- matrix(Inf, nrow = nr, ncol = n_trials)
  ok2 <- ok & !is.na(p$v)

  # Step single OU trajectory in R for trial i
  .sim_one <- function(v, k, B, A, s, bkind, Binf, tau, pw) {
    if (!is.finite(v) || !is.finite(B) || !is.finite(s) || s <= 0) return(Inf)
    kk <- if (is.finite(k) && k > 0) k else 0
    AA <- if (is.finite(A) && A > 0) A else 0
    X <- if (AA > 0) runif(1, 0, AA) else 0
    b0 <- B + AA
    b_t <- function(t) {
      if (bkind == 0L) return(b0)
      binf <- abs(Binf)
      if (bkind == 1L) binf + (b0 - binf) * exp(-(t / tau)^pw)
      else if (bkind == 2L) binf + (b0 - binf) * exp(-t / tau)
      else if (bkind == 3L) binf + (b0 - binf) * pmax(0, 1 - t / tau)
      else if (bkind == 4L) binf + (b0 - binf) / (1 + t / tau)
      else b0
    }
    b <- b_t(0)
    if (X >= b) return(0)

    phi <- exp(-kk * dt)
    m1 <- -expm1(-kk * dt)
    m2 <- -expm1(-2 * kk * dt)
    drift_gain <- if (kk > 1e-10) (m1 / kk) else dt
    var_val <- if (kk > 1e-10) (s * s * m2 / (2 * kk)) else (s * s * dt)
    sd_val <- sqrt(max(var_val, 0))
    inv_2var_bb <- 2 / (s * s * dt)

    t <- 0
    while (t < t_max) {
      X1 <- X * phi + v * drift_gain + sd_val * rnorm(1)
      b1 <- b_t(t + dt)
      t <- t + dt
      if (X1 >= b1) {
        d0 <- b - X; d1 <- b1 - X1
        frac <- d0 / max(d0 - d1, 1e-300)
        return(t - dt + min(max(frac, 0), 1) * dt)
      }
      pc <- exp(-(b - X) * (b1 - X1) * inv_2var_bb)
      if (runif(1) < pc) return(t - dt + runif(1) * dt)
      X <- X1
      b <- b1
    }
    Inf
  }

  if (any(ok2)) {
    idx <- which(ok2)
    hits <- numeric(length(idx))
    for (j in seq_along(idx)) {
      i <- idx[j]
      hits[j] <- .sim_one(p$v[i], p$k[i], p$B[i], p$A[i], p$s[i], p$bkind,
                          if (length(p$Binf) >= i) p$Binf[i] else 0,
                          if (length(p$tau) >= i) p$tau[i] else 0,
                          if (length(p$pw) >= i) p$pw[i] else 0)
    }
    dt_mat[ok2] <- hits
  }

  bad_col <- apply(dt_mat, 2, function(x) all(is.infinite(x)))
  R <- max.col(-t(dt_mat), ties.method = "first")
  pick <- cbind(R, seq_len(ncol(dt_mat)))
  rt <- matrix(p$t0, nrow = nr)[pick] + dt_mat[pick]
  R <- factor(levels(lR)[R], levels = levels(lR))
  R[bad_col] <- NA
  rt[bad_col] <- Inf
  out$R <- factor(levels(lR)[R], levels = levels(lR))
  out$rt <- rt
  out <- .apply_timed_guess_winner(out, levels(lR))
  out
}


#' The Racing Ornstein-Uhlenbeck Model (Leaky Accumulator)
#'
#' Model file to estimate a race between leaky accumulators, each an
#' Ornstein-Uhlenbeck process
#' \deqn{dX = (v - k X) dt + s dW,\quad X(0) \sim U(0, A),}
#' absorbed at \eqn{b = B + A}.
#'
#' Model files are almost exclusively used in `design()`.
#'
#' @details
#'
#' Default values are used for all parameters that are not explicitly listed in the `formula`
#' argument of `design()`. They can also be accessed with `ROU()$p_types`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default**   | **Mapping**          | **Interpretation**                                                |
#' |-----------|-----------|---------------|-----------|------------------|---------------------------------------------------------------|
#' | *v*       | log       | \[0, Inf\]      | log(1)    |                  | Accumulation rate (stimulus input)                             |
#' | *k*       | log       | \[0, Inf\]      | log(0)    |                  | Leak (1/s); *k* = 0 is the Wiener race (RDM)                   |
#' | *A*       | log       | \[0, Inf\]      | log(0)    |                  | Between-trial variation (range) in start point                 |
#' | *B*       | log       | \[0, Inf\]      | log(1)    | *b* = *B* + *A*      | Distance from *A* to *b* (response threshold)                  |
#' | *t0*      | log       | \[0, Inf\]      | log(0)    |                  | Non-decision time                                             |
#' | *s*       | log       | \[0, Inf\]      | log(1)    |                  | Within-trial standard deviation of the diffusion              |
#' | *pContaminant* | probit | \[0, 1\] | qnorm(0) | | Optional contamination probability handled by the data pipeline |
#'
#' The leak *k* has units of 1/time, so it is *not* absorbed by the scaling
#' convention that fixes *s* = 1; with *s* fixed, (*v*, *k*, *B*, *A*) are
#' identified. The parameterization is deliberately (*v*, *k*) rather than the
#' solver's (\eqn{\lambda}, \eqn{\theta}): the asymptote \eqn{\theta = v/k}
#' diverges as the leak vanishes, whereas the drift \eqn{v - kX} is perfectly
#' regular there. As a result *k* = 0 is an ordinary interior value at which the
#' model **is** the racing diffusion model (RDM), reached through the same
#' partial differential equation rather than by dispatching to the Wald
#' formulae. That makes `ROU` with `constants = c(k = log(0))` an independent
#' check on `RDM`, which it would not be if the code branched.
#'
#' The parameterization *b* = *B* + *A* ensures that the response threshold is
#' always higher than the between trial variation in start point. `A` is
#' non-negative, so the start point is always in \[0, *A*\].
#'
#' **This is not the leaky competing accumulator.** It is the leaky accumulator
#' of Smith and Ratcliff (2004, appendix): the accumulators race independently,
#' with no lateral inhibition, and activation is *not* rectified at zero — the
#' process may go negative and only upward threshold crossings count. Usher and
#' McClelland's (2001) LCA adds both of those features.
#'
#' Because there is no closed-form first-passage density, each distinct
#' parameter tuple costs one Fokker-Planck solve (about 1.7 ms at the default
#' resolution, roughly 45 times an analytic race kernel). The number of
#' solves is set by the number of distinct parameter rows in the design, not by
#' the number of trials, so cost is bounded by the design rather than the data.
#' A trend or covariate model makes every trial distinct and is correspondingly
#' much slower — permitted, but priced. Accuracy can be traded against speed
#' with `options(emc2.fpe_nx = )`, `options(emc2.fpe_dt = )`,
#' `options(emc2.fpe_grade = )` and `options(emc2.fpe_tgrade = )`.
#'
#' Smith, P. L., & Ratcliff, R. (2004). Psychology and neurobiology of simple
#' decisions. *Trends in Neurosciences, 27*(3), 161-168.
#'
#' Usher, M., & McClelland, J. L. (2001). The time course of perceptual choice:
#' the leaky, competing accumulator model. *Psychological Review, 108*(3), 550-592.
#'
#' @param boundary_collapse Character; the form of the decision boundary.
#'   `"fixed"` (the default) is a constant threshold. The collapsing forms start
#'   at \eqn{b_0 = B + A} and add `Binf` (the asymptotic boundary, measured from
#'   *zero* rather than from the top of the start-point range) and `tau` (the
#'   time scale of the collapse), plus `pw` (a shape exponent) for `"weibull"`:
#'   \itemize{
#'     \item `"exponential"`: \eqn{b(t) = b_\infty + (b_0 - b_\infty) e^{-t/\tau}}
#'     \item `"linear_additive"`: \eqn{b(t)} falls linearly from \eqn{b_0} to
#'       \eqn{b_\infty} over \eqn{[0, \tau]}, then holds (additive urgency)
#'     \item `"linear_multiplicative"`: \eqn{b(t) = b_\infty + (b_0 - b_\infty) / (1 + t/\tau)}
#'       (multiplicative urgency)
#'     \item `"weibull"`: \eqn{b(t) = b_\infty + (b_0 - b_\infty) e^{-(t/\tau)^{pw}}}
#'   }
#'   Because `Binf` is measured from zero, the boundary is free to fall into the
#'   start-point range \eqn{[0, A]} and, at `Binf = 0`, to collapse all the way
#'   to the start point -- a forced response. `Binf = 0` is an ordinary interior
#'   value, not a limit: the solver's domain is \eqn{[x_{lo}, b(t)]} with
#'   \eqn{x_{lo} < 0}, so it stays well conditioned there. A moving boundary forfeits the one-time factorisation of
#'   the solver's linear operator, so it is slower per solve; setting `Binf`
#'   equal to \eqn{B + A} makes the boundary constant again and recovers the
#'   fixed-bound cost exactly.
#'
#' @return A list defining the cognitive model
#' @examples
#' # A leaky accumulator design with the same structure as the RDM example.
#' ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
#' matchfun=function(d)d$S==d$lR
#' design_ROU <- design(data = forstmann,model=ROU,matchfun=matchfun,
#'                      formula=list(v~lM,k~1,B~E+lR,A~1,t0~1),
#'                      contrasts=list(v=list(lM=ADmat)),constants=c(s=log(1)))
#'
#' # The same race with an exponentially collapsing threshold.
#' design_ROUc <- design(data = forstmann,
#'                       model=function() ROU(boundary_collapse="exponential"),
#'                       matchfun=matchfun,
#'                       formula=list(v~lM,k~1,B~E+lR,A~1,t0~1,Binf~1,tau~1),
#'                       contrasts=list(v=list(lM=ADmat)),constants=c(s=log(1)))
#' @export

ROU <- function(boundary_collapse = c("fixed", "exponential", "linear_additive",
                                      "linear_multiplicative", "weibull")) {
  boundary_collapse <- match.arg(boundary_collapse)
  kind <- boundary_collapse

  p_types <- c("v" = log(1), "k" = log(0), "B" = log(1), "A" = log(0),
               "t0" = log(0), "s" = log(1))
  transform <- c(v = "exp", k = "exp", B = "exp", A = "exp", t0 = "exp",
                 s = "exp")
  minmax <- cbind(v = c(1e-3, Inf), k = c(0, Inf), B = c(0, Inf),
                  A = c(1e-4, Inf), t0 = c(0.05, Inf), s = c(0, Inf))
  exception <- c(A = 0, v = 0, k = 0)

  # Optional columns for the collapsing forms, in the order declared after
  # N_REQ in emc2col::rou (src/col_registry.h) and BEFORE pContaminant, exactly
  # as BAwL orders mG/mK/omega.  Only the columns the selected form actually
  # uses are added, so a fixed-bound ROU is unchanged by collapse existing.
  #
  # Binf is the asymptotic boundary measured from ZERO, not from the top of the
  # start-point range the way B is.  The bound may therefore descend into
  # [0, A]: a bound that meets the start point is a forced response, which the
  # model should be able to express.
  if (kind != "fixed") {
    p_types <- c(p_types, Binf = log(0.5), tau = log(1))
    transform <- c(transform, Binf = "exp", tau = "exp")
    # Binf = 0 is a legal interior value, not a limit to be approached: the
    # solver's domain is [x_lo, b(t)] with x_lo < 0, so b(t) -> 0 leaves it
    # perfectly well conditioned.  Verified smooth and monotone through
    # 1e-3 -> 1e-6 -> 0 for all three forms.  A bound that reaches zero has
    # met the start point, i.e. a forced response.
    minmax <- cbind(minmax, Binf = c(0, Inf), tau = c(1e-3, Inf))
    exception <- c(exception, Binf = 0)
    if (kind == "weibull") {
      p_types <- c(p_types, pw = log(1))
      transform <- c(transform, pw = "exp")
      minmax <- cbind(minmax, pw = c(1e-3, Inf))
    }
  }

  p_types <- c(p_types, pContaminant = qnorm(0))
  transform <- c(transform, pContaminant = "pnorm")
  minmax <- cbind(minmax, pContaminant = c(0.001, 0.999))
  exception <- c(exception, pContaminant = 0)

  list(
    type = "RACE",
    c_name = paste0("ROU", .ROU_SUFFIX[[kind]]),
    p_types = p_types,
    p_types_canonical = names(p_types)[names(p_types) != "pContaminant"],
    transform = list(func = transform),
    bound = list(minmax = minmax, exception = exception),
    # Trial dependent parameter transform
    Ttransform = function(pars, dadm) {
      pars <- cbind(pars, b = pars[, "B"] + pars[, "A"])
      pars
    },
    # Random function for racing accumulators
    rfun = function(data = NULL, pars) .rfun_ROU(data$lR, pars, ok = attr(pars, "ok"),
                                                 kind = kind),
    # Density function (PDF) for single accumulator
    dfun = function(rt, pars) dROU(rt, pars, kind = kind),
    # Probability function (CDF) for single accumulator
    pfun = function(rt, pars) pROU(rt, pars, kind = kind),
    # Race likelihood combining pfun and dfun
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      log_likelihood_race_missing(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}
