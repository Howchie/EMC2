# ============================================================================
# ROU: the racing Ornstein-Uhlenbeck process (leaky accumulator)
#
# Each accumulator is  dX = (v - k X) dt + s dW,  X(0) ~ U(0, A),  absorbed at
# b = B + A.  There is no closed form for the first-passage density, so the
# likelihood is a Fokker-Planck solve (src/fpe_solver.h) cached per distinct
# parameter tuple (src/fpe_race.h).
# ============================================================================

# Resolution knobs are exposed as options so accuracy can be adjusted without
# recompiling the package. They must match the grid used by the C++ likelihood.
.rou_grid <- function() {
  list(
    nx = getOption("emc2.fpe_nx", 512L),
    dt_target = getOption("emc2.fpe_dt", 4e-3),
    grade = getOption("emc2.fpe_grade", 8),
    tgrade = getOption("emc2.fpe_tgrade", 32)
  )
}

# Boundary-form codes, matching the C++ solver.
.ROU_BND <- c(fixed = 0L, weibull = 1L, exponential = 2L, linear_additive = 3L,
              linear_multiplicative = 4L)

# The suffix on c_name selects the boundary form in the C++ adapter.
.ROU_SUFFIX <- c(fixed = "", weibull = "_BWEIB", exponential = "_BEXP",
                 linear_additive = "_BLIN_ADD", linear_multiplicative = "_BLIN_MULT")

# Parameterisation codes, matching fperace::ROU_PAR_* in src/fpe_race.h, and the
# c_name infix that tells the C++ adapter which one is in play. The infix comes
# before the boundary suffix, so ROUCURV_BEXP is a curvature-parameterised model
# with an exponentially collapsing bound.
.ROU_PAR <- c(rate = 0L, curvature = 1L, equilibrium = 2L)
.ROU_PAR_INFIX <- c(rate = "", curvature = "CURV", equilibrium = "EQ")

# The three columns each parameterisation puts in place of (v, k, s).
.ROU_PAR_COLS <- list(rate = c("v", "k", "s"),
                      curvature = c("tstar", "k", "s"),
                      equilibrium = c("tk", "theta", "chi"))

.rou_cols <- function(pars, par = "rate") {
  n <- nrow(pars)
  cn <- dimnames(pars)[[2]]
  # The boundary form is inferred from which shape columns are present, so the
  # R-side dfun/pfun/rfun agree with the C++ likelihood without needing to be
  # told the model variant separately.
  bkind <- if (!("Binf" %in% cn) || !("tau" %in% cn)) .ROU_BND[["fixed"]]
           else if ("pw" %in% cn) .ROU_BND[["weibull"]]
           else .ROU_BND[["exponential"]]
  B <- pars[, "B"]
  A <- if ("A" %in% cn) pars[, "A"] else rep(0, n)
  if (par == "rate") {
    v <- pars[, "v"]
    k <- if ("k" %in% cn) pars[, "k"] else rep(0, n)
    s <- if ("s" %in% cn) pars[, "s"] else rep(1, n)
  } else {
    # The map lives in C++ (src/fpe_race.h) and is shared with the likelihood
    # kernels, so the R density cannot drift from the sampled one.
    nm <- .ROU_PAR_COLS[[par]]
    m <- rou_to_rate_vec(.ROU_PAR[[par]], pars[, nm[1]], pars[, nm[2]],
                         pars[, nm[3]], B, A)
    v <- m$v; k <- m$k; s <- m$s
  }
  list(
    v = v,
    k = k,
    B = B,
    A = A,
    t0 = pars[, "t0"],
    s = s,
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

.rou_pdf_cdf <- function(rt, pars, kind = NULL, par = "rate") {
  p <- .rou_with_kind(.rou_cols(pars, par), kind)
  g <- .rou_grid()
  # NA rates mark accumulators that are not in the race on this trial; the
  # solver has nothing to say about them, so keep them out and zero them after.
  # Non-finite mapped rates or shifts cannot produce a valid key.  Keep them
  # separate from a positive-infinite query so only a valid pROU(Inf) row gets
  # the proper OU limit below.
  bad_par <- !is.finite(p$v) | !is.finite(p$t0)
  bad_rt <- is.na(rt) | !is.finite(rt)
  bad <- bad_par | bad_rt
  rt_pos_inf <- is.infinite(rt) & rt > 0
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
    # With s > 0 and a finite upper boundary, the OU hits eventually with
    # probability one.  Preserve that limit for a direct pROU(Inf) query;
    # invalid rows and -Inf/NA remain at zero.
    out$cdf[bad] <- ifelse(!bad_par[bad] & rt_pos_inf[bad],
                           1, 0)
  }
  out
}

dROU <- function(rt, pars, kind = NULL, par = "rate")
  .rou_pdf_cdf(rt, pars, kind, par)$pdf

pROU <- function(rt, pars, kind = NULL, par = "rate")
  .rou_pdf_cdf(rt, pars, kind, par)$cdf

rROU <- function(lR, pars, ok = rep(TRUE, nrow(pars)), kind = NULL,
                 par = "rate",
                 dt = getOption("emc2.rou_sim_dt", 1e-3),
                 t_max = getOption("emc2.rou_sim_tmax", 30)) {
  p <- .rou_with_kind(.rou_cols(pars, par), kind)
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
  bad_col <- colSums(!is.infinite(dt_mat)) == 0L
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
                        par = "rate",
                        dt = getOption("emc2.rou_sim_dt", 1e-3),
                        t_max = getOption("emc2.rou_sim_tmax", 30)) {
  p <- .rou_with_kind(.rou_cols(pars, par), kind)
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

  bad_col <- colSums(!is.infinite(dt_mat)) == 0L
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
#'
#' @details
#'
#' Default values are used for all parameters that are not explicitly listed in the `formula`
#' argument of `design()`. They can also be accessed with `ROU()$p_types`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default**   | **Mapping**          | **Interpretation**                                                |
#' |-----------|-----------|---------------|-----------|------------------|---------------------------------------------------------------|
#' | *v*       | log       | \[0, Inf\]      | log(1)    |                  | Accumulation rate (stimulus input)                             |
#' | *k*       | log       | \[0, Inf\]      | log(0)    |                  | Leak, a rate in units of 1/time (**not** divided by *s*); *k* = 0 is the Wiener race (RDM) |
#' | *A*       | log       | \[0, Inf\]      | log(0)    |                  | Between-trial variation (range) in start point                 |
#' | *B*       | log       | \[0, Inf\]      | log(1)    | *b* = *B* + *A*      | Distance from *A* to *b* (response threshold)                  |
#' | *t0*      | log       | \[0, Inf\]      | log(0)    |                  | Non-decision time                                             |
#' | *s*       | log       | \[0, Inf\]      | log(1)    |                  | Within-trial standard deviation of the diffusion              |
#' | *pContaminant* | probit | \[0, 1\] | qnorm(0) | | Optional *omission* contaminant probability: mass at `rt = Inf` only, handled by the data pipeline |
#' | *pGuess* | probit | \[0, 1\] | qnorm(0) | | Optional uniform *guess* (outlier) probability, mixed into observed RT densities over the guess window |
#'
#' # Parameterizations
#'
#' `parameterization` changes which three parameters stand in for
#' \eqn{(v, k, s)}. The charts share one solver but cover restricted regions of
#' the SDE family; `B`, `A`, `t0` and the collapsing-bound parameters are common
#' to all three.
#'
#' `"rate"` (the default) is \eqn{(v, k, s)} as above.
#'
#' `"curvature"` replaces them with \eqn{(t_\star, k, s)}. `tstar` is the
#' reference crossing time of the deterministic path from \eqn{x = 0} to
#' \eqn{b = B + A}; `k` and `s` retain their physical units:
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Interpretation** |
#' |-----------|-----------|---------------|-----------|------------------------------------|
#' | *tstar*   | log       | \[0, Inf\]    | log(1)    | Time at which the deterministic mean path reaches *b* |
#' | *k*       | log       | \[0, Inf\]    | log(0)    | Physical leak rate; *k* = 0 is the Wiener race |
#' | *s*       | log       | \[0, Inf\]    | log(1)    | Diffusion standard deviation |
#'
#' The map is \eqn{v = bk/(1-e^{-k t_\star})}, with the smooth limit
#' \eqn{v = b/t_\star} at *k* = 0. This is a deterministic-crossing chart, so
#' `tstar ~ condition`, `k ~ 1`, and `s ~ 1` have their literal meanings.
#'
#' `"equilibrium"` replaces them with \eqn{(t_k, \theta, \chi)}:
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Interpretation** |
#' |-----------|-----------|---------------|-----------|------------------------------------|
#' | *tk*      | log       | \[0.001, Inf\]  | log(1)    | Leak time constant \eqn{t_k = 1/k} (the boundary-collapse time constant keeps the name *tau*) |
#' | *theta*   | log       | \[0.001, 20\]   | log(1)    | Physical OU equilibrium, \eqn{v/k} |
#' | *chi*     | log       | \[0.001, 10\]   | log(1)    | Physical noise over one relaxation, \eqn{s\sqrt{t_k}} |
#'
#' The map is \eqn{k = 1/t_k}, \eqn{v = \theta/t_k}, and
#' \eqn{s = \chi/\sqrt{t_k}}. Thus the equilibrium is below the upper boundary
#' when \eqn{\theta < B + A}; responses in that regime are noise-driven escapes
#' and effective omissions in a finite response window, not an intrinsic
#' never-response probability. The start point remains \eqn{U(0,A)}.
#'
#' # Scale identification
#'
#' Scaling the state shows the dynamics depend only on
#' \eqn{(v/s, k, B/s, A/s)} — four quantities from five parameters. The rate
#' parameterization resolves the redundancy by fixing `s` (usually
#' `constants = c(s = log(1))`). The alternatives leave the same state-scale
#' redundancy, so **fix `B` rather than `s`**, with a reference/intercept
#' `constants = c(B = log(1))`. In the equilibrium chart `theta` and `chi` are
#' physical-unit parameters; fixing the threshold reference keeps their scale
#' common when condition-specific threshold effects are added.
#' A bound manipulation is represented directly by `tstar ~ E` in the curvature
#' chart; `theta` and `tk` in the equilibrium chart change equilibrium and
#' relaxation, respectively.
#'
#' This is the independent leaky accumulator, not the leaky competing
#' accumulator (LCA): accumulators do not inhibit one another, and the process
#' is not rectified at zero. See [BOU()] for the two-boundary OU diffusion,
#' whose leak is defined relative to the starting point rather than zero.
#'
#'
#' Response-time compression is disabled for this model because the solver
#' evaluates the continuous-time density directly; binning would add flooring
#' bias without reducing the solve cost.
#'
#' Ratcliff, R., & Smith, P. L. (2004). A comparison of sequential sampling
#' models for two-choice reaction time. *Psychological Review, 111*(2), 333-367.
#'
#' Smith, P. L. (2000). Stochastic dynamic models of response time and accuracy:
#' A foundational primer. *Journal of Mathematical Psychology, 44*(3), 408-463.
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
#'   `Binf` is measured from zero, so the boundary may collapse into the
#'   start-point range. Setting `Binf = B + A` gives the fixed-boundary model.
#'
#' @param parameterization Character; which three parameters stand in for
#'   \eqn{(v, k, s)}. `"rate"` (the default) estimates them directly.
#'   `"curvature"` estimates \eqn{(t_\star, k, s)} — a reference crossing
#'   time, physical leak and diffusion — with \eqn{v = bk/(1-e^{-k t_\star})}
#'   and the smooth \eqn{k=0} limit.
#'   `"equilibrium"` estimates \eqn{(t_k, \theta, \chi)} — relaxation time,
#'   physical equilibrium \eqn{v/k}, and noise per relaxation \eqn{s\sqrt{t_k}}.
#'   The alternatives fix the state scale with `B` rather than with `s`; see
#'   Details.
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
#'
#' # The curvature parameterization: stimulus and speed emphasis move the
#' # reference crossing time, while physical leak and diffusion can be shared.
#' design_ROUcurv <- design(data = forstmann,
#'                          model=function() ROU(parameterization="curvature"),
#'                          matchfun=matchfun,
#'                          formula=list(tstar~lM+E,k~1,s~1,A~1,t0~1),
#'                          contrasts=list(tstar=list(lM=ADmat)),
#'                          constants=c(B=log(1)))
#'
#' # The equilibrium parameterization, for a deadline design in which theta is free
#' # to sit below the bound. Give the data a UC column to censor at it.
#' design_ROUeq <- design(data = forstmann,
#'                        model=function() ROU(parameterization="equilibrium"),
#'                        matchfun=matchfun,
#'                        formula=list(tk~1,theta~lM,chi~1,A~1,t0~1),
#'                        contrasts=list(theta=list(lM=ADmat)),
#'                        constants=c(B=log(1)))
#' @export

ROU <- function(boundary_collapse = c("fixed", "exponential", "linear_additive",
                                      "linear_multiplicative", "weibull"),
                parameterization = c("rate", "curvature", "equilibrium")) {
  boundary_collapse <- match.arg(boundary_collapse)
  parameterization <- match.arg(parameterization)
  kind <- boundary_collapse
  par <- parameterization

  if (par == "rate") {
    p_types <- c("v" = log(1), "k" = log(0), "B" = log(1), "A" = log(0),
                 "t0" = log(0), "s" = log(1))
    transform <- c(v = "exp", k = "exp", B = "exp", A = "exp", t0 = "exp",
                   s = "exp")
    minmax <- cbind(v = c(1e-3, Inf), k = c(0, Inf), B = c(0, Inf),
                    A = c(1e-4, Inf), t0 = c(0.05, Inf), s = c(0, Inf))
    exception <- c(A = 0, v = 0, k = 0)
  } else if (par == "curvature") {
    # Defaults map to the rate default: tstar = 1, k = 0, s = 1, B = 1, A = 0.
    p_types <- c("tstar" = log(1), "k" = log(0), "s" = log(1), "B" = log(1),
                 "A" = log(0), "t0" = log(0))
    transform <- c(tstar = "exp", k = "exp", s = "exp", B = "exp", A = "exp",
                   t0 = "exp")
    minmax <- cbind(tstar = c(0.05, Inf), k = c(0, Inf), s = c(0, Inf),
                    B = c(0, Inf), A = c(1e-4, Inf), t0 = c(0.05, Inf))
    exception <- c(A = 0, k = 0)
  } else {
    p_types <- c("tk" = log(1), "theta" = log(1), "chi" = log(1), "B" = log(1),
                 "A" = log(0), "t0" = log(0))
    transform <- c(tk = "exp", theta = "exp", chi = "exp", B = "exp", A = "exp",
                   t0 = "exp")
    minmax <- cbind(tk = c(0.01, Inf), theta = c(1e-3, 20), chi = c(1e-3, 10),
                    B = c(0, Inf), A = c(1e-4, Inf), t0 = c(0.05, Inf))
    exception <- c(A = 0)
  }

  # Collapsing forms append only the columns they use. Binf is measured from
  # zero, while B is the distance from the top of the start-point range.
  if (kind != "fixed") {
    p_types <- c(p_types, Binf = log(0.5), tau = log(1))
    transform <- c(transform, Binf = "exp", tau = "exp")
    # Binf = 0 is allowed; the solver handles a boundary that reaches the
    # start-point range.
    minmax <- cbind(minmax, Binf = c(0, Inf), tau = c(1e-3, Inf))
    exception <- c(exception, Binf = 0)
    if (kind == "weibull") {
      p_types <- c(p_types, pw = log(1))
      transform <- c(transform, pw = "exp")
      minmax <- cbind(minmax, pw = c(1e-3, Inf))
    }
  }

  # pContaminant (omission) and pGuess (uniform outlier); see add_nuisance_pars().
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- .nuis$p_types; transform <- .nuis$transform
  minmax <- .nuis$minmax; exception <- .nuis$exception

  list(
    type = "RACE",
    c_name = paste0("ROU", .ROU_PAR_INFIX[[par]], .ROU_SUFFIX[[kind]]),
    # The density is evaluated by a cached Fokker--Planck solve.
    compress_ok = FALSE,
    p_types = p_types,
    p_types_canonical = setdiff(names(p_types), .nuisance_par_names),
    transform = list(func = transform),
    bound = list(minmax = minmax, exception = exception),
    # Trial dependent parameter transform
    # b = B + A is the threshold absorbed by the kernel.
    Ttransform = function(pars, dadm) {
      pars <- cbind(pars, b = pars[, "B"] + pars[, "A"])
      pars
    },
    # Random function for racing accumulators
    rfun = function(data = NULL, pars) .rfun_ROU(data$lR, pars, ok = attr(pars, "ok"),
                                                 kind = kind, par = par),
    # Density function (PDF) for single accumulator
    dfun = function(rt, pars) dROU(rt, pars, kind = kind, par = par),
    # Probability function (CDF) for single accumulator
    pfun = function(rt, pars) pROU(rt, pars, kind = kind, par = par),
    # Race likelihood combining pfun and dfun
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      log_likelihood_race_missing(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}
