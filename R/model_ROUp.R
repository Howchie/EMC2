# ============================================================================
# ROUp: the racing Ornstein-Uhlenbeck process (leaky accumulator)
#
# Each accumulator is  dX = (v - k X) dt + s dW,  X(0) ~ U(0, A),  absorbed at
# b = B + A.  There is no closed form for the first-passage density, so the
# likelihood is a Fokker-Planck solve (src/fpe_solver.h) cached per distinct
# parameter tuple (src/fpe_race.h).
# ============================================================================

# Resolution knobs are exposed as options so accuracy can be adjusted without
# recompiling the package. They must match the grid used by the C++ likelihood.
.roup_grid <- function() {
  list(
    nx = getOption("emc2.fpe_nx", 512L),
    dt_target = getOption("emc2.fpe_dt", 4e-3),
    grade = getOption("emc2.fpe_grade", 8),
    tgrade = getOption("emc2.fpe_tgrade", 32)
  )
}

# Boundary-form codes, matching the C++ solver.
.ROUp_BND <- c(fixed = 0L, weibull = 1L, exponential = 2L, linear_additive = 3L,
              linear_multiplicative = 4L)

# The suffix on c_name selects the boundary form in the C++ adapter.
.ROUp_SUFFIX <- c(fixed = "", weibull = "_BWEIB", exponential = "_BEXP",
                 linear_additive = "_BLIN_ADD", linear_multiplicative = "_BLIN_MULT")

# Parameterisation codes, matching fperace::ROUP_PAR_* in src/fpe_race.h, and the
# c_name infix that tells the C++ adapter which transient chart is in play. The
# infix comes before the boundary suffix, so ROUpAREA_BEXP is the area chart with
# an exponentially collapsing bound.
.ROUp_PAR <- c(rate = 0L, area = 1L)
.ROUp_PAR_INFIX <- c(rate = "", area = "AREA")

# The transient columns supplied by each chart. The sustained channel (v_S,
# tau_S) remains present in both charts.
.ROUp_PAR_COLS <- list(rate = c("v_T", "tau_T"),
                      area = c("E_T", "tau_T"))

.roup_cols <- function(pars, par = "rate") {
  par <- match.arg(par, names(.ROUp_PAR))
  n <- nrow(pars)
  cn <- dimnames(pars)[[2]]
  bkind <- if (!("Binf" %in% cn) || !("tau" %in% cn)) .ROUp_BND[["fixed"]]
           else if ("pw" %in% cn) .ROUp_BND[["weibull"]]
           else .ROUp_BND[["exponential"]]
  B <- pars[, "B"]
  A <- if ("A" %in% cn) pars[, "A"] else rep(0, n)
  v_S <- pars[, "v_S"]
  v_T <- pars[, if (par == "area") "E_T" else "v_T"]
  tau_S <- pars[, "tau_S"]
  tau_T <- pars[, "tau_T"]
  k <- if ("k" %in% cn) pars[, "k"] else rep(0, n)
  s <- if ("s" %in% cn) pars[, "s"] else rep(1, n)
  
  list(
    v_S = v_S,
    v_T = v_T,
    tau_S = tau_S,
    tau_T = tau_T,
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

# .roup_cols() cannot distinguish "exponential" from "linear_additive" /
# "linear_multiplicative" -- they take the same columns -- so a linear-collapse
# model has to say so explicitly. Set by the model's dfun/pfun/rfun closures,
# which do know.
.roup_with_kind <- function(p, kind) {
  if (!is.null(kind) && p$bkind != 0L) p$bkind <- .ROUp_BND[[kind]]
  p
}

.roup_pdf_cdf <- function(rt, pars, kind = NULL, par = "rate") {
  p <- .roup_with_kind(.roup_cols(pars, par), kind)
  g <- .roup_grid()
  bad_par <- !is.finite(p$v_S) | !is.finite(p$t0)
  bad_rt <- is.na(rt) | !is.finite(rt)
  bad <- bad_par | bad_rt
  rt_pos_inf <- is.infinite(rt) & rt > 0
  if (any(bad)) {
    p$v_S[bad] <- 1; p$v_T[bad] <- 0; p$tau_S[bad] <- 1; p$tau_T[bad] <- 1;
    p$k[bad] <- 0; p$B[bad] <- 1; p$A[bad] <- 0; p$t0[bad] <- 0; p$s[bad] <- 1;
    rt[bad] <- 1
    if (p$bkind != 0L) {
      p$Binf[bad] <- 0.5; p$tau[bad] <- 1
      if (length(p$pw)) p$pw[bad] <- 1
    }
  }
  out <- droup_cpp(rt, p$v_S, p$v_T, p$tau_S, p$tau_T, p$k, p$B, p$A, p$t0, p$s,
                   as.integer(g$nx), g$dt_target, g$grade, g$tgrade,
                   p$bkind, p$Binf, p$tau, p$pw, .ROUp_PAR[[par]])
  if (any(bad)) {
    out$pdf[bad] <- 0
    out$cdf[bad] <- ifelse(!bad_par[bad] & rt_pos_inf[bad], 1, 0)
  }
  out
}

dROUp <- function(rt, pars, kind = NULL, par = "rate")
  .roup_pdf_cdf(rt, pars, kind, par)$pdf

pROUp <- function(rt, pars, kind = NULL, par = "rate")
  .roup_pdf_cdf(rt, pars, kind, par)$cdf

rROUp <- function(lR, pars, ok = rep(TRUE, nrow(pars)), kind = NULL,
                 par = "rate",
                 dt = getOption("emc2.roup_sim_dt", 1e-3),
                 t_max = getOption("emc2.roup_sim_tmax", 30)) {
  p <- .roup_with_kind(.roup_cols(pars, par), kind)
  nr <- length(levels(lR))
  bad <- rep(NA, length(lR) / nr)
  out <- data.frame(R = bad, rt = bad)

  dt_mat <- matrix(Inf, nrow = nr, ncol = nrow(pars) / nr)
  ok2 <- ok & !is.na(p$v_S)
  if (any(ok2)) {
    dt_mat[ok2] <- rroup_hit_times_cpp(p$v_S[ok2], p$v_T[ok2], p$tau_S[ok2], p$tau_T[ok2],
                                     p$k[ok2], p$B[ok2], p$A[ok2],
                                     p$s[ok2], dt, t_max, p$bkind,
                                     if (p$bkind == 0L) numeric(0) else p$Binf[ok2],
                                     if (p$bkind == 0L) numeric(0) else p$tau[ok2],
                                     if (length(p$pw)) p$pw[ok2] else numeric(0),
                                     .ROUp_PAR[[par]])
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

# Pure R reference simulator fallback for .rfun_ROUp when emc2.cpp_rfun is FALSE
.rfun_ROUp_R <- function(lR, pars, ok = rep(TRUE, nrow(pars)), kind = NULL,
                        par = "rate",
                        dt = getOption("emc2.roup_sim_dt", 1e-3),
                        t_max = getOption("emc2.roup_sim_tmax", 30)) {
  p <- .roup_with_kind(.roup_cols(pars, par), kind)
  nr <- length(levels(lR))
  bad <- rep(NA, length(lR) / nr)
  out <- data.frame(R = bad, rt = bad)

  n_trials <- nrow(pars) / nr
  dt_mat <- matrix(Inf, nrow = nr, ncol = n_trials)
  ok2 <- ok & !is.na(p$v_S)

  # Step single OU trajectory in R for trial i
  .sim_one <- function(v_S, v_T, tau_S, tau_T, k, B, A, s, bkind, Binf, tau, pw) {
    if (!is.finite(v_S) || !is.finite(v_T) || !is.finite(B) || !is.finite(s) || s <= 0) return(Inf)
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
      tm <- t + 0.5 * dt
      mu_S <- v_S * ifelse(tau_S > 1e-12, 1 - exp(-tm / tau_S), 1)
      mu_T <- v_T * ifelse(tau_T > 1e-12, (tm / tau_T) * exp(-tm / tau_T), 0)
      d_step <- (mu_S + mu_T) * drift_gain

      X1 <- X * phi + d_step + sd_val * rnorm(1)
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
    v_T_rate <- if (par == "area") {
      ifelse(is.finite(p$tau_T) & p$tau_T > 0, p$v_T / p$tau_T, 0)
    } else p$v_T
    for (j in seq_along(idx)) {
      i <- idx[j]
      hits[j] <- .sim_one(p$v_S[i], v_T_rate[i], p$tau_S[i], p$tau_T[i],
                          p$k[i], p$B[i], p$A[i], p$s[i], p$bkind,
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


#' The Racing Ornstein-Uhlenbeck Model with Time-Varying Pulse Drift (Smith 1995)
#'
#' Model file to estimate a race between leaky accumulators with time-varying
#' drift rates driven by sustained and transient visual channels (Smith, 1995).
#' Each accumulator follows an Ornstein-Uhlenbeck process:
#' \deqn{dX = (\mu(t) - k X) dt + s dW,\quad X(0) \sim U(0, A),}
#' absorbed at \eqn{b = B + A}, where the instantaneous drift rate \eqn{\mu(t)} is:
#' \deqn{\mu(t) = \mu_S(t) + \mu_T(t) = v_S (1 - e^{-t/\tau_S}) + v_T \left(\frac{t}{\tau_T}\right) e^{-t/\tau_T}.}
#' In the `area` parameterisation, `v_T` is replaced by \eqn{E_T/\tau_T},
#' where \eqn{E_T = v_T\tau_T} is the integrated transient evidence.
#'
#' @details
#'
#' Default values are used for all parameters that are not explicitly listed in the `formula`
#' argument of `design()`. They can also be accessed with `ROUp()$p_types`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default**   | **Mapping**          | **Interpretation**                                                |
#' |-----------|-----------|---------------|-----------|------------------|---------------------------------------------------------------|
#' | *v_S*     | log       | \[0, Inf\]      | log(1)    |                  | Sustained channel asymptotic drift rate                        |
#' | *v_T*     | log       | \[0, Inf\]      | log(1)    |                  | Transient channel peak drift rate factor                       |
#' | *E_T*     | log       | \[0, Inf\]      | log(1)    |                  | Integrated transient evidence (`area` only)                     |
#' | *tau_S*   | log       | \[0, Inf\]      | log(1)    |                  | Time constant of sustained low-pass filter                     |
#' | *tau_T*   | log       | \[0, Inf\]      | log(1)    |                  | Time constant of transient band-pass filter (peak at \eqn{\tau_T}) |
#' | *k*       | log       | \[0, Inf\]      | log(0)    |                  | Leak rate in units of 1/time; *k* = 0 is the Wiener race       |
#' | *A*       | log       | \[0, Inf\]      | log(0)    |                  | Between-trial variation (range) in start point                 |
#' | *B*       | log       | \[0, Inf\]      | log(1)    | *b* = *B* + *A*      | Distance from *A* to *b* (response threshold)                  |
#' | *t0*      | log       | \[0, Inf\]      | log(0)    |                  | Non-decision time                                             |
#' | *s*       | log       | \[0, Inf\]      | log(1)    |                  | Within-trial standard deviation of the diffusion              |
#' | *pContaminant* | probit | \[0, 1\] | qnorm(0) | | Optional *omission* contaminant probability: mass at `rt = Inf` only, handled by the data pipeline |
#' | *pGuess* | probit | \[0, 1\] | qnorm(0) | | Optional uniform *guess* (outlier) probability, mixed into observed RT densities over the guess window |
#'
#' # Scale identification
#'
#' Scaling the state shows the rate chart depends only on
#' \eqn{(v_S/s, v_T/s, k, B/s, A/s)} — quantities relative to *s*. In the `area`
#' chart, the corresponding transient coordinate is \eqn{E_T/s}; \eqn{\tau_T}
#' is unchanged. Fix `s` (usually `constants = c(s = log(1))`).
#'
#' Response-time compression is disabled for this model because the solver
#' evaluates the continuous-time density directly; binning would add flooring
#' bias without reducing the solve cost.
#'
#' # References
#'
#' Smith, P. L. (1995). Psychophysically principled models of visual simple
#' reaction time. *Psychological Review, 102*(4), 567-593.
#'
#' Smith, P. L. (1998). Attention and luminance in the detection of brief visual
#' stimuli. *Journal of Experimental Psychology: Human Perception and Performance, 24*(1), 105-133.
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
#' @param parameterization Character; transient parameterisation type. `"rate"`
#'   uses `(v_T, tau_T)`. `"area"` uses `(E_T, tau_T)`, where
#'   `E_T = v_T * tau_T`; the sustained `(v_S, tau_S)` channel is retained in
#'   both charts.
#'
#' @return A list defining the cognitive model
#' @examples
#' # A 2-choice pulse-drift race model
#' ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
#' matchfun <- function(d) d$S == d$lR
#' design_ROUp <- design(data = forstmann, model = ROUp, matchfun = matchfun,
#'                       formula = list(v_S ~ lM, v_T ~ 1, tau_S ~ 1, tau_T ~ 1,
#'                                      k ~ 1, B ~ E + lR, A ~ 1, t0 ~ 1),
#'                       contrasts = list(v_S = list(lM = ADmat)),
#'                       constants = c(s = log(1)))
#' @export

ROUp <- function(boundary_collapse = c("fixed", "exponential", "linear_additive",
                                      "linear_multiplicative", "weibull"),
                parameterization = c("rate", "area")) {
  boundary_collapse <- match.arg(boundary_collapse)
  parameterization <- match.arg(parameterization)
  kind <- boundary_collapse
  par <- parameterization

  if (par == "rate") {
    p_types <- c("v_S" = log(1), "v_T" = log(1), "tau_S" = log(1), "tau_T" = log(1),
                 "k" = log(0), "B" = log(1), "A" = log(0),
                 "t0" = log(0), "s" = log(1))
    transform <- c(v_S = "exp", v_T = "exp", tau_S = "exp", tau_T = "exp",
                   k = "exp", B = "exp", A = "exp", t0 = "exp", s = "exp")
    minmax <- cbind(v_S = c(1e-3, Inf), v_T = c(1e-3, Inf),
                    tau_S = c(1e-3, Inf), tau_T = c(1e-3, Inf),
                    k = c(0, Inf), B = c(0, Inf),
                    A = c(1e-4, Inf), t0 = c(0.05, Inf), s = c(0, Inf))
    exception <- c(A = 0, v_S = 0, v_T = 0, tau_S = 0, tau_T = 0, k = 0)
  } else if (par == "area") {
    p_types <- c("v_S" = log(1), "E_T" = log(1), "tau_S" = log(1), "tau_T" = log(1),
                 "k" = log(0), "B" = log(1), "A" = log(0),
                 "t0" = log(0), "s" = log(1))
    transform <- c(v_S = "exp", E_T = "exp", tau_S = "exp", tau_T = "exp",
                   k = "exp", B = "exp", A = "exp", t0 = "exp", s = "exp")
    minmax <- cbind(v_S = c(1e-3, Inf), E_T = c(1e-3, Inf),
                    tau_S = c(1e-3, Inf), tau_T = c(1e-3, Inf),
                    k = c(0, Inf), B = c(0, Inf),
                    A = c(1e-4, Inf), t0 = c(0.05, Inf), s = c(0, Inf))
    # tau_T remains strictly positive in this chart: E_T/tau_T is undefined
    # at a zero time constant, even when E_T is fixed at its zero exception.
    exception <- c(A = 0, v_S = 0, E_T = 0, tau_S = 0, k = 0)
  } else {
    stop("ROUp parameterization must be 'rate' or 'area'.")
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
    c_name = paste0("ROUp", .ROUp_PAR_INFIX[[par]], .ROUp_SUFFIX[[kind]]),
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
    rfun = function(data = NULL, pars) .rfun_ROUp(data$lR, pars, ok = attr(pars, "ok"),
                                                 kind = kind, par = par),
    # Density function (PDF) for single accumulator
    dfun = function(rt, pars) dROUp(rt, pars, kind = kind, par = par),
    # Probability function (CDF) for single accumulator
    pfun = function(rt, pars) pROUp(rt, pars, kind = kind, par = par),
    # Race likelihood combining pfun and dfun
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      log_likelihood_race_missing(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}
