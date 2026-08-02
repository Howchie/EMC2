
dLNR <- function(rt,pars){
  rt <- rt - pars[,"t0"]
  out <- numeric(length(rt))
  ok <- rt > 0
  ok[is.na(ok) | is.infinite(rt)] <- FALSE
  out[ok] <- stats::dlnorm(rt[ok],meanlog=pars[ok,"m"],sdlog=pars[ok,"s"])
  out
}

pLNR <- function(rt,pars){
  rt <- rt - pars[,"t0"]
  out <- numeric(length(rt))
  ok <- rt > 0
  ok[is.na(ok) | is.infinite(rt)] <- FALSE
  out[ok] <- stats::plnorm(rt[ok],meanlog=pars[ok,"m"],sdlog=pars[ok,"s"])
  out

}

rLNR <- function(lR,pars,p_types=c("m","s","t0"),ok=rep(TRUE,dim(pars)[1])){
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  nr <- length(levels(lR))
  dt <- matrix(Inf,ncol=nrow(pars)/nr,nrow=nr)
  t0 <- pars[,"t0"]
  pars <- pars[ok,]
  dt[ok] <- stats::rlnorm(dim(pars)[1],meanlog=pars[,"m"],sdlog=pars[,"s"])
  R <- max.col(-t(dt), ties.method='first')
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
  rt <- matrix(t0,nrow=nr)[pick] + dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  out <- cbind.data.frame(R=R,rt=rt)
  .apply_timed_guess_winner(out, levels(lR))
}


#' The Log-Normal Race Model
#'
#' Model file to estimate the Log-Normal Race Model (LNR) in EMC2.
#'
#' Model files are almost exclusively used in `design()`.
#'
#'
#' @details
#'
#' Default values are used for all parameters that are not explicitly listed in the `formula`
#' argument of `design()`.They can also be accessed with `LNR()$p_types`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default**   | **Mapping**                    | **Interpretation**            |
#'  |-----------|-----------|---------------|-----------|----------------------------|---------------------------|
#'  | *m*       | identity  | \[-Inf, Inf\]   | 1         |                            | Meanlog of the lognormal decision-time distribution |
#'  | *s*       | log       | \[0, Inf\]      | log(1)    |                            | SDlog of the lognormal decision-time distribution |
#'  | *t0*      | log       | \[0, Inf\]      | log(0)    |                            | Additive non-decision-time shift |
#'  | *pContaminant* | probit | \[0, 1\] | qnorm(0) | | Optional contamination probability handled by the data pipeline |
#'
#' Conditional on an accumulator's parameters, its decision time is
#' `T = t0 + Y`, where `log(Y) ~ Normal(m, s^2)`. Thus `m` and `s` are the
#' location and scale on the log-time axis, not the mean and SD on the raw
#' response-time axis. The raw decision-time mean is
#' `t0 + exp(m + s^2 / 2)`, before race selection. The distribution has support
#' above `t0`; the model's default lower bound for `t0` is 0.05 when parameters
#' are checked for simulation.
#'
#' Because the LNR is a race model, it has one independent lognormal
#' accumulator per response option. EMC2 automatically constructs a factor
#' representing the accumulators `lR` (i.e., the latent response) with level
#' names taken from the `R` column in the data. The observed density is the
#' winning accumulator density multiplied by the survivor probabilities of the
#' other accumulators.
#' In `design()`, `matchfun` can be used to automatically create a latent match
#' (`lM`) factor with levels `FALSE` (i.e., the stimulus does not match the accumulator)
#' and `TRUE` (i.e., the stimulus does match the accumulator). This is added internally
#' and can also be used in the model formula, typically for parameters related to
#' the rate of accumulation (see the example below).
#'
#' All model parameters are trial-dependent after the design formulas are
#' evaluated. `pContaminant` is a generic nuisance parameter and is not part of
#' the lognormal race distribution itself.
#'
#' Rouder, J. N., Province, J. M., Morey, R. D., Gomez, P., & Heathcote, A. (2015).
#' The lognormal race: A cognitive-process model of choice and latency with
#' desirable psychometric properties. *Psychometrika, 80*, 491-513.
#' https://doi.org/10.1007/s11336-013-9396-3
#'
#' @return A model list with all the necessary functions for EMC2 to sample
#' @examples
#' # When working with lM it is useful to design  an "average and difference"
#' # contrast matrix, which for binary responses has a simple canonical from:
#' ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
#' # We also define a match function for lM
#' matchfun=function(d)d$S==d$lR
#' # We now construct our design, with v ~ lM and the contrast for lM the ADmat.
#' design_LNRmE <- design(data = forstmann,model=LNR,matchfun=matchfun,
#'                        formula=list(m~lM + E,s~1,t0~1),
#'                        contrasts=list(m=list(lM=ADmat)))
#' # For all parameters that are not defined in the formula, default values are assumed
#' # (see Table above).
#' @export
#'


LNR <- function() {
  list(
    type="RACE",
    c_name = "LNR",
    p_types=c("m" = 1,"s" = log(1),"t0" = log(0), "pContaminant"=qnorm(0)),
    p_types_canonical = c("m", "s", "t0"),
    transform=list(func=c(m = "identity",s = "exp", t0 = "exp", pContaminant="pnorm")),
    bound=list(minmax=cbind(m=c(-Inf,Inf),s = c(0, Inf), t0=c(0.05,Inf)),pContaminant=c(0.001,0.999)),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) pars,
    # Random function for racing accumulators
    rfun=function(data=NULL,pars) rLNR(data$lR, pars, ok = attr(pars, "ok")),
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dLNR(rt,pars),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pLNR(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
      log_likelihood_race_missing(pars=pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}


# Response criteria in a counter model are counts.  `integer_K = TRUE` snaps the
# sampled criterion to the nearest integer (floored at 1) everywhere it is used
# -- density, distribution function, and simulation -- so the R and C++ paths
# stay in lockstep; the C++ kernels apply the same rule from
# ContextForRaceModels::pcounter_integer_K.
#
# Nearest, not ceiling: the snap makes the likelihood flat across the whole
# interval that maps to one count, so the sampled K is only identified up to
# that interval.  Rounding centres each interval on the integer it represents,
# which keeps the raw K unbiased for the criterion; ceiling would put the
# criterion at the interval's upper edge and bias every posterior downwards.
# floor(K + 0.5) rather than round() so that the R and C++ paths agree exactly
# on half-way values (R rounds half to even, C++ std::round does not).
.pcounter_K <- function(K, integer_K) {
  if (!integer_K) return(K)
  out <- floor(K + 0.5)
  out[!is.na(out) & out < 1] <- 1
  out
}

#' The Poisson Counter Race Model
#'
#' The counter model of Pike (1966), LaBerge (1962) and Townsend and Ashby
#' (1983), as specified in the Appendix of Ratcliff and Smith (2004). Each
#' response has an integer-valued evidence counter that accrues unit counts in
#' continuous time as a Poisson process with rate `alpha`, independently and in
#' parallel with the other counters, until it reaches its response criterion
#' `K`. The response is made by whichever counter reaches criterion first.
#'
#' Because the times between counts are exponential, the time for accumulator
#' `i` to collect `K_i` counts is Erlang (Gamma with integer shape), so the
#' joint density of a response by accumulator `i` at time `t` is
#'
#' \deqn{g_i(t) = \frac{(\alpha_i t)^{K_i - 1} \alpha_i e^{-\alpha_i t}}{(K_i - 1)!}
#'   \prod_{j \neq i} \sum_{n=0}^{K_j - 1} \frac{(\alpha_j t)^n}{n!} e^{-\alpha_j t},}
#'
#' which is Equations A10a and A10b of Ratcliff and Smith (2004) generalised
#' from two counters to the `R` accumulators of an EMC2 race. The first factor
#' is the Erlang density of the winner and the sum is the Erlang survivor
#' function of each loser, so this is the standard EMC2 race likelihood applied
#' to Erlang finish times.
#'
#' The model parameters are:
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Interpretation** |
#' |---|---|---|---|---|
#' | *alpha* | log | \[0, Inf\] | log(10) | Poisson accrual rate of the counter (counts per second). |
#' | *K* | log | \[0, Inf\] | log(5) | Response criterion; the number of counts needed to respond. |
#' | *t0* | log | \[0, Inf\] | log(0) | Non-decision time, a nonnegative shift of the finish-time distribution. |
#' | *pContaminant* | probit | \[0, 1\] | qnorm(0) | Optional contamination probability handled by the data pipeline. |
#'
#' The mean decision time of a counter is `K / alpha`, so `alpha` plays the role
#' of a drift rate and `K` the role of a threshold. Ratcliff and Smith held the
#' summed rate of the two counters constant and let the relative rate
#' `pi = alpha_a / (alpha_a + alpha_b)` carry the stimulus effect; that
#' constraint couples accumulators within a trial, so it is not imposed here.
#' Per-accumulator `alpha` is the unconstrained equivalent, and the sum can be
#' held constant through the design if wanted.
#'
#' @param integer_K Logical. If `FALSE` (the default) `K` is continuous, giving
#'   a fractional-counter generalisation with a smooth likelihood. If `TRUE`,
#'   `K` is snapped to the nearest integer (floored at 1) before entering the
#'   likelihood, so the model is the literal counter model.
#'
#'   Snapping makes the likelihood a step function of `K`, flat across the
#'   interval of `K` that maps to one count. The sampled `K` is therefore
#'   identified only up to that interval and its posterior will be roughly as
#'   wide as the interval, centred on the criterion; it is the snapped value,
#'   supplied as the derived parameter `K_int`, that is the estimate of the
#'   criterion. Pass `map = TRUE, add_recalculated = TRUE, remove_constants =
#'   FALSE` to [recovery()], [get_pars()] and friends to summarise `K_int`
#'   rather than the raw `K`. `remove_constants = FALSE` is needed because a
#'   criterion that is perfectly determined is constant across draws, and the
#'   default would drop it for precisely the fits that recovered it best.
#'
#' @return A model list compatible with `design()`.
#' @references
#' Ratcliff, R., & Smith, P. L. (2004). A comparison of sequential sampling
#' models for two-choice reaction time. *Psychological Review, 111*(2),
#' 333-367.
#'
#' Townsend, J. T., & Ashby, F. G. (1983). *The stochastic modeling of
#' elementary psychological processes*. Cambridge University Press.
#' @export
PCOUNTER <- function(integer_K = FALSE) {
  integer_K <- isTRUE(integer_K)

  dPCOUNTER <- function(rt, pars) {
    out <- numeric(length(rt))
    tt <- rt - pars[, "t0"]
    ok <- tt > 0 & is.finite(tt) & pars[, "alpha"] > 0 & pars[, "K"] > 0
    ok[is.na(ok)] <- FALSE
    out[ok] <- stats::dgamma(tt[ok], shape = .pcounter_K(pars[ok, "K"], integer_K),
                             rate = pars[ok, "alpha"])
    out
  }

  pPCOUNTER <- function(rt, pars) {
    out <- numeric(length(rt))
    tt <- rt - pars[, "t0"]
    ok <- tt > 0 & is.finite(tt) & pars[, "alpha"] > 0 & pars[, "K"] > 0
    ok[is.na(ok)] <- FALSE
    out[ok] <- stats::pgamma(tt[ok], shape = .pcounter_K(pars[ok, "K"], integer_K),
                             rate = pars[ok, "alpha"])
    out
  }

  rPCOUNTER <- function(lR, pars, p_types = c("alpha", "K", "t0"), ok = rep(TRUE, dim(pars)[1])) {
    if (!all(p_types %in% dimnames(pars)[[2]])) {
      stop("pars must have columns ", paste(p_types, collapse = " "))
    }
    nr <- length(levels(lR))
    dt <- matrix(Inf, ncol = nrow(pars) / nr, nrow = nr)
    idx_ok <- which(ok)
    pars <- pars[idx_ok, , drop = FALSE]
    alpha <- pars[, "alpha"]
    K <- pars[, "K"]
    t0 <- pars[, "t0"]
    ok_draw <- is.finite(alpha) & is.finite(K) & is.finite(t0) & alpha > 0 & K > 0
    if (any(ok_draw)) {
      dt[idx_ok[ok_draw]] <- stats::rgamma(sum(ok_draw),
                                           shape = .pcounter_K(K[ok_draw], integer_K),
                                           rate = alpha[ok_draw]) + t0[ok_draw]
    }
    R <- max.col(-t(dt), ties.method = "first")
    pick <- cbind(R, seq_len(ncol(dt)))
    rt <- dt[pick]
    R <- factor(levels(lR)[R], levels = levels(lR))
    out <- cbind.data.frame(R = R, rt = rt)
    .apply_timed_guess_winner(out, levels(lR))
  }

  list(
    type = "RACE",
    c_name = if (integer_K) "PCOUNTER_INTK" else "PCOUNTER",
    p_types = c("alpha" = log(10), "K" = log(5), "t0" = log(0), "pContaminant" = qnorm(0)),
    p_types_canonical = c("alpha", "K", "t0"),
    transform = list(func = c(alpha = "exp", K = "exp", t0 = "exp", pContaminant = "pnorm")),
    bound = list(
      minmax = cbind(alpha = c(1e-6, Inf), K = c(1e-6, Inf), t0 = c(0, Inf), pContaminant = c(0.001, 0.999)),
      exception = c(t0 = 0, pContaminant = 0)
    ),
    # In integer mode the sampled K diffuses freely inside the flat interval
    # that ceiling() maps to one criterion, so the raw K is not the criterion.
    # Expose the snapped count as a derived parameter for mapped_pars() and
    # posterior summaries.  Extra columns are appended after the p_types block,
    # so the alpha/K/t0 prefix the C++ kernels index positionally is unchanged.
    Ttransform = function(pars, dadm) {
      if (!integer_K) return(pars)
      cbind(pars, K_int = .pcounter_K(pars[, "K"], TRUE))
    },
    rfun = function(data = NULL, pars) rPCOUNTER(data$lR, pars, ok = attr(pars, "ok")),
    dfun = function(rt, pars) dPCOUNTER(rt, pars),
    pfun = function(rt, pars) pPCOUNTER(rt, pars),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      log_likelihood_race_missing(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}

# -----------------------------------------------------------------------------
# Refactored PCOUNTER implementation (see Math/new.md).
# The old implementation above is retained only as source-level compatibility
# history; this definition is the one returned when the package is loaded.
.pcounter_k_new <- function(k) {
  out <- floor(k + 0.5)
  out[is.finite(out) & out < 1] <- 1
  out
}

.pcounter_logsumexp_new <- function(x) {
  if (!length(x) || all(x == -Inf)) return(-Inf)
  m <- max(x)
  m + log(sum(exp(x - m)))
}

.pcounter_logdiffexp_new <- function(a, b) {
  if (!is.finite(a)) return(-Inf)
  if (!is.finite(b)) return(a)
  if (b >= a) return(-Inf)
  a + log1p(-exp(b - a))
}

.pcounter_stirling_logs_new <- function(nmax) {
  out <- matrix(-Inf, nrow = nmax + 1L, ncol = nmax + 1L)
  out[1L, 1L] <- 0
  if (nmax < 1L) return(out)
  for (n in seq_len(nmax)) {
    for (j in seq_len(n)) {
      a <- out[n, j]
      b <- if (j < n) log(n - 1) + out[n, j + 1L] else -Inf
      out[n + 1L, j + 1L] <- .pcounter_logsumexp_new(c(a, b))
    }
  }
  out
}

.pcounter_log_rising_new <- function(a, n) {
  if (n == 0L) return(0)
  sum(log(a + seq_len(n) - 1))
}

.pcounter_log_lm_new <- function(a, nu, sv) {
  if (sv <= 0) {
    logL <- -nu * a
    return(c(logL = logL, logM = log(nu) + logL))
  }
  shape <- nu^2 / sv^2
  rate <- nu / sv^2
  logL <- -shape * log1p(a / rate)
  c(logL = logL, logM = log(shape) - log(rate + a) + logL)
}

.pcounter_log_h_new <- function(n, t, nu, sv, gamma, stirling) {
  lm <- .pcounter_log_lm_new(t, nu, sv)
  if (sv <= 0) return(lm[["logL"]] + .pcounter_log_rising_new(nu / gamma, n))
  shape <- nu^2 / sv^2
  rate <- nu / sv^2
  terms <- vapply(0:n, function(j) {
    coeff <- stirling[n + 1L, j + 1L]
    if (!is.finite(coeff)) return(-Inf)
    coeff + .pcounter_log_rising_new(shape, j) - j * log(gamma) - j * log(rate + t)
  }, numeric(1))
  lm[["logL"]] + .pcounter_logsumexp_new(terms)
}

.pcounter_log_pi_phi_new <- function(n, t, nu, sv, gamma, stirling, gamma_zero) {
  if (gamma_zero) {
    if (sv <= 0) {
      logpi <- -nu * t + if (n == 0L) 0 else n * log(nu * t) - lgamma(n + 1)
      return(c(logpi = logpi, logphi = log(nu) + logpi))
    }
    shape <- nu^2 / sv^2
    rate <- nu / sv^2
    lm <- .pcounter_log_lm_new(t, nu, sv)
    logt <- if (t > 0) log(t) else -Inf
    logpi <- lm[["logL"]] + .pcounter_log_rising_new(shape, n) - lgamma(n + 1) +
      n * logt - n * log(rate + t)
    return(c(logpi = logpi, logphi = logpi + log(shape + n) - log(rate + t)))
  }
  q <- -expm1(-gamma * t)
  logq <- if (q > 0) log(q) else -Inf
  c(logpi = n * logq - lgamma(n + 1) + .pcounter_log_h_new(n, t, nu, sv, gamma, stirling),
    logphi = log(gamma) + n * logq - lgamma(n + 1) +
      .pcounter_log_h_new(n + 1L, t, nu, sv, gamma, stirling))
}

.pcounter_log_tail_new <- function(start, t, nu, sv, gamma, logr, stirling, gamma_zero, phi = FALSE) {
  ns <- start:(start + 64L)
  vals <- vapply(ns, function(n) {
    z <- .pcounter_log_pi_phi_new(n, t, nu, sv, gamma, stirling, gamma_zero)
    if (phi) z[["logphi"]] else z[["logpi"]]
  }, numeric(1))
  .pcounter_logsumexp_new(vals + ns * logr)
}

.pcounter_log_fixed_cdf_new <- function(start, t, nu, sv, gamma, stirling, gamma_zero) {
  ns <- start:(start + 64L)
  vals <- vapply(ns, function(n)
    .pcounter_log_pi_phi_new(n, t, nu, sv, gamma, stirling, gamma_zero)[["logpi"]], numeric(1))
  .pcounter_logsumexp_new(vals)
}

.pcounter_log_geom_cdf_new <- function(start, t, nu, sv, gamma, omega, stirling, gamma_zero) {
  lr <- log(omega) - log1p(omega)
  ns <- start:(start + 64L)
  vals <- vapply(ns, function(n) {
    lp <- .pcounter_log_pi_phi_new(n, t, nu, sv, gamma, stirling, gamma_zero)[["logpi"]]
    exponent <- n - start
    if (exponent == 0L) -Inf else lp + log(-expm1(exponent * lr))
  }, numeric(1))
  .pcounter_logsumexp_new(vals)
}

.pcounter_eval_one_new <- function(t, nu, sv, gamma, k, omega, stirling, eps = 1e-8) {
  t <- unname(t); nu <- unname(nu); sv <- unname(sv); gamma <- unname(gamma)
  k <- unname(k); omega <- unname(omega)
  if (is.nan(t) || !is.finite(nu) || !is.finite(sv) || !is.finite(gamma) ||
      !is.finite(k) || !is.finite(omega) || nu <= 0 || sv < 0 || gamma < 0 ||
      k <= 0 || omega < 0) return(c(logf = -Inf, logS = 0, logF = -Inf))
  if (is.infinite(t) && t > 0) return(c(logf = -Inf, logS = -Inf, logF = 0))
  if (!is.finite(t) || t <= 0) return(c(logf = -Inf, logS = 0, logF = -Inf))
  kk <- .pcounter_k_new(k)
  gamma_zero <- gamma < eps
  sv_zero <- sv < eps
  omega_zero <- omega < eps
  if (omega_zero) {
    logs <- vapply(0:(kk - 1L), function(n)
      .pcounter_log_pi_phi_new(n, t, nu, if (sv_zero) 0 else sv,
                               if (gamma_zero) 0 else gamma, stirling, gamma_zero)[["logpi"]], numeric(1))
    logS <- min(0, .pcounter_logsumexp_new(logs))
    z <- .pcounter_log_pi_phi_new(kk - 1L, t, nu, if (sv_zero) 0 else sv,
                                  if (gamma_zero) 0 else gamma, stirling, gamma_zero)
    logF <- if (logS > -1e-7) .pcounter_log_fixed_cdf_new(kk, t, nu,
      if (sv_zero) 0 else sv, if (gamma_zero) 0 else gamma, stirling, gamma_zero)
      else if (logS < 0) log(-expm1(logS)) else -Inf
    return(c(logf = z[["logphi"]], logS = logS, logF = logF))
  }
  lp <- -log1p(omega)
  lr <- log(omega) + lp
  p <- exp(lp)
  if (gamma_zero) {
    ar <- p * t
    logd <- 0
  } else {
    q <- -expm1(-gamma * t)
    logd <- log1p(-exp(lr) * q)
    ar <- t + logd / gamma
  }
  lm_ar <- .pcounter_log_lm_new(ar, nu, if (sv_zero) 0 else sv)
  low <- if (kk >= 2L) 0:(kk - 2L) else integer()
  vals <- if (length(low)) vapply(low, function(n)
    .pcounter_log_pi_phi_new(n, t, nu, if (sv_zero) 0 else sv,
                             if (gamma_zero) 0 else gamma, stirling, gamma_zero), numeric(2)) else
    matrix(numeric(), nrow = 2L, ncol = 0L)
  loglow <- if (length(low)) .pcounter_logsumexp_new(vals[1L, ]) else -Inf
  loglowr <- if (length(low)) .pcounter_logsumexp_new(vals[1L, ] + low * lr) else -Inf
  logtail <- .pcounter_logdiffexp_new(lm_ar[["logL"]], loglowr)
  if (length(low) && (loglowr > lm_ar[["logL"]] - 1e-7 || !is.finite(logtail)))
    logtail <- .pcounter_log_tail_new(kk - 1L, t, nu, if (sv_zero) 0 else sv,
                                       if (gamma_zero) 0 else gamma, lr, stirling,
                                       gamma_zero, phi = FALSE)
  logS <- min(0, .pcounter_logsumexp_new(c(loglow, (1 - kk) * lr + logtail)))
  logA <- lm_ar[["logM"]] - if (gamma_zero) 0 else logd
  loglowf <- if (length(low)) .pcounter_logsumexp_new(vals[2L, ] + low * lr) else -Inf
  logf <- lp + (1 - kk) * lr + .pcounter_logdiffexp_new(logA, loglowf)
  if (length(low) && (loglowf > logA - 1e-7 || !is.finite(logf))) {
    tailf <- .pcounter_log_tail_new(kk - 1L, t, nu, if (sv_zero) 0 else sv,
                                    if (gamma_zero) 0 else gamma, lr, stirling,
                                    gamma_zero, phi = TRUE)
    logf <- lp + (1 - kk) * lr + tailf
  }
  logF <- if (kk == 1L) log(-expm1(lm_ar[["logL"]])) else if (logS > -1e-7)
    .pcounter_log_geom_cdf_new(kk - 1L, t, nu, if (sv_zero) 0 else sv,
                               if (gamma_zero) 0 else gamma, omega, stirling, gamma_zero)
    else if (logS < 0) log(-expm1(logS)) else -Inf
  c(logf = logf, logS = logS, logF = logF)
}

#' The gamma-mixed, self-exciting Poisson Counter race model.
#'
#' `integer_K` is accepted for source compatibility with the former model and
#' ignored.  The new criterion `k` is integer-valued by definition.
#' @param integer_K Deprecated compatibility argument; ignored.
#' @return A model list compatible with `design()`.
#' @export
PCOUNTER <- function(integer_K = NULL) {
  req <- c("nu", "sv", "gamma", "k", "omega", "t0")
  prepare <- function(pars) {
    if (!all(req %in% colnames(pars))) stop("pars must have columns ", paste(req, collapse = " "))
    ks <- .pcounter_k_new(pars[, "k"])
    ks <- ks[is.finite(ks)]
    mk <- if (length(ks)) max(ks) else 1
    .pcounter_stirling_logs_new(max(65L, as.integer(mk) + 65L))
  }
  dPCOUNTER <- function(rt, pars) {
    st <- prepare(pars); out <- numeric(length(rt))
    for (i in seq_along(rt)) out[i] <- exp(.pcounter_eval_one_new(
      rt[i] - pars[i, "t0"], pars[i, "nu"], pars[i, "sv"], pars[i, "gamma"],
      pars[i, "k"], pars[i, "omega"], st)[["logf"]])
    out
  }
  pPCOUNTER <- function(rt, pars) {
    st <- prepare(pars); out <- numeric(length(rt))
    for (i in seq_along(rt)) {
      z <- .pcounter_eval_one_new(rt[i] - pars[i, "t0"], pars[i, "nu"], pars[i, "sv"],
                                  pars[i, "gamma"], pars[i, "k"], pars[i, "omega"], st)
      out[i] <- if (z[["logF"]] == 0) 1 else if (is.finite(z[["logF"]])) exp(z[["logF"]]) else 0
    }
    out
  }
  rPCOUNTER <- function(lR, pars, p_types = req, ok = rep(TRUE, nrow(pars))) {
    if (!all(p_types %in% colnames(pars))) stop("pars must have columns ", paste(p_types, collapse = " "))
    if (is.null(ok)) ok <- rep(TRUE, nrow(pars))
    nr <- length(levels(lR)); nd <- nrow(pars) / nr
    if (nr < 1L || nrow(pars) %% nr != 0L) stop("pars rows must be a multiple of accumulators")
    dt <- matrix(Inf, nrow = nr, ncol = nd)
    for (ii in which(ok)) {
      nu <- pars[ii, "nu"]; sv <- pars[ii, "sv"]; ga <- pars[ii, "gamma"]
      kk <- .pcounter_k_new(pars[ii, "k"]); om <- pars[ii, "omega"]
      if (!is.finite(nu) || !is.finite(sv) || !is.finite(ga) || !is.finite(kk) ||
          !is.finite(om) || !is.finite(pars[ii, "t0"]) || nu <= 0 || sv < 0 ||
          ga < 0 || kk < 1 || om < 0 || pars[ii, "t0"] < 0) next
      vv <- if (sv < 1e-8) nu else stats::rgamma(1L, nu^2 / sv^2, rate = nu / sv^2)
      if (!is.finite(vv) || vv <= 0) next
      excess <- if (om < 1e-8) 0L else stats::rgeom(1L, 1 / (1 + om))
      m <- kk + excess
      dt[ii] <- sum(stats::rexp(m, rate = vv + ga * seq.int(0, m - 1L))) + pars[ii, "t0"]
    }
    R <- max.col(-t(dt), ties.method = "first"); pick <- cbind(R, seq_len(nd))
    out <- cbind.data.frame(R = factor(levels(lR)[R], levels = levels(lR)), rt = dt[pick])
    .apply_timed_guess_winner(out, levels(lR))
  }
  list(
    type = "RACE", c_name = "PCOUNTER",
    p_types = c(nu = log(10), sv = log(0), gamma = log(0), k = log(5), omega = log(0),
                t0 = log(0), pContaminant = qnorm(0)),
    p_types_canonical = c("nu", "sv", "gamma", "k", "omega", "t0"),
    transform = list(func = c(nu = "exp", sv = "exp", gamma = "exp", k = "exp",
                              omega = "exp", t0 = "exp", pContaminant = "pnorm")),
    # Keep estimated variability/self-excitation parameters above the numerical
    # branch threshold. Exact zero remains available only as an explicit fixed
    # constant through the exceptions below.
    bound = list(minmax = cbind(nu = c(1e-6, Inf), sv = c(1e-4, Inf), gamma = c(1e-4, Inf),
                                k = c(1, Inf), omega = c(1e-4, Inf), t0 = c(0, Inf),
                                pContaminant = c(0.001, 0.999)),
                 exception = c(sv = 0, gamma = 0, omega = 0, t0 = 0, pContaminant = 0)),
    Ttransform = function(pars, dadm) pars,
    rfun = function(data = NULL, pars) rPCOUNTER(data$lR, pars, ok = attr(pars, "ok")),
    dfun = dPCOUNTER, pfun = pPCOUNTER,
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10))
      log_likelihood_race_missing(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
  )
}

#' The Ex-Gaussian Race Model
#'
#' Race model in which each accumulator has an ex-Gaussian finish-time
#' distribution. For accumulator `i`,
#' `T_i = mu_i + G_i + E_i`, with
#' `G_i ~ Normal(0, sigma_i^2)` and
#' `E_i ~ Exponential(rate = 1 / tau_i)`. Equivalently, `mu` is the Gaussian
#' location, `sigma` is the Gaussian SD, and `tau` is the mean of the
#' exponential tail. The observed response is the accumulator with the
#' smallest finish time.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Interpretation** |
#' |---|---|---|---|---|
#' | *mu* | log | \[0, Inf\] | log(.4) | Location of the Gaussian component. The current bounds restrict this location to be nonnegative. |
#' | *sigma* | log | \[0, Inf\] | log(.05) | SD of the Gaussian component. |
#' | *tau* | log | \[0, Inf\] | log(.1) | Mean of the exponential component; its rate is `1 / tau`. |
#' | *pContaminant* | probit | \[0, 1\] | qnorm(0) | Optional contamination probability handled by the data pipeline. |
#'
#' The ex-Gaussian distribution has support on the real line. Its mean is
#' `mu + tau` and its variance is `sigma^2 + tau^2`; these are descriptive
#' distributional summaries rather than separate model parameters. Unlike the
#' stop-signal ex-Gaussian models [SSEXG()] and [SSRDEX()], `REXG()` does not
#' apply a lower truncation bound. The `pContaminant` parameter is generic
#' nuisance infrastructure and is not part of the ex-Gaussian distribution.
#' EMC2 creates one accumulator per response level in `R` and evaluates the
#' race likelihood from the accumulator density and survivor functions.
#'
#' @return A model list compatible with `design()`.
#' @export
REXG <- function() {
  dREXG <- function(rt, pars) {
    out <- numeric(length(rt))
    ok <- is.finite(rt) & pars[, "sigma"] > 0 & pars[, "tau"] > 0
    ok[is.na(ok)] <- FALSE
    out[ok] <- dexGaussian(rt[ok], pars[ok, c("mu", "sigma", "tau"), drop = FALSE])
    out
  }

  pREXG <- function(rt, pars) {
    out <- numeric(length(rt))
    ok <- is.finite(rt) & pars[, "sigma"] > 0 & pars[, "tau"] > 0
    ok[is.na(ok)] <- FALSE
    out[ok] <- pexGaussian(rt[ok], pars[ok, c("mu", "sigma", "tau"), drop = FALSE])
    out[rt == Inf] <- 1
    out
  }

  rREXG <- function(lR, pars, p_types = c("mu", "sigma", "tau"),
                    ok = rep(TRUE, dim(pars)[1])) {
    if (!all(p_types %in% dimnames(pars)[[2]])) {
      stop("pars must have columns ", paste(p_types, collapse = " "))
    }
    nr <- length(levels(lR))
    dt <- matrix(Inf, ncol = nrow(pars) / nr, nrow = nr)
    idx_ok <- which(ok)
    pars <- pars[idx_ok, , drop = FALSE]
    mu <- pars[, "mu"]
    sigma <- pars[, "sigma"]
    tau <- pars[, "tau"]
    ok_draw <- is.finite(mu) & is.finite(sigma) & is.finite(tau) & sigma > 0 & tau > 0
    if (any(ok_draw)) {
      dt[idx_ok[ok_draw]] <-
        stats::rnorm(sum(ok_draw), mean = mu[ok_draw], sd = sigma[ok_draw]) +
        stats::rexp(sum(ok_draw), rate = 1 / tau[ok_draw])
    }
    R <- max.col(-t(dt), ties.method = "first")
    pick <- cbind(R, seq_len(ncol(dt)))
    rt <- dt[pick]
    R <- factor(levels(lR)[R], levels = levels(lR))
    out <- cbind.data.frame(R = R, rt = rt)
    .apply_timed_guess_winner(out, levels(lR))
  }

  list(
    type = "RACE",
    c_name = "REXG",
    p_types = c(mu = log(.4), sigma = log(.05), tau = log(.1), pContaminant = qnorm(0)),
    p_types_canonical = c("mu", "sigma", "tau"),
    transform = list(func = c(mu = "exp", sigma = "exp", tau = "exp",
                              pContaminant = "pnorm")),
    bound = list(
      minmax = cbind(mu = c(0, Inf), sigma = c(1e-6, Inf),
                     tau = c(1e-6, Inf), pContaminant = c(0.001, 0.999)),
      exception = c(pContaminant = 0)
    ),
    Ttransform = function(pars, dadm) pars,
    rfun = function(data = NULL, pars) rREXG(data$lR, pars, ok = attr(pars, "ok")),
    dfun = function(rt, pars) dREXG(rt, pars),
    pfun = function(rt, pars) pREXG(rt, pars),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      log_likelihood_race_missing(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}
