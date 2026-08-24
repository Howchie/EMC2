
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
  pick <- cbind(R,1:dim(dt)[2])
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
#'  | *pContaminant* | probit | \[0, 1\] | qnorm(0) | | Optional *omission* contaminant probability: mass at `rt = Inf` only, handled by the data pipeline |
#'  | *pGuess* | probit | \[0, 1\] | qnorm(0) | | Optional uniform *guess* (outlier) probability, mixed into observed RT densities over the guess window |
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
#' evaluated. `pContaminant` and `pGuess` are generic nuisance parameters and
#' are not part of the lognormal race distribution itself: `pContaminant` is a
#' Bernoulli *omission* rate that puts mass only at `rt = Inf`, while `pGuess`
#' mixes a uniform outlier density into observed RTs. Both are proportions among
#' *retained* trials -- they are applied after truncation renormalisation.
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
    p_types=c("m" = 1,"s" = log(1),"t0" = log(0), "pContaminant"=qnorm(0), "pGuess"=qnorm(0)),
    p_types_canonical = c("m", "s", "t0"),
    transform=list(func=c(m = "identity",s = "exp", t0 = "exp", pContaminant="pnorm", pGuess="pnorm")),
    # NOTE: pContaminant's bound sits outside cbind() (a long-standing slip, so
    # it is not actually a minmax column).  Left as-is deliberately: moving it
    # would change pContaminant's established behaviour.
    # pGuess needs the exception at 0: without it the [0.001, 0.999] minmax
    # clamps the default (pnorm(-Inf) = 0) up to 0.001, making the nuisance
    # parameter silently active.
    bound=list(minmax=cbind(m=c(-Inf,Inf),s = c(0, Inf), t0=c(0.05,Inf),pGuess=c(0.001,0.999)),
               exception=c(pGuess=0),pContaminant=c(0.001,0.999)),
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


# `k` is sampled as a non-negative offset.  The canonical count criterion is
# two events plus the nearest non-negative integer offset.  Keeping this
# conversion in one helper ensures that the R density, C++ likelihood, and
# simulator all use the same threshold convention.
.pcounter_k <- function(k) {
  out <- floor(k + 0.5)
  out[is.finite(out) & out < 0] <- 0
  out + 2
}

# The C++ implementation uses an int index for the exact finite-K recurrence.
# This is a machine-representability check, not a substantive model bound.
# Tail evaluation needs 67 additional row indices.
.pcounter_k_supported <- function(k) {
  is.finite(k) & k >= 0 &
    floor(k + 0.5) <= (.Machine$integer.max - 67)
}


#' The gamma-mixed, self-exciting Poisson Counter race model.
#'
#' Each accumulator is a pure-birth counter with event rate
#' `nu + gamma * n`. The trialwise input rate has mean `nu` and standard
#' deviation `sv`; the response threshold is
#' `K = 2 + floor(k + 1/2)` and has a geometric excess with mean `omega`.
#' The exact-zero values of `sv`, `gamma`, and `omega` select simpler nested
#' branches. The closed form evaluates its Stirling recurrence row by row,
#' avoiding a quadratic threshold-sized allocation.
#'
#' The proper `sv = 0`, `omega = 0`, `gamma > 0` member overlaps with the
#' non-defective boundary of FRQ: for integer `alpha = K`, `beta = nu/gamma`,
#' `h = 1`, and `lambda = gamma`, both have
#' `F(t) = I_{1-exp(-gamma*t)}(K, nu/gamma)`. FRQ is not otherwise redundant:
#' it supplies a finite-reservoir defect (`h < 1`) and threshold-shape
#' relaxation, while PCOUNTER supplies trialwise input variation,
#' self-excitation, and geometric threshold variation.
#'
#' @return A model list compatible with `design()`.
#' @export
PCOUNTER <- function() {
  req <- c("nu", "sv", "gamma", "k", "omega", "t0")
  check <- function(rt, pars) {
    if (!all(req %in% colnames(pars)))
      stop("pars must have columns ", paste(req, collapse = " "))
    if (nrow(pars) != length(rt))
      stop("PCOUNTER requires one parameter row per response time")
  }
  dPCOUNTER <- function(rt, pars) {
    check(rt, pars)
    dpcounter(rt, pars[, "nu"], pars[, "sv"], pars[, "gamma"], pars[, "k"],
              pars[, "omega"], pars[, "t0"])
  }
  pPCOUNTER <- function(rt, pars) {
    check(rt, pars)
    ppcounter(rt, pars[, "nu"], pars[, "sv"], pars[, "gamma"], pars[, "k"],
              pars[, "omega"], pars[, "t0"])
  }
  rPCOUNTER <- function(lR, pars, p_types = req, ok = rep(TRUE, nrow(pars))) {
    if (!all(p_types %in% colnames(pars))) stop("pars must have columns ", paste(p_types, collapse = " "))
    if (is.null(ok)) ok <- rep(TRUE, nrow(pars))
    if (length(ok) != nrow(pars)) stop("ok must have one value per parameter row")
    nr <- length(levels(lR))
    if (nr < 1L || nrow(pars) %% nr != 0L) stop("pars rows must be a multiple of accumulators")
    nd <- nrow(pars) / nr
    dt <- matrix(Inf, nrow = nr, ncol = nd)
    for (ii in which(ok)) {
      nu <- pars[ii, "nu"]; sv <- pars[ii, "sv"]; ga <- pars[ii, "gamma"]
      k_raw <- pars[ii, "k"]
      kk <- .pcounter_k(k_raw); om <- pars[ii, "omega"]
      if (!is.finite(nu) || !is.finite(sv) || !is.finite(ga) || !is.finite(k_raw) ||
          !is.finite(kk) || !is.finite(om) || !is.finite(pars[ii, "t0"]) ||
          !.pcounter_k_supported(k_raw) || nu <= 0 || sv < 0 || ga < 0 ||
          k_raw < 0 || om < 0 || pars[ii, "t0"] < 0) next
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
    p_types = c(nu = log(10), sv = log(0), gamma = log(0), k = log(3), omega = log(0),
                t0 = log(0), pContaminant = qnorm(0), pGuess = qnorm(0)),
    p_types_canonical = c("nu", "sv", "gamma", "k", "omega", "t0"),
    transform = list(func = c(nu = "exp", sv = "exp", gamma = "exp", k = "exp",
                              omega = "exp", t0 = "exp", pContaminant = "pnorm",
                              pGuess = "pnorm")),
    # Keep estimated variability/self-excitation parameters above the numerical
    # branch threshold. Exact zero remains available only as an explicit fixed
    # constant through the exceptions below.  `k` is a non-negative offset; the
    # canonical count criterion is K = 2 + floor(k + 0.5), so every accumulator needs
    # at least two events before responding.
    bound = list(minmax = cbind(nu = c(1e-6, Inf), sv = c(1e-4, Inf), gamma = c(1e-4, Inf),
                                # Finite thresholds are valid; the implementation
                                # checks only machine representability of its index.
                                k = c(0, Inf), omega = c(1e-4, Inf), t0 = c(0, Inf),
                                pContaminant = c(0.001, 0.999),
                                pGuess = c(0.001, 0.999)),
                 exception = c(k = 0, sv = 0, gamma = 0, omega = 0, t0 = 0, pContaminant = 0,
                               pGuess = 0)),
    Ttransform = function(pars, dadm) pars,
    rfun = function(data = NULL, pars) rPCOUNTER(data$lR, pars, ok = attr(pars, "ok")),
    dfun = dPCOUNTER, pfun = pPCOUNTER,
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10))
      log_likelihood_race_missing(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
  )
}

#' The Ex-Gaussian Race Model
#'
#' Race model in which each accumulator has a zero-truncated ex-Gaussian
#' process-time distribution. For accumulator `i`, let
#' `X*_i = mu_i + G_i + E_i`, with `G_i ~ Normal(0, sigma_i^2)` and
#' `E_i ~ Exponential(rate = 1 / tau_i)`. The process time is
#' `X_i = X*_i | X*_i > 0`, and the finish time is `T_i = t0_i + X_i`,
#' where `t0_i` is nonnegative. The observed response is the accumulator
#' with the smallest finish time.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Interpretation** |
#' |---|---|---|---|---|
#' | *mu* | log | \[0, Inf\] | log(.4) | Location of the underlying Gaussian component. |
#' | *sigma* | log | \[0, Inf\] | log(.05) | SD of the underlying Gaussian component. |
#' | *tau* | log | \[0, Inf\] | log(.1) | Mean of the underlying exponential component; its rate is `1 / tau`. |
#' | *t0* | log | \[0, Inf\] | log(0) | Additive non-decision-time shift. |
#' | *pContaminant* | probit | \[0, 1\] | qnorm(0) | Optional *omission* contaminant probability: mass at `rt = Inf` only, handled by the data pipeline |
#' | *pGuess* | probit | \[0, 1\] | qnorm(0) | Optional uniform *guess* (outlier) probability, mixed into observed RT densities over the guess window |
#'
#' The underlying, untruncated ex-Gaussian component has mean `mu + tau` and
#' variance `sigma^2 + tau^2`; these summaries describe `X*`, not the
#' zero-truncated process time `X` or the shifted finish time `T`.
#' `REXG()` conditions the process time on `X* > 0`, so the finish-time
#' support is `rt > t0`. For `rt <= t0`, the density and CDF are zero and the
#' survivor is one. Above `t0`, the density, CDF, and survivor include the
#' normalizer `1 - F_EXG(0)`. The `pContaminant` parameter is generic nuisance
#' infrastructure and is not part of the ex-Gaussian distribution.
#' EMC2 creates one accumulator per response level in `R` and evaluates the
#' race likelihood from the accumulator density and survivor functions.
#'
#' @return A model list compatible with `design()`.
#' @export
REXG <- function() {
  dREXG <- function(rt, pars) {
    out <- numeric(length(rt))
    t0 <- pars[, "t0"]
    ok <- is.finite(rt) & pars[, "sigma"] > 0 & pars[, "tau"] > 0 &
      is.finite(t0) & t0 >= 0
    ok[is.na(ok)] <- FALSE
    if (any(ok)) {
      pars_ok <- cbind(pars[ok, c("mu", "sigma", "tau"), drop = FALSE], exg_lb = 0)
      out[ok] <- dtexGaussian(rt[ok] - t0[ok], pars_ok)
    }
    out
  }
  pREXG <- function(rt, pars) {
    out <- numeric(length(rt))
    t0 <- pars[, "t0"]
    valid <- pars[, "sigma"] > 0 & pars[, "tau"] > 0 &
      is.finite(t0) & t0 >= 0
    valid[is.na(valid)] <- FALSE
    ok <- is.finite(rt) & valid
    if (any(ok)) {
      pars_ok <- cbind(pars[ok, c("mu", "sigma", "tau"), drop = FALSE], exg_lb = 0)
      out[ok] <- ptexGaussian(rt[ok] - t0[ok], pars_ok)
    }
    out[is.infinite(rt) & rt > 0 & valid] <- 1
    out
  }
  rREXG <- function(lR, pars, p_types = c("mu", "sigma", "tau", "t0"),
                    ok = rep(TRUE, dim(pars)[1])) {
    if (!all(p_types %in% dimnames(pars)[[2]]))
      stop("pars must have columns ", paste(p_types, collapse = " "))
    if (is.null(ok)) ok <- rep(TRUE, nrow(pars))
    if (length(ok) != nrow(pars))
      stop("ok must have one value per row of pars")
    ok <- as.logical(ok)
    ok[is.na(ok)] <- FALSE
    idx_ok <- which(ok)
    nr <- length(levels(lR))
    n_tri <- nrow(pars) / nr
    dt <- matrix(Inf, nrow = nr, ncol = n_tri)
    pars_ok <- pars[idx_ok, , drop = FALSE]
    mu_ok <- pars_ok[, "mu"]; sigma_ok <- pars_ok[, "sigma"]; tau_ok <- pars_ok[, "tau"]
    t0_ok <- pars_ok[, "t0"]
    ok_draw <- is.finite(mu_ok) & is.finite(sigma_ok) & is.finite(tau_ok) &
      is.finite(t0_ok) & sigma_ok > 0 & tau_ok > 0 & t0_ok >= 0
    if (any(ok_draw)) {
      dt[idx_ok[ok_draw]] <- rtexG(
        sum(ok_draw), mu = mu_ok[ok_draw], sigma = sigma_ok[ok_draw],
        tau = tau_ok[ok_draw], lb = rep(0, sum(ok_draw))
      )
    }
    t0 <- matrix(0, nrow = nr, ncol = n_tri)
    t0[idx_ok[ok_draw]] <- t0_ok[ok_draw]
    finish <- dt + t0
    R <- max.col(-t(finish), ties.method = "first")
    pick <- cbind(R, seq_len(n_tri))
    rt <- finish[pick]
    R <- factor(levels(lR)[R], levels = levels(lR))
    out <- cbind.data.frame(R = R, rt = rt)
    .apply_timed_guess_winner(out, levels(lR))
  }
  list(
    type = "RACE", c_name = "REXG",
    p_types = c(mu = log(.4), sigma = log(.05), tau = log(.1), t0 = log(0),
                pContaminant = qnorm(0), pGuess = qnorm(0)),
    p_types_canonical = c("mu", "sigma", "tau", "t0"),
    transform = list(func = c(mu = "exp", sigma = "exp", tau = "exp", t0 = "exp",
                              pContaminant = "pnorm", pGuess = "pnorm")),
    bound = list(
      minmax = cbind(mu = c(0, Inf), sigma = c(1e-6, Inf), tau = c(1e-6, Inf),
                     t0 = c(0, Inf), pContaminant = c(0.001, 0.999),
                     pGuess = c(0.001, 0.999)),
      exception = c(t0 = 0, pContaminant = 0, pGuess = 0)
    ),
    Ttransform = function(pars, dadm) pars,
    rfun = function(data = NULL, pars) rREXG(data$lR, pars, ok = attr(pars, "ok")),
    dfun = function(rt, pars) dREXG(rt, pars),
    pfun = function(rt, pars) pREXG(rt, pars),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10))
      log_likelihood_race_missing(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
  )
}
