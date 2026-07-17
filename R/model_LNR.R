
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


#' The Shifted-Gamma Race Model
#'
#' Race model in which each accumulator has a shifted Gamma finish-time
#' distribution. For accumulator `i`,
#' `T_i = shift_i + G_i`, with
#' `G_i ~ Gamma(shape_i, rate = lambda_i)`. The observed response is the
#' accumulator with the smallest finish time; the race likelihood combines
#' the winning density with the survivor probabilities of all losing
#' accumulators.
#'
#' The model parameters are:
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Interpretation** |
#' |---|---|---|---|---|
#' | *lambda* | log | \[0, Inf\] | log(1) | Gamma rate; the reciprocal is the scale parameter conditional on `shape`. |
#' | *shape* | log | \[0, Inf\] | log(1) | Gamma shape. `shape = 1` is exponential and `shape = 2` is Erlang-2. |
#' | *shift* | log | \[0, Inf\] | log(0) | Nonnegative additive shift of the finish-time distribution. |
#' | *pContaminant* | probit | \[0, 1\] | qnorm(0) | Optional contamination probability handled by the data pipeline. |
#'
#' All parameters are estimated on their indicated transformed scales. The
#' zero shift is an explicit exception to the positive log-scale bounds, so
#' `shift = 0` can be fixed in a design. Similarly, `shape` can be fixed to 1
#' or 2 through model constants to obtain exponential/Erlang special cases.
#' The finish-time distribution has support above `shift`, and its conditional
#' mean is `shift + shape / lambda`. EMC2 creates one accumulator per response
#' level in `R` and a latent accumulator factor `lR` for the race.
#'
#' @return A model list compatible with `design()`.
#' @export
RGAMMA <- function() {
  dRGAMMA <- function(rt, pars) {
    out <- numeric(length(rt))
    tt <- rt - pars[, "shift"]
    ok <- tt > 0 & is.finite(tt) & pars[, "lambda"] > 0 & pars[, "shape"] > 0
    ok[is.na(ok)] <- FALSE
    out[ok] <- stats::dgamma(tt[ok], shape = pars[ok, "shape"], rate = pars[ok, "lambda"])
    out
  }

  pRGAMMA <- function(rt, pars) {
    out <- numeric(length(rt))
    tt <- rt - pars[, "shift"]
    ok <- tt > 0 & is.finite(tt) & pars[, "lambda"] > 0 & pars[, "shape"] > 0
    ok[is.na(ok)] <- FALSE
    out[ok] <- stats::pgamma(tt[ok], shape = pars[ok, "shape"], rate = pars[ok, "lambda"])
    out
  }

  rRGAMMA <- function(lR, pars, p_types = c("lambda", "shape", "shift"), ok = rep(TRUE, dim(pars)[1])) {
    if (!all(p_types %in% dimnames(pars)[[2]])) {
      stop("pars must have columns ", paste(p_types, collapse = " "))
    }
    nr <- length(levels(lR))
    dt <- matrix(Inf, ncol = nrow(pars) / nr, nrow = nr)
    idx_ok <- which(ok)
    pars <- pars[idx_ok, , drop = FALSE]
    lam <- pars[, "lambda"]
    shp <- pars[, "shape"]
    shf <- pars[, "shift"]
    ok_draw <- is.finite(lam) & is.finite(shp) & is.finite(shf) & lam > 0 & shp > 0
    if (any(ok_draw)) {
      dt[idx_ok[ok_draw]] <- stats::rgamma(sum(ok_draw), shape = shp[ok_draw], rate = lam[ok_draw]) + shf[ok_draw]
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
    c_name = "RGAMMA",
    p_types = c("lambda" = log(1), "shape" = log(1), "shift" = log(0), "pContaminant" = qnorm(0)),
    p_types_canonical = c("lambda", "shape", "shift"),
    transform = list(func = c(lambda = "exp", shape = "exp", shift = "exp", pContaminant = "pnorm")),
    bound = list(
      minmax = cbind(lambda = c(1e-6, Inf), shape = c(1e-6, Inf), shift = c(0, Inf), pContaminant = c(0.001, 0.999)),
      exception = c(shift = 0, pContaminant = 0)
    ),
    Ttransform = function(pars, dadm) pars,
    rfun = function(data = NULL, pars) rRGAMMA(data$lR, pars, ok = attr(pars, "ok")),
    dfun = function(rt, pars) dRGAMMA(rt, pars),
    pfun = function(rt, pars) pRGAMMA(rt, pars),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      log_likelihood_race_missing(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
    }
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
