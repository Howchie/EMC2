# ============================================================================
# OPR: opponent Poisson race (Math/CounterAndQuorumModels.tex, sec:opponent)
#
# Each accumulator integrates X(u) = N+(u) - N-(u), two independent Poisson
# event streams with rates lambda_p and lambda_m, and finishes when it reaches
# the continuous criterion K.  The kernels (src/model_OPR.h) take the rates;
# the (v, rho) chart maps lambda_p = (1 + rho) v, lambda_m = rho v, and the
# (v, lambda0) chart lambda_p = v + lambda0, lambda_m = lambda0.
# ============================================================================

# Event rates per row from either chart.
.opr_rates <- function(pars) {
  cn <- colnames(pars)
  if (all(c("lambda_p", "lambda_m") %in% cn))
    return(list(a = pars[, "lambda_p"], b = pars[, "lambda_m"]))
  if (all(c("v", "lambda0") %in% cn))
    return(list(a = pars[, "v"] + pars[, "lambda0"], b = pars[, "lambda0"]))
  if (all(c("v", "rho") %in% cn))
    return(list(a = (1 + pars[, "rho"]) * pars[, "v"],
                b = pars[, "rho"] * pars[, "v"]))
  stop("OPR requires parameter columns v and rho, v and lambda0, or ",
       "lambda_p and lambda_m.")
}

.opr_check_cols <- function(pars) {
  missing <- setdiff(c("K", "t0"), colnames(pars))
  if (length(missing))
    stop("OPR requires parameter columns ", paste(missing, collapse = ", "))
  invisible(NULL)
}

dOPR <- function(rt, pars) {
  .opr_check_cols(pars)
  r <- .opr_rates(pars)
  dt <- rt - pars[, "t0"]
  ok <- is.finite(dt) & dt > 0 & is.finite(pars[, "t0"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- dopr(dt[ok], pars[ok, "K"], r$a[ok], r$b[ok])
    out[ok][is.na(out[ok])] <- 0
  }
  out
}

# rt = Inf is kept: the CDF there is the completion probability h.
pOPR <- function(rt, pars) {
  .opr_check_cols(pars)
  r <- .opr_rates(pars)
  dt <- rt - pars[, "t0"]
  ok <- !is.na(dt) & dt > 0 & is.finite(pars[, "t0"])
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- popr(dt[ok], pars[ok, "K"], r$a[ok], r$b[ok])
    out[ok][is.na(out[ok])] <- 0
  }
  out
}

# Reference simulator: inverts the analytic survivor per accumulator.
rOPR <- function(lR, pars, ok = rep(TRUE, nrow(pars))) {
  .opr_check_cols(pars)
  nr <- length(levels(lR))
  if (nr < 1L || nrow(pars) %% nr != 0L)
    stop("OPR: pars rows must be a multiple of accumulators")
  ok <- as.logical(ok)
  ok[is.na(ok)] <- FALSE
  r <- .opr_rates(pars)
  K <- pars[, "K"]
  t0 <- pars[, "t0"]
  dt <- rep(Inf, nrow(pars))
  for (i in which(ok & is.finite(t0) & is.finite(K) & K > 0 & r$a > 0 &
                    r$b >= 0)) {
    a <- r$a[i]; b <- r$b[i]
    if (a < b) {
      if (runif(1) >= (a / b)^K[i]) next
      tmp <- a; a <- b; b <- tmp
    }
    lv <- log(runif(1))
    g <- function(y) popr(exp(y), K[i], a, b, lower_tail = FALSE,
                          log_out = TRUE) - lv
    y0 <- log(K[i] / a)
    lo <- y0 - 1; while (g(lo) <= 0) lo <- lo - 1
    hi <- y0 + 1; while (g(hi) >= 0) hi <- hi + 1
    dt[i] <- exp(stats::uniroot(g, c(lo, hi), tol = 1e-12)$root)
  }
  n_tri <- nrow(pars) / nr
  finish <- matrix(dt + ifelse(is.finite(t0), t0, 0), nrow = nr, ncol = n_tri)
  R <- max.col(-t(finish), ties.method = "first")
  rt <- finish[cbind(R, seq_len(n_tri))]
  R <- factor(levels(lR)[R], levels = levels(lR))
  R[!is.finite(rt)] <- NA
  out <- cbind.data.frame(R = R, rt = rt)
  .apply_timed_guess_winner(out, levels(lR))
}

#' The Opponent Poisson Race Model (OPR)
#'
#' Race model in which each accumulator counts signed evidence events. For
#' accumulator `i`, positive and negative events arrive as independent Poisson
#' processes with rates `lambda_p` and `lambda_m`, and the evidence
#' `X(u) = N+(u) - N-(u)` starts at zero. The accumulator finishes at the
#' continuous criterion `K > 0` with the first-passage subdensity
#'
#' `f(u) = (K / u) (lambda_p / lambda_m)^(K / 2) exp(-(lambda_p + lambda_m) u) I_K(2 u sqrt(lambda_p lambda_m))`,
#'
#' with `I_K` the modified Bessel function of the first kind, and CDF
#'
#' `F(u) = H(2 lambda_p u; 2K, 2 lambda_m u) + (lambda_p / lambda_m)^K H(2 lambda_m u; 2K + 2, 2 lambda_p u)`,
#'
#' with `H(x; d, ncp)` the noncentral chi-squared CDF. At integer `K` this is
#' the first passage of the signed event count to `K`. The completion
#' probability is `h = F(Inf) = 1` when `lambda_p >= lambda_m` and
#' `(lambda_p / lambda_m)^K` otherwise; the remaining mass is an accumulator
#' that never finishes. `lambda_m = 0` gives the `Gamma(K, lambda_p)` first
#' passage. The finish time is `t0 + u`, and the observed response is the
#' accumulator with the smallest finish time.
#'
#' The default `parameterization = "ratio"` samples the net rate `v` and the
#' opposing-activity ratio `rho`, with `lambda_p = (1 + rho) v` and
#' `lambda_m = rho v`:
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *K* | log | \[1, Inf\] | log(3) | | Criterion (evidence units). |
#' | *v* | log | \[1e-4, Inf\] | log(10) | | Net event rate `lambda_p - lambda_m`. |
#' | *rho* | log | \[0, Inf\] | log(1) | | Negative event rate relative to `v`. |
#' | *t0* | log | \[0.05, Inf\] | log(0) | | Non-decision time. |
#'
#' `parameterization = "background"` samples the net rate `v` and a background
#' event rate `lambda0` shared by both streams, with `lambda_p = v + lambda0`
#' and `lambda_m = lambda0`:
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *K* | log | \[1, Inf\] | log(3) | | Criterion (evidence units). |
#' | *v* | log | \[1e-4, Inf\] | log(10) | | Net event rate `lambda_p - lambda_m`. |
#' | *lambda0* | log | \[0, Inf\] | log(10) | | Background event rate in each stream. |
#' | *t0* | log | \[0.05, Inf\] | log(0) | | Non-decision time. |
#'
#' `parameterization = "rates"` samples the event rates directly:
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *K* | log | \[1, Inf\] | log(3) | | Criterion (evidence units). |
#' | *lambda_p* | log | \[1e-4, Inf\] | log(20) | | Positive event rate. |
#' | *lambda_m* | log | \[0, Inf\] | log(10) | | Negative event rate. |
#' | *t0* | log | \[0.05, Inf\] | log(0) | | Non-decision time. |
#'
#' The `"ratio"` and `"background"` charts have `h = 1`. Under `"rates"` an accumulator with
#' `lambda_m > lambda_p` is defective.
#'
#' `mapped_pars()` adds `lambda_p` and `lambda_m` under `"ratio"` and
#' `"background"`, and `v`
#' (`lambda_p - lambda_m`) and `h` under `"rates"`.
#'
#' Optional nuisance parameters are `pContaminant`, an omission probability
#' outside the accumulators, and `pGuess`, a uniform-outlier probability.
#'
#' @param parameterization `"ratio"` (the default) samples *v* and *rho*;
#'   `"background"` samples *v* and *lambda0*; `"rates"` samples *lambda_p*
#'   and *lambda_m*.
#' @return A model list defining the OPR race model.
#' @examples
#' matchfun <- function(d) d$S == d$lR
#' design_OPR <- design(data = forstmann, model = OPR, matchfun = matchfun,
#'                      formula = list(K ~ E, v ~ lM, rho ~ 1, t0 ~ 1))
#' design_OPR_rates <- design(data = forstmann,
#'                            model = function() OPR(parameterization = "rates"),
#'                            matchfun = matchfun,
#'                            formula = list(K ~ E, lambda_p ~ lM,
#'                                           lambda_m ~ lM, t0 ~ 1))
#' design_OPR_bg <- design(data = forstmann,
#'                         model = function() OPR(parameterization = "background"),
#'                         matchfun = matchfun,
#'                         formula = list(K ~ E, v ~ lM, lambda0 ~ 1, t0 ~ 1))
#' @export
OPR <- function(parameterization = c("ratio", "background", "rates")) {
  parameterization <- match.arg(parameterization)
  rates <- parameterization == "rates"
  r1 <- if (rates) "lambda_p" else "v"
  r2 <- switch(parameterization, ratio = "rho", background = "lambda0",
               rates = "lambda_m")
  p_types <- c(K = log(3), r1 = if (rates) log(20) else log(10),
               r2 = switch(parameterization, ratio = log(1), log(10)),
               t0 = log(0))
  transform <- c(K = "exp", r1 = "exp", r2 = "exp", t0 = "exp")
  # K >= 1 keeps the onset density u^(K - 1) finite.  lambda_m = 0 (rho = 0,
  # lambda0 = 0) is the Gamma first passage, the continuous limit of the kernel.
  minmax <- cbind(K = c(1, Inf), r1 = c(1e-4, Inf), r2 = c(0, Inf),
                  t0 = c(0.05, Inf))
  exception <- c(t0 = 0)
  .rn <- function(x) { x[x == "r1"] <- r1; x[x == "r2"] <- r2; x }
  names(p_types) <- .rn(names(p_types))
  names(transform) <- .rn(names(transform))
  colnames(minmax) <- .rn(colnames(minmax))
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- .nuis$p_types; transform <- .nuis$transform
  minmax <- .nuis$minmax; exception <- .nuis$exception
  chart <- switch(parameterization, ratio = 0L, rates = 1L, background = 2L)

  list(
    type = "RACE",
    c_name = switch(parameterization, ratio = "OPR", rates = "OPR_RATES",
                    background = "OPR_BG"),
    p_types = p_types,
    p_types_canonical = c("K", r1, r2, "t0"),
    transform = list(func = transform),
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) {
      # Reporting only: the compiled likelihood maps the chart in the kernel.
      r <- .opr_rates(pars)
      add <- if (rates) {
        cbind(v = as.numeric(r$a - r$b),
              h = as.numeric(ifelse(r$a >= r$b, 1, (r$a / r$b)^pars[, "K"])))
      } else {
        cbind(lambda_p = as.numeric(r$a), lambda_m = as.numeric(r$b))
      }
      rownames(add) <- NULL
      cbind(pars, add)
    },
    rfun = function(data = NULL, pars) .rfun_OPR(data$lR, pars, chart,
                                                 ok = attr(pars, "ok")),
    dfun = function(rt, pars) dOPR(rt, pars),
    pfun = function(rt, pars) pOPR(rt, pars),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      stop("OPR: the likelihood is implemented in the compiled race path; ",
           "the R likelihood route is not supported.")
    }
  )
}
