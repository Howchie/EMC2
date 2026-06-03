#' The Ex-Gaussian Race Model
#'
#' Model file for a race model where each accumulator follows an ex-Gaussian
#' distribution with free mean, Gaussian scale, and exponential tail.
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
    cbind.data.frame(R = R, rt = rt)
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
