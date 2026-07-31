# The bounded Ornstein-Uhlenbeck model -- the DDM with leak.
#
# Wrappers only: the density, cdf and simulator all live in C++ (see
# src/bou_diffusion.cpp and src/fpe_bou.h), and the likelihood is the shared
# DDM one.  Nothing here reimplements the model, so the R-side density and the
# sampled likelihood cannot disagree.

# Boundary-form codes, matching fpe::FPE_BoundaryKind in src/fpe_models.h.
.BOU_BND <- c(fixed = 0L, weibull = 1L, exponential = 2L,
              linear_additive = 3L, linear_multiplicative = 4L)

# Pull the columns the C++ entry points want out of a parameter matrix, with the
# defaults a design may legitimately omit.
.bou_cols <- function(pars) {
  get1 <- function(nm, default) {
    if (nm %in% colnames(pars)) pars[, nm] else rep(default, nrow(pars))
  }
  list(v = get1("v", 1), a = get1("a", 1), Z = get1("Z", 0.5),
       sv = get1("sv", 0), SZ = get1("SZ", 0),
       t0 = get1("t0", 0), st0 = get1("st0", 0),
       s = get1("s", 1), beta = get1("beta", 0))
}

.bou_pdf_cdf <- function(rt, R, pars, want_cdf = TRUE) {
  p <- .bou_cols(pars)
  # lower is the first level of R and upper the second, exactly as for the DDM.
  Ri <- as.integer(R)
  bou_pdf_cdf_vec(rt, Ri, p$v, p$a, p$Z, p$sv, p$SZ, p$t0, p$st0, p$s, p$beta,
                  nx = getOption("emc2.bou_nx", 512L),
                  dt_target = getOption("emc2.bou_dt", 5e-4),
                  grade = getOption("emc2.bou_grade", 1),
                  tgrade = getOption("emc2.bou_tgrade", 32),
                  n_sv = getOption("emc2.bou_n_sv", 7L),
                  n_sz = getOption("emc2.bou_n_sz", 7L),
                  n_st0 = getOption("emc2.bou_n_st0", 7L),
                  want_cdf = want_cdf)
}

#' Density of the bounded OU model
#'
#' Defective density: the value returned for a response is that response's
#' contribution, so the two responses' densities integrate to 1 between them.
#' Matches the convention of [dDDM()].
#'
#' @param rt Numeric vector of response times
#' @param R Response factor; the first level is the lower boundary
#' @param pars Matrix of parameters, one row per element of `rt`
#' @return Numeric vector of densities
#' @keywords internal
#' @noRd
dBOU <- function(rt, R, pars) .bou_pdf_cdf(rt, R, pars, want_cdf = FALSE)$pdf

#' Distribution function of the bounded OU model
#'
#' Defective cdf, matching the convention of [pDDM()]: the two responses' values
#' sum to 1 in the limit, not each on its own.
#'
#' @inheritParams dBOU
#' @return Numeric vector of cumulative probabilities
#' @keywords internal
#' @noRd
pBOU <- function(rt, R, pars) .bou_pdf_cdf(rt, R, pars, want_cdf = TRUE)$cdf

#' Random generation for the bounded OU model
#'
#' Simulates by stepping the exact Ornstein-Uhlenbeck transition (not an Euler
#' discretisation of the drift), testing both absorbing boundaries each step and
#' applying a Brownian-bridge correction for paths that cross and return within
#' a step. Without that correction RTs are biased upwards and error rates
#' downwards, which would make the simulator useless as a check on the solver.
#'
#' @param R Response factor supplying the levels of the result
#' @param pars Matrix of parameters, one row per trial
#' @param ok Logical vector of rows to simulate
#' @param dt Integration step in seconds
#' @param t_max Horizon after which a trial is recorded as unfinished
#' @return A data frame with columns `R` and `rt`
#' @keywords internal
#' @noRd
rBOU <- function(R, pars, ok = rep(TRUE, length(R)), dt = 1e-4, t_max = 30) {
  bad <- rep(NA_real_, nrow(pars))
  out <- data.frame(R = rep(NA_integer_, nrow(pars)), rt = bad)
  if (any(ok)) {
    p <- .bou_cols(pars[ok, , drop = FALSE])
    sim <- rbou_cpp(sum(ok), p$v, p$a, p$Z, p$sv, p$SZ, p$t0, p$st0, p$s, p$beta,
                    dt = dt, t_max = t_max)
    out[ok, "R"] <- sim$R
    out[ok, "rt"] <- sim$rt
  }
  data.frame(R = factor(out$R, levels = c(1, 2), labels = levels(R)),
             rt = out$rt)
}

#' The Bounded Ornstein-Uhlenbeck Model
#'
#' The two-choice OU diffusion of Smith and Ratcliff (2004): the diffusion
#' decision model with a leak term. Evidence accumulates between two absorbing
#' boundaries, and the drift is pulled back toward the starting point at rate
#' `beta`. At `beta = 0` the model is exactly the Wiener diffusion of [DDM()],
#' which is how it is validated -- but it reaches that answer through a
#' Fokker-Planck solver rather than the DDM's series expansion.
#'
#' Model files are almost exclusively used in `design()`.
#'
#' @details
#'
#' Default values are used for all parameters that are not explicitly listed in
#' the `formula` argument of `design()`. They can also be accessed with
#' `BOU()$p_types`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *v* | - | \[-Inf, Inf\] | 1 | | Mean evidence-accumulation rate (drift rate) |
#' | *a* | log | \[0, Inf\] | log(1) | | Boundary separation |
#' | *beta* | log | \[0, Inf\] | log(0) | | Leak: decay toward the starting point. 0 gives the DDM |
#' | *t0* | log | \[0, Inf\] | log(0) | | Non-decision time |
#' | *s* | log | \[0, Inf\] | log(1) | | Within-trial standard deviation of drift rate |
#' | *Z* | probit | \[0, 1\] | qnorm(0.5) | *z* = *Z* x *a* | Relative start point (bias) |
#' | *SZ* | probit | \[0, 1\] | qnorm(0) | *sz* = 2 x *SZ* x min(*a* x *Z*, *a* x (1-*Z*)) | Relative between-trial variation in start point |
#' | *sv* | log | \[0, Inf\] | log(0) | | Between-trial standard deviation of drift rate |
#' | *st0* | log | \[0, Inf\] | log(0) | | Between-trial variation (range) in non-decision time |
#'
#' The parameterisation is [DDM()]'s with `beta` added, so a DDM design converts
#' to this model by adding `beta ~ 1`, and `beta` fixed to 0 recovers the DDM.
#'
#' Note that the decay pulls the process toward the *starting point*, following
#' Smith and Ratcliff's own choice (their footnote 2). Decay toward zero is not
#' an available alternative here, because zero is one of the response
#' boundaries.
#'
#' Estimation cost is dominated by across-trial variability rather than by the
#' diffusion: one Fokker-Planck march yields both response distributions, but
#' because the decay is anchored at the starting point, `sv` and `SZ` each add a
#' quadrature dimension whose nodes need their own march. With both switched on
#' the model costs roughly `n_sv * n_SZ` solves per distinct parameter set; node
#' counts are controlled by the options `emc2.bou_n_sv` and `emc2.bou_n_sz`.
#' `st0` is a convolution and is nearly free.
#'
#' Smith, P. L., & Ratcliff, R. (2004). A comparison of sequential sampling
#' models for two-choice reaction time. *Psychological Review, 111*(2), 333-367.
#' doi:10.1037/0033-295X.111.2.333
#'
#' @return A model list with all the necessary functions for EMC2 to sample
#' @examples
#' design_BOU <- design(data = forstmann, model = BOU,
#'                      formula = list(v ~ 0 + S, a ~ E, beta ~ 1, t0 ~ 1,
#'                                     s ~ 1, Z ~ 1, sv ~ 1, SZ ~ 1),
#'                      constants = c(s = log(1)))
#' @export
BOU <- function() {
  list(
    c_name = "BOU",
    type = "DDM",
    p_types = c("v" = 1, "a" = log(1), "beta" = log(0), "sv" = log(0),
                "t0" = log(0), "st0" = log(0), "s" = log(1),
                "Z" = qnorm(0.5), "SZ" = qnorm(0)),
    transform = list(func = c(v = "identity", a = "exp", beta = "exp",
                              sv = "exp", t0 = "exp", st0 = "exp", s = "exp",
                              Z = "pnorm", SZ = "pnorm")),
    bound = list(minmax = cbind(v = c(-20, 20), a = c(0, 10), beta = c(0, 50),
                                Z = c(.01, .99), t0 = c(0.05, Inf),
                                sv = c(.01, 10), s = c(0, Inf),
                                SZ = c(.01, .99), st0 = c(0, .5)),
                 exception = c(sv = 0, SZ = 0, st0 = 0, beta = 0)),
    # Same trial-dependent transform as the DDM, verbatim, so the two models are
    # directly comparable and a design converts between them by adding beta.
    Ttransform = function(pars, dadm) {
      pars[, "SZ"] <- 2 * pars[, "SZ"] * pmin(pars[, "Z"], 1 - pars[, "Z"])
      pars <- cbind(pars, z = pars[, "Z"] * pars[, "a"],
                    sz = pars[, "SZ"] * pars[, "a"])
      pars
    },
    rfun = function(data = NULL, pars) rBOU(data$R, pars, attr(pars, "ok")),
    dfun = function(rt, R, pars) dBOU(rt, R, pars),
    pfun = function(rt, R, pars) pBOU(rt, R, pars),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      log_likelihood_ddm(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}
