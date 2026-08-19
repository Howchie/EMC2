# Bounded Ornstein-Uhlenbeck model -- the DDM with leak.
#
# Density, CDF, and simulation wrappers call the C++ kernels.  Boundary-form
# codes and c_name suffixes match the C++ adapter in particle_ll.cpp.

# Boundary-form codes, matching fpe::FPE_BoundaryKind.
.BOU_BND <- c(fixed = 0L, weibull = 1L, exponential = 2L,
              linear_additive = 3L, linear_multiplicative = 4L)

.BOU_SUFFIX <- c(fixed = "", weibull = "_BWEIB", exponential = "_BEXP",
                 linear_additive = "_BLIN_ADD",
                 linear_multiplicative = "_BLIN_MULT")

# Pull the columns required by the C++ entry points.
.bou_cols <- function(pars) {
  get1 <- function(nm, default) {
    if (nm %in% colnames(pars)) pars[, nm] else rep(default, nrow(pars))
  }
  cn <- colnames(pars)
  bkind <- if (!("aInf" %in% cn) || !("tau" %in% cn)) .BOU_BND[["fixed"]]
           else if ("pw" %in% cn) .BOU_BND[["weibull"]]
           else .BOU_BND[["exponential"]]
  list(v = get1("v", 1), a = get1("a", 1), Z = get1("Z", 0.5),
       sv = get1("sv", 0), SZ = get1("SZ", 0),
       t0 = get1("t0", 0), st0 = get1("st0", 0),
       s = get1("s", 1), beta = get1("beta", 0),
       bkind = bkind,
       aInf = if (bkind == 0L) numeric(0) else pars[, "aInf"],
       tau = if (bkind == 0L) numeric(0) else pars[, "tau"],
       pw = if (bkind == 1L) pars[, "pw"] else numeric(0))
}

# The three linear/exponential forms take the same columns, so .bou_cols cannot
# tell them apart from the parameter matrix alone; `kind` carries the model's
# own choice through from BOU(boundary_collapse=) and wins when supplied.
.bou_kind_code <- function(p, kind = NULL) {
  if (is.null(kind)) p$bkind else .BOU_BND[[kind]]
}

.bou_pdf_cdf <- function(rt, R, pars, want_cdf = TRUE, kind = NULL) {
  p <- .bou_cols(pars)
  # lower is the first level of R and upper the second, exactly as for the DDM.
  Ri <- as.integer(R)
  bou_pdf_cdf_vec(rt, Ri, p$v, p$a, p$Z, p$sv, p$SZ, p$t0, p$st0, p$s, p$beta,
                  bkind = .bou_kind_code(p, kind),
                  aInf = p$aInf, tau = p$tau, pw = p$pw,
                  nx = getOption("emc2.bou_nx", 384L),
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
dBOU <- function(rt, R, pars, kind = NULL)
  .bou_pdf_cdf(rt, R, pars, want_cdf = FALSE, kind = kind)$pdf

#' Distribution function of the bounded OU model
#'
#' Defective cdf, matching the convention of [pDDM()]: the two responses' values
#' sum to 1 in the limit, not each on its own.
#'
#' @inheritParams dBOU
#' @return Numeric vector of cumulative probabilities
#' @keywords internal
#' @noRd
pBOU <- function(rt, R, pars, kind = NULL)
  .bou_pdf_cdf(rt, R, pars, want_cdf = TRUE, kind = kind)$cdf

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
rBOU <- function(R, pars, ok = rep(TRUE, length(R)), dt = 1e-4, t_max = 30,
                 kind = NULL) {
  bad <- rep(NA_real_, nrow(pars))
  out <- data.frame(R = rep(NA_integer_, nrow(pars)), rt = bad)
  if (any(ok)) {
    p <- .bou_cols(pars[ok, , drop = FALSE])
    sim <- rbou_cpp(sum(ok), p$v, p$a, p$Z, p$sv, p$SZ, p$t0, p$st0, p$s, p$beta,
                    bkind = .bou_kind_code(p, kind),
                    aInf = p$aInf, tau = p$tau, pw = p$pw,
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
#' When `boundary_collapse` is not `"fixed"` two further parameters appear,
#' `aInf` and `tau` (plus `pw` for `"weibull"`):
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Interpretation** |
#' |---|---|---|---|---|
#' | *aInf* | log | \[0, Inf\] | log(0.5) | Asymptotic separation, in the same units as *a* |
#' | *tau* | log | \[0, Inf\] | log(1) | Time scale of the collapse |
#' | *pw* | log | \[0, Inf\] | log(1) | Shape exponent (Weibull only) |
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
#' @param boundary_collapse Character; the form of the decision boundary. Both
#'   boundaries collapse SYMMETRICALLY toward the midpoint `a/2`, so the
#'   separation runs from `a` at `t = 0` down to `aInf`, with the halves meeting
#'   in the middle. `"fixed"` (the default) leaves the boundaries where they
#'   are; the other forms take `aInf` and `tau`, plus `pw` for `"weibull"`, and
#'   are the same four functional forms [ROU()] offers, written in terms of the
#'   separation \eqn{a(t)}:
#'   \itemize{
#'     \item `"exponential"`: \eqn{a(t) = a_\infty + (a_0 - a_\infty) e^{-t/\tau}}
#'     \item `"linear_additive"`: falls linearly from \eqn{a_0} to
#'       \eqn{a_\infty} over \eqn{[0, \tau]}, then holds
#'     \item `"linear_multiplicative"`: \eqn{a(t) = a_\infty + (a_0 - a_\infty) / (1 + t/\tau)}
#'     \item `"weibull"`: \eqn{a(t) = a_\infty + (a_0 - a_\infty) e^{-(t/\tau)^{pw}}}
#'   }
#'   Collapsing toward the midpoint is the only symmetric choice available:
#'   there is no fixed level to collapse toward, because zero and `a` are the
#'   two response boundaries. The separation is floored at 2% of `a`, so
#'   `aInf = 0` means "collapse essentially all the way" rather than exactly so;
#'   below that width the remaining survivor is absorbed within a millisecond
#'   and the time step, not the mesh, is what stops resolving it.
#' @return A model list with all the necessary functions for EMC2 to sample
#' @examples
#' design_BOU <- design(data = forstmann, model = BOU,
#'                      formula = list(v ~ 0 + S, a ~ E, beta ~ 1, t0 ~ 1,
#'                                     s ~ 1, Z ~ 1, sv ~ 1, SZ ~ 1),
#'                      constants = c(s = log(1)))
#'
#' # Leak and a symmetrically collapsing boundary together
#' design_BOUc <- design(data = forstmann,
#'                       model = function() BOU(boundary_collapse = "exponential"),
#'                       formula = list(v ~ 0 + S, a ~ E, beta ~ 1, t0 ~ 1,
#'                                      s ~ 1, Z ~ 1, aInf ~ 1, tau ~ 1),
#'                       constants = c(s = log(1)))
#' @export
BOU <- function(boundary_collapse = c("fixed", "exponential", "linear_additive",
                                      "linear_multiplicative", "weibull")) {
  boundary_collapse <- match.arg(boundary_collapse)
  kind <- boundary_collapse

  p_types <- c("v" = 1, "a" = log(1), "sv" = log(0), "t0" = log(0),
               "st0" = log(0), "s" = log(1), "Z" = qnorm(0.5),
               "SZ" = qnorm(0), "beta" = log(0))
  transform <- c(v = "identity", a = "exp", sv = "exp", t0 = "exp",
                 st0 = "exp", s = "exp", Z = "pnorm", SZ = "pnorm",
                 beta = "exp")
  minmax <- cbind(v = c(-20, 20), a = c(0, 10), sv = c(.01, 10),
                  t0 = c(0.05, Inf), st0 = c(0, .5), s = c(0, Inf),
                  Z = c(.01, .99), SZ = c(.01, .99), beta = c(0, 50))
  exception <- c(sv = 0, SZ = 0, st0 = 0, beta = 0)

  # Optional columns for the collapsing forms, in the order declared after
  # N_REQ in emc2col::bou (src/col_registry.h).  Only the columns the selected
  # form actually uses are added, so a fixed-bound BOU is unchanged by collapse
  # existing -- and, because they come AFTER the nine required ones, the
  # positional column contract the kernel relies on is untouched.
  if (kind != "fixed") {
    p_types <- c(p_types, aInf = log(0.5), tau = log(1))
    transform <- c(transform, aInf = "exp", tau = "exp")
    # aInf = 0 is a legal interior value: the solver floors the separation at 2%
    # of a, so the operator stays well conditioned all the way down.
    minmax <- cbind(minmax, aInf = c(0, 10), tau = c(1e-3, Inf))
    exception <- c(exception, aInf = 0)
    if (kind == "weibull") {
      p_types <- c(p_types, pw = log(1))
      transform <- c(transform, pw = "exp")
      minmax <- cbind(minmax, pw = c(1e-3, Inf))
    }
  }

  # pContaminant (omission) and pGuess (uniform outlier).  Appended LAST, after
  # the optional collapse columns, because the kernel indexes aInf/tau/pw
  # positionally from N_REQ; see add_nuisance_pars() and contaminant_mixture.h.
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- .nuis$p_types; transform <- .nuis$transform
  minmax <- .nuis$minmax; exception <- .nuis$exception

  list(
    c_name = paste0("BOU", .BOU_SUFFIX[[kind]]),
    type = "DDM",
    # Fokker-Planck solve cached per parameter tuple: rt only picks readout
    # points off a march that has already been paid for, so binning it saves
    # nothing and costs accuracy.  This also keeps make_emc's default
    # rt_resolution = 1/60 from evaluating the density at floored RTs.
    compress_ok = FALSE,
    # ORDER IS THE KERNEL'S COLUMN ORDER, not a stylistic choice: the C++ side
    # reads columns positionally and validate_col_prefix checks this against
    # emc2col::bou::spec() in src/col_registry.h.  It is DDM's order with beta
    # appended, which is what lets a DDM design convert by adding beta~1.
    p_types = p_types,
    p_types_canonical = setdiff(names(p_types), .nuisance_par_names),
    transform = list(func = transform),
    bound = list(minmax = minmax, exception = exception),
    # Same trial-dependent transform as the DDM, verbatim, so the two models are
    # directly comparable and a design converts between them by adding beta.
    Ttransform = function(pars, dadm) {
      pars[, "SZ"] <- 2 * pars[, "SZ"] * pmin(pars[, "Z"], 1 - pars[, "Z"])
      pars <- cbind(pars, z = pars[, "Z"] * pars[, "a"],
                    sz = pars[, "SZ"] * pars[, "a"])
      pars
    },
    rfun = function(data = NULL, pars) rBOU(data$R, pars, attr(pars, "ok"),
                                            kind = kind),
    dfun = function(rt, R, pars) dBOU(rt, R, pars, kind = kind),
    pfun = function(rt, R, pars) pBOU(rt, R, pars, kind = kind),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      log_likelihood_ddm(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}
