# Gompertz growth-process race

.GOM_BND <- c(fixed = 0L, weibull = 1L, exponential = 2L,
              linear_additive = 3L, linear_multiplicative = 4L)
.GOM_SUFFIX <- c(fixed = "", weibull = "_BWEIB", exponential = "_BEXP",
                 linear_additive = "_BLIN_ADD",
                 linear_multiplicative = "_BLIN_MULT")

.gom_cols <- function(pars, kind = NULL) {
  cn <- colnames(pars)
  bkind <- if (is.null(kind)) {
    if (!all(c("Binf", "tau") %in% cn)) .GOM_BND[["fixed"]]
    else if ("pw" %in% cn) .GOM_BND[["weibull"]]
    else .GOM_BND[["exponential"]]
  } else .GOM_BND[[kind]]
  list(
    alpha = pars[, "alpha"], beta = pars[, "beta"], K = pars[, "K"],
    B = pars[, "B"], A = pars[, "A"], t0 = pars[, "t0"],
    bkind = bkind,
    Binf = if (bkind == 0L) numeric(0) else pars[, "Binf"],
    tau = if (bkind == 0L) numeric(0) else pars[, "tau"],
    pw = if (bkind == 1L) pars[, "pw"] else numeric(0)
  )
}

.gom_pdf_cdf <- function(rt, pars, kind = NULL) {
  p <- .gom_cols(pars, kind)
  g <- .rou_grid()
  n <- length(rt)
  # The C++ helper groups identical physical parameter tuples and uses the
  # same FPE cache as the sampled likelihood.  Invalid/non-finite rows are
  # replaced only for the call; their output is set back to zero below.
  bad_par <- !is.finite(p$alpha) | !is.finite(p$beta) | !is.finite(p$K) |
             !is.finite(p$B) | !is.finite(p$A) | !is.finite(p$t0) |
             p$alpha <= 0 | p$beta <= 0 | p$K <= 0 |
             p$B < 0 | p$A < 0 | p$t0 < 0
  bad_rt <- is.na(rt) | !is.finite(rt) | (is.finite(rt) & rt <= p$t0)
  bad <- bad_par | bad_rt
  rt_pos_inf <- is.infinite(rt) & rt > 0
  rt_call <- rt
  if (any(bad)) rt_call[bad] <- p$t0[bad] + 1
  p$alpha[bad] <- 1; p$beta[bad] <- 1; p$K[bad] <- 5
  p$B[bad] <- 1; p$A[bad] <- 0; p$t0[bad] <- 0
  if (p$bkind != 0L) {
    p$Binf[bad] <- 1.5; p$tau[bad] <- 1
    if (length(p$pw)) p$pw[bad] <- 1
  }
  out <- gomp_pdf_cdf_vec(rt_call, p$alpha, p$beta, p$K, p$B, p$A, p$t0,
                          p$bkind, p$Binf, p$tau, p$pw,
                          as.integer(g$nx), g$dt_target, g$grade, g$tgrade)
  if (any(bad)) {
    out$pdf[bad] <- 0
    out$cdf[bad] <- ifelse(!bad_par[bad] & rt_pos_inf[bad], 1, 0)
  }
  out
}

dGOM <- function(rt, pars, kind = NULL)
  .gom_pdf_cdf(rt, pars, kind)$pdf

pGOM <- function(rt, pars, kind = NULL)
  .gom_pdf_cdf(rt, pars, kind)$cdf

.gomp_boundary <- function(t, b0, binf, tau, pw, kind) {
  if (kind == "fixed") return(b0)
  if (kind == "weibull") return(binf + (b0 - binf) * exp(-(t / tau)^pw))
  if (kind == "exponential") return(binf + (b0 - binf) * exp(-t / tau))
  if (kind == "linear_additive") return(binf + (b0 - binf) * pmax(0, 1 - t / tau))
  binf + (b0 - binf) / (1 + t / tau)
}

.gomp_hit_one <- function(alpha, beta, K, B, A, kind, Binf, tau, pw,
                          dt, t_max) {
  if (!all(is.finite(c(alpha, beta, K, B, A, dt, t_max))) ||
      alpha <= 0 || beta <= 0 || K <= 0 || B < 0 || A < 0 ||
      dt <= 0 || t_max <= 0) return(Inf)
  b0 <- 1 + B + A
  X <- 1 + if (A > 0) runif(1, 0, A) else 0
  b <- .gomp_boundary(0, b0, if (kind == "fixed") b0 else Binf, tau, pw, kind)
  if (X >= b) return(0)
  y <- log(X)
  v <- alpha * log(K) - 0.5 * beta^2
  phi <- exp(-alpha * dt)
  gain <- -expm1(-alpha * dt) / alpha
  var <- beta^2 * (-expm1(-2 * alpha * dt)) / (2 * alpha)
  sd <- sqrt(max(var, 0))
  inv2 <- 2 / (beta^2 * dt)
  t <- 0
  while (t < t_max) {
    y1 <- y * phi + v * gain + sd * rnorm(1)
    b1 <- .gomp_boundary(t + dt, b0, if (kind == "fixed") b0 else Binf,
                         tau, pw, kind)
    lb <- log(b); lb1 <- log(b1)
    t1 <- t + dt
    if (y1 >= lb1) {
      frac <- (lb - y) / max((lb - y) - (lb1 - y1), 1e-300)
      return(t + min(max(frac, 0), 1) * dt)
    }
    if (runif(1) < exp(-(lb - y) * (lb1 - y1) * inv2))
      return(t + runif(1) * dt)
    y <- y1; b <- b1; t <- t1
  }
  Inf
}

.rGOM_R <- function(lR, pars, ok = rep(TRUE, nrow(pars)), kind = "fixed",
                    dt = getOption("emc2.gom_sim_dt", 1e-3),
                    t_max = getOption("emc2.gom_sim_tmax", 30)) {
  if (!is.null(attr(pars, "ok"))) ok <- attr(pars, "ok")
  if (is.null(ok)) ok <- rep(TRUE, nrow(pars))
  nr <- length(levels(lR)); n_trials <- nrow(pars) / nr
  out <- data.frame(R = factor(rep(NA, n_trials), levels = levels(lR)),
                    rt = rep(Inf, n_trials))
  hits <- matrix(Inf, nr, n_trials)
  cn <- colnames(pars)
  has_binf <- "Binf" %in% cn
  has_tau <- "tau" %in% cn
  has_pw <- "pw" %in% cn
  for (i in which(ok)) {
    tr <- ((i - 1L) %% nr) + 1L
    j <- ((i - 1L) %/% nr) + 1L
    hits[tr, j] <- .gomp_hit_one(pars[i, "alpha"], pars[i, "beta"],
                                 pars[i, "K"], pars[i, "B"], pars[i, "A"],
                                 kind, if (has_binf) pars[i, "Binf"] else 0,
                                 if (has_tau) pars[i, "tau"] else 1,
                                 if (has_pw) pars[i, "pw"] else 1,
                                 dt, t_max)
  }
  win <- max.col(-t(hits), ties.method = "first")
  pick <- cbind(win, seq_len(n_trials))
  out$R <- factor(levels(lR)[win], levels = levels(lR))
  out$rt <- matrix(pars[, "t0"], nr, n_trials)[pick] + hits[pick]
  bad <- colSums(is.finite(hits)) == 0L
  out$R[bad] <- NA; out$rt[bad] <- Inf
  .apply_timed_guess_winner(out, levels(lR))
}

#' The Gompertz Growth Process race model
#'
#' Each accumulator follows `dX = alpha X log(K/X) dt + beta X dW`, starts at
#' `1 + Uniform(0, A)`, and responds at `1 + B + A`.  The likelihood is evaluated
#' by the cached OU Fokker--Planck solver after the exact `Y = log(X)` reduction.
#' The physical start range is retained when seeding the log-space solver.
#'
#' @param boundary_collapse Character; one of `"fixed"`, `"exponential"`,
#'   `"linear_additive"`, `"linear_multiplicative"`, or `"weibull"`.
#' @return A list defining an EMC2 race model.
#' @export
GOM <- function(boundary_collapse = c("fixed", "exponential", "linear_additive",
                                      "linear_multiplicative", "weibull")) {
  boundary_collapse <- match.arg(boundary_collapse)
  kind <- boundary_collapse
  p_types <- c(alpha = log(1), beta = log(1), K = log(5), B = log(1),
               A = log(0), t0 = log(0))
  transform <- c(alpha = "exp", beta = "exp", K = "exp", B = "exp",
                 A = "exp", t0 = "exp")
  minmax <- cbind(alpha = c(1e-4, Inf), beta = c(1e-4, Inf), K = c(1e-4, Inf),
                  B = c(0, Inf), A = c(1e-4, Inf), t0 = c(0.05, Inf))
  exception <- c(A = 0)
  if (kind != "fixed") {
    p_types <- c(p_types, Binf = log(1.5), tau = log(1))
    transform <- c(transform, Binf = "exp", tau = "exp")
    minmax <- cbind(minmax, Binf = c(1e-4, Inf), tau = c(1e-3, Inf))
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
    c_name = paste0("GOM", .GOM_SUFFIX[[kind]]),
    compress_ok = FALSE,
    p_types = p_types,
    p_types_canonical = setdiff(names(p_types), .nuisance_par_names),
    transform = list(func = transform),
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) cbind(pars, b = 1 + pars[, "B"] + pars[, "A"]),
    rfun = function(data = NULL, pars) .rGOM_R(data$lR, pars,
                                                ok = attr(pars, "ok"), kind = kind),
    dfun = function(rt, pars) dGOM(rt, pars, kind = kind),
    pfun = function(rt, pars) pGOM(rt, pars, kind = kind),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10))
      log_likelihood_race_missing(pars = pars, dadm = dadm, model = model,
                                   min_ll = min_ll)
  )
}

# Descriptive aliases; GOM is the canonical short model name used in c_name.
#' @export
GOMPERTZ <- GOM
#' @export
Gompertz <- GOM
