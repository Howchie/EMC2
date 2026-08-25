# Racing symmetric alpha-stable Lévy-flight model.

.rlf_grid <- function() {
  list(
    nx = getOption("emc2.rlf_nx", 160L),
    dt_target = getOption("emc2.rlf_dt", 1.6e-2),
    tgrade = getOption("emc2.rlf_tgrade", getOption("emc2.fpe_tgrade", 1)),
    adaptive = isTRUE(getOption("emc2.rlf_adaptive", FALSE)),
    explicit_inverse = isTRUE(
      getOption("emc2.rlf_explicit_inverse", TRUE)
    ),
    sparse_output = isTRUE(getOption("emc2.rlf_sparse_output", TRUE)),
    # Off by default: see the note on Grid::simd_batch in src/model_RLF.h.
    simd_batch = isTRUE(getOption("emc2.rlf_simd_batch", FALSE)),
    horizon_split = isTRUE(getOption("emc2.rlf_horizon_split", TRUE)),
    # Pair each solve with one at 1.25x the resolution and extrapolate; see
    # rlf_cache_solve in src/model_RLF.h.
    richardson = isTRUE(getOption("emc2.rlf_richardson", TRUE)),
    richardson_ratio = getOption("emc2.rlf_richardson_ratio", 1.25)
  )
}

.rlf_cols <- function(pars) {
  list(
    v = pars[, "v"],
    B = pars[, "B"],
    A = pars[, "A"],
    t0 = pars[, "t0"],
    s = pars[, "s"],
    alpha = pars[, "alpha"]
  )
}

.rlf_pdf_cdf <- function(rt, pars) {
  if (is.null(dim(pars))) {
    pars <- matrix(
      pars, nrow = length(rt), ncol = length(pars), byrow = TRUE,
      dimnames = list(NULL, names(pars))
    )
  } else if (nrow(pars) == 1L && length(rt) > 1L) {
    pars <- pars[rep(1L, length(rt)), , drop = FALSE]
  }
  if (nrow(pars) != length(rt)) {
    stop("pars must have one row per response time")
  }
  p <- .rlf_cols(pars)
  grid <- .rlf_grid()
  rlf_pdf_cdf_vec(
    rt, p$v, p$B, p$A, p$t0, p$s, p$alpha,
    as.integer(grid$nx), grid$dt_target, grid$tgrade,
    grid$adaptive, grid$explicit_inverse, grid$sparse_output,
    grid$simd_batch, grid$horizon_split, grid$richardson,
    grid$richardson_ratio
  )
}

dRLF <- function(rt, pars) .rlf_pdf_cdf(rt, pars)$pdf

pRLF <- function(rt, pars) .rlf_pdf_cdf(rt, pars)$cdf

rRLF <- function(lR, pars, ok = rep(TRUE, nrow(pars)),
                 dt = getOption("emc2.rlf_sim_dt", 1e-3),
                 t_max = getOption("emc2.rlf_sim_tmax", 30)) {
  required <- c("v", "B", "A", "t0", "s", "alpha")
  if (!all(required %in% colnames(pars))) {
    stop("pars must have columns ", paste(required, collapse = " "))
  }
  n_accumulators <- length(levels(lR))
  n_trials <- nrow(pars) / n_accumulators
  finish <- matrix(Inf, nrow = n_accumulators, ncol = n_trials)
  valid <- ok & is.finite(pars[, "v"])
  if (any(valid)) {
    finish[valid] <- rlf_hit_times_vec(
      pars[valid, "v"], pars[valid, "B"], pars[valid, "A"],
      pars[valid, "s"], pars[valid, "alpha"], dt, t_max
    )
  }

  all_infinite <- colSums(!is.infinite(finish)) == 0L
  winner <- max.col(-t(finish), ties.method = "first")
  pick <- cbind(winner, seq_len(ncol(finish)))
  rt <- matrix(pars[, "t0"], nrow = n_accumulators)[pick] + finish[pick]
  response <- factor(
    levels(lR)[winner], levels = levels(lR)
  )
  response[all_infinite] <- NA
  rt[all_infinite] <- Inf
  out <- data.frame(R = response, rt = rt)
  .apply_timed_guess_winner(out, levels(lR))
}

#' The Racing Lévy-Flight Model
#'
#' A race between independent evidence accumulators following
#' \deqn{dX_t = v\,dt + s\,dL_t^{(\alpha)},}
#' where the symmetric stable increment has characteristic function
#' \eqn{\exp[-s^\alpha dt |k|^\alpha / 2]}. Each accumulator starts uniformly
#' on \eqn{[0,A]} and finishes at \eqn{b=B+A}. At \eqn{\alpha=2} the process is
#' the racing diffusion model.
#'
#' The first-passage density is evaluated by a cached nonlocal
#' Fokker--Planck solver. Rows that share `(v/s, (B+A)/s, A/s, alpha)` share one
#' numerical march; `t0` only shifts response times and is not part of the
#' cache key.
#'
#' Numerical resolution is controlled with R options. `emc2.rlf_nx` is the
#' number of spatial intervals in the coarse grid (default 160).
#' `emc2.rlf_richardson` enables Richardson extrapolation by combining the
#' coarse grid with a fine grid whose size is set by
#' `emc2.rlf_richardson_ratio` (default 1.25, giving 160 and 200 intervals).
#' `emc2.rlf_dt` sets the target time step (default 0.016), and
#' `emc2.rlf_tgrade` controls time-grid grading (default 1). Keep these options
#' fixed throughout a fit. Lower values of \eqn{\alpha} generally require more
#' spatial resolution; increase `emc2.rlf_nx` when fitting especially
#' heavy-tailed data or when numerical accuracy in \eqn{\alpha} is critical.
#'
#' `emc2.rlf_horizon_split` (default `TRUE`) separates shorter and longer
#' response-time horizons so that a small number of long observations do not
#' coarsen the solver output for the rest of the data.
#'
#' Response times are never binned for this model. Because the likelihood is a
#' grid solve cached per parameter tuple, the cost is the time march out to
#' `max(rt)` and the individual response times only select readout points along
#' a march already paid for, so `rt_resolution` saves no computation while the
#' likelihood evaluates the density *at* the floored time and
#' \eqn{\hat{\alpha}} absorbs the mismatch. `RLF()` therefore declares
#' `compress_ok = FALSE` and `make_emc` forces `compress = FALSE` and
#' `rt_resolution = NULL`; nothing needs to be passed for this.
#'
#' Set `emc2.rlf_adaptive = TRUE` to enable the slower refinement checks used
#' by the standalone validation solver.
#'
#' @return A list defining an EMC2 race model.
#' @export
RLF <- function() {
  list(
    type = "RACE",
    c_name = "RLF",
    # The likelihood is a cached grid solve, so binning rt buys no speed (see
    # model_compress_ok) while costing accuracy in alpha.  make_emc() turns both
    # compression and rt_resolution off for this model.
    compress_ok = FALSE,
    p_types = c(
      v = log(1), B = log(1), A = log(0), t0 = log(0), s = log(1),
      alpha = qnorm(0.7), pContaminant = qnorm(0), pGuess = qnorm(0)
    ),
    p_types_canonical = c("v", "B", "A", "t0", "s", "alpha"),
    transform = list(
      func = c(
        v = "exp", B = "exp", A = "exp", t0 = "exp", s = "exp",
        alpha = "pnorm", pContaminant = "pnorm", pGuess = "pnorm"
      ),
      lower = c(alpha = 1),
      upper = c(alpha = 2)
    ),
    # Unlike the analytic race models, the parameters here set the solver's
    # discretisation as well as the process: the mesh is h = (b + extent)/nx, so
    # B near zero and alpha near its endpoints are not merely extreme, they are
    # numerically degenerate.  alpha is held off 2 because the heavy-tail term
    # that keeps the domain open scales as sin(pi*alpha/2) and vanishes there,
    # collapsing the domain onto b; B is floored for the same reason.  alpha = 2
    # is deliberately not an exception value: the pnorm transform returns
    # exactly 1 for sampled values past ~8.3, so the endpoint is reachable.
    bound = list(
      minmax = cbind(
        v = c(1e-3, Inf), B = c(1e-4, Inf), A = c(1e-4, Inf),
        t0 = c(0.05, Inf), s = c(0, Inf), alpha = c(1.01, 1.99),
        pContaminant = c(0.001, 0.999), pGuess = c(0.001, 0.999)
      ),
      exception = c(A = 0, pContaminant = 0, pGuess = 0)
    ),
    Ttransform = function(pars, dadm) {
      cbind(pars, b = pars[, "B"] + pars[, "A"])
    },
    rfun = function(data = NULL, pars) {
      rRLF(data$lR, pars, ok = attr(pars, "ok"))
    },
    dfun = function(rt, pars) dRLF(rt, pars),
    pfun = function(rt, pars) pRLF(rt, pars),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      log_likelihood_race_missing(
        pars = pars, dadm = dadm, model = model, min_ll = min_ll
      )
    }
  )
}
