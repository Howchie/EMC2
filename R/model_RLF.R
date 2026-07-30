# Racing symmetric alpha-stable Lévy-flight model.

.rlf_grid <- function() {
  list(
    nx = getOption("emc2.rlf_nx", 128L),
    dt_target = getOption("emc2.rlf_dt", 8e-3),
    tgrade = getOption("emc2.rlf_tgrade", getOption("emc2.fpe_tgrade", 1)),
    adaptive = isTRUE(getOption("emc2.rlf_adaptive", FALSE)),
    explicit_inverse = isTRUE(
      getOption("emc2.rlf_explicit_inverse", TRUE)
    ),
    sparse_output = isTRUE(getOption("emc2.rlf_sparse_output", TRUE)),
    # Off by default: see the note on Grid::simd_batch in src/model_RLF.h.
    simd_batch = isTRUE(getOption("emc2.rlf_simd_batch", FALSE)),
    horizon_split = isTRUE(getOption("emc2.rlf_horizon_split", TRUE))
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
    grid$simd_batch
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

  all_infinite <- apply(finish, 2L, function(x) all(is.infinite(x)))
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
#' Fokker--Planck solver that is second order in both space and time. Rows that
#' share `(v/s, (B+A)/s, A/s, alpha)` share one numerical march; `t0` only
#' shifts response times and is not part of the cache key. Likelihood
#' resolution is controlled by `emc2.rlf_nx`, `emc2.rlf_dt`, and
#' `emc2.rlf_tgrade`. The defaults use a single calibrated pass. Away from the
#' Brownian limit the spatial truncation dominates: below \eqn{\alpha = 2} the
#' error is unchanged across `emc2.rlf_dt` from 4e-3 to 8e-3, so accuracy there
#' is bought with `emc2.rlf_nx` (which converges at \eqn{O(h^2)}) rather than
#' with a smaller step. Strongly heavy-tailed fits (\eqn{\alpha} below about
#' 1.4) resolve the absorbing boundary least well and benefit most.
#'
#' The discretisation error is smooth in \eqn{\alpha} (residual roughness of
#' the profile log-likelihood is ~0.004 nats, far below anything a sampler
#' responds to), so it acts as a small bias rather than as noise. On a
#' recovery check at \eqn{\alpha = 1.6} with 2000 observations it displaces
#' \eqn{\hat{\alpha}} by -0.020 at `emc2.rlf_nx = 128`, against a posterior SD
#' of about 0.035; raising `emc2.rlf_nx` to 192 reduces the displacement to
#' -0.003. Because the bias is fixed while the posterior SD shrinks as
#' \eqn{1/\sqrt{n}}, raise `emc2.rlf_nx` when \eqn{\alpha} is of direct
#' inferential interest and \eqn{n} per participant is large.
#'
#' Set `emc2.rlf_adaptive = TRUE` for the slower refinement checks used by the
#' standalone validation solver.
#'
#' @return A list defining an EMC2 race model.
#' @export
RLF <- function() {
  list(
    type = "RACE",
    c_name = "RLF",
    p_types = c(
      v = log(1), B = log(1), A = log(0), t0 = log(0), s = log(1),
      alpha = qnorm(0.7), pContaminant = qnorm(0)
    ),
    p_types_canonical = c("v", "B", "A", "t0", "s", "alpha"),
    transform = list(
      func = c(
        v = "exp", B = "exp", A = "exp", t0 = "exp", s = "exp",
        alpha = "pnorm", pContaminant = "pnorm"
      ),
      lower = c(alpha = 1),
      upper = c(alpha = 2)
    ),
    bound = list(
      minmax = cbind(
        v = c(1e-3, Inf), B = c(0, Inf), A = c(1e-4, Inf),
        t0 = c(0.05, Inf), s = c(0, Inf), alpha = c(1, 2),
        pContaminant = c(0.001, 0.999)
      ),
      exception = c(A = 0, alpha = 2, pContaminant = 0)
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
