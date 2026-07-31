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
    horizon_split = isTRUE(getOption("emc2.rlf_horizon_split", TRUE)),
    # Pair each solve with one at 1.5x the resolution and extrapolate; see
    # rlf_cache_solve in src/model_RLF.h.
    richardson = isTRUE(getOption("emc2.rlf_richardson", TRUE)),
    richardson_ratio = getOption("emc2.rlf_richardson_ratio", 1.5)
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
#' `emc2.rlf_nx` sets the *coarse* member of the extrapolated pair, so the
#' default of 128 solves at 128 and 192 and combines them; the figures below
#' are all for the pair, not for a single solve at that resolution.
#'
#' How much resolution a fit needs depends strongly on \eqn{\alpha}. Scored as
#' mean \eqn{|\Delta \log f|} per trial over the central 96% of each
#' distribution against a converged reference, the `emc2.rlf_nx` at which the
#' pair holds 0.01 is about 146 at \eqn{\alpha = 1.1}, 99 at 1.3, 70 at 1.5 and
#' 48 from 1.7 up, so the default is generous above \eqn{\alpha \approx 1.5}
#' and marginal below 1.2. Lowering `emc2.rlf_nx` to 96 for data expected to be
#' near-Brownian roughly halves the cost, at a mean displacement of
#' \eqn{\hat{\alpha}} of 0.014 against 0.008 for the default (2000 trials,
#' \eqn{\alpha} from 1.1 to 1.9) -- still several times more accurate than a
#' single unextrapolated solve at 128, which displaces \eqn{\hat{\alpha}} by
#' 0.091. Do not go below 96: the profile log-likelihood becomes visibly jagged
#' in \eqn{\alpha} (RMS roughness about a local quadratic at
#' \eqn{\alpha = 1.7} is 0.13 nats at `emc2.rlf_nx` 128 and 96, 0.60 at 64 and
#' 1.54 at 48), which a sampler sees as noise.
#'
#' Set `emc2.rlf_nx` once per fit rather than trying to vary it with
#' \eqn{\alpha}: the discretisation error is a bias, and a bias that moves with
#' \eqn{\alpha} biases \eqn{\hat{\alpha}} itself.
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
#' A low \eqn{\alpha} does not announce itself with implausible response times.
#' In a race the loser truncates the winner's tail, so with human-plausible
#' parameters (`t0` 0.1, `A` 0, `s` 1, a second accumulator at 0.6 of the
#' drift) the probability of a response beyond 3 s stays below 0.3% across the
#' whole of \eqn{\alpha \in [1.1, 1.9]}, and the 99th percentile of the pooled
#' response times moves by less than 0.12 s.
#'
#' \eqn{\alpha} is nonetheless well identified, because it acts on choice and
#' on the correct/error contrast rather than on the pooled tail. At
#' \eqn{v = 2, B = 1} accuracy rises from 0.657 at \eqn{\alpha = 1.9} to 0.735
#' at \eqn{\alpha = 1.1}, and the two response distributions pull apart as
#' \eqn{\alpha} falls: at 1.9 they are near-identical, as a Wald race with
#' equal thresholds requires (correct 0.267/0.434/0.778 against error
#' 0.265/0.435/0.779 at the 10th, 50th and 90th percentiles), while at 1.1 the
#' errors gain a long tail and the correct responses lose one (correct
#' 0.337/0.544/0.767 against error 0.240/0.535/0.907). Pooling over responses
#' hides almost all of this, so judge \eqn{\alpha} from accuracy and from the
#' error distribution, not from the marginal response times.
#'
#' The corresponding discriminability against the nearest Wald race -- the best
#' \eqn{\alpha = 2} fit in Kullback-Leibler divergence, optimising over drifts,
#' threshold and `t0` -- is 0.145 nats per trial at \eqn{\alpha = 1.1}, 0.083 at
#' 1.3 and 0.042 at 1.5, so 3 nats of evidence accumulate in roughly 20, 35 and
#' 70 trials. Detection designs with a single accumulator and an upper censor at
#' 3 s are weaker but still workable: at \eqn{v = 0.8, B = 1} the censoring rate
#' moves from 9.6% to 12.1% and the conditional median from 0.83 s to 1.08 s
#' over the same range, for 0.053 nats per trial at \eqn{\alpha = 1.1} (about 57
#' trials). Do not lower `emc2.rlf_nx` on the assumption that the heavy-tailed
#' end is unobservable; it is observable, and the discretisation error is not
#' small relative to the effect being measured.
#'
#' At a fixed `emc2.rlf_nx` the error is smooth in \eqn{\alpha} and acts as a
#' small bias rather than as noise. On a recovery check at \eqn{\alpha = 1.6}
#' with 2000 observations it displaces \eqn{\hat{\alpha}} by -0.020 at
#' `emc2.rlf_nx = 128`, against a posterior SD of about 0.035; raising
#' `emc2.rlf_nx` to 192 reduces the displacement to -0.003. Because the bias is
#' fixed while the posterior SD shrinks as \eqn{1/\sqrt{n}}, raise
#' `emc2.rlf_nx` when \eqn{\alpha} is of direct inferential interest and
#' \eqn{n} per participant is large.
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
    # The likelihood is a cached grid solve, so binning rt buys no speed (see
    # model_compress_ok) while costing accuracy in alpha.  make_emc() turns both
    # compression and rt_resolution off for this model.
    compress_ok = FALSE,
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
