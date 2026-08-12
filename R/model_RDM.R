dRDM <- function(rt, pars)
# density for single accumulator
{
  out <- numeric(length(rt))
  ok <- rt > pars[, "t0"] & !pars[, "v"] < 0 # code handles rate zero case
  ok[is.na(ok)] <- FALSE
  if (any(dimnames(pars)[[2]] == "s")) { # rescale
    pars[ok, c("A", "B", "v")] <- pars[ok, c("A", "B", "v")] / pars[ok, "s"]
  }
  out[ok] <- dWald(rt[ok], v = pars[ok, "v"], B = pars[ok, "B"], A = pars[ok, "A"], t0 = pars[ok, "t0"])
  out
}


pRDM <- function(rt, pars)
# cumulative density for single accumulator
{
  out <- numeric(length(rt))
  ok <- rt > pars[, "t0"] & !pars[, "v"] < 0 # code handles rate zero case
  ok[is.na(ok)] <- FALSE
  if (any(dimnames(pars)[[2]] == "s")) { # rescale
    pars[ok, c("A", "B", "v")] <- pars[ok, c("A", "B", "v")] / pars[ok, "s"]
  }
  out[ok] <- pWald(rt[ok], v = pars[ok, "v"], B = pars[ok, "B"], A = pars[ok, "A"], t0 = pars[ok, "t0"])
  out
}

#### random

rWald <- function(n, B, v, A, s = 1, posdrift = TRUE)
# random function for single accumulator
{
  rwaldt <- function(n, k, l, s = 1, tiny = 1e-6) {
    # random sample of n from a Wald (or Inverse Gaussian)
    # k = criterion, l = rate, s = diffusion SD
    # about same speed as statmod rinvgauss

    rlevy <- function(n = 1, m = 0, c = 1) {
      if (any(c < 0)) stop("c must be positive")
      c / qnorm(1 - runif(n) / 2)^2 + m
    }

    s <- rep(s, length.out = n)
    flag <- l > abs(tiny)
    x <- rep(NA, times = n)

    x[!flag] <- rlevy(sum(!flag), 0, (k[!flag] / s[!flag])^2)
    mu <- k / l
    lambda <- (k / s)^2

    y <- rnorm(sum(flag))^2
    mu.0 <- mu[flag]
    lambda.0 <- lambda[flag]

    x.0 <- mu.0 + mu.0^2 * y / (2 * lambda.0) -
      sqrt(4 * mu.0 * lambda.0 * y + mu.0^2 * y^2) * mu.0 / (2 * lambda.0)

    z <- runif(length(x.0))
    test <- mu.0 / (mu.0 + x.0)
    x.0[z > test] <- mu.0[z > test]^2 / x.0[z > test]
    x[flag] <- x.0
    x[x < 0] <- max(x)
    x
  }
  
  out <- rep(Inf, n)
  s <- rep(s, length.out = n)
  neg <- rep(FALSE, n)
  if (posdrift) {
    pos <- v > 0
  } else {
    pos <- v >= 0
    neg <- v < 0
  }

  # positive (or zero) drift: standard inverse Gaussian
  # With posdrift=TRUE zero is excluded but for posdrift=FALSE it counts because mathematically it is an eventual guaranteed hit.
  npos <- sum(pos)
  if (npos > 0) {
    bs <- B[pos] + runif(npos, 0, A[pos])
    out[pos] <- rwaldt(npos, k = bs, l = v[pos], s = s[pos])
  }

  # negative drift with posdrift=FALSE: defective Wald via Bernoulli(p_hit)
  # Conditional FPT given hitting equals Wald with |v| (Girsanov / time-reversal)
  if (any(neg)) {
    nneg <- sum(neg)
    if (!posdrift) { # sample bernoulli hitting probability and use the absolute value of v for the finite finishes
      bs_neg <- B[neg] + runif(nneg, 0, A[neg])
      p_hit <- exp(2 * v[neg] * bs_neg / (s[neg]^2))  # v < 0, bs > 0 -> 0 < p_hit < 1
      hit <- as.logical(rbinom(nneg, 1, p_hit))
      if (any(hit)) {
        out[which(neg)[hit]] <- rwaldt(sum(hit), k = bs_neg[hit], l = abs(v[neg][hit]), s = s[neg][hit])
      }
    }
  }

  out
}

rRDM <- function(lR, pars, p_types=c("v", "B", "A", "t0"), ok=rep(TRUE, dim(pars)[1]))
                 # lR is an empty latent response factor lR with one level for each accumulator.
                 # pars is a matrix of corresponding parameter values named as in p_types
                 # pars must be sorted so accumulators and parameter for each trial are in
                 # contiguous rows. "s" parameter will be used but can be ommitted
                 #
                 # test
# pars=cbind(B=c(1,2),v=c(1,1),A=c(0,0),t0=c(.2,.2)); lR=factor(c(1,2))
{
  if (!all(p_types %in% dimnames(pars)[[2]])) {
    stop("pars must have columns ", paste(p_types, collapse = " "))
  }
  if (any(dimnames(pars)[[2]] == "s")) { # rescale
    pars[, c("A", "B", "v")] <- pars[, c("A", "B", "v")] / pars[, "s"]
  }
  pars[, "B"][pars[, "B"] < 0] <- 0 # Protection for negatives
  pars[, "A"][pars[, "A"] < 0] <- 0
  bad <- rep(NA, length(lR) / length(levels(lR)))
  out <- data.frame(R = bad, rt = bad)
  nr <- length(levels(lR))
  dt <- matrix(Inf, nrow = nr, ncol = nrow(pars) / nr)
  t0 <- pars[, "t0"]
  pars <- pars[ok, ]
  dt[ok] <- rWald(sum(ok), B = pars[, "B"], v = pars[, "v"], A = pars[, "A"])
  bad_col <- apply(dt, 2, function(x) all(is.infinite(x)))
  R <- max.col(-t(dt), ties.method = "first")
  pick <- cbind(R, 1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
  rt <- matrix(t0, nrow = nr)[pick] + dt[pick]
  R <- factor(levels(lR)[R], levels = levels(lR))
  R[bad_col] <- NA
  rt[bad_col] <- Inf
  out$R <- levels(lR)[R]
  out$R <- factor(out$R, levels = levels(lR))
  out$rt <- rt
  out <- .apply_timed_guess_winner(out, levels(lR))
  out
}

#' The Racing Diffusion Model
#'
#' Model file to estimate the Racing Diffusion Model (RDM), also known as the Racing Wald Model.
#'
#' Model files are almost exclusively used in `design()`.
#'
#' @details
#'
#' Default values are used for all parameters that are not explicitly listed in the `formula`
#' argument of `design()`.They can also be accessed with `RDM()$p_types`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default**   | **Mapping**          | **Interpretation**                                                |
#' |-----------|-----------|---------------|-----------|------------------|---------------------------------------------------------------|
#' | *v*       | log       | \[0, Inf\]      | log(1)    |                  | Evidence-accumulation rate (drift rate)                        |
#' | *A*       | log       | \[0, Inf\]      | log(0)    |                  | Between-trial variation (range) in start point                 |
#' | *B*       | log       | \[0, Inf\]      | log(1)    | *b* = *B* + *A*      | Distance from *A* to *b* (response threshold)                  |
#' | *t0*      | log       | \[0, Inf\]      | log(0)    |                  | Non-decision time                                             |
#' | *s*       | log       | \[0, Inf\]      | log(1)    |                  | Within-trial standard deviation of drift rate                 |
#' | *pContaminant* | probit | \[0, 1\] | qnorm(0) | | Optional *omission* contaminant probability: mass at `rt = Inf` only, handled by the data pipeline |
#' | *pGuess* | probit | \[0, 1\] | qnorm(0) | | Optional uniform *guess* (outlier) probability, mixed into observed RT densities over the guess window |
#'
#'
#' The core RDM parameters are estimated on the log scale. `pContaminant` is
#' estimated on the probit scale and is generic nuisance infrastructure rather
#' than an accumulator parameter.
#'
#' The parameterization *b* = *B* + *A* ensures that the response threshold is
#' always higher than the between trial variation in start point.
#'
#' Conventionally, `s` is fixed to 1 to satisfy scaling constraints.
#'
#' Because the RDM is a race model, it has one accumulator per response option.
#' EMC2 automatically constructs a factor representing the accumulators `lR` (i.e., the
#' latent response) with level names taken from the `R` column in the data.
#'
#' The `lR` factor is mainly used to allow for response bias, analogous to *Z* in the
#' DDM. For example, in the RDM, response thresholds are determined by the *B*
#' parameters, so `B~lR` allows for different thresholds for the accumulator
#' corresponding to "left" and "right" stimuli, for example, (e.g., a bias to respond left occurs
#' if the left threshold is less than the right threshold).
#'
#' For race models in general, the argument `matchfun` can be provided in `design()`.
#' One needs to supply a function that takes the `lR` factor (defined in the augmented data (d)
#' in the following function) and returns a logical defining the correct
#' response. In the example below, this is simply whether the `S` factor equals the
#' latent response factor: `matchfun=function(d)d$S==d$lR`. Using `matchfun` a latent match factor (`lM`) with
#' levels `FALSE` (i.e., the stimulus does not match the accumulator) and `TRUE`
#' (i.e., the stimulus does match the accumulator). This is added internally
#' and can also be used in model formula, typically for parameters related to
#' the rate of accumulation.
#'
#' Tillman, G., Van Zandt, T., & Logan, G. D. (2020). Sequential sampling models
#' without random between-trial variability: The racing diffusion model of speeded
#' decision making. *Psychonomic Bulletin & Review, 27*(5), 911-936.
#' https://doi.org/10.3758/s13423-020-01719-6
#'
#' @return A list defining the cognitive model
#' @examples

#' # When working with lM it is useful to design  an "average and difference"
#' # contrast matrix, which for binary responses has a simple canonical from:

#' ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
#' # We also define a match function for lM
#' matchfun=function(d)d$S==d$lR
#' # We now construct our design, with v ~ lM and the contrast for lM the ADmat.
#' design_RDMBE <- design(data = forstmann,model=RDM,matchfun=matchfun,
#'                        formula=list(v~lM,s~lM,B~E+lR,A~1,t0~1),
#'                        contrasts=list(v=list(lM=ADmat)),constants=c(s=log(1)))
#' # For all parameters that are not defined in the formula, default values are assumed
#' # (see Table above).
#' @export

RDM <- function() {
  list(
    type = "RACE",
    c_name = "RDM",
    p_types = c("v" = log(1), "B" = log(1), "A" = log(0), "t0" = log(0), "s" = log(1), "pContaminant" = qnorm(0), "pGuess" = qnorm(0)),
    p_types_canonical = c("v", "B", "A", "t0", "s"),
    transform = list(func = c(v = "exp", B = "exp", A = "exp", t0 = "exp", s = "exp", pContaminant = "pnorm", pGuess = "pnorm")),
    bound = list(
      minmax = cbind(v = c(1e-3, Inf), B = c(0, Inf), A = c(1e-4, Inf), t0 = c(0.05, Inf), s = c(0, Inf), pContaminant = c(0.001, 0.999), pGuess = c(0.001, 0.999)),
      exception = c(A = 0, v = 0, pContaminant = 0, pGuess = 0)
    ),
    # Trial dependent parameter transform
    Ttransform = function(pars, dadm) {
      pars <- cbind(pars, b = pars[, "B"] + pars[, "A"])
      pars
    },
    # Random function for racing accumulators
    rfun = function(data = NULL, pars) .rfun_RDM(data$lR, pars, ok = attr(pars, "ok")),
    # Density function (PDF) for single accumulator
    dfun = function(rt, pars) dRDM(rt, pars),
    # Probability function (CDF) for single accumulator
    pfun = function(rt, pars) pRDM(rt, pars),
    # Race likelihood combining pfun and dfun
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      log_likelihood_race_missing(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}

# ============================================================================
# RDMGBM: Racing Geometric Brownian Motion with start-point variability
# ============================================================================

dRDMGBM <- function(rt, pars, erlang = 1L) {
  if (is.null(dim(pars)) || (dim(pars)[1] == 1 & length(rt) > 1)) {
    original_names <- names(pars)
    if (is.null(original_names)) original_names <- colnames(pars)
    pars <- matrix(pars,
      nrow = length(rt), ncol = length(pars),
      dimnames = list(NULL, original_names), byrow = TRUE
    )
  }
  if (!("b" %in% colnames(pars)) && all(c("B", "A") %in% colnames(pars))) {
    pars <- cbind(pars, b = 1 + pars[, "B"] + pars[, "A"])
  }
  if (!all(c("lambda_g", "lambda_k") %in% colnames(pars))) {
    stop("RDMGBM requires parameter columns 'lambda_g' and 'lambda_k'.")
  }
  out <- rep(NaN, length(rt))
  erl <- (pars[, "lambda_g", drop = FALSE] > 0) | (pars[, "lambda_k", drop = FALSE] > 0)
  ok <- (rt > 0) & ((rt > pars[, "t0", drop = FALSE]) | erl) & !(pars[, "v", drop = FALSE] < 0)
  ok[is.na(ok)] <- FALSE
  if (any(ok)) {
    out[ok] <- dGBMspv(rt[ok],
      v = pars[ok, "v", drop = FALSE], b = pars[ok, "b", drop = FALSE],
      A = pars[ok, "A", drop = FALSE], t0 = pars[ok, "t0", drop = FALSE],
      s = pars[ok, "s", drop = FALSE], lambda_g = pars[ok, "lambda_g", drop = FALSE],
      lambda_k = pars[ok, "lambda_k", drop = FALSE],
      kill_shape = erlang, erlang_omega = .rdmswtn_erlang_omega(pars[ok, , drop = FALSE], erlang)
    )
  }
  out
}

pRDMGBM <- function(rt, pars, erlang = 1L) {
  if (is.null(dim(pars)) || (dim(pars)[1] == 1 & length(rt) > 1)) {
    original_names <- names(pars)
    if (is.null(original_names)) original_names <- colnames(pars)
    pars <- matrix(pars,
      nrow = length(rt), ncol = length(pars),
      dimnames = list(NULL, original_names), byrow = TRUE
    )
  }
  if (!("b" %in% colnames(pars)) && all(c("B", "A") %in% colnames(pars))) {
    pars <- cbind(pars, b = 1 + pars[, "B"] + pars[, "A"])
  }
  if (!all(c("lambda_g", "lambda_k") %in% colnames(pars))) {
    stop("RDMGBM requires parameter columns 'lambda_g' and 'lambda_k'.")
  }
  out <- rep(NaN, length(rt))
  erl <- (pars[, "lambda_g", drop = FALSE] > 0) | (pars[, "lambda_k", drop = FALSE] > 0)
  ok <- (rt > 0) & ((rt > pars[, "t0", drop = FALSE]) | erl) & !(pars[, "v", drop = FALSE] < 0)
  ok[is.na(ok)] <- FALSE
  if (any(ok)) {
    out[ok] <- pGBMspv(rt[ok],
      v = pars[ok, "v", drop = FALSE], b = pars[ok, "b", drop = FALSE],
      A = pars[ok, "A", drop = FALSE], t0 = pars[ok, "t0", drop = FALSE],
      s = pars[ok, "s", drop = FALSE], lambda_g = pars[ok, "lambda_g", drop = FALSE],
      lambda_k = pars[ok, "lambda_k", drop = FALSE],
      kill_shape = erlang, erlang_omega = .rdmswtn_erlang_omega(pars[ok, , drop = FALSE], erlang)
    )
  }
  out
}

rGBM <- function(n, b, v, A, s = 1) {
  out <- rep(Inf, n)
  if (n <= 0) {
    return(out)
  }
  if (n > 1 && all(length(b) == 1, length(v) == 1, length(A) == 1, length(s) == 1)) {
    b <- rep(b, n)
    v <- rep(v, n)
    A <- rep(A, n)
    s <- rep(s, n)
  }
  A[A < 0] <- 0
  x0 <- 1 + runif(n, 0, A)
  d <- log(b / x0)
  mu_log <- v - 0.5 * s^2
  ok <- is.finite(mu_log) & is.finite(d) & (d > 0)
  if (any(ok)) {
    out[ok] <- statmod::rinvgauss(sum(ok), mean = d[ok] / mu_log[ok], shape = d[ok]^2 / s[ok]^2)
  }
  out
}

rRDMGBM <- function(lR, pars, p_types = c("v", "b", "A", "t0", "s", "lambda_g", "lambda_k"),
                    ok = rep(TRUE, dim(pars)[1]), erlang_shape = 1L, erlang_type = "none") {
  if (!is.null(attr(pars, "ok"))) ok <- attr(pars, "ok")
  if (!("b" %in% dimnames(pars)[[2]]) && all(c("B", "A") %in% dimnames(pars)[[2]])) {
    pars <- cbind(pars, b = 1 + pars[, "B"] + pars[, "A"])
  }
  if (!all(c("lambda_g", "lambda_k") %in% colnames(pars))) {
    stop("RDMGBM requires parameter columns 'lambda_g' and 'lambda_k'.")
  }
  erlang_omega_all <- .rdmswtn_erlang_omega(pars, erlang_shape)
  required <- c("v", "b", "A", "t0", "s", "lambda_g", "lambda_k")
  if (!all(required %in% dimnames(pars)[[2]])) {
    stop("pars must have columns ", paste(required, collapse = " "))
  }
  pars[, "b"][pars[, "b"] < 1 + 1e-8] <- 1 + 1e-8
  pars[, "A"][pars[, "A"] < 0] <- 0
  bad <- rep(NA, length(lR) / length(levels(lR)))
  out <- data.frame(R = bad, rt = bad)
  nr <- length(levels(lR))
  n_trials <- nrow(pars) / nr
  if (length(ok) != nrow(pars)) stop("ok must have length nrow(pars).")
  dt <- matrix(Inf, nrow = nr, ncol = n_trials)
  t0 <- pars[, "t0"]

  guess  <- erlang_type == "local_guess" || erlang_type == "local_kill_guess"
  global <- erlang_type == "global_kill"
  local_kill <- erlang_type == "local_kill" || erlang_type == "local_kill_guess"

  # global_kill: draw one shared Erlang timer per trial
  if (global) {
    lambda_mat <- matrix(pars[, "lambda_k"], nrow = nr)
    if (any(apply(lambda_mat, 2, function(x) length(unique(x)) > 1)))
      stop("global_kill requires lambda_k to be constant across accumulators")
    lambda_trials <- lambda_mat[1, ]
    tk <- rep(Inf, n_trials)
    kill_ok <- !is.na(lambda_trials) & lambda_trials > 0
    if (any(kill_ok)) {
      shape_k <- if (as.integer(erlang_shape) == 3L) {
        ifelse(runif(sum(kill_ok)) <= matrix(erlang_omega_all, nrow = nr)[1, kill_ok], 1L, 2L)
      } else {
        as.integer(erlang_shape)
      }
      rate_k <- if (as.integer(erlang_shape) == 3L) {
        ifelse(shape_k == 2L, 2 * lambda_trials[kill_ok], lambda_trials[kill_ok])
      } else {
        lambda_trials[kill_ok]
      }
      tk[kill_ok] <- rgamma(sum(kill_ok), shape = shape_k, rate = rate_k)
    }
  }

  pars_ok <- pars[ok, , drop = FALSE]
  if (nrow(pars_ok) > 0) {
    # Local Erlang clocks are handled on the raw-time axis below, after t0 is
    # added to the evidence-accumulation finish times.
    k_vec <- rep(0, nrow(pars_ok))
    dt[ok] <- rGBM_killed(sum(ok),
      b = pars_ok[, "b"], v = pars_ok[, "v"], A = pars_ok[, "A"],
      s = pars_ok[, "s"], k = k_vec, erlang = erlang_shape,
      erlang_omega = erlang_omega_all[ok]
    )
  }
  # Put EAM on the same raw-time axis as Erlang clocks.
  dt <- dt + matrix(t0, nrow = nr)

  if (guess || local_kill) {
    tg_local <- matrix(Inf, nrow = nr, ncol = n_trials)
    tk_local <- matrix(Inf, nrow = nr, ncol = n_trials)

    if (guess) {
      lambda_g_local <- matrix(0, nrow = nr, ncol = n_trials)
      active_guess <- matrix(ok, nrow = nr) & (levels(lR) != "nogo")
      lambda_g_local[active_guess] <- pars[, "lambda_g"][active_guess]
      guess_ok_mat <- lambda_g_local > 0
      if (any(guess_ok_mat)) {
        shape_g <- if (as.integer(erlang_shape) == 3L) {
          ifelse(runif(sum(guess_ok_mat)) <= matrix(erlang_omega_all, nrow = nr)[guess_ok_mat], 1L, 2L)
        } else {
          as.integer(erlang_shape)
        }
        rate_g <- if (as.integer(erlang_shape) == 3L) {
          ifelse(shape_g == 2L, 2 * lambda_g_local[guess_ok_mat], lambda_g_local[guess_ok_mat])
        } else {
          lambda_g_local[guess_ok_mat]
        }
        tg_local[guess_ok_mat] <- rgamma(sum(guess_ok_mat), shape = shape_g,
                                         rate = rate_g)
      }
    }

    if (local_kill) {
      lambda_k_local <- matrix(0, nrow = nr, ncol = n_trials)
      lambda_k_local[matrix(ok, nrow = nr)] <- pars[ok, "lambda_k"]
      kill_ok_mat <- lambda_k_local > 0
      if (any(kill_ok_mat)) {
        shape_k <- if (as.integer(erlang_shape) == 3L) {
          ifelse(runif(sum(kill_ok_mat)) <= matrix(erlang_omega_all, nrow = nr)[kill_ok_mat], 1L, 2L)
        } else {
          as.integer(erlang_shape)
        }
        rate_k <- if (as.integer(erlang_shape) == 3L) {
          ifelse(shape_k == 2L, 2 * lambda_k_local[kill_ok_mat], lambda_k_local[kill_ok_mat])
        } else {
          lambda_k_local[kill_ok_mat]
        }
        tk_local[kill_ok_mat] <- rgamma(sum(kill_ok_mat), shape = shape_k,
                                        rate = rate_k)
      }
    }

    if (guess && local_kill) {
      guess_win = tg_local<dt & tg_local<tk_local
      dt_candidate <- pmin(dt, tg_local)
      dt <- ifelse(dt_candidate < tk_local, dt_candidate, Inf)
    } else if (guess) {
      guess_win = tg_local<dt
      dt <- pmin(dt, tg_local)
    } else {
      dt <- ifelse(dt < tk_local, dt, Inf)
    }
  }

  if (global) {
    # Global kill: shared timer fires → no response this trial
    is_killed <- tk < apply(dt, 2, min)
    if (any(is_killed)) {
      dt[, is_killed] <- Inf
    }
  }

  bad_col <- apply(dt, 2, function(x) all(is.infinite(x)))
  R <- apply(dt, 2, which.min)
  pick <- cbind(R, 1:dim(dt)[2])
  rt <- dt[pick]
  out$R <- levels(lR)[R]
  out$R <- factor(out$R, levels = levels(lR))
  out$rt <- rt
  out$R[bad_col] <- NA
  out$rt[bad_col] <- Inf
  out <- .apply_timed_guess_winner(out, levels(lR))
  if (guess) {
    out$isTime <- rep(NA, n_trials)
    out$isTime[!bad_col] <- guess_win[pick][!bad_col]
  }
  out
}

rGBM_killed <- function(n, b, v, A, s = 1, k = 0, erlang = 1L, erlang_omega = 1) {
  out <- rGBM(n, b = b, v = v, A = A, s = s)
  kill_idx <- which(k > 0)
  if (length(kill_idx) > 0) {
    shape <- if (as.integer(erlang) == 3L) {
      ifelse(runif(length(kill_idx)) <= erlang_omega[kill_idx], 1L, 2L)
    } else {
      as.integer(erlang)
    }
    rate <- if (as.integer(erlang) == 3L) {
      ifelse(shape == 2L, 2 * k[kill_idx], k[kill_idx])
    } else {
      k[kill_idx]
    }
    tk <- rgamma(length(kill_idx), shape = shape, rate = rate)
    out[kill_idx[tk <= out[kill_idx]]] <- Inf
  }
  out
}

# `drift_override` supplies a pre-drawn rate per element, which is how the
# correlated-draw models inject their equicorrelated draws; NULL keeps the
# ordinary independent per-element draw.
rSWTN <- function(n, b, v, A, sv, s = 1, k = 0, erlang = 1L, erlang_omega = 1,
                  posdrift = TRUE, drift_override = NULL) {
  if (n <= 0) return(numeric(0))
  b <- rep(b, length.out = n)
  v <- rep(v, length.out = n)
  A <- rep(A, length.out = n)
  sv <- rep(sv, length.out = n)
  s <- rep(s, length.out = n)
  k <- rep(k, length.out = n)
  erlang_omega <- rep(erlang_omega, length.out = n)
  out <- rep(Inf, n)
  # For sv > 0 and posdrift=TRUE, draw per-trial drifts from N(v, sv^2)
  # truncated at zero. Otherwise draw from the full normal and let rWald
  # handle defective negative-drift finite hits when posdrift=FALSE.
  # For sv == 0 and v < 0 with posdrift=FALSE: Bernoulli(p_hit) sampling in rWald.
  v_draw <- v
  if (!is.null(drift_override)) {
    v_draw <- rep(drift_override, length.out = n)
  } else {
    sample <- is.finite(sv) & sv > 1e-12
    if (any(sample)) {
      if (posdrift) {
        lo <- pnorm(0, mean = v[sample], sd = sv[sample])
        u <- lo + runif(sum(sample)) * (1 - lo)
        v_draw[sample] <- qnorm(u, mean = v[sample], sd = sv[sample])
      } else {
        v_draw[sample] <- rnorm(sum(sample), mean = v[sample], sd = sv[sample])
      }
    }
  }
  
  out <- rWald(n, B = b - A, v = v_draw, A = A, s = s, posdrift = posdrift)
  kill_idx <- which(k > 0)
  if (length(kill_idx) > 0) {
    shape <- if (as.integer(erlang) == 3L) {
      ifelse(runif(length(kill_idx)) <= erlang_omega[kill_idx], 1L, 2L)
    } else {
      as.integer(erlang)
    }
    rate <- if (as.integer(erlang) == 3L) {
      ifelse(shape == 2L, 2 * k[kill_idx], k[kill_idx])
    } else {
      k[kill_idx]
    }
    tk <- rgamma(length(kill_idx), shape = shape, rate = rate)
    out[kill_idx[tk <= out[kill_idx]]] <- Inf
  }
  out
}

#' RDMGBM Model
#'
#' Racing geometric Brownian motion (GBM) first-passage model with
#' start-point variability. Each accumulator starts at
#' `X(0) = 1 + U`, where `U ~ Uniform(0, A)`, and races to the boundary
#' `b = 1 + B + A`. On the log scale, the GBM has drift
#' `v - s^2 / 2` and diffusion SD `s`; `v` is the drift parameter on the
#' original GBM scale. The first-passage time is shifted by `t0`, and the
#' accumulator with the earliest finish wins the race.
#'
#' The model's user-facing parameters are:
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *v* | log | \[0, Inf\] | log(1) | | GBM drift parameter. |
#' | *B* | log | \[0, Inf\] | log(1) | *b* = 1 + *B* + *A* | Distance from the upper start-point range to the boundary, with the GBM baseline at 1. |
#' | *A* | log | \[0, Inf\] | log(0) | | Start-point range above 1. |
#' | *t0* | log | \[0, Inf\] | log(0) | | Non-decision time. |
#' | *s* | log | \[0, Inf\] | log(1) | | GBM diffusion SD. |
#' | *mG* | log | \[0, Inf\] | log(1) | *lambda_g* = *q* / *mG* | Mean of the optional guess clock. |
#' | *mK* | log | \[0, Inf\] | log(1) | *lambda_k* = *q* / *mK* | Mean of the optional kill clock. |
#' | *omega* | probit | \[0, 1\] | qnorm(.5) | | Erlang-1 mixture weight in mixed mode. |
#' | *pContaminant* | probit | \[0, 1\] | qnorm(0) | | Optional *omission* contaminant probability: mass at `rt = Inf` only, handled by the data pipeline |
#' | *pGuess* | probit | \[0, 1\] | qnorm(0) | | Optional uniform *guess* (outlier) probability, mixed into observed RT densities over the guess window |
#'
#' The internal `lambda_g` and `lambda_k` columns are created by the
#' `Ttransform`; they are rates, not parameters to include in a design
#' formula. Here `q = 1` for Erlang-1 and `q = 2` for Erlang-2. In mixed mode,
#' each clock is Erlang-1 with probability `omega` and Erlang-2 otherwise,
#' with the same mean `mG` or `mK` in either component.
#'
#' `erlang_type = "none"` gives the ordinary GBM race. `"local_kill"` adds an
#' independent kill clock to each accumulator, `"global_kill"` adds one kill
#' clock shared within a trial, `"local_guess"` adds local guess clocks, and
#' `"local_kill_guess"` adds both local clock types. A kill winner is an
#' omission; a guess winner follows the package's timed-guess convention. A
#' global kill requires the kill parameter to be constant across accumulators
#' within a trial. Mixed Erlang clocks are currently supported only for local
#' clock configurations, not `"global_kill"`.
#'
#' For finite-time simulation the model checks the GBM log-drift condition
#' `v > s^2 / 2` together with `s > 0`. The `v = 0` exception is retained for
#' parameter mapping, but it does not produce ordinary finite GBM first
#' passages. The model is a race model, so EMC2 constructs the accumulator
#' factor `lR` from the response levels in `R`.
#'
#' @param erlang_shape Integer `1` for exponential clocks, `2` for Erlang-2,
#'   or `"mixed"` for the Erlang-1/Erlang-2 mixture.
#' @param erlang_type Clock configuration: one of `"none"`, `"local_kill"`,
#'   `"global_kill"`, `"local_guess"`, or `"local_kill_guess"`.
#' @return A model list compatible with [design()].
#'
#' @export
#'
RDMGBM <- function(erlang_shape = 1L, erlang_type = "none") {
  erlang_type <- match.arg(erlang_type, c("none", "local_kill", "global_kill", "local_guess", "local_kill_guess"))
  erlang_mixed <- identical(erlang_shape, "mixed")
  if (erlang_mixed && erlang_type == "global_kill") {
    stop("RDMGBM erlang_shape = 'mixed' is currently implemented for local Erlang processes only.")
  }
  erlang_shape_cpp <- if (erlang_mixed) 3L else as.integer(erlang_shape)
  type_suffix <- switch(erlang_type,
    "none" = "",
    "local_kill" = "_LOCAL_KILL",
    "global_kill" = "_GLOBAL_KILL",
    "local_guess" = "_LOCAL_GUESS",
    "local_kill_guess" = "_LOCAL_KILL_GUESS"
  )
  base_name <- if (erlang_mixed) "RDMGBM_EMIX" else if (erlang_shape_cpp >= 2L) "RDMGBM_E2" else "RDMGBM"
  
  has_guess <- erlang_type %in% c("local_guess", "local_kill_guess")
  has_kill  <- erlang_type %in% c("local_kill", "global_kill", "local_kill_guess")
  
  p_types <- c(
    "v" = log(1), "B" = log(1), "A" = log(0), "t0" = log(0), "s" = log(1)
  )
  transform <- c(
    v = "exp", B = "exp", A = "exp", t0 = "exp", s = "exp"
  )
  minmax <- cbind(
    v = c(1e-3, Inf), B = c(0, Inf), A = c(0, Inf),
    t0 = c(0.05, Inf), s = c(0, Inf)
  )
  exception <- c(A = 0, v = 0)

  p_types  <- c(p_types,  mG = log(1))
  transform <- c(transform, mG = "exp")
  minmax   <- cbind(minmax, mG = c(1e-4, Inf))
  exception <- c(exception, mG = 0)

  p_types  <- c(p_types,  mK = log(1))
  transform <- c(transform, mK = "exp")
  minmax   <- cbind(minmax, mK = c(1e-4, Inf))
  exception <- c(exception, mK = 0)

  if (erlang_mixed) {
    p_types  <- c(p_types,  omega = qnorm(0.5))
    transform <- c(transform, omega = "pnorm")
    minmax   <- cbind(minmax, omega = c(0, 1))
    exception <- c(exception, omega = 0)
  }

  # pContaminant (omission) and pGuess (uniform outlier); see add_nuisance_pars().
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- .nuis$p_types; transform <- .nuis$transform
  minmax <- .nuis$minmax; exception <- .nuis$exception

  list(
    type = "RACE",
    c_name = paste0(base_name, type_suffix),
    p_types = p_types,
    p_types_canonical = c("v", "B", "A", "t0", "s"),
    transform = list(func = transform),
    bound = list(
      minmax = minmax,
      exception = exception,
      # Joint validity: positive log-space drift required for finite-time simulation.
      joint_ok = function(pars) {
        if (!all(c("v", "s") %in% colnames(pars))) return(rep(TRUE, nrow(pars)))
        is.finite(pars[, "v"]) & is.finite(pars[, "s"]) &
          pars[, "s"] > 0 & (pars[, "v"] > 0.5 * pars[, "s"]^2)
      }
    ),
    Ttransform = function(pars, dadm) {
      lambda_factor <- if (erlang_shape_cpp == 2L) 2 else 1
      n <- nrow(pars)
      mG_val <- pars[, "mG"]
      mK_val <- pars[, "mK"]
      lg <- if (has_guess) ifelse(mG_val <= 0, 0, lambda_factor / mG_val) else rep(0, n)
      lk <- if (has_kill)  ifelse(mK_val <= 0, 0, lambda_factor / mK_val) else rep(0, n)
      timed <- cbind(lambda_g = lg, lambda_k = lk)
      extra_drop <- c("v", "B", "A", "t0", "s", "mG", "mK", "omega")
      extra <- pars[, setdiff(colnames(pars), extra_drop), drop = FALSE]
      if (erlang_mixed) {
        pars <- cbind(
          pars[, c("v", "B", "A", "t0", "s"), drop = FALSE],
          timed,
          omega = pars[, "omega"],
          extra
        )
      } else {
        pars <- cbind(
          pars[, c("v", "B", "A", "t0", "s"), drop = FALSE],
          timed,
          extra
        )
      }
      pars <- cbind(pars, b = 1 + pars[, "B"] + pars[, "A"])
      pars
    },
    rfun = function(data = NULL, pars) {
      ok <- attr(pars, "ok")
      if (is.null(ok)) ok <- rep(TRUE, nrow(pars))
      rRDMGBM(data$lR, pars, ok = ok, erlang_shape = erlang_shape_cpp, erlang_type = erlang_type)
    },
    dfun = function(rt, pars) dRDMGBM(rt, pars, erlang = erlang_shape_cpp),
    pfun = function(rt, pars) pRDMGBM(rt, pars, erlang = erlang_shape_cpp),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      log_likelihood_race_missing(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}

# ============================================================================
# RDMSWTN: Racing Diffusion Model with Shifted Wald Truncated Normal
# Superset of RDM: sv=0,A=0 reduces to point Wald; sv=0 reduces to standard RDM
# ============================================================================

#' RDMSWTN Model
#'
#' Racing Diffusion Model with Shifted Wald Truncated Normal (SWTN)
#' accumulators. Conditional on a trial-specific drift draw `V`, each
#' accumulator has a Wald first-passage time with diffusion SD `s`, a
#' non-decision shift `t0`, and a uniformly varying first-passage distance.
#' The start-point/threshold parameterization is `b = B + A`, with the
#' corresponding Wald distance written as `B + U * A`, where
#' `U ~ Uniform(0, 1)`. The between-trial drift is
#' `V ~ Normal(v, sv^2)`; when `posdrift = TRUE`, this normal distribution is
#' truncated below at zero.
#'
#' The model parameters are:
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *v* | log | \[0, Inf\] | log(1) | | Mean between-trial drift rate. |
#' | *B* | log | \[0, Inf\] | log(1) | *b* = *B* + *A* | Baseline distance to the threshold. |
#' | *A* | log | \[0, Inf\] | log(0) | | Range of between-trial start-point/distance variability. |
#' | *t0* | log | \[0, Inf\] | log(0) | | Non-decision time. |
#' | *s* | log | \[0, Inf\] | log(1) | | Within-trial diffusion SD; conventionally fixed to 1 for scale identification. |
#' | *sv* | log | \[0, Inf\] | log(0) | | Between-trial SD of the drift rate. |
#' | *mG* | log | \[0, Inf\] | log(1) | *lambda_g* = *q* / *mG* | Mean of the optional guess clock. |
#' | *mK* | log | \[0, Inf\] | log(1) | *lambda_k* = *q* / *mK* | Mean of the optional kill clock. |
#' | *omega* | probit | \[0, 1\] | qnorm(.5) | | Erlang-1 mixture weight in mixed mode. |
#' | *pContaminant* | probit | \[0, 1\] | qnorm(0) | | Optional *omission* contaminant probability: mass at `rt = Inf` only, handled by the data pipeline |
#' | *pGuess* | probit | \[0, 1\] | qnorm(0) | | Optional uniform *guess* (outlier) probability, mixed into observed RT densities over the guess window |
#' | *rho* | scaled probit | \[-1, 1\] | qnorm(.5) | | Gaussian-copula correlation between the two participating finishing times; only when `correlated = TRUE`. |
#'
#' `erlang_shape = 1` uses exponential clocks and `erlang_shape = 2` uses
#' Erlang-2 clocks. In `erlang_shape = "mixed"`, each clock is Erlang-1 with
#' probability `omega` and Erlang-2 otherwise, with both components retaining
#' the same mean `mG` or `mK`. The internal `lambda_g` and `lambda_k` columns
#' are generated from these means by the `Ttransform`; they are not parameters
#' for `design()` formulas. `erlang_type` selects no clocks, local kill, global
#' kill, local guess, or local kill plus local guess, using the same values
#' `"none"`, `"local_kill"`, `"global_kill"`, `"local_guess"`, and
#' `"local_kill_guess"`. A global kill clock is shared across accumulators and
#' therefore requires the kill parameter to be constant within a trial. Mixed
#' Erlang mode is currently restricted to local clock configurations.
#'
#' Setting `sv = 0` reduces the model to the standard RDM parameterization;
#' setting both `sv = 0` and `A = 0` gives a point Wald accumulator (apart from
#' any optional clock). With `posdrift = FALSE`, negative drifts are allowed and
#' the finish-time distribution can be defective: an accumulator can have an
#' intrinsic infinite finish time. The compiled model name is suffixed with
#' `IO` in that case, and `v` is sampled on the natural scale
#' (`transform = "identity"`, bounds `(-Inf, Inf)`, default `1`) so that
#' negative mean rates are reachable. With `posdrift = TRUE`, the drift
#' distribution is truncated at zero, `v` is sampled on the log scale
#' (bounds `(1e-3, Inf)`), and the ordinary RDM is recovered when `sv = 0`.
#'
#' As a race model, RDMSWTN has one accumulator per response option. EMC2
#' constructs the latent accumulator factor `lR` from `R`, and the race
#' likelihood combines one winning density with the survivor probabilities of
#' all other accumulators. The optional `pContaminant` and `pGuess` parameters
#' are generic nuisance infrastructure and are not part of the SWTN
#' distribution: `pContaminant` is an *omission* rate contributing mass only at
#' `rt = Inf`, `pGuess` a uniform *outlier* density mixed into observed RTs.
#' Both are proportions among *retained* trials, being applied after truncation
#' renormalisation.
#'
#' With `correlated = TRUE`, exactly two active accumulator rows in a trial may
#' have the same signed, nonzero natural-scale `rho`; every other row must have
#' `rho = 0`. A binary race can therefore use `rho ~ 1`. In a race with two
#' coupled accumulators and additional independent accumulators, use a
#' participation factor, for example `rho ~ 0 + coupled`, and fix the
#' opted-out coefficient to zero. Zero or one nonzero row is an ordinary
#' independent RDMSWTN race.
#'
#' `correlate` picks which of two genuinely different data generating
#' processes `rho` describes.
#'
#' With `correlate = "times"` (the default) correlation is applied to
#' finishing-time ranks through a Gaussian copula, including when `sv > 0`; it
#' is not a correlation of drift draws or the Pearson correlation of observed
#' response times. For the participating pair, Kendall's tau is
#' `2 * asin(rho) / pi` and Spearman's rho is `6 * asin(rho / 2) / pi`.
#'
#' This copula is available under `posdrift = FALSE` provided `sv = 0` on the
#' correlated rows. The marginals are then defective: a negative mean rate
#' gives an accumulator with intrinsic never-finish mass, `F(Inf) = p < 1`.
#' The copula construction is unchanged and remains exact, because copula
#' uniforms at or above the marginal plateau represent the atom at infinity.
#' `rho` therefore couples the omissions as well as the finite finishing
#' times: the probability that neither of a coupled pair ever finishes is
#' `Phi2(qnorm(p1), qnorm(p2); rho)` rather than the product `(1-p1)(1-p2)`.
#'
#' With `correlate = "drifts"` `rho` is instead the correlation of the
#' *between-trial drift draws themselves*, exactly as in [BAwLcorr()]: the
#' rates are drawn from an equicorrelated normal (truncated to the positive
#' orthant when `posdrift = TRUE`) and the accumulators then race
#' independently on the drawn rates. This requires `sv > 0` on the correlated
#' rows, since at `sv = 0` there are no draws to correlate. Both models share
#' one implementation of the construction and of the shared-factor quadrature
#' that integrates it out; see `src/drift_factor.h`. A `matchfun` is required
#' either way, and the correct racer is the positive reference while the
#' incorrect racers carry the cell sign, so a pair's loadings multiply back to
#' `rho`.
#'
#' The two settings are not nested and neither is a limit of the other, so
#' `correlate` must be chosen deliberately; `"times"` is the default because
#' it is what earlier versions did.
#'
#' Correlated RDMSWTN currently supports only `erlang_type = "none"`.
#' [LogicalRulesRDMSWTN()] remains an independent model; correlated
#' logical-rule races are not supported.
#'
#' @param erlang_shape Integer `1` for exponential clocks, `2` for Erlang-2,
#'   or `"mixed"` for the Erlang-1/Erlang-2 mixture.
#' @param erlang_type Clock configuration: one of `"none"`, `"local_kill"`,
#'   `"global_kill"`, `"local_guess"`, or `"local_kill_guess"`.
#' @param posdrift Logical. If `TRUE` (default), truncate the between-trial
#'   normal drift distribution below at zero; if `FALSE`, use untruncated
#'   drifts and allow intrinsic omissions.
#' @param correlated Logical. If `TRUE`, add the correlation parameter `rho`
#'   and use the correlated race likelihood and simulator.
#' @param correlate What `rho` correlates when `correlated = TRUE`: `"times"`
#'   (default) for a Gaussian copula on the two participating finishing times,
#'   or `"drifts"` for correlated between-trial drift draws in the style of
#'   [BAwLcorr()]. `"drifts"` requires `sv > 0` on the correlated rows.
#' @return A model list compatible with [design()].
#'
#' @export
#'
RDMSWTN <- function(erlang_shape = 1L, erlang_type = "none", posdrift = TRUE,
                    correlated = FALSE, correlate = c("times", "drifts")) {
  erlang_type <- match.arg(erlang_type, c("none", "local_kill", "global_kill", "local_guess", "local_kill_guess"))
  correlate <- match.arg(correlate)
  drift_corr <- correlated && correlate == "drifts"
  if (correlated && erlang_type != "none") {
    stop("Correlated RDMSWTN does not support guess or kill clocks; use erlang_type = \"none\".")
  }
  erlang_mixed <- identical(erlang_shape, "mixed")
  if (erlang_mixed && erlang_type == "global_kill") {
    stop("RDMSWTN erlang_shape = 'mixed' is currently implemented for local Erlang processes only.")
  }
  erlang_shape_cpp <- if (erlang_mixed) 3L else as.integer(erlang_shape)
  type_suffix <- switch(erlang_type,
    "none" = "",
    "local_kill" = "_LOCAL_KILL",
    "global_kill" = "_GLOBAL_KILL",
    "local_guess" = "_LOCAL_GUESS",
    "local_kill_guess" = "_LOCAL_KILL_GUESS"
  )
  base_name <- paste0(
    if (erlang_mixed) "RDMSWTN_EMIX" else if (erlang_shape_cpp >= 2L) "RDMSWTN_E2" else "RDMSWTN",
    type_suffix
  )
  has_guess <- erlang_type %in% c("local_guess", "local_kill_guess")
  has_kill  <- erlang_type %in% c("local_kill", "global_kill", "local_kill_guess")
  # posdrift = FALSE is the unrestricted-drift (IO) model: negative mean rates
  # are legal and give a defective finishing time, so v must be sampled on the
  # natural scale.  The log/exp scale is kept for posdrift = TRUE, where v > 0
  # is a hard constraint of the truncated drift law.
  .v <- .rdmswtn_v_scale(posdrift)
  p_types <- c(
    "v" = .v$p_type, "B" = log(1), "A" = log(0), "t0" = log(0),
    "s" = log(1), "sv" = log(0)
  )
  transform <- c(
    v = .v$transform, B = "exp", A = "exp", t0 = "exp",
    s = "exp", sv = "exp"
  )
  minmax <- cbind(
    v = .v$minmax, B = c(0, Inf), A = c(0, Inf),
    t0 = c(0.05, Inf), s = c(0, Inf), sv = c(0, Inf)
  )
  exception <- c(A = 0, v = 0, sv = 0)
  p_types  <- c(p_types,  mG = log(1))
  transform <- c(transform, mG = "exp")
  minmax   <- cbind(minmax, mG = c(1e-4, Inf))
  exception <- c(exception, mG = 0)

  p_types  <- c(p_types,  mK = log(1))
  transform <- c(transform, mK = "exp")
  minmax   <- cbind(minmax, mK = c(1e-4, Inf))
  exception <- c(exception, mK = 0)

  if (erlang_mixed) {
    p_types  <- c(p_types,  omega = qnorm(0.5))
    transform <- c(transform, omega = "pnorm")
    minmax   <- cbind(minmax, omega = c(0, 1))
    exception <- c(exception, omega = 0)
  }

  # pContaminant (omission) and pGuess (uniform outlier); see add_nuisance_pars().
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- .nuis$p_types; transform <- .nuis$transform
  minmax <- .nuis$minmax; exception <- .nuis$exception
  if (correlated) {
    p_types <- c(p_types, rho = qnorm(0.5))
    transform <- c(transform, rho = "pnorm")
    # Bounds are checked strictly by the compiled mapper.  The scaled-probit
    # round trip used for the advertised +/- .99 safety limits can land
    # exactly on an endpoint, which would otherwise reject those valid draws
    # and create an artificial discontinuity at the edge of the rho domain.
    rho_bound_eps <- 2 * .Machine$double.eps
    minmax <- cbind(minmax,
                    rho = c(-.99 - rho_bound_eps, .99 + rho_bound_eps))
    exception <- c(exception, rho = 0)
  }
  transform_spec <- list(func = transform)
  if (correlated) {
    transform_spec$lower <- c(rho = -1)
    transform_spec$upper <- c(rho = 1)
  }
  list(
    type = "RACE",
    c_name = paste0(if (posdrift) base_name else paste0(base_name, "_IO"),
                    if (!correlated) "" else if (drift_corr) "_CORRD" else "_CORR"),
    correlated = correlated,
    correlation_type = if (!correlated) NULL else if (drift_corr)
      "rdmswtn_drift_factor" else "rdmswtn_gaussian_copula",
    p_types = p_types,
    p_types_canonical = c("v", "B", "A", "t0", "s", "sv"),
    transform = transform_spec,
    bound = list(
      minmax = minmax,
      exception = exception
    ),
    Ttransform = function(pars, dadm) {
      if (drift_corr) {
        pars <- .apply_drift_factor_rho(pars, dadm, sv = pars[, "sv"],
                                        model = "RDMSWTNcorr")
      } else if (correlated) {
        .validate_rdmswtn_corr_rows(pars[, "rho"], dadm = dadm,
                                    posdrift = posdrift, sv = pars[, "sv"])
      }
      lambda_factor <- if (erlang_shape_cpp == 2L) 2 else 1
      n <- nrow(pars)
      mG_val <- pars[, "mG"]
      mK_val <- pars[, "mK"]
      lg <- if (has_guess) ifelse(mG_val <= 0, 0, lambda_factor / mG_val) else rep(0, n)
      lk <- if (has_kill)  ifelse(mK_val <= 0, 0, lambda_factor / mK_val) else rep(0, n)
      timed <- cbind(lambda_g = lg, lambda_k = lk)
      extra_drop <- c("v", "B", "A", "t0", "s", "sv", "mG", "mK", "omega")
      extra <- pars[, setdiff(colnames(pars), extra_drop), drop = FALSE]
      if (erlang_mixed) {
        pars <- cbind(
          pars[, c("v", "B", "A", "t0", "s", "sv"), drop = FALSE],
          timed,
          omega = pars[, "omega"],
          extra
        )
      } else {
        pars <- cbind(
          pars[, c("v", "B", "A", "t0", "s", "sv"), drop = FALSE],
          timed,
          extra
        )
      }
      pars <- cbind(pars, b = pars[, "B"] + pars[, "A"])
      pars
    },
    rfun = function(data = NULL, pars) {
      ok <- attr(pars, "ok")
      if (is.null(ok)) ok <- rep(TRUE, nrow(pars))
      .rfun_RDMSWTN(data$lR, pars, ok = ok, erlang_shape = erlang_shape_cpp,
               erlang_type = erlang_type, posdrift = posdrift,
               correlated = correlated, correlate = correlate)
    },
    dfun = function(rt, pars) dRDMSWTN(rt, pars, erlang = erlang_shape_cpp, posdrift = posdrift),
    pfun = function(rt, pars) pRDMSWTN(rt, pars, erlang = erlang_shape_cpp, posdrift = posdrift),
    log_likelihood = if (correlated) {
      function(pars, dadm, model, min_ll = log(1e-10)) {
        stop("RDMSWTNcorr likelihood is implemented in the C++ race path; use fast_path=TRUE.")
      }
    } else {
      function(pars, dadm, model, min_ll = log(1e-10)) {
        log_likelihood_race_missing(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
      }
    }
  )
}

# Sampling scale for the RDMSWTN mean drift.  With posdrift = TRUE the
# truncated normal drift law requires v > 0, so v is sampled on the log scale.
# With posdrift = FALSE negative mean rates are part of the model (they give a
# defective accumulator with intrinsic never-finish mass), so v is sampled on
# the natural scale, matching the LBA/BAwL IO convention.
.rdmswtn_v_scale <- function(posdrift) {
  if (posdrift) {
    list(p_type = log(1), transform = "exp", minmax = c(1e-3, Inf))
  } else {
    list(p_type = 1, transform = "identity", minmax = c(-Inf, Inf))
  }
}

# The Gaussian copula couples finishing-time ranks, which stays well defined
# when posdrift = FALSE makes those finishing times defective: uniforms above
# the marginal plateau are the atom at infinity, so rho also correlates the
# intrinsic omissions.  sv > 0 is reserved on that route: a between-trial drift
# SD combined with unrestricted drifts is intended to become a
# correlated-*draws* model (as in BAwLcorr), which is a different data
# generating process, so it is rejected rather than silently given copula
# semantics.
.check_rdmswtn_corr_io_sv <- function(correlated, posdrift, sv, tol = 1e-12) {
  if (posdrift || !any(correlated)) return(invisible(TRUE))
  if (is.null(sv)) {
    stop("RDMSWTNcorr with posdrift = FALSE requires sv, which was not supplied for validation.")
  }
  sv <- rep_len(sv, length(correlated))
  if (any(correlated & (!is.finite(sv) | sv > tol))) {
    stop(
      "RDMSWTNcorr with posdrift = FALSE requires sv = 0 on the correlated ",
      "rows: the copula couples defective finishing times, while sv > 0 with ",
      "unrestricted drifts is reserved for a correlated-drift model."
    )
  }
  invisible(TRUE)
}

# Row-level encoding shared by every correlated-*draws* race model (BAwLcorr
# and RDMSWTN(correlate = "drifts")).  The sampled coefficient is one direct
# correlation per cell; the compiled kernel turns a row's signed rho into a
# loading sign(rho) * sv * sqrt(|rho|) on a single shared standard-normal
# factor, so the correct racer is made the positive reference and the
# incorrect racers carry the cell sign.  The product of any correct/incorrect
# pair of loadings is then exactly rho.  Independent rows (a PM racer, say)
# must be made structurally zero by the rho design and its constants.
#
# The remap is idempotent -- the C++ side applies the same role mapping from
# lM -- but it has to happen here as well so the R simulator, which reads the
# row-level rho directly, describes the same process as the likelihood.
.apply_drift_factor_rho <- function(pars, dadm, sv = NULL, model = "model") {
  if (is.null(dadm) || !"lM" %in% names(dadm) || nrow(pars) != nrow(dadm)) {
    stop(model, " with correlate = \"drifts\" requires a matchfun-generated ",
         "lM role indicator in the expanded data.")
  }
  rho_cell <- pars[, "rho"]
  if (any(!is.finite(rho_cell)) || any(abs(rho_cell) > 1)) {
    stop(model, " requires finite natural-scale rho values in [-1, 1].")
  }
  correct <- as.character(dadm$lM) == "TRUE"
  if (anyNA(correct)) {
    stop(model, " matchfun produced missing lM values.")
  }
  # Correlating the draws is only meaningful when there are draws to
  # correlate: at sv = 0 the drift is deterministic and rho has no effect at
  # all, which would leave it unidentified rather than merely uninformative.
  if (!is.null(sv)) {
    coupled <- abs(rho_cell) > 1e-12
    if (any(coupled & (!is.finite(sv) | sv <= 1e-12))) {
      stop(model, " with correlate = \"drifts\" requires sv > 0 on the ",
           "correlated rows: the correlation acts on the between-trial drift ",
           "draws, so it is undefined at sv = 0. Use correlate = \"times\" ",
           "for a finishing-time copula.")
    }
  }
  pars[, "rho"] <- ifelse(correct, abs(rho_cell), rho_cell)
  pars
}

.validate_rdmswtn_corr_rows <- function(rho, dadm = NULL, posdrift = TRUE,
                                        sv = NULL, tol = 1e-12) {
  if (any(!is.finite(rho)) || any(abs(rho) > 1)) {
    stop("RDMSWTNcorr requires finite natural-scale rho values in [-1, 1].")
  }
  if (is.null(dadm) || is.null(dadm$lR)) {
    .check_rdmswtn_corr_io_sv(abs(rho) > tol, posdrift, sv, tol)
    return(invisible(TRUE))
  }
  n_lR <- length(levels(dadm$lR))
  if (n_lR < 1L || length(rho) %% n_lR != 0L) {
    stop("RDMSWTNcorr could not determine accumulator rows within trials.")
  }
  active <- rep(TRUE, length(rho))
  if ("RACE" %in% names(dadm) && !is.null(attr(dadm, "RACE_mask"))) {
    race_mask <- attr(dadm, "RACE_mask")
    if (length(race_mask) == length(rho)) active <- race_mask
  } else if ("RACE" %in% names(dadm)) {
    for (j in seq_len(length(rho) / n_lR)) {
      rows <- ((j - 1L) * n_lR + 1L):(j * n_lR)
      n_acc <- suppressWarnings(
        as.integer(as.character(dadm$RACE[rows[1L]])))
      if (is.finite(n_acc)) active[rows] <- seq_len(n_lR) <= n_acc
    }
  }
  .check_rdmswtn_corr_io_sv(active & abs(rho) > tol, posdrift, sv, tol)
  for (j in seq_len(length(rho) / n_lR)) {
    rows <- ((j - 1L) * n_lR + 1L):(j * n_lR)
    values <- rho[rows][active[rows] & abs(rho[rows]) > tol]
    if (length(values) > 2L ||
        (length(values) == 2L && abs(values[1L] - values[2L]) > tol)) {
      stop(
        "RDMSWTNcorr requires exactly one directly specified pair: at most two ",
        "active rows may have the same signed nonzero rho, and all other rows ",
        "must have rho = 0."
      )
    }
  }
  invisible(TRUE)
}

#' Correlated RDMSWTN Gaussian-Copula Race
#'
#' Convenience constructor for `RDMSWTN(..., correlated = TRUE)`. With the
#' default `correlate = "times"` the two participating accumulator finishing
#' times retain their ordinary RDMSWTN marginals and are coupled by a Gaussian
#' copula; with `correlate = "drifts"` it is the between-trial drift draws that
#' are correlated, as in [BAwLcorr()]. See [RDMSWTN()] for the participation
#' design and the interpretation of `rho` under each.
#'
#' @param erlang_shape Retained for constructor compatibility. Only the
#'   no-clock model is currently supported.
#' @param erlang_type Must be `"none"`.
#' @param posdrift Logical. Nonzero `rho` is supported for both settings. With
#'   `correlate = "times"`, `posdrift = FALSE` additionally requires `sv = 0`
#'   on the correlated rows.
#' @param correlate `"times"` (default) for the finishing-time copula, or
#'   `"drifts"` for correlated drift draws, which requires `sv > 0`.
#' @return A model list compatible with [design()].
#' @export
RDMSWTNcorr <- function(erlang_shape = 1L, erlang_type = "none",
                        posdrift = TRUE,
                        correlate = c("times", "drifts")) {
  RDMSWTN(erlang_shape = erlang_shape, erlang_type = erlang_type,
          posdrift = posdrift, correlated = TRUE,
          correlate = match.arg(correlate))
}

#' Time-Changed RDMSWTN Race
#'
#' `RDMSWTN_TT()` applies a finite linear exhaustion clock to each ordinary
#' RDMSWTN accumulator. For decision time `x = t - t0`, operational time is
#' `q(x) = x - x^2 / (2 * tau)` on `0 < x < tau`. The clock has budget
#' `Q = tau / 2`; the density is zero at and beyond `t0 + tau`, while the CDF
#' remains fixed at the ordinary RDMSWTN CDF evaluated at `Q`. The remaining
#' probability is genuine omission mass. The analytic functions and both
#' simulators also accept `tau = Inf`, which is the exact identity-clock limit
#' and recovers the ordinary RDMSWTN process without an exhaustion plateau.
#'
#' The parameters are `v`, `B`, `A`, `t0`, `s`, `sv`, and `tau`, with
#' `b = B + A`. `tau` is the width of the decision-time support after `t0`.
#' These seven parameters use exponential transforms. Their defaults are
#' `log(1)` for `v`, `B`, `s`, and `tau`, and the boundary value `log(0)` for
#' `A`, `t0`, and `sv`. The drift, start-point variability, diffusion scale,
#' and `b = B + A` convention are otherwise identical to [RDMSWTN()].
#' The optional `pContaminant` (omission) and `pGuess` (uniform outlier)
#' parameters are handled by the generic data pipeline. With
#' `correlated = TRUE`, `rho` follows the same contract as [RDMSWTN()]:
#' `correlate = "times"` couples exactly two active finishing-time marginals
#' with a Gaussian copula (and under `posdrift = FALSE` additionally requires
#' `sv = 0` on the correlated rows), while `correlate = "drifts"` correlates
#' the between-trial drift draws and requires `sv > 0`.
#'
#' @param posdrift Logical. If `TRUE` (default), truncate the between-trial
#'   normal drift distribution below at zero. If `FALSE`, allow unrestricted
#'   drifts and their additional intrinsic omission mass.
#' @param correlated Logical. If `TRUE`, include the correlation parameter
#'   `rho` for one directly specified accumulator pair.
#' @param correlate `"times"` (default) for the finishing-time copula, or
#'   `"drifts"` for correlated drift draws. See [RDMSWTN()].
#' @return A model list compatible with [design()].
#' @export
RDMSWTN_TT <- function(posdrift = TRUE, correlated = FALSE,
                       correlate = c("times", "drifts")) {
  correlate <- match.arg(correlate)
  drift_corr <- correlated && correlate == "drifts"
  # See .rdmswtn_v_scale(): posdrift = FALSE samples v on the natural scale.
  .v <- .rdmswtn_v_scale(posdrift)
  p_types <- c(
    v = .v$p_type, B = log(1), A = log(0), t0 = log(0),
    s = log(1), sv = log(0), tau = log(1),
    pContaminant = qnorm(0), pGuess = qnorm(0)
  )
  transform <- c(
    v = .v$transform, B = "exp", A = "exp", t0 = "exp",
    s = "exp", sv = "exp", tau = "exp",
    pContaminant = "pnorm", pGuess = "pnorm"
  )
  minmax <- cbind(
    v = .v$minmax, B = c(0, Inf), A = c(0, Inf),
    t0 = c(0.05, Inf), s = c(0, Inf), sv = c(0, Inf),
    tau = c(1e-4, Inf), pContaminant = c(0.001, 0.999),
    pGuess = c(0.001, 0.999)
  )
  exception <- c(A = 0, v = 0, sv = 0, pContaminant = 0, pGuess = 0)
  if (correlated) {
    p_types <- c(p_types, rho = qnorm(0.5))
    transform <- c(transform, rho = "pnorm")
    rho_bound_eps <- 2 * .Machine$double.eps
    minmax <- cbind(
      minmax, rho = c(-.99 - rho_bound_eps, .99 + rho_bound_eps)
    )
    exception <- c(exception, rho = 0)
  }
  transform_spec <- list(func = transform)
  if (correlated) {
    transform_spec$lower <- c(rho = -1)
    transform_spec$upper <- c(rho = 1)
  }
  list(
    type = "RACE",
    c_name = paste0(
      if (posdrift) "RDMSWTN_TT" else "RDMSWTN_TT_IO",
      if (!correlated) "" else if (drift_corr) "_CORRD" else "_CORR"
    ),
    correlated = correlated,
    correlation_type = if (!correlated) NULL else if (drift_corr)
      "rdmswtn_drift_factor" else "rdmswtn_gaussian_copula",
    p_types = p_types,
    p_types_canonical = c("v", "B", "A", "t0", "s", "sv", "tau"),
    transform = transform_spec,
    bound = list(minmax = minmax, exception = exception),
    Ttransform = function(pars, dadm) {
      if (drift_corr) {
        pars <- .apply_drift_factor_rho(pars, dadm, sv = pars[, "sv"],
                                        model = "RDMSWTN_TTcorr")
      } else if (correlated) {
        .validate_rdmswtn_corr_rows(
          pars[, "rho"], dadm = dadm, posdrift = posdrift, sv = pars[, "sv"]
        )
      }
      canonical <- c("v", "B", "A", "t0", "s", "sv", "tau")
      extra <- pars[, setdiff(colnames(pars), canonical), drop = FALSE]
      pars <- cbind(pars[, canonical, drop = FALSE], extra)
      cbind(pars, b = pars[, "B"] + pars[, "A"])
    },
    rfun = function(data = NULL, pars) {
      ok <- attr(pars, "ok")
      if (is.null(ok)) ok <- rep(TRUE, nrow(pars))
      .rfun_RDMSWTN_TT(
        data$lR, pars, ok = ok, posdrift = posdrift,
        correlated = correlated, correlate = correlate
      )
    },
    dfun = function(rt, pars) {
      dRDMSWTN_TT(rt, pars, posdrift = posdrift)
    },
    pfun = function(rt, pars) {
      pRDMSWTN_TT(rt, pars, posdrift = posdrift)
    },
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      stop("RDMSWTN_TT likelihoods are implemented only in the compiled race path; use fast_path=TRUE.")
    }
  )
}

#' Correlated Time-Changed RDMSWTN Race
#'
#' Convenience constructor for `RDMSWTN_TT(correlated = TRUE)`.
#'
#' @param posdrift Logical. Active nonzero `rho` is supported for both
#'   settings. With `correlate = "times"`, `posdrift = FALSE` additionally
#'   requires `sv = 0` on the correlated rows.
#' @param correlate `"times"` (default) for the finishing-time copula, or
#'   `"drifts"` for correlated drift draws, which requires `sv > 0`.
#' @return A model list compatible with [design()].
#' @export
RDMSWTN_TTcorr <- function(posdrift = TRUE,
                           correlate = c("times", "drifts")) {
  RDMSWTN_TT(posdrift = posdrift, correlated = TRUE,
             correlate = match.arg(correlate))
}

.rdmswtn_tt_pars <- function(rt, pars) {
  if (is.null(dim(pars)) || (nrow(pars) == 1L && length(rt) > 1L)) {
    original_names <- names(pars)
    if (is.null(original_names)) original_names <- colnames(pars)
    pars <- matrix(
      pars, nrow = length(rt), ncol = length(pars),
      dimnames = list(NULL, original_names), byrow = TRUE
    )
  }
  if (length(rt) == 1L && nrow(pars) > 1L) rt <- rep(rt, nrow(pars))
  if (!"b" %in% colnames(pars) &&
      all(c("B", "A") %in% colnames(pars))) {
    pars <- cbind(pars, b = pars[, "B"] + pars[, "A"])
  }
  required <- c("v", "b", "A", "t0", "sv", "tau")
  if (!all(required %in% colnames(pars))) {
    stop("RDMSWTN_TT requires parameter columns ", paste(required, collapse = ", "), ".")
  }
  if (!"s" %in% colnames(pars)) pars <- cbind(pars, s = 1)
  list(rt = rt, pars = pars)
}

dRDMSWTN_TT <- function(rt, pars, posdrift = TRUE, log = FALSE) {
  input <- .rdmswtn_tt_pars(rt, pars)
  pars <- input$pars
  dRDMSWTN_TT_cpp(
    input$rt, pars[, "v"], pars[, "b"], pars[, "A"], pars[, "s"],
    pars[, "t0"], pars[, "sv"], pars[, "tau"],
    log_out = log, posdrift = posdrift
  )
}

pRDMSWTN_TT <- function(rt, pars, posdrift = TRUE, log.p = FALSE) {
  input <- .rdmswtn_tt_pars(rt, pars)
  pars <- input$pars
  pRDMSWTN_TT_cpp(
    input$rt, pars[, "v"], pars[, "b"], pars[, "A"], pars[, "s"],
    pars[, "t0"], pars[, "sv"], pars[, "tau"],
    log_out = log.p, posdrift = posdrift
  )
}

.rdmswtn_tt_qinv_R <- function(u, tau) {
  2 * u / (1 + sqrt(pmax(0, 1 - 2 * u / tau)))
}

rRDMSWTN_TT <- function(lR, pars, ok = rep(TRUE, nrow(pars)),
                        posdrift = TRUE, drift_override = NULL) {
  if (!is.null(attr(pars, "ok"))) ok <- attr(pars, "ok")
  nr <- length(levels(lR))
  if (nr < 1L || nrow(pars) %% nr != 0L) {
    stop("RDMSWTN_TT requires rows grouped by accumulator within trial.")
  }
  required <- c("v", "b", "A", "t0", "sv", "tau")
  if (!all(required %in% colnames(pars))) {
    stop("RDMSWTN_TT requires parameter columns ", paste(required, collapse = ", "), ".")
  }
  if (!"s" %in% colnames(pars)) pars <- cbind(pars, s = 1)
  finish <- rep(Inf, nrow(pars))
  active <- which(ok)
  if (length(active)) {
    operational <- rSWTN(
      length(active), b = pars[active, "b"], v = pars[active, "v"],
      A = pars[active, "A"], sv = pars[active, "sv"],
      s = pars[active, "s"], k = 0, posdrift = posdrift,
      drift_override = if (is.null(drift_override)) NULL else drift_override[active]
    )
    Q <- pars[active, "tau"] / 2
    respond <- is.finite(operational) & operational <= Q
    if (any(respond)) {
      rows <- active[respond]
      finish[rows] <- pars[rows, "t0"] +
        .rdmswtn_tt_qinv_R(operational[respond], pars[rows, "tau"])
    }
  }
  n_trials <- nrow(pars) / nr
  dt <- matrix(finish, nrow = nr)
  bad <- colSums(is.finite(dt)) == 0L
  response <- max.col(-t(dt), ties.method = "first")
  pick <- cbind(response, seq_len(n_trials))
  out <- data.frame(
    R = factor(levels(lR)[response], levels = levels(lR)),
    rt = dt[pick]
  )
  out$R[bad] <- NA
  out$rt[bad] <- Inf
  .apply_timed_guess_winner(out, levels(lR))
}

.qRDMSWTN_TT_operational <- function(u, Q, pars, posdrift = TRUE) {
  if (u <= 0) return(0)
  FQ <- prdmswtn(
    Q, pars[1L, "v"], pars[1L, "b"], pars[1L, "A"],
    pars[1L, "s"], 0, pars[1L, "sv"], 0, 0,
    posdrift = posdrift
  )
  if (u >= FQ) return(Q)
  stats::uniroot(
    function(x) {
      prdmswtn(
        x, pars[1L, "v"], pars[1L, "b"], pars[1L, "A"],
        pars[1L, "s"], 0, pars[1L, "sv"], 0, 0,
        posdrift = posdrift
      ) - u
    },
    c(0, Q), tol = 1e-10
  )$root
}

rRDMSWTN_TT_corr <- function(lR, pars, ok = rep(TRUE, nrow(pars)),
                             posdrift = TRUE) {
  if (!is.null(attr(pars, "ok"))) ok <- attr(pars, "ok")
  if (!"rho" %in% colnames(pars)) {
    stop("RDMSWTN_TTcorr requires parameter column 'rho'.")
  }
  if (!"s" %in% colnames(pars)) pars <- cbind(pars, s = 1)
  nr <- length(levels(lR))
  if (nr < 1L || nrow(pars) %% nr != 0L) {
    stop("RDMSWTN_TTcorr requires rows grouped by accumulator within trial.")
  }
  rho <- pars[, "rho"]
  # posdrift = TRUE here runs the shape checks only; the IO/sv rule is applied
  # below against the rows this call actually simulates.
  .validate_rdmswtn_corr_rows(rho, posdrift = TRUE)
  .check_rdmswtn_corr_io_sv(ok & abs(rho) > 1e-12, posdrift, pars[, "sv"])
  if (!any(ok & abs(rho) > 1e-12)) {
    return(rRDMSWTN_TT(lR, pars, ok = ok, posdrift = posdrift))
  }
  n_trials <- nrow(pars) / nr
  u <- runif(nrow(pars))
  for (tr in seq_len(n_trials)) {
    rows <- ((tr - 1L) * nr + 1L):(tr * nr)
    pair <- rows[ok[rows] & abs(rho[rows]) > 1e-12]
    if (length(pair) > 2L ||
        (length(pair) == 2L && abs(rho[pair[1L]] - rho[pair[2L]]) > 1e-12)) {
      stop("RDMSWTN_TTcorr requires at most two active rows with one signed nonzero rho.")
    }
    if (length(pair) == 2L) {
      z1 <- rnorm(1)
      z2 <- rho[pair[1L]] * z1 +
        sqrt(max(0, 1 - rho[pair[1L]]^2)) * rnorm(1)
      u[pair] <- pnorm(c(z1, z2))
    }
  }
  finish <- rep(Inf, nrow(pars))
  for (r in which(ok)) {
    Q <- pars[r, "tau"] / 2
    FQ <- pRDMSWTN_TT(Inf, pars[r, , drop = FALSE], posdrift = posdrift)
    if (u[r] <= FQ) {
      operational <- .qRDMSWTN_TT_operational(
        u[r], Q, pars[r, , drop = FALSE], posdrift = posdrift
      )
      finish[r] <- pars[r, "t0"] +
        .rdmswtn_tt_qinv_R(operational, pars[r, "tau"])
    }
  }
  dt <- matrix(finish, nrow = nr)
  bad <- colSums(is.finite(dt)) == 0L
  response <- max.col(-t(dt), ties.method = "first")
  pick <- cbind(response, seq_len(n_trials))
  out <- data.frame(
    R = factor(levels(lR)[response], levels = levels(lR)),
    rt = dt[pick]
  )
  out$R[bad] <- NA
  out$rt[bad] <- Inf
  .apply_timed_guess_winner(out, levels(lR))
}

#' RDMSWTN Logical Rules Model
#'
#' Use the current RDMSWTN accumulator distribution for the logical-rules
#' likelihood. The logical-rules evaluator is shared with the other race
#' models; adding the `LogicalRules` model-name marker routes the compiled
#' likelihood through the RDMSWTN adapter while retaining the ordinary
#' RDMSWTN parameter and clock contract.
#'
#' Logical rules use accumulator roles `A`, `B`, `n_A`, and `n_B` and the
#' `LogicalRule` data column. Supported rules are `OR`, `AND`, `XOR`, `ID`,
#' `OR_DETECTION_ANALYTIC`, and `OR_DETECTION_GNG`.
#'
#' @param erlang_shape Integer `1` for exponential clocks, `2` for Erlang-2,
#'   or `"mixed"` for the Erlang-1/Erlang-2 mixture.
#' @param erlang_type Clock configuration passed to [RDMSWTN()].
#' @param posdrift Logical. If `TRUE` (default), truncate the between-trial
#'   normal drift distribution below at zero; if `FALSE`, use untruncated
#'   drifts and allow intrinsic omissions.
#' @return A model list compatible with [design()].
#' @export
#'
LogicalRulesRDMSWTN <- function(erlang_shape = 1L, erlang_type = "none",
                                posdrift = TRUE) {
  out <- RDMSWTN(erlang_shape = erlang_shape,
                 erlang_type = erlang_type,
                 posdrift = posdrift)
  out$c_name <- paste0(out$c_name, "_LogicalRules")
  out$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
    stop("LogicalRulesRDMSWTN: R likelihood path not implemented. Use the compiled path.")
  }
  out
}

.rdmswtn_erlang_omega <- function(pars, erlang) {
  n <- nrow(pars)
  if (is.null(n)) n <- length(pars) > 0
  if (as.integer(erlang) <= 1L) return(rep(1, n))
  if (as.integer(erlang) == 2L) return(rep(0, n))
  if (!"omega" %in% colnames(pars)) {
    stop("mixed RDMSWTN Erlang mode requires parameter column 'omega'.")
  }
  pars[, "omega"]
}

dRDMSWTN <- function(rt, pars, erlang = 1L, posdrift = TRUE) {
  if (is.null(dim(pars)) || (dim(pars)[1] == 1 & length(rt) > 1)) {
    original_names <- names(pars)
    if (is.null(original_names)) original_names <- colnames(pars)
    pars <- matrix(pars,
      nrow = length(rt), ncol = length(pars),
      dimnames = list(NULL, original_names), byrow = TRUE
    )
  }
  if (length(rt) == 1 && !is.null(dim(pars)) && nrow(pars) > 1) {
    rt <- rep(rt, nrow(pars))
  }
  out <- rep(NaN, length(rt))
  if (!all(c("lambda_g", "lambda_k") %in% colnames(pars))) {
    stop("RDMSWTN requires parameter columns 'lambda_g' and 'lambda_k'.")
  }
  erl <- (pars[, "lambda_g"] > 0) | (pars[, "lambda_k"] > 0)
  ok <- (rt > 0) & ((rt > pars[, "t0"]) | erl)
  ok[is.na(ok)] <- FALSE
  if (any(ok)) {
    s_ok <- if ("s" %in% colnames(pars)) pars[ok, "s", drop = FALSE] else 1
    out[ok] <- dSWTNspv(rt[ok],
      v = pars[ok, "v", drop = FALSE], b = pars[ok, "b", drop = FALSE],
      A = pars[ok, "A", drop = FALSE],
      s = s_ok,
      t0 = pars[ok, "t0", drop = FALSE],
      sv = pars[ok, "sv", drop = FALSE], lambda_g = pars[ok, "lambda_g", drop = FALSE],
      lambda_k = pars[ok, "lambda_k", drop = FALSE],
      kill_shape = erlang, erlang_omega = .rdmswtn_erlang_omega(pars[ok, , drop = FALSE], erlang),
      posdrift = posdrift
    )
  }
  out
}

pRDMSWTN <- function(rt, pars, erlang = 1L, posdrift = TRUE) {
  if (is.null(dim(pars)) || (dim(pars)[1] == 1 & length(rt) > 1)) {
    original_names <- names(pars)
    if (is.null(original_names)) original_names <- colnames(pars)
    pars <- matrix(pars,
      nrow = length(rt), ncol = length(pars),
      dimnames = list(NULL, original_names), byrow = TRUE
    )
  }
  if (length(rt) == 1 && !is.null(dim(pars)) && nrow(pars) > 1) {
    rt <- rep(rt, nrow(pars))
  }
  out <- rep(NaN, length(rt))
  if (!all(c("lambda_g", "lambda_k") %in% colnames(pars))) {
    stop("RDMSWTN requires parameter columns 'lambda_g' and 'lambda_k'.")
  }
  erl <- (pars[, "lambda_g"] > 0) | (pars[, "lambda_k"] > 0)
  ok <- (rt > 0) & ((rt > pars[, "t0"]) | erl)
  ok[is.na(ok)] <- FALSE
  if (any(ok)) {
    s_ok <- if ("s" %in% colnames(pars)) pars[ok, "s", drop = FALSE] else 1
    out[ok] <- pSWTNspv(rt[ok],
      v = pars[ok, "v", drop = FALSE], b = pars[ok, "b", drop = FALSE],
      A = pars[ok, "A", drop = FALSE],
      s = s_ok,
      t0 = pars[ok, "t0", drop = FALSE],
      sv = pars[ok, "sv", drop = FALSE], lambda_g = pars[ok, "lambda_g", drop = FALSE],
      lambda_k = pars[ok, "lambda_k", drop = FALSE],
      kill_shape = erlang, erlang_omega = .rdmswtn_erlang_omega(pars[ok, , drop = FALSE], erlang),
      posdrift = posdrift
    )
  }
  out
}

rRDMSWTN <- function(lR, pars, p_types = c("v", "b", "A", "t0", "sv", "lambda_g", "lambda_k"),
                     ok = rep(TRUE, dim(pars)[1]), erlang_shape = 1L, erlang_type = "none",
                     posdrift = TRUE, drift_override = NULL) {
  if (!is.null(attr(pars, "ok"))) ok <- attr(pars, "ok")
  if (is.null(dim(pars)) || (dim(pars)[1] == 1 & length(lR) > 1)) {
    original_names <- names(pars)
    if (is.null(original_names)) original_names <- colnames(pars)
    pars <- matrix(pars,
      nrow = length(lR), ncol = length(pars),
      dimnames = list(NULL, original_names), byrow = TRUE
    )
  }
  if (!all(c("lambda_g", "lambda_k") %in% colnames(pars))) {
    stop("RDMSWTN requires parameter columns 'lambda_g' and 'lambda_k'.")
  }
  erlang_omega_all <- .rdmswtn_erlang_omega(pars, erlang_shape)
  if (!all(p_types %in% dimnames(pars)[[2]])) {
    stop("pars must have columns ", paste(p_types, collapse = " "))
  }
  if (!("s" %in% dimnames(pars)[[2]])) {
    pars <- cbind(pars, s = 1)
  }
  pars[, "b"][pars[, "b"] < 0] <- 0
  pars[, "A"][pars[, "A"] < 0] <- 0
  bad <- rep(NA, length(lR) / length(levels(lR)))
  out <- data.frame(R = bad, rt = bad)
  nr <- length(levels(lR))
  n_trials <- nrow(pars) / nr
  dt <- matrix(Inf, nrow = nr, ncol = n_trials)
  t0 <- pars[, "t0"]

  guess  <- erlang_type == "local_guess" || erlang_type == "local_kill_guess"
  global <- erlang_type == "global_kill"
  local_kill <- erlang_type == "local_kill" || erlang_type == "local_kill_guess"

  # global_kill: draw one shared Erlang timer per trial
  if (global) {
    lambda_mat <- matrix(pars[, "lambda_k"], nrow = nr)
    if (any(apply(lambda_mat, 2, function(x) length(unique(x)) > 1))) {
      stop("global_kill requires lambda_k to be constant across accumulators")
    }
    lambda_trials <- lambda_mat[1, ]
    tk <- rep(Inf, n_trials)
    kill_ok <- !is.na(lambda_trials) & lambda_trials > 0
    if (any(kill_ok)) {
      tk[kill_ok] <- rgamma(sum(kill_ok), shape = erlang_shape, rate = lambda_trials[kill_ok])
    }
  }

  pars_all <- pars
  pars <- pars[ok, , drop = FALSE]
  # Local Erlang clocks are handled on the raw-time axis below, after t0 is
  # added to the evidence-accumulation finish times.
  k_vec <- rep(0, nrow(pars))
  dt[ok] <- rSWTN(sum(ok),
    b = pars[, "b"], v = pars[, "v"], A = pars[, "A"], sv = pars[, "sv"],
    s = pars[, "s"],
    k = k_vec, erlang = erlang_shape, erlang_omega = erlang_omega_all[ok],
    posdrift = posdrift,
    drift_override = if (is.null(drift_override)) NULL else drift_override[ok]
  )
  # Put EAM on the same raw-time axis as Erlang clocks.
  dt <- dt + matrix(t0, nrow = nr)
  if (guess || local_kill) {
    tg_local <- matrix(Inf, nrow = nr, ncol = n_trials)
    tk_local <- matrix(Inf, nrow = nr, ncol = n_trials)

    if (guess) {
      lambda_g_local <- matrix(0, nrow = nr, ncol = n_trials)
      active_guess <- matrix(ok, nrow = nr) & (levels(lR) != "nogo")
      lambda_g_local[active_guess] <- pars_all[, "lambda_g"][active_guess]
      guess_ok_mat <- lambda_g_local > 0
      if (any(guess_ok_mat)) {
        shape_g <- if (as.integer(erlang_shape) == 3L) {
          ifelse(runif(sum(guess_ok_mat)) <= matrix(erlang_omega_all, nrow = nr)[guess_ok_mat], 1L, 2L)
        } else {
          as.integer(erlang_shape)
        }
        rate_g <- if (as.integer(erlang_shape) == 3L) {
          ifelse(shape_g == 2L, 2 * lambda_g_local[guess_ok_mat], lambda_g_local[guess_ok_mat])
        } else {
          lambda_g_local[guess_ok_mat]
        }
        tg_local[guess_ok_mat] <- rgamma(sum(guess_ok_mat), shape = shape_g,
                                         rate = rate_g)
      }
    }

    if (local_kill) {
      lambda_k_local <- matrix(0, nrow = nr, ncol = n_trials)
      lambda_k_local[matrix(ok, nrow = nr)] <- pars[, "lambda_k"]
      kill_ok_mat <- lambda_k_local > 0
      if (any(kill_ok_mat)) {
        shape_k <- if (as.integer(erlang_shape) == 3L) {
          ifelse(runif(sum(kill_ok_mat)) <= matrix(erlang_omega_all, nrow = nr)[kill_ok_mat], 1L, 2L)
        } else {
          as.integer(erlang_shape)
        }
        rate_k <- if (as.integer(erlang_shape) == 3L) {
          ifelse(shape_k == 2L, 2 * lambda_k_local[kill_ok_mat], lambda_k_local[kill_ok_mat])
        } else {
          lambda_k_local[kill_ok_mat]
        }
        tk_local[kill_ok_mat] <- rgamma(sum(kill_ok_mat), shape = shape_k,
                                        rate = rate_k)
      }
    }

    if (guess && local_kill) {
      guess_win = tg_local<dt & tg_local<tk_local
      dt_candidate <- pmin(dt, tg_local)
      dt <- ifelse(dt_candidate < tk_local, dt_candidate, Inf)
    } else if (guess) {
      guess_win = tg_local<dt
      dt <- pmin(dt, tg_local)
    } else {
      dt <- ifelse(dt < tk_local, dt, Inf)
    }
  }

  if (global) {
    # Global kill: shared timer fires → no response this trial.
    is_killed <- tk < apply(dt, 2, min)
    if (any(is_killed)) {
      dt[, is_killed] <- Inf
    }
  }

  bad_col <- apply(dt, 2, function(x) all(is.infinite(x)))
  R <- apply(dt, 2, which.min)
  pick <- cbind(R, 1:dim(dt)[2])
  rt <- dt[pick]
  out$R <- levels(lR)[R]
  out$R <- factor(out$R, levels = levels(lR))
  out$rt <- rt
  out$R[bad_col] <- NA
  out$rt[bad_col] <- Inf
  out <- .apply_timed_guess_winner(out, levels(lR))
  if (guess) {
    out$isTime <- rep(NA, n_trials)
    out$isTime[!bad_col] <- guess_win[pick][!bad_col]
  }
  out
}

.qRDMSWTN_row <- function(u, pars, posdrift = TRUE) {
  if (!(is.finite(u) && u > 0 && u < 1)) {
    stop("RDMSWTNcorr inverse CDF requires a probability strictly between zero and one.")
  }
  cdf <- function(t) {
    if (t <= pars[1L, "t0"]) return(-u)
    pRDMSWTN(t, pars, erlang = 1L, posdrift = posdrift) - u
  }
  if (!posdrift) {
    # Defective marginal: uniforms above the plateau F(Inf) = p_hit are the
    # atom at infinity, i.e. this accumulator never finishes.
    cap <- pRDMSWTN(Inf, pars, erlang = 1L, posdrift = FALSE)
    if (!is.finite(cap)) {
      stop("RDMSWTNcorr: non-finite marginal hit probability while inverting a finishing-time quantile.")
    }
    if (u >= cap) return(Inf)
  }
  lower <- 0
  upper <- max(1, pars[1L, "t0"] + 1)
  value <- cdf(upper)
  for (iter in seq_len(80L)) {
    if (is.finite(value) && value >= 0) break
    upper <- upper * 2
    if (!is.finite(upper) || upper > 1e12) {
      stop("RDMSWTNcorr could not bracket a finishing-time quantile; check the marginal parameters.")
    }
    value <- cdf(upper)
  }
  if (!is.finite(value) || value < 0) {
    stop("RDMSWTNcorr could not bracket a finishing-time quantile; check the marginal parameters.")
  }
  uniroot(cdf, c(lower, upper), tol = 1e-10)$root
}

rRDMSWTN_corr <- function(lR, pars, ok = rep(TRUE, nrow(pars)),
                          posdrift = TRUE) {
  if (!is.null(attr(pars, "ok"))) ok <- attr(pars, "ok")
  if (!"rho" %in% colnames(pars)) {
    stop("RDMSWTNcorr requires parameter column 'rho'.")
  }
  if (!all(c("lambda_g", "lambda_k") %in% colnames(pars)) ||
      any(pars[, "lambda_g"] != 0 | pars[, "lambda_k"] != 0)) {
    stop("RDMSWTNcorr does not support guess or kill clocks.")
  }
  nr <- length(levels(lR))
  if (nr < 1L || nrow(pars) %% nr != 0L) {
    stop("RDMSWTNcorr requires rows grouped by accumulator within trial.")
  }
  rho <- pars[, "rho"]
  # posdrift = TRUE here runs the shape checks only; the IO/sv rule is applied
  # below against the rows this call actually simulates.
  .validate_rdmswtn_corr_rows(rho, posdrift = TRUE)
  .check_rdmswtn_corr_io_sv(ok & abs(rho) > 1e-12, posdrift, pars[, "sv"])
  if (!any(ok & abs(rho) > 1e-12)) {
    return(rRDMSWTN(lR, pars, ok = ok, erlang_shape = 1L,
                    erlang_type = "none", posdrift = posdrift))
  }

  n_trials <- nrow(pars) / nr
  u <- runif(nrow(pars))
  for (tr in seq_len(n_trials)) {
    rows <- ((tr - 1L) * nr + 1L):(tr * nr)
    pair <- rows[ok[rows] & abs(rho[rows]) > 1e-12]
    if (length(pair) > 2L ||
        (length(pair) == 2L && abs(rho[pair[1L]] - rho[pair[2L]]) > 1e-12)) {
      stop("RDMSWTNcorr requires at most two active rows with the same signed nonzero rho in each trial.")
    }
    if (length(pair) == 2L) {
      z1 <- rnorm(1)
      z2 <- rho[pair[1L]] * z1 +
        sqrt(max(0, 1 - rho[pair[1L]]^2)) * rnorm(1)
      u[pair] <- pnorm(c(z1, z2))
    }
  }

  finish <- rep(Inf, nrow(pars))
  active <- which(ok)
  for (r in active) {
    finish[r] <- .qRDMSWTN_row(u[r], pars[r, , drop = FALSE],
                               posdrift = posdrift)
  }
  dt <- matrix(finish, nrow = nr)
  bad <- apply(dt, 2L, function(x) all(is.infinite(x)))
  response <- apply(dt, 2L, which.min)
  pick <- cbind(response, seq_len(n_trials))
  out <- data.frame(
    R = factor(levels(lR)[response], levels = levels(lR)),
    rt = dt[pick]
  )
  out$R[bad] <- NA
  out$rt[bad] <- Inf
  .apply_timed_guess_winner(out, levels(lR))
}
