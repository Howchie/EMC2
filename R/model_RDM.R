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

rWald <- function(n, B, v, A, posdrift = TRUE)
# random function for single accumulator
{
  rwaldt <- function(n, k, l, tiny = 1e-6) {
    # random sample of n from a Wald (or Inverse Gaussian)
    # k = criterion, l = rate, assumes sigma=1 Browninan motion
    # about same speed as statmod rinvgauss

    rlevy <- function(n = 1, m = 0, c = 1) {
      if (any(c < 0)) stop("c must be positive")
      c / qnorm(1 - runif(n) / 2)^2 + m
    }

    flag <- l > abs(tiny)
    x <- rep(NA, times = n)

    x[!flag] <- rlevy(sum(!flag), 0, k[!flag]^2)
    mu <- k / l
    lambda <- k^2

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
    out[pos] <- rwaldt(npos, k = bs, l = v[pos])
  }

  # negative drift with posdrift=FALSE: defective Wald via Bernoulli(p_hit)
  # Conditional FPT given hitting equals Wald with |v| (Girsanov / time-reversal)
  if (any(neg)) {
    nneg <- sum(neg)
    if (!posdrift) { # sample bernoulli hitting probability and use the absolute value of v for the finite finishes
      bs_neg <- B[neg] + runif(nneg, 0, A[neg])
      p_hit <- exp(2 * v[neg] * bs_neg)  # v < 0, bs > 0 → 0 < p_hit < 1
      hit <- as.logical(rbinom(nneg, 1, p_hit))
      if (any(hit)) {
        out[which(neg)[hit]] <- rwaldt(sum(hit), k = bs_neg[hit], l = abs(v[neg][hit]))
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
  R <- max.col(-t(dt), ties.method = "first")
  pick <- cbind(R, 1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
  rt <- matrix(t0, nrow = nr)[pick] + dt[pick]
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
#'
#'
#' All parameters are estimated on the log scale.
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
    p_types = c("v" = log(1), "B" = log(1), "A" = log(0), "t0" = log(0), "s" = log(1), "pContaminant" = qnorm(0)),
    p_types_canonical = c("v", "B", "A", "t0", "s"),
    transform = list(func = c(v = "exp", B = "exp", A = "exp", t0 = "exp", s = "exp", pContaminant = "pnorm")),
    bound = list(
      minmax = cbind(v = c(1e-3, Inf), B = c(0, Inf), A = c(1e-4, Inf), t0 = c(0.05, Inf), s = c(0, Inf), pContaminant = c(0.001, 0.999)),
      exception = c(A = 0, v = 0, pContaminant = 0)
    ),
    # Trial dependent parameter transform
    Ttransform = function(pars, dadm) {
      pars <- cbind(pars, b = pars[, "B"] + pars[, "A"])
      pars
    },
    # Random function for racing accumulators
    rfun = function(data = NULL, pars) rRDM(data$lR, pars, ok = attr(pars, "ok")),
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
      kill_shape = erlang
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
      kill_shape = erlang
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
      tk[kill_ok] <- rgamma(sum(kill_ok), shape = erlang_shape, rate = lambda_trials[kill_ok])
    }
  }

  pars_ok <- pars[ok, , drop = FALSE]
  if (nrow(pars_ok) > 0) {
    # In the combined local kill+guess case, handle both clocks jointly below.
    k_vec <- if (local_kill && !guess) pars_ok[, "lambda_k"] else rep(0, nrow(pars_ok))
    dt[ok] <- rGBM_killed(sum(ok),
      b = pars_ok[, "b"], v = pars_ok[, "v"], A = pars_ok[, "A"],
      s = pars_ok[, "s"], k = k_vec, erlang = erlang_shape
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
        tg_local[guess_ok_mat] <- rgamma(sum(guess_ok_mat), shape = erlang_shape,
                                         rate = lambda_g_local[guess_ok_mat])
      }
    }

    if (local_kill) {
      lambda_k_local <- matrix(0, nrow = nr, ncol = n_trials)
      lambda_k_local[matrix(ok, nrow = nr)] <- pars[ok, "lambda_k"]
      kill_ok_mat <- lambda_k_local > 0
      if (any(kill_ok_mat)) {
        tk_local[kill_ok_mat] <- rgamma(sum(kill_ok_mat), shape = erlang_shape,
                                        rate = lambda_k_local[kill_ok_mat])
      }
    }

    if (guess && local_kill) {
      dt_candidate <- pmin(dt, tg_local)
      dt <- ifelse(dt_candidate < tk_local, dt_candidate, Inf)
    } else if (guess) {
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
  out
}

rGBM_killed <- function(n, b, v, A, s = 1, k = 0, erlang = 1L) {
  out <- rGBM(n, b = b, v = v, A = A, s = s)
  kill_idx <- which(k > 0)
  if (length(kill_idx) > 0) {
    tk <- rgamma(length(kill_idx), shape = erlang, rate = k[kill_idx])
    out[kill_idx[tk <= out[kill_idx]]] <- Inf
  }
  out
}

rSWTN <- function(n, b, v, A, sv, k = 0, erlang = 1L, posdrift = TRUE) {
  if (n <= 0) return(numeric(0))
  b <- rep(b, length.out = n)
  v <- rep(v, length.out = n)
  A <- rep(A, length.out = n)
  sv <- rep(sv, length.out = n)
  k <- rep(k, length.out = n)
  out <- rep(Inf, n)
  # For sv > 0 and posdrift=TRUE, draw per-trial drifts from N(v, sv^2)
  # truncated at zero. Otherwise draw from the full normal and let rWald
  # handle defective negative-drift finite hits when posdrift=FALSE.
  # For sv == 0 and v < 0 with posdrift=FALSE: Bernoulli(p_hit) sampling in rWald.
  v_draw <- v
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
  
  out <- rWald(n, B = b - A, v = v_draw, A = A, posdrift = posdrift)
  kill_idx <- which(k > 0)
  if (length(kill_idx) > 0) {
    tk <- rgamma(length(kill_idx), shape = erlang, rate = k[kill_idx])
    out[kill_idx[tk <= out[kill_idx]]] <- Inf
  }
  out
}

#' RDMGBM Model
#'
#' Racing geometric Brownian first-passage model with start-point variability.
#' Supports exponential killing (omission / guessing) via the `lambda_g` and `lambda_k` parameters.
#'
#' @param erlang_shape integer shape of the killing process (1 = exponential, 2 = Erlang-2)
#' @param erlang_type string, one of "none", "local_kill", "global_kill", "local_guess", "local_kill_guess"
#'
#' @export
#'
RDMGBM <- function(erlang_shape = 1L, erlang_type = "none") {
  erlang_type <- match.arg(erlang_type, c("none", "local_kill", "global_kill", "local_guess", "local_kill_guess"))
  type_suffix <- switch(erlang_type,
    "none" = "",
    "local_kill" = "_LOCAL_KILL",
    "global_kill" = "_GLOBAL_KILL",
    "local_guess" = "_LOCAL_GUESS",
    "local_kill_guess" = "_LOCAL_KILL_GUESS"
  )
  list(
    type = "RACE",
    c_name = paste0(
      if (erlang_shape >= 2L) "RDMGBM_E2" else "RDMGBM",
      type_suffix
    ),
    p_types = c(
      "v" = log(1), "B" = log(1), "A" = log(0), "t0" = log(0),
      "s" = log(1), "lambda_g" = log(0), "lambda_k" = log(0), "pContaminant" = qnorm(0)
    ),
    p_types_canonical = c("v", "B", "A", "t0", "s"),
    transform = list(func = c(
      v = "exp", B = "exp", A = "exp", t0 = "exp",
      s = "exp", lambda_g = "exp", lambda_k = "exp", pContaminant = "pnorm"
    )),
    bound = list(
      minmax = cbind(
        v = c(1e-3, Inf), B = c(0, Inf), A = c(0, Inf),
        t0 = c(0.05, Inf), s = c(0, Inf), lambda_g = c(1e-4, Inf), lambda_k = c(1e-4, Inf),
        pContaminant = c(0.001, 0.999)
      ),
      exception = c(A = 0, v = 0, lambda_g = 0, lambda_k = 0, pContaminant = 0),
      # Joint validity: positive log-space drift required for finite-time simulation.
      joint_ok = function(pars) {
        if (!all(c("v", "s") %in% colnames(pars))) return(rep(TRUE, nrow(pars)))
        is.finite(pars[, "v"]) & is.finite(pars[, "s"]) &
          pars[, "s"] > 0 & (pars[, "v"] > 0.5 * pars[, "s"]^2)
      }
    ),
    Ttransform = function(pars, dadm) {
      pars <- cbind(pars, b = 1 + pars[, "B"] + pars[, "A"])
      pars
    },
    rfun = function(data = NULL, pars) {
      rRDMGBM(data$lR, pars, ok = attr(pars, "ok"), erlang_shape = erlang_shape, erlang_type = erlang_type)
    },
    dfun = function(rt, pars) dRDMGBM(rt, pars, erlang = erlang_shape),
    pfun = function(rt, pars) pRDMGBM(rt, pars, erlang = erlang_shape),
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
#' @details
#'
#' Racing Diffusion Model with Shifted Wald Truncated Normal (SWTN) accumulators.
#' Supports between-trial drift variability (sv) and start-point variability (A).
#' When sv=0 and A=0 the model reduces to a point Wald; when sv=0 it reduces
#' to the standard RDM.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Interpretation**                        |
#' |-----------|-----------|---------------|---------|---------------------------------------|
#' | *v*       | log       | \[0, Inf\]      | log(1)  | Mean drift rate                        |
#' | *B*       | log       | \[0, Inf\]      | log(1)  | Response threshold (before A offset)   |
#' | *A*       | log       | \[0, Inf\]      | log(0)  | Start-point variability range          |
#' | *t0*      | log       | \[0, Inf\]      | log(0)  | Non-decision time                      |
#' | *s*       | log       | \[0, Inf\]      | log(1)  | Within-trial SD of drift rate          |
#' | *sv*      | log       | \[0, Inf\]      | log(0)  | Between-trial SD of drift rate         |
#' | *lambda_g*| log       | \[0, Inf\]      | log(0)  | Exponential guessing rate              |
#' | *lambda_k*| log       | \[0, Inf\]      | log(0)  | Exponential killing rate               |
#'
#' @param erlang_shape integer shape of the killing process (1=exponential, 2=Erlang-2)
#' @param erlang_type string, one of "none", "local_kill", "global_kill", "local_guess", "local_kill_guess"
#'
#' @return a list of parameters
#'
#' @export
#'
RDMSWTN <- function(erlang_shape = 1L, erlang_type = "none", posdrift = TRUE) {
  erlang_type <- match.arg(erlang_type, c("none", "local_kill", "global_kill", "local_guess", "local_kill_guess"))
  type_suffix <- switch(erlang_type,
    "none" = "",
    "local_kill" = "_LOCAL_KILL",
    "global_kill" = "_GLOBAL_KILL",
    "local_guess" = "_LOCAL_GUESS",
    "local_kill_guess" = "_LOCAL_KILL_GUESS"
  )
  base_name <- paste0(if (erlang_shape >= 2L) "RDMSWTN_E2" else "RDMSWTN", type_suffix)
  list(
    type = "RACE",
    c_name = if (posdrift) base_name else paste0(base_name, "_IO"),
    p_types = c(
      "v" = 1, "B" = log(1), "A" = log(0), "t0" = log(0),
      "s" = log(1), "sv" = log(0), "lambda_g" = log(0), "lambda_k" = log(0), "pContaminant" = qnorm(0)
    ),
    p_types_canonical = c("v", "B", "A", "t0", "s", "sv"),
    transform = list(func = c(
      v = "identity", B = "exp", A = "exp", t0 = "exp",
      s = "exp", sv = "exp", lambda_g = "exp", lambda_k = "exp", pContaminant = "pnorm"
    )),
    bound = list(
      minmax = cbind(
        v = c(-Inf, Inf), B = c(0, Inf), A = c(0, Inf),
        t0 = c(0.05, Inf), s = c(0, Inf), sv = c(0, Inf),
        lambda_g = c(1e-4, Inf), lambda_k = c(1e-4, Inf),
        pContaminant = c(0.001, 0.999)
      ),
      exception = c(A = 0, v = 0, sv = 0, lambda_g = 0, lambda_k = 0, pContaminant = 0)
    ),
    Ttransform = function(pars, dadm) {
      pars <- cbind(pars, b = pars[, "B"] + pars[, "A"])
      pars
    },
    rfun = function(data = NULL, pars) rRDMSWTN(data$lR, pars, ok = attr(pars, "ok"),
                                                erlang_shape = erlang_shape, erlang_type = erlang_type,
                                                posdrift = posdrift),
    dfun = function(rt, pars) dRDMSWTN(rt, pars, erlang = erlang_shape, posdrift = posdrift),
    pfun = function(rt, pars) pRDMSWTN(rt, pars, erlang = erlang_shape, posdrift = posdrift),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      log_likelihood_race_missing(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
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
    if (any(dimnames(pars)[[2]] == "s")) {
      pars_ok <- pars[ok, , drop = FALSE]
      pars_ok[, c("A", "b", "v", "sv")] <- pars_ok[, c("A", "b", "v", "sv")] / pars_ok[, "s"]
      pars[ok, ] <- pars_ok
    }
    out[ok] <- dSWTNspv(rt[ok],
      v = pars[ok, "v", drop = FALSE], b = pars[ok, "b", drop = FALSE],
      A = pars[ok, "A", drop = FALSE],
      s = if ("s" %in% colnames(pars)) pars[ok, "s", drop = FALSE] else 1,
      t0 = pars[ok, "t0", drop = FALSE],
      sv = pars[ok, "sv", drop = FALSE], lambda_g = pars[ok, "lambda_g", drop = FALSE],
      lambda_k = pars[ok, "lambda_k", drop = FALSE],
      kill_shape = erlang, posdrift = posdrift
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
    if (any(dimnames(pars)[[2]] == "s")) {
      pars_ok <- pars[ok, , drop = FALSE]
      pars_ok[, c("A", "b", "v", "sv")] <- pars_ok[, c("A", "b", "v", "sv")] / pars_ok[, "s"]
      pars[ok, ] <- pars_ok
    }
    out[ok] <- pSWTNspv(rt[ok],
      v = pars[ok, "v", drop = FALSE], b = pars[ok, "b", drop = FALSE],
      A = pars[ok, "A", drop = FALSE],
      s = if ("s" %in% colnames(pars)) pars[ok, "s", drop = FALSE] else 1,
      t0 = pars[ok, "t0", drop = FALSE],
      sv = pars[ok, "sv", drop = FALSE], lambda_g = pars[ok, "lambda_g", drop = FALSE],
      lambda_k = pars[ok, "lambda_k", drop = FALSE],
      kill_shape = erlang, posdrift = posdrift
    )
  }
  out
}

rRDMSWTN <- function(lR, pars, p_types = c("v", "b", "A", "t0", "sv", "lambda_g", "lambda_k"),
                     ok = rep(TRUE, dim(pars)[1]), erlang_shape = 1L, erlang_type = "none",
                     posdrift = TRUE) {
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
  if (!all(p_types %in% dimnames(pars)[[2]])) {
    stop("pars must have columns ", paste(p_types, collapse = " "))
  }
  if (any(dimnames(pars)[[2]] == "s")) {
    pars[, c("A", "b", "v", "sv")] <- pars[, c("A", "b", "v", "sv")] / pars[, "s"]
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
  # In the combined local kill+guess case, handle both clocks jointly below.
  k_vec <- if (local_kill && !guess) pars[, "lambda_k"] else rep(0, nrow(pars))
  dt[ok] <- rSWTN(sum(ok),
    b = pars[, "b"], v = pars[, "v"], A = pars[, "A"], sv = pars[, "sv"],
    k = k_vec, erlang = erlang_shape, posdrift = posdrift
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
        tg_local[guess_ok_mat] <- rgamma(sum(guess_ok_mat), shape = erlang_shape,
                                         rate = lambda_g_local[guess_ok_mat])
      }
    }

    if (local_kill) {
      lambda_k_local <- matrix(0, nrow = nr, ncol = n_trials)
      lambda_k_local[matrix(ok, nrow = nr)] <- pars[, "lambda_k"]
      kill_ok_mat <- lambda_k_local > 0
      if (any(kill_ok_mat)) {
        tk_local[kill_ok_mat] <- rgamma(sum(kill_ok_mat), shape = erlang_shape,
                                        rate = lambda_k_local[kill_ok_mat])
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
