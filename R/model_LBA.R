.lba_dfun <- function (rt, pars, posdrift = TRUE)
  # posdrift = truncated positive normal rates
{
  dt <- rt - pars[,"t0"]
  ok <- (dt>0) & (pars[,"b"] >= pars[,"A"])
  ok[is.na(ok) | !is.finite(dt)] <- FALSE
  out <- numeric(length(dt))
  # dlba is the k = 0 BAwL member evaluated with the LBA normalizer floor,
  # so this R path matches the C++ likelihood kernels exactly.
  out[ok] <- dlba(t = dt[ok], A = pars[ok,"A"], b = pars[ok,"b"],
                  v = pars[ok,"v"], sv = pars[ok,"sv"],
                  posdrift = posdrift)
  out
}

.lba_pfun <- function (rt, pars, posdrift = TRUE)
  # posdrift = truncated positive normal rates
{
  dt <- rt - pars[,"t0"]
  ok <- (dt>0) & (pars[,"b"] >= pars[,"A"])
  ok[is.na(ok) | !is.finite(dt)] <- FALSE
  out <- numeric(length(dt))
  out[ok] <- plba(t = dt[ok], A = pars[ok,"A"], b = pars[ok,"b"],
                  v = pars[ok,"v"], sv = pars[ok,"sv"],
                  posdrift = posdrift)
  is_inf <- is.infinite(rt) & rt > 0 & (pars[,"b"] >= pars[,"A"])
  is_inf[is.na(is_inf)] <- FALSE
  if (any(is_inf)) {
    out[is_inf] <- if (posdrift) 1 else
      pnorm(0, mean = pars[is_inf,"v"], sd = pars[is_inf,"sv"], lower.tail = FALSE)
  }
  out
}


.lba_rfun <- function(lR,pars,p_types=c("v","sv","b","A","t0"),
                 ok=rep(TRUE,length(lR)),posdrift = TRUE)
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in
  # contiguous rows.
{
  bad <- rep(NA, length(lR)/length(levels(lR)))
  out <- data.frame(R = bad, rt = bad)
  nr <- length(levels(lR))
  dt <- matrix(Inf,nrow=nr,ncol=nrow(pars)/nr)
  t0 <- pars[,"t0"]
  pars <- pars[ok,]
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  dt[ok] <- (pars[,"b"]-pars[,"A"]*runif(dim(pars)[1]))/
    msm::rtnorm(dim(pars)[1],pars[,"v"],pars[,"sv"],ifelse(posdrift,0,-Inf))
  dt[dt<0] <- Inf
  bad <- colSums(is.infinite(dt)) == nrow(dt)
  R <- max.col(-t(dt), ties.method='first')
  pick <- cbind(R,1:dim(dt)[2])
  rt <- matrix(t0,nrow=nr)[pick] + dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  ok <- matrix(ok,nrow=length(levels(lR)))[1,]
  out$R[ok] <- levels(lR)[R][ok]
  out$R <- factor(out$R,levels=levels(lR))
  out$rt[ok] <- rt[ok]
  out$R[ok & bad] <- NA
  out <- .apply_timed_guess_winner(out, levels(lR))
  out
}

# Backwards compatibility for emc objects saved before the C++ rfun wrapper
# was introduced. Their serialized model closures call rLBA() directly;
# keep this internal name so loading such a fit uses the current simulator.
rLBA <- function(lR, pars, p_types = c("v", "sv", "b", "A", "t0"),
                 ok = rep(TRUE, length(lR)), posdrift = TRUE) {
  .rfun_LBA(lR, pars, ok = ok, posdrift = posdrift)
}

#### Model functions ----

#' The Linear Ballistic Accumulator model
#'
#' Model file to estimate the Linear Ballistic Accumulator (LBA) in EMC2.
#'
#' Model files are almost exclusively used in `design()`.
#'
#' @details
#'
#' Default values are used for all parameters that are not explicitly listed in the `formula`
#' argument of `design()`.They can also be accessed with `LBA()$p_types`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *v*       | identity  | \[-Inf, Inf\] | 1         |                            | Mean evidence-accumulation rate                                              |
#' | *A*       | log       | \[0, Inf\]    | log(0)    |                            | Between-trial variation (range) in start point                     |
#' | *B*       | log       | \[0, Inf\]    | log(1)    | *b* = *B*+*A*              | Distance from *A* to *b* (response threshold)                                       |
#' | *t0*      | log       | \[0, Inf\]    | log(0)    |                            | Non-decision time                                         |
#' | *sv*      | log       | \[0, Inf\]    | log(1)    |                            | Between-trial variation in evidence-accumulation rate                      |
#'
#'
#' All core LBA parameters are estimated on the log scale, except for the drift
#' rate which is estimated on the real line. Optional fitting parameters
#' `pContaminant` and `pGuess` are estimated on the probit scale; they are the
#' omission and uniform-outlier probabilities, respectively.
#'
#' Conventionally, `sv` is fixed to 1 to satisfy scaling constraints.
#'
#' The *b* = *B* + *A* parameterization ensures that the response threshold is always higher than the between trial variation in start point of the drift rate.
#'
#' Because the LBA is a race model, it has one accumulator per response option.
#' EMC2 automatically constructs a factor representing the accumulators `lR` (i.e., the
#' latent response) with level names taken from the `R` column in the data.
#'
#' The `lR` factor is mainly used to allow for response bias, analogous to `Z` in the
#' DDM. For example, in the LBA, response thresholds are determined by the *B*
#' parameters, so `B~lR` allows for different thresholds for the accumulator
#' corresponding to left and right stimuli (e.g., a bias to respond left occurs
#' if the left threshold is less than the right threshold).
#' For race models, the `design()` argument `matchfun` can be provided, a
#' function that takes the `lR` factor (defined in the augmented data (d)
#' in the following function) and returns a logical defining the correct response.
#' In the example below, the match is simply such that the `S` factor equals the
#' latent response factor: `matchfun=function(d)d$S==d$lR`. Then `matchfun` is
#' used to automatically create a latent match (`lM`) factor with
#' levels `FALSE` (i.e., the stimulus does not match the accumulator) and `TRUE`
#' (i.e., the stimulus does match the accumulator). This is added internally
#' and can also be used in model formula, typically for parameters related to
#' the rate of accumulation.
#'
#' Brown, S. D., & Heathcote, A. (2008). The simplest complete model of choice response time: Linear ballistic accumulation.
#' *Cognitive Psychology, 57*(3), 153-178. https://doi.org/10.1016/j.cogpsych.2007.12.002
#'
#' @param posdrift Logical. If TRUE (default), drift rates are truncated to be positive.
#' @return A model list with all the necessary functions for EMC2 to sample
#' @examples
#' # When working with lM it is useful to design  an "average and difference"
#' # contrast matrix, which for binary responses has a simple canonical from:

#' ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
#' # We also define a match function for lM
#' matchfun=function(d)d$S==d$lR
#' # We now construct our design, with v ~ lM and the contrast for lM the ADmat.
#' design_LBABE <- design(data = forstmann,model=LBA,matchfun=matchfun,
#'                        formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
#'                        contrasts=list(v=list(lM=ADmat)),constants=c(sv=log(1)))
#' # For all parameters that are not defined in the formula, default values are assumed
#' # (see Table above).
#' @export
#'

LBA <- function(posdrift=TRUE){
  list(
    type="RACE",
    c_name = ifelse(posdrift,"LBA","LBAIO"),
    # The C++ likelihood adapter evaluates LBA through the shared BAwL
    # kernels with k = 0 and clock means mG = mK = 0 on the natural scale.
    # Those optional BAwL parameters are deliberately not part of the LBA
    # parameter vector; if represented on the transformed scale, their off
    # value would be log(0) = -Inf.
    # p_vector transform, sets sv as a scaling parameter
    p_types=c("v" = 1,"sv" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0), "pContaminant"=qnorm(0), "pGuess"=qnorm(0)),
    p_types_canonical = c("v", "sv", "B", "A", "t0"),
    transform=list(func=c(v = "identity",sv = "exp", B = "exp", A = "exp",t0 = "exp",pContaminant="pnorm",pGuess="pnorm")),
    bound=list(minmax=cbind(v=c(-Inf,Inf),sv = c(1e-4, Inf), A=c(1e-4,Inf),B=c(1e-4,Inf),t0=c(0.05,Inf),pContaminant=c(0.001,0.999),pGuess=c(0.001,0.999)),
               exception=c(A=0,pContaminant=0,pGuess=0)),
    # Transform to natural scale
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
      pars
    },
    # Random function for racing accumulator
    rfun=ifelse(posdrift,function(data,pars) .rfun_LBA(data$lR,pars,ok = attr(pars, "ok"),posdrift=TRUE),
                function(data,pars) .rfun_LBA(data$lR,pars,ok = attr(pars, "ok"),posdrift=FALSE)),
    # Density function (PDF) for single accumulator
    dfun=ifelse(posdrift,function(rt,pars) .lba_dfun(rt,pars,posdrift=TRUE),
                function(rt,pars) .lba_dfun(rt,pars,posdrift=FALSE)),
    # Probability function (CDF) for single accumulator
    pfun=ifelse(posdrift,function(rt,pars) .lba_pfun(rt,pars,posdrift=TRUE),
                function(rt,pars) .lba_pfun(rt,pars,posdrift=FALSE)),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
      log_likelihood_race_missing(pars=pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}

#' LBA Logical Rules Model
#'
#' This is an LBA-based model for tasks in which the observable response is
#' assembled from two latent target-detection subraces. It is not an ordinary
#' race over the levels of `R`: the compiled likelihood combines the finish
#' times of the accumulator roles according to the `LogicalRule` column in the
#' data.
#'
#' The ordinary logical-rules model uses the accumulator roles `A`, `B`,
#' `n_A`, and `n_B`. `A` and `B` are target accumulators, while `n_A` and
#' `n_B` are their corresponding nontarget accumulators. The supported choice
#' rules are:
#'
#' * `OR`: respond `yes` when either target subrace is positive, otherwise
#'   respond `no`.
#' * `AND`: respond `yes` only when both target subraces are positive.
#' * `XOR`: respond `yes` when exactly one target subrace is positive.
#' * `ID`: report which targets are positive using `NN`, `AN`, `NB`, or `AB`.
#'
#' Two further rules target detection paradigms. `OR_DETECTION_ANALYTIC` uses
#' only the active target accumulators `A` and `B`, and reports the first
#' active detector to finish; its stimulus column (`S`, `stimulus`, or
#' `condition`) must identify `NN`, `AN`, `NB`, or `AB` on every trial (missing
#' values are treated as `NN`). `OR_DETECTION_GNG` is the go/no-go version of
#' the ordinary `OR` task: it uses the same four `A`, `n_A`, `B`, `n_B`
#' subraces, an overt "go" response is the first target subrace to win (coded
#' `yes`), and the "no" outcome — both subraces resolving to their nontarget —
#' is instead a withheld response coded as `rt = Inf` (with a missing `R`).
#' Because it is a four-horse race like `OR`, the drifts are driven by the
#' design (not by gating detectors on the stimulus), so false alarms remain
#' possible on `NN` trials.
#'
#' The LBA parameters are:
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *v* | identity | \[-Inf, Inf\] | 1 | | Mean drift rate for each accumulator. |
#' | *sv* | log | \[0, Inf\] | log(1) | | Between-trial SD of drift rate; conventionally fixed to 1 for scale identification. |
#' | *B* | log | \[0, Inf\] | log(1) | *b* = *B* + *A* | Distance from the upper end of the start-point range to the threshold. |
#' | *A* | log | \[0, Inf\] | log(0) | | Width of the uniform start-point range. |
#' | *t0* | log | \[0, Inf\] | log(0) | | Non-decision time. |
#'
#' Thus each accumulator starts at `A * U`, with `U ~ Uniform(0, 1)`, and
#' reaches the threshold `b = B + A`. With `posdrift = TRUE`, the drift draw
#' is a positive normal draw; with `posdrift = FALSE`, the untruncated normal
#' is used and negative-drift accumulators can produce intrinsic omissions.
#' All of the parameters above are trial-dependent after the design formulas
#' have been evaluated. The optional capacity extension adds `kappa` on the
#' identity scale (default `0`) and `tau` on the log scale (default `log(0)`).
#' They are active when `capacity = TRUE`.
#'
#' With `capacity = TRUE`, an `AB` trial receives one latent `Z ~ N(0, 1)` and
#' the target drift draws are
#' `V_i = v_i + kappa + tau * Z + epsilon_i`, for `i = A, B`, with
#' independent `epsilon_i ~ N(0, sv_i^2)`. The logical-rule calculation is
#' performed conditional on this same `Z` and then integrated over `Z`; the
#' factor is not integrated separately for the two subraces. Nontarget
#' accumulators do not load on the factor. Capacity has no effect on
#' `AN`, `NB`, or `NN` trials, and `kappa = 0, tau = 0` recovers the ordinary
#' logical-rules likelihood exactly. `kappa` and `tau` must be shared by the
#' `A` and `B` rows within a trial. In positive-drift mode, the active target
#' draws are conditioned jointly to be positive.
#' The optional fitting parameter `pContaminant` is the omission probability.
#'
#' The likelihood is implemented in the compiled logical-rules fast path. The
#' R `dfun` and `pfun` entries retain the single-accumulator LBA functions for
#' model integration and diagnostics, but the model-level R likelihood is not
#' available as a fallback.
#'
#' @param posdrift Logical. If `TRUE` (default), use positive-truncated LBA
#'   drift rates; if `FALSE`, use untruncated normal drift rates and append
#'   `IO` to the compiled model name.
#' @param fast_path Logical argument retained for compatibility. The current
#'   logical-rules likelihood is always dispatched through the compiled path.
#' @param capacity Logical. If `TRUE`, include the `kappa` and `tau` parameters
#'   and use the shared-capacity likelihood on redundant-target `AB` trials.
#' @return A model list defining the logical-rules race model.
#' @export
#'
LogicalRulesLBA <- function(posdrift = TRUE, fast_path=TRUE, capacity = FALSE){
  p_types <- c("v" = 1,"sv" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0), "pContaminant"=qnorm(0))
  transform <- c(v = "identity",sv = "exp", B = "exp", A = "exp",t0 = "exp", pContaminant="pnorm")
  minmax <- cbind(v=c(-Inf,Inf),sv = c(1e-4, Inf), A=c(1e-4,Inf),B=c(0,Inf),t0=c(0.05,Inf), pContaminant=c(0.001,0.999))
  exception <- c(A=0, pContaminant=0)
  if (capacity) {
    # kappa/tau are appended after the emc2col::lba kernel prefix so the raw
    # batch kernels keep their positional contract.  tau = 0 (from the
    # default log(0)) and kappa = 0 (identity scale) are exactly representable, so
    # the ordinary route is recovered bit-for-bit when capacity is off.
    p_types <- c(p_types, kappa = 0, tau = log(0))
    transform <- c(transform, kappa = "identity", tau = "exp")
    minmax <- cbind(minmax, kappa = c(-Inf, Inf), tau = c(1e-4, Inf))
    exception <- c(exception, tau = 0)
  }
  list(
    type="RACE",
    # Note: `calc_ll_oo()` infers LBA `posdrift` from whether `c_name` contains "IO".
    # Keep this in sync so likelihood and simulation use the same setting.
    # The compiled likelihood detects the capacity variant from the presence
    # of the kappa/tau parameter columns, so the c_name is shared.
    c_name = paste0("LBA_LogicalRules",ifelse(posdrift,"","IO")),
    capacity = capacity,
    # p_vector transform, sets sv as a scaling parameter
    p_types=p_types,
    transform=list(func=transform),
    bound=list(minmax=minmax,
               exception=exception),
    # Transform to natural scale
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
      pars
    },
    # Random function for racing accumulator
    rfun=function(data,pars) .rfun_LBA(data$lR,pars,ok = attr(pars, "ok"),posdrift=ifelse(posdrift,TRUE,FALSE)),
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) .lba_dfun(rt,pars,posdrift=ifelse(posdrift,TRUE,FALSE)),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) .lba_pfun(rt,pars,posdrift=ifelse(posdrift,TRUE,FALSE)),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
      stop("LogicalRulesLBA: R likelihood path not implemented. Use fast_path=TRUE (the default).")
    }
  )
}

#### BAwL (Ballistic Accumulator with Leak) ----

# 0 = normal launch strength (v, sv); 1 = lognormal launch strength
# (mu, sigma); 2 = median-parameterised continuous split-lognormal launch
# (mu, sigma, delta); 3 = Weibull launch (shape, scale). Must match
# BAWL_LAUNCH_* in src/model_LBA.h, the
# `launch` argument of the compiled kernels, and the value the adapter derives
# from the `_LOGN`/`_SPLIT` c_name suffix.  Deliberately the same convention as
# BAwD.

# The launch pair occupies the same leading kernel columns; split-lognormal
# adds delta immediately after (mu, sigma).

.bawl_check_cols <- function(pars, launch) {
  need <- .ba_par_names(launch)
  missing <- setdiff(need, colnames(pars))
  if (length(missing))
    stop("BAwL requires parameter columns ", paste(missing, collapse = ", "))
  need
}

dBAwL <- function(rt, pars, posdrift = TRUE, erlang = 1L, guess = FALSE,
                  launch = 0L) {
  if (!all(c("lambda_g", "lambda_k") %in% colnames(pars))) {
    stop("BAwL requires parameter columns 'lambda_g' and 'lambda_k'.")
  }
  nm  <- .bawl_check_cols(pars, launch)
  dt  <- rt - pars[, "t0"]
  erl <- (pars[, "lambda_g"] > 0) | (pars[, "lambda_k"] > 0)
  ok  <- (rt > 0) & ((dt > 0) | erl) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok) | !is.finite(dt)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- dkilledleakyba(
      t = rt[ok], v = pars[ok, nm[1]], b = pars[ok, "b"], A = pars[ok, "A"],
      sv = pars[ok, nm[2]], t0 = pars[ok, "t0"], k = pars[ok, "k"],
      lambda_g = pars[ok, "lambda_g"], lambda_k = pars[ok, "lambda_k"],
      posdrift = posdrift, log_out = FALSE,
      kill_shape = as.integer(erlang), guess = guess,
      erlang_omega = .rdmswtn_erlang_omega(pars[ok, , drop = FALSE], erlang),
      launch = as.integer(launch),
      delta = if (launch == 2L) pars[ok, "delta"] else 0
    )
  }
  out
}

pBAwL <- function(rt, pars, posdrift = TRUE, erlang = 1L, guess = FALSE,
                  launch = 0L) {
  if (!all(c("lambda_g", "lambda_k") %in% colnames(pars))) {
    stop("BAwL requires parameter columns 'lambda_g' and 'lambda_k'.")
  }
  nm  <- .bawl_check_cols(pars, launch)
  dt  <- rt - pars[, "t0"]
  erl <- (pars[, "lambda_g"] > 0) | (pars[, "lambda_k"] > 0)
  ok  <- (rt > 0) & ((dt > 0) | erl) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- pkilledleakyba(
      t = rt[ok], v = pars[ok, nm[1]], b = pars[ok, "b"], A = pars[ok, "A"],
      sv = pars[ok, nm[2]], t0 = pars[ok, "t0"], k = pars[ok, "k"],
      lambda_g = pars[ok, "lambda_g"], lambda_k = pars[ok, "lambda_k"],
      posdrift = posdrift, log_out = FALSE,
      kill_shape = as.integer(erlang), guess = guess,
      erlang_omega = .rdmswtn_erlang_omega(pars[ok, , drop = FALSE], erlang),
      launch = as.integer(launch),
      delta = if (launch == 2L) pars[ok, "delta"] else 0
    )
  }
  out
}

rBAwL <- function(lR, pars, ok = rep(TRUE, length(lR)),
                  p_types = NULL,
                  posdrift = TRUE, eps = 1e-10, erlang = 1L, guess = FALSE, global = FALSE,
                  .drifts = NULL, launch = 0L) {
  if (!all(c("lambda_g", "lambda_k") %in% colnames(pars))) {
    stop("BAwL requires parameter columns 'lambda_g' and 'lambda_k'.")
  }
  nm <- .bawl_check_cols(pars, launch)
  if (is.null(p_types))
    p_types <- c(nm, "b", "A", "t0", "k", "lambda_g", "lambda_k")
  erlang_omega_all <- .rdmswtn_erlang_omega(pars, erlang)
  bad  <- rep(NA, length(lR) / length(levels(lR)))
  out  <- data.frame(R = bad, rt = bad)
  nr   <- length(levels(lR))
  n_trials <- nrow(pars) / nr
  dt   <- matrix(Inf, nrow = nr, ncol = n_trials)
  t0   <- pars[, "t0"]

  if (global) {
    lambda_mat <- matrix(pars[, "lambda_k"], nrow = nr)
    if (any(apply(lambda_mat, 2, function(x) length(unique(x)) > 1))) {
      stop("global=TRUE requires lambda_k to be constant across accumulators (lambda_k must not vary by any accumulator-level factor)")
    }
    lambda_trials <- lambda_mat[1, ]
    tk_global <- rep(Inf, n_trials)
    kill_ok <- !is.na(lambda_trials) & lambda_trials > 0
    if (any(kill_ok)) {
      shape_k <- if (as.integer(erlang) == 3L) {
        ifelse(runif(sum(kill_ok)) <= matrix(erlang_omega_all, nrow = nr)[1, kill_ok], 1L, 2L)
      } else {
        as.integer(erlang)
      }
      rate_k <- if (as.integer(erlang) == 3L) {
        ifelse(shape_k == 2L, 2 * lambda_trials[kill_ok], lambda_trials[kill_ok])
      } else {
        lambda_trials[kill_ok]
      }
      tk_global[kill_ok] <- rgamma(sum(kill_ok), shape = shape_k, rate = rate_k)
    }
  }

  pars_all <- pars
  pars <- pars[ok, , drop = FALSE]
  ok_idx <- which(ok)
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ", paste(p_types, collapse = " "))
  lower  <- if (posdrift) 0 else -Inf
  if (is.null(.drifts)) {
    # A lognormal launch strength is positive by construction, so posdrift
    # never applies to it.
    drifts <- if (launch == 3L) {
      rweibull(nrow(pars), pars[, nm[1]], pars[, nm[2]])
    } else if (launch == 1L) {
      rlnorm(nrow(pars), pars[, nm[1]], pars[, nm[2]])
    } else if (launch == 2L) {
      .bawd_split_rlnorm(pars[, "mu"], pars[, "sigma"], pars[, "delta"])
    } else {
      msm::rtnorm(nrow(pars), mean = pars[, nm[1]], sd = pars[, nm[2]], lower = lower)
    }
  } else {
    if (length(.drifts) != nrow(pars_all))
      stop(".drifts must have one value per row of the original parameter matrix.")
    drifts <- .drifts[ok]
  }

  small_k <- pars[, "k"] < eps
  if (any(small_k))
    dt[ok_idx[small_k]] <- (pars[small_k, "b"] - pars[small_k, "A"] * runif(sum(small_k))) /
                           drifts[small_k]

  big_k <- !small_k & (drifts > pars[, "k"] * pars[, "b"])
  if (any(big_k)) {
    num        <- drifts[big_k] - pars[big_k, "k"] * pars[big_k, "b"]
    den        <- drifts[big_k] - pars[big_k, "k"] * pars[big_k, "A"] * runif(sum(big_k))
    ratio      <- num / den
    ratio      <- pmin(pmax(ratio, .Machine$double.xmin), 1 - 1e-15)
    dt[ok_idx[big_k]] <- (-1 / pars[big_k, "k"]) * log(ratio)
  }
  dt[dt < 0] <- Inf
  # Put EAM on the same raw-time axis as Erlang clocks.
  dt <- dt + matrix(t0, nrow = nr)

  if (guess || !global) {
    tg_local <- matrix(Inf, nrow = nr, ncol = n_trials)
    tk_local <- matrix(Inf, nrow = nr, ncol = n_trials)

    if (!global) {
      lambda_k_local <- matrix(0, nrow = nr, ncol = n_trials)
      lambda_k_local[matrix(ok, nrow = nr)] <- pars[, "lambda_k"]
      kill_ok_mat <- lambda_k_local > 0
      if (any(kill_ok_mat)) {
        shape_k <- if (as.integer(erlang) == 3L) {
          ifelse(runif(sum(kill_ok_mat)) <= matrix(erlang_omega_all, nrow = nr)[kill_ok_mat], 1L, 2L)
        } else {
          as.integer(erlang)
        }
        rate_k <- if (as.integer(erlang) == 3L) {
          ifelse(shape_k == 2L, 2 * lambda_k_local[kill_ok_mat], lambda_k_local[kill_ok_mat])
        } else {
          lambda_k_local[kill_ok_mat]
        }
        tk_local[kill_ok_mat] <- rgamma(sum(kill_ok_mat), shape = shape_k, rate = rate_k)
      }
    }

    lambda_g_local <- matrix(0, nrow = nr, ncol = n_trials)
    active_guess <- matrix(ok, nrow = nr) & (levels(lR) != "nogo")
    if (guess) {
      lambda_g_local[active_guess] <- pars_all[, "lambda_g"][active_guess]
      guess_ok_mat <- lambda_g_local > 0
      if (any(guess_ok_mat)) {
        shape_g <- if (as.integer(erlang) == 3L) {
          ifelse(runif(sum(guess_ok_mat)) <= matrix(erlang_omega_all, nrow = nr)[guess_ok_mat], 1L, 2L)
        } else {
          as.integer(erlang)
        }
        rate_g <- if (as.integer(erlang) == 3L) {
          ifelse(shape_g == 2L, 2 * lambda_g_local[guess_ok_mat], lambda_g_local[guess_ok_mat])
        } else {
          lambda_g_local[guess_ok_mat]
        }
        tg_local[guess_ok_mat] <- rgamma(sum(guess_ok_mat), shape = shape_g, rate = rate_g)
      }
    }

    if (guess && !global) {
      dt_candidate <- pmin(dt, tg_local)
      dt <- ifelse(dt_candidate < tk_local, dt_candidate, Inf)
    } else if (guess) {
      dt <- pmin(dt, tg_local)
    } else {
      dt <- ifelse(dt < tk_local, dt, Inf)
    }
  }

  if (global) {
    is_killed <- tk_global < apply(dt, 2, min)
    if (any(is_killed)) dt[, is_killed] <- Inf
  }
  bad_col <- colSums(!is.infinite(dt)) == 0L
  R   <- apply(dt, 2, which.min)
  pick <- cbind(R, 1:dim(dt)[2])
  rt   <- dt[pick]
  R    <- factor(levels(lR)[R], levels = levels(lR))
  R[bad_col]  <- NA
  rt[bad_col] <- Inf
  ok  <- matrix(ok, nrow = length(levels(lR)))[1, ]
  out$R[ok]   <- levels(lR)[R][ok]
  out$R       <- factor(out$R, levels = levels(lR))
  out$rt[ok]  <- rt[ok]
  out <- .apply_timed_guess_winner(out, levels(lR))
  out
}

# One-factor correlated BAwL reference simulator.  The low-level `pars` input
# carries the signed row-level factor-variance share consumed by the scalar
# kernels.  BAwLcorr's Ttransform maps a cell-level rho onto these row values
# when the expanded data contain the lM role indicator.
rBAwL_corr <- function(lR, pars, ok = rep(TRUE, nrow(pars)),
                       posdrift = TRUE, eps = 1e-10, erlang = 1L,
                       guess = FALSE, global = FALSE) {
  if (!"rho" %in% colnames(pars))
    stop("Correlated BAwL requires parameter column 'rho'.")
  if (any(!is.finite(pars[, "rho"]) | abs(pars[, "rho"]) > 1))
    stop("Correlated BAwL rho values must be finite and lie in [-1, 1].")
  nr <- length(levels(lR))
  if (nr <= 0 || nrow(pars) %% nr != 0)
    stop("Correlated BAwL requires rows grouped by accumulator within trial.")
  n_trials <- nrow(pars) / nr
  rho <- pars[, "rho"]
  drifts <- rep(Inf, nrow(pars))

  # For posdrift=TRUE the model is a jointly truncated Gaussian drift vector:
  # draw the shared factor and all active untruncated drifts together, then
  # reject the complete trial if any active drift is non-positive.  This keeps
  # the factor distribution correctly reweighted by the probability that all
  # active drifts are positive.
  max_iter <- 100000L
  for (j in seq_len(n_trials)) {
    rows <- ((j - 1L) * nr + 1L):(j * nr)
    active <- rows[ok[rows]]
    if (!length(active)) next
    accepted <- FALSE
    for (iter in seq_len(max_iter)) {
      z <- rnorm(1)
      magnitude <- abs(rho[active])
      loading <- sqrt(magnitude)
      residual <- pmax(sqrt(pmax(0, 1 - magnitude)), 1e-12)
      mu <- pars[active, "v"] + sign(rho[active]) * pars[active, "sv"] * loading * z
      draw <- rnorm(length(active), mean = mu, sd = pars[active, "sv"] * residual)
      if (!posdrift || all(draw > 0)) {
        drifts[active] <- draw
        accepted <- TRUE
        break
      }
    }
    if (!accepted)
      stop("Correlated BAwL jointly positive drift rejection exceeded ",
           max_iter, " attempts; check that the drift means are not far below zero.")
  }
  rBAwL(lR, pars, ok = ok, posdrift = posdrift, eps = eps,
        erlang = erlang, guess = guess, global = global, .drifts = drifts)
}

#' The Ballistic Accumulator with Leak (BAwL)
#'
#' A race model in which each accumulator follows a leaky evidence trajectory
#' and can race against optional guess and kill clocks. For accumulator `i`,
#' the start point is `A * U_i`, where `U_i ~ Uniform(0, 1)`, the threshold is
#' `b = B + A`, and the drift is a trialwise draw from the launch distribution
#' selected by `drift_distribution`: a normal draw with mean `v_i` and SD
#' `sv_i` (the default), or a lognormal draw with `log V_i ~ N(mu_i, sigma_i^2)`.
#' With leak rate `k`, the evidence trajectory is
#' `x_i(t) = x_i(0) * exp(-k * t) + v_i / k * (1 - exp(-k * t))` for `k > 0`;
#' its `k = 0` limit is exactly the ordinary LBA trajectory with the same
#' launch distribution: positive-truncated normal when `posdrift = TRUE`, or
#' unrestricted normal when `posdrift = FALSE`. If the leaky asymptote cannot
#' reach the threshold, that accumulator has no finite hit. For `k > 0`,
#' finite hits have no common hard right endpoint: launch values approaching
#' the threshold condition from above can produce arbitrarily late responses.
#' This is distinct from BAwD, whose decaying drive plus positive clearance can
#' impose a finite endpoint.
#'
#' The model uses the following parameterization. The `B` parameter is the
#' distance from the upper end of the start-point range to the threshold, so
#' `b = B + A` is always at least `A`. Parameters are mapped to the natural
#' scale before trial-dependent transforms are applied.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *v* | identity | \[-Inf, Inf\] | 1 | | Mean evidence-accumulation rate. |
#' | *sv* | log | \[0, Inf\] | log(1) | | Between-trial SD of the drift rate; conventionally fixed to 1 for scale identification. |
#' | *B* | log | \[0, Inf\] | log(1) | *b* = *B* + *A* | Distance from the upper start-point range to the threshold. |
#' | *A* | log | \[0, Inf\] | log(0) | | Start-point range. |
#' | *t0* | log | \[0, Inf\] | log(0) | | Non-decision time. |
#' | *k* | log | \[0, Inf\] | log(0) | | Leak rate; `k = 0` is the LBA limit. |
#'
#' With `drift_distribution = "lognormal"`, `v` and `sv` are replaced by `mu`
#' (identity, default 0) and `sigma` (log, default `log(1)`), giving
#' `log V ~ N(mu, sigma^2)`.  With `"splitlognormal"`, `mu` is the exact
#' median of `log V`, `delta` is an unbounded real parameter, and the
#' continuous split-normal widths are
#' `sigma_L = sigma * exp(delta/2)` and `sigma_R = sigma * exp(-delta/2)`;
#' the split point is derived from the median condition.  `delta = 0`
#' reduces exactly to the lognormal launch.  Both lognormal launches are
#' positive by construction, so `posdrift` does not apply.  The `correlated`
#' one-factor path — which decomposes the *Gaussian* drift vector — is not
#' available for either lognormal variant.
#' With `drift_distribution = "weibull"`, `V ~ Weibull(shape, scale)` with
#' positive shape and scale. Weibull launches are also incompatible with the
#' Gaussian `correlated` path.
#'
#' **Fixing the evidence scale.** The evidence axis is defined only up to a
#' scale: `(V, b, A) -> (cV, cb, cA)` leaves every crossing time unchanged, so
#' exactly one parameter must be fixed. With `drift_distribution = "normal"`
#' the usual `constants = c(sv = log(1))` does it. With
#' `drift_distribution = "lognormal"` it does **not**: `sigma` is
#' dimensionless, and on the sampled scale the rescaling is the single
#' direction `mu -> mu + log c`, `B -> B + log c`, `A -> A + log c`, leaving
#' `sigma`, `k`, `t0` and every contrast coefficient untouched. Fixing *any
#' one* of `mu`'s intercept, `B`'s intercept, or `A` therefore identifies the
#' scale; `constants = c(mu = 0)` and `constants = c(B = log(1))` are both
#' valid and differ only by reparameterization. Because the rescaling shifts
#' only the intercept, condition effects on the fixed parameter survive: with
#' `B ~ E`, `constants = c(B = log(1))` leaves `B_Eneutral` and `B_Eaccuracy`
#' free and identified. Omitting the constraint produces a ridge in the
#' posterior rather than an error.
#'
#' Here `q = 1` for Erlang-1 clocks and `q = 2` for Erlang-2 clocks. In
#' `erlang_shape = "mixed"`, the clock is Erlang-1 with probability `omega` and
#' Erlang-2 otherwise, with both components having the same mean `mG` or `mK`.
#' The internal rates `lambda_g` and `lambda_k` are generated by the
#' `Ttransform`; they are not user-facing model parameters and should not be
#' put in a `design()` formula.
#' The optional clock means `mG` and `mK` use log/exp transforms with default
#' `log(1)` when their corresponding clocks are active; mixed mode additionally
#' uses optional `omega` on the probit scale with default `qnorm(.5)`. Optional
#' fitting parameters are `pContaminant`, the omission probability,
#' and `pGuess`, the uniform-outlier probability.
#'
#' `erlang_type` selects which clocks are active. `"none"` has no clocks;
#' `"local_kill"` adds an independent kill clock to each accumulator;
#' `"global_kill"` adds one kill clock shared by all accumulators in a trial;
#' `"local_guess"` adds local guess clocks; and `"local_kill_guess"` adds both
#' local guess and local kill clocks. A kill clock that wins produces an
#' omission. A guess clock that wins is represented using the package's timed
#' guess convention; if a `time` accumulator is present, its win is converted
#' to a sampled response and marked with `isTime`. Global kill requires the
#' `mK`/`lambda_k` value to be constant across accumulators within a trial.
#'
#' With `posdrift = TRUE` (the default), the drift draws are truncated to be
#' positive. With `posdrift = FALSE`, the normal draws are untruncated and the
#' model can have intrinsic omissions; the compiled model name is suffixed with
#' `IO`. The BAwL likelihood is a race likelihood over the accumulator levels
#' in `lR`, which EMC2 constructs from the response levels or fixed accumulator
#' roles.
#'
#' If `correlated = TRUE`, `rho` is a trial/cell-level parameter and a
#' `matchfun`-generated `lM` factor is required. The correct racer is the
#' positive reference and the incorrect racer receives the sign of `rho`; this
#' produces the requested pairwise correlation while allowing structurally
#' independent racers to be assigned `rho = 0`. The positive-drift correlated
#' simulator jointly conditions the active drift vector on all active drifts
#' being positive, so `rho` refers to the underlying untruncated Gaussian
#' draws, not to their marginally truncated correlation.
#'
#' @param posdrift Logical. If `TRUE` (default), drift rates are truncated at
#'   zero; if `FALSE`, they are sampled from the untruncated normal.
#' @param erlang_shape Integer `1` for exponential clocks, `2` for Erlang-2,
#'   or `"mixed"` for an Erlang-1/Erlang-2 mixture controlled by `omega`.
#' @param erlang_type Clock configuration: one of `"none"`, `"local_kill"`,
#'   `"global_kill"`, `"local_guess"`, or `"local_kill_guess"`.
#' @param correlated Logical. If `TRUE`, add the correlated-drift parameter
#'   `rho` and use the correlated BAwL race simulator and likelihood path.
#'   Only available for `drift_distribution = "normal"`.
#' @param drift_distribution Distribution of the trialwise launch strength:
#'   `"normal"` (the default) for `V ~ N(v, sv^2)`, `"lognormal"` for
#'   `log V ~ N(mu, sigma^2)`, or `"splitlognormal"` for a continuous
#'   split-normal `log V` with widths
#'   `sigma_L = sigma * exp(delta/2)` and `sigma_R = sigma * exp(-delta/2)`.
#'   In the split variant `mu` is the exact median and `delta` is unbounded;
#'   `delta = 0` is exactly the lognormal launch. `"weibull"` uses
#'   `V ~ Weibull(shape, scale)` on the positive launch scale.
#' @return A model list defining the BAwL race model.
#' @examples
#' # A lognormal-launch BAwL. mu's intercept is fixed to identify the evidence
#' # scale, which sigma no longer does; B and the lM effect stay free.
#' ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
#' matchfun <- function(d) d$S == d$lR
#' design_BAwL_logn <- design(
#'   data = forstmann, model = function() BAwL(drift_distribution = "lognormal"),
#'   matchfun = matchfun,
#'   formula = list(mu ~ lM, sigma ~ 1, B ~ E, A ~ 1, t0 ~ 1, k ~ 1),
#'   contrasts = list(mu = list(lM = ADmat)),
#'   constants = c(mu = 0))
#'
#' @export
BAwL <- function(posdrift = TRUE, erlang_shape = 1L,
                 erlang_type = c("none", "local_kill", "global_kill", "local_guess", "local_kill_guess"),
                 correlated = FALSE,
                 drift_distribution = c("normal", "lognormal", "splitlognormal", "weibull")) {
  erlang_type <- match.arg(erlang_type)
  drift_distribution <- match.arg(drift_distribution)
  launch <- .ba_launch_code(drift_distribution, "BAwL")
  lognormal <- launch %in% c(1L, 2L)
  weibull <- launch == 3L
  splitlognormal <- launch == 2L
  if ((lognormal || weibull) && !isTRUE(posdrift)) {
    stop("BAwL: posdrift only applies to drift_distribution = \"normal\"; a ",
         "lognormal and Weibull launch strengths are positive by construction.")
  }
  if ((lognormal || weibull) && correlated) {
    # The correlated path is a one-factor decomposition of the Gaussian drift
    # vector (drift_factor.h); a lognormal analogue is a different model, not
    # a swapped marginal, so refuse rather than silently ignore one of them.
    stop("BAwL: correlated = TRUE is only implemented for ",
         "drift_distribution = \"normal\".")
  }
  erlang_mixed <- identical(erlang_shape, "mixed")
  erlang_shape_cpp <- if (erlang_mixed) 3L else as.integer(erlang_shape)

  has_guess <- erlang_type %in% c("local_guess", "local_kill_guess")
  has_kill  <- erlang_type %in% c("local_kill", "global_kill", "local_kill_guess")

  # "_SPLIT" follows "_LOGN"; neither lognormal variant is combined with IO.
  base_name <- paste0(ifelse(posdrift, "BAwL", "BAwLIO"),
                    if (weibull) "_WEIB" else if (lognormal) "_LOGN" else "",
                    if (splitlognormal) "_SPLIT" else "",
                    if (erlang_mixed) "_EMIX" else if (erlang_shape_cpp >= 2L) "_E2" else "")
  type_suffix <- if (erlang_type == "local_guess") "_LOCAL_GUESS"
                    else if (erlang_type == "local_kill_guess") "_LOCAL_KILL_GUESS"
                    else if (erlang_type == "local_kill") "_LOCAL_KILL"
                    else if (erlang_type == "global_kill") "_GLOBAL_KILL"
                    else ""

  # The launch pair occupies the leading two kernel columns either way; only
  if (weibull) {
    p_types <- c(shape = log(1), scale = log(1))
    transform <- c(shape = "exp", scale = "exp")
    minmax <- cbind(shape = c(1e-4, Inf), scale = c(1e-4, Inf))
  } else if (lognormal) {
    p_types <- c("mu" = 0, "sigma" = log(1))
    transform <- c(mu = "identity", sigma = "exp")
    minmax <- cbind(mu = c(-Inf, Inf), sigma = c(1e-4, Inf))
    if (splitlognormal) {
      p_types <- c(p_types, delta = 0)
      transform <- c(transform, delta = "identity")
      minmax <- cbind(minmax, delta = c(-Inf, Inf))
    }
  } else {
    p_types <- c("v" = 1, "sv" = log(1))
    transform <- c(v = "identity", sv = "exp")
    minmax <- cbind(v = c(-Inf, Inf), sv = c(1e-4, Inf))
  }
  p_types <- c(p_types, "B" = log(1), "A" = log(0), "t0" = log(0),
               "k" = log(0))
  transform <- c(transform, B = "exp", A = "exp", t0 = "exp", k = "exp")
  minmax <- cbind(minmax, A  = c(1e-4, Inf), B  = c(1e-4, Inf),
                     t0 = c(0.05, Inf), k  = c(1e-4, Inf))
  exception <- c(A = 0, k = 0)
  launch_pars <- .ba_par_names(launch)
  
  p_types <- c(p_types, mG = log(1))
  transform <- c(transform, mG = "exp")
  minmax <- cbind(minmax, mG = c(1e-4, Inf))
  exception <- c(exception, mG = 0)

  p_types <- c(p_types, mK = log(1))
  transform <- c(transform, mK = "exp")
  minmax <- cbind(minmax, mK = c(1e-4, Inf))
  exception <- c(exception, mK = 0)

  if (erlang_mixed) {
    p_types <- c(p_types, omega = qnorm(0.5))
    transform <- c(transform, omega = "pnorm")
    minmax <- cbind(minmax, omega = c(0, 1))
    exception <- c(exception, omega = 0)
  }

  # pContaminant (omission) and pGuess (uniform outlier); see add_nuisance_pars().
  .nuis <- add_nuisance_pars(p_types, transform, minmax, exception)
  p_types <- .nuis$p_types; transform <- .nuis$transform
  minmax <- .nuis$minmax; exception <- .nuis$exception

  if (correlated) {
    # rho is a direct cell-level race correlation.  BAwLcorr's Ttransform
    # derives the signed row-level variance shares from it when lM is present.
    # The endpoints are not reachable, which avoids zero conditional drift SDs
    # in the scalar BAwL kernels.
    p_types <- c(p_types, rho = qnorm(0.5))
    transform <- c(transform, rho = "pnorm")
    minmax <- cbind(minmax, rho = c(-.99, .99))
    exception <- c(exception, rho = 0)
  }

  transform_spec <- list(func = transform)
  if (correlated) {
    transform_spec$lower <- c(rho = -1)
    transform_spec$upper <- c(rho = 1)
  }

  list(
    type   = "RACE",
    c_name = paste0(base_name, type_suffix, if (correlated) "_CORR" else ""),
    correlated = correlated,
    drift_distribution = drift_distribution,
    p_types = p_types,
    p_types_canonical = c(launch_pars, "B", "A", "t0", "k"),
    transform = transform_spec,
    bound = list(
      minmax = minmax,
      exception = exception),
    Ttransform = function(pars, dadm) {
      if (correlated) {
        # Shared with RDMSWTN(correlate = "drifts"); see
        # .apply_drift_factor_rho() in R/model_RDM.R for the encoding.
        pars <- .apply_drift_factor_rho(pars, dadm, sv = pars[, "sv"],
                                        model = "BAwLcorr")
      }
      lambda_factor <- if (erlang_shape_cpp == 2L) 2 else 1
      n <- nrow(pars)
      mG_val <- pars[, "mG"]
      mK_val <- pars[, "mK"]
      lg <- if (has_guess) ifelse(mG_val <= 0, 0, lambda_factor / mG_val) else rep(0, n)
      lk <- if (has_kill)  ifelse(mK_val <= 0, 0, lambda_factor / mK_val) else rep(0, n)
      timed <- cbind(lambda_g = lg, lambda_k = lk)
      # The kernel column order is the p_types prefix, so the launch pair must
      # lead whichever names it uses.
      lead <- c(launch_pars, "B", "A", "t0", "k")
      extra_drop <- c(lead, "mG", "mK", "omega")
      extra <- pars[, setdiff(colnames(pars), extra_drop), drop = FALSE]
      if (erlang_mixed) {
        pars <- cbind(
          pars[, lead, drop = FALSE],
          timed,
          omega = pars[, "omega"],
          extra
        )
      } else {
        pars <- cbind(
          pars[, lead, drop = FALSE],
          timed,
          extra
        )
      }
      cbind(pars, b = pars[, "B"] + pars[, "A"])
    },
    rfun = function(data, pars) {
      if (correlated) {
        .rfun_BAwL_corr(data$lR, pars, ok = attr(pars, "ok"), posdrift = posdrift,
                        erlang = erlang_shape_cpp, guess = has_guess,
                        global = erlang_type == "global_kill")
      } else {
        .rfun_BAwL(data$lR, pars, ok = attr(pars, "ok"), posdrift = posdrift,
                   erlang = erlang_shape_cpp, guess = has_guess,
                   global = erlang_type == "global_kill", launch = launch)
      }
    },
    dfun = function(rt, pars) dBAwL(rt, pars, posdrift = posdrift,  erlang = erlang_shape_cpp,
                                      guess = has_guess, launch = launch),
    pfun = function(rt, pars) pBAwL(rt, pars, posdrift = posdrift,  erlang = erlang_shape_cpp,
                                      guess = has_guess, launch = launch),
    log_likelihood = if (correlated) {
      function(pars, dadm, model, min_ll = log(1e-10)) {
        stop("BAwLcorr likelihood is implemented in the C++ race path; use fast_path=TRUE.")
      }
    } else {
      function(pars, dadm, model, min_ll = log(1e-10)) {
        log_likelihood_race_missing(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
      }
    }
  )
}

#' Correlated Ballistic Accumulator with Leak
#'
#' This is the `correlated = TRUE` form of [BAwL()]. It adds a one-factor
#' correlation to the Gaussian drift draws while retaining the BAwL leak,
#' start-point, non-decision-time, clock, and positive-drift parameterization.
#' In an expanded race design with `lM`, `rho` is the pairwise correlation
#' between the correct and incorrect racers' *underlying untruncated* Gaussian
#' drifts for each trial type. The correct racer is the positive reference and
#' the incorrect racer receives the cell sign.
#'
#' `rho` is a direct natural-scale correlation in `[-1, 1]`, represented by a
#' scaled probit transform. It must be shared by all rows of a trial. A
#' `matchfun` is therefore required so that the generated `lM` factor can
#' identify the correct racer. Structurally independent racers, such as a
#' PM/false-alarm row, should be assigned `rho = 0` using a participation
#' factor. For example:
#'
#' ```r
#' coupled <- function(d)
#'   factor(d$lR != "pm", c(FALSE, TRUE), c("no", "yes"))
#' # In the design: rho ~ 0 + coupled,
#' # constants = c(rho_coupledno = 0)
#' ```
#'
#' With `posdrift = TRUE`, BAwLcorr conditions the joint correlated Gaussian
#' drift vector on all active drifts being positive. Thus `rho` is the
#' correlation of the underlying untruncated Gaussian, while the marginal
#' drift draws retain standard positive-drift BAwL semantics. The marginal
#' truncated correlation is therefore not equal to `rho` in general.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *v* | identity | \[-Inf, Inf\] | 1 | | Mean normal launch strength. |
#' | *sv* | log | \[0, Inf\] | log(1) | | SD of normal launch strength. |
#' | *shape* | log | \[0, Inf\] | log(1) | | Weibull shape (Weibull launch only). |
#' | *scale* | log | \[0, Inf\] | log(1) | | Weibull scale (Weibull launch only). |
#' | *B* | log | \[0, Inf\] | log(1) | *b* = *B* + *A* | Distance from the upper start-point range to the threshold. |
#' | *A* | log | \[0, Inf\] | log(0) | | Start-point range. |
#' | *t0* | log | \[0, Inf\] | log(0) | | Non-decision time. |
#' | *k* | log | \[0, Inf\] | log(0) | | Leak rate. |
#'
#' All BAwL parameters are described in [BAwL()], including `v`, `sv`, `B`,
#' `A`, `t0`, and `k`. Optional fitting parameters are `pContaminant`, the
#' omission probability, and `pGuess`, the uniform-outlier probability; `rho`
#' is the optional correlation parameter.
#'
#' @param posdrift Logical. If `TRUE` (default), drift rates are jointly
#'   conditioned to be positive; if `FALSE`, use untruncated normal drifts.
#' @param erlang_shape Integer `1` for exponential clocks, `2` for Erlang-2,
#'   or `"mixed"` for an Erlang-1/Erlang-2 mixture.
#' @param erlang_type Clock configuration, as in [BAwL()].
#' @return A model list defining the correlated BAwL race model.
#' @export
BAwLcorr <- function(posdrift = TRUE, erlang_shape = 1L,
                     erlang_type = c("none", "local_kill", "global_kill", "local_guess", "local_kill_guess")) {
  BAwL(posdrift = posdrift, erlang_shape = erlang_shape,
       erlang_type = erlang_type, correlated = TRUE)
}
