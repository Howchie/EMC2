.lba_dfun <- function (rt, pars, posdrift = TRUE)
  # posdrift = truncated positive normal rates
  # robust slower, deals with extreme rate values
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
  # robust slower, deals with extreme rate values
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
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
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
#' | **Parameter** | **Transform** | **Natural scale** | **Default**   | **Mapping**                    | **Interpretation**                                            |
#' |-----------|-----------|---------------|-----------|----------------------------|-----------------------------------------------------------|
#' | *v*       | -         | \[-Inf, Inf\] | 1         |                            | Mean evidence-accumulation rate                                              |
#' | *A*       | log       | \[0, Inf\]    | log(0)    |                            | Between-trial variation (range) in start point                     |
#' | *B*       | log       | \[0, Inf\]    | log(1)    | *b* = *B*+*A*              | Distance from *A* to *b* (response threshold)                                       |
#' | *t0*      | log       | \[0, Inf\]    | log(0)    |                            | Non-decision time                                         |
#' | *sv*      | log       | \[0, Inf\]    | log(1)    |                            | Between-trial variation in evidence-accumulation rate                      |
#' | *pContaminant* | probit | \[0, 1\] | qnorm(0) | | Optional contamination probability handled by the data pipeline |
#'
#'
#' All core LBA parameters are estimated on the log scale, except for the drift
#' rate which is estimated on the real line. `pContaminant` is estimated on the
#' probit scale and is generic nuisance infrastructure rather than an LBA
#' accumulator parameter.
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
    p_types=c("v" = 1,"sv" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0), "pContaminant"=qnorm(0)),
    p_types_canonical = c("v", "sv", "B", "A", "t0"),
    transform=list(func=c(v = "identity",sv = "exp", B = "exp", A = "exp",t0 = "exp",pContaminant="pnorm")),
    bound=list(minmax=cbind(v=c(-Inf,Inf),sv = c(1e-4, Inf), A=c(1e-4,Inf),B=c(1e-4,Inf),t0=c(0.05,Inf),pContaminant=c(0.001,0.999)),
               exception=c(A=0,pContaminant=0)),
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
#' have been evaluated. The capacity extension adds the following parameters:
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Interpretation** |
#' |---|---|---|---|---|
#' | *kappa* | identity | \[-Inf, Inf\] | 0 | Additive mean capacity shift shared by both active targets in an `AB` trial. |
#' | *tau* | log | \[0, Inf\] | log(0) | Between-trial SD of the shared additive drift shift in an `AB` trial. |
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

dBAwL <- function(rt, pars, posdrift = TRUE, erlang = 1L, guess = FALSE) {
  if (!all(c("lambda_g", "lambda_k") %in% colnames(pars))) {
    stop("BAwL requires parameter columns 'lambda_g' and 'lambda_k'.")
  }
  dt  <- rt - pars[, "t0"]
  erl <- (pars[, "lambda_g"] > 0) | (pars[, "lambda_k"] > 0)
  ok  <- (rt > 0) & ((dt > 0) | erl) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok) | !is.finite(dt)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- dkilledleakyba(
      t = rt[ok], v = pars[ok, "v"], b = pars[ok, "b"], A = pars[ok, "A"],
      sv = pars[ok, "sv"], t0 = pars[ok, "t0"], k = pars[ok, "k"],
      lambda_g = pars[ok, "lambda_g"], lambda_k = pars[ok, "lambda_k"],
      posdrift = posdrift, log_out = FALSE,
      kill_shape = as.integer(erlang), guess = guess,
      erlang_omega = .rdmswtn_erlang_omega(pars[ok, , drop = FALSE], erlang)
    )
  }
  out
}

pBAwL <- function(rt, pars, posdrift = TRUE, erlang = 1L, guess = FALSE) {
  if (!all(c("lambda_g", "lambda_k") %in% colnames(pars))) {
    stop("BAwL requires parameter columns 'lambda_g' and 'lambda_k'.")
  }
  dt  <- rt - pars[, "t0"]
  erl <- (pars[, "lambda_g"] > 0) | (pars[, "lambda_k"] > 0)
  ok  <- (rt > 0) & ((dt > 0) | erl) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- pkilledleakyba(
      t = rt[ok], v = pars[ok, "v"], b = pars[ok, "b"], A = pars[ok, "A"],
      sv = pars[ok, "sv"], t0 = pars[ok, "t0"], k = pars[ok, "k"],
      lambda_g = pars[ok, "lambda_g"], lambda_k = pars[ok, "lambda_k"],
      posdrift = posdrift, log_out = FALSE,
      kill_shape = as.integer(erlang), guess = guess,
      erlang_omega = .rdmswtn_erlang_omega(pars[ok, , drop = FALSE], erlang)
    )
  }
  out
}

rBAwL <- function(lR, pars, ok = rep(TRUE, length(lR)),
                  p_types = c("v", "sv", "b", "A", "t0", "k", "lambda_g", "lambda_k"),
                  posdrift = TRUE, eps = 1e-10, erlang = 1L, guess = FALSE, global = FALSE,
                  .drifts = NULL) {
  if (!all(c("lambda_g", "lambda_k") %in% colnames(pars))) {
    stop("BAwL requires parameter columns 'lambda_g' and 'lambda_k'.")
  }
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
    drifts <- msm::rtnorm(nrow(pars), mean = pars[, "v"], sd = pars[, "sv"], lower = lower)
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
  bad_col <- apply(dt, 2, function(x) all(is.infinite(x)))
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
#' `b = B + A`, and the drift is a normal draw with mean `v_i` and SD `sv_i`.
#' With leak rate `k`, the evidence trajectory is
#' `x_i(t) = x_i(0) * exp(-k * t) + v_i / k * (1 - exp(-k * t))` for `k > 0`;
#' its `k = 0` limit is the ordinary ballistic trajectory. If the leaky
#' asymptote cannot reach the threshold, that accumulator has no finite hit.
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
#' | *mG* | log | \[0, Inf\] | log(1) | *lambda_g* = *q* / *mG* | Mean of the optional guess clock. |
#' | *mK* | log | \[0, Inf\] | log(1) | *lambda_k* = *q* / *mK* | Mean of the optional kill clock. |
#' | *omega* | probit | \[0, 1\] | qnorm(.5) | | Probability of the Erlang-1 component in mixed mode. |
#' | *pContaminant* | probit | \[0, 1\] | qnorm(0) | | Optional contamination probability handled by the data pipeline. |
#' | *rho* | scaled probit | \[-1, 1\] | qnorm(.5) | | Direct cell-level correlation of the underlying Gaussian drifts; only when `correlated = TRUE`. |
#'
#' Here `q = 1` for Erlang-1 clocks and `q = 2` for Erlang-2 clocks. In
#' `erlang_shape = "mixed"`, the clock is Erlang-1 with probability `omega` and
#' Erlang-2 otherwise, with both components having the same mean `mG` or `mK`.
#' The internal rates `lambda_g` and `lambda_k` are generated by the
#' `Ttransform`; they are not user-facing model parameters and should not be
#' put in a `design()` formula.
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
#' @return A model list defining the BAwL race model.
#'
#' @export
BAwL <- function(posdrift = TRUE, erlang_shape = 1L,
                 erlang_type = c("none", "local_kill", "global_kill", "local_guess", "local_kill_guess"),
                 correlated = FALSE) {
  erlang_type <- match.arg(erlang_type)
  erlang_mixed <- identical(erlang_shape, "mixed")
  erlang_shape_cpp <- if (erlang_mixed) 3L else as.integer(erlang_shape)
  
  has_guess <- erlang_type %in% c("local_guess", "local_kill_guess")
  has_kill  <- erlang_type %in% c("local_kill", "global_kill", "local_kill_guess")
  
  base_name <- paste0(ifelse(posdrift, "BAwL", "BAwLIO"),
                    if (erlang_mixed) "_EMIX" else if (erlang_shape_cpp >= 2L) "_E2" else "")
  type_suffix <- if (erlang_type == "local_guess") "_LOCAL_GUESS"
                    else if (erlang_type == "local_kill_guess") "_LOCAL_KILL_GUESS"
                    else if (erlang_type == "local_kill") "_LOCAL_KILL"
                    else if (erlang_type == "global_kill") "_GLOBAL_KILL"
                    else ""

  p_types <- c("v"  = 1, "sv" = log(1), "B" = log(1), "A" = log(0),
                "t0" = log(0), "k" = log(0))
  transform <- c(v = "identity", sv = "exp", B = "exp",
                              A = "exp", t0 = "exp", k = "exp")
  minmax <- cbind(v  = c(-Inf, Inf), sv = c(1e-4, Inf),
                     A  = c(1e-4, Inf), B  = c(1e-4, Inf),
                     t0 = c(0.05, Inf), k  = c(1e-4, Inf))
  exception <- c(A = 0, k = 0)
  
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

  p_types <- c(p_types, pContaminant = qnorm(0))
  transform <- c(transform, pContaminant = "pnorm")
  minmax <- cbind(minmax, pContaminant = c(0.001, 0.999))
  exception <- c(exception, pContaminant = 0)

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
    p_types = p_types,
    p_types_canonical = c("v", "sv", "B", "A", "t0", "k"),
    transform = transform_spec,
    bound = list(
      minmax = minmax,
      exception = exception),
    Ttransform = function(pars, dadm) {
      if (correlated) {
        if (is.null(dadm) || !"lM" %in% names(dadm) ||
            nrow(pars) != nrow(dadm)) {
          stop("BAwLcorr requires a matchfun-generated lM role indicator in the expanded data.")
        }
        # The sampled coefficient is one direct correlation for the cell.
        # Correct is the positive reference racer; incorrect receives the cell
        # sign.  Independent rows (for example PM) must be made structurally
        # zero by the rho design and its constants.
        rho_cell <- pars[, "rho"]
        correct <- as.character(dadm$lM) == "TRUE"
        if (anyNA(correct))
          stop("BAwLcorr matchfun produced missing lM values.")
        pars[, "rho"] <- ifelse(
          correct, abs(rho_cell), sign(rho_cell) * abs(rho_cell)
        )
      }
      lambda_factor <- if (erlang_shape_cpp == 2L) 2 else 1
      n <- nrow(pars)
      mG_val <- pars[, "mG"]
      mK_val <- pars[, "mK"]
      lg <- if (has_guess) ifelse(mG_val <= 0, 0, lambda_factor / mG_val) else rep(0, n)
      lk <- if (has_kill)  ifelse(mK_val <= 0, 0, lambda_factor / mK_val) else rep(0, n)
      timed <- cbind(lambda_g = lg, lambda_k = lk)
      extra_drop <- c("v", "sv", "B", "A", "t0", "k", "mG", "mK", "omega")
      extra <- pars[, setdiff(colnames(pars), extra_drop), drop = FALSE]
      if (erlang_mixed) {
        pars <- cbind(
          pars[, c("v", "sv", "B", "A", "t0", "k"), drop = FALSE],
          timed,
          omega = pars[, "omega"],
          extra
        )
      } else {
        pars <- cbind(
          pars[, c("v", "sv", "B", "A", "t0", "k"), drop = FALSE],
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
                   global = erlang_type == "global_kill")
      }
    },
    dfun = function(rt, pars) dBAwL(rt, pars, posdrift = posdrift,  erlang = erlang_shape_cpp,
                                      guess = has_guess),
    pfun = function(rt, pars) pBAwL(rt, pars, posdrift = posdrift,  erlang = erlang_shape_cpp,
                                      guess = has_guess),
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
#' All BAwL parameters are described in [BAwL()], including `v`, `sv`, `B`,
#' `A`, `t0`, `k`, `mG`, `mK`, `omega`, and `pContaminant`. This constructor
#' always sets `correlated = TRUE` internally.
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
