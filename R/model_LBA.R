#### Standard LBA ----
# # Moved to C+_ in model_LBA.cpp
#
# pnormP <- function (x, mean = 0, sd = 1, lower.tail = TRUE)
#       ifelse(abs(x) < 7, pnorm(x, mean = mean, sd = sd, lower.tail = lower.tail),
#              ifelse(x < 0, 0, 1))
#
# dnormP <- function (x, mean = 0, sd = 1)
#       ifelse(abs(x) < 7, dnorm(x, mean = mean, sd = sd), 0)
#
#
# dlba_norm <- function (dt,A,b,v,sv,posdrift=TRUE,robust=FALSE)
#     # like dlba_norm_core but t0 dealt with outside (removed from dt)
# {
#
#
#     if (robust) {
#       pnorm1 <- pnormP
#       dnorm1 <- dnormP
#     } else {
#       pnorm1 <- pnorm
#       dnorm1 <- dnorm
#     }
#
#     if (posdrift)
#       denom <- pmax(pnorm1(v/sv), 1e-10) else
#         denom <- rep(1, length(t))
#
#     A_small <- A < 1e-10
#     if (any(A_small)) {
#       out <- numeric(length(dt))
#       out[A_small] <- pmax(0, ((b[A_small]/dt[A_small]^2) *
#                                  dnorm1(b[A_small]/dt[A_small],v[A_small], sd = sv[A_small]))/denom[A_small])
#       zs <- dt[!A_small] * sv[!A_small]
#       zu <- dt[!A_small] * v[!A_small]
#       chiminuszu <- b[!A_small] - zu
#       chizu <- chiminuszu/zs
#       chizumax <- (chiminuszu - A[!A_small])/zs
#       out[!A_small] <- pmax(0, (v[!A_small] * (pnorm1(chizu) -
#                                                  pnorm1(chizumax)) + sv[!A_small] * (dnorm1(chizumax) -
#                                                                                        dnorm1(chizu)))/(A[!A_small] * denom[!A_small]))
#       return(out)
#     } else {
#       zs <- dt * sv
#       zu <- dt * v
#       chiminuszu <- b - zu
#       chizu <- chiminuszu/zs
#       chizumax <- (chiminuszu - A)/zs
#       return(pmax(0, (v * (pnorm1(chizu) - pnorm1(chizumax)) +
#                         sv * (dnorm1(chizumax) - dnorm1(chizu)))/(A * denom)))
#     }
# }
#
#
# plba_norm <- function (dt,A,b,v,sv,posdrift=TRUE,robust=FALSE)
#     # like plba_norm_core but t0 dealt with outside (removed from dt)
# {
#
#     if (robust) {
#       pnorm1 <- pnormP
#       dnorm1 <- dnormP
#     } else {
#       pnorm1 <- pnorm
#       dnorm1 <- dnorm
#     }
#     if (posdrift)
#       denom <- pmax(pnorm1(v/sv), 1e-10) else
#         denom <- 1
#     A_small <- A < 1e-10
#     if (any(A_small)) {
#       out <- numeric(length(dt))
#       out[A_small] <- pmin(1, pmax(0, (pnorm1(b[A_small]/dt[A_small],
#                                               mean = v[A_small], sd = sv[A_small],
#                                               lower.tail = FALSE))/denom[A_small]))
#       zs <- dt[!A_small] * sv[!A_small]
#       zu <- dt[!A_small] * v[!A_small]
#       chiminuszu <- b[!A_small] - zu
#       xx <- chiminuszu - A[!A_small]
#       chizu <- chiminuszu/zs
#       chizumax <- xx/zs
#       tmp1 <- zs * (dnorm1(chizumax) - dnorm1(chizu))
#       tmp2 <- xx * pnorm1(chizumax) - chiminuszu * pnorm1(chizu)
#       out[!A_small] <- pmin(pmax(0, (1 + (tmp1 + tmp2)/A[!A_small])/denom[!A_small]),1)
#       return(out)
#     } else {
#       zs <- dt * sv
#       zu <- dt * v
#       chiminuszu <- b - zu
#       xx <- chiminuszu - A
#       chizu <- chiminuszu/zs
#       chizumax <- xx/zs
#       tmp1 <- zs * (dnorm1(chizumax) - dnorm1(chizu))
#       tmp2 <- xx * pnorm1(chizumax) - chiminuszu * pnorm1(chizu)
#       return(pmin(pmax(0, (1 + (tmp1 + tmp2)/A)/denom), 1))
#     }
# }
#
#

dLBA <- function (rt, pars, posdrift = TRUE)
  # posdrift = truncated positive normal rates
  # robust slower, deals with extreme rate values
{
  dt <- rt - pars[,"t0"]
  ok <- (dt>0) & (pars[,"b"] >= pars[,"A"])
  ok[is.na(ok) | !is.finite(dt)] <- FALSE
  out <- numeric(length(dt))
  out[ok] <- dlba(t = dt[ok], A = pars[ok,"A"], b = pars[ok,"b"],
                         v = pars[ok,"v"], sv = pars[ok,"sv"],
                         posdrift = posdrift)
  out
}

pLBA <- function (rt, pars, posdrift = TRUE)
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
  out
}


rLBA <- function(lR,pars,p_types=c("v","sv","b","A","t0"),
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
  out
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
#'
#'
#' All parameters are estimated on the log scale, except for the drift rate which is estimated on the real line.
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
    # p_vector transform, sets sv as a scaling parameter
    p_types=c("v" = 1,"sv" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0), "pContaminant"=qnorm(0)),
    transform=list(func=c(v = "identity",sv = "exp", B = "exp", A = "exp",t0 = "exp",pContaminant="pnorm")),
    bound=list(minmax=cbind(v=c(-Inf,Inf),sv = c(0, Inf), A=c(1e-4,Inf),B=c(1e-4,Inf),t0=c(0.05,Inf),pContaminant=c(0.001,0.999)),
               exception=c(A=0,pContaminant=0)),
    # Transform to natural scale
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
      pars
    },
    # Random function for racing accumulator
    rfun=ifelse(posdrift,function(data,pars) rLBA(data$lR,pars,ok = attr(pars, "ok"),posdrift=TRUE),
                function(data,pars) rLBA(data$lR,pars,ok = attr(pars, "ok"),posdrift=FALSE)),
    # Density function (PDF) for single accumulator
    dfun=ifelse(posdrift,function(rt,pars) dLBA(rt,pars,posdrift=TRUE),
                function(rt,pars) dLBA(rt,pars,posdrift=FALSE)),
    # Probability function (CDF) for single accumulator
    pfun=ifelse(posdrift,function(rt,pars) pLBA(rt,pars,posdrift=TRUE),
                function(rt,pars) pLBA(rt,pars,posdrift=FALSE)),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
      log_likelihood_race_missing(pars=pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}

#' LBA Logical Rules Model
#' @export
#'
LogicalRulesLBA <- function(posdrift = TRUE, fast_path=TRUE){
  list(
    type="RACE",
    # Note: `calc_ll()` infers LBA `posdrift` from whether `c_name` contains "IO".
    # Keep this in sync so likelihood and simulation use the same setting.
    c_name = paste0("LBA_LogicalRules",ifelse(posdrift,"","IO")),
    # p_vector transform, sets sv as a scaling parameter
    p_types=c("v" = 1,"sv" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0)),
    transform=list(func=c(v = "identity",sv = "exp", B = "exp", A = "exp",t0 = "exp")),
    bound=list(minmax=cbind(v=c(-Inf,Inf),sv = c(0, Inf), A=c(1e-4,Inf),B=c(0,Inf),t0=c(0.05,Inf)),
               exception=c(A=0)),
    # Transform to natural scale
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
      pars
    },
    # Random function for racing accumulator
    rfun=function(data,pars) rLBA(data$lR,pars,ok = attr(pars, "ok"),posdrift=ifelse(posdrift,TRUE,FALSE)),
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dLBA(rt,pars,posdrift=ifelse(posdrift,TRUE,FALSE)),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pLBA(rt,pars,posdrift=ifelse(posdrift,TRUE,FALSE)),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
      stop("LogicalRulesLBA: R likelihood path not implemented. Use fast_path=TRUE (the default).")
    }
  )
}

#' LBA Redundant Target Race Model
#' @export
#' 
RedundantTargetLBA <- function(posdrift = TRUE){
  list(
    type = "RACE",
    c_name = paste0("LBA_RedundantTarget", ifelse(posdrift, "", "IO")),
    p_types = c("v" = 1, "sv" = log(1), "B" = log(1), "A" = log(0), "t0" = log(0), "pContaminant" = qnorm(0)),
    transform = list(func = c(v = "identity", sv = "exp", B = "exp", A = "exp", t0 = "exp", pContaminant = "pnorm")),
    bound = list(
      minmax = cbind(v = c(-Inf, Inf), sv = c(0, Inf), A = c(1e-4, Inf), B = c(1e-4, Inf), t0 = c(0.05, Inf), pContaminant = c(0.001, 0.999)),
      exception = c(A = 0, pContaminant = 0)
    ),
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars, b = pars[, "B"] + pars[, "A"])
      pars
    },
    rfun = ifelse(posdrift,
      function(data,pars) rLBA(data$lR,pars,ok = attr(pars, "ok"),posdrift = TRUE),
      function(data,pars) rLBA(data$lR,pars,ok = attr(pars, "ok"),posdrift = FALSE)),
    dfun = ifelse(posdrift,
      function(rt,pars) dLBA(rt,pars,posdrift = TRUE),
      function(rt,pars) dLBA(rt,pars,posdrift = FALSE)),
    pfun = ifelse(posdrift,
      function(rt,pars) pLBA(rt,pars,posdrift = TRUE),
      function(rt,pars) pLBA(rt,pars,posdrift = FALSE)),
    log_likelihood = function(pars,dadm,model,min_ll = log(1e-10)){
      stop("RedundantTargetLBA: R likelihood path not implemented. Use the C++ path (default).")
    }
  )
}

#### BAwL (Ballistic Accumulator with Leak) ----

# Erlang-n kill survival: exp(-lambda*t) * sum_{m=0}^{n-1} (lambda*t)^m / m!
# n=1: exponential; n=2: Erlang-2.
erlang_surv <- function(t, lambda, n = 1L) {
  s <- exp(-lambda * t)
  if (n >= 2L) s <- s * (1 + lambda * t)
  s
}

dBAwL <- function(rt, pars, posdrift = TRUE, erlang = 1L) {
  dt  <- rt - pars[, "t0"]
  ok  <- (dt > 0) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok) | !is.finite(dt)] <- FALSE
  out <- numeric(length(dt))
  base <- dleakyba(t  = dt[ok],
                   A  = pars[ok, "A"],
                   b  = pars[ok, "b"],
                   v  = pars[ok, "v"],
                   sv = pars[ok, "sv"],
                   k  = pars[ok, "k"],
                   posdrift = posdrift)
  out[ok] <- base * erlang_surv(dt[ok], pars[ok, "lambda"], erlang)
  out
}

pBAwL <- function(rt, pars, posdrift = TRUE, erlang = 1L) {
  dt  <- rt - pars[, "t0"]
  ok  <- (dt > 0) & (pars[, "b"] >= pars[, "A"])
  ok[is.na(ok) | !is.finite(dt)] <- FALSE
  out <- numeric(length(dt))
  if (any(ok)) {
    out[ok] <- vapply(which(ok), function(i) {
      stats::integrate(
        f = function(x) {
          dleakyba(
            t = x, A = pars[i, "A"], b = pars[i, "b"], v = pars[i, "v"],
            sv = pars[i, "sv"], k = pars[i, "k"], posdrift = posdrift
          ) * erlang_surv(x, pars[i, "lambda"], erlang)
        },
        lower = 0, upper = dt[i], subdivisions = 200L, rel.tol = 1e-6
      )$value
    }, numeric(1))
  }
  out
}

rBAwL <- function(lR, pars, ok = rep(TRUE, length(lR)),
                  p_types = c("v", "sv", "b", "A", "t0", "k", "lambda"),
                  posdrift = TRUE, eps = 1e-10, erlang = 1L) {
  bad  <- rep(NA, length(lR) / length(levels(lR)))
  out  <- data.frame(R = bad, rt = bad)
  nr   <- length(levels(lR))
  dt   <- matrix(Inf, nrow = nr, ncol = nrow(pars) / nr)
  t0   <- pars[, "t0"]
  pars <- pars[ok, ]
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ", paste(p_types, collapse = " "))
  lower  <- if (posdrift) 0 else -Inf
  drifts <- msm::rtnorm(nrow(pars), mean = pars[, "v"], sd = pars[, "sv"], lower = lower)

  small_k <- pars[, "k"] < eps
  if (any(small_k))
    dt[small_k] <- (pars[small_k, "b"] - pars[small_k, "A"] * runif(sum(small_k))) /
                   drifts[small_k]

  big_k <- !small_k & (drifts > pars[, "k"] * pars[, "b"])
  if (any(big_k)) {
    num        <- drifts[big_k] - pars[big_k, "k"] * pars[big_k, "b"]
    den        <- drifts[big_k] - pars[big_k, "k"] * pars[big_k, "A"] * runif(sum(big_k))
    ratio      <- num / den
    ratio      <- pmin(pmax(ratio, .Machine$double.xmin), 1 - 1e-15)
    dt[big_k]  <- (-1 / pars[big_k, "k"]) * log(ratio)
  }
  dt[dt < 0] <- Inf
  lam_pos <- pars[, "lambda"] > 0
  tk <- ifelse(lam_pos, rexp(nrow(pars), pars[, "lambda"]), Inf)
  if (erlang >= 2L && any(lam_pos))
    tk[lam_pos] <- tk[lam_pos] + rexp(sum(lam_pos), pars[lam_pos, "lambda"])
  dt <- ifelse(dt < tk, dt, Inf)
  bad_col <- apply(dt, 2, function(x) all(is.infinite(x)))
  R   <- apply(dt, 2, which.min)
  pick <- cbind(R, 1:dim(dt)[2])
  rt   <- matrix(t0, nrow = nr)[pick] + dt[pick]
  R    <- factor(levels(lR)[R], levels = levels(lR))
  R[bad_col]  <- NA
  rt[bad_col] <- Inf
  ok  <- matrix(ok, nrow = length(levels(lR)))[1, ]
  out$R[ok]   <- levels(lR)[R][ok]
  out$R       <- factor(out$R, levels = levels(lR))
  out$rt[ok]  <- rt[ok]
  out
}

#' The Ballistic Accumulator with Leak (BAwL)
#'
#' Race model where each accumulator follows a leaky-integration trajectory and
#' races against an independent exponential killing clock.
#' Setting k=0 recovers the standard LBA decision dynamics; setting lambda=0
#' removes killing.
#'
#' Parameters as in LBA plus:
#' * k: leak rate (log scale, [0, Inf); default log(0) = effectively 0).
#' * lambda: killing/expiry rate (log scale, [0, Inf); default log(0) = 0).
#'
#' @export
BAwL <- function(posdrift = TRUE, erlang = 1L) {
  list(
    type   = "RACE",
    c_name = paste0(ifelse(posdrift, "BAwL", "BAwLIO"), if (erlang >= 2L) "_E2" else ""),
    p_types = c("v"  = 1, "sv" = log(1), "B" = log(1), "A" = log(0),
                "t0" = log(0), "k" = log(0), "lambda" = log(0), "pContaminant" = qnorm(0)),
    transform = list(func = c(v = "identity", sv = "exp", B = "exp",
                              A = "exp", t0 = "exp", k = "exp", lambda = "exp",
                              pContaminant = "pnorm")),
    bound = list(
      minmax = cbind(v  = c(-Inf, Inf), sv = c(0, Inf),
                     A  = c(1e-4, Inf), B  = c(1e-4, Inf),
                     t0 = c(0.05, Inf), k  = c(1e-4, Inf), lambda = c(1e-4, Inf),
                     pContaminant = c(0.001, 0.999)),
      exception = c(A = 0, pContaminant = 0, k = 0, lambda = 0)),
    Ttransform = function(pars, dadm) {
      cbind(pars, b = pars[, "B"] + pars[, "A"])
    },
    rfun = if (posdrift)
             function(data, pars) rBAwL(data$lR, pars, ok = attr(pars, "ok"), posdrift = TRUE,  erlang = erlang)
           else
             function(data, pars) rBAwL(data$lR, pars, ok = attr(pars, "ok"), posdrift = FALSE, erlang = erlang),
    dfun = if (posdrift)
             function(rt, pars) dBAwL(rt, pars, posdrift = TRUE,  erlang = erlang)
           else
             function(rt, pars) dBAwL(rt, pars, posdrift = FALSE, erlang = erlang),
    pfun = if (posdrift)
             function(rt, pars) pBAwL(rt, pars, posdrift = TRUE,  erlang = erlang)
           else
             function(rt, pars) pBAwL(rt, pars, posdrift = FALSE, erlang = erlang),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      log_likelihood_race_missing(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}
