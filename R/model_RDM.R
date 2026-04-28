
dRDM <- function(rt,pars)
  # density for single accumulator
{
  out <- numeric(length(rt))
  ok <- rt > pars[,"t0"] & !pars[,"v"] < 0  # code handles rate zero case
  ok[is.na(ok)] <- FALSE
  if (any(dimnames(pars)[[2]]=="s")) # rescale
    pars[ok,c("A","B","v")] <- pars[ok,c("A","B","v")]/pars[ok,"s"]
  out[ok] <- dWald(rt[ok],v=pars[ok,"v"],B=pars[ok,"B"],A=pars[ok,"A"],t0=pars[ok,"t0"])
  out
}


pRDM <- function(rt,pars)
  # cumulative density for single accumulator
{
  out <- numeric(length(rt))
  ok <- rt > pars[,"t0"] & !pars[,"v"] < 0  # code handles rate zero case
  ok[is.na(ok)] <- FALSE
  if (any(dimnames(pars)[[2]]=="s")) # rescale
    pars[ok,c("A","B","v")] <- pars[ok,c("A","B","v")]/pars[ok,"s"]
  out[ok] <- pWald(rt[ok],v=pars[ok,"v"],B=pars[ok,"B"],A=pars[ok,"A"],t0=pars[ok,"t0"])
  out
}

#### random

rWald <- function(n,B,v,A)
  # random function for single accumulator
{

  rwaldt <- function(n,k,l,tiny=1e-6) {
    # random sample of n from a Wald (or Inverse Gaussian)
    # k = criterion, l = rate, assumes sigma=1 Browninan motion
    # about same speed as statmod rinvgauss

    rlevy <- function(n=1, m=0, c=1) {
      if (any(c<0)) stop("c must be positive")
      c/qnorm(1-runif(n)/2)^2+m
    }

    flag <- l>tiny
    x <- rep(NA,times=n)

    x[!flag] <- rlevy(sum(!flag),0,k[!flag]^2)
    mu <- k/l
    lambda <- k^2

    y <- rnorm(sum(flag))^2
    mu.0 <- mu[flag]
    lambda.0 <- lambda[flag]

    x.0 <- mu.0 + mu.0^2*y/(2*lambda.0) -
      sqrt(4*mu.0*lambda.0*y + mu.0^2*y^2)*mu.0/(2*lambda.0)

    z <- runif(length(x.0))
    test <- mu.0/(mu.0+x.0)
    x.0[z>test] <- mu.0[z>test]^2/x.0[z>test]
    x[flag] <- x.0
    x[x<0] <- max(x)
    x
  }

  # Kluge to return Inf for negative rates
  out <- numeric(n)
  ok <- !v<0
  nok <- sum(ok)
  bs <- B[ok]+runif(nok,0,A[ok])
  out[ok] <- rwaldt(nok,k=bs,l=v[ok])
  out[!ok] <- Inf
  out
}

rRDM <- function(lR,pars,p_types=c("v","B","A","t0"),ok=rep(TRUE,dim(pars)[1]))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in
  # contiguous rows. "s" parameter will be used but can be ommitted
  #
  # test
  # pars=cbind(B=c(1,2),v=c(1,1),A=c(0,0),t0=c(.2,.2)); lR=factor(c(1,2))
{
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  if (any(dimnames(pars)[[2]]=="s")) # rescale
    pars[,c("A","B","v")] <- pars[,c("A","B","v")]/pars[,"s"]
  pars[,"B"][pars[,"B"]<0] <- 0 # Protection for negatives
  pars[,"A"][pars[,"A"]<0] <- 0
  bad <- rep(NA, length(lR)/length(levels(lR)))
  out <- data.frame(R = bad, rt = bad)
  nr <- length(levels(lR))
  dt <- matrix(Inf,nrow=nr,ncol=nrow(pars)/nr)
  t0 <- pars[,"t0"]
  pars <- pars[ok,]
  dt[ok] <- rWald(sum(ok),B=pars[,"B"],v=pars[,"v"],A=pars[,"A"])
  R <- max.col(-t(dt), ties.method='first')
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
  rt <- matrix(t0,nrow=nr)[pick] + dt[pick]
  out$R <- levels(lR)[R]
  out$R <- factor(out$R,levels=levels(lR))
  out$rt <- rt
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

RDM <- function(){
  list(
    type="RACE",
    c_name = "RDM",
    p_types=c("v" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),"s" = log(1), "pContaminant"=qnorm(0)),
    transform=list(func=c(v = "exp", B = "exp", A = "exp",t0 = "exp", s = "exp",pContaminant="pnorm")),
    bound=list(minmax=cbind(v=c(1e-3,Inf), B=c(0,Inf), A=c(1e-4,Inf),t0=c(0.05,Inf), s=c(0,Inf),pContaminant=c(0.001,0.999)),
               exception=c(A=0, v=0,pContaminant=0)),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
      pars
    },
    # Random function for racing accumulators
    rfun=function(data=NULL,pars)  rRDM(data$lR,pars,ok=attr(pars, "ok")),
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dRDM(rt,pars),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pRDM(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
      log_likelihood_race_missing(pars=pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}

# ============================================================================
# RDMGBM: Racing Geometric Brownian Motion with start-point variability
# ============================================================================

dRDMGBM <- function(rt, pars) {
  if (is.null(dim(pars)) || (dim(pars)[1] == 1 & length(rt) > 1)) {
    original_names <- names(pars); if (is.null(original_names)) original_names <- colnames(pars)
    pars <- matrix(pars, nrow=length(rt), ncol=length(pars),
                   dimnames=list(NULL, original_names), byrow=TRUE)
  }
  if (!("b" %in% colnames(pars)) && all(c("B", "A") %in% colnames(pars))) {
    pars <- cbind(pars, b=1 + pars[, "B"] + pars[, "A"])
  }
  out <- rep(NaN, length(rt))
  ok <- rt > pars[, "t0", drop=FALSE] & !pars[, "v", drop=FALSE] < 0
  ok[is.na(ok)] <- FALSE
  if (any(ok)) {
    out[ok] <- dGBMspv(rt[ok], v=pars[ok, "v", drop=FALSE], b=pars[ok, "b", drop=FALSE],
                       A=pars[ok, "A", drop=FALSE], t0=pars[ok, "t0", drop=FALSE],
                       s=pars[ok, "s", drop=FALSE])
  }
  out
}

pRDMGBM <- function(rt, pars) {
  if (is.null(dim(pars)) || (dim(pars)[1] == 1 & length(rt) > 1)) {
    original_names <- names(pars); if (is.null(original_names)) original_names <- colnames(pars)
    pars <- matrix(pars, nrow=length(rt), ncol=length(pars),
                   dimnames=list(NULL, original_names), byrow=TRUE)
  }
  if (!("b" %in% colnames(pars)) && all(c("B", "A") %in% colnames(pars))) {
    pars <- cbind(pars, b=1 + pars[, "B"] + pars[, "A"])
  }
  out <- rep(NaN, length(rt))
  ok <- rt > pars[, "t0", drop=FALSE] & !pars[, "v", drop=FALSE] < 0
  ok[is.na(ok)] <- FALSE
  if (any(ok)) {
    out[ok] <- pGBMspv(rt[ok], v=pars[ok, "v", drop=FALSE], b=pars[ok, "b", drop=FALSE],
                       A=pars[ok, "A", drop=FALSE], t0=pars[ok, "t0", drop=FALSE],
                       s=pars[ok, "s", drop=FALSE])
  }
  out
}

rGBM <- function(n, b, v, A, s = 1) {
  out <- rep(Inf, n)
  if (n <= 0) return(out)
  if (n > 1 && all(length(b) == 1, length(v) == 1, length(A) == 1, length(s) == 1)) {
    b <- rep(b, n); v <- rep(v, n); A <- rep(A, n); s <- rep(s, n)
  }
  A[A < 0] <- 0
  x0 <- 1 + runif(n, 0, A)
  d <- log(b / x0)
  mu_log <- v - 0.5 * s^2
  ok <- is.finite(mu_log) & is.finite(d) & (d > 0)
  if (any(ok)) {
    out[ok] <- statmod::rinvgauss(sum(ok), mean=d[ok] / mu_log[ok], shape=d[ok]^2 / s[ok]^2)
  }
  out
}

rRDMGBM <- function(lR, pars, p_types=c("v", "B", "A", "t0", "s"), ok=rep(TRUE, dim(pars)[1])) {
  if (!("b" %in% dimnames(pars)[[2]]) && all(c("B", "A") %in% dimnames(pars)[[2]])) {
    pars <- cbind(pars, b=1 + pars[, "B"] + pars[, "A"])
  }
  required <- c("v", "b", "A", "t0", "s")
  if (!all(required %in% dimnames(pars)[[2]]))
    stop("pars must have columns ", paste(required, collapse = " "))
  pars[, "b"][pars[, "b"] < 1 + 1e-8] <- 1 + 1e-8
  pars[, "A"][pars[, "A"] < 0] <- 0
  bad <- rep(NA, length(lR) / length(levels(lR)))
  out <- data.frame(R = bad, rt = bad)
  nr <- length(levels(lR))
  n_trials <- nrow(pars) / nr
  if (length(ok) != nrow(pars)) stop("ok must have length nrow(pars).")
  trial_ok <- colSums(matrix(ok, nrow = nr)) == nr
  dt <- matrix(NA_real_, nrow = nr, ncol = n_trials)
  t0 <- pars[, "t0"]
  pars_ok <- pars[ok, , drop = FALSE]
  if (nrow(pars_ok) > 0) {
    dt[ok] <- rGBM(sum(ok), b = pars_ok[, "b"], v = pars_ok[, "v"], A = pars_ok[, "A"], s = pars_ok[, "s"])
  }
  if (any(trial_ok)) {
    dt_valid <- dt[, trial_ok, drop = FALSE]
    R <- max.col(-t(dt_valid), ties.method = "first")
    pick <- cbind(R, which(trial_ok))
    rt <- matrix(t0, nrow = nr)[pick] + dt[pick]
    out$R[trial_ok] <- levels(lR)[R]
    out$rt[trial_ok] <- rt
  }
  out$R <- factor(out$R, levels = levels(lR))
  out
}

#' RDMGBM Model
#'
#' Racing geometric Brownian first-passage model with start-point variability.
#' Equivalent parameterization to RDMSWTN without `sv`.
#'
#' @export
#' 
RDMGBM <- function() {
  list(
    type = "RACE",
    c_name = "RDMGBM",
    p_types = c("v"=log(1), "B"=log(1), "A"=log(0), "t0"=log(0),
                "s"=log(1), "pContaminant"=qnorm(0)),
    transform = list(func = c(v="exp", B="exp", A="exp", t0="exp",
                               s="exp", pContaminant="pnorm")),
    bound = list(
      minmax = cbind(v=c(1e-3, Inf), B=c(0, Inf), A=c(0, Inf),
                     t0=c(0.05, Inf), s=c(0, Inf), pContaminant=c(0.001, 0.999)),
      exception = c(A=0, v=0, pContaminant=0),
      # Joint validity for GBM first-passage parameterization.
      # Enforces positive effective log-drift so simulators do not emit NA RTs.
      joint_ok = function(pars) {
        if (!all(c("v", "s") %in% colnames(pars))) return(rep(TRUE, nrow(pars)))
        is.finite(pars[, "v"]) & is.finite(pars[, "s"]) &
          pars[, "s"] > 0 & (pars[, "v"] > 0.5 * pars[, "s"]^2)
      }
    ),
    Ttransform = function(pars, dadm) {
      pars <- cbind(pars, b=1 + pars[, "B"] + pars[, "A"])
      pars
    },
    rfun = function(data=NULL, pars) rRDMGBM(data$lR, pars, ok=attr(pars, "ok")),
    dfun = function(rt, pars) dRDMGBM(rt, pars),
    pfun = function(rt, pars) pRDMGBM(rt, pars),
    log_likelihood = function(pars, dadm, model, min_ll=log(1e-10)) {
      log_likelihood_race_missing(pars=pars, dadm=dadm, model=model, min_ll=min_ll)
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
#' | *lambda*  | log       | \[0, Inf\]      | log(0)  | Exponential killing rate               |
#' 
#' @return a list of parameters
#' 
#' @export
#' 
RDMSWTN <- function(erlang = 1L, guess = FALSE){
  list(
    type="RACE",
    c_name = paste0(if (erlang >= 2L) "RDMSWTN_E2" else "RDMSWTN", if (guess) "_GUESS" else ""),
    p_types=c("v"=log(1), "B"=log(1), "A"=log(0), "t0"=log(0),
              "s"=log(1), "sv"=log(0), "lambda"=log(0), "pContaminant"=qnorm(0)),
    transform=list(func=c(v="exp", B="exp", A="exp", t0="exp",
                          s="exp", sv="exp", lambda="exp", pContaminant="pnorm")),
    bound=list(minmax=cbind(v=c(1e-3,Inf), B=c(0,Inf), A=c(0,Inf),
                            t0=c(0.05,Inf), s=c(0,Inf), sv=c(0,Inf), lambda=c(1e-4,Inf),
                            pContaminant=c(0.001,0.999)),
               exception=c(A=0, v=0, sv=0, lambda=0, pContaminant=0)),
    Ttransform = function(pars, dadm) {
      pars <- cbind(pars, b=pars[,"B"] + pars[,"A"])
      pars
    },
    rfun=function(data=NULL, pars) rRDMSWTN(data$lR, pars, ok=attr(pars, "ok"), erlang=erlang, guess=guess),
    dfun=function(rt, pars) dRDMSWTN(rt, pars, erlang=erlang),
    pfun=function(rt, pars) pRDMSWTN(rt, pars, erlang=erlang),
    log_likelihood=function(pars, dadm, model, min_ll=log(1e-10)){
      log_likelihood_race_missing(pars=pars, dadm=dadm, model=model, min_ll=min_ll)
    }
  )
}

dRDMSWTN <- function(rt, pars, erlang = 1L) {
  if (is.null(dim(pars)) || (dim(pars)[1]==1 & length(rt)>1)) {
    original_names <- names(pars); if (is.null(original_names)) original_names <- colnames(pars)
    pars <- matrix(pars, nrow=length(rt), ncol=length(pars),
                   dimnames=list(NULL, original_names), byrow=TRUE)
  }
  out <- rep(NaN, length(rt))
  ok <- rt > pars[,"t0",drop=FALSE] & !pars[,"v",drop=FALSE] < 0
  ok[is.na(ok)] <- FALSE
  if (any(ok)) {
    if (!("lambda" %in% colnames(pars))) stop("RDMSWTN requires parameter column 'lambda'.")
    if (any(dimnames(pars)[[2]] == "s")) {
      pars_ok <- pars[ok,,drop=FALSE]
      pars_ok[,c("A","b","v","sv")] <- pars_ok[,c("A","b","v","sv")] / pars_ok[,"s"]
      pars[ok,] <- pars_ok
    }
    out[ok] <- dSWTNspv(rt[ok], v=pars[ok,"v",drop=FALSE], b=pars[ok,"b",drop=FALSE],
                        A=pars[ok,"A",drop=FALSE], t0=pars[ok,"t0",drop=FALSE],
                        sv=pars[ok,"sv",drop=FALSE], lambda=pars[ok,"lambda",drop=FALSE],
                        kill_shape=erlang)
  }
  out
}

pRDMSWTN <- function(rt, pars, erlang = 1L) {
  if (is.null(dim(pars)) || (dim(pars)[1]==1 & length(rt)>1)) {
    original_names <- names(pars); if (is.null(original_names)) original_names <- colnames(pars)
    pars <- matrix(pars, nrow=length(rt), ncol=length(pars),
                   dimnames=list(NULL, original_names), byrow=TRUE)
  }
  out <- rep(NaN, length(rt))
  ok <- rt > pars[,"t0",drop=FALSE] & !pars[,"v",drop=FALSE] < 0
  ok[is.na(ok)] <- FALSE
  if (any(ok)) {
    if (!("lambda" %in% colnames(pars))) stop("RDMSWTN requires parameter column 'lambda'.")
    if (any(dimnames(pars)[[2]] == "s")) {
      pars_ok <- pars[ok,,drop=FALSE]
      pars_ok[,c("A","b","v","sv")] <- pars_ok[,c("A","b","v","sv")] / pars_ok[,"s"]
      pars[ok,] <- pars_ok
    }
    out[ok] <- pSWTNspv(rt[ok], v=pars[ok,"v",drop=FALSE], b=pars[ok,"b",drop=FALSE],
                        A=pars[ok,"A",drop=FALSE], t0=pars[ok,"t0",drop=FALSE],
                        sv=pars[ok,"sv",drop=FALSE], lambda=pars[ok,"lambda",drop=FALSE],
                        kill_shape=erlang)
  }
  out
}

rRDMSWTN <- function(lR, pars, p_types=c("v","b","A","t0","sv","lambda"),
                     ok=rep(TRUE, dim(pars)[1]), erlang=1L, guess=FALSE) {
  if (!is.null(attr(pars, "ok"))) ok <- attr(pars, "ok")
  if (is.null(dim(pars)) || (dim(pars)[1]==1 & length(lR)>1)) {
    original_names <- names(pars); if (is.null(original_names)) original_names <- colnames(pars)
    pars <- matrix(pars, nrow=length(lR), ncol=length(pars),
                   dimnames=list(NULL, original_names), byrow=TRUE)
  }
  if (!("lambda" %in% colnames(pars))) stop("RDMSWTN requires parameter column 'lambda'.")
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ", paste(p_types, collapse=" "))
  if (any(dimnames(pars)[[2]] == "s"))
    pars[,c("A","b","v","sv")] <- pars[,c("A","b","v","sv")] / pars[,"s"]
  pars[,"b"][pars[,"b"] < 0] <- 0
  pars[,"A"][pars[,"A"] < 0] <- 0
  bad <- rep(NA, length(lR)/length(levels(lR)))
  out <- data.frame(R=bad, rt=bad)
  nr <- length(levels(lR))
  dt <- matrix(Inf, nrow=nr, ncol=nrow(pars)/nr)
  t0 <- pars[,"t0"]
  # If guessing, drawn ONE kill process per trial
  if (guess) {
    n_trials <- nrow(pars)/nr
    # lambda is constant across accumulators for a trial in RDMSWTN
    lambda_trials <- matrix(pars[,"lambda"], nrow=nr)[1,]
    tk <- rexp(n_trials, rate=lambda_trials)
    if (erlang >= 2L) tk <- tk + rexp(n_trials, rate=lambda_trials)
    # Simulator rSWTN without killing (k=0)
    pars[,"lambda"] <- 0
  }
  pars <- pars[ok,,drop=FALSE]
  dt[ok] <- rSWTN(sum(ok), b=pars[,"b"], v=pars[,"v"], A=pars[,"A"], sv=pars[,"sv"],
                  k=pars[,"lambda"], erlang=erlang)
  
  if (guess) {
    # Response is a guess if tk < any(dt)
    is_guess <- tk < apply(dt, 2, min)
    if (any(is_guess)) {
      dt[, is_guess] <- Inf
      ok_mat <- matrix(ok, nrow=nr)
      # For each guess trial, pick one of the active, non-nogo accumulators
      for (trial_idx in which(is_guess)) {
        active_indices <- which(ok_mat[, trial_idx] & levels(lR) != "nogo")
        if (length(active_indices) > 0) {
          # Use sample safely
          winner_acc <- if (length(active_indices) > 1) sample(active_indices, 1) else active_indices
          dt[winner_acc, trial_idx] <- tk[trial_idx]
        }
      }
    }
  }

  bad_col <- apply(dt, 2, function(x) all(is.infinite(x)))
  R <- apply(dt, 2, which.min)
  pick <- cbind(R, 1:dim(dt)[2])
  rt <- matrix(t0, nrow=nr)[pick] + dt[pick]
  out$R <- levels(lR)[R]
  out$R <- factor(out$R, levels=levels(lR))
  out$rt <- rt
  out$R[bad_col] <- NA
  out$rt[bad_col] <- Inf
  out
}

rSWTN <- function(n, b, v, A, sv, s=1, k=0, erlang=1L) {
  out <- numeric(n)
  if (n > 1 & all(length(A)==1, length(v)==1, length(b)==1, length(sv)==1, length(k)==1)) {
    A  <- rep(A, n)
    b  <- rep(b, n)
    v  <- rep(v, n)
    sv <- rep(sv, n)
    k  <- rep(k, n)
  }
  b <- ifelse(A==0, b, runif(n, b-A, b))
  l <- ifelse(sv==0, v, msm::rtnorm(n, mean=v, sd=sv, lower=0, upper=Inf))
  ok <- !l < 0
  nok <- sum(ok)
  out[ok] <- statmod::rinvgauss(nok, mean=(b[ok]/s) / (l[ok]/s), shape=(b[ok]/s)^2)
  kill_ok <- ok & (k > 0)
  if (any(kill_ok)) {
    nk <- sum(kill_ok)
    t_kill <- rexp(nk, rate=k[kill_ok])
    if (erlang >= 2L) t_kill <- t_kill + rexp(nk, rate=k[kill_ok])  # Erlang-2 = sum of 2 Exp
    kill_idx <- which(kill_ok)
    out[kill_idx[t_kill <= out[kill_idx]]] <- Inf
  }
  out[!ok] <- Inf
  out
}
