# n-choice uniformly varying start point (0-A) Wald
#    race, with t0, v, A, b (boundary) parameterization. 
# ZH added variant with optional between-trial variability following a truncated normal distribution (Steingrover et al 2021 extended to the full Racing framework and paired with start-point variability)

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
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10))
      log_likelihood_race(pars=pars, dadm = dadm, model = model, min_ll = min_ll)
  )
}

## Wrapper functions call C++ functions for likelihood estimation
dRDM <- function(rt,pars)
  # density for single accumulator
{
  if (is.null(dim(pars)) || (dim(pars)[1]==1 & length(rt)>1) ) { # Check if pars is a vector
    original_names <- names(pars); if (is.null(original_names)) {original_names = colnames(pars)}
    pars <- matrix(pars, nrow = length(rt), ncol=length(pars), dimnames = list(NULL, original_names),byrow=TRUE)
  }
  out=rep(NaN, length(rt))
  # Applying drop=FALSE for subsetting pars for robustness, though t0 and v are single columns
  ok <- rt > pars[,"t0",drop=FALSE] & !pars[,"v",drop=FALSE] < 0  # code handles rate zero case
  ok[is.na(ok)] <- FALSE
  if (any(ok)){
    if (any(dimnames(pars)[[2]]=="s")) { # rescale
      # Ensure pars[ok,] remains a matrix even if sum(ok)==1
      pars_ok <- pars[ok,,drop=FALSE]
      pars_ok[,c("A","B","v")] <- pars_ok[,c("A","B","v")]/pars_ok[,"s"]
      pars[ok,] <- pars_ok
    }
    out[ok] <- dRDM_c(rt[ok],v=pars[ok,"v",drop=FALSE],B=pars[ok,"B",drop=FALSE],A=pars[ok,"A",drop=FALSE],t0=pars[ok,"t0",drop=FALSE])
  }
  out
}

pRDM <- function(rt,pars)
  # cumulative density for single accumulator
{
  if (is.null(dim(pars)) || (dim(pars)[1]==1 & length(rt)>1) ) { # Check if pars is a vector
    original_names <- names(pars); if (is.null(original_names)) {original_names = colnames(pars)}
    pars <- matrix(pars, nrow = length(rt), ncol=length(pars), dimnames = list(NULL, original_names),byrow=TRUE)
  }
  out=rep(NaN, length(rt))
  ok <- rt > pars[,"t0",drop=FALSE] & !pars[,"v",drop=FALSE] < 0  # code handles rate zero case
  ok[is.na(ok)] <- FALSE
  if (any(ok)){
    if (any(dimnames(pars)[[2]]=="s")) { # rescale
      pars_ok <- pars[ok,,drop=FALSE]
      pars_ok[,c("A","B","v")] <- pars_ok[,c("A","B","v")]/pars_ok[,"s"]
      pars[ok,] <- pars_ok
    }
    out[ok] <- pRDM_c(rt[ok],v=pars[ok,"v",drop=FALSE],B=pars[ok,"B",drop=FALSE],A=pars[ok,"A",drop=FALSE],t0=pars[ok,"t0",drop=FALSE])
  }
  out
}

#### random
rWald <- function(n,B,v,A)
  # random function for single accumulator
{
  if (n>1 & all(length(A)==1,length(v)==1,length(B)==1)) {
    A=rep(A,n)
    B=rep(B,n)
    v=rep(v,n)
  }
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
  if (is.null(dim(pars)) || (dim(pars)[1]==1 & length(lR)>1) ) { # Check if pars is a vector
    original_names <- names(pars); if (is.null(original_names)) {original_names = colnames(pars)}
    pars <- matrix(pars, nrow = length(lR), ncol=length(pars), dimnames = list(NULL, original_names),byrow=TRUE)
  }
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
  pars <- pars[ok,,drop=FALSE]
  dt[ok] <- rWald(sum(ok),B=pars[,"B"],v=pars[,"v"],A=pars[,"A"])
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
  rt <- matrix(t0,nrow=nr)[pick] + dt[pick]
  out$R <- levels(lR)[R]
  out$R <- factor(out$R,levels=levels(lR))
  out$rt <- rt
  out
}

#' The Racing Diffusion Model with Trial-Varying Drift (DSWTN)
#'
#' Model file to estimate a Racing Diffusion Model where each accumulator
#' follows a Wald distribution with trial-varying drift rates, modeled
#' by the DSWTN (Distribution Scalar Wald Truncated Normal) distribution.
#'
#' @details
#' This model replaces the standard RDM's single drift rate `v` and start-point
#' variability `A` with parameters for a distribution of drift rates:
#' `mu_drift` (mean of the drift distribution) and `sigma_drift_sq` (variance
#' of the drift distribution).
#'
#' Default values are used for all parameters that are not explicitly listed in the `formula`
#' argument of `design()`.
#'
#'#' Default values are used for all parameters that are not explicitly listed in the `formula`
#' argument of `design()`.They can also be accessed with `RDM()$p_types`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default**   | **Mapping**          | **Interpretation**                                                |
#' |-----------|-----------|--------------  -|-----------|------------------|---------------------------------------------------------------|
#' | *v*       | log       | \[0, Inf\]      | log(1)    |                  | Evidence-accumulation rate (drift rate)                        |
#' | *zA*       | qnorm     | \[0, 1\]        | log(0)    |                  | Between-trial variation (range) in start point as a proportion of b (threshold)                 |
#' | *b*       | log       | \[0, Inf\]      | log(1)    |                  | Response threshold                  |
#' | *t0*      | log       | \[0, Inf\]      | log(0)    |                  | Non-decision time                                             |
#' | *s*       | log       | \[0, Inf\]      | log(1)    |                  | Within-trial standard deviation of drift rate                 |
#' | *cv*      | log       | \[0, Inf\]      | log(0.1)  |                  | Standard deviation of the trial-by-trial drift rate distribution, expressed as a coefficient of variation in terms of the corresponding drift rate   |
#'
#' All parameters are typically estimated on the log scale. `A` (start-point variability range)
#' and `sv` (drift rate between trial variance) can be zero.
#' * If `zA=0` and `cv=0`, the model reduces to a simple Wald distribution with threshold `B` and drift `v`.
#' * If `cv=0` (and `zA > 0`), the model reduces to the standard RDM (with start-point variability `A` and fixed drift `v`).
#' * If `zA=0` (and `cv > 0`), the model reduces to a SWTN with a fixed threshold `B`.
#' `b` should generally be positive, especially if `zA=0`.
#'
#' Like the standard RDM, this is a race model with one SWTN-spv accumulator per response option.

#' @export
#' 
RDMSWTN <- function(){
  list(
    type="RACE",
    c_name = "RDMSWTN",
    p_types=c("v" = log(1),"b" = log(1),"zA" = qnorm(0),"t0" = log(0),"s" = log(1),"cv" = qnorm(0), "pContaminant"=qnorm(0)),
    transform=list(func=c(v = "exp", b = "exp", zA = "pnorm",t0 = "exp", s = "exp", cv="pnorm",pContaminant="pnorm")),
    bound=list(minmax=cbind(v=c(1e-3,Inf), b=c(0,Inf), zA=c(0.01,0.99),t0=c(0.05,Inf), s=c(0,Inf), cv=c(0.01,0.99),pContaminant=c(0.001,0.999)),
               exception=c(zA=0, v=0, cv=0,pContaminant=0)),
    # Trial dependent parameter transform. sv is sampled as a coefficient of variance and transformed to standard deviation of drift, tying its magnitude to the mean_drift.
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,A=pars[,"b"] * pars[,"zA"],sv=pars[,"cv"]*pars[,"v"])
      pars
    },
    # Random function for racing accumulators
    rfun=function(data=NULL,pars)  rRDMSWTN(data$lR,pars,ok=attr(pars, "ok")),
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dRDMSWTN(rt,pars),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pRDMSWTN(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10))
      log_likelihood_race(pars=pars, dadm = dadm, model = model, min_ll = min_ll)
  )
}

## The equivalents to dWald and pWald are implemented in C++ and then wrap to their equivalent likelihoods
dRDMSWTN <- function(rt,pars)
  # density for single accumulator
{
  if (is.null(dim(pars)) || (dim(pars)[1]==1 & length(rt)>1) ) { # Check if pars is a vector
    original_names <- names(pars); if (is.null(original_names)) {original_names = colnames(pars)}
    pars <- matrix(pars, nrow = length(rt), ncol=length(pars), dimnames = list(NULL, original_names),byrow=TRUE)
  }
  out=rep(NaN, length(rt))
  # Applying drop=FALSE for subsetting pars for robustness, though t0 and v are single columns
  ok <- rt > pars[,"t0",drop=FALSE] & !pars[,"v",drop=FALSE] < 0  # code handles rate zero case
  ok[is.na(ok)] <- FALSE
  if (any(ok)){
    if (any(dimnames(pars)[[2]]=="s")) { # rescale
      # Ensure pars[ok,] remains a matrix even if sum(ok)==1
      pars_ok <- pars[ok,,drop=FALSE]
      pars_ok[,c("A","b","v","sv")] <- pars_ok[,c("A","b","v","sv")]/pars_ok[,"s"]
      pars[ok,] <- pars_ok
    }
    # dSWTN is a C++ wrapper which takes vectorized input and feeds drdmswtn (the sequential likelihood function in C++)
    out[ok] <- dSWTNspv(rt[ok],v=pars[ok,"v",drop=FALSE],b=pars[ok,"b",drop=FALSE],A=pars[ok,"A",drop=FALSE],t0=pars[ok,"t0",drop=FALSE],sv=pars[ok,"sv",drop=FALSE])
  }
  out
}

pRDMSWTN <- function(rt,pars)
  # cumulative density for single accumulator
{
  if (is.null(dim(pars)) || (dim(pars)[1]==1 & length(rt)>1) ) { # Check if pars is a vector
    original_names <- names(pars); if (is.null(original_names)) {original_names = colnames(pars)}
    pars <- matrix(pars, nrow = length(rt), ncol=length(pars), dimnames = list(NULL, original_names),byrow=TRUE)
  }
  out=rep(NaN, length(rt))
  ok <- rt > pars[,"t0",drop=FALSE] & !pars[,"v",drop=FALSE] < 0  # code handles rate zero case
  ok[is.na(ok)] <- FALSE
  if (any(ok)){
    if (any(dimnames(pars)[[2]]=="s")) { # rescale
      pars_ok <- pars[ok,,drop=FALSE]
      pars_ok[,c("A","b","v","sv")] <- pars_ok[,c("A","b","v","sv")]/pars_ok[,"s"]
      pars[ok,] <- pars_ok
    }
    # pSWTN is a C++ wrapper which takes vectorized input and feeds prdmswtn (the sequential likelihood function in C++)
    out[ok] <- pSWTNspv(rt[ok],v=pars[ok,"v",drop=FALSE],b=pars[ok,"b",drop=FALSE],A=pars[ok,"A",drop=FALSE],t0=pars[ok,"t0",drop=FALSE],sv=pars[ok,"sv",drop=FALSE])
  }
  out
}


rRDMSWTN <- function(lR,pars,p_types=c("v","b","A","t0","sv"),ok=rep(TRUE,dim(pars)[1])) 
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in
  # contiguous rows. "s" parameter will be used but can be ommitted
  #
  # test
  # pars=cbind(B=c(1,2),v=c(1,1),A=c(0,0),t0=c(.2,.2)); lR=factor(c(1,2))
{
  if (!is.null(attr(pars,"ok"))) {ok=attr(pars,"ok")}
  if (is.null(dim(pars)) || (dim(pars)[1]==1 & length(lR)>1) ) { # Check if pars is a vector
    original_names <- names(pars); if (is.null(original_names)) {original_names = colnames(pars)}
    pars <- matrix(pars, nrow = length(lR), ncol=length(pars), dimnames = list(NULL, original_names),byrow=TRUE)
  }
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  if (any(dimnames(pars)[[2]]=="s")) # rescale
    pars[,c("A","b","v","sv")] <- pars[,c("A","b","v","sv")]/pars[,"s"]
  pars[,"b"][pars[,"b"]<0] <- 0 # Protection for negatives
  pars[,"A"][pars[,"A"]<0] <- 0
  bad <- rep(NA, length(lR)/length(levels(lR)))
  out <- data.frame(R = bad, rt = bad)
  nr <- length(levels(lR))
  dt <- matrix(Inf,nrow=nr,ncol=nrow(pars)/nr)
  t0 <- pars[,"t0"]
  pars <- pars[ok,,drop=FALSE]
  dt[ok] <- rSWTN(sum(ok),b=pars[,"b"],v=pars[,"v"],A=pars[,"A"],sv=pars[,"sv"])
  
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
  rt <- matrix(t0,nrow=nr)[pick] + dt[pick]
  out$R <- levels(lR)[R]
  out$R <- factor(out$R,levels=levels(lR))
  out$rt <- rt
  out
}

rSWTN <- function(n,b,v,A,sv,s=1)
  # random function for single accumulator
{
  
  out=numeric(n)
  if (n>1 & all(length(A)==1,length(v)==1,length(b)==1,length(sv==1))) {
    A=rep(A,n)
    b=rep(b,n)
    v=rep(v,n)
    sv=rep(sv,n)
  }
  b = ifelse(A==0,b,runif(n,b-A, b)) # adjust for spv U[0,A]
  l = ifelse(sv==0,v,msm::rtnorm(n,mean=v,sd=sv,lower=0,upper=Inf)) # between trial variability
  
  ok <- !l<0
  nok <- sum(ok)
  out[ok] <- statmod::rinvgauss(nok, mean = (b / s) / (l / s), shape = (b / s)^2)
  out[!ok] <- Inf
  out  
  
}

#' @export
#'

LogicalRulesRDMSWTN <- function(){
  list(
    type="RACE",
    c_name = "RDMSWTN_LogicalRules",
    # Trial dependent parameter transform. sv is sampled as a coefficient of variance and transformed to standard deviation of drift, tying its magnitude to the mean_drift.
    p_types=c("v" = log(1),"b" = log(1),"zA" = qnorm(0),"t0" = log(0),"s" = log(1),"cv" = qnorm(0), "p"=qnorm(1),"q"=qnorm(0.5), "r"=qnorm(1)),
    transform=list(func=c(v = "exp", b = "exp", zA = "pnorm",t0 = "exp", s = "exp", cv="pnorm",p="pnorm",q="pnorm", r="pnorm")),
    bound=list(minmax=cbind(v=c(1e-3,Inf), b=c(0,Inf), zA=c(0.01,0.99),t0=c(0.05,Inf), s=c(0,Inf), cv=c(0.01,0.99),p=c(0.01,0.99),q=c(0.01,0.99),r=c(0.01,0.99)),
               exception=c(zA=0, v=0, cv=0,p=1,r=0)),
    # Transform to natural scale
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,A=pars[,"b"] * pars[,"zA"],sv=pars[,"cv"]*pars[,"v"])
      pars
    },
    # Random function for racing accumulators
    rfun=function(data=NULL,pars)  rRDMSWTN(data$lR,pars,ok=attr(pars, "ok")),
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dRDMSWTN(rt,pars),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pRDMSWTN(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
      log_likelihood_redundant_target_race(pars=pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}
