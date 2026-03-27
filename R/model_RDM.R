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
RDM <- function(){
  list(
    type="RACE",
    c_name = "RDM",
    p_types=c("v" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),"s" = log(1),"sv" = log(0), "pContaminant"=qnorm(0)),
    transform=list(func=c(v = "exp", B = "exp", A = "exp",t0 = "exp", s = "exp", sv="exp",pContaminant="pnorm")),
    bound=list(minmax=cbind(v=c(1e-3,Inf), B=c(0,Inf), A=c(0,Inf),t0=c(0.05,Inf), s=c(0,Inf), sv=c(0,Inf),pContaminant=c(0.001,0.999)),
               exception=c(A=0, v=0, sv=0,pContaminant=0)),
    # Trial dependent parameter transform.
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
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
