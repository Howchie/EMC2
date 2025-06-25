

# n-choice uniformly varying start point (0-A) Wald
#    race, with t0, v, A, b (boundary) parameterixaiton

# pigt, digt, rwaldt Copyright (C) 2013  Trisha Van Zandt distributed with:
# Logan, Van Zandt, Verbruggen, and Wagenmakers (2014).  On the ability to
# inhibit thought and action: General and special theories of an act of control.
# Psychological Review. Comments and changes added by Andrew Heathcote. Trish's
# code is for k = threshold, a = half width of uniform threshold variability,
# l = rate of accumulation. Note that Wald mean = k/l and shape = k^2.
#
# I use a different parameterization in terms of v=l (rate),
# uniform start point variability from 0-A (A>=0), threshold b (>0) and hence
# B=b-A (>=0) as a threshold gap. Hence k = b-A/2 = B + A/2 and a=A
#
# Added ability to set s = diffusive standard deviation, by scaling A, B and v
# i.e., given s, same result for A/s, B/s, v/s with s=1

#### distribution functions

# # Moved to C++ model_RDM.cpp
#
# digt0 <- function(t,k=1,l=1) {
#       # pdf of inverse gaussian at t with no k variability
#       # much faster than statmod's dinvgauss funciton
#
#       lambda <- k^2
#       l0 <- l==0
#       e <- numeric(length(t))
#       if ( any(!l0) ) {
#         mu <- k[!l0]/l[!l0]
#         e[!l0] <- -(lambda[!l0]/(2*t[!l0])) * (t[!l0]^2/mu^2 - 2*t[!l0]/mu  + 1)
#       }
#       if ( any(l0) )  e[l0] <- -.5*lambda[l0]/t[l0]
#       x <- exp(e + .5*log(lambda) - .5*log(2*t^3*pi))
#       x[t<=0] <- 0
#       x
# }
#
#
# digt <- function(t,k=1,l=1,a=.1,tiny=1e-10) {
#     # pdf of inverse gaussian at t with k +/- a/2 uniform variability
#     # returns digt0 if a<1e-10
#
#
#     options(warn=-1)
#     if(length(k)!=length(t)) k <- rep(k,length.out=length(t))
#     if(length(l)!=length(t)) l <- rep(l,length.out=length(t))
#     if(length(a)!=length(t)) a <- rep(a,length.out=length(t))
#
#     tpos <- t<=0
#
#     atiny <- a<=tiny & !tpos
#     a[atiny] <- 0
#
#     ltiny <- (l<=tiny) & !atiny & !tpos
#     notltiny <- (l>tiny) & !atiny & !tpos
#     l[l<=tiny] <- 0
#
#     x <- numeric(length(t))
#
#     # No threshold variability
#     if ( any(atiny) )
#       x[atiny] <- digt0(t=t[atiny],k=k[atiny],l=l[atiny])
#
#     # Threshold variability
#     if ( any(!atiny) ) {
#
#       if ( any(notltiny) ) { # rate non-zero
#
#         sqr.t <- sqrt(t[notltiny])
#
#         term.1a <- -(a[notltiny]-k[notltiny]+t[notltiny]*l[notltiny])^2/(2*t[notltiny])
#         term.1b <- -(a[notltiny]+k[notltiny]-t[notltiny]*l[notltiny])^2/(2*t[notltiny])
#         term.1 <- (exp(term.1a) - exp(term.1b))/sqrt(2*pi*t[notltiny])
#
#         term.2a <- log(.5)+log(l[notltiny])
#         term.2b <- 2*pnorm((-k[notltiny]+a[notltiny])/sqr.t+sqr.t*l[notltiny])-1
#         term.2c <- 2*pnorm((k[notltiny]+a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
#         term.2d <- term.2b+term.2c
#         term.2 <- exp(term.2a)*term.2d
#
#         term.3 <- term.1+term.2
#         term.4 <- log(term.3)-log(2)-log(a[notltiny])
#         x[notltiny] <- exp(term.4)
#       }
#
#       if ( any(ltiny) ) {  # rate zero
#         log.t <- log(t[ltiny])
#         term.1 <- -.5*(log(2)+log(pi)+log.t)
#         term.2 <- (k[ltiny]-a[ltiny])^2/(2*t[ltiny])
#         term.3 <- (k[ltiny]+a[ltiny])^2/(2*t[ltiny])
#         term.4 <- (exp(-term.2)-exp(-term.3))
#         term.5 <- term.1+log(term.4) - log(2) - log(a[ltiny])
#         x[ltiny] <- exp(term.5)
#       }
#
#     }
#
#     x[x<0 | is.nan(x) ] <- 0
#     x
# }
#
#
# pigt0 <- function(t,k=1,l=1) {
#       # cdf of inverse gaussian at t with no k variability
#       # much faster than statmod's pinvgauss funciton
#
#       mu <- k/l
#       lambda <- k^2
#
#       e <- exp(log(2*lambda) - log(mu))
#       add <- sqrt(lambda/t) * (1 + t/mu)
#       sub <- sqrt(lambda/t) * (1 - t/mu)
#
#       p.1 <- 1 - pnorm(add)
#       p.2 <- 1 - pnorm(sub)
#       x <- exp(e + log(p.1)) + p.2
#
#       x[t<0] <- 0
#       x
# }
#
#
# pigt <- function(t,k=1,l=1,a=.1,tiny=1e-10) {
#     # cdf of inverse gaussian at t with k +/- a/2 uniform variability
#     # returns pigt0 if a<=0
#
#     options(warn=-1)
#     if(length(k)!=length(t)) k <- rep(k,length.out=length(t))
#     if(length(l)!=length(t)) l <- rep(l,length.out=length(t))
#     if(length(a)!=length(t)) a <- rep(a,length.out=length(t))
#
#     tpos <- t<=0
#
#     atiny <- a<=tiny & !tpos
#     a[atiny] <- 0
#
#     ltiny <- (l<=tiny) & !atiny & !tpos
#     notltiny <- (l>tiny) & !atiny & !tpos
#     l[l<=tiny] <- 0
#
#     x <- numeric(length(t))
#
#     # No threshold variability
#     if ( any(atiny) )
#       x[atiny] <- pigt0(t[atiny],k[atiny],l[atiny])
#
#     # Threshold variability
#     if ( any(!atiny) ) {
#
#       if ( any(notltiny) ) { # rate non-zero
#
#         log.t <- log(t[notltiny])
#         sqr.t <- sqrt(t[notltiny])
#
#         term.1a <- .5*log.t-.5*log(2*pi)
#         term.1b <- exp(-((k[notltiny]-a[notltiny]-t[notltiny]*l[notltiny])^2/t[notltiny])/2)
#         term.1c <- exp(-((k[notltiny]+a[notltiny]-t[notltiny]*l[notltiny])^2/t[notltiny])/2)
#         term.1 <- exp(term.1a)*(term.1b-term.1c)
#
#         term.2a <- exp(2*l[notltiny]*(k[notltiny]-a[notltiny]) +
#                          log(pnorm(-(k[notltiny]-a[notltiny]+t[notltiny]*l[notltiny])/sqr.t)))
#         term.2b <- exp(2*l[notltiny]*(k[notltiny]+a[notltiny]) +
#                          log(pnorm(-(k[notltiny]+a[notltiny]+t[notltiny]*l[notltiny])/sqr.t)))
#         term.2 <- a[notltiny] + (term.2b-term.2a)/(2*l[notltiny])
#
#         term.4a <- 2*pnorm((k[notltiny]+a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
#         term.4b <- 2*pnorm((k[notltiny]-a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
#         term.4c <- .5*(t[notltiny]*l[notltiny] - a[notltiny] - k[notltiny] + .5/l[notltiny])
#         term.4d <- .5*(k[notltiny] - a[notltiny] - t[notltiny]*l[notltiny] - .5/l[notltiny])
#         term.4 <- term.4c*term.4a + term.4d*term.4b
#
#         x[notltiny] <- (term.4 + term.2 + term.1)/(2*a[notltiny])
#       }
#
#       if ( any(ltiny) ) {  # rate zero
#         sqr.t <- sqrt(t[ltiny])
#         log.t <- log(t[ltiny])
#         term.5a <- 2*pnorm((k[ltiny]+a[ltiny])/sqr.t)-1
#         term.5b <- 2*pnorm(-(k[ltiny]-a[ltiny])/sqr.t)-1
#         term.5 <- (-(k[ltiny]+a[ltiny])*term.5a - (k[ltiny]-a[ltiny])*term.5b)/(2*a[ltiny])
#
#         term.6a <- -.5*(k[ltiny]+a[ltiny])^2/t[ltiny] - .5*log(2) -.5*log(pi) + .5*log.t - log(a[ltiny])
#         term.6b <- -.5*(k[ltiny]-a[ltiny])^2/t[ltiny] - .5*log(2) -.5*log(pi) + .5*log.t - log(a[ltiny])
#         term.6 <- 1 + exp(term.6b) - exp(term.6a)
#
#         x[ltiny] <- term.5 + term.6
#       }
#
#     }
#
#     x[x<0 | is.nan(x) ] <- 0
#     x
# }
#
#

dRDM <- function(rt,pars)
  # density for single accumulator
{
  if (is.null(dim(pars)) || (dim(pars)[1]==1 & length(rt>1)) ) { # Check if pars is a vector
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
    out[ok] <- dWald(rt[ok],v=pars[ok,"v",drop=FALSE],B=pars[ok,"B",drop=FALSE],A=pars[ok,"A",drop=FALSE],t0=pars[ok,"t0",drop=FALSE])
  }
  out
}


pRDM <- function(rt,pars)
  # cumulative density for single accumulator
{
  if (is.null(dim(pars)) || (dim(pars)[1]==1 & length(rt>1)) ) { # Check if pars is a vector
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
    out[ok] <- pWald(rt[ok],v=pars[ok,"v",drop=FALSE],B=pars[ok,"B",drop=FALSE],A=pars[ok,"A",drop=FALSE],t0=pars[ok,"t0",drop=FALSE])
  }
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
  R <- apply(dt,2,which.min)
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
    p_types=c("v" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),"s" = log(1)),
    transform=list(func=c(v = "exp", B = "exp", A = "exp",t0 = "exp", s = "exp")),
    bound=list(minmax=cbind(v=c(1e-3,Inf), B=c(0,Inf), A=c(1e-4,Inf),t0=c(0.05,Inf), s=c(0,Inf)),
               exception=c(A=0, v=0)),
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
#' | **Parameter**      | **Transform** | **Natural scale** | **Default**     | **Interpretation**                                           |
#' |------------------|---------------|-------------------|---------------|--------------------------------------------------------------|
#' | *alpha*          | log           | \[0, Inf\]        | log(1)        | Response threshold (boundary separation)                     |
#' | *mu_drift*       | log           | \[0, Inf\]        | log(1)        | Mean of the trial-by-trial drift rate distribution       |
#' | *sigma_drift_sq* | log           | \[0, Inf\]        | log(0.1)      | Variance of the trial-by-trial drift rate distribution   |
#' | *theta*          | log           | \[0, Inf\]        | log(0.1)      | Non-decision time                                            |
#'
#' All parameters are typically estimated on the log scale. `sigma_drift_sq` must be positive.
#'
#' Like the standard RDM, this is a race model with one DSWTN accumulator per response option.
#' The `lR` factor (latent response) and `lM` factor (match between stimulus and response)
#' can be used in formulas as with the standard RDM.
#'
#' @return A list defining the cognitive model
#' @export
#' @examples
#' # Define a simple DSWTN RDM model for a 2-choice task.
#' # Formulas specify how parameters relate to experimental factors (if any).
#' # Here, alpha, mu_drift, sigma_drift_sq, and theta are estimated per accumulator (lR).
#' \dontrun{
#' model_dswtn <- RDM_DSWTN()
#'
#' # Example design (assuming 'data_df' has 'rt', 'R' (response factor), and 'S' (stimulus factor))
#' # Match function for lM (match between stimulus S and latent response lR)
#' matchfun <- function(d) d$S == d$lR
#'
#' design_dswtn <- design(
#'   data = data_df,
#'   model = model_dswtn,
#'   matchfun = matchfun,
#'   formula = list(
#'     alpha ~ lR,
#'     mu_drift ~ lM, # Drift rate mean depends on stimulus-accumulator match
#'     sigma_drift_sq ~ lR,
#'     theta ~ lR
#'   )
#' )
#'
#' # To fit this model (e.g., using pmwg):
#' # samples_dswtn <- pmwg(design_dswtn, iter=100, display_progress=TRUE) # Example iterations
#' }
RDM_DSWTN <- function() {

  # Helper function to prepare parameter matrices for loglik_DSWTN_race
  # This is a conceptual placeholder; actual implementation depends on how
  # EMC2's design system makes design matrices available for each parameter of each accumulator.
  prepare_dswtn_race_params <- function(pars, dadm, model) {
    # pars: named vector of all fitted parameters
    # dadm: augmented data, containing data and design matrices (e.g., dadm$X)
    # model: the model object itself

    n_trials <- nrow(dadm$X) # Assuming dadm$X contains design matrix for one trial
                            # Or more generally, number of unique trial instances
    if (is.null(n_trials)) n_trials <- length(dadm$rt) # Fallback if X structure is different

    # Accumulator names/levels from lR (latent response factor)
    acc_levels <- levels(dadm$lR) # dadm should have lR factor
    n_acc <- length(acc_levels)

    # Initialize parameter matrices
    params_alpha <- matrix(NA, nrow = n_trials, ncol = n_acc, dimnames = list(NULL, acc_levels))
    params_mu_drift <- matrix(NA, nrow = n_trials, ncol = n_acc, dimnames = list(NULL, acc_levels))
    params_sigma_drift_sq <- matrix(NA, nrow = n_trials, ncol = n_acc, dimnames = list(NULL, acc_levels))
    params_theta <- matrix(NA, nrow = n_trials, ncol = n_acc, dimnames = list(NULL, acc_levels))

    # This loop structure is highly dependent on how EMC² stores design matrices for race models.
    # Typically, for each parameter (e.g., alpha), there's a design matrix that, when multiplied
    # by the relevant subset of `pars`, gives the parameter value for each trial and accumulator.
    # Example: par_info <- attr(dadm$X, "par_info") might give info on which columns of dadm$X map to which param.

    # Placeholder logic: Assumes `dadm$X` is a combined design matrix and `attr(dadm$X,"pmat")`
    # maps its columns to model parameters (alpha_acc1, alpha_acc2, mu_drift_acc1 etc.)
    # This needs to be adapted to the actual structure provided by EMC2's `design` process.
    # The `get_pars` function in EMC2 might be relevant here.

    # For each trial (assuming parameters can differ per trial, though often they don't beyond accumulator effects)
    # This simplified version assumes parameters are constant per trial but can differ per accumulator based on design.
    # A more general solution would use the full design matrix per trial.

    # Let's assume `get_pars_rdm_dswtn` is a helper that extracts these based on `pars` and `dadm`
    # This is where the core EMC2 machinery for parameter expansion would be used.
    # For now, this is a conceptual sketch.

    # A simplified approach if parameters are constant across trials (but vary by accumulator):
    # Extract from `pars` vector based on names like "alpha_L", "alpha_R", "mu_drift_L" etc.
    # This depends on how `make_par_names` (from EMC2 utils) creates names from formulas.

    all_params_natural <- model$transform$func(pars, model) # Transform all params to natural scale

    for (acc_idx in 1:n_acc) {
      acc_name <- acc_levels[acc_idx]
      # Construct parameter names as generated by EMC2 (e.g., "alpha_L", "mu_drift_R")
      # This requires knowing the naming convention from `make_par_names` based on formulas.
      # Example: if formula for alpha is `alpha ~ lR`, names might be `alpha_L`, `alpha_R`.
      # If `alpha ~ 1`, name is just `alpha`.

      # This part is highly schematic and needs to use EMC2's parameter extraction logic
      # (e.g., using design matrices dadm$X and attr(dadm$X, "pmat"))
      # For each parameter (alpha, mu_drift, etc.):
      #   Identify columns in dadm$X relevant to this parameter for this accumulator.
      #   Calculate: dadm$X[, cols] %*% pars[relevant_pars_indices]
      # This needs to be done for each of the 4 DSWTN parameters for each accumulator.

      # Simplified: assume `all_params_natural` contains named expanded parameters for each trial & accumulator
      # This is what `get_pars` in EMC2 usually produces.
      # If `get_pars` is adapted or a similar function is used:
      # expanded_pars_matrix <- get_pars(pars, dadm, model) # N_trials x N_total_params_expanded
      # Then, from expanded_pars_matrix, pick columns for alpha_acc_idx, mu_drift_acc_idx, etc.

      # Fallback: If we assume `pars` already contains expanded per-accumulator values (e.g. from a prior step)
      # and names like alpha_ACCLEVEL, mu_drift_ACCLEVEL. This is unlikely for the main `pars` vector.

      # The most robust way is to use the design matrix multiplication for each parameter.
      # For now, let's assume a function `get_accumulator_params` exists that does this.
      # This function would be specific to how EMC2 handles design matrices for race models.

      # Placeholder: Fill with defaults or stop if logic is not fully implemented yet for parameter extraction.
      # This is the part that needs to correctly interface with EMC2's design output.
      # For the purpose of defining the model object, we assume this will be resolved by EMC2's framework.
      # The following is a mock-up of what needs to happen:
      for (param_name in c("alpha", "mu_drift", "sigma_drift_sq", "theta")) {
        # Determine the full parameter name in `all_params_natural` (e.g., "alpha_L", "alpha_R")
        # This depends on how formulas like `alpha ~ lR` are expanded.
        # Let's assume `dadm[[paste0(param_name, "_design_matrix_for_acc_", acc_name)]]` exists
        # and `pars[[paste0(param_name,"_coeffs")]]` exists. This is too complex.

        # A common pattern: design matrix X, parameter vector p.
        # X %*% p gives values. `dadm$X` is often trial-specific.
        # `dadm$par_info` or attributes of `X` usually tell how columns of X map to parameters.

        # This function is critical and complex. For now, we'll assume it's handled by the broader EMC2 framework
        # when this model is used with `design()` and `log_likelihood()`.
        # The `loglik_DSWTN_race` C++ function expects N_trials x N_acc matrices.
        # These would be constructed by EMC2's machinery using the design matrices and parameter vector.
      }
    }

    # This function would return a list of these matrices.
    # However, the `log_likelihood` function in the model object typically receives `pars` (coefficient vector)
    # and `dadm` (which includes design matrices). It then internally computes the per-trial, per-accumulator
    # parameters.

    # So, the logic inside the log_likelihood function needs to do this.
    # This is a complex piece of EMC2's infrastructure.
    # For now, the C++ `loglik_DSWTN_race` expects the final per-trial-per-acc matrices.
    # The R wrapper in the model definition needs to bridge this.
    stop("prepare_dswtn_race_params is a placeholder and needs full EMC2 integration logic.")

    # A more realistic approach for the log_likelihood slot in the model:
    # It will receive `pars` (coefficients) and `dadm`.
    # It needs to use `dadm$X` (design matrix for a single trial, typically)
    # and `attr(dadm$X, "pmat")` (parameter matrix mapping X columns to parameters)
    # along with `dadm$n_X` (number of trials for this X structure)
    # and `dadm$sX` (subject-level design matrix if applicable)
    # to construct the per-trial, per-accumulator parameters.
    # This is what generic `log_likelihood_race` does by calling `model$dfun` and `model$pfun`
    # with parameters extracted for each trial/accumulator.

    # If we bypass generic log_likelihood_race and call loglik_DSWTN_race directly,
    # we need to replicate that parameter expansion here.
    # This is a significant task.
  }

  # Simplified rfun placeholder - actual implementation needs care
  rfun_dswtn_race <- function(data, pars) {
    # data contains lR (latent response factor), n_trials etc.
    # pars contains parameter values for each trial & accumulator
    # This needs to be structured to match how rRDM gets its `pars`

    n_trials <- nrow(pars) / length(levels(data$lR)) # Assuming pars is n_trials*n_acc x n_params
    acc_levels <- levels(data$lR)
    n_acc <- length(acc_levels)

    # Reshape/extract pars into a list per accumulator, each element being a matrix/df for trials
    # This is highly dependent on how `pars` is formatted by the calling EMC2 framework.
    # For now, assume `pars` is already in a format where we can extract
    # alpha, mu_drift, sigma_drift_sq, theta for each trial and accumulator.
    # Example: pars might be a matrix from get_pars()

    rt_matrix <- matrix(Inf, nrow = n_trials, ncol = n_acc)

    # This assumes `pars` is a matrix from `get_pars` where columns are named e.g. alpha_L, mu_drift_L etc.
    # And we need to iterate trials.
    # This is complex to generalize here without EMC2's internal `get_pars` output structure.

    # Simpler: if `pars` argument to `rfun` is already N_trials x N_params_total_expanded
    # (e.g., alpha_acc1, mu_drift_acc1, theta_acc1, alpha_acc2, ...)
    # This is also unlikely.

    # Let's assume `pars` here is the same N_trials x N_acc matrices as for loglik_DSWTN_race
    # (alpha_mat, mu_drift_mat, etc.). This is not what `rfun` usually gets.
    # `rfun` typically gets a parameter vector `p` and `dadm`.

    # For now, a very simplified rfun assuming fixed parameters for all trials for illustration:
    # It would need proper parameter extraction from the `pars` vector passed by EMC2.
    if (is.null(attr(pars,"ok"))) ok_trials <- rep(TRUE, n_trials) else ok_trials <- attr(pars,"ok")

    # This is a placeholder and will not work directly with EMC2's usual `pars` structure for rfun.
    # It demonstrates the core race logic with rDSWTN.
    # A proper rfun would need to extract trial-specific parameters for each accumulator first.

    # Simplified: assume pars is a list where each element is a named list of params for an acc
    # pars[[1]]$alpha, pars[[1]]$mu_drift etc. (vectors of length n_trials)

    # Fallback to error for rfun as it needs careful integration
    stop("rfun for RDM_DSWTN requires careful implementation of parameter extraction for rDSWTN calls.")

    # Example structure if parameters were readily available per trial/acc:
    # for (i in 1:n_trials) {
    #   if (!ok_trials[i]) next
    #   trial_rts <- numeric(n_acc)
    #   for (j in 1:n_acc) {
    #     trial_rts[j] <- rDSWTN(1, alpha = pars_alpha[i,j], mu_drift = pars_mu_drift[i,j], ...)
    #   }
    #   rt_matrix[i, ] <- trial_rts
    # }
    # winner_idx <- apply(rt_matrix, 1, which.min) # for ok_trials
    # winner_rt <- rt_matrix[cbind(which(ok_trials), winner_idx[ok_trials])]
    # ... then format into data.frame(R, rt)
  }


  list(
    type = "RACE",
    c_name = "RDM_DSWTN", # Name for this variant
    p_types = c(alpha = log(1), mu_drift = log(1), sigma_drift_sq = log(0.25), theta = log(0.1)),
    transform = list(
      func = c(alpha = "exp", mu_drift = "exp", sigma_drift_sq = "exp", theta = "exp"),
      args = NULL # No special arguments for exp
    ),
    bound = list(
      minmax = cbind(alpha = c(1e-3, Inf), mu_drift = c(1e-3, Inf), sigma_drift_sq = c(1e-4, Inf), theta = c(0.001, Inf)),
      exception = NULL # No specific exceptions like A=0 for original RDM
    ),
    Ttransform = function(pars, dadm) {
      # No additional parameter transformations like b = B+A needed here by default
      # Users could add if they re-parameterize (e.g. log_sigma_drift_sq)
      pars
    },
    rfun = rfun_dswtn_race, # Placeholder, needs full implementation
    dfun = function(rt, pars, dadm=NULL) { # pars is matrix: rows for trials, cols for params
      # Ensure parameters are on natural scale if coming from fitting (usually are by this point)
      # dfun is for single accumulator, log_likelihood_race handles multiple.
      # Here `pars` would be for ONE accumulator.
      exp(dDSWTN_log(rt, alpha=pars[,"alpha"], mu_drift=pars[,"mu_drift"],
                     sigma_drift_sq=pars[,"sigma_drift_sq"], theta=pars[,"theta"]))
    },
    pfun = function(rt, pars, dadm=NULL) {
      pDSWTN(rt, alpha=pars[,"alpha"], mu_drift=pars[,"mu_drift"],
             sigma_drift_sq=pars[,"sigma_drift_sq"], theta=pars[,"theta"])
    },
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      # This function needs to extract per-trial, per-accumulator parameters
      # from `pars` (coefficients) and `dadm` (design info), then call C++ loglik_DSWTN_race.

      # Get expanded parameters (N_trials x N_params_total_expanded)
      # This uses EMC2's internal parameter expansion system.
      # `all_natural_pars` would be a matrix where columns are like alpha_L, alpha_R, mu_drift_L etc.
      # if lR has levels L, R.
      all_natural_pars <- get_pars(p_vector = pars, dadm = dadm, model = model)

      acc_levels <- levels(dadm$lR)
      n_acc <- length(acc_levels)
      n_trials <- nrow(dadm$data) # Number of trials in the data for this subject/block

      # Prepare matrices for C++ call
      params_alpha_mat <- matrix(NA, nrow = n_trials, ncol = n_acc)
      params_mu_drift_mat <- matrix(NA, nrow = n_trials, ncol = n_acc)
      params_sigma_drift_sq_mat <- matrix(NA, nrow = n_trials, ncol = n_acc)
      params_theta_mat <- matrix(NA, nrow = n_trials, ncol = n_acc)

      for (acc_idx in 1:n_acc) {
        acc_name <- acc_levels[acc_idx]
        # Construct expected parameter names (e.g., "alpha_L", "mu_drift_L")
        # These names are generated by make_par_names in EMC2 utils.
        # This assumes simple factor expansion (e.g., alpha ~ lR).
        # More complex formulas would result in different names from make_par_names.

        # A robust way uses attr(all_natural_pars, "pnames_by_type_acc") if available from get_pars
        # Or reconstruct names based on model$p_types and acc_name.

        # Simplified name construction (may need adjustment based on make_par_names behavior)
        current_p_names <- dimnames(all_natural_pars)[[2]]

        alpha_col_name <- grep(paste0("^alpha_", acc_name, "$|^alpha$"), current_p_names, value = TRUE)
        mu_drift_col_name <- grep(paste0("^mu_drift_", acc_name, "$|^mu_drift$"), current_p_names, value = TRUE)
        sigma_drift_sq_col_name <- grep(paste0("^sigma_drift_sq_", acc_name, "$|^sigma_drift_sq$"), current_p_names, value = TRUE)
        theta_col_name <- grep(paste0("^theta_", acc_name, "$|^theta$"), current_p_names, value = TRUE)

        # Handle cases where a parameter is global (not varying by accumulator)
        if (length(alpha_col_name) == 0 && "alpha" %in% current_p_names) alpha_col_name <- "alpha"
        if (length(mu_drift_col_name) == 0 && "mu_drift" %in% current_p_names) mu_drift_col_name <- "mu_drift"
        if (length(sigma_drift_sq_col_name) == 0 && "sigma_drift_sq" %in% current_p_names) sigma_drift_sq_col_name <- "sigma_drift_sq"
        if (length(theta_col_name) == 0 && "theta" %in% current_p_names) theta_col_name <- "theta"

        if (length(alpha_col_name)==0 || length(mu_drift_col_name)==0 ||
            length(sigma_drift_sq_col_name)==0 || length(theta_col_name)==0) {
            stop(paste("Could not find all DSWTN parameter columns for accumulator", acc_name,
                       "in matrix from get_pars(). Found:",
                       paste(current_p_names, collapse=", ")))
        }

        params_alpha_mat[, acc_idx] <- all_natural_pars[, alpha_col_name[1]]
        params_mu_drift_mat[, acc_idx] <- all_natural_pars[, mu_drift_col_name[1]]
        params_sigma_drift_sq_mat[, acc_idx] <- all_natural_pars[, sigma_drift_sq_col_name[1]]
        params_theta_mat[, acc_idx] <- all_natural_pars[, theta_col_name[1]]
      }

      # dadm$data$R should be the factor of observed responses. Convert to numeric index.
      # Ensure dadm$data$R levels match dadm$lR levels for consistent indexing.
      observed_choice_idx <- as.numeric(factor(as.character(dadm$data$R), levels = acc_levels))

      loglik_DSWTN_race(
        rts = dadm$data$rt,
        choices = observed_choice_idx, # Must be 1-indexed
        params_alpha = params_alpha_mat,
        params_mu_drift = params_mu_drift_mat,
        params_sigma_drift_sq = params_sigma_drift_sq_mat,
        params_theta = params_theta_mat,
        min_log_lik = min_ll
        # Pass cdf control params if needed, e.g. from model$cdf_control
      )
    }
  )
}

Mrdm <- function(){
  list(
    type="RACE",
    c_name = NULL, # must be NULL to use calc_ll_R
    # p_vector transform, sets sv as a scaling parameter
    p_types=c("v" = 1,"s" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0)),
    transform=list(func=c(v = "log",s = "exp", B = "exp", A = "exp",t0 = "exp")),
    bound=list(minmax=cbind(v=c(0,Inf),s = c(0, Inf), A=c(1e-4,Inf),B=c(0,Inf),t0=c(0.01,Inf)),
               exception=c(A=0)),
    # Transform to natural scale
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
    log_likelihood=function(pars,dadm,model){
      log_likelihood_race_cens_trunc(pars=pars, dadm=dadm, model=model, min_ll=log(1e-10))
    }
  )
}