#### Staircase ----
#' Assign Stop-Signal Delays (SSDs) to trials
#'
#' @description
#' This function assigns stop-signal delays (SSDs) to trials based on specified probabilities.
#'
#' @param d A data frame containing trial data with a response factor column 'lR'
#' @param SSD A vector of stop-signal delays to be assigned to trials
#' @param pSSD A vector of probabilities for each SSD value. If length is one less than SSD,
#'             the remaining probability is calculated automatically. Default is 0.25.
#'
#' @return A vector of SSDs with the same length as the number of rows in 'd'.
#'         Trials without a stop signal are assigned Inf.
SSD_function <- function(d,SSD=NA,pSSD=.25) {
  if (sum(pSSD)>1) stop("pSSD sum cannot exceed 1.")
  if (length(pSSD)==length(SSD)-1) pSSD <- c(pSSD,1-sum(pSSD))
  if (length(pSSD)!=length(SSD))
    stop("pSSD must be the same length or 1 less than the length of SSD")
  n_acc <- length(levels(d$lR))
  n_trial <- nrow(d)
  out <- rep(Inf,n_trial)
  trials <- c(1:n_trial)
  for (i in 1:length(pSSD)) {
    pick <- sample(trials,floor(pSSD[i]*n_trial))
    out[pick] <- SSD[i]
    trials <- trials[!(trials %in% pick)]
  }
  return(out)
}

staircase_function <- function(dts,staircase) {
  ns <- ncol(dts)
  SSD <- sR <- srt <- numeric()
  SSD[1] <- staircase$SSD0
  for (i in 1:ns) {
    if (SSD[i]<staircase$stairmin) SSD[i] <- staircase$stairmin
    if (SSD[i]>staircase$stairmax) SSD[i] <- staircase$stairmax
    dts[1,i] <- dts[1,i] + SSD[i]
    if (all(is.infinite(dts[-1,i]))) Ri <- 1 else # GF or GF & TF
      Ri <- which.min(dts[,i])
    if (Ri==1) {
      sR[i] <- srt[i] <- NA
      if (i<ns) SSD[i+1] <- round(SSD[i] + staircase$stairstep,3)
    } else {
      sR[i] <- Ri-1
      srt[i] <- min(dts[-1,i])
    if (i<ns) SSD[i+1] <- round(SSD[i] - staircase$stairstep,3)
    }
  }
  list(sR=sR,srt=srt,SSD=SSD)
}


check_staircase <- function(staircase){
  if (!is.list(staircase)){
    staircase <- list(SSD0=.25,stairstep=.05,stairmin=0,stairmax=Inf)
  }
  return(staircase)
}




ST_staircase_function <- function(dts,staircase) {
  ns <- ncol(dts)
  SSD <- sR <- srt <- numeric()
  SSD[1] <- staircase$SSD0
  if (is.null(staircase$accST))  # Indices for accumulator with SSD
    iSSD <- 1 else               # Only stop accumulator
    iSSD <- c(1,staircase$accST) # NB: accST refers to rows in dts
  for (i in 1:ns) {
    if (SSD[i]<staircase$stairmin) SSD[i] <- staircase$stairmin
    if (SSD[i]>staircase$stairmax) SSD[i] <- staircase$stairmax
    dts[iSSD,i] <- dts[iSSD,i] + SSD[i]
    if (all(is.infinite(dts[-1,i]))) { # GF (no ST) or GF & TF
      sR[i] <- srt[i] <- NA
      SSD[i+1] <- round(SSD[i] + staircase$stairstep,3)
    } else {
      Ri <- which.min(dts[,i])
      if (Ri==1) {
        if (!is.null(staircase$accST)) { # ST present
          sR[i] <- staircase$accST[which.min(dts[staircase$accST,i])]-1
          srt[i] <- dts[sR[i]+1,i]
        } else sR[i] <- srt[i] <- NA
        if (i<ns) SSD[i+1] <- round(SSD[i] + staircase$stairstep,3)
      } else {
        sR[i] <- Ri-1
        srt[i] <- dts[Ri,i]
        if (i<ns) {
          if (Ri %in% iSSD) # Stop-triggered response
            SSD[i+1] <- round(SSD[i] + staircase$stairstep,3) else
            SSD[i+1] <- round(SSD[i] - staircase$stairstep,3)
        }
      }
    }
  }
  list(sR=sR,srt=srt,SSD=SSD)
}

#### Single exGaussian functions ----

# Following functions moved to C++ model_SS_EXG.cpp

# pEXG <- function (q, mu = 5, sigma = 1, tau = 1, lower_tail = TRUE, log_p = FALSE)
#     # ex-Gaussian cumulative density
#     # Modified from gamlss.dist to make cdf in tau > 0.05 * sigma case robust,
#     # and robust to -Inf and Inf inputs, returns NA for bad sigma or tau and
#     # robust against small sigma cases.
#   {
#     if (sigma <= 0) return(rep(NA,length(q)))
#       if (tau <= 0) return(rep(NA,length(q)))
#
#       # if (sigma < 0.05*tau)
#       if (sigma < 1e-4)
#        return(pexp(q-mu,1/tau,log.p=log_p,lower.tail=lower_tail)) # shfited exponential
#
#     ly <- length(q)
#     sigma <- rep(sigma, length = ly)
#     mu <- rep(mu, length = ly)
#     tau <- rep(tau, length = ly)
#     index <- seq(along = q)
#     z <- q - mu - ((sigma^2)/tau)
#     cdf <- ifelse(is.finite(q),
#       ifelse(tau > 0.05 * sigma,
#          pnorm((q - mu)/sigma) - exp(log(pnorm(z/sigma)) + ((mu + (sigma^2/tau))^2 -
#            (mu^2) -  2 * q * ((sigma^2)/tau))/(2 * sigma^2)),
#          pnorm(q, mean = mu, sd = sigma)),
#         ifelse(q<0,0,1)
#     )
#     if (lower_tail == TRUE)
#       cdf <- cdf
#     else cdf <- 1 - cdf
#     if (log_p == FALSE)
#       cdf <- cdf
#     else cdf <- log(cdf)
#     cdf
#   }
#
# dEXG <- function (x, mu = 5, sigma = 1, tau = 1, log = FALSE)
#   # ex-Gaussian density
#   # gamlss.dist function, but returns NA for bad sigma or tau, and
#   # robust against small sigma cases.
# {
#     if (sigma <= 0) return(rep(NA,length(x)))
#     if (tau <= 0) return(rep(NA,length(x)))
#
#     # if (sigma < 0.05*tau)
#     if (sigma < 1e-4)
#       return(dexp(x-mu,1/tau,log=log)) # shfited exponential
#
#     ly <- length(x)
#     sigma <- rep(sigma, length = ly)
#     mu <- rep(mu, length = ly)
#     tau <- rep(tau, length = ly)
#     z <- x - mu - ((sigma^2)/tau)
#     logfy <- ifelse(tau > 0.05 * sigma,
#       -log(tau) - (z + (sigma^2/(2 *  tau)))/tau + log(pnorm(z/sigma)),
#       dnorm(x, mean = mu, sd = sigma, log = TRUE))
#     if (log == FALSE)
#       fy <- exp(logfy)
#     else fy <- logfy
#     fy
# }
#

dexGaussian <- function(rt,pars)
  # exGaussian pdf (returns normal or exponential for small tau/sigma)
{
  isexp <- pars[,"sigma"] < 1e-4 # shifted exponential
  rt[isexp] <- dexp(rt[isexp]-pars[isexp,"mu"],1/pars[isexp,"tau"])
  isnorm <- !isexp & pars[,"tau"] < 0.05 * pars[,"sigma"] # normal
  rt[isnorm] <- dnorm(rt[isnorm], mean = pars[isnorm,"mu"],
                                  sd = pars[isnorm,"sigma"])
  isexg <- !(isexp | isnorm)
  if (any(isexg)) {
    s2 <- pars[isexg,"sigma"]^2
    z <- rt[isexg] - pars[isexg,"mu"] - (s2/pars[isexg,"tau"])
    rt[isexg] <- exp(
      log(pnorm(z/pars[isexg,"sigma"])) -
        log(pars[isexg,"tau"]) -
        (z + (s2/(2 *  pars[isexg,"tau"])))/pars[isexg,"tau"]
    )
  }
  rt
}

pexGaussian <- function(rt,pars)
  # exGaussian cdf (returns normal or exponential for small tau/sigma)
{
  isexp <- pars[,"sigma"] < 1e-4 # shifted exponential
  rt[isexp] <- pexp(rt[isexp]-pars[isexp,"mu"],1/pars[isexp,"tau"])
  isnorm <- !isexp & pars[,"tau"] < 0.05 * pars[,"sigma"] # normal
  rt[isnorm] <- pnorm(rt[isnorm], mean = pars[isnorm,"mu"],
                                  sd = pars[isnorm,"sigma"])
  isexg <- !(isexp | isnorm)
  if (any(isexg)) {
    s2 <- pars[isexg,"sigma"]^2
    z <- rt[isexg] - pars[isexg,"mu"] - (s2/pars[isexg,"tau"])
    rt[isexg] <-
      pnorm((rt[isexg] - pars[isexg,"mu"])/pars[isexg,"sigma"]) -
      exp(log(pnorm(z/pars[isexg,"sigma"])) +
        ((pars[isexg,"mu"] + (s2/pars[isexg,"tau"]))^2 - (pars[isexg,"mu"]^2) -
          2 * rt[isexg] * (s2/pars[isexg,"tau"]))/(2 * s2))
  }
  rt
}

# Go cdf/pdf (strips out NAs and calls d/pexGaussian)
# Is stripping out rt NAs really necessary?

dexGaussianG <- function(rt,pars)
{
  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- dexGaussian(rt[ok],pars[ok,,drop=FALSE])
  out
}

pexGaussianG <- function(rt,pars)
{
  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- pexGaussian(rt[ok],pars[ok,,drop=FALSE])
  out
}


#### Stop Single ExGaussian ----

# Stop cdf/pdf, subtracts SSD and uses muS/sigmaS/tauS by renaming

dexGaussianS <- function(rt,pars)
{
  rt <- rt - pars[,"SSD"]
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="muS"] <- "mu"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="sigmaS"] <- "sigma"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="tauS"] <- "tau"
  dexGaussian(rt,pars)
}


pexGaussianS <- function(rt,pars)
{
  rt <- rt - pars[,"SSD"]
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="muS"] <- "mu"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="sigmaS"] <- "sigma"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="tauS"] <- "tau"
  pexGaussian(rt,pars)
}


#### ExG Race function ----

# Following functions moved to C++ model_SS_EXG.cpp

# dEXGrace <- function(dt,mu,sigma,tau)
#   # Generates defective PDF for win by first runner, dt (decison time) is
#   # a matrix with length(mu) rows, one row for each runner, and one column
#   # for each decision time for which a defective density value will be
#   # returned.
# {
#   dt[1,] <- dEXG(dt[1,],mu[1],sigma[1],tau[1])
#   if (length(mu)>1) for (i in 2:length(mu))
#     dt[1,] <- dt[1,]*pEXG(dt[i,],mu[i],sigma[i],tau[i],lower_tail=FALSE)
#   dt[1,]
# }

#
#
# stopfn_exg <- function(t,mu,sigma,tau,SSD)
#   # Used by my.integrate, t = vector of times, SSD is a scalar stop-signal delay.
# {
#   dt <- matrix(rep(t+SSD,each=length(mu)),nrow=length(mu))
#   dt[1,] <- dt[1,]-SSD
#   dEXGrace(dt,mu,sigma,tau)
# }

#### ExGaussian random ----

rexG <- function(n,mu,sigma,tau)
  rnorm(n,mean=mu,sd=sigma) + rexp(n,rate=1/tau)

rexGaussian <- function(lR,pars,p_types=c("mu","sigma","tau"),
                        ok=rep(TRUE,dim(pars)[1]))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in
  # contiguous rows.
  #
  # test
  # pars=cbind(mu=c(.5,.6),sigma=c(.1,.1),tau=c(.2,.2)); lR=factor(c(1))
{
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  dt <- matrix(rexG(dim(pars)[1],pars[,"mu"],pars[,"sigma"],pars[,"tau"]),
               nrow=length(levels(lR)))
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  rt <- dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  cbind.data.frame(R=R,rt=rt)
}


#### EXG Stop signal random -----

# FIX ME ok NOT IMPLEMENTED

rSSexGaussian <- function(data,pars,ok=rep(TRUE,dim(pars)[1]))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars must contain an SSD column and an lI column indicating if an
  # accumulator is triggered by the stop signal (ST = stop-triggered).

  # NB1: Go failures will only apply to accumulators where lI = TRUE
  #      and can still have a stop-triggered response on a go-failure trial.
{
  lR <- data$lR
  nacc <- length(levels(lR))   # Does not include stop runner
  ntrials <- dim(pars)[1]/nacc # Number of trials to simulate
  is1 <- lR==levels(lR)[1]     # First go accumulator
  acc <- 1:nacc                # choice (go and ST) accumulator index
  allSSD <- data$SSD[is1]
  # stop-triggered racers
  isST <- pars[,"lI"]==1              # Boolean for all pars
  accST <- acc[pars[1:nacc,"lI"]==1]  # Index of ST accumulator, from 1st trial

  # Setup race among nacc+1 (includes stop accumulator)
  # Default Inf finishing time so if not changed always looses
  dt <- matrix(Inf,nrow=nacc+1,ncol=ntrials)

  # Go failure trials (rep for each accumulator)
  isgf <- rep(pars[is1,"gf"] > runif(ntrials),each=nacc)
  # Go accumulators that don't fail (removing ST)
  isGO <- !isgf & !isST
  ngo <- sum(isGO)

  # Fill in go accumulators
  if (any(isGO)) dt[-1,][isGO] <- rexG(ngo,pars[isGO,"mu"],pars[isGO,"sigma"],
                                       pars[isGO,"tau"])

  # pick out stop trials and races with SSD that is not Inf (i.e., finite or
  # NA, the latter so a staircase can be filled in)
  isStrial <- !is.infinite(pars[is1,"SSD"])
  isSrace <- rep(isStrial,each=nacc)

  # pick out stop trials that are triggered
  isT <- pars[is1,"tf"][isStrial] < runif(sum(isStrial))

  # Logical to pick stop-triggered accumulators that are triggered
  isSTT <- logical(ntrials*nacc) # Default false

  # Pick out stop-triggered accumulators that are triggered
  # NB: can occur on gf trials
  isSTT[isSrace][rep(isT,each=nacc) & isST[isSrace]] <- TRUE
  nst <- sum(isSTT) # Number of triggered ST accumulators

  # Fill in stop-triggered accumulators
  if (any(isSTT)) dt[-1,][isSTT] <- rexG(nst,pars[isSTT,"mu"],
                                      pars[isSTT,"sigma"],pars[isSTT,"tau"])

  # pick out triggered stop racers
  isTS <- logical(ntrials)
  isTS[isStrial][isT] <- TRUE
  ns <- sum(isTS)

  # Fill in stop accumulators
  if (any(isTS)) dt[1,isTS] <- rexG(ns,pars[is1,"muS"][isTS],
                                 pars[is1,"sigmaS"][isTS],pars[is1,"tauS"][isTS])

  # staircase algorithm
  pstair <- is.na(pars[,"SSD"])
  stair <- pstair[is1]
  if (any(stair)) {
    if (is.null(attr(data,"staircase")))
      stop("When SSD has NAs a staircase list must be supplied!")
    staircase <- attr(data,"staircase")
    if (length(accST)>0)
      staircase$accST <- 1+accST
    if (is.null(attr(staircase,"staircase_function"))) {
      if (length(accST)>0)
        stop("Do not use default starircase function with stop-triggered accumulators")
      attr(staircase,"staircase_function") <- staircase_function
    }

    allR <- allrt <- numeric(ncol(dt))  # to store unified results
    dts <- dt[,stair,drop=F]

    # Non-staircase trials
    dt <- dt[,!stair,drop=F]
    spars <- pars[pstair,,drop=F]
    pars <- pars[!pstair,,drop=F]
    isTS <- isTS[!stair]
    isSTT <- isSTT[!pstair]
    isST <- isST[!pstair]
    ntrials <- sum(!stair)
    is1 <- is1[!pstair]
    stair_res <- attr(staircase,"staircase_function")(dts,staircase)
    allR[stair] <- stair_res$sR
    allrt[stair] <- stair_res$srt
  }


  # All SSD already filled in so return R and rt

  # Add SSD to triggered stop accumulators
  if (any(isTS)) dt[1,isTS] <- dt[1,isTS] + pars[is1,"SSD"][isTS]

  # Add SSD to stop-triggered accumulators that are triggered
  if (any(isSTT)) dt[-1,][isSTT] <- dt[-1,][isSTT] + pars[isSTT,"SSD"]

  # R <- factor(rep(NA,ntrials),levels=levels(lR))
  R <- rt <- rep(NA,ntrials)

  # All accumulators Inf (usually when both gf and tf)
  allinf <- apply(dt,2,function(x)all(is.infinite(x)))

  # get winner of stop and go where there is a race
  r <- c(0, acc)[apply(dt[,!allinf,drop=FALSE],2,which.min)]

  # stop wins
  stopwins <- r==0

  # First fill in cases where stop looses, can include wins by ST
  if (any(!stopwins)) {
    rgo <- r[!stopwins]
    R[!allinf][!stopwins] <- rgo
    pick <- cbind(rgo,c(1:sum(!stopwins))) # Matrix to pick winner
    rt[!allinf][!stopwins] <-
      dt[-1,!allinf,drop=FALSE][,!stopwins,drop=FALSE][pick]
  }

  # then if stop wins and kills go, find ST winner
  if (any(isST) & any(stopwins)) {
    # stop triggered accumulators that are racing
    rst <- dt[-1,!allinf,drop=FALSE][accST,stopwins,drop=FALSE]
    # stop-triggered winners
    rtw <- apply(rst,2,which.min)
    # index for stop-triggered
    R[!allinf][stopwins] <- accST[rtw]
    pick <- cbind(rtw,1:ncol(rst))
    rt[!allinf][stopwins] <- rst[pick]
  }

  rt[is.na(R)] <- NA

  if (any(stair)) {
    allrt[!stair] <- rt
    allR[!stair] <- R
    allSSD[stair] <- stair_res$SSD
    out <- cbind.data.frame(R=factor(allR,levels=1:nacc,labels=levels(lR)),
                            rt=allrt, SSD = allSSD)
    return(out)
  }
  cbind.data.frame(R=factor(R,levels=1:nacc,labels=levels(lR)),rt=rt)
}


#### ExG stop probability (no stop triggered) ----

pstopEXG <- function(parstop,n_acc,upper=Inf,
                     gpars=c("mu","sigma","tau"),spars=c("muS","sigmaS","tauS"))
{
  sindex <- seq(1,nrow(parstop),by=n_acc)  # Stop accumulator index
  ps <- parstop[sindex,spars,drop=FALSE]   # Stop accumulator parameters
  SSDs <- parstop[sindex,"SSD",drop=FALSE] # SSDs
  ntrials <- length(SSDs)
  if (length(upper)==1) upper <- rep(upper,length.out=ntrials)
  pgo <- array(parstop[,gpars],dim=c(n_acc,ntrials,length(gpars)),
               dimnames=list(NULL,NULL,gpars))
  cells <- apply(
    cbind(SSDs,ps,upper,matrix(as.vector(aperm(pgo,c(2,1,3))),nrow=ntrials))
  ,1,paste,collapse="")
  # cells <- character(ntrials)
  # for (i in 1:ntrials)
  #   cells[i] <- paste(SSDs[i],ps[i,],pgo[,i,],upper[i],collapse="")
  uniq <- !duplicated(cells)
  ups <- sapply(1:sum(uniq),function(i){
    my.integrate(f=stopfn_exg,lower=-Inf,SSD=SSDs[i],upper=upper[i],
                           mu=c(ps[i,"muS"],pgo[,i,"mu"]),
                           sigma=c(ps[i,"sigmaS"],pgo[,i,"sigma"]),
                           tau=c(ps[i,"tauS"],pgo[,i,"tau"]))
  })
  ups[as.numeric(factor(cells,levels=cells[uniq]))]
}

# #### Stop probability stop triggered ----
#
#
# stopfn_exgST <- function(t,mu,sigma,tau,SSD,st=1)
#   # Used by my.integrate, t = vector of times, SSD is a scalar stop-signal delay.
#   # st is a vector of indices for the stop and stop-triggered accumulators
# {
#   dt <- matrix(rep(t+SSD,each=length(mu)),nrow=length(mu))
#   dt[st,] <- dt[st,]-SSD
#   dEXGrace(dt,mu,sigma,tau)
# }
#
#
# pstopEXGST <- function(parstop,n_acc,upper=Inf,st=1,
#                      gpars=c("mu","sigma","tau"),spars=c("muS","sigmaS","tauS"))
# {
#   sindex <- seq(1,nrow(parstop),by=n_acc)
#   ps <- parstop[sindex,spars,drop=FALSE]
#   SSDs <- parstop[sindex,"SSD",drop=FALSE]
#   ntrials <- length(SSDs)
#   if (length(upper)==1) upper <- rep(upper,length.out=ntrials)
#   pgo <- array(parstop[,gpars],dim=c(n_acc,ntrials,length(gpars)),
#                dimnames=list(NULL,NULL,gpars))
#   cells <- apply(cbind(SSDs,ps,upper,
#     matrix(as.vector(aperm(pgo,c(2,1,3))),nrow=ntrials)),1,paste,collapse="")
#   # cells <- character(ntrials)
#   # for (i in 1:ntrials)
#   #   cells[i] <- paste(SSDs[i],ps[i,],pgo[,i,],upper[i],collapse="")
#   uniq <- !duplicated(cells)
#   ups <- sapply(1:sum(uniq),function(i){
#     my.integrate(f=stopfn_exgST,lower=-Inf,SSD=SSDs[i],upper=upper[i],
#                            mu=c(ps[i,"muS"],pgo[,i,"mu"]),
#                            sigma=c(ps[i,"sigmaS"],pgo[,i,"sigma"]),
#                            tau=c(ps[i,"tauS"],pgo[,i,"tau"]),st=st)
#   })
#   ups[as.numeric(factor(cells,levels=cells[uniq]))]
# }



#### SSexG Model list ----

#' Stop-signal ex-Gaussian race
#'
#' Model file to estimate the ex-Gaussian race model for Stop-Signal data
#'
#' Model files are almost exclusively used in `design()`.
#'
#' @details
#'
#' Default values are used for all parameters that are not explicitly listed in the `formula`
#' argument of `design()`.They can also be accessed with `SSexG()$p_types`.
#'
#' @details
#' | **Parameter** | **Transform** | **Natural scale** | **Default**   | **Mapping**                    | **Interpretation**                                            |
#' |-----------|-----------|---------------|-----------|----------------------------|-----------------------------------------------------------|
#' | *mu*       | log         | \[0, Inf\]     | log(.4)         |                            | mu parameter of ex-Gaussian go finishing time distribution              |
#' | *sigma*       | log       | \[0, Inf\]        | log(.05)    |                            | sigma parameter of ex-Gaussian go finishing time distribution                                        |
#' | *tau*      | log       | \[0, Inf\]        | log(.1)    |                            | tau parameter of ex-Gaussian go finishing time distribution                                          |
#' | *muS*       | log       | \[0, Inf\]        | log(.3)    |                            | mu parameter of ex-Gaussian stopping finishing time distribution           |
#' | *sigmaS*       | log    | \[0, Inf\]        | log(.025)|                   | sigma parameter of ex-Gaussian stopping finishing time distribution                              |
#' | *tauS*      | log    | \[0, Inf\]        | log(.05)  |  | tau parameter of ex-Gaussian stopping finishing time distribution       |
#' | *tf*      | probit       | \[0, 1\]        | qnorm(0)    |                            | Trigger failure probability           |
#' | *gf*     | probit       | \[0, 1\]        | qnorm(0)    |                            | Go failure probability    |
#'
#'
#'
#' @return A model list with all the necessary functions to sample
#' @export
SSexG <- function() {
  list(
    type="RACE",
    p_types=c(mu=log(.4),sigma=log(.05),tau=log(.1),
              muS=log(.3),sigmaS=log(.025),tauS=log(.05),tf=qnorm(0),gf=qnorm(0)),
    transform=list(func=c( mu = "exp",  sigma = "exp",  tau = "exp",
                          muS = "exp", sigmaS = "exp", tauS = "exp",
                          tf="pnorm",gf="pnorm")),
    bound=list(minmax=cbind( mu=c(0,Inf),  sigma=c(1e-4,Inf),  tau=c(1e-4,Inf),
                            muS=c(0,Inf), sigmaS=c(1e-4,Inf), tauS=c(1e-4,Inf),
                            tf=c(.001,.999),gf=c(.001,.999)),
                            exception=c(tf=0,gf=0)),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      # if (any(names(dadm)=="SSD")) pars <- cbind(pars,SSD=dadm$SSD) else
      #   pars <- cbind(pars,SSD=rep(Inf,dim(pars)[1]))
      pars <- cbind(pars,SSD=dadm$SSD)
      pars <- cbind(pars,lI=as.numeric(dadm$lI))  # Only necessary for data generation.
      pars
    },
    # Density function (PDF) for single go racer
    dfunG=function(rt,pars) dexGaussianG(rt,pars),
    # Probability function (CDF) for single go racer
    pfunG=function(rt,pars) pexGaussianG(rt,pars),
    # Density function (PDF) for single stop racer
    dfunS=function(rt,pars)
      dexGaussianS(rt,pars[,c("muS","sigmaS","tauS","SSD"),drop=FALSE]),
    # Probability function (CDF) for single stop racer
    pfunS=function(rt,pars)
      pexGaussianS(rt,pars[,c("muS","sigmaS","tauS","SSD"),drop=FALSE]),
    # Stop probability integral
    sfun=function(pars,n_acc,upper=Inf) pstopEXG(pars,n_acc,upper=upper),
    # Random function for SS race
    rfun=function(data=NULL,pars) {
      rSSexGaussian(data,pars,ok=attr(pars, "ok"))
    },
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10))
      log_likelihood_race_ss(pars, dadm, model, min_ll = min_ll),
    # SS likelihood with deadline censoring/truncation support (rt==Inf uses UC)
    log_likelihood_cens_trunc=function(pars,dadm,model,min_ll=log(1e-10), normalise_trunc=FALSE)
      log_likelihood_race_ss_cens_trunc(pars, dadm, model, min_ll = min_ll, normalise_trunc = normalise_trunc)
  )
}


#####################  RDEX ----

#### RDEX SS random ----

rSShybrid <- function(data,pars,ok=rep(TRUE,dim(pars)[1]))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars must contain an SSD column and an lI column indicating if an
  # accumulator is triggered by the stop signal (ST = stop-triggered).

  # NB1: Go failures will only apply to accumulators where lI = TRUE
  #      and can still have a stop-triggered response on a go-failure trial.
{
  lR <- data$lR
  pars[,c("A","B","v")] <- pars[,c("A","B","v")]/pars[ok,"s"]

  nacc <- length(levels(lR))   # Does not include stop runner
  ntrials <- dim(pars)[1]/nacc # Number of trials to simulate
  is1 <- lR==levels(lR)[1]     # First go accumulator
  acc <- 1:nacc                # choice (go and ST) accumulator index

  # stop-triggered racers
  isST <- pars[,"lI"]==1              # Boolean for all pars
  accST <- acc[pars[1:nacc,"lI"]==1]  # Index of ST accumulator, from 1st trial

  # Setup race among nacc+1 (includes stop accumulator)
  # Default Inf finishing time so if not changed always looses
  dt <- matrix(Inf,nrow=nacc+1,ncol=ntrials)

  # Go failure trials (rep for each accumulator)
  isgf <- rep(pars[is1,"gf"] > runif(ntrials),each=nacc)
  # Go accumulators that don't fail (removing ST)
  isGO <- !isgf & !isST
  ngo <- sum(isGO)

  # Fill in go accumulators
  if (any(isGO)) dt[-1,][isGO] <-
    pars[isGO,"t0"] + rWald(ngo,pars[isGO,"B"],pars[isGO,"v"],pars[isGO,"A"])

  # pick out stop trials and races with SSD that is not Inf (i.e., finite or
  # NA, the latter so a staircase can be filled in)
  isStrial <- !is.infinite(pars[is1,"SSD"])
  isSrace <- rep(isStrial,each=nacc)

  # pick out stop trials that are triggered
  isT <- pars[is1,"tf"][isStrial] < runif(sum(isStrial))

  # Logical to pick stop-triggered accumulators that are triggered
  isSTT <- logical(ntrials*nacc) # Default false

  # Pick out stop-triggered accumulators that are triggered
  # NB: can occur on gf trials
  isSTT[isSrace][rep(isT,each=nacc) & isST[isSrace]] <- TRUE
  nst <- sum(isSTT)

  # Fill in stop-triggered accumulators
  if (any(isSTT)) dt[-1,][isSTT] <-
    pars[isSTT,"t0"] + rWald(ngo,pars[isSTT,"B"],pars[isSTT,"v"],pars[isSTT,"A"])

  # pick out triggered stop racers
  isTS <- logical(ntrials)
  isTS[isStrial][isT] <- TRUE
  ns <- sum(isTS)

  # Fill in stop accumulators
  if (any(isTS)) dt[1,isTS] <- rexG(ns,pars[is1,"muS"][isTS],
    pars[is1,"sigmaS"][isTS],pars[is1,"tauS"][isTS])

  # staircase algorithm
  pstair <- is.na(pars[,"SSD"])
  stair <- pstair[is1]
  if (any(stair)) {
    if (is.null(attr(pars,"staircase")))
      stop("When SSD has NAs a staircase list must be supplied!")

    staircase <- attr(pars,"staircase")
    data <- data[is1,]
    if (length(accST)>0)
      staircase$accST <- 1+accST
    if (is.null(attr(staircase,"staircase_function")))
      attr(staircase,"staircase_function") <- staircase_function

    allR <- allrt <- numeric(ncol(dt))  # to store unified results
    dts <- dt[,stair,drop=F]

    # Non-staircase trials
    dt <- dt[,!stair,drop=F]
    spars <- pars[pstair,,drop=F]
    pars <- pars[!pstair,,drop=F]
    isTS <- isTS[!stair]
    isSTT <- isSTT[!pstair]
    isST <- isST[!pstair]
    ntrials <- sum(!stair)
    is1 <- is1[!pstair]
    data <- data[stair,]

    stair_res <- attr(staircase,"staircase_function")(dts,staircase)
    allR[stair] <- stair_res$sR
    allrt[stair] <- stair_res$srt
  }


  # All SSD already filled in so return R and rt

  # Add SSD to triggered stop accumulators
  if (any(isTS)) dt[1,isTS] <- dt[1,isTS] + pars[is1,"SSD"][isTS]

  # Add SSD to stop-triggered accumulators that are triggered
  if (any(isSTT)) dt[-1,][isSTT] <- dt[-1,][isSTT] + pars[isSTT,"SSD"]

  # R <- factor(rep(NA,ntrials),levels=levels(lR))
  R <- rt <- rep(NA,ntrials)

  # All accumulators Inf (usually when both gf and tf)
  allinf <- apply(dt,2,function(x)all(is.infinite(x)))

  # get winner of stop and go where there is a race
  r <- c(0, acc)[apply(dt[,!allinf,drop=FALSE],2,which.min)]

  # stop wins
  stopwins <- r==0

  # First fill in cases where stop looses, can include wins by ST
  if (any(!stopwins)) {
    rgo <- r[!stopwins]
    R[!allinf][!stopwins] <- rgo
    pick <- cbind(rgo,c(1:sum(!stopwins))) # Matrix to pick winner
    rt[!allinf][!stopwins] <-
      dt[-1,!allinf,drop=FALSE][,!stopwins,drop=FALSE][pick]
  }

  # then if stop wins and kills go, find ST winner
  if (any(isST) & any(stopwins)) {
    # stop triggered accumulators that are racing
    rst <- dt[-1,!allinf,drop=FALSE][accST,stopwins,drop=FALSE]
    # stop-triggered winners
    rtw <- apply(rst,2,which.min)
    # index for stop-triggered
    R[!allinf][stopwins] <- accST[rtw]
    pick <- cbind(rtw,1:ncol(rst))
    rt[!allinf][stopwins] <- rst[pick]
  }

  rt[is.na(R)] <- NA

  if (any(stair)) {
    allrt[!stair] <- rt
    allR[!stair] <- R
    allSSD <- NA
    allSSD[stair] <- stair_res$SSD
    out <- cbind.data.frame(R=factor(allR,levels=1:nacc,labels=levels(lR)),
                            rt=allrt, SSD = allSSD)
    return(out)
  }
  cbind.data.frame(R=factor(R,levels=1:nacc,labels=levels(lR)),rt=rt)
}

#### RDEX stop probability ----

# # NB: these functions are in Rcpp
# dWald_RDEX
# pWald_RDEX
# EMC2:::stopfn_rdex


pstopHybrid <- function(parstop,n_acc,upper=Inf,
  gpars=c("v","B","t0","A"),spars=c("muS","sigmaS","tauS"))
{
  sindex <- seq(1,nrow(parstop),by=n_acc)  # Stop accumulator index
  ps <- parstop[sindex,spars,drop=FALSE]   # Stop accumulator parameters
  SSDs <- parstop[sindex,"SSD",drop=FALSE] # SSDs
  ntrials <- length(SSDs)
  if (length(upper)==1) upper <- rep(upper,length.out=ntrials)
  pgo <- array(parstop[,gpars],dim=c(n_acc,ntrials,length(gpars)),
               dimnames=list(NULL,NULL,gpars))
  cells <- apply(
    cbind(SSDs,ps,upper,matrix(as.vector(aperm(pgo,c(2,1,3))),nrow=ntrials))
  ,1,paste,collapse="")
  uniq <- !duplicated(cells)
  ups <- sapply(1:sum(uniq),function(i){
    my.integrate(f=stopfn_rdex,lower=0,upper=upper[i],
      mu=ps[i,"muS"],sigma=ps[i,"sigmaS"],tau=ps[i,"tauS"],
      v=pgo[,i,"v"],B=pgo[,i,"B"],A=pgo[,i,"A"],t0=pgo[,i,"t0"],
      SSD=SSDs[i],n_acc=n_acc)
  })
  ups[as.numeric(factor(cells,levels=cells[uniq]))]
}


#### RDEX model list ----

#' Stop-signal Hybrid (RDM go, ExGaussian stop) race
#'
#' Model file to estimate the Hybrid race model for Stop-Signal data
#'
#' Model files are almost exclusively used in `design()`.
#'
#' @details
#'
#' Default values are used for all parameters that are not explicitly listed in the `formula`
#' argument of `design()`.They can also be accessed with `SShybrid()$p_types`.
#'
#' @return A model list with all the necessary functions to sample
#' @export
SShybrid <- function() {
  list(
    type="RACE",
    p_types=c("v" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),"s" = log(1),
              muS=log(.3),sigmaS=log(.025),tauS=log(.05),tf=qnorm(0),gf=qnorm(0)),
    transform=list(func=c(v = "exp", B = "exp", A = "exp",t0 = "exp", s = "exp",
                          muS = "exp", sigmaS = "exp", tauS = "exp",
                          tf="pnorm",gf="pnorm")),
    bound=list(minmax=cbind(v=c(1e-3,Inf), B=c(0,Inf), A=c(1e-4,Inf),t0=c(0.05,Inf),
                            s=c(0,Inf),muS=c(0,Inf), sigmaS=c(1e-4,Inf), tauS=c(1e-4,Inf),
                            tf=c(.001,.999),gf=c(.001,.999)),
               exception=c(A=0, v=0,tf=0,gf=0)),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
      pars <- cbind(pars,SSD=dadm$SSD)
      pars <- cbind(pars,lI=as.numeric(dadm$lI))  # Only necessary for data generation.
      pars
    },
    # Density function (PDF) for single go racer
    dfunG=function(rt,pars) dRDM(rt,pars),
    # Probability function (CDF) for single go racer
    pfunG=function(rt,pars) pRDM(rt,pars),
    # Density function (PDF) for single stop racer
    dfunS=function(rt,pars) dexGaussianS(rt,
      pars[,c("muS","sigmaS","tauS","SSD"),drop=FALSE]),
    # Probability function (CDF) for single stop racer
    pfunS=function(rt,pars) pexGaussianS(rt,
      pars[,c("muS","sigmaS","tauS","SSD"),drop=FALSE]),
    # Stop probability integral
    sfun=function(pars,n_acc,upper=Inf) pstopHybrid(pars,n_acc,upper=upper),
    # Random function for SS race
    rfun=function(data=NULL,pars) {
      rSShybrid(data,pars,ok=attr(pars, "ok"))
    },
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10))
      log_likelihood_race_ss(pars, dadm, model, min_ll = min_ll),
    # SS likelihood with deadline censoring/truncation support (rt==Inf uses UC)
    log_likelihood_cens_trunc=function(pars,dadm,model,min_ll=log(1e-10), normalise_trunc=FALSE)
      log_likelihood_race_ss_cens_trunc(pars, dadm, model, min_ll = min_ll, normalise_trunc = normalise_trunc)
  )
}


####  Stop-signal ----

my.integrate <- function(...,upper=Inf,big=10)
  # Avoids bug in integrate upper=Inf that uses only 1  subdivision
  # Use of  big=10 is arbitrary ...
{
  out <- try(integrate(...,upper=upper),silent=TRUE)
  if (inherits(out, "try-error")) 0 else
  {
    if (upper==Inf & out$subdivisions==1)
    {
      out <- try(integrate(...,upper=big),silent=TRUE)
      if (inherits(out, "try-error")) 0 else
      {
        if (out$subdivisions==1) 0 else out$value
      }
    } else out$value
  }
}

log_likelihood_race_ss <- function(pars,dadm,model,min_ll=log(1e-10))
{
  # All bad?
  if (is.null(attr(pars,"ok")))
    ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")
    if (!any(ok)) return(min_ll*length(attr(dadm, "expand")))

    # Nomenclature for indices:
    # "is" = logical, "isp" pars/dadm index,
    # "t" trials index, "ist" logical on trial index
    # "n_" number of integer

    # Counts
    n_acc <- length(levels(dadm$lR))                   # total number of accumulators
    n_trials <- nrow(dadm)/n_acc                       # number of trials
    n_accG <- sum(as.numeric(dadm[1:n_acc,"lI"])==2)   # go accumulators
    n_accST <- sum(as.numeric(dadm[1:n_acc,"lI"])==1)

    # Likelihood for all trials and for ok trials
    allLL <- rep(min_ll,n_trials)
    # used to put results into allLL[allok]
    allok <- ok[dadm$lR==levels(dadm$lR)[1]]

    # Remove trials not ok
    pars <- pars[ok,,drop=FALSE]
    dadm <- dadm[ok,,drop=FALSE]

    # Only keep OK trials for stop trial computations
    n_trials <- nrow(dadm)/n_acc # number of trials
    trials <- 1:n_trials         # trial number

    # Booleans for ok pars/dadm
    isp1 <- dadm$lR==levels(dadm$lR)[1]      # 1st accumulator rows
    ispGOacc <- dadm$lI==levels(dadm$lI)[2]  # Go accumulator rows
    ispStop <- is.finite(dadm$SSD)           # Stop-trial rows

    # Failure parameters for each trial
    gf <- pars[isp1,"gf"]
    tf <- pars[isp1,"tf"]

    # Go response names
    GoR <- as.character(dadm[1:n_acc,"lR"][dadm[1:n_acc,"lI"]==2])

    # No response
    ispNR <- is.na(dadm$R)
    if ( any(ispNR) ) {  # Note by definition no ST present

      # Go failures
      ispgoNR <- ispNR & !ispStop            # No response and no stop signal
      tgoNR <- c(1:sum(isp1))[ispgoNR[isp1]] # trial number
      if (any(ispgoNR))
        allLL[allok][tgoNR] <- log(gf[tgoNR])

      # Stop trial with no response and no ST accumulator
      ispstopNR <- ispNR & ispStop & (n_accST == 0)
      if ( any(ispstopNR) ) { # Stop trial probability
        # Non-response and stop trial & go/stop accumulator parameters
        pStop <- pmin(1,pmax(0,  # protection to keep in 0-1
                             model$sfun(
                               pars[ispNR & ispStop & ispGOacc,,drop=FALSE],n_acc=n_accG)
        ))
        # Fill in stop-trial non-response probabilities, either 1) go failure
        # 2) not go failure and not trigger failure and stop wins
        tstopNR <- trials[ispstopNR[isp1]]  # trial number
        allLL[allok][tstopNR] <- log(gf[tstopNR] + (1-gf[tstopNR])*(1-tf[tstopNR])*pStop)
      }
    }

    # Response made
    if (any(!ispNR)) {
      # Only keep response trials for further computation
      allr <- !ispNR[isp1] # used to put back into allLL[allok]
      pars <- pars[!ispNR,,drop=FALSE]
      dadm <- dadm[!ispNR,,drop=FALSE]
      n_trials <- nrow(dadm)/n_acc # number of trials
      trials <- 1:n_trials
      ptrials <- rep(trials,each=n_acc) # trial number for pars/dadm

      isp1 <- dadm$lR==levels(dadm$lR)[1]      # 1st accumulator rows
      ispGOacc <- dadm$lI==levels(dadm$lI)[2]  # Go accumulator rows
      ispStop <- is.finite(dadm$SSD)           # stop-trial rows with a response
      gf <- pars[isp1,"gf"]
      tf <- pars[isp1,"tf"]

      # likelihoods for trials with a response (eventually, intermediately probabilities)
      like <- numeric(n_trials)
      lds <- numeric(nrow(dadm)) # log density and survivor, used for both go and stop trials

      # Go trials with response
      if (any(!ispStop)) {
        ispGOwin <-  !ispStop & dadm$winner # Winner go accumulator rows
        tGO <- ptrials[ispGOwin]  # Go trials
        # Winner density
        like[tGO] <- log(model$dfunG(
          rt=dadm$rt[ispGOwin],pars=pars[ispGOwin,,drop=FALSE]))
        if (n_accG >1) {  # Looser survivor go accumulator(s)
          ispGOloss <- !ispStop & !dadm$winner & ispGOacc # Looser go accumulator rows
          n_expected <- (n_accG - 1) * length(tGO)
          if (sum(ispGOloss) != n_expected) {
            cdf_loss <- model$pfunG(rt = dadm$rt[ispGOloss], pars = pars[ispGOloss, , drop = FALSE])
            cdf_loss[is.na(cdf_loss)] <- 1
            cdf_loss <- pmin(1, pmax(0, cdf_loss))
            ll_loss <- log(pmax(1 - cdf_loss, 0))
            loss_sum <- rowsum(ll_loss, group = ptrials[ispGOloss], reorder = FALSE)
            like[as.integer(rownames(loss_sum))] <- like[as.integer(rownames(loss_sum))] + loss_sum[, 1]
          } else {
            like[tGO] <- like[tGO] + apply(matrix(log(1 - model$pfunG(
              rt = dadm$rt[ispGOloss], pars = pars[ispGOloss, , drop = FALSE]
            )), nrow = n_accG - 1), 2, sum)
          }
        }
        # Transform back to densities to include go failure
        like[tGO] <- (1-gf[tGO])*exp(like[tGO])
      }

      # Stop trials with a response, occurs if
      # 1) All triggered
      #   a) Stop does not beat go before rt, produces go and ST
      #   b) Stop beats go before rt, produces ST only
      # 2) go failure, stop triggered, produces ST only
      # 3) go triggered, stop failure, produces go only
      # NB: cant have both go and stop failure as then no response

      if (any(ispStop)) {                       # Stop trials
        ispSGO <- ispStop & (dadm$R %in% GoR)   # Go responses
        ispSST <- ispStop & !(dadm$R %in% GoR)  # ST responses
        # Go beats stop and ST (if any)
        if (any(ispSGO)) {
          ispGOwin <-  ispSGO & dadm$winner # Winner go accumulator rows
          tGO <- ptrials[ispGOwin]          # Go trials
          like[tGO] <- log(model$dfunG(
            rt=dadm$rt[ispGOwin],pars=pars[ispGOwin,,drop=FALSE]))
          if (n_accG > 1) {  # Looser survivor gp accumulators
            ispGOloss <- ispSGO & !dadm$winner & ispGOacc
            like[tGO] <- like[tGO] + apply(matrix(log(1-model$pfunG(
              rt=dadm$rt[ispGOloss],pars=pars[ispGOloss,,drop=FALSE])),
              nrow=n_accG-1),2,sum)
          }
          # trigger stop, add in stop survivor
          ts <- like[tGO] + log(1-model$pfunS(
            rt=dadm$rt[ispGOwin],pars=pars[ispGOwin,,drop=FALSE]))
          # ST loosers
          if (n_accST == 0) stl <- 0 else {
            ispSTloss <- ispSGO & !ispGOacc
            stl <- apply(matrix(log(1-model$pfunG(
              rt=dadm$rt[ispSTloss]-pars[ispSTloss,"SSD"], # correct for SSD
              pars=pars[ispSTloss,,drop=FALSE])),
              nrow=n_accST),2,sum)
          }

          # Transform back to densities to include failures
          like[tGO] <- (1-gf[tGO])*(tf[tGO]*exp(like[tGO]) +
                                      (1-tf[tGO])*exp(ts+stl))
        }

        # ST WINS (never tf)
        if (any(ispSST)) { # Go triggered 1a) with (1-pStop), all race,
          #              1b) with pStop, only ST races
          # Go failure    2) Only ST races
          # ST winner rows
          ispSSTwin <-  dadm$winner &  ispSST
          # Stop probability on ST win trials, only needs go accumulators
          # model$sfun integrates over stop *duration* (SSRT); convert absolute RT to duration (RT-SSD).
          pStop <- model$sfun(pars[ispSST & ispGOacc,,drop=FALSE],n_acc=n_accG,
                              upper=dadm$rt[ispSSTwin] - dadm$SSD[ispSSTwin]) # pStop before observed RT
          # ST win ll
          tST <- ptrials[ispSSTwin]
          like[tST] <- log(model$dfunG(
            rt=dadm$rt[ispSSTwin]-dadm$SSD[ispSSTwin], # correct ST racers for SSD delay
            pars=pars[ispSSTwin,,drop=FALSE]))
          # ST looser survivors
          if (n_accST > 1) {  # Survivor for looser for ST accumulator(s)
            ispSSTloss <-  !dadm$winner &  ispSST & !ispGOacc
            llST <-  log(1-model$pfunG(
              rt=dadm$rt[ispSSTloss]-dadm$SSD[ispSSTloss],
              pars=pars[ispSSTloss,,drop=FALSE]))
            if (n_accST == 2) # Could remove branch, maybe faster as no matrix sum?
              like[tST] <- like[tST] + llST else
                like[tST] <- like[tST] + apply(matrix(llST,nrow=n_accST-1),2,sum)
          }
          # Go looser survivor
          ispSGloss <- ispSST & ispGOacc
          llG <- apply(matrix(log(1-model$pfunG(
            rt=dadm$rt[ispSGloss],pars=pars[ispSGloss,,drop=FALSE])),
            nrow=n_accG),2,sum)
          like[tST] <- (1-tf[tST])*(                 # Never trigger failure
            gf[tST]*exp(like[tST]) +             # Case 2, gf only ST race
              (1-gf[tST])*(pStop*exp(like[tST]) +      # Case 1b, no gf, stop beats go, ST race
                             (1-pStop)*exp(like[tST]+llG))) # Case 1a, no  gf, all race (no stop win)
        }
      }
      allLL[allok][allr] <- log(like)
    }

    allLL[is.na(allLL)|is.nan(allLL)] <- min_ll
    allLL <- pmax(min_ll,allLL)

    sum(allLL[attr(dadm,"expand")])
}


log_likelihood_race_ss_cens_trunc <- function(pars, dadm, model, min_ll = log(1e-10),
                                              normalise_trunc = FALSE) {
  if (is.null(attr(pars, "ok"))) {
    ok <- !logical(dim(pars)[1])
  } else {
    ok <- attr(pars, "ok")
  }
  if (!any(ok)) return(min_ll * length(attr(dadm, "expand")))

  n_acc <- length(levels(dadm$lR))
  n_trials_all <- nrow(dadm) / n_acc
  n_accG <- sum(as.numeric(dadm[1:n_acc, "lI"]) == 2)
  n_accST <- sum(as.numeric(dadm[1:n_acc, "lI"]) == 1)

  if (!("UC" %in% names(dadm)) && any(is.infinite(dadm$rt))) {
    stop("Censored SS likelihood requires dadm to have UC when rt contains Inf.")
  }
  if (normalise_trunc && !all(c("LT", "UT") %in% names(dadm))) {
    stop("normalise_trunc=TRUE requires dadm to include LT and UT columns.")
  }
  if ((normalise_trunc || any(is.infinite(dadm$rt))) && n_accST > 0) {
    if (normalise_trunc) {
      stop("normalise_trunc=TRUE not yet supported when stop-triggered accumulators are present (n_accST > 0).")
    }
  }

  allLL <- rep(min_ll, n_trials_all)
  allok_trial <- ok[dadm$lR == levels(dadm$lR)[1]]

  pars <- pars[ok, , drop = FALSE]
  dadm <- dadm[ok, , drop = FALSE]

  n_trials <- nrow(dadm) / n_acc
  trials <- 1:n_trials

  isp1 <- dadm$lR == levels(dadm$lR)[1]
  ispGOacc <- dadm$lI == levels(dadm$lI)[2]
  ispStop <- is.finite(dadm$SSD)

  gf <- pars[isp1, "gf"]
  tf <- pars[isp1, "tf"]

  GoR <- as.character(dadm[1:n_acc, "lR"][dadm[1:n_acc, "lI"] == 2])

  .prod_surv_go_min <- function(rt_upper, idx_go_rows) {
    if (length(rt_upper) == 0) return(numeric(0))
    rt_rep <- rep(rt_upper, each = n_accG)
    cdf <- model$pfunG(rt = rt_rep, pars = pars[idx_go_rows, , drop = FALSE])
    cdf[is.na(cdf)] <- 1
    cdf <- pmin(1, pmax(0, cdf))
    surv <- 1 - cdf
    surv_mat <- matrix(surv, nrow = n_accG)
    as.numeric(apply(surv_mat, 2, prod))
  }

  .prod_surv_ST_min <- function(rt_upper, idx_st_rows, ssd_upper) {
    if (length(rt_upper) == 0) return(numeric(0))
    if (n_accST == 0) return(rep(1, length(rt_upper)))
    # ST racers are shifted by SSD; evaluate their survivor at (rt_upper - ssd_upper)
    t_eff <- rt_upper - ssd_upper
    # If deadline precedes trigger, ST cannot respond yet
    t_eff[!is.finite(t_eff) | t_eff <= 0] <- 0
    t_rep <- rep(t_eff, each = n_accST)
    cdf <- model$pfunG(rt = t_rep, pars = pars[idx_st_rows, , drop = FALSE])
    cdf[is.na(cdf)] <- 1
    cdf <- pmin(1, pmax(0, cdf))
    surv <- 1 - cdf
    surv_mat <- matrix(surv, nrow = n_accST)
    as.numeric(apply(surv_mat, 2, prod))
  }

  .trunc_Z_stoptrial <- local({
    cache_env <- new.env(parent = emptyenv())
    function(trial_idx, lt, ut, idx_go_rows, idx_trial_rows) {
      if (lt == 0 && is.infinite(ut)) return(1)
      key <- paste0(
        formatC(lt, digits = 12, width = 0, format = "fg"), "_",
        formatC(ut, digits = 12, width = 0, format = "fg"), "_",
        paste(formatC(as.numeric(pars[idx_trial_rows, , drop = FALSE]), digits = 12, width = 0, format = "fg"),
              collapse = "_")
      )
      hit <- cache_env[[key]]
      if (!is.null(hit)) return(hit)

      go_pars_mat <- pars[idx_go_rows, , drop = FALSE]
      stop_pars <- pars[idx_trial_rows[1], , drop = FALSE]
      gf_j <- stop_pars[1, "gf"]
      tf_j <- stop_pars[1, "tf"]

      integrand <- function(t) {
        vapply(t, function(tt) {
          t_rep <- rep(tt, n_accG)
          pdfs <- model$dfunG(rt = t_rep, pars = go_pars_mat)
          cdfs <- model$pfunG(rt = t_rep, pars = go_pars_mat)
          pdfs[is.na(pdfs)] <- 0
          cdfs[is.na(cdfs)] <- 1
          cdfs <- pmin(1, pmax(0, cdfs))
          survs <- 1 - cdfs
          f_min <- 0
          for (k in seq_len(n_accG)) {
            f_min <- f_min + pdfs[k] * prod(survs[-k])
          }
          S_stop <- 1 - model$pfunS(rt = tt, pars = stop_pars)
          f_min * (tf_j + (1 - tf_j) * S_stop)
        }, numeric(1))
      }

      out <- try(integrate(integrand, lower = lt, upper = ut), silent = TRUE)
      Z <- if (inherits(out, "try-error")) 0 else out$value
      Z <- (1 - gf_j) * max(Z, .Machine$double.eps)
      cache_env[[key]] <- Z
      Z
    }
  })

  ispNR <- is.na(dadm$R)

  # No-response trials (R is NA). rt may be NA (intrinsic omission) or Inf (deadline-censored).
  isNR_trial <- is.na(dadm$R[isp1])
  if (any(isNR_trial)) {
    t_NR <- trials[isNR_trial]
    rt_NR <- dadm$rt[isp1][isNR_trial]
    SSD_NR <- dadm$SSD[isp1][isNR_trial]
    is_deadline <- is.infinite(rt_NR) & rt_NR > 0

    # Intrinsic omissions: preserve existing SS assumption.
    if (any(!is_deadline)) {
      t_intr <- t_NR[!is_deadline]
      is_stop_intr <- is.finite(SSD_NR[!is_deadline])
      allLL[allok_trial][t_intr[!is_stop_intr]] <- log(gf[t_intr[!is_stop_intr]])

      if (any(is_stop_intr)) {
        t_stop_intr <- t_intr[is_stop_intr]
        mask_trials <- rep(FALSE, n_trials)
        mask_trials[t_stop_intr] <- TRUE
        idx_go <- which(rep(mask_trials, each = n_acc) & ispGOacc)
        pStop <- pmin(1, pmax(0, model$sfun(pars[idx_go, , drop = FALSE], n_acc = n_accG)))
        allLL[allok_trial][t_stop_intr] <- log(gf[t_stop_intr] + (1 - gf[t_stop_intr]) * (1 - tf[t_stop_intr]) * pStop)
      }
    }

    # Deadline-censored no-response: probability no response by UC.
    if (any(is_deadline)) {
      t_dead <- t_NR[is_deadline]
      UC_dead <- dadm$UC[isp1][isNR_trial][is_deadline]
      SSD_dead <- SSD_NR[is_deadline]
      is_stop_dead <- is.finite(SSD_dead)

      mask_trials_dead <- rep(FALSE, n_trials)
      mask_trials_dead[t_dead] <- TRUE
      idx_go_dead_all <- which(rep(mask_trials_dead, each = n_acc) & ispGOacc)
      idx_st_dead_all <- which(rep(mask_trials_dead, each = n_acc) & (dadm$lI == levels(dadm$lI)[1]))
      S_go_UC_all <- .prod_surv_go_min(UC_dead, idx_go_dead_all)
      S_st_UC_all <- .prod_surv_ST_min(UC_dead, idx_st_dead_all, SSD_dead)
      # Stop survivor evaluated at trial time UC
      stop_pars_dead <- pars[isp1, , drop = FALSE][t_dead, , drop = FALSE]
      S_stop_UC_all <- 1 - model$pfunS(rt = UC_dead, pars = stop_pars_dead)
      S_stop_UC_all[is.na(S_stop_UC_all)] <- 0
      S_stop_UC_all <- pmin(1, pmax(0, S_stop_UC_all))

      # Go trials: gf + (1-gf)*S_go(UC)
      t_go <- t_dead[!is_stop_dead]
      if (length(t_go) > 0) {
        p_nr_go <- gf[t_go] + (1 - gf[t_go]) * S_go_UC_all[!is_stop_dead]
        allLL[allok_trial][t_go] <- log(pmax(p_nr_go, exp(min_ll)))
      }

      # Stop trials: deadline-censored NR by UC.
      # With stop-triggered accumulators, multiply by ST survivor by UC (ST is shifted by SSD).
      # pStop(UC) is the stop-win CDF up to UC in the model's convention (via sfun upper=UC).
      t_stop <- t_dead[is_stop_dead]
      if (length(t_stop) > 0) {
        UC_stop <- UC_dead[is_stop_dead]
        mask_trials_stop <- rep(FALSE, n_trials)
        mask_trials_stop[t_stop] <- TRUE
        idx_go_stop <- which(rep(mask_trials_stop, each = n_acc) & ispGOacc)
        # model$sfun integrates over stop *duration* (t) with go evaluated at t+SSD,
        # so an absolute deadline UC must be converted to an upper bound on duration (UC-SSD).
        upper_stop <- UC_stop - SSD_dead[is_stop_dead]
        pStop_UC <- pmin(1, pmax(0, model$sfun(pars[idx_go_stop, , drop = FALSE], n_acc = n_accG, upper = upper_stop)))

        S_go_stop <- S_go_UC_all[is_stop_dead]
        S_st_stop <- S_st_UC_all[is_stop_dead]
        S_stop_stop <- S_stop_UC_all[is_stop_dead]

        # If trigger failure: no stop process and ST not triggered => go-only.
        p_nr_if_tf <- gf[t_stop] + (1 - gf[t_stop]) * S_go_stop
        # If stop triggers: ST is active; NR requires ST survive, and either stop wins by UC or both go and stop survive by UC.
        p_nr_if_trig <- S_st_stop * (gf[t_stop] + (1 - gf[t_stop]) * (pStop_UC + S_go_stop * S_stop_stop))
        p_nr_stop <- tf[t_stop] * p_nr_if_tf + (1 - tf[t_stop]) * p_nr_if_trig
        allLL[allok_trial][t_stop] <- log(pmax(p_nr_stop, exp(min_ll)))
      }
    }
  }

  # Response trials: reuse SS likelihood and optionally apply truncation normalisation.
  if (any(!ispNR)) {
    allr <- !ispNR[isp1]
    pars_r <- pars[!ispNR, , drop = FALSE]
    dadm_r <- dadm[!ispNR, , drop = FALSE]

    n_trials_r <- nrow(dadm_r) / n_acc
    trials_r <- 1:n_trials_r
    ptrials <- rep(trials_r, each = n_acc)

    isp1_r <- dadm_r$lR == levels(dadm_r$lR)[1]
    ispGOacc_r <- dadm_r$lI == levels(dadm_r$lI)[2]
    ispStop_r <- is.finite(dadm_r$SSD)
    gf_r <- pars_r[isp1_r, "gf"]
    tf_r <- pars_r[isp1_r, "tf"]

    like <- numeric(n_trials_r)

    # Go trials with response
    if (any(!ispStop_r)) {
      ispGOwin <- !ispStop_r & dadm_r$winner
      tGO <- ptrials[ispGOwin]
      like[tGO] <- log(model$dfunG(rt = dadm_r$rt[ispGOwin], pars = pars_r[ispGOwin, , drop = FALSE]))
      if (n_accG > 1) {
        ispGOloss <- !ispStop_r & !dadm_r$winner & ispGOacc_r
        n_expected <- (n_accG - 1) * length(tGO)
        if (sum(ispGOloss) != n_expected) {
          cdf_loss <- model$pfunG(rt = dadm_r$rt[ispGOloss], pars = pars_r[ispGOloss, , drop = FALSE])
          cdf_loss[is.na(cdf_loss)] <- 1
          cdf_loss <- pmin(1, pmax(0, cdf_loss))
          ll_loss <- log(pmax(1 - cdf_loss, 0))
          loss_sum <- rowsum(ll_loss, group = ptrials[ispGOloss], reorder = FALSE)
          like[as.integer(rownames(loss_sum))] <- like[as.integer(rownames(loss_sum))] + loss_sum[, 1]
        } else {
        like[tGO] <- like[tGO] + apply(
          matrix(log(1 - model$pfunG(rt = dadm_r$rt[ispGOloss], pars = pars_r[ispGOloss, , drop = FALSE])),
                 nrow = n_accG - 1),
          2, sum
        )
        }
      }
      like[tGO] <- (1 - gf_r[tGO]) * exp(like[tGO])
    }

    # Stop trials with a response (n_accST==0 only)
    if (any(ispStop_r)) {
      ispSGO <- ispStop_r & (dadm_r$R %in% GoR)
      if (any(ispSGO)) {
        ispGOwin <- ispSGO & dadm_r$winner
        tGO <- ptrials[ispGOwin]
        like[tGO] <- log(model$dfunG(rt = dadm_r$rt[ispGOwin], pars = pars_r[ispGOwin, , drop = FALSE]))
        if (n_accG > 1) {
          ispGOloss <- ispSGO & !dadm_r$winner & ispGOacc_r
          n_expected <- (n_accG - 1) * length(tGO)
          if (sum(ispGOloss) != n_expected) {
            cdf_loss <- model$pfunG(rt = dadm_r$rt[ispGOloss], pars = pars_r[ispGOloss, , drop = FALSE])
            cdf_loss[is.na(cdf_loss)] <- 1
            cdf_loss <- pmin(1, pmax(0, cdf_loss))
            ll_loss <- log(pmax(1 - cdf_loss, 0))
            loss_sum <- rowsum(ll_loss, group = ptrials[ispGOloss], reorder = FALSE)
            like[as.integer(rownames(loss_sum))] <- like[as.integer(rownames(loss_sum))] + loss_sum[, 1]
          } else {
            like[tGO] <- like[tGO] + apply(
              matrix(log(1 - model$pfunG(rt = dadm_r$rt[ispGOloss], pars = pars_r[ispGOloss, , drop = FALSE])),
                nrow = n_accG - 1
              ),
              2, sum
            )
          }
          ts <- like[tGO] + log(1 - model$pfunS(rt = dadm_r$rt[ispGOwin], pars = pars_r[ispGOwin, , drop = FALSE]))
          like[tGO] <- (1 - gf_r[tGO]) * (tf_r[tGO] * exp(like[tGO]) + (1 - tf_r[tGO]) * exp(ts))
        }
      }
    }

    if (normalise_trunc) {
      LT1 <- dadm_r$LT[isp1_r]
      UT1 <- dadm_r$UT[isp1_r]
      for (j in seq_len(n_trials_r)) {
        if (LT1[j] == 0 && is.infinite(UT1[j])) next
        is_stop_j <- is.finite(dadm_r$SSD[which(isp1_r)[j]])
        start_row <- (j - 1) * n_acc + 1
        idx_trial_rows <- start_row:(start_row + n_acc - 1)
        idx_go_rows <- idx_trial_rows[ispGOacc_r[idx_trial_rows]]
        if (!is_stop_j) {
          S_lt <- .prod_surv_go_min(LT1[j], idx_go_rows)
          S_ut <- .prod_surv_go_min(UT1[j], idx_go_rows)
          Z <- (1 - gf_r[j]) * max(S_lt - S_ut, .Machine$double.eps)
        } else {
          Z <- .trunc_Z_stoptrial(j, LT1[j], UT1[j], idx_go_rows, idx_trial_rows)
        }
        like[j] <- like[j] / Z
      }
    }

    like[!is.finite(like)] <- exp(min_ll)
    like <- pmax(exp(min_ll), like)
    allLL[allok_trial][allr] <- log(like)
  }

  allLL[is.na(allLL) | is.nan(allLL)] <- min_ll
  allLL <- pmax(min_ll, allLL)
  sum(allLL[attr(dadm, "expand")])
}
