#### Standard LBA ----


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

dBAwL <- function (rt, pars, posdrift = TRUE)
  # posdrift = truncated positive normal rates
  # robust slower, deals with extreme rate values
{
  dt <- rt - pars[,"t0"]
  ok <- (dt>0) & (pars[,"b"] >= pars[,"A"])
  ok[is.na(ok) | !is.finite(dt)] <- FALSE
  out <- numeric(length(dt))
  out[ok] <- dleakyba(t = dt[ok], A = pars[ok,"A"], b = pars[ok,"b"],
                         v = pars[ok,"v"], sv = pars[ok,"sv"], k = pars[ok, "k"],
                         posdrift = posdrift)
  out
}

pBAwL <- function (rt, pars, posdrift = TRUE)
  # posdrift = truncated positive normal rates
  # robust slower, deals with extreme rate values
{
  dt <- rt - pars[,"t0"]
  ok <- (dt>0) & (pars[,"b"] >= pars[,"A"])
  ok[is.na(ok) | !is.finite(dt)] <- FALSE
  out <- numeric(length(dt))
  out[ok] <- pleakyba(t = dt[ok], A = pars[ok,"A"], b = pars[ok,"b"],
                         v = pars[ok,"v"], sv = pars[ok,"sv"], k = pars[ok, "k"],
                         posdrift = posdrift)
  out
}


rLBA <- function(lR,pars,posdrift = TRUE,p_types=c("v","sv","b","A","t0"),
                 ok=rep(TRUE,length(lR)))
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
  bad <- apply(dt,2,function(x){all(is.infinite(x))})
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
  rt <- matrix(t0,nrow=nr)[pick] + dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  R[bad] <- NA
  rt[bad] <- Inf
  ok <- matrix(ok,nrow=length(levels(lR)))[1,]
  out$R[ok] <- levels(lR)[R][ok]
  out$R <- factor(out$R,levels=levels(lR))
  out$rt[ok] <- rt[ok]
  out
}

rBAwL <- function(lR,pars,ok=rep(TRUE,length(lR)),p_types=c("v","sv","b","A","t0","k"), posdrift = TRUE, eps = 1e-10) {
  bad <- rep(NA, length(lR)/length(levels(lR)))
  out <- data.frame(R = bad, rt = bad)
  nr <- length(levels(lR))
  dt <- matrix(Inf,nrow=nr,ncol=nrow(pars)/nr)
  t0 <- pars[,"t0"]
  pars <- pars[ok,]
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  lower <- if (posdrift) 0 else -Inf
  drifts <- msm::rtnorm(nrow(pars), mean = pars[,"v"], sd = pars[,"sv"], lower = lower)
  
  # For k<eps it's the LBA
  small_k = pars[,"k"]<eps
  if (any(small_k)) {
    dt[small_k] <- (pars[small_k,"b"]-pars[small_k,"A"]*runif(sum(small_k)))/
      drifts[small_k]
  }
  # leaky: ever-cross iff v > k*b
  big_k = pars[,"k"]>=eps & drifts > pars[,"k"] * pars[,"b"]
 
  if (any(big_k)) {
    num <- drifts[big_k] - pars[big_k,"k"] * pars[big_k,"b"]
    den <- drifts[big_k] - pars[big_k,"k"] * pars[big_k,"A"]*runif(sum(big_k))
    ratio <- num / den  # in (0,1)
    ratio <- pmin(pmax(ratio, .Machine$double.xmin), 1 - 1e-15)
    dt[big_k] <- (-1 / pars[big_k,"k"]) * log(ratio)
  }
  dt[dt<0] <- Inf
  bad <- apply(dt,2,function(x){all(is.infinite(x))})
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
  rt <- matrix(t0,nrow=nr)[pick] + dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  R[bad] <- NA
  rt[bad] <- Inf
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
    rfun=function(data,pars) rLBA(data$lR,pars,ok = attr(pars, "ok"),ifelse(posdrift,TRUE,FALSE)),
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dLBA(rt,pars,ifelse(posdrift,TRUE,FALSE)),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pLBA(rt,pars,ifelse(posdrift,TRUE,FALSE)),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
      log_likelihood_race_cens_trunc(pars=pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}

#' The Ballistic Accumulator Model with Leakage (modified race-variant of Brown & Heathcote, 2005)
#' @export
#'
BAwL <- function(posdrift=TRUE){
  list(
    type="RACE",
    c_name = ifelse(posdrift,"BAwL","BAwLIO"),
    # p_vector transform, sets sv as a scaling parameter
    p_types=c("v" = 1,"sv" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0), "k" = log(0),"pContaminant"=qnorm(0)),
    transform=list(func=c(v = "identity",sv = "exp", B = "exp", A = "exp",t0 = "exp",k = "exp",pContaminant="pnorm")),
    bound=list(minmax=cbind(v=c(-Inf,Inf),sv = c(0, Inf), A=c(1e-4,Inf),B=c(1e-4,Inf),t0=c(0.05,Inf),k=c(1e-4,Inf),pContaminant=c(0.001,0.999)),
               exception=c(A=0,pContaminant=0,k=0)),
    # Transform to natural scale
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
      pars
    },
    # Random function for racing accumulator
    rfun=function(data,pars) rBAwL(data$lR,pars,ok = attr(pars, "ok"),posdrift=posdrift,),
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dLBA(rt,pars,ifelse(posdrift,TRUE,FALSE)),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pLBA(rt,pars,ifelse(posdrift,TRUE,FALSE)),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
      log_likelihood_race_cens_trunc(pars=pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}


#' LBAGNG
#' LBA model with Go/No-Go racers
#' @export

LBAGNG <- function(){
  list(
    type="GNG",
    c_name = "LBAGNG",
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
    rfun=function(data,pars) rLBAGNG(data,pars,posdrift=TRUE,ok = attr(pars, "ok")),
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dLBA(rt,pars,posdrift = TRUE),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pLBA(rt,pars,posdrift = TRUE),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model){
      log_likelihood_race_cens_trunc(pars=pars, dadm=dadm, model=model, min_ll=log(1e-10))
    }
  )
}

rLBAGNG <- function(data,pars,p_types=c("v","sv","b","A","t0"),posdrift = TRUE,
                 ok=rep(TRUE,length(lR)))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in
  # contiguous rows.
{
  if (!"nogo"%in%levels(data$lR)) {
    stop("Go/No-Go models must have a nogo level for Response")
  }
  if (any(is.infinite(data$UC))) {
    stop("UC must be finite for Go/No-Go Models")
  }
  lR=data$lR
  D = data[data$lR==levels(data$lR)[1],"UC"]
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
  bad <- apply(dt,2,function(x){all(is.infinite(x))})
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
  rt <- matrix(t0,nrow=nr)[pick] + dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  R[bad] <- NA
  rt[bad] <- Inf
  ok <- matrix(ok,nrow=length(levels(lR)))[1,]
  out$R[ok] <- levels(lR)[R][ok]
  out$R <- factor(out$R,levels=levels(lR))
  out$rt[ok] <- rt[ok]
  out$R[out$rt>D]="nogo"
  out$rt[out$R=="nogo"]=Inf
  out
}

#' @export
#'

LogicalRulesLBA <- function(posdrift = TRUE, fast_path=TRUE, spline=FALSE){
  list(
    type="RACE",
    # Note: `calc_ll()` infers LBA `posdrift` from whether `c_name` contains "IO".
    # Keep this in sync so likelihood and simulation use the same setting.
    c_name = paste0("LBA_LogicalRules",ifelse(posdrift,"","IO"),ifelse(spline,"_spline",ifelse(fast_path,"_quad","_adapt"))),
    # p_vector transform, sets sv as a scaling parameter
    p_types=c("v" = 1,"sv" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0), "p"=qnorm(1),"q"=qnorm(0.5)),
    transform=list(func=c(v = "identity",sv = "exp", B = "exp", A = "exp",t0 = "exp",p="pnorm",q="pnorm")),
    bound=list(minmax=cbind(v=c(-Inf,Inf),sv = c(0, Inf), A=c(1e-4,Inf),B=c(0,Inf),t0=c(0.05,Inf),p=c(0.01,0.99),q=c(0.01,0.99)),
               exception=c(A=0,p=1,q=1)),
    # Transform to natural scale
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
      pars
    },
    # Random function for racing accumulator
    rfun=function(data,pars) rLBA(data$lR,pars,ok = attr(pars, "ok"),posdrift=ifelse(posdrift,TRUE,FALSE),),
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dLBA(rt,pars,posdrift=ifelse(posdrift,TRUE,FALSE)),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pLBA(rt,pars,posdrift=ifelse(posdrift,TRUE,FALSE)),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
      log_likelihood_redundant_target_race(pars=pars, dadm = dadm, model = model, min_ll = min_ll)
    }	
  )
}
