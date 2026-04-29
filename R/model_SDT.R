
pPROBIT <- function(lt,ut,pars)
  # probability between lt and ut
{
  pnorm(ut,mean=pars[,"mean"],sd=pars[,"sd"]) - pnorm(lt,mean=pars[,"mean"],sd=pars[,"sd"])
}


rPROBIT <- function(lR,pars,p_types=c("mean","sd","threshold"),lt=-Inf)
  # lR is an empty latent response factor lR with one level for response.
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in
  # contiguous rows.

{
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  nr <- length(levels(lR)) # Number of responses
  n <- dim(pars)[1]/nr     # Number of simulated trials
  first <- seq(1,dim(pars)[1]-nr+1,length.out=n)  # pick out mean and sd
  threshold <- matrix(pars[,"threshold"],nrow=nr) # format thresholds
  threshold[dim(threshold)[1],] <- Inf
  pmat <- rbind(rnorm(n,pars[first,"mean"],pars[first,"sd"]), # sample normal
                rep(lt,dim(threshold)[2]),threshold)          # lt ... ut
  pmat[dim(pmat)[1],] <- Inf
  R <- factor(apply(pmat,2,function(x){.bincode(x[1],x[-1])}),
              levels=1:length(levels(lR)),labels=levels(lR))
  cbind.data.frame(R=R,rt=NA)
}

#' Gaussian Signal Detection Theory Model for Binary Responses
#'
#' Discrete binary choice based on continuous Gaussian latent, with no rt (rt
#' must be set to NA in data).
#'
#' Model parameters are:
#'    mean (unbounded)
#'    sd (log scale) and
#'    threshold (unbounded).
#'
#' For identifiability in one condition two parameters must be fixed
#' (conventionally mean=0 and sd = 1). When used with data that records only
#' accuracy (so reponse bias cannot be evaluated) a single threshold must be
#' assumed and fixed (e.g., threshold = 0).
#'
#' At present this model is not fully implemented in C, but as its likelihood
#' requires only pnorm evaluation it is quite fast.
#'
#' @return A model list with all the necessary functions to sample
#' @examples
#' dprobit <- design(Rlevels = c("left","right"),
#'            factors=list(subjects=1,S=c("left","right")),
#'            formula=list(mean ~ 0+S, sd ~ 1,threshold ~ 1),
#'            matchfun=function(d)d$S==d$lR,
#'            constants=c(sd=log(1),threshold=0),
#'            model=SDT)
#'
#' p_vector <- sampled_pars(dprobit)
#' @export

SDT <- function(){
  list(
  type="SDT",
  p_types=c("mean" = 0,"sd" = log(1),"threshold" = 0),
  # Trial dependent parameter transform
  transform=list(func=c(mean = "identity",sd = "exp",threshold="identity")),
  bound=list(minmax=cbind(mean=c(-Inf,Inf),sd = c(0, Inf), threshold=c(-Inf,Inf))),
  Ttransform = function(pars,dadm) {
    pars
  },
  # Random function for discrete choices
  rfun=function(data=NULL,pars) {
    rPROBIT(data$lR,pars)
  },
  # probability of choice between lower and upper thresholds (lt & ut)
  pfun=function(lt,ut,pars) pPROBIT(lt,ut,pars),
  # quantile function, p = probability, used in making linear ROCs
  qfun=function(p) qnorm(p),
  # Likelihood, lb is lower bound threshold for first response
  log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
  log_likelihood_sdt(pars=pars, dadm = dadm, model = model, min_ll = min_ll, lb=-Inf)
  })
  }

  #' Hierarchical Unequal-Variance Signal Detection (hUVSD) Model for Binary Responses
  #'
  #' Binary choice SDT model parameterized with sensitivity (d) and bias (c),
  #' following Lages (2024).
  #'
  #' Model parameters are:
  #'    d (sensitivity, distance between signal and noise means)
  #'    c (bias, deviation from the midpoint between means)
  #'    sd (standard deviation of the signal distribution, noise SD fixed at 1)
  #'
  #' @return A model list with all the necessary functions to sample
  #' @export
  hUVSD <- function(){
  list(
  type="SDT",
  c_name="hUVSD",
  p_types=c("d" = 0,"c" = 0,"sd" = log(1)),
  # Trial dependent parameter transform
  transform=list(func=c(d = "identity",c = "identity",sd = "exp")),
  bound=list(minmax=cbind(d=c(-Inf,Inf),c = c(-Inf,Inf), sd=c(0,Inf))),
  Ttransform = function(pars,dadm) {
  pars
  },
  # Random function for discrete choices
  rfun=function(data=NULL,pars) {
  # Assume level 2 is Signal and Yes
  is_signal <- data$S == levels(data$S)[2]
  m <- ifelse(is_signal, 0.5 * pars[,"d"], -0.5 * pars[,"d"])
  s <- ifelse(is_signal, pars[,"sd"], 1.0)
  p_yes <- pnorm((m - pars[,"c"]) / s)
  R <- factor(ifelse(runif(length(p_yes)) < p_yes, levels(data$R)[2], levels(data$R)[1]),
              levels=levels(data$R))
  cbind.data.frame(R=R,rt=NA)
  },
  # quantile function, p = probability, used in making linear ROCs
  qfun=function(p) qnorm(p),
  # Likelihood
  log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
  # R implementation for completeness/fallback
  is_signal <- dadm$S == levels(dadm$S)[2]
  chosen_yes <- dadm$R == levels(dadm$R)[2]
  mu <- ifelse(is_signal, 0.5 * pars[,"d"], -0.5 * pars[,"d"])
  s <- ifelse(is_signal, pars[,"sd"], 1.0)
  z <- (mu - pars[,"c"]) / s
  p_yes <- pnorm(z)
  p <- ifelse(chosen_yes, p_yes, 1-p_yes)
  # Handle expand attribute for unique trials
  p_uniq <- p[dadm$winner]
  sum(log(pmax(exp(min_ll), p_uniq[attr(dadm,"expand")])))
  })
  }

