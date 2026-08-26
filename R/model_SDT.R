
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
#' Discrete binary choice based on a continuous Gaussian latent variable. No
#' response time is modeled; `rt` must be `NA` in the data. For a binary
#' response, the model draws `X ~ Normal(mean, sd^2)` and compares it with the
#' criterion `threshold`: the first response level is selected below the
#' criterion and the second response level above it.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *mean* | identity | \[-Inf, Inf\] | 0 | | Mean of the Gaussian latent variable. |
#' | *sd* | log | \[0, Inf\] | log(1) | | Gaussian SD. |
#' | *threshold* | identity | \[-Inf, Inf\] | 0 | | Response criterion; for binary responses, the cutoff between the two response levels. |
#'
#' The model is invariant to a common location and scale transformation, so in
#' one condition two parameters must be fixed for identification, conventionally
#' `mean = 0` and `sd = 1`. When only accuracy is observed and response bias
#' cannot be estimated, the criterion must also be fixed (for example,
#' `threshold = 0`). The likelihood is the Gaussian probability between the
#' lower and upper response thresholds; the model's `qfun` is `qnorm` for ROC
#' construction.
#'
#' At present this model is not fully implemented in C, but as its likelihood
#' requires only pnorm evaluation it is quite fast.
#' The factor order matters: the first level of `R` is the lower response and
#' the second is the upper response, while `lR` supplies the latent response
#' thresholds used during simulation.
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
#' Binary-choice signal-detection model parameterized by sensitivity `d`,
#' criterion bias `c`, and the signal-distribution SD `sd`. No response time is
#' modeled; `rt` is always `NA`. The noise SD is fixed to 1 for scale
#' identification. With the current factor-level convention, the first level
#' of `S` is noise and the second is signal, while the first level of `R` is
#' "no" and the second is "yes".
#'
#' The latent distributions are
#'
#' * noise: `X | noise ~ Normal(-d / 2, 1^2)`;
#' * signal: `X | signal ~ Normal(d / 2, sd^2)`.
#'
#' A yes response is generated when `X > c`, so
#' `P(yes | noise) = Phi(-d / 2 - c)` and
#' `P(yes | signal) = Phi((d / 2 - c) / sd)`. Here `d` and `c` are estimated
#' on the identity scale and `sd` is estimated on the log scale.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Mapping** | **Interpretation** |
#' |---|---|---|---|---|---|
#' | *d* | identity | \[-Inf, Inf\] | 0 | | Distance between signal and noise means. |
#' | *c* | identity | \[-Inf, Inf\] | 0 | | Criterion relative to the midpoint of the two means. |
#' | *sd* | log | \[0, Inf\] | log(1) | | Signal-distribution SD; noise SD is fixed at 1. |
#'
#' The model is hierarchical in the sense that the same subject-level
#' parameter vector can be used across signal/noise conditions; it does not
#' add a response-time hierarchy by itself. The R likelihood path is used
#' because the probabilities require only `pnorm` evaluations.
#'
#' @return A model list with all the necessary functions for sampling.
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
