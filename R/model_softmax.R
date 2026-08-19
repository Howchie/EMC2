
pSOFTMAX <- function(trials, pars)
  # Probability for each option in a trial.
{
  V <- pars[,'x']*pars[,'beta']
  Vmax <- ave(V, trials, FUN = max)
  expV <- exp(V - Vmax)

  denom <- ave(expV, trials, FUN = sum)

  expV / denom
}


log_likelihood_softmax <- function(pars, dadm, model, min_ll = log(1e-10)) {
  if (nrow(pars) != nrow(dadm)) {
    stop("`pars` and `dadm` must have the same number of rows.")
  }
  if (!all(c("x", "beta") %in% colnames(pars))) {
    stop("Softmax parameters must include `x` and `beta`.")
  }
  if (!"lR" %in% names(dadm) || !is.factor(dadm$lR)) {
    stop("Softmax data must contain a factor `lR` column.")
  }

  n_acc <- nlevels(dadm$lR)
  if (n_acc < 1L || nrow(pars) %% n_acc != 0L) {
    stop("Softmax data must contain complete accumulator blocks per trial.")
  }
  trial_id <- rep(seq_len(nrow(pars) / n_acc), each = n_acc)
  prob <- pSOFTMAX(trial_id, pars)
  winner <- if ("winner" %in% names(dadm)) dadm$winner else {
    if (!"R" %in% names(dadm)) stop("Softmax data must contain `R` or `winner`.")
    !is.na(dadm$R) & as.character(dadm$lR) == as.character(dadm$R)
  }
  winner <- !is.na(winner) & winner
  if (!any(winner)) return(min_ll)

  out <- sum(log(prob[winner]))
  if (!is.finite(out)) min_ll else max(out, min_ll)
}


rSOFTMAX <- function(lR,pars,p_types=c("x","beta"))

{
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  nr <- length(levels(lR)) # Number of responses
  n <- dim(pars)[1]/nr     # Number of simulated trials

  # trial index for each accumulator row
  trial_id <- rep(seq_len(n), each=nr)
  p <- pSOFTMAX(trial_id, pars)

  u <- runif(n)
  u_expanded <- u[trial_id]   # match accumulator rows

  # Grouped cumulative sums select the first option crossing the draw.
  cp <- cumsum(p)

  idx <- c(TRUE, diff(trial_id) != 0)
  cp[idx] <- p[idx]          # restart cumsum at trial start

  is_hit <- cp >= u_expanded
  hit_mat <- matrix(is_hit, nrow = nr, ncol = n)
  chosen <- max.col(t(hit_mat), ties.method = "first")
  chosen[colSums(hit_mat) == 0] <- nr # fallback (rare)

  # --- build response factor ---------------------------------------------------
  R_levels <- levels(lR)
  R <- factor(R_levels[chosen], levels=R_levels)

  # Return the sampled response and no response time.
  data.frame(R = R, rt = NA_real_)
}

#' Softmax Discrete Choice Model
#'
#' A multinomial (or binary) discrete-choice model using the softmax
#' (Luce/Logit) rule. The probability of choosing option \eqn{i} in trial
#' \eqn{t} is
#'
#' \deqn{
#'   P(R_t = i \mid x_{ti}, \beta_t)
#'      = \frac{\exp(\beta_t \, x_{ti})}
#'             {\sum_{j} \exp(\beta_t \, x_{tj})}.
#' }
#'
#' The model operates on *long-format* accumulator data, where each trial
#' consists of a block of rows—one per response option (accumulator). No
#' response time is modeled; \code{rt} must be \code{NA} in the input data.
#'
#' @section Parameters:
#' The softmax model uses the following parameters:
#'
#' \describe{
#'   \item{\code{x}}{A continuous latent Q-value or score for each response
#'        option. Unbounded.}
#'
#'   \item{\code{beta}}{Inverse temperature (precision) parameter. Must be
#'        non-negative. Larger values produce more deterministic choice.}
#' }
#'
#' Returned sampled data consist of a response factor \code{R} and an
#' \code{rt} column always set to \code{NA}.
#'
#' @return A model list with all the necessary functions to sample
#'
#' @examples
#' soft <- design(
#'   Rlevels = c("left","right"),
#'   factors = list(subjects = 1, S = c("left","right")),
#'   formula = list(x ~ 0 + S, beta ~ 1),
#'   matchfun = function(d) d$S == d$lR,
#'   constants = c(beta = log(1)),
#'   model = softmax
#' )
#'
#' # Sample parameter vector
#' p_vector <- sampled_pars(soft)
#'
#' @export
softmax <- function(){
  list(
    type="SOFTMAX",
    c_name = "SOFTMAX",
    p_types=c("x" = 0,"beta" = log(1)),
    # Trial dependent parameter transform
    transform=list(func=c(x = "identity", beta = "exp")),
    bound=list(minmax=cbind(x=c(-Inf,Inf),beta = c(0, Inf))),
    Ttransform = function(pars,dadm) {
      pars
    },
    # Random function for discrete choices
    rfun=function(data=NULL,pars) {
      rSOFTMAX(data$lR,pars)
    },
    # probability of choice between lower and upper thresholds (lt & ut)
    pfun=function(trials,pars) pSOFTMAX(trials,pars),
    # Likelihood, lb is lower bound threshold for first response
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
      log_likelihood_softmax(pars=pars, dadm = dadm, model = model, min_ll = min_ll)
    })
}
