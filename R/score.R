# Out-of-sample scoring of held-out data --------------------------------------
#
# compare(WAIC = TRUE, LOO = TRUE) estimates out-of-sample fit by reweighting
# the fitted posterior: every trial in the emc object contributed to that
# posterior, and PSIS approximates what would have happened had it not.  The
# functions here score data the fit genuinely never saw, so nothing is
# approximated -- but the caller has to do the split themselves and fit on the
# training half first.
#
# Two quantities, selected by `type`, because held-out trials and held-out
# subjects are not the same prediction problem:
#
#   "trial"   the held-out rows belong to subjects that were in the fit, so
#             their alpha posteriors exist and the score is conditional:
#             elpd_i = log (1/S) sum_s p(y_i | alpha_j^(s)).
#
#   "subject" the held-out rows belong to people the fit never saw, so there is
#             no alpha for them and the subject likelihood must be marginalised
#             over the group level, exactly as .marg_ll_matrix() does for
#             leave-one-subject-out LOO:
#             elpd_j = log (1/S) sum_s [ (1/K) sum_k p(y_j | alpha^(k)),
#                                        alpha^(k) ~ N(theta_mu^(s), theta_var^(s)) ].
#
#   "calibrated"
#             also held-out subjects, but each one's trials are split: alpha_j is
#             estimated from the calibration half with the group level *frozen*
#             at the training posterior (refit_alpha()), and the remaining trials
#             are scored conditionally on that alpha.  "subject" asks whether the
#             fitted population predicts these people; "calibrated" asks whether
#             the model can accommodate them once it is allowed to learn who they
#             are, which is the question a marginal score cannot separate from
#             "were they typical of the training sample".

# Build a dadm for data the fit never saw, using the emc object's own design.
#
# attr(design, "data") holds the *training* data, but design_model() takes data
# as an explicit argument and only mapped_pars() reads that attribute, so there
# is nothing to swap: the design travels as-is and the held-out frame is passed
# in.  What does have to be repeated is the preprocessing make_emc() applies
# before it calls design_model() (see R/fitting.R), in particular copying the
# design's censoring/truncation bounds onto the data.
.held_out_design <- function(emc) {
  if (!is(emc, "emc")) stop("`emc` must be an emc object")
  design <- get_design(emc)
  if (length(design) > 1L || is.list(emc[[1]]$model))
    stop("Held-out scoring is not implemented for joint models")
  des <- design[[1]]
  if (!is.null(attr(des, "custom_ll")))
    stop("Held-out scoring is not available for designs with a custom likelihood")
  des
}

.held_out_data <- function(data) {
  if (is.data.frame(data)) {
    # fine
  } else if (is.list(data) && length(data) == 1L && is.data.frame(data[[1]])) {
    data <- data[[1]]
  } else {
    stop("`data` must be a data frame (or a list of length 1 containing one)")
  }
  if (is.null(data$subjects)) stop("`data` must have a `subjects` column")
  if (nrow(data) == 0) stop("`data` has no rows")
  data
}

# The held-out parameterisation has to be the fitted one, or the frozen group
# level and the fitted alphas do not refer to the same quantities.
.check_held_out_pars <- function(fit_pars, new_pars) {
  if (identical(as.character(new_pars), as.character(fit_pars))) return(invisible(NULL))
  only_fit <- setdiff(fit_pars, new_pars)
  only_new <- setdiff(new_pars, fit_pars)
  stop("Held-out design matrix does not match the fitted parameter vector.\n",
       if (length(only_fit)) paste0("  only in the fit: ",
                                    paste(only_fit, collapse = ", "), "\n") else "",
       if (length(only_new)) paste0("  only in the held-out data: ",
                                    paste(only_new, collapse = ", "), "\n") else "",
       "This usually means the fit was made with use_data = TRUE (the default), ",
       "which drops design columns that the training data never populated. ",
       "Refit with make_emc(..., use_data = FALSE) so the parameterisation is ",
       "fixed by the design rather than by the split.")
}

.held_out_dadm <- function(emc, data, rt_resolution = 1/60, verbose = FALSE) {
  des <- .held_out_design(emc)
  model <- emc[[1]]$model
  data <- .held_out_data(data)

  # Same preprocessing make_emc() does before design_model().
  data$subjects <- factor(as.character(data$subjects))
  data <- data[order(data$subjects), , drop = FALSE]
  # add_trials() renumbers `trials` within the held-out frame, destroying any
  # link back to the parent data set, so keep the row names first: they are the
  # only stable identity a pooled score can use to check that folds are
  # disjoint.  order() is stable and dm_list() splits on the same (sorted)
  # factor levels, so these stay aligned with the likelihood columns.
  held_out_ids <- rownames(data)
  data <- add_trials(data)
  for (nm in c("LT", "LC", "UT", "UC")) {
    if (!any(names(data) == nm) && !is.null(des$TC[[nm]])) data[[nm]] <- des$TC[[nm]]
  }

  # compress = FALSE keeps one likelihood column per held-out trial;
  # drop_unobserved = FALSE keeps the parameterisation design-determined rather
  # than letting the held-out cells decide which columns exist.
  dadm <- design_model(data = data, design = des, model = model,
                       compress = FALSE, rt_resolution = rt_resolution,
                       drop_unobserved = FALSE, verbose = verbose)

  .check_held_out_pars(emc[[1]]$par_names, attr(dadm, "sampled_p_names"))
  attr(dadm, "held_out_ids") <- held_out_ids
  dadm
}

# rt_resolution is consumed by make_emc() and never stored, so a caller can
# silently score at a different binning than they fitted at, which evaluates the
# held-out densities at differently rounded rts.  Compare against the training
# data instead of guessing.
.check_held_out_rt_resolution <- function(emc, rt_resolution) {
  model <- emc[[1]]$model
  # design_model() forces rt_resolution to NULL for grid-solve likelihoods
  # (R/design.R), so there is nothing the caller can get wrong for those.
  if (!is.null(model) && !model_compress_ok(model)) return(invisible(NULL))
  train <- emc[[1]]$data
  if (!is.list(train) || length(train) == 0) return(invisible(NULL))
  rt <- unlist(lapply(train, function(d) if (is.null(d$rt)) numeric(0) else d$rt),
               use.names = FALSE)
  rt <- rt[is.finite(rt) & rt > 0]
  if (length(rt) == 0) return(invisible(NULL))
  on_grid <- function(res) {
    if (is.null(res)) return(NA)
    q <- rt / res
    all(abs(q - round(q)) < 1e-8)
  }
  if (is.null(rt_resolution)) {
    if (isTRUE(on_grid(1/60)))
      warning("Scoring with rt_resolution = NULL, but the training rts all lie on a ",
              "1/60 grid, so the fit was probably made with rt_resolution = 1/60. ",
              "Pass the same value used for the fit.", call. = FALSE)
  } else if (isFALSE(on_grid(rt_resolution))) {
    warning("Scoring with rt_resolution = ", format(rt_resolution),
            ", but the training rts are not on that grid, so the fit used a ",
            "different (or no) binning. Pass the same value used for the fit.",
            call. = FALSE)
  }
  invisible(NULL)
}

# Frozen group level ----------------------------------------------------------
#
# The group level enters an alpha-only ("single") sampler as a fixed Gaussian
# prior, so freezing it is a matter of summarising p(theta_mu, theta_var | y_train)
# by two moments.  Which two depends on what the prior is meant to represent:
#
#   group_uncertainty = TRUE   the *predictive* distribution of a new subject's
#                              alpha, marginal over the group-level posterior.
#                              By the law of total covariance its first two
#                              moments are exactly mean(theta_mu) and
#                              mean(theta_var) + cov(theta_mu), so this is not an
#                              approximation of the moments, only a Gaussian
#                              approximation of a mixture of Gaussians.
#
#   group_uncertainty = FALSE  the plug-in / empirical-Bayes prior
#                              N(mean(theta_mu), mean(theta_var)), which treats
#                              the group level as known.  Narrower, so it shrinks
#                              the new alphas harder.
.frozen_group_prior <- function(emc, stage = "sample", filter = 0,
                                group_uncertainty = TRUE) {
  theta_mu  <- get_pars(emc, selection = "mu", stage = stage, filter = filter,
                        merge_chains = TRUE, return_mcmc = FALSE)
  theta_var <- get_pars(emc, selection = "Sigma", stage = stage, filter = filter,
                        merge_chains = TRUE, return_mcmc = FALSE)
  if (is.null(theta_mu) || is.null(theta_var))
    stop("Freezing the group level needs a hierarchical fit: this emc object has ",
         "no group-level parameters. (A type = \"single\" fit has nothing to freeze.)")
  mu  <- rowMeans(theta_mu)
  var <- apply(theta_var, c(1, 2), mean)
  if (group_uncertainty) var <- var + stats::var(t(theta_mu))
  dimnames(var) <- list(names(mu), names(mu))
  list(theta_mu_mean = mu, theta_mu_var = var,
       group_uncertainty = group_uncertainty, n_iter = ncol(theta_mu))
}

# Symmetrise and, if needed, nudge to positive definite: the averaged covariance
# is PD by construction, but cov(theta_mu) added on top can still trip the
# Cholesky in the proposal machinery at very low iteration counts.
.frozen_prior_pd <- function(var) {
  var <- (var + t(var)) / 2
  if (!inherits(try(chol(var), silent = TRUE), "try-error")) return(var)
  ev <- eigen(var, symmetric = TRUE)
  nm <- dimnames(var)
  var <- ev$vectors %*% diag(pmax(ev$values, max(ev$values) * 1e-8),
                             nrow = length(ev$values)) %*% t(ev$vectors)
  dimnames(var) <- nm
  var
}

#' Estimate Subject Parameters With the Group Level Frozen
#'
#' Fits subject-level parameters (\code{alpha}) for data the model was not
#' fitted to, holding the group level fixed at the posterior of an existing fit.
#' Each subject in \code{data} is estimated independently, with the training
#' fit's group level acting as a fixed prior rather than as something the new
#' data can change.
#'
#' This answers a question neither the posterior predictive nor
#' \code{score_held_out(type = "subject")} can: both of those integrate
#' \code{alpha} over the fitted population, so a held-out subject who is atypical
#' (say, with a much higher omission rate than anyone in the training sample)
#' scores badly whether or not the model's mechanisms could have produced them.
#' Freezing the group level and estimating \code{alpha} separates the two: if the
#' refitted subject's predictions match their data, the mechanism can accommodate
#' them and it was the population distribution that did not.
#'
#' Implemented as a \code{type = "single"} fit (see \code{\link{make_emc}}) whose
#' fixed prior is taken from \code{emc}, so subjects are independent of each
#' other and the group level cannot move. The prior is
#' \code{N(mean(theta_mu), mean(Sigma) + cov(theta_mu))} by default: the first two
#' moments of the predictive distribution of a new subject's \code{alpha},
#' marginal over the group-level posterior. With
#' \code{group_uncertainty = FALSE} it is the plug-in \code{N(mean(theta_mu),
#' mean(Sigma))}, which treats the group level as known and shrinks the new
#' subjects harder.
#'
#' The returned object is an ordinary emc object, so \code{\link{predict}},
#' \code{\link{plot_pars}}, \code{\link{plot_fit}} and friends work on it.
#' \code{\link{alpha_group_distance}} reports how far each refitted subject
#' landed in the tail of the frozen group level, which is the other half of the
#' answer: reproducing an atypical subject only at \code{alpha} four SD out of
#' the population says something quite different from reproducing them at one.
#'
#' Note that predictions from this fit are \emph{not} out-of-sample for the data
#' used to estimate \code{alpha}. To score it, hold some of each subject's trials
#' back -- which is what \code{score_held_out(type = "calibrated")} does.
#'
#' @param emc A hierarchical emc object fitted to the training data.
#' @param data A data frame of new subjects' data, in the same format as the data
#'   the model was fitted to. Subjects need not be new, but nothing stops the
#'   same person appearing in both.
#' @param stage Character. The sampling stage of \code{emc} to freeze the group
#'   level from, defaults to \code{"sample"}.
#' @param filter Integer or vector. Iterations of \code{emc} to remove, passed to
#'   \code{get_pars}.
#' @param group_uncertainty Boolean. If \code{TRUE} (the default), widen the
#'   frozen prior by the posterior uncertainty in \code{theta_mu}. See details.
#' @param rt_resolution The rt binning to build the new data model with; should
#'   match what \code{emc} was fitted with.
#' @param n_chains Integer. Number of chains for the refit.
#' @param iter Integer. Sampling iterations, passed to \code{\link{fit}}.
#' @param do_fit Boolean. If \code{FALSE}, return the unfitted emc object (with
#'   the frozen prior already installed) instead of running the sampler.
#' @param verbose Boolean. Passed to \code{\link{fit}}.
#' @param ... Further arguments for \code{\link{fit}}, e.g. \code{cores_per_chain},
#'   \code{cores_for_chains}, \code{stop_criteria}.
#' @return An emc object of \code{type = "single"} containing one \code{alpha}
#'   posterior per subject in \code{data}, whose \code{prior} is the frozen group
#'   level.
#' @examples \donttest{
#' # new_fit <- refit_alpha(emc, held_out_data)
#' # alpha_group_distance(new_fit)
#' # plot_fit(new_fit)
#' }
#' @export
refit_alpha <- function(emc, data, stage = "sample", filter = 0,
                        group_uncertainty = TRUE, rt_resolution = 1/60,
                        n_chains = 3, iter = 1000, do_fit = TRUE,
                        verbose = TRUE, ...) {
  des  <- .held_out_design(emc)
  data <- .held_out_data(data)
  frozen <- .frozen_group_prior(emc, stage = stage, filter = filter,
                                group_uncertainty = group_uncertainty)
  .check_held_out_rt_resolution(emc, rt_resolution)

  # use_data = FALSE so the parameterisation comes from the design rather than
  # from whichever cells these particular subjects happen to have filled.
  build <- function() make_emc(data, des, model = list(emc[[1]]$model),
                               type = "single", n_chains = n_chains,
                               rt_resolution = rt_resolution, use_data = FALSE)
  new_emc <- if (verbose) build() else suppressMessages(build())

  .check_held_out_pars(emc[[1]]$par_names, new_emc[[1]]$par_names)
  idx <- match(new_emc[[1]]$par_names, names(frozen$theta_mu_mean))
  if (anyNA(idx))
    stop("Group-level parameters do not cover the held-out parameter vector: ",
         paste(new_emc[[1]]$par_names[is.na(idx)], collapse = ", "))
  mu  <- frozen$theta_mu_mean[idx]
  var <- .frozen_prior_pd(frozen$theta_mu_var[idx, idx, drop = FALSE])

  for (i in seq_along(new_emc)) {
    new_emc[[i]]$prior$theta_mu_mean <- mu
    new_emc[[i]]$prior$theta_mu_var  <- var
    attr(new_emc[[i]]$prior, "frozen_group") <- TRUE
  }
  if (!do_fit) return(new_emc)
  fit(new_emc, iter = iter, verbose = verbose, ...)
}

#' Distance of Refitted Subjects From the Frozen Group Level
#'
#' For an emc object returned by \code{\link{refit_alpha}}, reports how unusual
#' each subject's estimated \code{alpha} is under the group level it was frozen
#' against: the per-parameter z score of the posterior mean, and the Mahalanobis
#' distance of the whole vector.
#'
#' This is the companion to a posterior predictive check on a refitted subject.
#' Reproducing an atypical subject's data is only evidence that the model
#' accommodates them if it does not require an \code{alpha} the fitted population
#' effectively excludes; the distances say which of the two happened.
#'
#' Under the frozen group level, a subject drawn from the population has squared
#' Mahalanobis distance distributed as chi-squared with \code{n_pars} degrees of
#' freedom, so \code{p} is the probability that a randomly drawn subject would be
#' at least this far out. It is a rough reference (the frozen prior is a Gaussian
#' summary, and posterior means are shrunk towards it) and not a test.
#'
#' @param emc An emc object from \code{\link{refit_alpha}}.
#' @param stage Character. Sampling stage to summarise, defaults to \code{"sample"}.
#' @param filter Integer or vector. Iterations to remove.
#' @return A data frame with one row per subject: the Mahalanobis distance
#'   \code{d}, its chi-squared tail probability \code{p}, and one z-score column
#'   per parameter.
#' @export
alpha_group_distance <- function(emc, stage = "sample", filter = 0) {
  prior <- emc[[1]]$prior
  if (!isTRUE(attr(prior, "frozen_group")))
    stop("`emc` does not carry a frozen group level; it should be the output of ",
         "refit_alpha().")
  mu  <- prior$theta_mu_mean
  var <- prior$theta_mu_var
  alpha <- get_pars(emc, selection = "alpha", stage = stage, filter = filter,
                    by_subject = TRUE, merge_chains = TRUE)
  means <- t(vapply(alpha, function(a) colMeans(do.call(rbind, a)),
                    numeric(length(mu))))
  means <- means[, names(mu), drop = FALSE]
  centred <- sweep(means, 2, mu, "-")
  d <- sqrt(rowSums((centred %*% solve(var)) * centred))
  z <- sweep(centred, 2, sqrt(diag(var)), "/")
  out <- data.frame(subjects = rownames(means), d = d,
                    p = stats::pchisq(d^2, df = length(mu), lower.tail = FALSE),
                    stringsAsFactors = FALSE)
  out <- cbind(out, as.data.frame(z))
  rownames(out) <- NULL
  out
}

# Split each subject's trials into a calibration half (used to estimate alpha
# with the group level frozen) and a scored half.  Chronological by default:
# "we have seen n trials from this person, how well do we predict the rest" is
# the situation the score is meant to stand in for.
.calibration_split <- function(data, n_calibration, order = c("first", "random")) {
  order <- match.arg(order)
  if (!is.numeric(n_calibration) || length(n_calibration) != 1 || n_calibration <= 0)
    stop("`n_calibration` must be a single positive number: a proportion of each ",
         "subject's trials if < 1, otherwise a number of trials")
  by_sub <- split(seq_len(nrow(data)), factor(as.character(data$subjects)))
  cal <- lapply(names(by_sub), function(s) {
    rows <- by_sub[[s]]
    n <- if (n_calibration < 1) max(1L, floor(n_calibration * length(rows))) else
      as.integer(n_calibration)
    if (n >= length(rows))
      stop("Subject ", s, " has ", length(rows), " trials, which leaves nothing ",
           "to score after ", n, " calibration trials. Lower `n_calibration`.")
    if (order == "random") rows <- sample(rows)
    rows[seq_len(n)]
  })
  cal <- sort(unlist(cal, use.names = FALSE))
  list(calibration = data[cal, , drop = FALSE],
       scored = data[setdiff(seq_len(nrow(data)), cal), , drop = FALSE])
}

# Conditional path: held-out trials from subjects that were in the fit.
.score_trials <- function(emc, dadm_list, stage, filter, cores) {
  model <- emc[[1]]$model
  alpha <- get_pars(emc, selection = "alpha", stage = stage, filter = filter,
                    by_subject = TRUE, merge_chains = TRUE)
  subs <- names(dadm_list)
  unknown <- setdiff(subs, names(alpha))
  if (length(unknown) > 0)
    stop("These subjects are not in the fit: ", paste(unknown, collapse = ", "),
         ".\nUse type = \"subject\" to score subjects the model has never seen.")
  ll_list <- auto_mclapply(subs, function(s) {
    proposals <- do.call(rbind, alpha[[s]])          # [n_iter x n_pars]
    calc_ll_pw(proposals, dadm_list[[s]], model)     # [n_iter x n_trials_s]
  }, mc.cores = cores)
  names(ll_list) <- subs
  n_trials <- vapply(ll_list, ncol, integer(1))
  list(ll_mat = do.call(cbind, ll_list),
       point_subject = rep(subs, n_trials))
}

# Marginal path: held-out subjects the fit has no alpha for.
.score_subjects <- function(emc, dadm_list, stage, filter, K, cores) {
  model <- emc[[1]]$model
  theta_mu  <- get_pars(emc, selection = "mu", stage = stage, filter = filter,
                        merge_chains = TRUE, return_mcmc = FALSE)
  theta_var <- get_pars(emc, selection = "Sigma", stage = stage, filter = filter,
                        merge_chains = TRUE, return_mcmc = FALSE)
  if (is.null(theta_mu) || is.null(theta_var))
    stop("type = \"subject\" needs a hierarchical fit: this emc object has no ",
         "group-level parameters to draw new subjects from.")

  n_iter <- ncol(theta_mu)
  subs   <- names(dadm_list)
  log_K  <- log(K)

  # Draw the proposals once per iteration and share them across subjects, so
  # every subject is scored against the same group-level draws (and a fit with
  # many iterations does not pay for the draw n_subjects times).
  props <- lapply(seq_len(n_iter), function(iter) {
    tryCatch(MASS::mvrnorm(K, theta_mu[, iter], theta_var[, , iter]),
             error = function(e) NULL)   # non-PD Sigma; this draw is skipped
  })
  n_bad <- sum(vapply(props, is.null, logical(1)))
  if (n_bad == n_iter)
    stop("Could not draw subject parameters from any posterior sample ",
         "(theta_var was never positive definite).")

  ll_list <- auto_mclapply(subs, function(s) {
    dadm_s <- dadm_list[[s]]
    out <- rep(NA_real_, n_iter)
    for (iter in seq_len(n_iter)) {
      if (is.null(props[[iter]])) next
      subj_lls <- rowSums(calc_ll_pw(props[[iter]], dadm_s, model))   # [K]
      out[iter] <- matrixStats::logSumExp(subj_lls) - log_K
    }
    out
  }, mc.cores = cores)

  ll_mat <- do.call(cbind, ll_list)                  # [n_iter x n_subjects]
  colnames(ll_mat) <- subs
  list(ll_mat = ll_mat, point_subject = subs, n_bad = n_bad)
}

#' Score Held-Out Data
#'
#' Computes the expected log pointwise predictive density (elpd) of data the
#' model was **not** fitted to. Unlike the WAIC and PSIS-LOO columns of
#' \code{\link{compare}}, which reweight the fitted posterior to approximate
#' out-of-sample fit, this evaluates genuinely new data: split the data
#' yourself, run \code{\link{fit}} on the training half, and pass the held-out
#' half here.
#'
#' Two prediction problems are supported, selected by \code{type}.
#'
#' \code{"trial"} scores held-out trials from subjects that \emph{were} in the
#' fit. Their subject-level parameters are known, so each trial is scored
#' conditionally on the posterior of that subject's \code{alpha}.
#'
#' \code{"subject"} scores subjects the fit has never seen. There is no
#' \code{alpha} for them, so each subject's likelihood is marginalised over the
#' group level: for every posterior draw of \code{(theta_mu, theta_var)},
#' \code{K} subject parameter vectors are drawn and the likelihood averaged.
#' This requires a hierarchical fit.
#'
#' \code{"calibrated"} also scores unseen subjects, but lets the model learn who
#' they are first: each subject's trials are split, \code{\link{refit_alpha}}
#' estimates their \code{alpha} from the calibration part with the group level
#' frozen at the training posterior, and the remaining trials are scored
#' conditionally on that \code{alpha}. The two subject-level scores answer
#' different questions -- \code{"subject"} asks whether the fitted population
#' predicts these people, \code{"calibrated"} asks whether the model can predict
#' them once it has seen some of their data, which is the question that survives
#' deliberately holding out atypical subjects. The gap between them is the cost
#' of not knowing who the subject is.
#'
#' Two things to get right before the numbers mean anything:
#'
#' \code{rt_resolution} is not recorded on the emc object, so it must be given
#' here and must match the value used for the fit (\code{\link{make_emc}}
#' defaults to \code{1/60}); scoring at a different binning evaluates the
#' held-out densities at differently rounded response times. A mismatch against
#' the training data is warned about but cannot always be detected.
#'
#' The fit should be made with \code{make_emc(..., use_data = FALSE)}. With the
#' default \code{use_data = TRUE}, design-matrix columns that the training data
#' never populated are dropped, so the parameter vector depends on the split and
#' need not match the held-out data's. This function checks and errors rather
#' than scoring against a mismatched parameterisation.
#'
#' @param emc An emc object fitted to the training data.
#' @param data A data frame of held-out data, in the same format as the data the
#'   model was fitted to.
#' @param type \code{"trial"} for held-out trials from fitted subjects,
#'   \code{"subject"} for entirely new subjects. See details.
#' @param rt_resolution The rt binning used when the model was fitted. Must
#'   match \code{make_emc}'s \code{rt_resolution} (default \code{1/60}); pass
#'   \code{NULL} if the fit used none.
#' @param stage Character. The sampling stage to draw from, defaults to
#'   \code{"sample"}.
#' @param filter Integer or vector. Iterations to remove, passed to
#'   \code{get_pars}.
#' @param K Integer. Number of subject parameter draws per posterior sample for
#'   \code{type = "subject"}. Defaults to 200.
#' @param n_calibration For \code{type = "calibrated"}, how many of each
#'   subject's trials are used to estimate their \code{alpha}: a proportion if
#'   less than 1 (the default, \code{0.5}), otherwise a number of trials. The
#'   rest are scored.
#' @param calibration_order For \code{type = "calibrated"}, whether the
#'   calibration trials are the \code{"first"} ones (the default, matching the
#'   applied situation of having seen a subject's first block) or a
#'   \code{"random"} subset (which removes any confound with time on task, e.g.
#'   practice or fatigue).
#' @param refit For \code{type = "calibrated"}, a list of further arguments for
#'   \code{\link{refit_alpha}}, e.g.
#'   \code{list(iter = 500, group_uncertainty = FALSE, cores_for_chains = 3)}.
#' @param cores Integer. Number of cores to parallelise subjects across.
#' @param verbose Boolean. Passed to \code{design_model} when building the
#'   held-out data model.
#' @return An object of class \code{emc.score}: a list with the total
#'   \code{elpd} and its standard error, the deviance-scaled \code{ic}
#'   (\code{-2 * elpd}, comparable with the DIC/BPIC/WAIC columns of
#'   \code{compare}), the \code{pointwise} elpd of each scored point, and the
#'   subject each point belongs to. For \code{type = "calibrated"} it also
#'   carries the frozen-group \code{refit} emc object.
#' @examples \donttest{
#' # split the data, fit the training half, score the rest
#' # dat <- get_data(samples_LNR)
#' # train <- dat[dat$trials <= 50, ]; test <- dat[dat$trials > 50, ]
#' # emc <- fit(make_emc(train, design_LNR, use_data = FALSE))
#' # score_held_out(emc, test, type = "trial")
#' }
#' @export
score_held_out <- function(emc, data, type = c("trial", "subject", "calibrated"),
                           rt_resolution = 1/60, stage = "sample", filter = 0,
                           K = 200, n_calibration = 0.5,
                           calibration_order = c("first", "random"),
                           refit = list(), cores = 1, verbose = FALSE) {
  type <- match.arg(type)
  calibration_order <- match.arg(calibration_order)

  # The calibrated path scores with a *different* emc object -- the frozen-group
  # refit -- but builds its data model from the original design, so the held-out
  # dadm is constructed identically to the other two paths.
  refit_emc <- NULL
  n_cal <- NULL
  if (type == "calibrated") {
    data <- .held_out_data(data)
    parts <- .calibration_split(data, n_calibration, calibration_order)
    n_cal <- table(factor(as.character(parts$calibration$subjects)))
    refit_emc <- do.call(refit_alpha,
                         c(list(emc = emc, data = parts$calibration, stage = stage,
                                filter = filter, rt_resolution = rt_resolution,
                                verbose = verbose), refit))
    data <- parts$scored
  }

  dadm <- .held_out_dadm(emc, data, rt_resolution = rt_resolution, verbose = verbose)
  .check_held_out_rt_resolution(emc, rt_resolution)
  dadm_list <- dm_list(dadm)

  res <- if (type == "subject") {
    .score_subjects(emc, dadm_list, stage = stage, filter = filter, K = K,
                    cores = cores)
  } else {
    # Conditional on alpha either way; for "calibrated" the alphas come from the
    # refit, whose own stage/filter are the refit's, not the training fit's.
    .score_trials(if (type == "calibrated") refit_emc else emc, dadm_list,
                  stage = if (type == "calibrated") "sample" else stage,
                  filter = if (type == "calibrated") 0 else filter,
                  cores = cores)
  }
  ll_mat <- res$ll_mat
  n_iter <- nrow(ll_mat)
  if (anyNA(ll_mat) && type != "subject")
    warning("Held-out log-likelihood contains NAs; they are dropped from the ",
            "posterior average.", call. = FALSE)

  # elpd_i = log mean_s p(y_i | theta_s), averaged over the retained draws.
  pointwise <- apply(ll_mat, 2, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_real_)
    matrixStats::logSumExp(x) - log(length(x))
  })
  n_points <- length(pointwise)
  elpd <- sum(pointwise, na.rm = TRUE)
  se   <- if (n_points > 1) sqrt(n_points) * stats::sd(pointwise, na.rm = TRUE) else NA_real_

  # Identity of each scored point, so combine_scores() can verify that folds do
  # not overlap: original row names for trials, subject names for subjects.
  point_id <- if (type == "subject") {
    res$point_subject
  } else {
    ids <- attr(dadm, "held_out_ids")
    if (length(ids) == n_points) ids else NULL
  }

  out <- list(elpd = elpd, se_elpd = se, ic = -2 * elpd,
              pointwise = pointwise, subject = res$point_subject,
              point_id = point_id,
              type = type, n_points = n_points, n_iter = n_iter,
              n_subjects = length(dadm_list), K = if (type == "subject") K else NA_integer_,
              n_skipped_draws = if (is.null(res$n_bad)) 0L else res$n_bad,
              n_calibration = n_cal, refit = refit_emc,
              n_folds = 1L)
  class(out) <- "emc.score"
  out
}

#' Combine Held-Out Scores Across Folds
#'
#' Pools the \code{\link{score_held_out}} results of several cross-validation
#' folds into a single score. Each fold must come from a model fitted without
#' the data that fold scores; this function does no fitting and no scoring, it
#' only aggregates results you already have.
#'
#' The cross-validated elpd is the sum of the fold elpds, but the standard error
#' is \emph{not} a combination of the per-fold standard errors: it is computed
#' from the pooled vector of pointwise values, \code{sqrt(n) * sd(pointwise)}
#' over all \code{n} scored points. Per-fold standard errors are close to
#' meaningless for \code{type = "subject"}, where a fold contributes only as many
#' points as it held out subjects.
#'
#' Folds must be disjoint or the same data is counted twice. This is checked
#' using the identity of each scored point (subject names for
#' \code{type = "subject"}, the held-out data's row names for
#' \code{type = "trial"}, which requires the folds to have been subset from a
#' common data frame with unique row names).
#'
#' @param scores A list of \code{emc.score} objects, all of the same
#'   \code{type}. Objects returned by this function can themselves be combined.
#' @param check_disjoint Boolean. If \code{TRUE} (the default), error when the
#'   folds score the same point more than once.
#' @return An \code{emc.score} object covering all folds, with \code{n_folds}
#'   set and a \code{folds} data frame giving each fold's contribution.
#' @examples \donttest{
#' # fit each fold on the other subjects, score the held-out ones, then pool
#' # scores <- Map(function(emc, held) score_held_out(emc, held, type = "subject"),
#' #               fold_fits, fold_data)
#' # combine_scores(scores)
#' }
#' @export
combine_scores <- function(scores, check_disjoint = TRUE) {
  if (is(scores, "emc.score")) scores <- list(scores)
  if (!is.list(scores) || length(scores) == 0)
    stop("`scores` must be a non-empty list of emc.score objects")
  ok <- vapply(scores, function(x) is(x, "emc.score"), logical(1))
  if (!all(ok))
    stop("Element(s) ", paste(which(!ok), collapse = ", "),
         " of `scores` are not emc.score objects (as returned by score_held_out)")

  types <- unique(vapply(scores, function(x) x$type, character(1)))
  if (length(types) > 1)
    stop("Cannot combine scores of different types: ",
         paste(types, collapse = " and "),
         ". Held-out trials and held-out subjects are different quantities.")
  type <- types

  pointwise <- unlist(lapply(scores, function(x) unname(x$pointwise)), use.names = FALSE)
  subject   <- unlist(lapply(scores, function(x) as.character(x$subject)), use.names = FALSE)
  ids <- lapply(scores, function(x) x$point_id)
  have_ids <- !vapply(ids, is.null, logical(1))

  if (check_disjoint) {
    if (!all(have_ids)) {
      warning("Fold(s) ", paste(which(!have_ids), collapse = ", "),
              " carry no point identities, so overlap between folds could not ",
              "be checked. Make sure the folds are disjoint.", call. = FALSE)
    } else {
      all_ids <- unlist(ids, use.names = FALSE)
      dup <- unique(all_ids[duplicated(all_ids)])
      if (length(dup) > 0) {
        what <- if (type == "subject") "subject(s)" else "row(s) of the data"
        stop("Folds overlap: ", what, " ", paste(utils::head(dup, 5), collapse = ", "),
             if (length(dup) > 5) paste0(" and ", length(dup) - 5, " more") else "",
             " are scored in more than one fold, which would count them twice.",
             if (type != "subject")
               "\n(For trial-level scores this test uses row names, so it also fires if the folds were not subset from a common data frame.)"
             else "")
      }
    }
  }

  n_points <- length(pointwise)
  elpd <- sum(pointwise, na.rm = TRUE)
  se   <- if (n_points > 1) sqrt(n_points) * stats::sd(pointwise, na.rm = TRUE) else NA_real_

  folds <- data.frame(
    fold = if (is.null(names(scores))) seq_along(scores) else names(scores),
    n_points = vapply(scores, function(x) as.integer(x$n_points), integer(1)),
    n_subjects = vapply(scores, function(x) as.integer(x$n_subjects), integer(1)),
    elpd = vapply(scores, function(x) x$elpd, numeric(1)),
    stringsAsFactors = FALSE
  )

  out <- list(elpd = elpd, se_elpd = se, ic = -2 * elpd,
              pointwise = stats::setNames(pointwise,
                                          if (all(have_ids)) unlist(ids, use.names = FALSE) else NULL),
              subject = subject,
              point_id = if (all(have_ids)) unlist(ids, use.names = FALSE) else NULL,
              type = type, n_points = n_points,
              n_iter = max(vapply(scores, function(x) as.integer(x$n_iter), integer(1))),
              n_subjects = length(unique(subject)),
              K = scores[[1]]$K,
              n_skipped_draws = sum(vapply(scores, function(x)
                as.integer(x$n_skipped_draws), integer(1))),
              n_folds = sum(vapply(scores, function(x)
                as.integer(x$n_folds %||% 1L), integer(1))),
              folds = folds)
  class(out) <- "emc.score"
  out
}

#' Score Cross-Validation Folds
#'
#' Convenience wrapper for the case where the fold fits already exist: scores
#' each fold's held-out data against that fold's model and pools the result with
#' \code{\link{combine_scores}}. No fitting happens here -- \code{emcs[[i]]} must
#' already have been fitted to everything except \code{data[[i]]}.
#'
#' With one model per subject-fold and \code{type = "subject"} this is K-fold
#' cross-validation at the subject level (leave-one-subject-out when each fold
#' holds out one person). Fit every fold with
#' \code{make_emc(..., use_data = FALSE)} so all folds share one parameter
#' vector; otherwise folds can drop different design columns and
#' \code{score_held_out} will (correctly) refuse to score them.
#'
#' @param emcs A list of fitted emc objects, one per fold.
#' @param data A list of held-out data frames, the same length as \code{emcs},
#'   where \code{data[[i]]} is the data \code{emcs[[i]]} was \emph{not} fitted to.
#' @param check_disjoint Boolean, passed to \code{\link{combine_scores}}.
#' @param ... Passed to \code{\link{score_held_out}} (\code{type},
#'   \code{rt_resolution}, \code{K}, \code{cores}, ...). \code{rt_resolution}
#'   must match what the folds were fitted with.
#' @return An \code{emc.score} object covering all folds.
#' @examples \donttest{
#' # score_folds(fold_fits, fold_data, type = "subject", rt_resolution = 1/60)
#' }
#' @export
score_folds <- function(emcs, data, check_disjoint = TRUE, ...) {
  if (is(emcs, "emc")) emcs <- list(emcs)
  if (is.data.frame(data)) data <- list(data)
  if (length(emcs) != length(data))
    stop("`emcs` and `data` must be the same length: one held-out data set per ",
         "fold fit (got ", length(emcs), " and ", length(data), ")")
  scores <- vector("list", length(emcs))
  for (i in seq_along(emcs)) {
    scores[[i]] <- tryCatch(score_held_out(emcs[[i]], data[[i]], ...),
                            error = function(e)
                              stop("Fold ", i, ": ", conditionMessage(e), call. = FALSE))
  }
  names(scores) <- names(emcs)
  combine_scores(scores, check_disjoint = check_disjoint)
}

#' @export
print.emc.score <- function(x, digits = 1, ...) {
  what <- switch(x$type,
                 trial = "held-out trials",
                 subject = "held-out subjects",
                 calibrated = "held-out trials, new subjects calibrated")
  n_folds <- x$n_folds %||% 1L
  cat("Held-out score (", what,
      if (n_folds > 1) paste0(", ", n_folds, " folds") else "", ")\n", sep = "")
  cat("  points scored : ", x$n_points, " from ", x$n_subjects,
      " subject(s), ", x$n_iter, " posterior draws\n", sep = "")
  if (x$type == "subject")
    cat("  marginalised over ", x$K, " subject draws per posterior sample\n", sep = "")
  if (x$type == "calibrated" && !is.null(x$n_calibration))
    cat("  alpha calibrated on ",
        if (length(unique(x$n_calibration)) == 1) as.integer(x$n_calibration[1]) else
          paste0(min(x$n_calibration), "-", max(x$n_calibration)),
        " trials per subject, group level frozen\n", sep = "")
  if (x$n_skipped_draws > 0)
    cat("  skipped draws : ", x$n_skipped_draws, " (non-PD covariance)\n", sep = "")
  cat("  elpd          : ", round(x$elpd, digits),
      if (is.na(x$se_elpd)) "" else paste0(" (SE ", round(x$se_elpd, digits), ")"),
      "\n", sep = "")
  cat("  -2*elpd       : ", round(x$ic, digits), "\n", sep = "")
  if (!is.null(x$folds)) {
    cat("\n")
    folds <- x$folds
    folds$elpd <- round(folds$elpd, digits)
    print(folds, row.names = FALSE)
  }
  invisible(x)
}
