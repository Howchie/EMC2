# Last observation carried forward
# Replaces NA values with the last non-NA value
# Mimics zoo::na.locf behavior (vectorized for speed)
na_locf <- function(x, na.rm = FALSE) {
  if (length(x) == 0) return(x)

  # Check if all values are NA
  na_mask <- is.na(x)
  if (all(na_mask)) {
    if (na.rm) {
      return(x[0])  # Return empty vector
    } else {
      return(x)  # Return as is
    }
  }

  # Vectorized approach: create an index that tracks the last non-NA position
  # For each position, we need the index of the last non-NA value up to that point
  idx <- seq_along(x)
  idx[na_mask] <- NA  # Set NA positions to NA in index
  idx <- cummax(ifelse(na_mask, 0, idx))  # Cumulative max gives us last non-NA index

  # Replace values: use the index to look up the last non-NA value
  # For positions where idx is 0 (leading NAs), keep as NA
  result <- x
  non_zero <- idx > 0
  result[non_zero] <- x[idx[non_zero]]

  # If na.rm=TRUE, remove leading NAs
  if (na.rm) {
    first_non_na <- which(!na_mask)[1]
    if (!is.na(first_non_na)) {
      result <- result[first_non_na:length(result)]
    }
  }

  return(result)
}

.floor_to_rt_resolution <- function(x, rt_resolution=1/60) {
  if (is.null(rt_resolution)) return(x)
  finite <- is.finite(x)
  if (!any(finite)) return(x)

  # Avoid dropping exact-bin values one bin lower due to binary rounding, e.g.
  # floor(0.7 / 0.05) == 13 instead of 14 on typical platforms.
  q <- x[finite] / rt_resolution
  tol <- 16 * .Machine$double.eps * pmax(1, abs(q))
  x[finite] <- floor(q + tol) * rt_resolution
  x
}

# Whether a model tolerates having its response times binned.
#
# Binning rt (rt_resolution) is a lossy compression: the likelihood then
# evaluates the density AT the bin edge instead of integrating over the bin.
# For closed-form models that is worth it, because collapsing duplicate rows is
# what makes the likelihood cheap.  For models whose likelihood is a grid solve
# cached per parameter tuple (RLF, ROU) it is pure loss: the cost is the time
# march out to max(rt), and the rt values only select readout points along a
# solve that has already been paid for, so binning saves nothing and biases the
# parameter that absorbs the bin-edge error.  Such a model declares
# compress_ok = FALSE; anything that says nothing is compressible.
model_compress_ok <- function(model) {
  ok <- model()$compress_ok
  is.null(ok) || isTRUE(ok)
}

# Names of the contaminant nuisance parameters, in the order they are appended.
# Both are *trailing* p_types: design() turns any p_type absent from the formula
# into a constant at its default, and emc2col::validate_col_prefix (src/col_registry.h)
# only checks the canonical prefix of the column order, so appending here is free.
.nuisance_par_names <- c("pContaminant", "pGuess")

# Append the contaminant nuisance parameters to a model's parameter machinery.
#
#   pContaminant  Bernoulli *omission* rate -- mass at rt == +Inf only.
#   pGuess        uniform-outlier ("guess") rate -- a flat density over the
#                 guess window, contributing to observed RT densities.
#
# See src/contaminant_mixture.h for how the two combine (nested weights).  Both
# default to qnorm(0) == -Inf, i.e. probability 0, so a model that gains them
# is unchanged until the user puts them in a formula.
#
# `p_types`/`transform`/`exception` are named vectors and `minmax` a 2-row
# matrix, matching the shape each model constructor already builds.  Returns the
# four updated objects in a list.
add_nuisance_pars <- function(p_types, transform, minmax, exception = NULL,
                              which = .nuisance_par_names) {
  for (nm in which) {
    if (nm %in% names(p_types)) next
    p_types[[nm]] <- stats::qnorm(0)
    transform[[nm]] <- "pnorm"
    minmax <- cbind(minmax, c(0.001, 0.999))
    colnames(minmax)[ncol(minmax)] <- nm
    exception[[nm]] <- 0
  }
  list(p_types = p_types, transform = transform, minmax = minmax,
       exception = exception)
}
# Name of the operational-time warp parameter (Math/ballistic-time.md).  Like
# the nuisance parameters this is a *trailing* p_type kept out of
# p_types_canonical: eta = 0 is the exact parent model, so a design that never
# mentions it is unchanged and should not be warned about.
.time_warp_par_name <- "eta"

# Append the operational-time warp parameter to a ballistic model's parameter
# machinery.  Call this immediately BEFORE add_nuisance_pars() so the column
# order stays <model parameters>, eta, pContaminant, pGuess.  The kernels
# resolve eta by name, so the position is cosmetic.
add_time_warp_par <- function(p_types, transform, minmax, exception = NULL) {
  nm <- .time_warp_par_name
  if (!(nm %in% names(p_types))) {
    p_types[[nm]] <- 0
    transform[[nm]] <- "identity"
    # Unbounded, like v/mu/delta.  eta lives on the whole real line and the
    # numerics are total there; EMC2 bounds encode support, never numerical
    # convenience.
    minmax <- cbind(minmax, c(-Inf, Inf))
    colnames(minmax)[ncol(minmax)] <- nm
  }
  list(p_types = p_types, transform = transform, minmax = minmax,
       exception = exception)
}

# --- the warp itself, vectorised; this mirrors emc2tw:: in src/time_warp.h.
# eta == 0 short-circuits to the identity so that a design carrying a constant
# eta = 0 reproduces the parent model bit for bit.
.tw_active <- function(eta) (!is.na(eta) & eta != 0) | is.nan(eta)

.tw_fwd <- function(u, eta) {                       # s = c_eta(u)
  eta <- rep_len(eta, length(u))
  out <- u
  act <- .tw_active(eta) & !is.na(u) & is.finite(u) & u > 0
  if (!any(act)) return(out)
  om <- exp(eta[act])
  nan_eta <- is.nan(eta[act])
  l1p <- log1p(u[act])
  x <- om * l1p
  s <- ifelse(om < 1e-12, l1p, expm1(x) / om)   # omega -> 0 limit
  s[nan_eta] <- NaN
  s[(!nan_eta & !is.finite(om)) | x > 709] <- Inf
  out[act] <- s
  out
}

.tw_log_jac <- function(u, eta) {                 # log c'_eta(u)
  eta <- rep_len(eta, length(u))
  out <- numeric(length(u))
  act <- .tw_active(eta) & !is.na(u) & is.finite(u) & u > 0
  out[act] <- (exp(eta[act]) - 1) * log1p(u[act])
  out
}

.tw_jac <- function(u, eta) exp(.tw_log_jac(u, eta))

.tw_inv <- function(s, eta) {                     # u = c_eta^{-1}(s)
  eta <- rep_len(eta, length(s))
  out <- s
  act <- .tw_active(eta) & !is.na(s) & is.finite(s) & s > 0
  if (!any(act)) return(out)
  om <- exp(eta[act])
  ss <- s[act]
  os <- om * ss
  # log1p(om*s), with the overflow-safe log(om) + log(s) fallback.
  y <- ifelse(is.finite(os), log1p(os), log(om) + log(ss)) / om
  nan_eta <- is.nan(om)
  u <- ifelse(om < 1e-12, expm1(pmin(ss, 709)), expm1(y)) # omega -> 0 limit
  u[nan_eta] <- NaN
  u[!is.finite(om) & !nan_eta] <- 0              # omega -> Inf; NaN propagates
  u[is.finite(om) & y > 709] <- Inf
  u[om < 1e-12 & ss > 709] <- Inf
  out[act] <- u
  out
}

# eta column of a mapped-parameter matrix; zeros when the model has no warp.
.tw_eta <- function(pars) {
  if (!is.null(colnames(pars)) && .time_warp_par_name %in% colnames(pars))
    pars[, .time_warp_par_name]
  else rep(0, NROW(pars))
}


.apply_timed_guess_winner <- function(out, lR_levels) {
  if (is.null(out$R) || !("time" %in% lR_levels)) return(out)
  
  r_chr <- as.character(out$R)
  is_time <- !is.na(r_chr) & r_chr == "time"
  if (!any(is_time)) {
    out$R <- factor(r_chr, levels = lR_levels)
    return(out)
  }
  out$isTime = ifelse(is_time,TRUE,FALSE)
  guess_levels <- lR_levels[!lR_levels %in% c("time", "nogo")]
  if (length(guess_levels) == 0L) {
    out$R <- factor(r_chr, levels = lR_levels)
    return(out)
  }

  r_chr[is_time] <- sample(guess_levels, sum(is_time), replace = TRUE)
  out$R <- factor(r_chr, levels = lR_levels)
  out
}

.emc2_ll_cache_version <- 3L

.is_valid_ll_cache <- function(dadm, n_trials, n_lR, has_RACE_col) {
  if (!is.data.frame(dadm)) return(FALSE)
  if (!identical(attr(dadm, "emc2_ll_cache_version"), .emc2_ll_cache_version)) return(FALSE)
  if (is.null(attr(dadm, "emc2_all_finite_trials"))) return(FALSE)
  all_finite <- attr(dadm, "emc2_all_finite_trials")
  if (!is.logical(all_finite) || length(all_finite) != 1L || is.na(all_finite)) return(FALSE)

  if (has_RACE_col) {
    race_nacc <- attr(dadm, "RACE_nacc_by_row")
    race_mask <- attr(dadm, "RACE_mask")
    if (is.null(race_nacc) || is.null(race_mask)) return(FALSE)
    if (length(race_nacc) != n_trials || length(race_mask) != n_trials) return(FALSE)
  }

  if (isTRUE(all_finite)) return(TRUE)
  if (n_lR <= 0L || (n_trials %% n_lR) != 0L) return(FALSE)
  n_unique_trials <- n_trials %/% n_lR

  finite_rt_mask <- attr(dadm, "finite_rt_mask")
  finite_rt_unique <- attr(dadm, "finite_rt_unique_trial_indices")
  other_unique <- attr(dadm, "other_unique_trial_indices")
  active_nogo_trial_mask <- attr(dadm, "active_nogo_trial_mask")
  if (is.null(finite_rt_mask) || is.null(finite_rt_unique) || is.null(other_unique) ||
      is.null(active_nogo_trial_mask)) return(FALSE)
  if (length(finite_rt_mask) != n_trials) return(FALSE)
  if ((length(finite_rt_unique) + length(other_unique)) != n_unique_trials) return(FALSE)
  if (!is.logical(active_nogo_trial_mask) || length(active_nogo_trial_mask) != n_unique_trials) return(FALSE)

  valid_idx <- function(x) {
    if (!is.numeric(x) && !is.integer(x)) return(FALSE)
    if (length(x) == 0L) return(TRUE)
    all(is.finite(x)) && all(x == as.integer(x)) &&
      all(x >= 0L) && all(x <= (n_unique_trials - 1L))
  }
  if (!valid_idx(finite_rt_unique) || !valid_idx(other_unique)) return(FALSE)
  finite_rt_unique <- as.integer(finite_rt_unique)
  other_unique <- as.integer(other_unique)
  if (length(unique(finite_rt_unique)) != length(finite_rt_unique)) return(FALSE)
  if (length(unique(other_unique)) != length(other_unique)) return(FALSE)
  if (length(intersect(finite_rt_unique, other_unique)) > 0L) return(FALSE)
  if (!setequal(c(finite_rt_unique, other_unique), 0:(n_unique_trials - 1L))) return(FALSE)

  # Validate cache partition against current dadm content to guard against stale
  # attributes.  Every likelihood call runs this, so it has to be O(n) and
  # vectorised.  The per-trial loop it replaces cost 0.41 s per call on 20,000
  # censored trials -- more than the RLF grid solve it was guarding -- because
  # `j0 %in% finite_set` rescans a 20,000-element vector once per iteration, and
  # `R_idx[start_row]` dispatches [.factor 20,000 times.  Only data sets with a
  # non-finite trial reach here at all: `all_finite` returns above.
  rts <- dadm[["rt"]]
  R_idx <- dadm[["R"]]
  if (is.null(rts) || is.null(R_idx)) return(FALSE)
  lR_codes <- as.integer(dadm[["lR"]])
  nogo_code <- match("nogo", levels(dadm[["lR"]]))
  race_nacc <- if (has_RACE_col) attr(dadm, "RACE_nacc_by_row") else NULL
  race_mask <- if (has_RACE_col) attr(dadm, "RACE_mask") else NULL
  if (has_RACE_col && (is.null(race_nacc) || length(race_nacc) != n_trials)) return(FALSE)

  start_row <- seq.int(1L, n_trials, by = n_lR)
  nacc <- if (has_RACE_col) race_nacc[start_row] else rep.int(n_lR, n_unique_trials)
  if (any(!is.finite(nacc) | nacc < 1L | nacc > n_lR)) return(FALSE)

  # Rows are trial-major blocks of n_lR, so a trial index and a within-trial
  # accumulator index address every row without a loop.
  trial_of_row <- rep.int(seq_len(n_unique_trials), rep.int(n_lR, n_unique_trials))
  acc_of_row <- rep.int(seq_len(n_lR), n_unique_trials)
  in_race <- acc_of_row <= nacc[trial_of_row]

  # The disjointness and covering checks above already prove the two index sets
  # partition the trials, so membership in `other_unique` is exactly the
  # complement of membership in `finite_rt_unique` and one flag vector serves
  # both of the loop's `%in%` tests.
  finite_flag <- logical(n_unique_trials)
  finite_flag[finite_rt_unique + 1L] <- TRUE
  rt_start <- unname(rts[start_row])
  finite_trial <- is.finite(rt_start) & rt_start > 0 & !is.na(R_idx[start_row])
  if (!identical(unname(finite_trial), finite_flag)) return(FALSE)
  # A finite trial's own rows must be flagged; the loop asserted nothing about
  # the rest, so neither does this.
  if (!all(finite_rt_mask[in_race & finite_flag[trial_of_row]])) return(FALSE)

  if (is.na(nogo_code)) {
    if (any(active_nogo_trial_mask)) return(FALSE)
  } else {
    active <- in_race
    if (has_RACE_col) active <- active & race_mask
    is_nogo <- active & !is.na(lR_codes) & lR_codes == nogo_code
    have_nogo <- colSums(matrix(is_nogo, nrow = n_lR)) > 0
    if (!identical(unname(have_nogo), active_nogo_trial_mask)) return(FALSE)
  }
  TRUE
}

.cache_ll_data_attrs <- function(dadm, force_rebuild = FALSE) {
  if (!is.data.frame(dadm)) return(dadm)

  cols <- names(dadm)
  if (!all(c("lR", "rt", "R") %in% cols)) {
    attr(dadm, "emc2_ll_cache_version") <- .emc2_ll_cache_version
    return(dadm)
  }

  n_trials <- nrow(dadm)
  lR <- dadm[["lR"]]
  lR_codes <- as.integer(lR)
  n_lR <- length(unique(lR_codes))

  has_RACE_col <- "RACE" %in% cols && is.factor(dadm[["RACE"]])
  if (!isTRUE(force_rebuild) &&
      .is_valid_ll_cache(dadm, n_trials = n_trials, n_lR = n_lR, has_RACE_col = has_RACE_col)) {
    return(dadm)
  }

  if (isTRUE(force_rebuild)) {
    attr(dadm, "emc2_all_finite_trials") <- NULL
    attr(dadm, "finite_rt_mask") <- NULL
    attr(dadm, "finite_rt_unique_trial_indices") <- NULL
    attr(dadm, "other_unique_trial_indices") <- NULL
    attr(dadm, "active_nogo_trial_mask") <- NULL
    attr(dadm, "RACE_nacc_by_row") <- NULL
    attr(dadm, "RACE_mask") <- NULL
  }
  if (has_RACE_col) {
    race_idx <- dadm[["RACE"]]
    race_levels <- levels(race_idx)
    nacc_by_level <- suppressWarnings(as.integer(race_levels))
    if (anyNA(nacc_by_level)) stop("RACE column levels must be integer-valued (e.g., '2', '3').")

    race_nacc_by_row <- rep.int(as.integer(n_lR), n_trials)
    race_mask <- rep.int(TRUE, n_trials)
    race_codes <- as.integer(race_idx)
    for (i in seq_len(n_trials)) {
      code <- race_codes[i]
      if (is.na(code)) next
      nacc <- nacc_by_level[code]
      race_nacc_by_row[i] <- nacc
      lR_i <- lR_codes[i]
      if (!is.na(lR_i) && lR_i > nacc) race_mask[i] <- FALSE
    }
    attr(dadm, "RACE_nacc_by_row") <- race_nacc_by_row
    attr(dadm, "RACE_mask") <- race_mask
  }

  all_finite_trials <- TRUE
  if (n_trials > 0L) {
    if (n_lR <= 0L || (n_trials %% n_lR) != 0L) {
      all_finite_trials <- FALSE
    } else {
      start_idx <- seq.int(1L, n_trials, by = n_lR)
      rts <- dadm[["rt"]][start_idx]
      R_idx <- dadm[["R"]][start_idx]
      ok <- is.finite(rts) & rts > 0 & !is.na(R_idx)
      LT <- if ("LT" %in% cols) dadm[["LT"]][start_idx] else rep(0, length(start_idx))
      UT <- if ("UT" %in% cols) dadm[["UT"]][start_idx] else rep(Inf, length(start_idx))
      ok <- ok & (LT == 0) & is.infinite(UT)
      if (!all(ok)) all_finite_trials <- FALSE
    }
  }
  attr(dadm, "emc2_all_finite_trials") <- all_finite_trials

  if (!all_finite_trials && n_trials > 0L && n_lR > 0L && (n_trials %% n_lR) == 0L) {
    n_unique_trials <- n_trials %/% n_lR
    finite_rt_mask <- rep.int(FALSE, n_trials)
    finite_rt_unique_trial_indices <- integer(0)
    other_unique_trial_indices <- integer(0)
    active_nogo_trial_mask <- rep.int(FALSE, n_unique_trials)

    start_idx <- seq.int(1L, n_trials, by = n_lR)
    rts_dadm <- dadm[["rt"]]
    R_idxs_dadm <- dadm[["R"]]
    race_nacc_by_row <- if (has_RACE_col) attr(dadm, "RACE_nacc_by_row") else NULL
    race_mask <- if (has_RACE_col) attr(dadm, "RACE_mask") else NULL
    nogo_code <- match("nogo", levels(lR))

    for (j0 in 0:(n_unique_trials - 1L)) {
      start_row_idx <- start_idx[j0 + 1L]
      rt_j <- rts_dadm[start_row_idx]
      R_j <- R_idxs_dadm[start_row_idx]
      n_lR_j <- if (has_RACE_col && length(race_nacc_by_row) == n_trials) {
        race_nacc_by_row[start_row_idx]
      } else {
        n_lR
      }
      rows <- start_row_idx + 0:(n_lR_j - 1L)
      if (!is.na(nogo_code) && length(rows) > 0L) {
        active_rows <- if (has_RACE_col) rows[which(race_mask[rows])] else rows
        if (length(active_rows) > 0L) {
          active_nogo_trial_mask[j0 + 1L] <- any(lR_codes[active_rows] == nogo_code, na.rm = TRUE)
        }
      }
      if (is.finite(rt_j) && rt_j > 0 && !is.na(R_j)) {
        finite_rt_unique_trial_indices <- c(finite_rt_unique_trial_indices, j0)
        finite_rt_mask[rows] <- TRUE
      } else {
        other_unique_trial_indices <- c(other_unique_trial_indices, j0)
      }
    }

    attr(dadm, "finite_rt_mask") <- finite_rt_mask
    attr(dadm, "finite_rt_unique_trial_indices") <- finite_rt_unique_trial_indices
    attr(dadm, "other_unique_trial_indices") <- other_unique_trial_indices
    attr(dadm, "active_nogo_trial_mask") <- active_nogo_trial_mask
  }

  attr(dadm, "emc2_ll_cache_version") <- .emc2_ll_cache_version
  dadm
}

#
# augment <- function(s,da,design)
#   # Adds attributes to augmented data
#   # learn: empty array for Q values with dim = choice alternative (low,high) x
#   #   stimulus x trials (max across stimuli)
#   # index: look up (row number) for stimuli in da, a matrix dim  = max trials x
#   #   stimulus matrix (rows for each choice alternative contiguous)
# {
#   if (!is.null(design$adapt$stimulus)) {
#     targets <- design$adapt$stimulus$targets
#     par <- design$adapt$stimulus$output_name
#     maxn <- max(sapply(dimnames(targets)[[1]],function(x){table(da[da$subjects==s,x])}))
#     # da index x stimulus
#     out <- sapply(targets[1,],getIndex,cname=dimnames(targets)[[1]][1],
#                   da=da[da$subjects==s,],maxn=maxn)
#     stimulus <- list(index=out)
#     # accumulator x stimulus x trials
#     stimulus$learn <- array(NA,dim=c(dim(targets),maxn/dim(targets)[1]),
#                             dimnames=list(rownames(targets),targets[1,],NULL))
#     stimulus$targets <- targets
#     stimulus$par <- par
#   } # add other types here
#   list(stimulus=stimulus)
# }
#
# getIndex <- function(typei,cname,da,maxn) {
#   out <- which(da[,cname]==typei)
#   c(out,rep(NA,maxn-length(out)))
# }

# ============================================================================
# Shared Ballistic Accumulator (BA) Launch Helpers
# ============================================================================

.ba_launch_code <- function(drift_distribution, caller = "Model") {
  switch(drift_distribution,
         normal = 0L,
         lognormal = 1L,
         splitlognormal = 2L,
         weibull = 3L,
         stop("Unknown ", caller, " drift_distribution: ", drift_distribution))
}

# Public names of the launch pair, in kernel column order.  The plain
# lognormal launch (1) is sampled on the natural scale as the arithmetic mean
# and the coefficient of variation; the C++ kernels convert to (mu, sigma) via
# launch_lognormal_pair() in src/wald_functions.h.  The split-lognormal launch
# (2) is NOT reparameterized -- mu stays the exact median.  BTAwL has not been
# converted yet and passes meancv = FALSE to keep its (mu, sigma) columns.
.ba_par_names <- function(launch, meancv = TRUE) {
  if (launch == 1L) if (meancv) c("mean", "cv") else c("mu", "sigma")
  else if (launch == 2L) c("mu", "sigma", "delta")
  else if (launch == 3L) c("shape", "mean")
  else c("v", "sv")
}

# Weibull launch parameters use the shape and arithmetic mean.  R's
# `rweibull()` and the C++ kernels use the usual scale internally, so keep the
# conversion in one place.  Working on the log scale avoids overflowing the
# gamma function for the small shapes allowed by the model bounds.
.weibull_scale_from_mean <- function(shape, mean) {
  exp(log(mean) - lgamma(1 + 1 / shape))
}

.rweibull_mean <- function(n, shape, mean) {
  rweibull(n, shape, .weibull_scale_from_mean(shape, mean))
}
