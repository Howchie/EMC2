# Omission decomposition for fitted race-model emc objects.

.omx_key <- function(x) {
  if (!length(x)) "<none>" else paste(sort(x), collapse = "|")
}

.omx_setcol <- function(cols, value) {
  force(cols)
  force(value)
  function(pars) {
    hit <- intersect(cols, colnames(pars))
    if (length(hit)) pars[, hit] <- value
    pars
  }
}

# Model pfuns are closures over constructor options.  Cloning the closure
# environment lets the decomposition toggle the positive-drift switch without
# reconstructing the model (and therefore without losing gamma, rho, eta, or
# launch variants).
.omx_pfun <- function(model_list, posdrift = NULL) {
  f <- model_list$pfun
  if (!is.function(f)) return(NULL)
  if (is.null(posdrift)) return(f)
  c_name <- if (is.null(model_list$c_name)) "" else model_list$c_name
  # LBA() selects a literal branch with ifelse(), so its pfun does not retain
  # the constructor flag in the closure environment.
  if (grepl("^LBA(IO)?$", c_name)) {
    return(function(rt, pars) .lba_pfun(rt, pars, posdrift = isTRUE(posdrift)))
  }
  env <- environment(f)
  if (!is.environment(env) || !exists("posdrift", env, inherits = FALSE)) {
    return(NULL)
  }
  out <- f
  env2 <- list2env(as.list(env, all.names = TRUE), parent = parent.env(env))
  env2$posdrift <- isTRUE(posdrift)
  environment(out) <- env2
  out
}

.omx_family <- function(c_name) {
  if (grepl("^BAwL", c_name)) "BAwL"
  else if (grepl("^BAwDp", c_name)) "BAwDp"
  else if (grepl("^BAwDD", c_name)) "BAwDD"
  else if (grepl("^BAwD", c_name)) "BAwD"
  else if (grepl("^BAwF", c_name)) "BAwF"
  else if (grepl("^BAwR", c_name)) "BAwR"
  else if (grepl("^BTAwL", c_name)) "BTAwL"
  else if (grepl("^RDMSWTN", c_name)) "RDMSWTN"
  else if (grepl("^RDM", c_name)) "RDM"
  else if (grepl("^LBA", c_name)) "LBA"
  else "other"
}

.omx_mechanisms <- function(model_list, pars) {
  c_name <- if (is.null(model_list$c_name)) "" else model_list$c_name
  family <- .omx_family(c_name)
  live <- function(col, f = function(x) any(x > 0, na.rm = TRUE)) {
    col %in% colnames(pars) && isTRUE(f(pars[, col]))
  }
  m <- list()
  add <- function(name, desc, null = identity, posdrift = FALSE) {
    m[[name]] <<- list(name = name, desc = desc, null = null,
                       posdrift = posdrift)
  }

  # A kill clock is an omission mechanism in every family that exposes it.
  if (live("lambda_k")) {
    add("kill", "kill clock fires before threshold",
        null = .omx_setcol("lambda_k", 0))
  }

  # BAwL's leak defines the intrinsic endpoint defect; keep that mass in the
  # named asymptotic residual rather than registering a duplicate switch.

  # BAwD has two distinct sources of defective (never-finish) mass whenever
  # clearance is present.  Launches with V < ell never have positive
  # accumulation (the k = 0 limit), while launches with V > ell can rise but
  # still miss the threshold when the decaying drive turns around.  Mark the
  # residual for a direct split: the decomposition evaluates the k = 0 limit
  # and subtracts it from the full defective mass, rather than treating the
  # difference as a Shapley attribution.  Other ballistic families do not have
  # an independent ell clearance scale, so this split is specific to BAwD.
  if (identical(family, "BAwD")) {
    if (live("ell")) {
      attr(m, "split") <- list(
        name = "dead_launch",
        desc = "launch strength below clearance ell: never rises"
      )
      attr(m, "residual") <- c(
        name = "asymptotic_subthreshold",
        desc = "launch exceeds clearance but the turnaround remains below threshold"
      )
    }
  }

  # The IO constructors expose an untruncated normal drift.  Re-evaluating
  # their own pfun with positive drifts isolates negative accumulation.
  current_pos <- .omx_pfun(model_list, NULL)
  env <- if (is.function(current_pos)) environment(current_pos) else NULL
  pos_now <- if (is.environment(env) &&
                 exists("posdrift", envir = env, inherits = FALSE))
    get("posdrift", envir = env, inherits = FALSE) else NULL
  if (identical(pos_now, FALSE) && !is.null(.omx_pfun(model_list, TRUE))) {
    add("negative_drift", "untruncated non-positive accumulation",
        posdrift = TRUE)
  }

  # This residual is deliberately named: it covers finite-peak, finite-drive,
  # reservoir, or other model-specific sub-threshold mass without an opaque
  # residual bucket.
  if (is.null(attr(m, "residual"))) {
    attr(m, "residual") <- c(
      name = "asymptotic_subthreshold",
      desc = if (identical(family, "BAwL"))
        "launch strength below leak threshold (V < k b)"
      else
        "accumulation remains below threshold without a registered switch"
    )
  }
  m
}

.omx_alpha <- function(emc, n_post = 50,
                       stat = c("random", "mean", "median")) {
  stat <- match.arg(stat)
  if (stat == "random") {
    valid_n <- is.numeric(n_post) && length(n_post) == 1L &&
      is.finite(n_post) && n_post >= 1L &&
      n_post <= .Machine$integer.max && n_post == floor(n_post)
    if (!isTRUE(valid_n)) stop("n_post must be a positive integer")
  }
  raw <- get_pars(
    emc, selection = "alpha", merge_chains = TRUE,
    by_subject = TRUE, remove_constants = FALSE, remove_dup = FALSE,
    return_mcmc = FALSE
  )
  if (is.list(raw) && length(raw) == 1L) raw <- raw[[1L]]
  d <- dim(raw)
  if (is.null(d) || length(d) < 2L)
    stop("alpha posterior samples are not available; run fit() first")
  if (length(d) == 2L) {
    # A one-draw/single-subject object can be returned without a third
    # dimension by older emc objects.
    raw <- array(raw, dim = c(d[1L], 1L, d[2L]),
                 dimnames = list(dimnames(raw)[[1L]], NULL, NULL))
    d <- dim(raw)
  }
  pnames <- dimnames(raw)[[1L]]
  snames <- dimnames(raw)[[2L]]
  n_draw <- d[3L]
  as_matrix <- function(i) {
    z <- raw[, , i, drop = FALSE][, , 1L]
    out <- t(z)
    dimnames(out) <- list(snames, pnames)
    out
  }
  if (stat != "random") {
    out <- apply(raw, c(2L, 1L), stat)
    out <- matrix(out, nrow = d[2L], ncol = d[1L],
                  dimnames = list(snames, pnames))
    return(list(out))
  }
  lapply(seq_len(as.integer(n_post)), function(i)
    as_matrix(sample.int(n_draw, 1L)))
}

.omx_map <- function(alpha, dadm, model_list) {
  pars <- get_pars_oo(alpha, dadm, function() model_list)
  pars <- model_list$Ttransform(pars, dadm)
  fix_bound(pars, model_list$bound, dadm$lR)
}

.omx_prep <- function(emc, data = NULL) {
  design <- get_design(emc)
  if (length(design) != 1L)
    stop("decompose_omissions() handles one design at a time")
  design <- design[[1L]]
  model <- design$model
  model_list <- model()
  if (!identical(model_list$type, "RACE"))
    stop("decompose_omissions() is defined for race models; got type '",
         model_list$type, "'")
  c_name <- if (is.null(model_list$c_name)) "" else model_list$c_name
  if (isTRUE(model_list$correlated) ||
      grepl("LogicalRules", c_name, fixed = TRUE))
    stop("decompose_omissions() currently requires independent accumulator races")
  if (is.null(data)) data <- get_data(emc)
  if (is.list(data) && !is.data.frame(data)) data <- data[[1L]]
  if (!is.data.frame(data)) stop("data must be a data frame")
  if (!all(c("R", "subjects") %in% names(data)))
    stop("data must contain R and subjects columns")
  data <- data[order(data$subjects), , drop = FALSE]

  expanded <- add_accumulators(
    add_trials(data), design$matchfun, simulate = TRUE,
    type = model_list$type, Fcovariates = design$Fcovariates,
    fixed_accumulator_roles = design$fixed_accumulator_roles
  )
  dadm <- design_model(
    expanded, design, model, add_acc = FALSE, compress = FALSE,
    verbose = FALSE, rt_check = FALSE, rt_resolution = NULL
  )
  n_acc <- length(levels(dadm$lR))
  if (!n_acc || nrow(dadm) %% n_acc)
    stop("could not identify complete accumulator trials")
  expected <- rep(seq_len(n_acc), nrow(dadm) / n_acc)
  if (!identical(as.integer(dadm$lR), expected))
    stop("unexpected accumulator ordering in expanded design")
  first <- dadm$lR == levels(dadm$lR)[1L]
  trial <- dadm[first, , drop = FALSE]

  tc <- design$TC
  if (is.null(tc)) tc <- list()
  getcol <- function(nm, default) {
    if (!is.null(dadm[[nm]])) return(as.numeric(trial[[nm]]))
    if (!is.null(tc[[nm]])) return(rep_len(as.numeric(tc[[nm]]), nrow(trial)))
    rep(default, nrow(trial))
  }
  UC <- getcol("UC", Inf)
  LC <- getcol("LC", 0)
  LT <- getcol("LT", 0)
  UT <- getcol("UT", Inf)
  if (any(LC > LT, na.rm = TRUE)) {
    stop("lower censoring (LC > LT) is not supported by decompose_omissions()")
  }
  if (any(UC > UT, na.rm = TRUE)) {
    stop("upper censoring (UC > UT) is not supported by decompose_omissions()")
  }
  UCresp <- if (!is.null(dadm$UCresponse)) {
    as.logical(trial$UCresponse)
  } else {
    rep_len(isTRUE(tc$UCresponse), nrow(trial))
  }
  UCresp[is.na(UCresp)] <- FALSE

  ref <- .omx_alpha(emc, n_post = 1L, stat = "mean")[[1L]]
  pars_ref <- .omx_map(ref, dadm, model_list)
  pC <- if ("pContaminant" %in% colnames(pars_ref))
    pars_ref[first, "pContaminant"] else rep(0, nrow(trial))

  list(design = design, model = model, model_list = model_list,
       data = data, dadm = dadm, trial = trial, n_acc = n_acc, first = first,
       UC = UC, LC = LC, LT = LT, UT = UT, UCresp = UCresp,
       pars_ref = pars_ref, pars_ref_pC = pC)
}

.omx_S <- function(t_trial, pars, pfun, n_acc) {
  n_trial <- length(t_trial)
  if (!n_trial) return(numeric())
  p <- pfun(rep(t_trial, each = n_acc), pars)
  p <- as.numeric(p)
  p[!is.finite(p)] <- 0
  p <- pmin(pmax(p, 0), 1)
  out <- apply(matrix(1 - p, nrow = n_acc), 2L, prod)
  pmin(pmax(as.numeric(out), 0), 1)
}

.omx_global_kill <- function(model_list) {
  c_name <- if (is.null(model_list$c_name)) "" else model_list$c_name
  isTRUE(grepl("_GLOBAL_KILL", c_name, fixed = TRUE))
}

.omx_erlang_shape <- function(model_list) {
  env <- if (is.function(model_list$pfun)) environment(model_list$pfun) else NULL
  if (is.environment(env) && exists("erlang_shape_cpp", env, inherits = FALSE))
    return(as.integer(get("erlang_shape_cpp", env, inherits = FALSE)))
  c_name <- if (is.null(model_list$c_name)) "" else model_list$c_name
  if (grepl("_E2", c_name, fixed = TRUE)) 2L else 1L
}

.omx_clock_survival <- function(t, lambda, shape, omega = 1) {
  if (!is.finite(lambda) || lambda <= 0) return(1)
  if (!is.finite(t)) return(0)
  if (shape == 1L) return(exp(-lambda * t))
  if (shape == 2L) {
    x <- lambda * t
    return(if (is.finite(x)) exp(-x) * (1 + x) else 0)
  }
  x <- lambda * t
  y <- 2 * x
  if (!is.finite(x) || !is.finite(y)) return(0)
  omega * exp(-lambda * t) +
    (1 - omega) * exp(-y) * (1 + y)
}

.omx_clock_density <- function(t, lambda, shape, omega = 1) {
  if (!is.finite(lambda) || lambda <= 0 || !is.finite(t)) return(0)
  x <- lambda * t
  if (!is.finite(x)) return(0)
  if (shape == 1L) return(lambda * exp(-x))
  if (shape == 2L) return(lambda * x * exp(-x))
  y <- 2 * x
  if (!is.finite(y)) return(0)
  omega * lambda * exp(-x) + (1 - omega) * 2 * lambda * y * exp(-y)
}

# A global kill clock is shared by all accumulators, so the product of the
# independently-killed marginal survivors is not the race survivor.  Average
# over the clock explicitly: S(t) = S0(t) P(K > t) + integral_0^t fK(u) S0(u) du.
.omx_S_global_kill <- function(t_trial, pars, pfun, n_acc, shape) {
  n_trial <- length(t_trial)
  if (!n_trial) return(numeric())
  if (!"lambda_k" %in% colnames(pars))
    stop("global kill decomposition requires a lambda_k column")
  base_pars <- pars
  base_pars[, "lambda_k"] <- 0
  out <- numeric(n_trial)
  for (j in seq_len(n_trial)) {
    rows <- ((j - 1L) * n_acc + 1L):(j * n_acc)
    pj <- base_pars[rows, , drop = FALSE]
    lambda <- pars[rows[1L], "lambda_k"]
    if (any(abs(pars[rows, "lambda_k"] - lambda) >
            1e-10 * pmax(1, abs(lambda))))
      stop("global kill requires lambda_k to be constant across accumulators")
    omega <- if (shape == 3L && "omega" %in% colnames(pars))
      pars[rows[1L], "omega"] else 1
    s0 <- function(u) .omx_S(u, pj, pfun, n_acc)[1L]
    if (!is.finite(lambda) || lambda <= 0) {
      out[j] <- s0(t_trial[j])
      next
    }
    integrand <- function(u) vapply(u, function(z)
      .omx_clock_density(z, lambda, shape, omega) * s0(z), numeric(1L))
    integral <- try(stats::integrate(integrand, lower = 0,
                                     upper = t_trial[j], rel.tol = 1e-6),
                    silent = TRUE)
    if (inherits(integral, "try-error") || !is.finite(integral$value)) {
      out[j] <- NA_real_
    } else {
      out[j] <- .omx_clock_survival(t_trial[j], lambda, shape, omega) *
        s0(t_trial[j]) + integral$value
    }
  }
  pmin(pmax(out, 0), 1)
}

.omx_off_state <- function(model_list, mech, pars) {
  p <- pars
  f <- model_list$pfun
  if (length(mech)) for (nm in names(mech)) {
    mm <- mech[[nm]]
    p <- mm$null(p)
    if (isTRUE(mm$posdrift)) {
      f2 <- .omx_pfun(model_list, TRUE)
      if (!is.null(f2)) f <- f2
    }
  }
  list(pars = p, pfun = f)
}

.omx_survivor_inf <- function(model_list, pars, pfun, n_trial, n_acc) {
  global_kill <- .omx_global_kill(model_list)
  erlang_shape <- .omx_erlang_shape(model_list)
  if (global_kill)
    .omx_S_global_kill(rep(Inf, n_trial), pars, pfun, n_acc, erlang_shape)
  else
    .omx_S(rep(Inf, n_trial), pars, pfun, n_acc)
}

.omx_Vfun <- function(model_list, mech, pars, n_trial, n_acc) {
  function(set) {
    off <- .omx_off_state(model_list, mech[setdiff(names(mech), set)], pars)
    .omx_survivor_inf(model_list, off$pars, off$pfun, n_trial, n_acc)
  }
}

# For BAwD, split the residual never-finish mass after all registered
# mechanisms have been ablated.  At k = 0, a trial can finish iff at least one
# launch exceeds ell, so the survivor is precisely the all-accumulator
# V < ell ("dead launch") event.  The difference from the full k survivor is
# the V > ell turnaround-asymptote event.
.omx_split_inf <- function(model_list, mech, pars, n_trial, n_acc) {
  split <- attr(mech, "split")
  if (is.null(split)) return(NULL)
  off <- .omx_off_state(model_list, mech, pars)
  s_none <- .omx_survivor_inf(model_list, off$pars, off$pfun,
                              n_trial, n_acc)
  p_dead <- .omx_setcol("k", 0)(off$pars)
  s_dead <- .omx_survivor_inf(model_list, p_dead, off$pfun,
                              n_trial, n_acc)
  list(dead_launch = pmin(pmax(s_dead, 0), 1),
       asymptotic_subthreshold = pmin(pmax(s_none - s_dead, 0), 1),
       residual = s_none)
}

.omx_shapley <- function(mech, Vfun) {
  nms <- names(mech)
  M <- length(nms)
  subsets <- lapply(0:(2^M - 1L), function(z) {
    nms[bitwAnd(z, 2L^(seq_len(M) - 1L)) > 0]
  })
  V <- lapply(subsets, Vfun)
  names(V) <- vapply(subsets, .omx_key, character(1L))
  phi <- setNames(vector("list", M), nms)
  for (m in nms) {
    others <- setdiff(nms, m)
    acc <- 0
    for (z in 0:(2^length(others) - 1L)) {
      T <- others[bitwAnd(z, 2L^(seq_along(others) - 1L)) > 0]
      w <- factorial(length(T)) * factorial(M - length(T) - 1L) /
        factorial(M)
      acc <- acc + w * (V[[.omx_key(c(T, m))]] - V[[.omx_key(T)]])
    }
    phi[[m]] <- acc
  }
  list(phi = phi, V = V, V_none = V[[.omx_key(character())]],
       V_all = V[[.omx_key(nms)]])
}

.omx_ablate_all <- function(prep, mech) {
  .omx_Vfun(prep$model_list, mech, prep$pars_ref,
            nrow(prep$trial), prep$n_acc)(character())
}

#' Which omission mechanisms are live in a fitted emc object.
#'
#' @param emc A fitted single-design race-model emc object.
#' @param data Optional data frame on which to evaluate the decomposition.
#' @return A data frame containing only live omission mechanisms, with
#'   `mechanism` and `description` columns.
#' @export
omission_mechanisms <- function(emc, data = NULL) {
  prep <- .omx_prep(emc, data)
  mech <- .omx_mechanisms(prep$model_list, prep$pars_ref)
  rows <- list()
  add_row <- function(name, is_live, description) {
    if (isTRUE(is_live)) rows[[length(rows) + 1L]] <<- data.frame(
      mechanism = name, description = description,
      stringsAsFactors = FALSE
    )
  }
  add_row("contaminant", any(prep$pars_ref_pC > 0, na.rm = TRUE),
          "pContaminant: omission before the accumulator race")
  if (length(mech)) for (mm in mech)
    add_row(mm$name, TRUE, mm$desc)
  residual <- attr(mech, "residual")
  split <- attr(mech, "split")
  if (is.null(split)) {
    add_row(residual[["name"]],
            any(.omx_ablate_all(prep, mech) > 1e-10, na.rm = TRUE),
            residual[["desc"]])
  } else {
    split_ref <- .omx_split_inf(prep$model_list, mech, prep$pars_ref,
                                nrow(prep$trial), prep$n_acc)
    add_row(split$name,
            any(split_ref$dead_launch > 1e-10, na.rm = TRUE), split$desc)
    add_row(residual[["name"]],
            any(split_ref$asymptotic_subthreshold > 1e-10, na.rm = TRUE),
            residual[["desc"]])
  }
  censor_window <- !prep$UCresp & is.finite(prep$UC) & prep$UC <= prep$UT
  censor_ref <- if (any(censor_window, na.rm = TRUE)) {
    s_ref <- function(t) {
      if (.omx_global_kill(prep$model_list))
        .omx_S_global_kill(t, prep$pars_ref, prep$model_list$pfun,
                           prep$n_acc, .omx_erlang_shape(prep$model_list))
      else
        .omx_S(t, prep$pars_ref, prep$model_list$pfun, prep$n_acc)
    }
    ifelse(censor_window, pmax(s_ref(prep$UC) - s_ref(prep$UT), 0), 0)
  } else numeric(nrow(prep$trial))
  add_row(
    "censor_slow",
    any(censor_ref > 1e-10, na.rm = TRUE),
    "response after UC and before UT (upper-censor omission)"
  )
  out <- if (length(rows)) do.call(rbind, rows) else
    data.frame(mechanism = character(), description = character(),
               stringsAsFactors = FALSE)
  attr(out, "c_name") <- prep$model_list$c_name
  out
}

#' Decompose expected omissions by generating mechanism.
#'
#' The returned probabilities are conditional on a trial being retained after
#' lower/upper truncation.  For BAwL, `asymptotic_subthreshold` is the intrinsic
#' mass with launch strength below the leak threshold (`V < k b`).
#' `censor_slow` is an upper-censor omission.  For BAwD with
#' positive clearance, that residual is split into `dead_launch` (all launches
#' below clearance, `V < ell`) and the turnaround-asymptote remainder
#' (`V > ell` but still never reaching threshold).
#'
#' @param emc A fitted single-design race-model emc object.
#' @param factors Character data columns used to form cells, or NULL.
#' @param n_post Number of posterior draws for `stat = "random"`.
#' @param stat Posterior summary: `"random"`, `"mean"`, or `"median"`.
#' @param by_subject Split cells by subject when TRUE.
#' @param data Optional data frame on which to evaluate the decomposition.
#' @param probs Quantiles to report in `summary`.
#' @return A `omission_decomposition` object with draws, summary, long,
#'   observed, coalitions, and mechanisms components.
#' @export
decompose_omissions <- function(emc, factors = NULL, n_post = 50,
                                stat = c("random", "mean", "median"),
                                by_subject = FALSE, data = NULL,
                                probs = c(.025, .5, .975)) {
  stat <- match.arg(stat)
  if (length(probs) < 1L || any(!is.finite(probs)) || any(probs < 0 | probs > 1))
    stop("probs must lie in [0, 1]")
  prep <- .omx_prep(emc, data)
  model_list <- prep$model_list
  dadm <- prep$dadm
  n_acc <- prep$n_acc
  first <- prep$first
  mech <- .omx_mechanisms(model_list, prep$pars_ref)
  residual_name <- attr(mech, "residual")[["name"]]

  if (!is.null(factors)) {
    miss <- setdiff(factors, names(prep$trial))
    if (length(miss)) stop("factors not in data: ", paste(miss, collapse = ", "))
  }
  cell_cols <- unique(c(if (isTRUE(by_subject)) "subjects", factors))
  cell <- if (length(cell_cols)) {
    interaction(prep$trial[, cell_cols, drop = FALSE], sep = " ", drop = TRUE)
  } else {
    factor(rep("all", nrow(prep$trial)))
  }

  alphas <- .omx_alpha(emc, n_post = n_post, stat = stat)
  if (stat != "random") n_post <- 1L else n_post <- length(alphas)
  UC <- prep$UC
  UT <- prep$UT
  LT <- prep$LT
  global_kill <- .omx_global_kill(model_list)
  erlang_shape <- .omx_erlang_shape(model_list)

  draws <- vector("list", n_post)
  coal_draws <- vector("list", n_post)
  for (i in seq_len(n_post)) {
    pars <- .omx_map(alphas[[i]], dadm, model_list)
    pC <- if ("pContaminant" %in% colnames(pars))
      pars[first, "pContaminant"] else rep(0, nrow(prep$trial))
    S <- function(t, p = pars, f = model_list$pfun) {
      if (global_kill)
        .omx_S_global_kill(t, p, f, n_acc, erlang_shape)
      else
        .omx_S(t, p, f, n_acc)
    }
    Vfun <- .omx_Vfun(model_list, mech, pars, nrow(prep$trial), n_acc)

    S_inf <- S(rep(Inf, nrow(prep$trial)))
    S_LT <- S(LT)
    S_UT <- S(UT)
    S_UC <- S(prep$UC)
    # An upper-censored response is retained as a response; an upper-censored
    # omission contributes only when UCresponse is FALSE.
    censor <- ifelse(!prep$UCresp & is.finite(prep$UC),
                     pmax(S_UC - S_UT, 0), 0)
    response_edge <- ifelse(!prep$UCresp & is.finite(prep$UC),
                            pmin(prep$UC, UT), UT)
    responded <- pmax(S_LT - S(response_edge), 0)

    sh <- if (length(mech)) .omx_shapley(mech, Vfun) else
      list(phi = list(), V = list(), V_none = S_inf, V_all = S_inf)
    comp <- list(contaminant = pC)
    if (length(sh$phi)) for (nm in names(sh$phi))
      comp[[nm]] <- sh$phi[[nm]]
    split_values <- .omx_split_inf(model_list, mech, pars,
                                    nrow(prep$trial), n_acc)
    if (is.null(split_values)) {
      comp[[residual_name]] <- sh$V_none
    } else {
      # The Shapley game allocates only the mechanisms registered in `mech`;
      # split the remaining intrinsic mass into dead launches and the
      # turnaround asymptote.  This keeps the two entries additive and avoids
      # attributing their nested relationship as an interaction.
      split_name <- attr(mech, "split")[["name"]]
      comp[[split_name]] <- split_values$dead_launch
      comp[[residual_name]] <- split_values$asymptotic_subthreshold
    }
    comp$censor_slow <- censor

    # The race process is truncated and normalised first; pContaminant is a
    # mixture weight among retained trials (see contaminant_mixture.h).
    process_Z <- pmax(S_LT - S_UT, 0) + S_inf
    process_Z[process_Z <= 0] <- NA_real_
    process_names <- setdiff(names(comp), "contaminant")
    comp[process_names] <- lapply(comp[process_names],
                                  function(x) x / process_Z)
    comp[process_names] <- lapply(comp[process_names],
                                  function(x) (1 - pC) * x)
    responded <- (1 - pC) * responded / process_Z

    df <- data.frame(draw = i, cell = cell, as.data.frame(comp),
                     responded = responded, check.names = FALSE)
    num <- setdiff(names(df), c("draw", "cell"))
    agg <- aggregate(df[num], by = list(cell = df$cell), FUN = mean,
                     na.rm = TRUE)
    agg$draw <- i
    draws[[i]] <- agg

    if (length(sh$V)) {
      cm <- do.call(cbind, lapply(sh$V, function(v) {
        tapply((1 - pC) * v / process_Z, cell, mean, na.rm = TRUE)
      }))
      cm <- as.matrix(cm)
      if (is.null(dimnames(cm))) {
        cm <- matrix(cm, nrow = length(levels(cell)),
                     dimnames = list(levels(cell), names(sh$V)))
      }
      coal_draws[[i]] <- cm
    }
  }

  draws <- do.call(rbind, draws)
  component_cols <- setdiff(names(draws), c("draw", "cell", "responded"))
  draws$omission_total <- rowSums(draws[, component_cols, drop = FALSE])

  value_cols <- setdiff(names(draws), c("draw", "cell"))
  long <- do.call(rbind, lapply(value_cols, function(v) data.frame(
    draw = draws$draw, cell = draws$cell, component = v, p = draws[[v]],
    check.names = FALSE
  )))
  long$component <- factor(long$component, levels = value_cols)
  groups <- split(long, list(long$cell, long$component), drop = TRUE)
  summary <- do.call(rbind, lapply(groups, function(d) {
    q <- quantile(d$p, probs = probs, na.rm = TRUE, names = FALSE)
    out <- data.frame(cell = d$cell[1L], component = d$component[1L],
                      mean = mean(d$p, na.rm = TRUE), check.names = FALSE)
    qnames <- paste0(probs * 100, "%")
    for (j in seq_along(q)) out[[qnames[j]]] <- q[j]
    out
  }))
  summary <- summary[order(summary$cell, summary$component), , drop = FALSE]

  obs_cell <- if (length(cell_cols)) {
    interaction(prep$data[, cell_cols, drop = FALSE], sep = " ", drop = TRUE)
  } else factor(rep("all", nrow(prep$data)))
  observed <- data.frame(
    cell = levels(obs_cell),
    omission_observed = as.numeric(tapply(is.na(prep$data$R), obs_cell, mean)),
    n = as.numeric(table(obs_cell))
  )

  coalitions <- if (length(coal_draws)) {
    good <- coal_draws[!vapply(coal_draws, is.null, logical(1L))]
    if (length(good)) {
      cm <- Reduce(`+`, good) / length(good)
      colnames(cm)[colnames(cm) == .omx_key(character())] <- "none"
      colnames(cm)[colnames(cm) == .omx_key(names(mech))] <- "all"
      cm
    } else NULL
  } else NULL

  structure(list(
    draws = draws, summary = summary, long = long, observed = observed,
    coalitions = coalitions, mechanisms = omission_mechanisms(emc, data),
    cell_cols = cell_cols, c_name = model_list$c_name
  ), class = "omission_decomposition")
}

#' @export
print.omission_decomposition <- function(x, ...) {
  cat("Omission decomposition for", x$c_name, "\n")
  cat("cells:", if (length(x$cell_cols)) paste(x$cell_cols, collapse = " x ")
      else "(pooled)", "\n\n")
  if (nrow(x$mechanisms)) print(x$mechanisms, row.names = FALSE)
  cat("\n")
  live <- if (nrow(x$mechanisms)) as.character(x$mechanisms$mechanism) else
    character()
  breakdown <- x$summary[
    as.character(x$summary$component) %in% c(live, "omission_total"), ,
    drop = FALSE
  ]
  print(breakdown, row.names = FALSE, digits = 3)
  if (!is.null(x$coalitions) && ncol(x$coalitions) > 2L) {
    cat("\nnever-finish mass with registered mechanisms switched on:\n")
    print(x$coalitions, digits = 3)
  }
  cat("\nobserved:\n")
  print(x$observed, row.names = FALSE, digits = 3)
  invisible(x)
}

#' Stacked bars for an omission decomposition.
#' @param x An `omission_decomposition` object.
#' @param use Summary column to plot, usually `"mean"`.
#' @param ... Additional arguments passed to `barplot`.
#' @return The matrix passed to `barplot`, invisibly.
#' @export
plot_omission_decomposition <- function(x, use = "mean", ...) {
  s <- x$summary[!x$summary$component %in%
                   c("responded", "omission_total"), , drop = FALSE]
  if (!nrow(s)) return(invisible(matrix(numeric(), 0, 0)))
  keep <- tapply(s[[use]], s$component, function(z) any(abs(z) > 1e-10))
  s <- s[as.character(s$component) %in% names(keep)[keep], , drop = FALSE]
  if (!nrow(s)) return(invisible(matrix(numeric(), 0, 0)))
  m <- tapply(s[[use]], list(s$component, s$cell), function(z) z[1L])
  m <- as.matrix(m)
  m[is.na(m)] <- 0
  cols <- grDevices::hcl.colors(nrow(m), "Set 2")
  op <- par(mar = c(5, 4, 4, 8) + .1)
  on.exit(par(op), add = TRUE)
  ymax <- max(c(colSums(m), x$observed$omission_observed), na.rm = TRUE)
  bp <- graphics::barplot(m, col = cols, ylab = "P(omission)", las = 2,
                ylim = c(0, if (ymax > 0) ymax * 1.15 else 1), ...)
  obs <- x$observed$omission_observed[match(colnames(m), x$observed$cell)]
  points(bp, obs, pch = 4, lwd = 2)
  legend(par("usr")[2], par("usr")[4], xpd = NA, bty = "n", cex = .8,
         legend = c(rownames(m), "observed"),
         fill = c(cols, NA), border = c(rep("black", nrow(m)), NA),
         pch = c(rep(NA, nrow(m)), 4))
  invisible(m)
}
