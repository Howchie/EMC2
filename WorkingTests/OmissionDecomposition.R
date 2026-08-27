## ---------------------------------------------------------------------------
## decompose_omissions(): mechanism-wise decomposition of predicted omissions
##
## Every EMC2 race model produces omissions (R = NA) through at most three
## structurally different routes:
##
##   1. contaminant   pContaminant: a Bernoulli omission applied before anything
##                    else, i.e. mass at rt = Inf that the process never saw.
##   2. never_finish  no accumulator EVER crosses threshold.  This is the
##                    model-dependent part and is split further (see below).
##   3. censor_slow   an accumulator does finish, but after UC ("design
##                    omissions").  Only an omission when UCresponse is FALSE.
##
## Lower censoring is NOT handled: these designs have LT but no LC, so fast
## responses are truncated away rather than turned into omissions.  A finite
## LC > LT would create a fourth route and is refused rather than ignored.
##
## Everything is computed exactly from the model's own pfun, so it agrees with
## the likelihood by construction; nothing is simulated.  Writing
##
##     S(t) = prod_r (1 - F_r(t))          (no accumulator has finished by t)
##
## and pC for pContaminant, then, per trial,
##
##     P(retained)    = pC + (1 - pC) * S(LT)      [rt < LT is truncated away]
##     never_finish   = (1 - pC) * S(Inf)
##     censor_slow    = (1 - pC) * (S(UC) - S(Inf))
##     responded      = (1 - pC) * (S(LT) - S(UC))
##
## all divided by P(retained), which is what makes them comparable to the
## omission rate you measure in the fitted data.  S(Inf) < 1 is exactly the
## defective mass the EMC2 pfuns return at rt = Inf.
##
## The never_finish mass is split into named mechanisms by ABLATION: each
## mechanism is a switch that can be turned off (k = 0 removes BAwL leak,
## lambda_k = 0 removes the kill clock, posdrift = TRUE removes negative drift
## draws, ell = 0 removes BAwD clearance, ...).  Re-evaluating S(Inf) with any
## subset of the switches on gives a coalitional value function, and the
## mechanisms are credited their Shapley values, which (a) sum exactly to the
## total never-finish mass and (b) are order independent, which matters because
## the mechanisms interact (leak and a kill clock both eat the same trials).
## Whatever survives with EVERY switch off is reported under the residual name,
## "unexplained" unless the family gives it a meaning -- that is also what an
## unregistered model gets, so the totals are always right even if the split is
## not available.
##
## CAVEAT, and read $coalitions before interpreting a split: Shapley divides
## interaction between mechanisms evenly.  When no mechanism produces much
## never-finish mass ON ITS OWN, the reported split is that convention rather
## than two separable causes, and the function warns.  Where mechanisms nest
## (ablating one gives a strict SUBSET of the failures) the split is exact and
## needs no convention; that is why BAwD registers one mechanism plus a named
## residual rather than "decay" and "clearance" as a pair.
##
## A mechanism is only reported if it is actually live in the fitted parameters
## (k > 0 somewhere, mK finite, posdrift = FALSE, ...), which is the
## "recognises which components are valid" part.
##
## Usage
##   source("OmissionDecomposition.R")
##   d <- decompose_omissions(emc, factors = c("L", "S"), n_post = 100)
##   d$summary                       # median + 95% CI per cell per component
##   plot_omission_decomposition(d)  # stacked bars, observed rate overlaid
##   omission_mechanisms(emc)        # just the mechanism audit
## ---------------------------------------------------------------------------

## --- model introspection ----------------------------------------------------

# TRUE when the model's own drift draws are truncated positive.  The IO
# ("intrinsic omission") suffix is EMC2's marker for the untruncated variant.
.omx_posdrift <- function(c_name) !grepl("IO", c_name)

.omx_erlang <- function(c_name) {
  if (grepl("EMIX", c_name)) 3L else if (grepl("_E2", c_name)) 2L else 1L
}

.omx_family <- function(c_name) {
  if (grepl("BAwD", c_name)) "BAwD"
  else if (grepl("BAwL", c_name)) "BAwL"
  else if (grepl("RDMSWTN|RDMGBM", c_name)) "RDMSWTN"
  else if (grepl("^RDM", c_name)) "RDM"
  else if (grepl("LBA", c_name)) "LBA"
  else NA_character_
}

# A pfun for this model with posdrift forced to `posdrift`.  Returns NULL when
# the family has no such switch, in which case negative drift is not a
# separable mechanism and the mass lands in "unexplained".
.omx_pfun <- function(model_list, posdrift) {
  c_name <- model_list$c_name
  if (identical(posdrift, .omx_posdrift(c_name))) return(model_list$pfun)
  erl <- .omx_erlang(c_name)
  guess <- grepl("GUESS", c_name)
  switch(.omx_family(c_name),
    BAwL = function(rt, pars) EMC2:::pBAwL(rt, pars, posdrift = posdrift,
                                           erlang = erl, guess = guess),
    BAwD = function(rt, pars) EMC2:::pBAwD(
      rt, pars, launch = if (grepl("LOGN", c_name)) 1L else 0L,
      posdrift = posdrift),
    RDMSWTN = function(rt, pars) EMC2:::pRDMSWTN(rt, pars, erlang = erl,
                                                 posdrift = posdrift),
    LBA = function(rt, pars) EMC2:::.lba_pfun(rt, pars, posdrift = posdrift),
    NULL)
}

.omx_setcol <- function(cols, value) {
  force(cols); force(value)
  function(pars) {
    hit <- intersect(cols, colnames(pars))
    if (length(hit)) pars[, hit] <- value
    pars
  }
}

# The mechanism registry.  Each entry: name, a description, `null` (how to turn
# the mechanism OFF in the post-Ttransform parameter matrix) and `posdrift`
# (TRUE when turning it off means calling pfun with posdrift = TRUE instead).
# `pars` is a mapped parameter matrix, used only to decide what is live.
#
# The mass that survives with every switch off is reported under the name in
# attr(, "residual").  For most families that is genuinely "unexplained", but
# where it has a meaning it gets named -- see the BAwD note below.
.omx_mechanisms <- function(model_list, pars) {
  c_name <- model_list$c_name
  fam <- .omx_family(c_name)
  live <- function(col, f = function(x) any(x > 0, na.rm = TRUE)) {
    col %in% colnames(pars) && f(pars[, col])
  }
  m <- list()
  residual <- c(name = "unexplained",
                desc = "never-finish mass left with every mechanism switched off")
  add <- function(name, desc, null = identity, posdrift = FALSE)
    m[[name]] <<- list(name = name, desc = desc, null = null, posdrift = posdrift)

  # Generic across every race model that carries a kill clock.
  if (live("lambda_k"))
    add("kill", "kill clock fires before threshold (mK)",
        null = .omx_setcol("lambda_k", 0))

  if (identical(fam, "BAwL") && live("k"))
    add("leak", "leak asymptote v/k below threshold",
        null = .omx_setcol("k", 0))

  # BAwD.  Decay and clearance are NOT two separable causes: an omission needs
  # the peak of z + (V/k)(1 - e^{-kt}) - ell*t to fall short of b, which with
  # k > 0 and ell > 0 is one event, and ablating either switch alone leaves
  # (empirically) almost no mass -- a Shapley split would be reporting a pure
  # interaction as if it were two main effects.  What IS separable is nested and
  # convention-free: X(t) is decreasing in k pointwise, so the trials that fail
  # even at k = 0 (V <= ell: the accumulator never rises at all) are a SUBSET of
  # the failures.  So split total = dead_launch + drive_decay exactly.
  # Fixing ell (the usual scale convention for the lognormal launch) does not
  # bias this: both V <= ell and V/k <= b - z are dimensionless, so the split is
  # invariant under (V, ell, b, A) -> c*(V, ell, b, A), verified to 1e-15.
  if (identical(fam, "BAwD")) {
    if (live("k"))
      add("drive_decay", "rose, but the decaying drive peaked below threshold",
          null = .omx_setcol("k", 0))
    if (live("ell"))
      residual <- c(name = "dead_launch",
                    desc = "launch strength below clearance ell: never rises")
  }

  # Untruncated drift draws: only separable when the family has a posdrift
  # switch we can flip.
  if (!.omx_posdrift(c_name) && !is.null(.omx_pfun(model_list, TRUE)))
    add("negative_drift", "non-positive drift draw (posdrift = FALSE)",
        posdrift = TRUE)

  attr(m, "residual") <- residual
  m
}

## --- exported: mechanism audit ---------------------------------------------

#' Which omission mechanisms are live in a fitted emc object.
omission_mechanisms <- function(emc, data = NULL) {
  prep <- .omx_prep(emc, data)
  mech <- .omx_mechanisms(prep$model_list, prep$pars_ref)
  rows <- list(data.frame(
    mechanism = "contaminant", live = any(prep$pars_ref_pC > 0),
    description = "pContaminant: Bernoulli omission before the race"))
  for (mm in mech) rows[[length(rows) + 1L]] <-
    data.frame(mechanism = mm$name, live = TRUE, description = mm$desc)
  # The residual is live only if switching everything off still leaves mass.
  res <- attr(mech, "residual")
  rows[[length(rows) + 1L]] <- data.frame(
    mechanism = res[["name"]],
    live = any(.omx_ablate_all(prep, mech) > 1e-10),
    description = res[["desc"]])
  rows[[length(rows) + 1L]] <- data.frame(
    mechanism = "censor_slow", live = any(is.finite(prep$UC)),
    description = "response after UC (design omission)")
  out <- do.call(rbind, rows)
  attr(out, "c_name") <- prep$model_list$c_name
  out
}

## --- internals: design + parameter mapping ----------------------------------

.omx_prep <- function(emc, data = NULL) {
  design <- get_design(emc)
  if (length(design) > 1L)
    stop("decompose_omissions() handles single-design (non-joint) fits")
  design <- design[[1]]
  model <- design$model
  model_list <- model()
  if (!identical(model_list$type, "RACE"))
    stop("decompose_omissions() is defined for race models; got type '",
         model_list$type, "'")
  if (is.null(data)) data <- get_data(emc)
  if (is.list(data) && !is.data.frame(data)) data <- data[[1]]
  data <- data[order(data$subjects), ]

  dadm <- EMC2:::design_model(
    EMC2:::add_accumulators(
      EMC2:::add_trials(data), design$matchfun, simulate = TRUE,
      type = model_list$type, Fcovariates = design$Fcovariates,
      fixed_accumulator_roles = design$fixed_accumulator_roles),
    design, model, add_acc = FALSE, compress = FALSE, verbose = FALSE,
    rt_check = FALSE)

  n_acc <- length(levels(dadm$lR))
  # Everything below reshapes the parameter matrix with matrix(nrow = n_acc),
  # which is only valid because dadm cycles lR fastest within trial.
  if (!identical(as.integer(dadm$lR), rep(seq_len(n_acc), nrow(dadm) / n_acc)))
    stop("unexpected accumulator ordering in the expanded design")
  first <- dadm$lR == levels(dadm$lR)[1]
  trial <- dadm[first, , drop = FALSE]

  # Censoring/truncation bounds: trial columns win, then the design's TC, then
  # the no-censoring defaults make_missing() would apply.
  tc <- attr(design, "TC")
  getcol <- function(nm, default) {
    if (!is.null(dadm[[nm]])) return(as.numeric(trial[[nm]]))
    if (!is.null(tc[[nm]]) && length(tc[[nm]]) == 1L)
      return(rep(as.numeric(tc[[nm]]), nrow(trial)))
    rep(default, nrow(trial))
  }
  UC <- getcol("UC", Inf); LC <- getcol("LC", 0)
  LT <- getcol("LT", 0);   UT <- getcol("UT", Inf)
  # Upper censoring only makes an omission when the response is discarded with
  # the RT.
  UCresp <- if (!is.null(dadm$UCresponse)) as.logical(trial$UCresponse) else
    rep(isTRUE(tc$UCresponse), nrow(trial))
  # Lower censoring is out of scope (see the header): with LC <= LT a fast
  # response is truncated away, never recorded as an omission.
  if (any(LC > LT, na.rm = TRUE))
    stop("lower censoring (LC > LT) present: decompose_omissions() only ",
         "handles LT, so fast-censored omissions would be missed")

  # A reference parameter matrix (posterior means) used only to decide which
  # mechanisms are live.
  ref <- .omx_alpha(emc, n_post = 1, stat = "mean")[[1]]
  pars_ref <- .omx_map(ref, dadm, model_list)
  pC <- if ("pContaminant" %in% colnames(pars_ref))
    pars_ref[first, "pContaminant"] else rep(0, nrow(trial))

  list(design = design, model = model, model_list = model_list, data = data,
       dadm = dadm, trial = trial, n_acc = n_acc, first = first,
       UC = UC, LT = LT, UT = UT, UCresp = UCresp,
       pars_ref = pars_ref, pars_ref_pC = pC)
}

.omx_alpha <- function(emc, n_post = 50, stat = "random") {
  samps <- get_pars(emc, selection = "alpha", merge_chains = TRUE,
                    by_subject = TRUE)
  if (stat != "random") {
    p <- do.call(rbind, lapply(samps, function(x) apply(x[[1]], 2, stat)))
    return(list(p))
  }
  lapply(seq_len(n_post), function(i)
    do.call(rbind, lapply(samps, function(x)
      x[[1]][sample(seq_len(nrow(x[[1]])), 1), ])))
}

# alpha -> per-accumulator natural-scale parameters, exactly as make_data() does.
.omx_map <- function(alpha, dadm, model_list) {
  pars <- EMC2:::get_pars_oo(alpha, dadm, function() model_list)
  pars <- model_list$Ttransform(pars, dadm)
  EMC2:::fix_bound(pars, model_list$bound, dadm$lR)
}

## --- internals: survivor + Shapley -----------------------------------------

# S(t) = P(no accumulator has finished by t), one value per trial.
.omx_S <- function(t_trial, pars, pfun, n_acc) {
  n_trial <- length(t_trial)
  rt <- rep(t_trial, each = n_acc)
  p <- pfun(rt, pars)
  p[!is.finite(p)] <- 0
  surv <- matrix(1 - pmin(pmax(p, 0), 1), nrow = n_acc)
  pmin(pmax(apply(surv, 2, prod), 0), 1)
}

# Value function for the ablation game: never-finish mass with the mechanisms
# in `set` switched ON and every other registered mechanism switched OFF.
.omx_Vfun <- function(model_list, mech, pars, n_trial, n_acc) {
  function(set) {
    p <- pars; f <- model_list$pfun
    for (nm in setdiff(names(mech), set)) {
      mm <- mech[[nm]]
      p <- mm$null(p)
      if (isTRUE(mm$posdrift)) {
        f2 <- .omx_pfun(model_list, TRUE)
        if (!is.null(f2)) f <- f2
      }
    }
    .omx_S(rep(Inf, n_trial), p, f, n_acc)
  }
}

# Never-finish mass with every mechanism off, at the reference parameters.
.omx_ablate_all <- function(prep, mech) {
  .omx_Vfun(prep$model_list, mech, prep$pars_ref, nrow(prep$trial),
            prep$n_acc)(character(0))
}

# "" is not usable as a list name in R, so the empty coalition gets a sentinel.
.omx_key <- function(s) if (!length(s)) "<none>" else paste(sort(s), collapse = "|")

# Shapley split of V(all) - V(none) over the mechanisms in `mech`.
.omx_shapley <- function(mech, Vfun) {
  nms <- names(mech)
  M <- length(nms)
  subsets <- lapply(0:(2^M - 1), function(z)
    nms[bitwAnd(z, 2L^(seq_len(M) - 1L)) > 0])
  V <- list()
  for (s in subsets) V[[.omx_key(s)]] <- Vfun(s)
  phi <- setNames(vector("list", M), nms)
  for (m in nms) {
    others <- setdiff(nms, m)
    acc <- 0
    for (z in 0:(2^length(others) - 1)) {
      T_ <- others[bitwAnd(z, 2L^(seq_along(others) - 1L)) > 0]
      w <- factorial(length(T_)) * factorial(M - length(T_) - 1) / factorial(M)
      acc <- acc + w * (V[[.omx_key(c(T_, m))]] - V[[.omx_key(T_)]])
    }
    phi[[m]] <- acc
  }
  # V is returned as well: when a mechanism ON ITS OWN explains almost none of
  # V_all, the Shapley numbers are splitting an interaction rather than
  # reporting two causes, and should be read as a convention.  See $coalitions.
  list(phi = phi, V = V, V_none = V[[.omx_key(character(0))]],
       V_all = V[[.omx_key(nms)]])
}

## --- exported: the decomposition -------------------------------------------

#' Decompose predicted omissions by generating mechanism.
#'
#' @param emc a fitted emc object (single design, race model).
#' @param factors character vector of data columns to split by.  NULL pools.
#' @param n_post number of posterior draws (ignored unless stat = "random").
#' @param stat "random" (default), "mean" or "median" over the posterior.
#' @param by_subject also split by subject.
#' @param data optional data to evaluate on; defaults to the fitted data.
#' @param probs credible interval for the summary.
#' @return list(draws=, summary=, mechanisms=, observed=)
decompose_omissions <- function(emc, factors = NULL, n_post = 50,
                                stat = c("random", "mean", "median"),
                                by_subject = FALSE, data = NULL,
                                probs = c(.025, .5, .975)) {
  stat <- match.arg(stat)
  prep <- .omx_prep(emc, data)
  model_list <- prep$model_list
  dadm <- prep$dadm; n_acc <- prep$n_acc; first <- prep$first
  mech <- .omx_mechanisms(model_list, prep$pars_ref)
  residual_name <- attr(mech, "residual")[["name"]]

  if (!is.null(factors)) {
    miss <- setdiff(factors, names(prep$trial))
    if (length(miss)) stop("factors not in data: ", paste(miss, collapse = ", "))
  }
  cell_cols <- c(if (by_subject) "subjects", factors)
  cell <- if (length(cell_cols))
    interaction(prep$trial[, cell_cols, drop = FALSE], sep = " ", drop = TRUE) else
      factor(rep("all", nrow(prep$trial)))

  alphas <- .omx_alpha(emc, n_post = n_post, stat = stat)
  if (stat != "random") n_post <- 1L

  UC <- prep$UC; LT <- prep$LT
  # Censoring that keeps the response is not an omission.
  UC_om <- ifelse(prep$UCresp, Inf, UC)

  out <- vector("list", n_post)
  coal <- vector("list", n_post)
  for (i in seq_len(n_post)) {
    pars <- .omx_map(alphas[[i]], dadm, model_list)
    ok <- attr(pars, "ok")
    if (!is.null(ok) && mean(!ok) > 0) warning(
      "draw ", i, ": ", round(100 * mean(!ok), 2), "% of rows out of bounds")
    pC <- if ("pContaminant" %in% colnames(pars)) pars[first, "pContaminant"] else
      rep(0, nrow(prep$trial))

    S <- function(t, p = pars, f = model_list$pfun) .omx_S(t, p, f, n_acc)
    Vfun <- .omx_Vfun(model_list, mech, pars, nrow(prep$trial), n_acc)

    S_inf <- S(rep(Inf, nrow(prep$trial)))
    S_UC  <- S(UC_om)
    S_LT  <- S(LT)

    sh <- if (length(mech)) .omx_shapley(mech, Vfun) else
      list(phi = list(), V = list(), V_none = S_inf, V_all = S_inf)

    comp <- list(contaminant = pC)
    for (nm in names(sh$phi)) comp[[nm]] <- (1 - pC) * sh$phi[[nm]]
    comp[[residual_name]] <- (1 - pC) * sh$V_none
    comp$censor_slow <- (1 - pC) * pmax(S_UC - S_inf, 0)
    responded <- (1 - pC) * pmax(S_LT - S_UC, 0)

    # Condition on being retained: trials with rt < LT are truncated away.
    retained <- pC + (1 - pC) * S_LT
    retained[retained <= 0] <- NA_real_
    comp <- lapply(comp, function(x) x / retained)
    responded <- responded / retained

    df <- data.frame(draw = i, cell = cell,
                     as.data.frame(comp), responded = responded,
                     check.names = FALSE)
    num <- setdiff(names(df), c("draw", "cell"))
    agg <- aggregate(df[num], by = list(cell = df$cell), FUN = mean,
                     na.rm = TRUE)
    agg$draw <- i
    out[[i]] <- agg

    # Raw coalition values, so a split that is really an interaction is visible.
    if (length(sh$V)) {
      cm <- sapply(sh$V, function(v)
        tapply((1 - pC) * v / retained, cell, mean, na.rm = TRUE))
      if (is.null(dim(cm)))  # single cell: sapply drops to a vector
        cm <- matrix(cm, nrow = 1L, dimnames = list(levels(cell), names(sh$V)))
      coal[[i]] <- cm
    }
  }
  draws <- do.call(rbind, out)
  draws$omission_total <- rowSums(
    draws[, setdiff(names(draws), c("draw", "cell", "responded")), drop = FALSE])

  # Long form + summary.
  value_cols <- setdiff(names(draws), c("draw", "cell"))
  long <- do.call(rbind, lapply(value_cols, function(v)
    data.frame(draw = draws$draw, cell = draws$cell, component = v,
               p = draws[[v]])))
  long$component <- factor(long$component, levels = value_cols)
  summ <- do.call(rbind, lapply(split(long, list(long$cell, long$component),
                                      drop = TRUE), function(d) {
    q <- quantile(d$p, probs = probs, na.rm = TRUE)
    data.frame(cell = d$cell[1], component = d$component[1], mean = mean(d$p),
               t(q), check.names = FALSE, row.names = NULL)
  }))
  summ <- summ[order(summ$cell, summ$component), ]

  # Observed omission rate per cell, for calibration.
  obs_cell <- if (length(cell_cols))
    interaction(prep$data[, cell_cols, drop = FALSE], sep = " ", drop = TRUE) else
      factor(rep("all", nrow(prep$data)))
  observed <- data.frame(
    cell = levels(obs_cell),
    omission_observed = as.numeric(tapply(is.na(prep$data$R), obs_cell, mean)),
    n = as.numeric(table(obs_cell)))

  # Coalition values averaged over draws: the never-finish mass with only the
  # named mechanisms switched on.  If no single mechanism gets near the "all"
  # column, the Shapley split above is dividing an interaction, not two causes.
  coal <- if (length(mech)) {
    cm <- Reduce(`+`, coal) / length(coal)
    colnames(cm)[colnames(cm) == .omx_key(character(0))] <- "none"
    colnames(cm)[colnames(cm) == .omx_key(names(mech))] <- "all"
    cm
  } else NULL
  if (!is.null(coal) && length(mech) > 1L) {
    solo <- rowSums(coal[, names(mech), drop = FALSE] - coal[, "none"])
    if (any(coal[, "all"] - coal[, "none"] > 1e-8 &
            solo < 0.5 * (coal[, "all"] - coal[, "none"])))
      warning("mechanisms interact: no single mechanism reproduces much of the ",
              "never-finish mass, so the split is a Shapley convention -- see ",
              "$coalitions")
  }

  structure(list(draws = draws, summary = summ, long = long,
                 observed = observed, coalitions = coal,
                 mechanisms = omission_mechanisms(emc, data),
                 cell_cols = cell_cols, c_name = model_list$c_name),
            class = "omission_decomposition")
}

print.omission_decomposition <- function(x, ...) {
  cat("Omission decomposition for", x$c_name, "\n")
  cat("cells:", if (length(x$cell_cols)) paste(x$cell_cols, collapse = " x ") else
    "(pooled)", "\n\n")
  print(x$mechanisms, row.names = FALSE)
  cat("\n")
  s <- x$summary[x$summary$component != "responded", ]
  print(s, row.names = FALSE, digits = 3)
  if (!is.null(x$coalitions) && ncol(x$coalitions) > 2L) {
    cat("\nnever-finish mass with only these mechanisms switched on:\n")
    print(x$coalitions, digits = 3)
  }
  cat("\nobserved:\n")
  print(x$observed, row.names = FALSE, digits = 3)
  invisible(x)
}

#' Stacked bars of the omission decomposition, observed rate overlaid.
plot_omission_decomposition <- function(x, use = "mean", ...) {
  s <- x$summary[!x$summary$component %in% c("responded", "omission_total"), ]
  s <- s[!tapply(s[[use]], s$component, function(z) all(abs(z) < 1e-10))[
    as.character(s$component)], ]
  s$component <- droplevels(s$component)
  m <- tapply(s[[use]], list(s$component, s$cell), function(z) z[1])
  m[is.na(m)] <- 0
  cols <- grDevices::hcl.colors(nrow(m), "Set 2")
  op <- par(mar = c(5, 4, 4, 8) + .1); on.exit(par(op))
  bp <- barplot(m, col = cols, ylab = "P(omission)", las = 2,
                ylim = c(0, max(colSums(m), x$observed$omission_observed,
                                na.rm = TRUE) * 1.15), ...)
  obs <- x$observed$omission_observed[match(colnames(m), x$observed$cell)]
  points(bp, obs, pch = 4, lwd = 2)
  legend(par("usr")[2], par("usr")[4], xpd = NA, bty = "n", cex = .8,
         legend = c(rownames(m), "observed"),
         fill = c(cols, NA), border = c(rep("black", nrow(m)), NA),
         pch = c(rep(NA, nrow(m)), 4))
  invisible(m)
}
