# ---------------------------------------------------------------------------
# Where does a likelihood call spend its time?
#
#   Rscript WorkingTests/bench_likelihood_prologue.R
#   EMC_LIB=/path/to/lib Rscript WorkingTests/bench_likelihood_prologue.R
#
# Reports three things per benchmark:
#   * the whole call            EMC2:::calc_ll_manager()
#   * the parameter prologue    EMC2:::pt_prologue_oo()  -- refill from
#     the particle row, map the designs, transform, check bounds; i.e. every
#     step calc_ll_oo takes before it reaches a model kernel
#   * the remainder, which is kernel arithmetic and summation
#
# The prologue costs what it costs because parameter columns carry only
# DESIGN-CELL information while the prologue runs at TRIAL resolution. Watch
# the "prologue" share: intercept-only parameters should be nearly free (they
# take a compute-once + std::fill path in ParamTable/transform_utils), so a
# share creeping back above ~10% means that path has been lost.
#
# TRAPS, both of which make this benchmark measure the wrong thing:
#   * Do NOT benchmark at sampled_pars()' zeros. t0 = exp(0) = 1 s exceeds
#     every RT, every row floors at min_ll and the call is ~3x too fast for a
#     reason that has nothing to do with the code under test.
#   * library(EMC2) silently picks up a stale system install. The resolved
#     namespace path is printed below -- check it.
# ---------------------------------------------------------------------------

lib <- Sys.getenv("EMC_LIB")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(EMC2))
cat("EMC2 loaded from:", dirname(getNamespaceInfo("EMC2", "path")), "\n")

NT   <- as.integer(Sys.getenv("NT", "2000"))
NP   <- as.integer(Sys.getenv("NP", "100"))
REPS <- as.integer(Sys.getenv("REPS", "20"))

cpu <- function(expr, reps = REPS) {
  expr <- substitute(expr); env <- parent.frame()
  eval(expr, env)                                     # warm up
  t <- system.time(for (i in seq_len(reps)) eval(expr, env))
  t[["user.self"]] / reps
}

# --- fixtures --------------------------------------------------------------
make_fixture <- function(n_trials, censor_frac = 0) {
  set.seed(1)
  d <- forstmann[forstmann$subjects == levels(forstmann$subjects)[1L], ]
  d <- d[rep(seq_len(nrow(d)), length.out = n_trials), , drop = FALSE]
  d$subjects <- factor("s1")
  d$rt <- d$rt + runif(nrow(d), 0, 0.3)   # defeat rt compression
  if (censor_frac > 0) {
    cens <- sample(nrow(d), floor(censor_frac * nrow(d)))
    d$rt[cens] <- Inf
    d$R[cens] <- NA
    d$UC <- 2.0
  }
  rownames(d) <- NULL
  des <- design(data = d, model = LBA, matchfun = function(x) x$S == x$lR,
                formula = list(v ~ lM, B ~ E, A ~ 1, t0 ~ 1),
                constants = c(sv = log(1)))
  emc <- make_emc(d, des, n_chains = 1, compress = TRUE,
                  rt_resolution = NULL, type = "single")
  dadm <- emc[[1]]$data[[1]]

  # Realistic parameter values -- see the trap note above.
  centre <- sampled_pars(des); centre[] <- 0
  centre["v"] <- 2; centre["v_lMTRUE"] <- 1
  centre["B"] <- log(1); centre["A"] <- log(0.3); centre["t0"] <- log(0.15)
  if ("B_Eneutral"  %in% names(centre)) centre["B_Eneutral"]  <- 0.1
  if ("B_Eaccuracy" %in% names(centre)) centre["B_Eaccuracy"] <- 0.2
  prop <- matrix(rnorm(NP * length(centre), 0, 0.05), NP, length(centre),
                 dimnames = list(NULL, names(centre)))
  prop <- sweep(prop, 2, as.numeric(centre), "+")
  list(dadm = dadm, model = des$model, prop = prop, des = des)
}

report <- function(label, fx) {
  dadm  <- fx$dadm
  model <- fx$model()
  p_types <- names(model$p_types)
  designs <- lapply(p_types, function(p) attr(dadm, "designs")[[p]])
  names(designs) <- p_types
  constants <- attr(dadm, "constants"); if (is.null(constants)) constants <- NA

  t_all <- cpu(EMC2:::calc_ll_manager(fx$prop, dadm, fx$model, r_cores = 1))
  t_pro <- cpu(EMC2:::pt_prologue_oo(
    fx$prop, dadm, constants = constants, designs = designs,
    bounds = model$bound, transforms = model$transform,
    pretransforms = model$pre_transform, p_types = p_types,
    trend = model$trend))

  cat(sprintf("\n%s  (%d dadm rows, %d particles, all_finite=%s)\n",
              label, nrow(dadm), nrow(fx$prop),
              isTRUE(attr(dadm, "emc2_all_finite_trials"))))
  cat(sprintf("  whole call     %8.2f ms\n", 1000 * t_all))
  cat(sprintf("  prologue       %8.2f ms   %5.1f %%\n",
              1000 * t_pro, 100 * t_pro / t_all))
  cat(sprintf("  kernel + sum   %8.2f ms   %5.1f %%\n",
              1000 * (t_all - t_pro), 100 * (t_all - t_pro) / t_all))
  cat(sprintf("  first lls      %s\n",
              paste(sprintf("%.4f",
                    head(EMC2:::calc_ll_manager(fx$prop, dadm, fx$model, r_cores = 1), 3)),
                    collapse = "  ")))
  invisible(t_all)
}

report("all-finite RTs", make_fixture(NT))
report("15% omissions, UC = 2.0 (mixed path)", make_fixture(NT, censor_frac = 0.15))

# --- fixed vs per-particle cost --------------------------------------------
# make_pt_mapper (pretransform, ParamTable construction, design-plan build) is
# paid once per call. If it ever becomes a large share, caching it across calls
# would be worth its staleness risk; at the time of writing it is ~1%.
fx <- make_fixture(NT)
cat("\ncost vs particle count (all-finite):\n")
for (np in c(1L, 2L, 5L, 10L, 25L, NP)) {
  p <- fx$prop[seq_len(np), , drop = FALSE]
  t <- cpu(EMC2:::calc_ll_manager(p, fx$dadm, fx$model, r_cores = 1),
           reps = max(3L, as.integer(REPS * NP / np / 4)))
  cat(sprintf("  %4d particles  %8.2f ms  (%.3f ms/particle)\n",
              np, 1000 * t, 1000 * t / np))
}
