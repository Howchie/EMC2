#!/usr/bin/env Rscript

# Benchmark for the correlated BAwL/LBA likelihood routes
# (bawl_corr_exact_kernel_plan.md).  Fixed seeds; each scenario evaluates
# calc_ll_oo over EMC2_BAWLCORR_BENCH_N unique trials, reporting absolute
# time, the ratio to the ordinary (independent) LBA/BAwL likelihood, and the
# route/node counters from src/bawl_corr_counters.h.
#
# Run from the package root:
#   Rscript tools/benchmark-bawl-corr-exact.R
# Env:
#   EMC2_BAWLCORR_BENCH_N     unique trials per scenario (default 20000)
#   EMC2_BAWLCORR_BENCH_REPS  timing repetitions (default 3)
#   EMC2_BAWLCORR_BENCH_LIB   if set, library path holding an installed EMC2
#                             (otherwise devtools::load_all() of this tree)

suppressPackageStartupMessages({
  lib <- Sys.getenv("EMC2_BAWLCORR_BENCH_LIB", "")
  if (nzchar(lib)) {
    library(EMC2, lib.loc = lib)
    cat(sprintf("EMC2 from %s\n", dirname(system.file(package = "EMC2"))))
  } else {
    devtools::load_all(quiet = TRUE)
    cat("EMC2 via devtools::load_all()\n")
  }
})

n_unique <- as.integer(Sys.getenv("EMC2_BAWLCORR_BENCH_N", "20000"))
reps <- as.integer(Sys.getenv("EMC2_BAWLCORR_BENCH_REPS", "3"))
stopifnot(is.finite(n_unique), n_unique >= 10L, is.finite(reps), reps >= 1L)
Sys.setenv(EMC2_BAWLCORR_COUNTERS = "1")

set.seed(20260717)

# --- data builders ---------------------------------------------------------

two_racer_data <- function(n, near_t0 = FALSE) {
  rt <- if (near_t0) .1 + runif(n, .02, .08) else .2 + rlnorm(n, log(.45), .4)
  data.frame(
    subjects = factor(1),
    S = factor(rep("correct", n), levels = c("correct", "error")),
    R = factor(ifelse(runif(n) < .7, "correct", "error"),
               levels = c("correct", "error")),
    rt = rt
  )
}

three_racer_data <- function(n, race23 = TRUE) {
  race <- if (race23) sample(2:3, n, replace = TRUE) else rep(3L, n)
  resp <- ifelse(runif(n) < .7, "correct", "error")
  pm_win <- race == 3L & runif(n) < .12
  resp[pm_win] <- "pm"
  out <- data.frame(
    subjects = factor(1),
    S = factor(rep("correct", n), levels = c("correct", "error", "pm")),
    R = factor(resp, levels = c("correct", "error", "pm")),
    rt = .2 + rlnorm(n, log(.45), .4)
  )
  if (race23) out$RACE <- factor(race, levels = 2:3)
  out
}

# --- design/context builders (mirrors tests/testthat/test-bawl-correlated.R)

coupled_fun <- function(d) factor(d$lR != "pm", levels = c(FALSE, TRUE),
                                  labels = c("no", "yes"))

make_ctx <- function(dat, model, formula, functions = NULL, constants = NULL) {
  des <- design(
    data = dat,
    Rlevels = levels(dat$R),
    matchfun = function(d) as.character(d$S) == as.character(d$lR),
    functions = functions,
    formula = formula,
    constants = constants,
    model = model,
    report_p_vector = FALSE
  )
  emc <- make_emc(dat, des, type = "single", n_chains = 1,
                  compress = FALSE, verbose = FALSE, rt_resolution = NULL)
  list(emc = emc[[1]], design = des)
}

context_args <- function(ctx, p) {
  dadm <- ctx$emc$data[[1]]
  model <- ctx$emc$model()
  designs <- lapply(names(model$p_types), function(nm) {
    dm <- attr(dadm, "designs")[[nm]]
    dm[attr(dm, "expand"), , drop = FALSE]
  })
  names(designs) <- names(model$p_types)
  list(
    particle_matrix = matrix(p, nrow = 1, dimnames = list(NULL, names(p))),
    data = dadm,
    constants = attr(dadm, "constants"),
    designs = designs,
    type = model$c_name,
    bounds = model$bound,
    transforms = model$transform,
    pretransforms = model$pre_transform,
    p_types = names(model$p_types),
    min_ll = log(1e-10),
    trend = model$trend
  )
}

set_values <- function(p, rho = NULL, v = 1, sv = .8, B = 1.4, A = .2,
                       t0 = .1, k = 0, mG = 1, mK = 1) {
  natural <- c(v = v, sv = sv, B = B, A = A, t0 = t0, k = k, mG = mG, mK = mK)
  for (nm in names(natural)) {
    hit <- names(p) == nm
    p[hit] <- if (nm == "v") natural[[nm]] else log(natural[[nm]])
  }
  if (!is.null(rho)) {
    rho_cols <- grep("^rho", names(p), value = TRUE)
    p[rho_cols] <- qnorm((rho + 1) / 2)
  }
  p
}

# --- scenario table --------------------------------------------------------

base_formula <- list(v ~ 1, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1,
                     mG ~ 1, mK ~ 1)
rho1_formula <- c(base_formula, rho ~ 1)
rho_coupled_formula <- c(base_formula, rho ~ 0 + coupled)

scenarios <- list()
add_scenario <- function(name, build) {
  scenarios[[name]] <<- build
}

add_scenario("plain_lba_2", function(n) {
  ctx <- make_ctx(two_racer_data(n), BAwL(posdrift = TRUE), base_formula)
  p <- set_values(sampled_pars(ctx$design, doMap = FALSE))
  list(ctx = ctx, p = p)
})

for (rho in c(.5, .8, .95)) {
  local({
    r <- rho
    add_scenario(sprintf("corr_lba_2_rho%02d", round(100 * r)), function(n) {
      ctx <- make_ctx(two_racer_data(n), BAwLcorr(posdrift = TRUE), rho1_formula)
      p <- set_values(sampled_pars(ctx$design, doMap = FALSE), rho = r)
      list(ctx = ctx, p = p)
    })
  })
}

add_scenario("corr_lba_2_rho08_unrestricted", function(n) {
  ctx <- make_ctx(two_racer_data(n), BAwLcorr(posdrift = FALSE), rho1_formula)
  p <- set_values(sampled_pars(ctx$design, doMap = FALSE), rho = .8)
  list(ctx = ctx, p = p)
})

add_scenario("corr_lba_2_rho08_near_t0", function(n) {
  ctx <- make_ctx(two_racer_data(n, near_t0 = TRUE),
                  BAwLcorr(posdrift = TRUE), rho1_formula)
  # Small threshold/start so densities at rt - t0 in [.02, .08] stay above
  # the min_ll floor; the point is stressing the narrow-integrand regime,
  # not timing floored trials.
  p <- set_values(sampled_pars(ctx$design, doMap = FALSE), rho = .8,
                  v = 1.5, sv = 1, B = .25, A = .15)
  list(ctx = ctx, p = p)
})

add_scenario("leak_bawl_2_rho08", function(n) {
  ctx <- make_ctx(two_racer_data(n), BAwLcorr(posdrift = TRUE), rho1_formula)
  p <- set_values(sampled_pars(ctx$design, doMap = FALSE), rho = .8, k = .2)
  list(ctx = ctx, p = p)
})

add_scenario("plain_bawl_pm_race23", function(n) {
  ctx <- make_ctx(three_racer_data(n), BAwL(posdrift = TRUE), base_formula)
  p <- set_values(sampled_pars(ctx$design, doMap = FALSE), k = .2)
  list(ctx = ctx, p = p)
})

add_scenario("corr_lba_pm_race23_rho08", function(n) {
  ctx <- make_ctx(three_racer_data(n), BAwLcorr(posdrift = TRUE),
                  rho_coupled_formula,
                  functions = list(coupled = coupled_fun),
                  constants = c(rho_coupledno = 0))
  p <- set_values(sampled_pars(ctx$design, doMap = FALSE), rho = .8)
  list(ctx = ctx, p = p)
})

add_scenario("leak_bawl_pm_race23_rho08", function(n) {
  ctx <- make_ctx(three_racer_data(n), BAwLcorr(posdrift = TRUE),
                  rho_coupled_formula,
                  functions = list(coupled = coupled_fun),
                  constants = c(rho_coupledno = 0))
  p <- set_values(sampled_pars(ctx$design, doMap = FALSE), rho = .8, k = .2)
  list(ctx = ctx, p = p)
})

add_scenario("forced_gh_3loaded_rho08", function(n) {
  ctx <- make_ctx(three_racer_data(n, race23 = FALSE),
                  BAwLcorr(posdrift = TRUE), rho1_formula)
  p <- set_values(sampled_pars(ctx$design, doMap = FALSE), rho = .8)
  list(ctx = ctx, p = p)
})

add_scenario("forced_generic_clock_rho08", function(n) {
  ctx <- make_ctx(two_racer_data(n),
                  BAwLcorr(posdrift = TRUE, erlang_type = "local_kill_guess"),
                  rho1_formula)
  p <- set_values(sampled_pars(ctx$design, doMap = FALSE), rho = .8,
                  mG = .8, mK = 1.2)
  list(ctx = ctx, p = p)
})

# --- run -------------------------------------------------------------------

counter_names <- names(EMC2:::bawl_corr_counter_values())
rows <- list()
counter_rows <- list()

for (name in names(scenarios)) {
  built <- scenarios[[name]](n_unique)
  args <- context_args(built$ctx, built$p)
  # Warm-up (also validates finiteness) and per-call counters.
  EMC2:::bawl_corr_counters_reset()
  ll <- as.numeric(do.call(EMC2:::calc_ll_oo, args))
  counters <- unlist(EMC2:::bawl_corr_counter_values())
  if (!is.finite(ll)) warning(sprintf("scenario %s: non-finite ll", name))
  times <- vapply(seq_len(reps), function(i)
    system.time(do.call(EMC2:::calc_ll_oo, args))[["elapsed"]], numeric(1))
  rows[[name]] <- data.frame(
    scenario = name, ll = ll, median_seconds = median(times),
    stringsAsFactors = FALSE
  )
  counter_rows[[name]] <- counters
  rm(built, args)
  gc(verbose = FALSE)
}

result <- do.call(rbind, rows)
ref_time <- result$median_seconds[result$scenario == "plain_lba_2"]
result$ratio_to_plain_lba <- result$median_seconds / ref_time
print(result, row.names = FALSE, digits = 4)

cat("\nRoute/node counters per single calc_ll_oo call:\n")
counters <- do.call(rbind, counter_rows)
nonzero <- counters[, colSums(counters) > 0, drop = FALSE]
print(nonzero)
cat(sprintf("\nn_unique_trials = %d, reps = %d\n", n_unique, reps))
