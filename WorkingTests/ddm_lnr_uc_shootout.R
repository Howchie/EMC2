#!/usr/bin/env Rscript

usage <- function(status = 0) {
  cat(
    "Usage:\n",
    "  Rscript WorkingTests/ddm_lnr_uc_shootout.R --lib LIB [--label NAME] [--lib LIB2 --label NAME2] [options]\n\n",
    "Options:\n",
    "  --lib PATH         Installed-library path containing EMC2. Repeat for shootouts.\n",
    "  --label NAME       Label for the preceding --lib. Defaults to basename(PATH).\n",
    "  --particles N      Number of particles per likelihood call. Default: 10000.\n",
    "  --reps N           Timing repetitions per workload. Default: 7.\n",
    "  --out PATH         Write combined CSV to PATH. Default: stdout.\n",
    "  --worker           Internal mode used to isolate compiled package builds.\n",
    sep = ""
  )
  quit(status = status)
}

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args || "-h" %in% args) usage()

take_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop(flag, " needs a value")
  args[[idx + 1]]
}

take_all <- function(flag) {
  idx <- which(args == flag)
  if (!length(idx)) return(character())
  if (any(idx == length(args))) stop(flag, " needs a value")
  args[idx + 1]
}

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) return(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
  normalizePath("WorkingTests/ddm_lnr_uc_shootout.R", mustWork = TRUE)
}

quiet_value <- function(expr) {
  value <- NULL
  capture.output(value <- suppressMessages(suppressWarnings(force(expr))))
  value
}

calc_ll_context <- function(emc, n_particles, seed = 12345) {
  model <- emc[[1]]$model()
  p_types <- names(model$p_types)
  dadm <- emc[[1]]$data[[1]]

  designs <- list()
  for (p in p_types) {
    designs[[p]] <- attr(dadm, "designs")[[p]][attr(attr(dadm, "designs")[[p]], "expand"), , drop = FALSE]
  }
  constants <- attr(dadm, "constants")
  if (is.null(constants)) constants <- NA

  set.seed(seed)
  p_mat <- matrix(rnorm(n_particles * length(p_types), sd = 0.25), ncol = length(p_types))
  colnames(p_mat) <- p_types

  old <- EMC2:::calc_ll(p_mat, dadm, constants = constants, designs = designs,
                        type = model$c_name, model$bound, model$transform,
                        model$pre_transform, p_types = p_types,
                        min_ll = log(1e-10), model$trend)
  new <- EMC2:::calc_ll_oo(p_mat, dadm, constants = constants, designs = designs,
                           type = model$c_name, model$bound, model$transform,
                           model$pre_transform, p_types = p_types,
                           min_ll = log(1e-10), model$trend)

  list(
    old = old, new = new, p_mat = p_mat, dadm = dadm, constants = constants,
    designs = designs, model = model, p_types = p_types
  )
}

time_reps <- function(ctx, reps) {
  gc()
  timings <- vapply(seq_len(reps), function(i) {
    system.time(
      EMC2:::calc_ll_oo(ctx$p_mat, ctx$dadm, constants = ctx$constants,
                        designs = ctx$designs, type = ctx$model$c_name,
                        bounds = ctx$model$bound, transforms = ctx$model$transform,
                        pretransforms = ctx$model$pre_transform,
                        p_types = ctx$p_types, min_ll = log(1e-10),
                        trend = ctx$model$trend)
    )[["elapsed"]]
  }, numeric(1))
  c(min_sec = min(timings), median_sec = stats::median(timings))
}

benchmark_one_build <- function(lib, label, particles, reps) {
  .libPaths(c(normalizePath(lib, mustWork = TRUE), .libPaths()))
  library(EMC2)

  ADmat <- matrix(c(-1 / 2, 1 / 2), ncol = 1, dimnames = list(NULL, "d"))
  matchfun <- function(d) d$S == d$lR
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2], ]
  dat$subjects <- droplevels(dat$subjects)

  design_DDM <- quiet_value(design(
    data = dat, model = DDM,
    formula = list(v ~ 1, a ~ 1, t0 ~ 1, s ~ 1, Z ~ 1, sv ~ 1, SZ ~ 1, st0 ~ 1),
    constants = c(s = 1)
  ))
  design_LNR <- quiet_value(design(
    data = dat, model = LNR, matchfun = matchfun,
    formula = list(m ~ lM, s ~ 1, t0 ~ 1),
    contrasts = list(m = list(lM = ADmat))
  ))

  rt_resolution <- 1/60
  make_case <- function(model_name, design_obj, uc) {
    data_uc <- if (is.infinite(uc)) dat else make_missing(dat, UC = uc, verbose = FALSE, rt_resolution = rt_resolution)
    quiet_value(make_emc(data_uc, design_obj, rt_resolution = rt_resolution, n_chains = 2, compress = TRUE))
  }

  cases <- list(
    DDM = list(design = design_DDM, ucs = c(base = Inf, mild = 1.05, moderate = 0.7)),
    LNR = list(design = design_LNR, ucs = c(base = Inf, mild = 1.05, moderate = 0.7))
  )

  rows <- list()
  for (model_name in names(cases)) {
    base_min <- NA_real_
    for (case_name in names(cases[[model_name]]$ucs)) {
      uc <- unname(cases[[model_name]]$ucs[[case_name]])
      emc <- make_case(model_name, cases[[model_name]]$design, uc)
      ctx <- calc_ll_context(emc, particles)
      timing <- time_reps(ctx, reps)
      if (case_name == "base") base_min <- timing[["min_sec"]]

      rows[[length(rows) + 1]] <- data.frame(
        build = label,
        lib = normalizePath(lib, mustWork = TRUE),
        model = model_name,
        case = case_name,
        UC = uc,
        particles = particles,
        reps = reps,
        input_nonfinite = if (is.infinite(uc)) 0L else sum(!is.finite(make_missing(dat, UC = uc, verbose = FALSE, rt_resolution = rt_resolution)$rt)),
        dadm_rows = nrow(ctx$dadm),
        dadm_nonfinite = sum(!is.finite(ctx$dadm$rt)),
        all_finite_attr = identical(attr(ctx$dadm, "emc2_all_finite_trials"), TRUE),
        old_new_max_abs = max(abs(ctx$old - ctx$new), na.rm = TRUE),
        min_sec = unname(timing[["min_sec"]]),
        median_sec = unname(timing[["median_sec"]]),
        ratio_to_base_min = unname(timing[["min_sec"]] / base_min),
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
}

worker <- "--worker" %in% args
particles <- as.integer(take_arg("--particles", "10000"))
reps <- as.integer(take_arg("--reps", "7"))
out <- take_arg("--out", NULL)

if (worker) {
  lib <- take_arg("--lib")
  label <- take_arg("--label", basename(normalizePath(lib, mustWork = TRUE)))
  result <- benchmark_one_build(lib, label, particles, reps)
  if (is.null(out)) {
    utils::write.csv(result, row.names = FALSE)
  } else {
    utils::write.csv(result, out, row.names = FALSE)
  }
  quit(status = 0)
}

libs <- take_all("--lib")
if (!length(libs)) usage(status = 1)
labels <- take_all("--label")
if (length(labels) < length(libs)) {
  labels <- c(labels, basename(normalizePath(libs[(length(labels) + 1):length(libs)], mustWork = TRUE)))
}
labels <- labels[seq_along(libs)]

tmp_files <- file.path(tempdir(), paste0("emc2_uc_shootout_", seq_along(libs), ".csv"))
for (i in seq_along(libs)) {
  cmd_args <- c(
    "--vanilla", script_path(), "--worker",
    "--lib", libs[[i]], "--label", labels[[i]],
    "--particles", as.character(particles),
    "--reps", as.character(reps),
    "--out", tmp_files[[i]]
  )
  status <- system2(file.path(R.home("bin"), "Rscript"), cmd_args)
  if (!identical(status, 0L)) stop("Benchmark worker failed for ", labels[[i]], " at ", libs[[i]])
}

combined <- do.call(rbind, lapply(tmp_files, utils::read.csv, stringsAsFactors = FALSE))
if (is.null(out)) {
  utils::write.csv(combined, row.names = FALSE)
} else {
  utils::write.csv(combined, out, row.names = FALSE)
  message("Wrote ", normalizePath(out, mustWork = FALSE))
}
