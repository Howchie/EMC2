# Tests for recover_sbc() and the shared SBC assembly helpers.
# These fabricate the intermediate rep_<i>.rds / prior_samples.rds files that
# run_sbc() writes, so no model fitting is needed: the assembly logic is
# exercised directly and deterministically.

# --- helpers to fabricate intermediate files --------------------------------

make_single_rep <- function(par_names, seed) {
  set.seed(seed)
  p <- length(par_names)
  list(
    rank     = setNames(runif(p), par_names),
    med      = setNames(rnorm(p), par_names),
    bias     = setNames(rnorm(p), par_names),
    coverage = setNames(runif(p) > 0.1, par_names)
  )
}

make_hier_rep <- function(par_names, var_names, seed) {
  set.seed(seed)
  list(
    rank_mu_row   = setNames(runif(length(par_names)), NULL),
    rank_var_row  = setNames(runif(length(var_names)), NULL),
    var_col_names = var_names,
    rand_effects  = matrix(rnorm(length(par_names) * 2), nrow = 2,
                           dimnames = list(NULL, par_names))
  )
}

write_single_run <- function(dir, par_names, reps_idx, replicates,
                             write_prior = TRUE) {
  temp_dir <- file.path(dir, "run_temp")
  dir.create(temp_dir, showWarnings = FALSE)
  for (i in reps_idx)
    saveRDS(make_single_rep(par_names, i), file.path(temp_dir, paste0("rep_", i, ".rds")))
  if (write_prior) {
    prior_alpha <- matrix(rnorm(replicates * length(par_names)),
                          nrow = replicates, dimnames = list(NULL, par_names))
    saveRDS(prior_alpha, file.path(temp_dir, "prior_samples.rds"))
  }
  file.path(dir, "run.RData")
}

write_hier_run <- function(dir, par_names, var_names, reps_idx, replicates,
                           write_prior = TRUE) {
  temp_dir <- file.path(dir, "run_temp")
  dir.create(temp_dir, showWarnings = FALSE)
  for (i in reps_idx)
    saveRDS(make_hier_rep(par_names, var_names, i), file.path(temp_dir, paste0("rep_", i, ".rds")))
  if (write_prior) {
    prior_mu  <- matrix(rnorm(length(par_names) * replicates), nrow = length(par_names),
                        dimnames = list(par_names, NULL))
    prior_var <- array(rnorm(length(par_names)^2 * replicates),
                       dim = c(length(par_names), length(par_names), replicates))
    saveRDS(list(prior_mu = prior_mu, prior_var = prior_var),
            file.path(temp_dir, "prior_samples.rds"))
  }
  file.path(dir, "run.RData")
}

# --- single ------------------------------------------------------------------

test_that("recover_sbc reassembles a complete single run", {
  dir <- tempfile("sbc_"); dir.create(dir)
  par_names <- c("m", "m_lMTRUE", "s", "t0")
  fileName <- write_single_run(dir, par_names, reps_idx = 1:5, replicates = 5)

  out_path <- file.path(dir, "out.RData")
  SBC <- recover_sbc(fileName, type = "single", pathName = out_path, verbose = FALSE)
  expect_named(SBC, c("rank", "med", "bias", "coverage"))
  expect_equal(nrow(SBC$rank$alpha), 5)
  expect_equal(colnames(SBC$rank$alpha), par_names)
  expect_equal(attr(SBC, "recovered_reps"), 1:5)

  # rows match what split_list_to_dfs would produce directly from the reps
  ref <- split_list_to_dfs(lapply(1:5, function(i) make_single_rep(par_names, i)))
  expect_equal(unname(SBC$rank$alpha), unname(ref$rank$alpha))

  # saved file round-trips in run_sbc's single format (SBC + prior_alpha)
  expect_true(file.exists(out_path))
  env <- new.env(); load(out_path, envir = env)
  expect_true(all(c("SBC", "prior_alpha") %in% ls(env)))
  expect_equal(env$SBC$rank$alpha, SBC$rank$alpha)
})

test_that("recover_sbc reassembles a partial single run and reports missing", {
  dir <- tempfile("sbc_"); dir.create(dir)
  par_names <- c("m", "s", "t0")
  fileName <- write_single_run(dir, par_names, reps_idx = c(1, 2, 5), replicates = 10)

  expect_message(SBC <- recover_sbc(fileName, type = "single",
                                    pathName = file.path(dir, "out.RData")),
                 "Recovered 3 single replicate")
  expect_equal(nrow(SBC$rank$alpha), 3)
  expect_equal(attr(SBC, "recovered_reps"), c(1, 2, 5))
})

test_that("recover_sbc auto-detects single and works without a prior file", {
  dir <- tempfile("sbc_"); dir.create(dir)
  par_names <- c("m", "s")
  fileName <- write_single_run(dir, par_names, reps_idx = 1:4, replicates = 4,
                               write_prior = FALSE)

  # inspect-only message when prior draws are absent
  expect_message(SBC <- recover_sbc(fileName, pathName = file.path(dir, "out.RData")),
                 "Inspect-only")
  expect_equal(nrow(SBC$rank$alpha), 4)
  expect_equal(attr(SBC, "recovered_reps"), 1:4)
})

test_that("recover_sbc skips failed single replicates", {
  dir <- tempfile("sbc_"); dir.create(dir)
  par_names <- c("m", "s")
  temp_dir <- file.path(dir, "run_temp")
  dir.create(temp_dir)
  saveRDS(make_single_rep(par_names, 1), file.path(temp_dir, "rep_1.rds"))
  saveRDS(list(rank = NULL, med = NULL, bias = NULL, coverage = NULL, failed = TRUE),
          file.path(temp_dir, "rep_2.rds"))
  saveRDS(make_single_rep(par_names, 3), file.path(temp_dir, "rep_3.rds"))

  SBC <- recover_sbc(file.path(dir, "run.RData"), type = "single",
                     pathName = file.path(dir, "out.RData"), verbose = FALSE)
  expect_equal(nrow(SBC$rank$alpha), 2)
  expect_equal(attr(SBC, "recovered_reps"), c(1, 3))
})

# --- hierarchical ------------------------------------------------------------

test_that("recover_sbc reassembles a complete hierarchical run", {
  dir <- tempfile("sbc_"); dir.create(dir)
  par_names <- c("v", "a", "t0")
  var_names <- c("v.v", "a.a", "t0.t0")
  fileName <- write_hier_run(dir, par_names, var_names, reps_idx = 1:6, replicates = 6)

  out_path <- file.path(dir, "out.RData")
  out <- recover_sbc(fileName, type = "hierarchical", pathName = out_path, verbose = FALSE)
  expect_named(out, c("rank", "prior", "rand_effects"))
  expect_equal(nrow(out$rank$mu), 6)
  expect_equal(colnames(out$rank$mu), par_names)         # from prior_mu rownames
  expect_equal(colnames(out$rank$var), var_names)        # from rep var_col_names
  expect_false(is.null(out$prior$mu))
  expect_length(out$rand_effects, 6)
  expect_equal(attr(out, "recovered_reps"), 1:6)

  # saved file round-trips in run_sbc's hierarchical format (SBC_temp)
  expect_true(file.exists(out_path))
  env <- new.env(); load(out_path, envir = env)
  expect_true("SBC_temp" %in% ls(env))
  expect_equal(env$SBC_temp$rank$mu, out$rank$mu)
})

test_that("recover_sbc recovers hierarchical ranks without prior (prior NULL + warning)", {
  dir <- tempfile("sbc_"); dir.create(dir)
  par_names <- c("v", "a")
  var_names <- c("v.v", "a.a")
  fileName <- write_hier_run(dir, par_names, var_names, reps_idx = 1:3, replicates = 3,
                             write_prior = FALSE)

  expect_warning(out <- recover_sbc(fileName, type = "hierarchical",
                                    pathName = file.path(dir, "out.RData"), verbose = FALSE),
                 "Prior draws not found")
  expect_equal(nrow(out$rank$mu), 3)
  expect_equal(colnames(out$rank$var), var_names)
  expect_null(out$prior$mu)
})

# --- errors ------------------------------------------------------------------

test_that("recover_sbc errors when no rep files exist", {
  dir <- tempfile("sbc_"); dir.create(dir)
  expect_error(recover_sbc(file.path(dir, "run.RData"),
                           pathName = file.path(dir, "out.RData"), verbose = FALSE),
               "No replicate files found")
})
