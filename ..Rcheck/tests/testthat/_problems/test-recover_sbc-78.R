# Extracted from test-recover_sbc.R:78

# prequel ----------------------------------------------------------------------
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
  temp_dir
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
  temp_dir
}

# test -------------------------------------------------------------------------
dir <- tempfile("sbc_")
dir.create(dir)
par_names <- c("m", "m_lMTRUE", "s", "t0")
tdir <- write_single_run(dir, par_names, reps_idx = 1:5, replicates = 5)
save_to <- file.path(dir, "out.RData")
SBC <- recover_sbc(tdir, par_names, fileName = save_to, type = "single", verbose = FALSE)
expect_named(SBC, c("rank", "med", "bias", "coverage"))
expect_equal(nrow(SBC$rank$alpha), 5)
expect_equal(colnames(SBC$rank$alpha), par_names)
expect_equal(attr(SBC, "recovered_reps"), 1:5)
ref <- split_list_to_dfs(lapply(1:5, function(i) make_single_rep(par_names, i)))
