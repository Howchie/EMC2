# The C++ ensemble kernels (src/ensemble_move.cpp) against their R reference:
# the pool draws, the stacked pools' proposal densities, and the Gibbs rounds
# of every hierarchy they cover. Same seed, same random numbers in the same
# order, so the selections and the RNG state after the call are identical and
# the group state agrees to rounding.

ens_cpp_sampler <- function(type, p, n, par_groups = NULL) {
  pn <- paste0("v", seq_len(p))
  s <- list(type = type, n_pars = p, n_subjects = n, par_names = pn,
            nuisance = rep(FALSE, p), marginalised_idx = rep(FALSE, p),
            model = NULL, marginalise = NULL)
  s <- if (type == "diagonal-gamma") EMC2:::add_info_diag_gamma(s, prior = NULL) else
    EMC2:::add_info_standard(s, prior = NULL, par_groups = par_groups)
  s$samples <- list(theta_mu = matrix(0, p, 1), theta_var = array(diag(p), c(p, p, 1)),
                    last_theta_var_inv = diag(p), a_half = matrix(1, p, 1),
                    alpha = array(0, c(p, n, 1)), subj_ll = matrix(0, n, 1),
                    stage = "sample", idx = 1)
  s
}

# The pool proposal and its draws as .emc_ensemble_q / .emc_ensemble_draw_q
# define them, one subject after another.
ens_draw_pools_ref <- function(sampler, chains_mu, chains_var) {
  cfg <- EMC2:::.emc_ensemble_options()
  n <- sampler$n_subjects
  q <- lapply(seq_len(n), function(s)
    EMC2:::.emc_ensemble_q(as.numeric(chains_mu[[s]]), as.matrix(chains_var[[s]]), cfg))
  lq <- vector("list", n)
  rows <- lapply(seq_len(n), function(s) {
    fresh <- EMC2:::.emc_ensemble_draw_q(cfg$pool - 1L, q[[s]], cfg)
    lq[[s]] <<- EMC2:::.emc_ensemble_log_q_m2(attr(fresh, "m2"), q[[s]], cfg)
    attr(fresh, "m2") <- NULL
    colnames(fresh) <- sampler$par_names
    fresh
  })
  list(q = q, rows = rows, lq = lq)
}

ens_cpp_case <- function(type, p, n, par_groups = NULL, seed = 1) {
  set.seed(seed)
  s <- ens_cpp_sampler(type, p, n, par_groups)
  alpha <- matrix(rnorm(p * n, 0, 0.7), p, n, dimnames = list(s$par_names, NULL))
  chains_mu <- lapply(seq_len(n), function(k) alpha[, k] + rnorm(p, 0, 0.1))
  chains_var <- lapply(seq_len(n), function(k) {
    A <- matrix(rnorm(p * p, 0, 0.3), p); crossprod(A) + diag(0.2, p)
  })
  pars <- EMC2:::gibbs_step(s, alpha, type)
  list(s = s, alpha = alpha, chains_mu = chains_mu, chains_var = chains_var, pars = pars)
}

test_that("pool draws from C++ match the R reference draw for draw", {
  withr::local_options(list(emc2.ensemble_pool = 12L))
  cs <- ens_cpp_case("standard", 4, 5)
  set.seed(11)
  pools <- EMC2:::.emc_ensemble_draw_pools(cs$s, "sample", cs$chains_mu, cs$chains_var)
  seed_cpp <- .Random.seed
  set.seed(11)
  ref <- ens_draw_pools_ref(cs$s, cs$chains_mu, cs$chains_var)
  expect_identical(.Random.seed, seed_cpp)
  for (k in seq_len(cs$s$n_subjects)) {
    expect_equal(pools$rows[[k]], ref$rows[[k]], tolerance = 1e-12)
    expect_equal(pools$lq[[k]], ref$lq[[k]], tolerance = 1e-10)
    expect_equal(pools$L[, , k], ref$q[[k]]$L, tolerance = 1e-12, ignore_attr = TRUE)
    expect_equal(pools$q[[k]]$logdet, ref$q[[k]]$logdet, tolerance = 1e-12)
  }
  # The stack carries the current state's density from the same proposal.
  ll <- rnorm(cs$s$n_subjects)
  ll_new <- matrix(rnorm(11 * cs$s$n_subjects), 11)
  k0 <- sample.int(12L, cs$s$n_subjects, replace = TRUE)
  cfg <- EMC2:::.emc_ensemble_options()
  st <- EMC2:::emc_ensemble_stack_cpp(pools$rows, pools$lq, ll_new, rbind(cs$alpha, ll), k0,
                                      pools$L, pools$mu, pools$logdet, pools$c,
                                      cfg$tail_weight, cfg$tail_df)
  for (k in seq_len(cs$s$n_subjects)) {
    cur <- (k - 1L) * 12L + k0[k]
    expect_equal(st$pool_stack[cur, ], cs$alpha[, k], ignore_attr = TRUE)
    expect_equal(st$base_w[cur], ll[k] - EMC2:::.emc_ensemble_log_q(t(cs$alpha[, k]), ref$q[[k]], cfg),
                 tolerance = 1e-10)
    others <- st$pool_stack[(k - 1L) * 12L + seq_len(12L)[-k0[k]], ]
    expect_equal(others, ref$rows[[k]], tolerance = 1e-12, ignore_attr = TRUE)
  }
})

test_that("C++ Gibbs rounds reproduce the R rounds for every covered hierarchy", {
  withr::local_options(list(emc2.ensemble_pool = 10L, emc2.ensemble_rounds = 4L))
  cases <- list(list("standard", 3, 6, NULL), list("standard", 4, 7, c(1, 1, 2, 2)),
                list("standard", 4, 5, c(1, 2, 1, 3)), list("diagonal-gamma", 3, 6, NULL))
  for (cse in cases) {
    cs <- ens_cpp_case(cse[[1]], cse[[2]], cse[[3]], cse[[4]])
    s <- cs$s; p <- s$n_pars; n <- s$n_subjects; M <- 10L
    expect_true(EMC2:::.emc_ensemble_rounds_in_cpp(s))
    pools <- EMC2:::.emc_ensemble_draw_pools(s, "sample", cs$chains_mu, cs$chains_var)
    cfg <- EMC2:::.emc_ensemble_options()
    ll_new <- matrix(rnorm((M - 1L) * n, -2, 1), M - 1L)
    k0 <- sample.int(M, n, replace = TRUE)
    st <- EMC2:::emc_ensemble_stack_cpp(pools$rows, pools$lq, ll_new, rbind(cs$alpha, rnorm(n, -2, 1)),
                                        k0, pools$L, pools$mu, pools$logdet, pools$c,
                                        cfg$tail_weight, cfg$tail_df)
    stack <- st$pool_stack; colnames(stack) <- s$par_names
    set.seed(7)
    r_ref <- EMC2:::.emc_ensemble_rounds_r(s, cs$pars, stack, st$base_w, k0, M, cfg)
    seed_ref <- .Random.seed
    set.seed(7)
    r_cpp <- EMC2:::.emc_ensemble_rounds_cpp(s, cs$pars, stack, st$base_w, k0, M, cfg)
    label <- paste(cse[[1]], paste(cse[[4]], collapse = ""))
    expect_identical(.Random.seed, seed_ref, label = paste(label, "rng"))
    expect_identical(r_cpp$selected, r_ref$selected, label = paste(label, "selection"))
    expect_equal(r_cpp$wess, r_ref$wess, tolerance = 1e-8, label = paste(label, "wess"))
    expect_equal(r_cpp$pars$tmu, r_ref$pars$tmu, tolerance = 1e-8, label = paste(label, "tmu"))
    expect_equal(r_cpp$pars$tvar, r_ref$pars$tvar, tolerance = 1e-8, ignore_attr = TRUE)
    expect_equal(r_cpp$pars$tvinv, r_ref$pars$tvinv, tolerance = 1e-8, ignore_attr = TRUE)
    expect_identical(names(r_cpp$pars), names(r_ref$pars))
    if (cse[[1]] == "standard") {
      expect_equal(r_cpp$pars$a_half, r_ref$pars$a_half, tolerance = 1e-8)
      expect_equal(r_cpp$pars$subj_mu, r_ref$pars$subj_mu, tolerance = 1e-8, ignore_attr = TRUE)
    } else {
      expect_identical(names(r_cpp$pars$tmu), names(r_ref$pars$tmu))
    }
    # Control: a different number of rounds consumes different random numbers.
    set.seed(7)
    cfg2 <- cfg; cfg2$rounds <- cfg$rounds + 1L
    EMC2:::.emc_ensemble_rounds_cpp(s, cs$pars, stack, st$base_w, k0, M, cfg2)
    expect_false(identical(.Random.seed, seed_ref))
  }
})

test_that("the R rounds remain the path for hierarchies the C++ rounds do not cover", {
  s <- list(type = "factor")
  expect_false(EMC2:::.emc_ensemble_rounds_in_cpp(s))
  s <- list(type = "standard", gd = list(matrix(1, 3, 1)))
  expect_false(EMC2:::.emc_ensemble_rounds_in_cpp(s))
  s <- list(type = "standard", gd = NULL, group_designs = list(a = matrix(1, 3, 1)))
  expect_false(EMC2:::.emc_ensemble_rounds_in_cpp(s))
})

