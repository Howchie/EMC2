# Equivalence tests for the proposal-factorisation cache and the compressed
# design pass-through.  Both are pure optimisations: they must not change what
# the sampler computes, only how often it computes it.

test_that(".chol_factor reproduces the inline chol/backsolve it replaced", {
  set.seed(11)
  p <- 7
  A <- crossprod(matrix(rnorm(p * p), p))
  f <- EMC2:::.chol_factor(A, p)

  expect_equal(f$R, chol(A))
  expect_equal(f$rooti, backsolve(chol(A), diag(p)))
  expect_equal(f$sum_log_diag, sum(log(diag(backsolve(chol(A), diag(p))))))
})

test_that(".chol_factor keeps the escalating-ridge fallback for singular input", {
  p <- 4
  # Exactly singular: chol() fails, the 1e-6 ridge must rescue it.
  A <- matrix(1, p, p)
  f <- EMC2:::.chol_factor(A, p)
  expect_equal(f$R, chol(A + diag(1e-6, p)))

  # Not even the ridges help (negative diagonal) -> diagonal fallback.
  B <- diag(c(-1, -1, -1, -1))
  g <- EMC2:::.chol_factor(B, p)
  expect_equal(g$R, diag(sqrt(rep(1e-4, p)), p))

  expect_null(EMC2:::.chol_factor(matrix(numeric(0), 0, 0), 0))
})

test_that("epsilon scaling factors out of the cached factorisation", {
  # The whole point of caching: an adapting epsilon rescales R and rooti
  # analytically instead of forcing a fresh chol()/backsolve().
  set.seed(12)
  p <- 6
  A <- crossprod(matrix(rnorm(p * p), p))
  eps <- 0.37
  f <- EMC2:::.chol_factor(A, p)

  R_direct <- chol(A) * eps
  rooti_direct <- backsolve(R_direct, diag(p))
  lc_direct <- sum(log(diag(rooti_direct))) - 0.5 * p * log(2 * pi)

  R_cached <- f$R * eps
  rooti_cached <- f$rooti / eps
  lc_cached <- f$sum_log_diag - p * log(eps) - 0.5 * p * log(2 * pi)

  expect_equal(R_cached, R_direct)
  expect_equal(rooti_cached, rooti_direct)
  expect_equal(lc_cached, lc_direct)

  # And the density these feed is the same as the unfactorised reference.
  x <- matrix(rnorm(3 * p), 3, p)
  mu <- rnorm(p)
  expect_equal(
    as.vector(EMC2:::fast_dmvnorm_rooti(x = x, mean = mu, rooti = rooti_cached,
                                        log_const = lc_cached)),
    as.vector(mvtnorm::dmvnorm(x, mu, A * eps^2, log = TRUE))
  )
})

test_that("component index and marginalise helpers agree with the sampler", {
  components <- c(1L, 1L, 2L, 2L, 2L)
  par_names <- c("a", "b", "c", "d", "t0")

  expect_equal(EMC2:::.marginal_par_idx(par_names, NULL), rep(FALSE, 5))
  marg <- list(param = "t0", mu = 0, sigma = 0.5)
  expect_equal(EMC2:::.marginal_par_idx(par_names, marg),
               c(FALSE, FALSE, FALSE, FALSE, TRUE))

  idx_list <- EMC2:::.component_idx_list(components, EMC2:::.marginal_par_idx(par_names, marg))
  expect_equal(idx_list[[1]], c(TRUE, TRUE, FALSE, FALSE, FALSE))
  # t0 is held out of component 2's proposal block
  expect_equal(idx_list[[2]], c(FALSE, FALSE, TRUE, TRUE, FALSE))

  V <- diag(5)
  Vm <- EMC2:::.apply_marginal_group_var(V, c(FALSE, FALSE, FALSE, FALSE, TRUE), marg)
  expect_equal(Vm[5, 5], 0.25)
  expect_equal(Vm[5, 1:4], rep(0, 4))
  expect_equal(Vm[1:4, 1:4], diag(4))
})

test_that("cached subject factorisations match direct factorisation", {
  set.seed(13)
  p <- 5
  chains_var <- crossprod(matrix(rnorm(p * p), p))
  eff_var <- crossprod(matrix(rnorm(p * p), p))
  idx_list <- EMC2:::.component_idx_list(rep(1L, p), rep(FALSE, p))

  cache <- EMC2:::build_subject_chol_cache(chains_var, eff_var, idx_list)
  expect_equal(cache$chains[[1]]$R, chol(chains_var))
  expect_equal(cache$eff[[1]]$R, chol(eff_var))
  expect_identical(cache$chains_ref, chains_var)

  # preburn has no chains_var/eff_var at all
  expect_null(EMC2:::build_subject_chol_cache(NULL, NULL, idx_list))

  gcache <- EMC2:::build_group_chol_cache(chains_var, idx_list)
  expect_equal(gcache$f[[1]]$R, chol(chains_var))
  expect_identical(gcache$ref, chains_var)
})

test_that("compressed designs give bit-identical likelihoods to expanded ones", {
  # calc_ll_manager now hands the C++ mapper compressed design matrices plus
  # their "expand" attribute instead of materialising full-length copies.
  # ParamTable::map_from_designs must produce exactly the same parameters.
  skip_on_os("windows")
  ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
  dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2], ]
  dat$subjects <- droplevels(dat$subjects)
  des <- design(data = dat, model = LNR, matchfun = function(d) d$S == d$lR,
                formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                contrasts = list(m = list(lM = ADmat)))
  emc <- make_emc(dat, des, rt_resolution = 0.05, n_chains = 1)
  dadm <- emc[[1]]$data[[1]]

  set.seed(14)
  pars <- sampled_pars(des, doMap = FALSE)
  props <- matrix(rep(0, length(pars)), nrow = 20, ncol = length(pars), byrow = TRUE,
                  dimnames = list(NULL, names(pars)))
  props <- props + matrix(rnorm(length(props), 0, 0.3), nrow = nrow(props))

  m <- emc[[1]]$model()
  const <- attr(dadm, "constants"); if (is.null(const)) const <- NA
  call_with <- function(designs) {
    EMC2:::calc_ll_oo(props, dadm, constants = const, designs = designs,
                      type = m$c_name, bounds = m$bound, transforms = m$transform,
                      pretransforms = m$pre_transform, p_types = names(m$p_types),
                      min_ll = log(1e-10), trend = m$trend, marginalise = NULL)
  }
  ll_expanded <- call_with(EMC2:::.oo_expanded_designs(dadm, expand = TRUE))
  ll_compressed <- call_with(EMC2:::.oo_expanded_designs(dadm, expand = FALSE))

  expect_identical(ll_compressed, ll_expanded)
  expect_true(all(is.finite(ll_expanded)))
})
