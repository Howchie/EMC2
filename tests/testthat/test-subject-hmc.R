# Subject-conditional HMC (R/subject_hmc.R).
local_rng_guard()

hmc_banana <- function(b = 1.5, s = 0.4) {
  list(ll = function(X) stats::dnorm(X[, 2], b * (X[, 1]^2 - 1), s, log = TRUE),
       mu = c(0, 0), V = diag(c(1, 4)),
       H = matrix(c(4 * b^2 / s^2 * 0.2, 0, 0, 1 / s^2), 2), b = b, s = s)
}

# Exact draws from N(a; mu, V) N(a2; b (a1^2 - 1), s^2): a1 by rejection, then a2 | a1.
hmc_banana_exact <- function(tg, n) {
  out <- numeric(0)
  vy <- tg$V[2, 2] + tg$s^2
  while (length(out) < n) {
    x <- stats::rnorm(4 * n)
    k <- stats::dnorm(tg$b * (x^2 - 1), 0, sqrt(vy)) / stats::dnorm(0, 0, sqrt(vy))
    out <- c(out, x[stats::runif(length(x)) < k])
  }
  x <- out[seq_len(n)]
  prec <- 1 / tg$V[2, 2] + 1 / tg$s^2
  cbind(x, stats::rnorm(n, (tg$b * (x^2 - 1) / tg$s^2) / prec, sqrt(1 / prec)))
}

test_that("the finite-difference force is the whitened gradient of the target", {
  local_test_rng(1)
  p <- 4
  A <- matrix(stats::rnorm(p * p), p)
  PL <- solve(crossprod(A) / p + diag(0.3, p))
  m <- stats::rnorm(p); mu <- stats::rnorm(p); V <- diag(stats::runif(p, 0.5, 2))
  ll <- function(X) { D <- sweep(X, 2, m); -0.5 * rowSums((D %*% PL) * D) }
  fr <- EMC2:::.emc_hmc_frame(PL, mu, V)
  expect_equal(fr$S %*% t(fr$S), solve(PL + solve(V)), tolerance = 1e-10)
  z <- stats::rnorm(p)
  e <- EMC2:::.emc_hmc_eval(z, fr, ll, 1e-4, -Inf)
  a <- mu + drop(fr$S %*% z)
  exact <- drop(crossprod(fr$S, -PL %*% (a - m) - solve(V, a - mu)))
  # A forward difference: exact up to h times the whitened curvature.
  expect_lt(max(abs(e$g - exact)), 1e-3)
  expect_equal(e$rows, p + 1L)
  # A known centre likelihood is not re-evaluated.
  expect_equal(EMC2:::.emc_hmc_eval(z, fr, ll, 1e-4, -Inf, ll_centre = e$ll)$rows, p)
})

test_that("leapfrog with the finite-difference force is reversible", {
  tg <- hmc_banana()
  fr <- EMC2:::.emc_hmc_frame(tg$H, tg$mu, tg$V)
  leap <- function(z0, v0, eps, L) {
    cur <- EMC2:::.emc_hmc_eval(z0, fr, tg$ll, 1e-5, -Inf)
    v <- v0 + 0.5 * eps * cur$g; z <- z0
    for (l in seq_len(L)) {
      z <- z + eps * v
      new <- EMC2:::.emc_hmc_eval(z, fr, tg$ll, 1e-5, -Inf)
      v <- v + (if (l < L) eps else 0.5 * eps) * new$g
    }
    list(z = z, v = v)
  }
  z0 <- c(0.3, -1.2); v0 <- c(0.7, 0.4)
  f <- leap(z0, v0, 0.13, 12)
  b <- leap(f$z, -f$v, 0.13, 12)
  expect_lt(max(abs(b$z - z0)), 1e-9)
  expect_lt(max(abs(b$v + v0)), 1e-9)
})

test_that("energy error shrinks as eps^2 at a fixed integration time", {
  tg <- hmc_banana()
  fr <- EMC2:::.emc_hmc_frame(tg$H, tg$mu, tg$V)
  med <- vapply(c(0.2, 0.1, 0.05), function(eps) {
    local_test_rng(3)
    stats::median(replicate(100, {
      a0 <- c(stats::rnorm(1, 0, 0.7), stats::rnorm(1))
      abs(EMC2:::.emc_hmc_transition(a0, tg$ll(matrix(a0, 1)), fr, tg$ll, eps, 1e-5, -Inf)$energy)
    }))
  }, numeric(1))
  slope <- unname(coef(stats::lm(log(med) ~ log(c(0.2, 0.1, 0.05))))[2])
  expect_gt(slope, 1.7)
  expect_lt(slope, 2.3)
})

test_that("a transition draws the same number of variates whatever happens", {
  tg <- hmc_banana()
  fr <- EMC2:::.emc_hmc_frame(tg$H, tg$mu, tg$V)
  fl <- -1e6
  llt <- function(X) ifelse(X[, 1] > 0, tg$ll(X), fl - 1)
  after <- function(a0, eps) {
    local_test_rng(5)
    tr <- EMC2:::.emc_hmc_transition(a0, llt(matrix(a0, 1)), fr, llt, eps, 1e-5, fl)
    list(tr = tr, next_u = stats::runif(1))
  }
  local_test_rng(5); stats::runif(1); stats::rnorm(2); stats::runif(1); ref <- stats::runif(1)
  edge <- after(c(1e-6, 0), 0.5)   # a probe outside the bound: aborts at once
  expect_identical(edge$tr$status, "abort")
  expect_identical(edge$tr$alpha, c(1e-6, 0))
  expect_identical(edge$next_u, ref)
  ok <- after(c(0.8, 0), 0.3)
  expect_identical(ok$next_u, ref)
})

test_that("the kernel leaves a curved target invariant and a broken one does not", {
  skip_model_validation()
  tg <- hmc_banana()
  fr <- EMC2:::.emc_hmc_frame(tg$H, tg$mu, tg$V)
  apply_kernel <- function(X, eps, broken) t(apply(X, 1L, function(a0) {
    for (k in 1:3) {
      l0 <- tg$ll(matrix(a0, 1))
      if (!broken) {
        a0 <- EMC2:::.emc_hmc_transition(a0, l0, fr, tg$ll, eps, 1e-5, -Inf)$alpha
      } else {
        # Acceptance without the kinetic energy: not a valid correction.
        e <- eps * stats::runif(1, 0.9, 1.1); L <- ceiling((pi / 2) / e)
        v0 <- stats::rnorm(2); u <- stats::runif(1)
        cur <- EMC2:::.emc_hmc_eval(drop(fr$R %*% (a0 - fr$mu)), fr, tg$ll, 1e-5, -Inf, l0)
        v <- v0 + 0.5 * e * cur$g; z <- cur$z
        for (l in seq_len(L)) {
          z <- z + e * v
          new <- EMC2:::.emc_hmc_eval(z, fr, tg$ll, 1e-5, -Inf)
          v <- v + (if (l < L) e else 0.5 * e) * new$g
        }
        if (u < min(1, exp(new$lp - cur$lp))) a0 <- new$a
      }
    }
    a0
  }))
  shift <- function(eps, broken) {
    local_test_rng(11)
    X0 <- hmc_banana_exact(tg, 4000)
    X1 <- apply_kernel(X0, eps, broken)
    stat <- function(X) c(mean(X[, 1]^2), mean(X[, 2]), mean(X[, 2] > 1), mean(X[, 1]^2 * X[, 2]))
    d <- vapply(1:200, function(b) { i <- sample(nrow(X0), replace = TRUE); stat(X1[i, ]) - stat(X0[i, ]) },
                numeric(4))
    max(abs((stat(X1) - stat(X0)) / apply(d, 1L, stats::sd)))
  }
  expect_lt(shift(0.5, FALSE), 3.5)
  expect_lt(shift(1.2, FALSE), 3.5)
  expect_gt(shift(1.2, TRUE), 3.5)
})

test_that("unsupported hierarchies take no HMC and a round without rates draws nothing", {
  s <- list(type = "factor", nuisance = rep(FALSE, 2), n_pars = 2, par_names = c("a", "b"),
            data = list(data.frame(x = 1)))
  expect_false(EMC2:::.emc_hmc_support(s)$ok)
  s$type <- "standard"
  expect_true(EMC2:::.emc_hmc_support(s)$ok)
  # Marginalised and pooled samplers store complete subject vectors and their
  # likelihoods, so HMC moves them like any other.
  s$marginalise <- list(param = "a")
  expect_true(EMC2:::.emc_hmc_support(s)$ok)
  s$marginalise <- NULL; s$pooled <- list(value = c(c = 0))
  expect_true(EMC2:::.emc_hmc_support(s)$ok)
  s$pooled <- NULL
  local_test_rng(9)
  streams <- list(.Random.seed, .Random.seed)
  before <- get(".Random.seed", envir = globalenv())
  props <- rbind(matrix(1:4, 2), c(-3, -4))
  pars <- list(tmu = c(a = 0, b = 0), tvar = diag(2))
  plan <- list(NULL, list(prec = diag(2), eps = 0.5, rate = 0, h = 1e-4, floor = -Inf))
  out <- EMC2:::.emc_hmc_round(plan, props, pars, "standard", streams, NULL, NULL,
                               par_names = c("a", "b"))
  expect_identical(get(".Random.seed", envir = globalenv()), before)
  expect_identical(out$streams, streams)
  expect_identical(out$proposals, props)
  expect_true(all(vapply(out$tasks, is.null, logical(1))))
  # A subject whose joint draw with the particle step declined HMC draws
  # nothing more and takes no transition.
  plan[[1]] <- plan[[2]] <- list(prec = diag(2), eps = 0.5, rate = 0.5, h = 1e-4, floor = -Inf)
  out <- EMC2:::.emc_hmc_round(plan, props, pars, "standard", streams, NULL, NULL,
                               par_names = c("a", "b"), pick = c(FALSE, FALSE))
  expect_identical(out$streams, streams)
  expect_identical(out$proposals, props)
  expect_true(all(vapply(out$tasks, is.null, logical(1))))
  expect_equal(EMC2:::.emc_hmc_plan_rates(plan, 2), c(0.5, 0.5))
  expect_equal(EMC2:::.emc_hmc_plan_rates(NULL, 2), c(0, 0))
})

hmc_small_fit <- function(route = "static", stages = c(20, 50), seed = 17, hmc = TRUE,
                          rate = 0.5) {
  withr::local_options(emc2.subject_hmc = hmc, emc2.subject_hmc_rate = rate,
                       emc2.worker_queue = route == "queued")
  if (route == "serial") local_mocked_bindings(.emc_wpool_start = function(...) NULL,
                                               .package = "EMC2")
  dat <- forstmann[forstmann$subjects %in% levels(forstmann$subjects)[1:4], ]
  dat$subjects <- droplevels(dat$subjects)
  utils::capture.output(des <- suppressMessages(design(
    data = dat, model = LNR, matchfun = function(d) d$S == d$lR,
    formula = list(m ~ lM, s ~ 1, t0 ~ 1))))
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  emc <- suppressMessages(make_emc(dat, des, n_chains = 1, compress = TRUE))
  sc <- function(i) list(iter = i, max_gd = Inf, min_unique = 0, min_es = 0)
  emc <- suppressMessages(run_emc(emc, stage = "preburn", cores_for_chains = 1,
    cores_per_chain = 2, step_size = stages[1], max_tries = 1, verbose = FALSE,
    stop_criteria = sc(stages[1])))
  if (length(stages) > 1L)
    emc <- suppressMessages(run_emc(emc, stage = "burn", cores_for_chains = 1,
      cores_per_chain = 2, step_size = 25, max_tries = 1, verbose = FALSE,
      stop_criteria = sc(stages[2])))
  emc
}

test_that("a subject's warmup stops after a fixed number of transitions", {
  ps <- list(eps = 0.5, rate = EMC2:::.EMC_HMC_WARMUP_RATE, warmup = TRUE, n_warm = 0L)
  n <- EMC2:::.EMC_HMC_WARMUP_MOVES
  for (k in seq_len(n - 1L)) ps <- EMC2:::.emc_hmc_warmup_step(ps, list(acc_prob = 0.95))
  expect_true(ps$warmup)
  expect_gt(ps$rate, 0)
  # Acceptance above target lengthens the step.
  expect_gt(ps$eps, 0.5)
  ps <- EMC2:::.emc_hmc_warmup_step(ps, list(acc_prob = NA_real_))
  expect_false(ps$warmup)
  expect_equal(ps$rate, 0)
})

test_that("the first burn block moves every eligible subject while its step size adapts", {
  skip_on_os("windows")
  skip_on_cran()
  local_rng_guard()
  fit <- hmc_small_fit("static", stages = c(20, 45), rate = NULL)
  h <- fit[[1]]$subject_hmc
  first <- fit[[1]]$alloc$history[[1]]
  eligible <- vapply(h$subjects, function(x) isTRUE(x$eligible), logical(1))
  expect_true(any(eligible))
  expect_identical(first$stage, "burn")
  expect_equal(first$rate[eligible], rep(EMC2:::.EMC_HMC_WARMUP_RATE, sum(eligible)))
  # No separate calibration: the moves themselves set the step size.
  expect_equal(first$calib_cpu, 0)
  expect_gt(h$totals$attempts, 0L)
  p <- fit[[1]]$n_pars
  for (s in which(eligible)) {
    ev <- h$evidence[[s]]
    expect_gt(ev$n, 0L)
    expect_identical(ev$key, EMC2:::.emc_hmc_key(h$subjects[[s]]$eps))
    expect_false(isTRUE(all.equal(h$subjects[[s]]$eps, p^(-1 / 4))))
  }
})

test_that("serial, static and queued HMC rounds give the same chain", {
  skip_on_os("windows")
  skip_on_cran()
  local_rng_guard()
  # Keep the HMC kernel fixed across worker routes. The production allocator
  # prices shear and collapsed moves from measured CPU time, so comparing full
  # adaptive fits can select different move plans even with identical seeds.
  withr::local_options(emc2.group_move = "none", emc2.collapsed_move = FALSE)
  # The ensemble's between-block plan prices pool members by wall time too.
  local_mocked_bindings(.emc_ensemble_plan_update = function(sampler) sampler, .package = "EMC2")
  static <- hmc_small_fit("static")
  h <- static[[1]]$subject_hmc
  expect_true(h$supported)
  expect_gt(h$totals$attempts, 0L)
  expect_gt(h$totals$accepted, 0L)
  for (route in c("serial", "queued")) {
    other <- hmc_small_fit(route)
    expect_identical(other[[1]]$samples$alpha, static[[1]]$samples$alpha, label = route)
    expect_identical(other[[1]]$samples$subj_ll, static[[1]]$samples$subj_ll, label = route)
    expect_identical(other[[1]]$rng, static[[1]]$rng, label = route)
    expect_identical(other[[1]]$subject_hmc$totals[c("attempts", "accepted", "aborts", "rows")],
                     h$totals[c("attempts", "accepted", "aborts", "rows")], label = route)
  }
})

test_that("a transition is a function of its input and RNG state", {
  skip_on_os("windows")
  skip_on_cran()
  local_rng_guard()
  emc <- hmc_small_fit("static")
  ch <- emc[[1]]
  withr::local_options(emc2.subject_hmc = TRUE, emc2.subject_hmc_rate = 1)
  plan <- EMC2:::.emc_hmc_worker_plan(ch)
  s <- which(!vapply(plan, is.null, logical(1)))[1]
  skip_if(is.na(s), "no eligible subject")
  prior <- EMC2:::.emc_burn_subject_prior(ch, s)
  run <- function() {
    assign(".Random.seed", ch$rng$subjects[[s]], envir = globalenv())
    tk <- EMC2:::.emc_hmc_run(ch$samples$alpha[, s, ch$samples$idx],
                              ch$samples$subj_ll[s, ch$samples$idx], ch$data[[s]], ch$model,
                              prior$mu, prior$var, plan[[s]], 1L, ch$par_names)
    tk$cpu <- NULL
    list(tk = tk, seed = get(".Random.seed", envir = globalenv()))
  }
  first <- run(); again <- run()
  expect_identical(first, again)
  expect_true(first$tk$attempted)
  # The returned likelihood is the likelihood of the returned state.
  ll <- EMC2:::calc_ll_manager(matrix(first$tk$alpha, 1, dimnames = list(NULL, ch$par_names)),
                               dadm = ch$data[[s]], model = ch$model)
  expect_equal(as.numeric(ll), unname(first$tk$ll))
})

test_that("a split run continues exactly and an old checkpoint starts without HMC", {
  skip_on_os("windows")
  skip_on_cran()
  local_rng_guard()
  # The continuation assertion is about saved HMC state and RNG streams.
  # Keep CPU-priced group moves off so wall-time noise cannot change the
  # transition plan between the whole and split runs.
  withr::local_options(emc2.group_move = "none", emc2.collapsed_move = FALSE)
  whole <- hmc_small_fit("static", stages = c(20, 50))
  half <- hmc_small_fit("static", stages = c(20, 25))
  half <- EMC2:::restore_duplicates(EMC2:::strip_duplicates(half))
  withr::local_options(emc2.subject_hmc = TRUE, emc2.subject_hmc_rate = 0.5)
  rest <- suppressMessages(run_emc(half, stage = "burn", cores_for_chains = 1,
    cores_per_chain = 2, step_size = 25, max_tries = 1, verbose = FALSE,
    stop_criteria = list(iter = 25, max_gd = Inf, min_unique = 0, min_es = 0)))
  expect_identical(rest[[1]]$samples$alpha, whole[[1]]$samples$alpha)
  expect_identical(rest[[1]]$subject_hmc$totals[c("attempts", "accepted")],
                   whole[[1]]$subject_hmc$totals[c("attempts", "accepted")])
  # Without HMC metadata a block runs with HMC off until a rebuild initialises it.
  old <- half
  old[[1]]$subject_hmc <- NULL
  expect_null(EMC2:::.emc_hmc_worker_plan(old[[1]]))
})

test_that("the HMC switch and its tuning are recorded production options", {
  withr::local_options(list(emc2.subject_hmc = TRUE, emc2.subject_hmc_h = 1e-5))
  rec <- EMC2:::.emc_record_production_options(list(list(n_pars = 2L)))
  expect_true(rec$options$emc2.subject_hmc)
  expect_identical(rec$options$emc2.subject_hmc_h, 1e-5)
  withr::local_options(list(emc2.subject_hmc = FALSE))
  again <- EMC2:::.emc_record_production_options(rec$emc)
  expect_true(again$options$emc2.subject_hmc)
})

test_that("step sizes stay in their bounds", {
  ev <- list(n = 10L, accp = rep(0.05, 10))
  sj <- list(eps = 0.03, n_tuned = 0L)
  expect_equal(EMC2:::.emc_hmc_update_eps(sj, ev)$eps, 0.02)
  expect_equal(EMC2:::.emc_hmc_clip_eps(5), 2)
})
