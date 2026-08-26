testthat::test_that("continuous split-normal median geometry is exact", {
  skip_if_not(exists(".bawd_split_params", envir = asNamespace("EMC2"), inherits = FALSE),
             "split launch helpers are unavailable")
  split_params <- getFromNamespace(".bawd_split_params", "EMC2")
  mu <- c(-1, 0, 1)
  sigma <- c(.3, 1, 2)
  delta <- c(-3, 0, 4)
  sp <- split_params(mu, sigma, delta)
  cdf <- function(y, i) {
    if (y <= sp$c[i])
      2 * sp$sigma_L[i] / (sp$sigma_L[i] + sp$sigma_R[i]) *
        pnorm((y - sp$c[i]) / sp$sigma_L[i])
    else
      2 * (sp$sigma_L[i] / (sp$sigma_L[i] + sp$sigma_R[i]) * .5 +
             sp$sigma_R[i] / (sp$sigma_L[i] + sp$sigma_R[i]) *
             (pnorm((y - sp$c[i]) / sp$sigma_R[i]) - .5))
  }
  expect_equal(vapply(seq_along(mu), function(i) cdf(mu[i], i), numeric(1)),
               rep(.5, length(mu)), tolerance = 5e-15)
  expect_equal(sp$c[2], mu[2], tolerance = 0)
  expect_equal(sp$sigma_L * sp$sigma_R, sigma^2, tolerance = 1e-14)
})

testthat::test_that("split constructors expose canonical launch columns and suffixes", {
  specs <- list(
    BAwD = BAwD(drift_distribution = "splitlognormal"),
    BAwDp = BAwDp(drift_distribution = "splitlognormal"),
    BAwF = BAwF(drift_distribution = "splitlognormal"),
    BAwR = BAwR(drift_distribution = "splitlognormal"),
    BAwL = BAwL(drift_distribution = "splitlognormal"),
    BTAwL = BTAwL(drift_distribution = "splitlognormal"),
    BTAwLTransient = BTAwLTransient(drift_distribution = "splitlognormal"),
    BTAwLSustained = BTAwLSustained(drift_distribution = "splitlognormal")
  )
  expect_true(all(vapply(specs, function(x) grepl("_LOGN_SPLIT", x$c_name, fixed = TRUE), logical(1))))
  expect_true(all(vapply(specs, function(x) identical(unname(x$p_types_canonical[1:3]),
                                                       c("mu", "sigma", "delta")), logical(1))))
  expect_true(all(vapply(specs, function(x) identical(unname(x$transform$func["delta"]), "identity"), logical(1))))
  expect_true(all(vapply(specs, function(x) identical(unname(x$bound$minmax[ , "delta"]), c(-Inf, Inf)), logical(1))))
})

testthat::test_that("delta zero preserves the ordinary lognormal scalar path", {
  compiled <- c("dbawd", "pbawd", "dbawdp", "pbawdp", "dbawf", "pbawf",
                "dbawr", "pbawr", "dBTAwLTransient", "pBTAwLTransient",
                "dBTAwL", "pBTAwL")
  available <- all(vapply(compiled, function(x) exists(x, envir = asNamespace("EMC2"), inherits = FALSE), logical(1)))
  skip_if_not(available, "compiled split-launch bindings are unavailable")
  expect_equal(EMC2:::dbawd(1, 0, 1, 0, 1, 0, 1, launch = 1),
               EMC2:::dbawd(1, 0, 1, 0, 1, 0, 1, launch = 2, delta = 0), tolerance = 1e-12)
  expect_equal(EMC2:::pbawd(1, 0, 1, 0, 1, 0, 1, launch = 1),
               EMC2:::pbawd(1, 0, 1, 0, 1, 0, 1, launch = 2, delta = 0), tolerance = 1e-12)
  expect_equal(EMC2:::dbawf(1, 0, 1, 0, 1, 0.5, 2, launch = 1),
               EMC2:::dbawf(1, 0, 1, 0, 1, 0.5, 2, launch = 2, delta = 0), tolerance = 1e-12)
  expect_equal(EMC2:::pbawf(1, 0, 1, 0, 1, 0.5, 2, launch = 1),
               EMC2:::pbawf(1, 0, 1, 0, 1, 0.5, 2, launch = 2, delta = 0), tolerance = 1e-12)
  expect_equal(EMC2:::dbawr(1, 0, 1, 0, 1, 0.5, 1, launch = 1),
               EMC2:::dbawr(1, 0, 1, 0, 1, 0.5, 1, launch = 2, delta = 0), tolerance = 1e-12)
  expect_equal(EMC2:::pbawr(1, 0, 1, 0, 1, 0.5, 1, launch = 1),
               EMC2:::pbawr(1, 0, 1, 0, 1, 0.5, 1, launch = 2, delta = 0), tolerance = 1e-12)
})

testthat::test_that("delta reaches point-start kernels and the pure-R fallback sampler", {
  compiled <- c("pbawf", ".bawd_split_rlnorm")
  available <- all(vapply(compiled, function(x)
    exists(x, envir = asNamespace("EMC2"), inherits = FALSE), logical(1)))
  skip_if_not(available, "compiled split-launch bindings are unavailable")
  tv <- c(0.35, 0.6, 1.0)
  expect_equal(EMC2:::pbawf(tv, 0, 1.5, 0, 1, 2, launch = 2L),
               EMC2:::pbawf(tv, 0, 1.5, 0, 1, 2, launch = 1L), tolerance = 1e-12)
  expect_true(!isTRUE(all.equal(
    EMC2:::pbawf(tv, 0, 1.5, 0, 1, 2, launch = 2L, delta = 0.4),
    EMC2:::pbawf(tv, 0, 1.5, 0, 1, 2, launch = 2L))))
  set.seed(11)
  y <- EMC2:::.bawd_split_rlnorm(rep(0, 5000), rep(1, 5000), rep(0.5, 5000))
  expect_true(all(is.finite(y)))
  expect_lt(abs(median(y) - 1), 0.05)  # median(V) = exp(mu)
})

testthat::test_that("BTAwL split designs accept documented pi endpoints", {
  d <- BTAwL(drift_distribution = "splitlognormal")
  skip_if_not("pi" %in% colnames(d$bound$minmax), "full BTAwL design carries no pi")
  mm <- d$bound$minmax
  lo <- mm[1, ]; hi <- mm[2, ]
  mid <- ifelse(is.finite(lo) & is.finite(hi), (lo + hi) / 2,
                ifelse(is.finite(lo), pmin(lo + 0.5, 1),
                       pmax(pmin(1, hi - 0.5), 0)))
  mk <- function(pi_val) {
    p <- as.data.frame(as.list(mid))
    p[["pi"]] <- pi_val
    p[, colnames(mm), drop = FALSE]
  }
  expect_true(all(EMC2:::do_bound(mk(0), d$bound)))
  expect_true(all(EMC2:::do_bound(mk(1), d$bound)))
  b_no_pi <- d$bound
  b_no_pi$exception <- b_no_pi$exception[names(b_no_pi$exception) != "pi"]
  expect_false(all(EMC2:::do_bound(mk(0), b_no_pi)))
})

testthat::test_that("ballistic model constructors work with default arguments", {
  expect_silent(m_bawl <- BAwL())
  expect_silent(m_bawd <- BAwD())
  expect_silent(m_bawdp <- BAwDp())
  expect_silent(m_bawf <- BAwF())
  expect_silent(m_bawr <- BAwR())
  expect_silent(m_btawl <- BTAwL())
  expect_equal(m_bawl$type, "RACE")
  expect_equal(m_bawl$c_name, "BAwL")
  expect_equal(m_bawd$type, "RACE")
  expect_equal(m_bawd$c_name, "BAwD_LOGN")
  expect_equal(m_bawdp$type, "RACE")
  expect_equal(m_bawdp$c_name, "BAwDp_LOGN")
})

testthat::test_that("BAwF keeps the nuisance columns available to design", {
  expect_true(all(c("pContaminant", "pGuess") %in% names(BAwF()$p_types)))
  expect_silent(design(
    factors = list(subjects = 1, S = 1), Rlevels = 1,
    formula = list(mu ~ 1, sigma ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1,
                   pContaminant ~ 1),
    model = BAwF, report_p_vector = FALSE
  ))
})

testthat::test_that("split lognormal CDF is non-decreasing and well-behaved", {
  t_grid <- seq(0.2, 3.0, by = 0.2)
  for (del in c(-0.8, 0, 0.8)) {
    # BAwF CDF
    p_f <- EMC2:::pbawf(t_grid, 0, 1.2, 0.3, 1, 0.5, 2, launch = 2L, delta = del)
    expect_true(all(is.finite(p_f)))
    expect_true(all(diff(p_f) >= -1e-12))
    expect_true(all(p_f >= 0 & p_f <= 1))

    # BAwD CDF
    p_d <- EMC2:::pbawd(t_grid, 0, 1.2, 0.3, 1, 0.5, 2, launch = 2L, delta = del)
    expect_true(all(is.finite(p_d)))
    expect_true(all(diff(p_d) >= -1e-12))
    expect_true(all(p_d >= 0 & p_d <= 1))

    # BAwR CDF
    p_r <- EMC2:::pbawr(t_grid, 0, 1.2, 0.3, 1, 0.5, 1, launch = 2L, delta = del)
    expect_true(all(is.finite(p_r)))
    expect_true(all(diff(p_r) >= -1e-12))
    expect_true(all(p_r >= 0 & p_r <= 1))
  }
})

testthat::test_that("split launch wrappers accept row-wise delta vectors", {
  compiled <- c("dbawd", "pbawd", "dbawf", "pbawf", "dbawr", "pbawr",
                "dbtawl_transient", "pbtawl_transient",
                "dbtawl_local_race", "pbtawl_local_race",
                "btawl_transient_log_surv_vec", "btawl_local_race_log_surv_vec")
  skip_if_not(all(vapply(compiled, function(x)
    exists(x, envir = asNamespace("EMC2"), inherits = FALSE), logical(1))),
    "compiled split-launch bindings are unavailable")
  n <- 5
  t <- seq(0.2, 2, length.out = n)
  delta <- seq(-0.8, 0.8, length.out = n)
  A <- rep(0.3, n); b <- rep(1, n); mu <- rep(0, n); sigma <- rep(1, n)
  k <- rep(0.5, n)
  expect_length(EMC2:::dbawd(t, A, b, mu, sigma, k, rep(0.5, n),
                             launch = 2, delta = delta), n)
  expect_length(EMC2:::pbawd(t, A, b, mu, sigma, k, rep(0.5, n),
                             launch = 2, delta = delta), n)
  expect_length(EMC2:::dbawf(t, A, b, mu, sigma, k,
                             launch = 2, delta = delta), n)
  expect_length(EMC2:::pbawf(t, A, b, mu, sigma, k,
                             launch = 2, delta = delta), n)
  expect_length(EMC2:::dbawr(t, A, b, mu, sigma, rep(0.5, n), rep(1, n),
                             launch = 2, delta = delta), n)
  expect_length(EMC2:::pbawr(t, A, b, mu, sigma, rep(0.5, n), rep(1, n),
                             launch = 2, delta = delta), n)

  pars <- data.frame(mu = mu, sigma = sigma, delta = delta, b = b, A = A,
                     t0 = rep(0.1, n), k = k, tau = rep(1, n),
                     tau_s = rep(1, n), tau_t = rep(1, n), pi = rep(0.5, n))
  expect_length(EMC2:::pBTAwLTransient(t, pars, launch = 2), n)
  expect_length(EMC2:::pBTAwL(t, pars, launch = 2), n)
  expect_length(EMC2:::pBTAwLSustained(t, pars, launch = 2), n)
  expect_length(EMC2:::btawl_transient_log_surv_vec(
    t, A, b, mu, sigma, k, rep(1, n), launch = 2, delta = delta), n)
  expect_length(EMC2:::btawl_local_race_log_surv_vec(
    t, A, b, mu, sigma, k, rep(1, n), rep(1, n), rep(0.5, n),
    launch = 2, delta = delta), n)
})

testthat::test_that("large delta produces valid finite splitlognormal parameters without intermediate overflow", {
  skip_if_not(exists(".bawd_split_params", envir = asNamespace("EMC2"), inherits = FALSE),
             "split launch helpers are unavailable")
  split_params <- getFromNamespace(".bawd_split_params", "EMC2")
  
  # delta = -1417.5 causes 3*sigma_R and 4*sigma_R to individually overflow to Inf
  # during intermediate steps (since sigma_R = 6.4e307), but sigma_R itself is finite.
  # The simplified math avoids Inf/Inf = NaN.
  sp <- split_params(mu = 0, sigma = 1, delta = -1417.5)
  
  expect_true(is.finite(sp$sigma_R))
  expect_true(is.finite(sp$c))
  expect_true(!is.nan(sp$c))
})

testthat::test_that("delta zero preserves BAwL lognormal path", {
  skip_if_not(exists("dbawl", envir = asNamespace("EMC2"), inherits = FALSE),
              "dbawl unavailable")
  expect_equal(EMC2:::dbawl(1, 0, 1, 0, 1, 0.5, launch = 1),
               EMC2:::dbawl(1, 0, 1, 0, 1, 0.5, launch = 2, delta = 0), tolerance = 1e-12)
  expect_equal(EMC2:::pbawl(1, 0, 1, 0, 1, 0.5, launch = 1),
               EMC2:::pbawl(1, 0, 1, 0, 1, 0.5, launch = 2, delta = 0), tolerance = 1e-12)
})
