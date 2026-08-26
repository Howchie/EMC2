# Distributional cross-checks between the C++ rfun kernels (src/model_rng.cpp)
# and the R reference rfuns they replace, per rfun_port_plan.md section 4.1.
# The two paths are distributionally, not stream-, equivalent -- draws differ
# trial by trial at a fixed seed -- so we compare response proportions,
# omission rates and winner-RT distributions (KS statistic) rather than exact
# values.

# tol = 0.05 keeps the false-positive rate low at n ~ 8000 (asymptotic 0.01
# critical value for two ~4000-sample groups is ~0.03) while still catching
# real distributional mismatches, which in testing produced KS > 0.15.
ks_ok <- function(a, b, tol = 0.05) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < 20 || length(b) < 20) return(TRUE)
  suppressWarnings(ks.test(a, b)$statistic) < tol
}

# Absolute-difference proportion comparison (expect_equal's default relative
# tolerance is too tight for small proportions like omission rates).
prop_close <- function(a, b, tol = 0.02) expect_lt(abs(a - b), tol)

n <- 8000

test_that("LBA cpp kernel matches the R LBA reference distributionally", {
  skip_on_cran()
  lR <- factor(rep(c("left", "right"), n), levels = c("left", "right"))
  pars <- cbind(v = rep(c(1, 0.3), n), sv = 1, b = 1.5, A = 0.3, t0 = 0.2)
  ok <- rep(TRUE, nrow(pars))

  set.seed(1); r <- EMC2:::.lba_rfun(lR, pars, ok = ok, posdrift = TRUE)
  set.seed(2); cpp <- EMC2:::rlba_cpp(pars, levels(lR), ok, TRUE)

  prop_close(mean(r$R == "left", na.rm = TRUE), mean(cpp$R == 1, na.rm = TRUE))
  prop_close(mean(is.na(r$R)), mean(is.na(cpp$R)))
  expect_true(ks_ok(r$rt[r$R == "left"], cpp$rt[cpp$R == 1]))
})

test_that("LBA cpp kernel matches the R LBA reference for posdrift = FALSE (defective/LBAIO)", {
  skip_on_cran()
  lR <- factor(rep(c("left", "right"), n), levels = c("left", "right"))
  pars <- cbind(v = rep(c(-0.5, 0.3), n), sv = 1, b = 1.5, A = 0.3, t0 = 0.2)
  ok <- rep(TRUE, nrow(pars))

  set.seed(1); r <- EMC2:::.lba_rfun(lR, pars, ok = ok, posdrift = FALSE)
  set.seed(2); cpp <- EMC2:::rlba_cpp(pars, levels(lR), ok, FALSE)

  prop_close(mean(is.na(r$R)), mean(is.na(cpp$R)))
  expect_true(ks_ok(r$rt[r$R == "right"], cpp$rt[cpp$R == 2]))
})

test_that("RDM cpp kernel matches R rRDM distributionally", {
  skip_on_cran()
  lR <- factor(rep(c("left", "right"), n), levels = c("left", "right"))
  pars <- cbind(v = rep(c(1, 0.3), n), B = 1, A = 0.3, t0 = 0.2, s = 1)
  ok <- rep(TRUE, nrow(pars))

  set.seed(1); r <- EMC2:::rRDM(lR, pars, ok = ok)
  set.seed(2); cpp <- EMC2:::rrdm_cpp(pars, levels(lR), ok)

  prop_close(mean(r$R == "left", na.rm = TRUE), mean(cpp$R == 1, na.rm = TRUE))
  expect_true(ks_ok(r$rt[r$R == "left"], cpp$rt[cpp$R == 1]))
})

test_that("RDM no-finish trials are coded as omissions with rt = Inf", {
  skip_on_cran()

  n_edge <- 1000
  lR <- factor(rep(c("left", "right"), n_edge), levels = c("left", "right"))
  ok <- rep(TRUE, length(lR))
  pars <- cbind(
    v = c(rep(0, n_edge), rep(c(1, 0.3), length.out = n_edge)),
    B = 1, A = 0.3, t0 = 0.2, s = 1
  )

  set.seed(101)
  r <- EMC2:::rRDM(lR, pars, ok = ok)
  set.seed(102)
  cpp <- EMC2:::rrdm_cpp(pars, levels(lR), ok)

  expect_true(any(is.na(r$R)))
  expect_true(any(is.na(cpp$R)))
  expect_true(all(is.infinite(r$rt[is.na(r$R)])))
  expect_true(all(is.na(r$R[is.infinite(r$rt)])))
  expect_true(all(is.infinite(cpp$rt[is.na(cpp$R)])))
  expect_true(all(is.na(cpp$R[is.infinite(cpp$rt)])))
})

test_that("BAwL cpp kernel matches R rBAwL: local guess+kill, shape 1", {
  skip_on_cran()
  lR <- factor(rep(c("left", "right"), n), levels = c("left", "right"))
  pars <- cbind(v = rep(c(1, 0.3), n), sv = 1, b = 1.5, A = 0.3, t0 = 0.2,
               k = 0.5, lambda_g = 0.3, lambda_k = 0.2)
  ok <- rep(TRUE, nrow(pars))

  set.seed(1); r <- EMC2:::rBAwL(lR, pars, ok = ok, posdrift = TRUE, erlang = 1L, guess = TRUE, global = FALSE)
  set.seed(2); cpp <- EMC2:::rbawl_cpp(pars, levels(lR), ok, TRUE, 1L, TRUE, FALSE)

  prop_close(mean(is.na(r$R)), mean(is.na(cpp$R)))
  expect_true(ks_ok(r$rt[r$R == "left"], cpp$rt[cpp$R == 1]))
})

test_that("BAwL cpp kernel matches R rBAwL: global kill, shape 2 (E2)", {
  skip_on_cran()
  lR <- factor(rep(c("left", "right"), n), levels = c("left", "right"))
  pars <- cbind(v = rep(c(1, 0.3), n), sv = 1, b = 1.5, A = 0.3, t0 = 0.2,
               k = 0.5, lambda_g = 0, lambda_k = 0.3)
  ok <- rep(TRUE, nrow(pars))

  set.seed(1); r <- EMC2:::rBAwL(lR, pars, ok = ok, posdrift = TRUE, erlang = 2L, guess = FALSE, global = TRUE)
  set.seed(2); cpp <- EMC2:::rbawl_cpp(pars, levels(lR), ok, TRUE, 2L, FALSE, TRUE)

  prop_close(mean(is.na(r$R)), mean(is.na(cpp$R)))
  expect_true(ks_ok(r$rt[r$R == "left"], cpp$rt[cpp$R == 1]))
})

test_that("cpp race kernels code omissions as rt = Inf, matching R rfuns", {
  skip_on_cran()

  n_edge <- 1000
  lR <- factor(rep(c("left", "right"), n_edge), levels = c("left", "right"))
  ok <- rep(TRUE, length(lR))

  pars_lba <- cbind(
    v = c(rep(-1, n_edge), rep(c(1, 0.5), length.out = n_edge)),
    sv = 1e-8, b = 1.5, A = 0.3, t0 = 0.2
  )
  set.seed(201)
    r_lba <- EMC2:::.lba_rfun(lR, pars_lba, ok = ok, posdrift = FALSE)
  set.seed(202)
  cpp_lba <- EMC2:::rlba_cpp(pars_lba, levels(lR), ok, FALSE)
  expect_true(any(is.na(r_lba$R)))
  expect_true(any(is.na(cpp_lba$R)))
  expect_true(all(is.infinite(r_lba$rt[is.na(r_lba$R)])))
  expect_true(all(is.infinite(cpp_lba$rt[is.na(cpp_lba$R)])))

  pars_bawl <- cbind(
    v = rep(1, length(lR)), sv = 1e-8, b = 1.5, A = 0.3, t0 = 0.2,
    k = 0.5, lambda_g = 0, lambda_k = 100
  )
  set.seed(203)
  r_bawl <- EMC2:::rBAwL(lR, pars_bawl, ok = ok, posdrift = TRUE,
                         erlang = 1L, guess = FALSE, global = FALSE)
  set.seed(204)
  cpp_bawl <- EMC2:::rbawl_cpp(pars_bawl, levels(lR), ok, TRUE, 1L, FALSE, FALSE)
  expect_true(any(is.na(r_bawl$R)))
  expect_true(any(is.na(cpp_bawl$R)))
  expect_true(all(is.infinite(r_bawl$rt[is.na(r_bawl$R)])))
  expect_true(all(is.infinite(cpp_bawl$rt[is.na(cpp_bawl$R)])))

  pars_rdmswtn <- cbind(
    v = rep(1, length(lR)), b = 1.5, A = 0.3, t0 = 0.2, sv = 1e-8,
    lambda_g = 0, lambda_k = 100, s = 1
  )
  set.seed(205)
  r_rdmswtn <- EMC2:::rRDMSWTN(lR, pars_rdmswtn, ok = ok, erlang_shape = 1L,
                               erlang_type = "local_kill", posdrift = TRUE)
  set.seed(206)
  cpp_rdmswtn <- EMC2:::rrdmswtn_cpp(pars_rdmswtn, levels(lR), ok, 1L, "local_kill", TRUE)
  expect_true(any(is.na(r_rdmswtn$R)))
  expect_true(any(is.na(cpp_rdmswtn$R)))
  expect_true(all(is.infinite(r_rdmswtn$rt[is.na(r_rdmswtn$R)])))
  expect_true(all(is.infinite(cpp_rdmswtn$rt[is.na(cpp_rdmswtn$R)])))
})

test_that("BAwL cpp kernel matches R rBAwL: k -> 0 reduces to LBA formula", {
  skip_on_cran()
  lR <- factor(rep(c("left", "right"), n), levels = c("left", "right"))
  pars_bawl <- cbind(v = rep(c(1, 0.3), n), sv = 1, b = 1.5, A = 0.3, t0 = 0.2,
                     k = 1e-12, lambda_g = 0, lambda_k = 0)
  pars_lba <- cbind(v = rep(c(1, 0.3), n), sv = 1, b = 1.5, A = 0.3, t0 = 0.2)
  ok <- rep(TRUE, n * 2)

  set.seed(3); cpp_bawl <- EMC2:::rbawl_cpp(pars_bawl, levels(lR), ok, TRUE, 1L, FALSE, FALSE)
  set.seed(4); cpp_lba <- EMC2:::rlba_cpp(pars_lba, levels(lR), ok, TRUE)

  prop_close(mean(cpp_bawl$R == 1, na.rm = TRUE), mean(cpp_lba$R == 1, na.rm = TRUE))
  expect_true(ks_ok(cpp_bawl$rt[cpp_bawl$R == 1], cpp_lba$rt[cpp_lba$R == 1]))
})

test_that("BTAwL cpp simulator covers transient, sustained and local-race launches", {
  skip_on_cran()
  n_btawl <- 1200L
  lR <- factor(rep(c("left", "right"), n_btawl), levels = c("left", "right"))
  ok <- rep(TRUE, length(lR))

  for (launch in 0:2) {
    launch_cols <- if (launch == 0L) {
      cbind(v = rep(c(3.0, 2.5), n_btawl), sv = 0.35)
    } else if (launch == 1L) {
      cbind(mu = rep(c(log(3.0), log(2.5)), n_btawl), sigma = 0.25)
    } else {
      cbind(mu = rep(c(log(3.0), log(2.5)), n_btawl), sigma = 0.25,
            delta = 0.4)
    }
    pars_t <- cbind(launch_cols, b = 1.2, A = 0.2, t0 = 0.1, k = 0.5,
                    tau = 0.8)
    pars_s <- cbind(launch_cols, b = 1.2, A = 0.2, t0 = 0.1, k = 0.5,
                    tau_s = 0.9)
    pars_f <- cbind(launch_cols, b = 1.2, A = 0.2, t0 = 0.1, k = 0.5,
                    tau_s = 0.9, tau_t = 0.8, pi = 0.4)

    mode_pars_list <- list(pars_t, pars_s, pars_f)
    for (mode in 0:2) {
      mode_pars <- mode_pars_list[[mode + 1L]]
      set.seed(100 + 10 * launch + mode)
      out <- EMC2:::rbta_wl_cpp(as.matrix(mode_pars), levels(lR), ok,
                                mode, TRUE, launch)
      expect_equal(length(out$R), n_btawl)
      expect_equal(length(out$rt), n_btawl)
      expect_true(any(is.finite(out$rt)))
      expect_true(all(is.na(out$R) == is.infinite(out$rt)))
    }
  }

  pars_tau <- cbind(v = rep(3, length(lR)), sv = 0.35, b = 1.2, A = 0.2,
                    t0 = 0.1, k = 0.5, tau = 2.0)
  set.seed(901)
  tau_out <- EMC2:::rbta_wl_cpp(pars_tau, levels(lR), ok, 0L, TRUE, 0L)
  expect_true(any(is.finite(tau_out$rt)))
})

test_that("BTAwL constructors dispatch their default rfun to C++", {
  skip_on_cran()
  old_opt <- getOption("emc2.cpp_rfun")
  on.exit(options(emc2.cpp_rfun = old_opt), add = TRUE)
  options(emc2.cpp_rfun = TRUE)

  lR <- factor(rep(c("left", "right"), 80), levels = c("left", "right"))
  cases <- list(
    BTAwLTransient(drift_distribution = "normal"),
    BTAwLSustained(drift_distribution = "normal"),
    BTAwL(drift_distribution = "normal")
  )
  raw <- list(
    cbind(v = rep(c(3, 2.5), 80), sv = .35, B = 1.0, A = .2,
          t0 = .1, k = .5, tau = 2),
    cbind(v = rep(c(3, 2.5), 80), sv = .35, B = 1.0, A = .2,
          t0 = .1, k = .5, tau_s = .9),
    cbind(v = rep(c(3, 2.5), 80), sv = .35, B = 1.0, A = .2,
          t0 = .1, k = .5, tau_s = .9, tau_t = .8, pi = .4)
  )

  for (i in seq_along(cases)) {
    pars <- cases[[i]]$Ttransform(raw[[i]], NULL)
    attr(pars, "ok") <- rep(TRUE, nrow(pars))
    out <- cases[[i]]$rfun(list(lR = lR), pars)
    expect_equal(nrow(out), 80L)
    expect_equal(length(out$R), length(out$rt))
    expect_true(all(is.na(out$R) == is.infinite(out$rt)))
  }

  pars <- cases[[3]]$Ttransform(raw[[3]], NULL)
  out <- EMC2:::rBTAwL(lR, pars, ok = rep(TRUE, nrow(pars)), launch = 0L)
  expect_equal(nrow(out), 80L)
})

test_that("BTAwL C++ draws match the R reference distributions", {
  skip_on_cran()
  n_btawl <- 3000L
  lR <- factor(rep(c("left", "right"), n_btawl), levels = c("left", "right"))
  ok <- rep(TRUE, length(lR))
  cases <- list(
    list(mode = 0L, ref = EMC2:::.rBTAwLTransient_R,
         pars = cbind(v = rep(c(5, 4), n_btawl), sv = .35, b = 1.2,
                      A = .2, t0 = .1, k = .5, tau = .8)),
    list(mode = 1L, ref = EMC2:::.rBTAwLSustained_R,
         pars = cbind(v = rep(c(5, 4), n_btawl), sv = .35, b = 1.2,
                      A = .2, t0 = .1, k = .5, tau_s = .9)),
    list(mode = 2L, ref = EMC2:::.rBTAwL_R,
         pars = cbind(v = rep(c(5, 4), n_btawl), sv = .35, b = 1.2,
                      A = .2, t0 = .1, k = .5, tau_s = .9,
                      tau_t = .8, pi = .4))
  )

  for (case in cases) {
    set.seed(1000 + case$mode)
    r <- case$ref(lR, case$pars, ok = ok, posdrift = TRUE, launch = 0L)
    set.seed(2000 + case$mode)
    cpp <- EMC2:::rbta_wl_cpp(case$pars, levels(lR), ok, case$mode,
                              TRUE, 0L)
    prop_close(mean(is.na(r$R)), mean(is.na(cpp$R)))
    prop_close(mean(r$R == "left", na.rm = TRUE),
               mean(cpp$R == 1, na.rm = TRUE))
    expect_true(ks_ok(r$rt[r$R == "left"], cpp$rt[cpp$R == 1]))
  }
})

test_that("RDMSWTN cpp kernel matches R rRDMSWTN: none", {
  skip_on_cran()
  lR <- factor(rep(c("left", "right"), n), levels = c("left", "right"))
  pars <- cbind(v = rep(c(1, 0.3), n), b = 1.5, A = 0.3, t0 = 0.2, sv = 0.3,
               lambda_g = 0, lambda_k = 0, s = 1)
  ok <- rep(TRUE, nrow(pars))

  set.seed(1); r <- EMC2:::rRDMSWTN(lR, pars, ok = ok, erlang_shape = 1L, erlang_type = "none", posdrift = TRUE)
  set.seed(2); cpp <- EMC2:::rrdmswtn_cpp(pars, levels(lR), ok, 1L, "none", TRUE)

  prop_close(mean(r$R == "left", na.rm = TRUE), mean(cpp$R == 1, na.rm = TRUE))
  expect_true(ks_ok(r$rt[r$R == "left"], cpp$rt[cpp$R == 1]))
})

test_that("RDMSWTN cpp kernel matches R rRDMSWTN: local_kill_guess, mixed Erlang", {
  skip_on_cran()
  lR <- factor(rep(c("left", "right"), n), levels = c("left", "right"))
  pars <- cbind(v = rep(c(1, 0.3), n), b = 1.5, A = 0.3, t0 = 0.2, sv = 0.3,
               lambda_g = 0.3, lambda_k = 0.2, s = 1, omega = 0.6)
  ok <- rep(TRUE, nrow(pars))

  set.seed(1); r <- EMC2:::rRDMSWTN(lR, pars, ok = ok, erlang_shape = 3L,
                                     erlang_type = "local_kill_guess", posdrift = TRUE)
  set.seed(2); cpp <- EMC2:::rrdmswtn_cpp(pars, levels(lR), ok, 3L, "local_kill_guess", TRUE)

  prop_close(mean(is.na(r$R)), mean(is.na(cpp$R)))
  expect_true(ks_ok(r$rt[r$R == "left"], cpp$rt[cpp$R == 1]))
  prop_close(mean(r$isTime, na.rm = TRUE), mean(cpp$isTime, na.rm = TRUE), tol = 0.05)
})

test_that("cpp kernels resample the time-level winner like .apply_timed_guess_winner", {
  skip_on_cran()
  lR <- factor(rep(c("left", "right", "time"), n), levels = c("left", "right", "time"))
  pars <- cbind(v = rep(c(1, 0.5, 3), n), sv = 1, b = 1.5, A = 0.3, t0 = 0.2)
  ok <- rep(TRUE, nrow(pars))

  set.seed(301)
  cpp <- EMC2:::rlba_cpp(pars, levels(lR), ok, TRUE)
  expect_false(any(cpp$R == 3, na.rm = TRUE))  # "time" level never a final response
  expect_true(any(cpp$isTime))
})

test_that("emc2.cpp_rfun option gates the fast path via make_data()", {
  skip_on_cran()
  matchfun <- function(d) d$S == d$lR
  design_lba <- design(data = forstmann, model = LBA, matchfun = matchfun,
                       formula = list(v ~ lM, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
                       contrasts = list(v = list(lM = matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d")))),
                       constants = c(sv = log(1)))
  p_vector <- sampled_pars(design_lba)
  p_vector[] <- 0
  p_vector["v"] <- 1; p_vector["v_lMd"] <- 1.5; p_vector["t0"] <- log(0.2)

  old_opt <- getOption("emc2.cpp_rfun")
  on.exit(options(emc2.cpp_rfun = old_opt), add = TRUE)

  options(emc2.cpp_rfun = TRUE)
  set.seed(1)
  dat_cpp <- make_data(p_vector, design = design_lba, model = LBA, n_trials = 1000)

  options(emc2.cpp_rfun = FALSE)
  set.seed(1)
  dat_r <- make_data(p_vector, design = design_lba, model = LBA, n_trials = 1000)

  expect_equal(nrow(dat_cpp), nrow(dat_r))
  expect_equal(mean(dat_cpp$R == "left"), mean(dat_r$R == "left"), tolerance = 0.05)
})
