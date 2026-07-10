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

test_that("LBA cpp kernel matches R rLBA distributionally", {
  skip_on_cran()
  lR <- factor(rep(c("left", "right"), n), levels = c("left", "right"))
  pars <- cbind(v = rep(c(1, 0.3), n), sv = 1, b = 1.5, A = 0.3, t0 = 0.2)
  ok <- rep(TRUE, nrow(pars))

  set.seed(1); r <- EMC2:::rLBA(lR, pars, ok = ok, posdrift = TRUE)
  set.seed(2); cpp <- EMC2:::rlba_cpp(pars, levels(lR), ok, TRUE)

  prop_close(mean(r$R == "left", na.rm = TRUE), mean(cpp$R == 1, na.rm = TRUE))
  prop_close(mean(is.na(r$R)), mean(is.na(cpp$R)))
  expect_true(ks_ok(r$rt[r$R == "left"], cpp$rt[cpp$R == 1]))
})

test_that("LBA cpp kernel matches R rLBA for posdrift = FALSE (defective/LBAIO)", {
  skip_on_cran()
  lR <- factor(rep(c("left", "right"), n), levels = c("left", "right"))
  pars <- cbind(v = rep(c(-0.5, 0.3), n), sv = 1, b = 1.5, A = 0.3, t0 = 0.2)
  ok <- rep(TRUE, nrow(pars))

  set.seed(1); r <- EMC2:::rLBA(lR, pars, ok = ok, posdrift = FALSE)
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

  lR <- factor(rep(c("left", "right"), 4), levels = c("left", "right"))
  ok <- rep(TRUE, length(lR))
  pars <- cbind(
    v = c(0, 0, 0, 0, 1, 0.3, 1, 0.3),
    B = 1, A = 0.3, t0 = 0.2, s = 1
  )

  r <- EMC2:::rRDM(lR, pars, ok = ok)
  cpp <- EMC2:::rrdm_cpp(pars, levels(lR), ok)

  expect_true(any(is.infinite(r$rt)))
  expect_true(all(is.na(r$R[is.infinite(r$rt)])))
  expect_true(all(is.infinite(cpp$rt[is.infinite(r$rt)])))
  expect_true(all(is.na(cpp$R[is.infinite(r$rt)])))
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

  lR <- factor(rep(c("left", "right"), 4), levels = c("left", "right"))
  ok <- rep(TRUE, length(lR))

  pars_lba <- cbind(
    v = c(-1, -1, -1, -1, 1, 0.5, 1, 0.5),
    sv = 1e-8, b = 1.5, A = 0.3, t0 = 0.2
  )
  r_lba <- EMC2:::rLBA(lR, pars_lba, ok = ok, posdrift = FALSE)
  cpp_lba <- EMC2:::rlba_cpp(pars_lba, levels(lR), ok, FALSE)
  expect_true(any(is.na(r_lba$R)))
  expect_true(all(is.infinite(r_lba$rt[is.na(r_lba$R)])))
  expect_true(all(is.infinite(cpp_lba$rt[is.na(cpp_lba$R)])))

  pars_bawl <- cbind(
    v = rep(1, length(lR)), sv = 1e-8, b = 1.5, A = 0.3, t0 = 0.2,
    k = 0.5, lambda_g = 0, lambda_k = 100
  )
  r_bawl <- EMC2:::rBAwL(lR, pars_bawl, ok = ok, posdrift = TRUE,
                         erlang = 1L, guess = FALSE, global = FALSE)
  cpp_bawl <- EMC2:::rbawl_cpp(pars_bawl, levels(lR), ok, TRUE, 1L, FALSE, FALSE)
  expect_true(any(is.na(r_bawl$R)))
  expect_true(all(is.infinite(r_bawl$rt[is.na(r_bawl$R)])))
  expect_true(all(is.infinite(cpp_bawl$rt[is.na(cpp_bawl$R)])))

  pars_rdmswtn <- cbind(
    v = rep(1, length(lR)), b = 1.5, A = 0.3, t0 = 0.2, sv = 1e-8,
    lambda_g = 0, lambda_k = 100, s = 1
  )
  r_rdmswtn <- EMC2:::rRDMSWTN(lR, pars_rdmswtn, ok = ok, erlang_shape = 1L,
                               erlang_type = "local_kill", posdrift = TRUE)
  cpp_rdmswtn <- EMC2:::rrdmswtn_cpp(pars_rdmswtn, levels(lR), ok, 1L, "local_kill", TRUE)
  expect_true(any(is.na(r_rdmswtn$R)))
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
  dat_cpp <- make_data(p_vector, design = design_lba, model = LBA, n_trials = 20)

  options(emc2.cpp_rfun = FALSE)
  set.seed(1)
  dat_r <- make_data(p_vector, design = design_lba, model = LBA, n_trials = 20)

  expect_equal(nrow(dat_cpp), nrow(dat_r))
  expect_equal(mean(dat_cpp$R == "left"), mean(dat_r$R == "left"), tolerance = 0.15)
})
