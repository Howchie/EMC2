# FRQ(rate_sharing = "shared"): one Gamma rate multiplier per trial.

frqsh_cols <- c("K", "v", "cv", "h", "t0", "plo", "clock")

expect_abs <- function(a, b, tol) expect_lt(max(abs(a - b)), tol)

frqsh_pars <- function(n, h = 1, plo = 0, clock = 0, cv = 0.4) {
  cbind(K = runif(n, 2, 7), v = runif(n, 3, 9), cv = cv, h = h,
        t0 = runif(n, 0.1, 0.3), plo = plo, clock = clock)
}

# Joint density of winner w at t by integration over the shared multiplier:
# given W = w each accumulator is an independent counter with clock w phi(u).
frqsh_ref_cause <- function(t, P, w_i, accel = FALSE) {
  n <- nrow(P)
  beta <- 1 / P[1, "cv"]^2
  g <- EMC2:::frq_rate(P[, "K"], P[, "v"], P[, "cv"] * P[, "v"], P[, "h"],
                       P[, "plo"], if (accel) rep(1, n) else P[, "clock"], rep(0, n))
  phi_d <- function(i, t) {
    u <- t - P[i, "t0"]
    if (u <= 0) return(c(0, 0))
    a <- if (accel) P[i, "clock"] else 0
    x <- if (a > 0) expm1(a * u) / a else u
    J <- if (a > 0) exp(a * u) else 1
    c2 <- if (accel) 1 else P[i, "clock"]^2
    lam <- g[i, "lambda"]; p <- g[i, "p"]
    if (c2 == 0) { G <- -expm1(-lam * x); gd <- lam * exp(-lam * x) }
    else { G <- 1 - (1 + c2 * lam * x)^(-1 / c2); gd <- lam * (1 + c2 * lam * x)^(-1 / c2 - 1) }
    q <- p * G
    c(q / (1 - q), p * gd / (1 - q)^2 * J)
  }
  sapply(t, function(tt) {
    pd <- sapply(seq_len(n), function(i) phi_d(i, tt))
    f <- function(w) sapply(w, function(ww) {
      out <- dgamma(ww, beta)
      for (i in seq_len(n)) {
        y <- ww * pd[1, i]; pl <- P[i, "plo"]
        out <- out * if (i == w_i)
          ww * pd[2, i] * ((1 - pl) * dgamma(y, P[i, "K"]) + pl * dexp(y))
        else (1 - pl) * pgamma(y, P[i, "K"], lower.tail = FALSE) + pl * exp(-y)
      }
      out
    })
    log(integrate(f, 0, Inf, rel.tol = 1e-12, subdivisions = 5000)$value)
  })
}

test_that("the shared-rate expectation matches a direct integral", {
  ex <- EMC2:::frq_shared_expect_cpp
  ref <- function(b, K, lpsi, plo) {
    ss <- seq(-60, log(b) + 4, length.out = 200001)
    out <- b * ss - exp(ss) - lgamma(b)
    for (l in seq_along(K)) {
      y <- exp(ss + lpsi[l])
      lq <- pgamma(y, K[l], lower.tail = FALSE, log.p = TRUE)
      out <- out + if (plo[l] > 0) log((1 - plo[l]) * exp(lq) + plo[l] * exp(-y)) else lq
    }
    m <- max(out); w <- exp(out - m)
    m + log(diff(ss[1:2]) * (sum(w) - (w[1] + w[length(w)]) / 2))
  }
  set.seed(2)
  for (i in 1:40) {
    m <- 1 + (i %% 3)
    b <- runif(1, 1.5, 60); K <- runif(m, 1, 9); lpsi <- runif(m, -4, 2.5)
    plo <- ifelse(runif(m) < 0.4, runif(m, 0, 0.3), 0)
    r <- ref(b, K, lpsi, plo)
    expect_abs(ex(b, K, lpsi, plo), r, 1e-8)
    if (all(plo == 0)) expect_abs(ex(b, K, lpsi, plo, 1L), r, 1e-8)
  }
  # One factor is the beta tail.
  expect_equal(ex(48.99, 7.29, 2.43, 0),
               pbeta(1 / (1 + exp(2.43)), 48.99, 7.29, log.p = TRUE), tolerance = 1e-12)
})

test_that("the joint density matches integration over the shared rate", {
  skip_on_cran()
  cause <- EMC2:::frq_shared_cause_cpp
  set.seed(3)
  tt <- c(0.35, 0.5, 0.8, 1.3, 2.5)
  cfgs <- list(list(n = 2), list(n = 2, plo = 0.2),
               list(n = 2, h = 0.9, clock = 0.5), list(n = 3, plo = 0.15, h = 0.95))
  for (cfg in cfgs) {
    P <- do.call(frqsh_pars, cfg)
    for (w in seq_len(cfg$n))
      expect_abs(cause(tt, P, w), frqsh_ref_cause(tt, P, w), 1e-9)
  }
  P <- frqsh_pars(2, cv = 0.5); P[, "clock"] <- c(0.8, 1.5)
  for (w in 1:2)
    expect_abs(cause(tt, P, w, TRUE), frqsh_ref_cause(tt, P, w, TRUE), 1e-9)
})

test_that("the joint survivor is the atom plus the integrated causes", {
  skip_on_cran()
  cause <- EMC2:::frq_shared_cause_cpp; surv <- EMC2:::frq_shared_surv_cpp
  set.seed(4)
  for (cfg in list(list(n = 2, h = 0.9, plo = 0.1), list(n = 3, h = 0.97),
                   list(n = 2, h = 0.9, cv = 1.3))) {
    P <- do.call(frqsh_pars, cfg); t1 <- 0.6
    tail <- sum(sapply(seq_len(cfg$n), function(w)
      integrate(function(t) exp(cause(t, P, w)), t1, Inf, rel.tol = 1e-11)$value))
    expect_equal(exp(surv(t1, P)), exp(surv(Inf, P)) + tail, tolerance = 1e-8)
  }
})

test_that("the shared simulator agrees with the joint density", {
  skip_on_cran()
  cause <- EMC2:::frq_shared_cause_cpp; surv <- EMC2:::frq_shared_surv_cpp
  set.seed(5)
  P <- frqsh_pars(3, h = 0.93, clock = 0.7, cv = 0.5)
  nsim <- 1e5
  lR <- factor(rep(c("a", "b", "c"), nsim), levels = c("a", "b", "c"))
  sim <- EMC2:::rFRQ_shared(lR, P[rep(1:3, nsim), ])
  for (w in 1:3) for (q in c(0.6, 1.4)) {
    pred <- integrate(function(t) exp(cause(t, P, w)), 0, q, rel.tol = 1e-10)$value
    obs <- mean(sim$R == levels(lR)[w] & sim$rt <= q, na.rm = TRUE)
    expect_lt(abs(obs - pred), 4 * sqrt(pred * (1 - pred) / nsim) + 1e-4)
  }
  expect_lt(abs(mean(is.na(sim$R)) - exp(surv(Inf, P))), 0.002)
})

frqsh_ll <- function(emc, p, pw = FALSE) {
  model <- emc[[1]]$model(); p_types <- names(model$p_types); dadm <- emc[[1]]$data[[1]]
  designs <- list()
  for (nm in p_types) designs[[nm]] <- attr(dadm, "designs")[[nm]][
    attr(attr(dadm, "designs")[[nm]], "expand"), , drop = FALSE]
  f <- if (pw) EMC2:::calc_ll_oo_pw else EMC2:::calc_ll_oo
  f(matrix(p, nrow = 1, dimnames = list(NULL, names(p))), dadm,
    constants = attr(dadm, "constants"), designs = designs, type = model$c_name,
    bounds = model$bound, transforms = model$transform,
    pretransforms = model$pre_transform, p_types = p_types,
    min_ll = log(1e-10), trend = model$trend)
}

frqsh_fixture <- function(n = 100, sharing = "shared", TC = NULL, dat = NULL,
                          constants = NULL) {
  matchfun <- function(d) d$S == d$lR
  if (is.null(dat)) {
    dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1], ]
    dat$subjects <- droplevels(dat$subjects)
    dat <- dat[seq_len(n), ]
  }
  ADmat <- matrix(c(-1 / 2, 1 / 2), ncol = 1, dimnames = list(NULL, "d"))
  formula <- list(K ~ 1, v ~ lM, h ~ lM, t0 ~ 1, plo ~ 1, clock ~ 1)
  if (is.null(constants)) formula <- c(formula, list(cv ~ 1))
  des <- design(data = dat, model = function() FRQ(parameterization = "cv",
                                                   rate_sharing = sharing),
                matchfun = matchfun, formula = formula, constants = constants,
                contrasts = list(h = list(lM = ADmat), v = list(lM = ADmat)), TC = TC)
  p <- sampled_pars(des); p[] <- 0
  p["K"] <- log(1.5); p["h"] <- qnorm(0.95); p["h_lMd"] <- 0.6
  p["v"] <- log(4); p["v_lMd"] <- 0.8; p["t0"] <- log(0.2)
  p["plo"] <- qnorm(0.1); p["clock"] <- log(0.5)
  if ("cv" %in% names(p)) p["cv"] <- log(0.6)
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
                                   compress = FALSE, rt_resolution = NULL))
  list(dat = dat, des = des, p = p, emc = emc)
}

frqsh_ref_ll <- function(e, p, LT = 0, UT = Inf, keep_atom = TRUE) {
  cause <- EMC2:::frq_shared_cause_cpp; surv <- EMC2:::frq_shared_surv_cpp
  dadm <- e[[1]]$data[[1]]; pars <- EMC2:::get_pars_matrix_oo(p, dadm, e[[1]]$model())
  n_lR <- nlevels(dadm$lR)
  tr <- rep(seq_len(nrow(dadm) / n_lR), each = n_lR)
  win <- as.logical(dadm$winner)
  sum(sapply(split(seq_len(nrow(dadm)), tr), function(ix) {
    P <- pars[ix, frqsh_cols, drop = FALSE]; rt <- dadm$rt[ix][1]
    num <- if (!is.finite(rt)) surv(Inf, P) else cause(rt, P, which(win[ix]))
    if (LT > 0 || is.finite(UT))
      num <- num - log(exp(surv(LT, P)) - exp(surv(UT, P)) + keep_atom * exp(surv(Inf, P)))
    max(num, log(1e-10))
  }))
}

test_that("the constructor exposes the shared contract", {
  m <- FRQ(parameterization = "cv", rate_sharing = "shared")
  expect_identical(m$c_name, "FRQ_CV_SHARED")
  expect_false("omega" %in% names(m$p_types))
  expect_identical(FRQ(parameterization = "accel", rate_sharing = "shared")$c_name,
                   "FRQ_ACCEL_SHARED")
  expect_error(FRQ(rate_sharing = "shared"), "parameterization")
})

test_that("the compiled shared likelihood matches the trial reference", {
  skip_on_cran()
  fx <- frqsh_fixture()
  expect_equal(frqsh_ll(fx$emc, fx$p), frqsh_ref_ll(fx$emc, fx$p), tolerance = 1e-8)
  expect_equal(sum(frqsh_ll(fx$emc, fx$p, TRUE)), frqsh_ll(fx$emc, fx$p), tolerance = 1e-10)

  LT <- 0.25; UT <- 1.2
  d2 <- fx$dat; d2$rt[1:5] <- Inf; d2$R[1:5] <- NA
  d2 <- d2[(d2$rt >= LT & d2$rt <= UT) | is.infinite(d2$rt), ]
  d2$LT <- LT; d2$UT <- UT
  e2 <- suppressMessages(make_emc(d2, fx$des, type = "single", n_chains = 1,
                                  compress = FALSE, rt_resolution = NULL))
  expect_equal(frqsh_ll(e2, fx$p), frqsh_ref_ll(e2, fx$p, LT, UT), tolerance = 1e-8)

  d3 <- d2[is.finite(d2$rt), ]
  f3 <- frqsh_fixture(dat = d3, TC = list(filter_defective = TRUE))
  expect_equal(frqsh_ll(f3$emc, fx$p), frqsh_ref_ll(f3$emc, fx$p, LT, UT, FALSE),
               tolerance = 1e-8)
})

test_that("cv = 0 is exactly the independent race", {
  skip_on_cran()
  fs <- frqsh_fixture(constants = c(cv = log(0)))
  fi <- frqsh_fixture(sharing = "independent", constants = c(cv = log(0)))
  expect_equal(frqsh_ll(fs$emc, fs$p), frqsh_ll(fi$emc, fs$p), tolerance = 1e-12)
})

test_that("cv must be shared by the accumulators of a trial", {
  skip_on_cran()
  fx <- frqsh_fixture()
  dat <- fx$dat
  matchfun <- function(d) d$S == d$lR
  des <- design(data = dat, model = function() FRQ(parameterization = "cv",
                                                   rate_sharing = "shared"),
                matchfun = matchfun, formula = list(K ~ 1, v ~ 1, cv ~ lM, h ~ 1, t0 ~ 1))
  e <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
                                 compress = FALSE, rt_resolution = NULL))
  p <- sampled_pars(des); p[] <- 0; p["t0"] <- log(0.2); p["cv"] <- log(0.5)
  p["cv_lMTRUE"] <- 0.3
  expect_error(frqsh_ll(e, p), "cv must be equal")
})
