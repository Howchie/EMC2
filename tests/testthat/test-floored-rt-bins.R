# Floored response times: design_model() stores a point inside each
# rt_resolution bin and the race likelihood averages the density over the bin
# (src/rt_bins.h: midpoint away from onset, integration near onset).

bin_design <- function(model, formula) {
  design(factors = list(subjects = 1, S = c("left", "right")),
         Rlevels = c("left", "right"),
         matchfun = function(d) d$S == d$lR,
         formula = formula, model = model, verbose = FALSE)
}

bin_dadm <- function(dat, des, rt_resolution, compress = FALSE) {
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1,
                                   rt_resolution = rt_resolution, compress = compress))
  emc[[1]]$data[[1]]
}

bin_ll <- function(p, dadm, des) {
  as.numeric(EMC2:::calc_ll_manager(matrix(p, 1, dimnames = list(NULL, names(p))),
                                    dadm, des$model))
}

gauss_legendre <- function(n) {
  b <- seq_len(n - 1) / sqrt(4 * seq_len(n - 1)^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(1:(n - 1), 2:n)] <- b
  J[cbind(2:n, 1:(n - 1))] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = 2 * e$vectors[1, ]^2)
}

# log of the bin-averaged race density of each trial, integrated in
# z = sqrt(t - t0) over the part of the bin above t0.
bin_reference <- function(dat, des, p, h, n_nodes = 24) {
  g <- gauss_legendre(n_nodes)
  t0 <- exp(p[["t0"]])
  vapply(seq_len(nrow(dat)), function(i) {
    r <- h * floor(dat$rt[i] / h + 1e-9)
    lo <- max(r, t0); hi <- r + h
    zl <- sqrt(lo - t0); zh <- sqrt(hi - t0)
    z <- (zh - zl) / 2 * g$x + (zh + zl) / 2
    w <- (zh - zl) / 2 * g$w
    d <- bin_dadm(dat[i, ], des, NULL)
    f <- vapply(t0 + z^2, function(t) { d$rt[] <- t; exp(bin_ll(p, d, des)) }, 0)
    log(sum(w * 2 * z * f) / h)
  }, 0)
}

rdm_case <- list(model = RDM, formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1),
                 p = c(v = log(2), v_lMTRUE = log(1.6), B = 0, A = log(0.3),
                       t0 = log(0.25), s = 0))
lba_case <- list(model = LBA, formula = list(v ~ lM, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
                 p = c(v = 1, v_lMTRUE = 0.8, sv = 0, B = log(0.5), A = log(0.4),
                       t0 = log(0.25)))
lnr_case <- list(model = LNR, formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                 p = c(m = -0.8, m_lMTRUE = -0.4, s = log(0.4), t0 = log(0.2)))

test_that("design_model stores bin midpoints and records the bin width", {
  des <- bin_design(rdm_case$model, rdm_case$formula)
  set.seed(1)
  dat <- make_data(rdm_case$p, des, n_trials = 20, rt_resolution = NULL)
  h <- 0.05
  dadm <- bin_dadm(dat, des, h)
  expect_equal(attr(dadm, "rt_resolution"), h)
  q <- dadm$rt / h
  expect_true(all(abs(q - floor(q) - 0.5) < 1e-9))
  # Re-binning a binned rt is a no-op, so simulated data floored by make_data
  # land on the same midpoints.
  expect_equal(EMC2:::.rt_bin_midpoint(dadm$rt, h), dadm$rt)
  expect_null(attr(bin_dadm(dat, des, NULL), "rt_resolution"))
})

test_that("bin midpoints stay inside the censoring and truncation bounds", {
  x <- c(0.40, 1.47, 1.49)
  mid <- EMC2:::.rt_bin_midpoint(x, 0.05, lower = 0.42, upper = 1.48)
  expect_equal(mid, c((0.42 + 0.45) / 2, (1.45 + 1.48) / 2, (1.45 + 1.48) / 2))
  expect_true(all(mid >= 0.42 & mid < 1.48))
})

test_that("the binned race likelihood matches the bin-averaged density", {
  h <- 1 / 60
  for (case in list(rdm_case, lba_case, lnr_case)) {
    des <- bin_design(case$model, case$formula)
    set.seed(2)
    dat <- make_data(case$p, des, n_trials = 6, rt_resolution = NULL)
    ref <- bin_reference(dat, des, case$p, h)
    got <- vapply(seq_len(nrow(dat)), function(i) bin_ll(case$p, bin_dadm(dat[i, ], des, h), des), 0)
    edge <- vapply(seq_len(nrow(dat)), function(i)
      bin_ll(case$p, bin_dadm(transform(dat[i, ], rt = h * floor(rt / h + 1e-9)), des, NULL), des), 0)
    expect_lt(max(abs(got - ref)), 0.02)
    expect_lt(sum(abs(got - ref)), sum(abs(edge - ref)) / 5)
  }
})

test_that("t0 inside the earliest bin softens the likelihood instead of walling it", {
  h <- 0.05
  des <- bin_design(lba_case$model, lba_case$formula)
  set.seed(3)
  dat <- make_data(lba_case$p, des, n_trials = 10, rt_resolution = NULL)
  dadm <- bin_dadm(dat, des, h)
  lo <- h * floor(min(dat$rt) / h + 1e-9)
  ll_at <- function(t0) { p <- lba_case$p; p[["t0"]] <- log(t0); bin_ll(p, dadm, des) }
  # Continuous where t0 crosses the bin's lower edge.
  expect_equal(ll_at(lo - 1e-10), ll_at(lo + 1e-10), tolerance = 1e-6)
  # Above the stored midpoint the earliest trial still has probability, and
  # the likelihood falls as less of its bin remains feasible.
  above <- vapply(lo + c(0.6, 0.75, 0.9) * h, ll_at, 0)
  expect_true(all(is.finite(above)))
  expect_true(all(diff(above) < 0))
})

test_that("compression and the censored-data path apply the same bin rule", {
  h <- 1 / 60
  des <- bin_design(rdm_case$model, rdm_case$formula)
  set.seed(4)
  dat <- make_data(rdm_case$p, des, n_trials = 40, rt_resolution = NULL)
  p <- rdm_case$p
  p[["t0"]] <- log(h * floor(min(dat$rt) / h + 1e-9) + 0.4 * h)
  expect_equal(bin_ll(p, bin_dadm(dat, des, h, compress = TRUE), des),
               bin_ll(p, bin_dadm(dat, des, h), des), tolerance = 1e-10)
  # A censored trial sends the data down the general race path; its own
  # contribution does not depend on binning.
  cens <- dat[1, ]
  cens$rt <- Inf
  both <- rbind(dat, cens)
  both$trials <- seq_len(nrow(both))
  expect_equal(bin_ll(p, bin_dadm(both, des, h), des) - bin_ll(p, bin_dadm(dat, des, h), des),
               bin_ll(p, bin_dadm(both, des, NULL), des) - bin_ll(p, bin_dadm(dat, des, NULL), des),
               tolerance = 1e-8)
})

test_that("every subject's data keeps the bin width", {
  des <- design(factors = list(subjects = 1:2, S = c("left", "right")),
                Rlevels = c("left", "right"), matchfun = function(d) d$S == d$lR,
                formula = rdm_case$formula, model = RDM, verbose = FALSE)
  set.seed(5)
  dat <- make_data(rdm_case$p, des, n_trials = 10)
  emc <- suppressMessages(make_emc(dat, des, n_chains = 1, rt_resolution = 1 / 60))
  for (d in emc[[1]]$data) expect_equal(attr(d, "rt_resolution"), 1 / 60)
})

bin_trial <- function(rt = 0.44, response = "left") {
  data.frame(subjects = factor(1),
             S = factor("left", levels = c("left", "right")),
             R = factor(response, levels = c("left", "right")), rt = rt)
}

# Independent analytical reference: no likelihood floors inside the integral.
lnr_bin_ref <- function(p, lo, hi, h, response = "left") {
  t0 <- exp(p[["t0"]]); s <- exp(p[["s"]])
  m <- p[["m"]] + c(left = p[["m_lMTRUE"]], right = 0)
  f <- function(t) {
    x <- t - t0
    d <- vapply(c("left", "right"), function(r)
      dlnorm(x, m[[r]], s) * plnorm(x, m[[setdiff(names(m), r)]], s,
                                   lower.tail = FALSE), numeric(length(t)))
    if (is.na(response)) rowSums(d) else d[, response]
  }
  if (hi <= max(lo, t0)) return(log(1e-10))
  max(log(1e-10), log(integrate(f, max(lo, t0), hi, rel.tol = 1e-9)$value / h))
}

test_that("boundary bins agree with unfloored analytical densities", {
  des <- bin_design(LNR, lnr_case$formula)
  for (h in c(1/60, 0.05)) {
    dat <- bin_trial()
    lo <- h * floor(dat$rt / h); hi <- lo + h
    for (t0 in c(lo + 0.1 * h, lo + 0.5 * h, lo + 0.9 * h, hi)) {
      p <- lnr_case$p; p[["t0"]] <- log(t0)
      expect_equal(bin_ll(p, bin_dadm(dat, des, h), des),
                   lnr_bin_ref(p, lo, hi, h), tolerance = 1e-5)
    }
  }
})

test_that("clipped bins weight both process and uniform-guess densities", {
  des <- bin_design(LNR, c(lnr_case$formula, list(pGuess ~ 1)))
  p <- c(lnr_case$p, pGuess = qnorm(0.3)); p[["t0"]] <- log(0.405)
  for (bound in c("LC", "UC", "LT", "UT")) {
    dat <- bin_trial()
    lower <- bound %in% c("LC", "LT")
    dat[[bound]] <- if (lower) 0.42 else 0.445
    dadm <- bin_dadm(dat, des, 0.05)
    lo <- if (lower) 0.42 else 0.4
    hi <- if (lower) 0.45 else 0.445
    # The unbinned likelihood supplies the existing truncation normaliser.
    ub <- bin_dadm(dat, des, NULL)
    f <- function(t) vapply(t, function(tt) {
      ub$rt[] <- tt; exp(bin_ll(p, ub, des))
    }, 0)
    expect_equal(bin_ll(p, dadm, des),
      log(integrate(f, lo, hi, rel.tol = 1e-8)$value / 0.05), tolerance = 1e-5)
  }
})

test_that("upper grid endpoints remain inside the retained window", {
  des <- bin_design(LNR, lnr_case$formula)
  for (bound in c("UC", "UT")) {
    dat <- bin_trial(0.45); dat[[bound]] <- 0.45
    dadm <- bin_dadm(dat, des, 0.05)
    expect_equal(dadm$rt, rep(0.425, 2))
    expect_true(all(dadm$rt < dadm[[bound]]))
  }
  expect_error(EMC2:::.rt_bin_midpoint(0.45, 0.05, lower = 0.45, upper = 0.45),
               "no retained RT bin")
})

test_that("unknown responses equal the sum of generating response bins", {
  des <- bin_design(LNR, lnr_case$formula)
  p <- lnr_case$p; p[["m"]] <- -3.4; p[["t0"]] <- log(0.43)
  known <- vapply(c("left", "right"), function(r)
    exp(bin_ll(p, bin_dadm(bin_trial(response = r), des, 0.05), des)), 0)
  expect_equal(bin_ll(p, bin_dadm(bin_trial(response = NA_character_), des, 0.05), des),
               log(sum(known)), tolerance = 1e-7)
})

test_that("correlated bin likelihoods approach their independent limit", {
  des <- bin_design(LNRcorr, c(lnr_case$formula, list(rho ~ 1)))
  p <- c(lnr_case$p, rho = 0); p[["m"]] <- -3.4; p[["t0"]] <- log(0.43)
  dadm <- bin_dadm(bin_trial(), des, 0.05)
  independent <- bin_ll(p, dadm, des)
  for (rho in c(-1e-8, 1e-8)) {
    p[["rho"]] <- qnorm((rho + 1) / 2)
    expect_equal(bin_ll(p, dadm, des), independent, tolerance = 1e-6)
  }
})

test_that("shared-gamma counter bins agree with the integrated R density", {
  des <- bin_design(PCOUNTERcorr,
    list(nu ~ 1, c ~ 1, gamma ~ 1, k ~ 1, omega ~ 1, t0 ~ 1, rho ~ 1))
  p <- c(nu = log(12), c = log(0.3), gamma = -Inf, k = log(2),
         omega = -Inf, t0 = log(0.1), rho = qnorm(0.4))
  pars <- matrix(rep(c(12, 0.3, 0, 2, 0, 0.1, 0.4), 2), nrow = 2, byrow = TRUE,
    dimnames = list(NULL, c("nu", "c", "gamma", "k", "omega", "t0", "rho")))
  expected <- log(integrate(function(t) dPCOUNTERcorr(t, pars, 1),
                             0.1, 0.15, rel.tol = 1e-9)$value / 0.05)
  expect_equal(bin_ll(p, bin_dadm(bin_trial(0.14), des, 0.05), des),
               expected, tolerance = 1e-6)
})

test_that("drift and copula races retain their boundary-bin probability", {
  cases <- list(
    list(model = BAwLcorr,
      formula = list(v ~ 1, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1,
                     mG ~ 1, mK ~ 1, rho ~ 1),
      p = c(v = 1, sv = log(0.4), B = log(0.012), A = log(0.002),
            t0 = log(0.43), k = -Inf, mG = -Inf, mK = -Inf, rho = 0)),
    list(model = RDMSWTNcorr,
      formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1,
                     mG ~ 1, mK ~ 1, rho ~ 1),
      p = c(v = log(10), B = log(0.1), A = log(0.02), t0 = log(0.43),
            s = 0, sv = -Inf, mG = -Inf, mK = -Inf, rho = 0)))
  for (case in cases) {
    des <- bin_design(case$model, case$formula)
    dadm <- bin_dadm(bin_trial(), des, 0.05)
    p <- case$p
    independent <- bin_ll(p, dadm, des)
    expect_gt(independent, log(1e-10))
    p[["rho"]] <- qnorm((1 + 1e-8) / 2)
    expect_equal(bin_ll(p, dadm, des), independent, tolerance = 1e-6)
    p[["rho"]] <- qnorm((1 + 0.35) / 2)
    ub <- bin_dadm(bin_trial(), des, NULL)
    f <- function(t) vapply(t, function(tt) {
      ub$rt[] <- tt; exp(bin_ll(p, ub, des))
    }, 0)
    expect_equal(bin_ll(p, dadm, des),
      log(integrate(f, 0.43, 0.45, rel.tol = 1e-7)$value / 0.05),
      tolerance = 1e-5)
  }
})

test_that("timed bins retain responses before the response accumulator onset", {
  des <- design(factors = list(subjects = 1, S = "stim"),
    Rlevels = c("left", "right", "time"), model = LNR,
    formula = list(m ~ 0 + lR, s ~ 1, t0 ~ 0 + lR), verbose = FALSE,
    report_p_vector = FALSE)
  dat <- data.frame(subjects = factor(1), S = factor("stim"),
    R = factor("left", levels = c("left", "right", "time")), rt = 0.44)
  p <- sampled_pars(des, doMap = FALSE)
  p[c("m_lRleft", "m_lRright", "m_lRtime")] <- c(-1, -1, log(0.01))
  p[["s"]] <- log(0.4)
  p[c("t0_lRleft", "t0_lRright", "t0_lRtime")] <- log(c(0.46, 0.46, 0.43))
  expected <- log((plnorm(0.02, log(0.01), 0.4) / 2) / 0.05)
  dadm <- bin_dadm(dat, des, 0.05)
  expect_equal(bin_ll(p, dadm, des), expected, tolerance = 1e-6)
  P <- matrix(p, 1, dimnames = list(NULL, names(p)))
  expect_equal(sum(EMC2:::calc_ll_pw(P, dadm, des$model)), bin_ll(p, dadm, des),
               tolerance = 1e-8)
})

test_that("ordinary bins retain the single midpoint evaluation", {
  des <- bin_design(LNR, lnr_case$formula)
  dat <- bin_trial()
  dadm <- bin_dadm(dat, des, 1/60)
  ub <- bin_dadm(transform(dat, rt = dadm$rt[1]), des, NULL)
  expect_equal(bin_ll(lnr_case$p, dadm, des), bin_ll(lnr_case$p, ub, des),
               tolerance = 1e-12)
})

test_that("the boundary fallback turns on continuously", {
  des <- bin_design(LNR, lnr_case$formula)
  dadm <- bin_dadm(bin_trial(), des, 0.05)
  p <- lnr_case$p
  ll <- function(t0) { p[["t0"]] <- log(t0); bin_ll(p, dadm, des) }
  # The bin starts at .4; one width before it is .35.
  expect_equal(ll(0.35 - 1e-10), ll(0.35 + 1e-10), tolerance = 1e-7)
})

test_that("marginalisation includes feasible sub-millisecond bin support", {
  des <- bin_design(LNR, lnr_case$formula)
  p <- lnr_case$p; p[["m"]] <- log(0.0002); p[["m_lMTRUE"]] <- 0
  dadm <- bin_dadm(bin_trial(0.44), des, 0.001)
  P <- matrix(p, 1, dimnames = list(NULL, names(p)))
  grid <- EMC2:::compute_marginal_grid(P, dadm, des$model,
    list(param = "t0", mu = log(0.4405), sigma = 0.001, n_nodes = 48, n_scan = 15))
  expect_gt(max(exp(grid$nodes), na.rm = TRUE), 0.4405)
  expect_lt(max(exp(grid$nodes), na.rm = TRUE), 0.441)
})

test_that("uniform guesses remove the marginal t0 ceiling", {
  des <- bin_design(LNR, c(lnr_case$formula, list(pGuess ~ 1)))
  p <- c(lnr_case$p, pGuess = qnorm(0.5))
  dadm <- bin_dadm(bin_trial(), des, 0.05)
  P <- matrix(p, 1, dimnames = list(NULL, names(p)))
  grid <- EMC2:::compute_marginal_grid(P, dadm, des$model,
    list(param = "t0", mu = log(0.8), sigma = 0.1, n_nodes = 48))
  expect_gt(max(exp(grid$nodes), na.rm = TRUE), 0.8)
  expect_true(is.finite(grid$ll))
})

bin_erlang_cases <- list(
  list(model = RDMSWTN,
    formula = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1, mG ~ 1, mK ~ 1),
    p = c(v = log(10), B = log(0.1), A = log(0.02), t0 = log(0.43),
          s = 0, sv = -Inf, mG = log(0.3), mK = log(0.5))),
  list(model = BAwL,
    formula = list(v ~ 1, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1, mG ~ 1, mK ~ 1),
    p = c(v = 1, sv = log(0.4), B = log(0.012), A = log(0.002), t0 = log(0.43),
          k = -Inf, mG = log(0.3), mK = log(0.5))))

test_that("Erlang races integrate across their internal decision onset", {
  for (case in bin_erlang_cases) for (shape in list(1L, 2L, "mixed")) {
    for (type in c("local_guess", "local_kill_guess", "local_kill", "global_kill")) {
      if (identical(shape, "mixed") && type == "global_kill") next
      p <- case$p
      formula <- case$formula
      if (type == "local_guess") p[["mK"]] <- -Inf
      if (identical(shape, "mixed")) {
        formula <- c(formula, list(omega ~ 1)); p <- c(p, omega = qnorm(0.35))
      }
      des <- bin_design(case$model(erlang_shape = shape, erlang_type = type), formula)
      dadm <- bin_dadm(bin_trial(), des, 0.05)
      ub <- bin_dadm(bin_trial(), des, NULL)
      f <- function(t) vapply(t, function(tt) {
        ub$rt[] <- tt; exp(bin_ll(p, ub, des))
      }, 0)
      # Split the independent reference at t0: guess clocks use raw time,
      # while the decision component uses t - t0.
      expected <- log((integrate(f, 0.4, 0.43, rel.tol = 1e-8)$value +
                       integrate(f, 0.43, 0.45, rel.tol = 1e-8)$value) / 0.05)
      expect_equal(bin_ll(p, dadm, des), expected, tolerance = 1e-5,
                   info = paste(des$model()$c_name, shape, type))
      if (type %in% c("local_guess", "local_kill_guess")) {
        p[["t0"]] <- log(0.44)
        expected <- log((integrate(f, 0.4, 0.44, rel.tol = 1e-8)$value +
                         integrate(f, 0.44, 0.45, rel.tol = 1e-8)$value) / 0.05)
        expect_equal(bin_ll(p, dadm, des), expected, tolerance = 1e-5)
        ll <- function(t0) { p[["t0"]] <- log(t0); bin_ll(p, dadm, des) }
        for (edge in c(0.45, 0.5))
          expect_equal(ll(edge - 1e-10), ll(edge + 1e-10), tolerance = 1e-7)
        p[["t0"]] <- log(0.6)
        expected <- log(integrate(f, 0, 0.05, rel.tol = 1e-8)$value / 0.05)
        expect_equal(bin_ll(p, bin_dadm(bin_trial(0.01), des, 0.05), des),
                     expected, tolerance = 1e-6)
      }
    }
  }
})

test_that("Erlang guess bins keep their true-zero support before t0", {
  for (case in bin_erlang_cases) for (shape in list(1L, 2L, "mixed")) {
    p <- case$p; p[["mK"]] <- -Inf; p[["t0"]] <- log(0.8)
    formula <- case$formula
    if (identical(shape, "mixed")) {
      formula <- c(formula, list(omega ~ 1)); p <- c(p, omega = qnorm(0.35))
    }
    des <- bin_design(case$model(erlang_shape = shape, erlang_type = "local_guess"), formula)
    f <- function(t) {
      if (identical(shape, "mixed")) {
        pdf <- 0.35 * dgamma(t, 1, rate = 1 / 0.3) + 0.65 * dgamma(t, 2, rate = 2 / 0.3)
        surv <- 0.35 * pgamma(t, 1, rate = 1 / 0.3, lower.tail = FALSE) +
                0.65 * pgamma(t, 2, rate = 2 / 0.3, lower.tail = FALSE)
      } else {
        pdf <- dgamma(t, shape, rate = shape / 0.3)
        surv <- pgamma(t, shape, rate = shape / 0.3, lower.tail = FALSE)
      }
      pdf * surv
    }
    expected <- log(integrate(f, 0, 0.05, rel.tol = 1e-9)$value / 0.05)
    expect_equal(bin_ll(p, bin_dadm(bin_trial(0.01), des, 0.05), des),
                 expected, tolerance = 1e-6)
    expect_equal(bin_ll(p, bin_dadm(bin_trial(0.01, NA_character_), des, 0.05), des),
                 expected + log(2), tolerance = 1e-6)
  }
})

test_that("unequal Erlang clocks do not use the identical-clock shortcut", {
  for (case in bin_erlang_cases) {
    formula <- case$formula
    formula[[which(vapply(formula, function(f) as.character(f[[2]]) == "mG", TRUE))]] <- mG ~ lM
    p <- c(case$p, mG_lMTRUE = log(2))
    p[["t0"]] <- log(0.8); p[["mK"]] <- -Inf
    des <- bin_design(case$model(erlang_shape = 2L, erlang_type = "local_guess"), formula)
    expected <- log(integrate(function(t)
      dgamma(t, 2, rate = 2 / 0.6) * pgamma(t, 2, rate = 2 / 0.3, lower.tail = FALSE),
      0, 0.05, rel.tol = 1e-9)$value / 0.05)
    expect_equal(bin_ll(p, bin_dadm(bin_trial(0.01), des, 0.05), des),
                 expected, tolerance = 1e-6)
  }
})
