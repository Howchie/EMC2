# Racing Ornstein-Uhlenbeck with Smith (1995) Pulse Drift (ROUp)
# src/fpe_race.h, src/model_ROUp.h, R/model_ROUp.R.

roup_pars <- function(t, v_S, v_T, tau_S, tau_T, k, B, A, t0 = 0, s = 1) {
  n <- length(t)
  cbind(v_S = rep(v_S, n), v_T = rep(v_T, n), tau_S = rep(tau_S, n),
        tau_T = rep(tau_T, n), k = rep(k, n), B = rep(B, n), A = rep(A, n),
        t0 = rep(t0, n), s = rep(s, n))
}

tq <- seq(0.05, 2.0, by = 0.05)

test_that("ROUp reproduces ROU in the limit tau_S -> 0 and v_T = 0", {
  p_rou <- cbind(v = 1.5, k = 0.5, B = 1.0, A = 0.5, t0 = 0, s = 1.0)[rep(1, length(tq)), ]
  p_roup <- cbind(v_S = 1.5, v_T = 0.0, tau_S = 1e-15, tau_T = 1.0, k = 0.5,
                  B = 1.0, A = 0.5, t0 = 0, s = 1.0)[rep(1, length(tq)), ]

  pdf_rou <- EMC2:::dROU(tq, p_rou)
  pdf_roup <- EMC2:::dROUp(tq, p_roup)
  cdf_rou <- EMC2:::pROU(tq, p_rou)
  cdf_roup <- EMC2:::pROUp(tq, p_roup)

  expect_lt(max(abs(pdf_rou - pdf_roup)), 5e-4)
  expect_lt(max(abs(cdf_rou - cdf_roup)), 5e-4)
})

test_that("s is scaled out exactly in ROUp", {
  p1 <- roup_pars(tq, v_S = 1.5, v_T = 2.0, tau_S = 0.5, tau_T = 0.2, k = 0.5, B = 1.0, A = 0.5, s = 1.0)
  p2 <- roup_pars(tq, v_S = 3.0, v_T = 4.0, tau_S = 0.5, tau_T = 0.2, k = 0.5, B = 2.0, A = 1.0, s = 2.0)

  expect_identical(EMC2:::dROUp(tq, p1), EMC2:::dROUp(tq, p2))
  expect_identical(EMC2:::pROUp(tq, p1), EMC2:::pROUp(tq, p2))
})

test_that("one solve serves every row sharing a parameter tuple in ROUp", {
  n <- 300
  t_seq <- seq(0.1, 2.0, length.out = n)
  p <- roup_pars(t_seq, v_S = 1.5, v_T = 2.0, tau_S = 0.5, tau_T = 0.2, k = 0.5, B = 1.0, A = 0.5, s = 1.0)

  res <- EMC2:::droup_cpp(t_seq, p[, "v_S"], p[, "v_T"], p[, "tau_S"], p[, "tau_T"],
                          p[, "k"], p[, "B"], p[, "A"], p[, "t0"], p[, "s"])
  expect_identical(res$n_solves, 1L)
})

test_that("C++ particle likelihood matches R-side likelihood on fixed and collapsing bounds", {
  data <- data.frame(
    subjects = factor(rep(1, 10)),
    trials = 1:10,
    S = factor(rep("r1", 10), levels = c("r1", "r2")),
    R = factor(rep("r1", 10), levels = c("r1", "r2")),
    rt = rep(0.6, 10)
  )
  matchfun <- function(d) d$S == d$lR

  # Fixed boundary
  design_fixed <- design(data = data, model = ROUp, matchfun = matchfun,
                         formula = list(v_S ~ 1, v_T ~ 1, tau_S ~ 1, tau_T ~ 1,
                                        k ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
                         constants = c(s = log(1), A = log(0)))
  p_vec <- sampled_pars(design_fixed)
  p_vec[] <- c(log(1.5), log(2.5), log(0.3), log(0.1), log(0.5), log(1.0), log(0.2))
  dadm <- design_model(data, design_fixed)
  props <- matrix(p_vec, nrow = 1, dimnames = list(NULL, names(p_vec)))

  ll_r <- calc_ll_R(p_vec, design_fixed$model(), dadm)
  ll_cpp <- calc_ll_manager(props, dadm = dadm, model = design_fixed$model)
  expect_equal(ll_r, as.numeric(ll_cpp), tolerance = 1e-10)

  # Collapsing boundaries
  for (bnd in c("exponential", "linear_additive", "linear_multiplicative", "weibull")) {
    form <- if (bnd == "weibull") {
      list(v_S ~ 1, v_T ~ 1, tau_S ~ 1, tau_T ~ 1,
           k ~ 1, B ~ 1, A ~ 1, t0 ~ 1, Binf ~ 1, tau ~ 1, pw ~ 1)
    } else {
      list(v_S ~ 1, v_T ~ 1, tau_S ~ 1, tau_T ~ 1,
           k ~ 1, B ~ 1, A ~ 1, t0 ~ 1, Binf ~ 1, tau ~ 1)
    }
    design_c <- design(data = data,
                       model = function() ROUp(boundary_collapse = bnd),
                       matchfun = matchfun,
                       formula = form,
                       constants = c(s = log(1), A = log(0)))
    p_vec_c <- sampled_pars(design_c)
    p_vec_c[] <- c(log(1.5), log(2.5), log(0.3), log(0.1), log(0.5), log(1.0), log(0.2), log(0.4), log(0.8),
                   if (bnd == "weibull") log(1.2) else NULL)
    dadm_c <- design_model(data, design_c)
    props_c <- matrix(p_vec_c, nrow = 1, dimnames = list(NULL, names(p_vec_c)))

    ll_r_c <- calc_ll_R(p_vec_c, design_c$model(), dadm_c)
    ll_cpp_c <- calc_ll_manager(props_c, dadm = dadm_c, model = design_c$model)
    expect_equal(ll_r_c, as.numeric(ll_cpp_c), tolerance = 1e-10)
  }
})

test_that("make_data works with ROUp designs", {
  data <- data.frame(
    subjects = factor(rep(1, 10)),
    trials = 1:10,
    S = factor(rep("r1", 10), levels = c("r1", "r2")),
    R = factor(rep("r1", 10), levels = c("r1", "r2")),
    rt = rep(0.6, 10)
  )
  matchfun <- function(d) d$S == d$lR
  design_m <- design(data = data, model = ROUp, matchfun = matchfun,
                     formula = list(v_S ~ 1, v_T ~ 1, tau_S ~ 1, tau_T ~ 1,
                                    k ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
                     constants = c(s = log(1), A = log(0)))
  p_vec <- sampled_pars(design_m)
  p_vec[] <- c(log(1.5), log(2.5), log(0.3), log(0.1), log(0.5), log(1.0), log(0.2))

  sim_data <- make_data(p_vec, design_m, n_trials = 20)
  expect_s3_class(sim_data, "data.frame")
  expect_equal(nrow(sim_data), 40)
  expect_true(all(is.finite(sim_data$rt)))
  expect_true(all(sim_data$R %in% c("r1", "r2")))
})

test_that("ROUp simulator matches analytical CDF within sampling error", {
  skip_on_cran()
  set.seed(42)
  N <- 20000
  p <- roup_pars(seq_len(N), v_S = 1.0, v_T = 2.0, tau_S = 0.5, tau_T = 0.2,
                 k = 0.5, B = 2.0, A = 0, t0 = 0.2, s = 1.0)
  rts <- EMC2:::rROUp(lR = factor(rep("A", N)), pars = p, ok = rep(TRUE, N))

  t_eval <- seq(0.4, 2.0, by = 0.2)
  emp_cdf <- vapply(t_eval, function(t) mean(rts$rt <= t & !is.na(rts$rt)), numeric(1))
  anal_cdf <- EMC2:::pROUp(t_eval, p[seq_along(t_eval), ])

  # Kolmogorov-Smirnov / binomial tolerance: max deviation < 0.015 for N=20000
  expect_lt(max(abs(emp_cdf - anal_cdf)), 0.015)
})
