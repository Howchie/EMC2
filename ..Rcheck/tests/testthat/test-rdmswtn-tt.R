library(EMC2)

rdmswtn_tt_formula <- function(correlated = FALSE, rho_formula = rho ~ 1) {
  out <- list(
    v ~ 0 + lR, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1, tau ~ 1
  )
  if (correlated) out <- c(out, list(rho_formula))
  out
}

make_rdmswtn_tt_context <- function(
    data, model = RDMSWTN_TT(), rho_formula = rho ~ 1,
    functions = NULL, constants = NULL, compress = FALSE) {
  constants <- c(s = 0, constants)
  des <- design(
    data = data, Rlevels = levels(data$R),
    formula = rdmswtn_tt_formula(isTRUE(model$correlated), rho_formula),
    functions = functions, constants = constants, model = model,
    report_p_vector = FALSE
  )
  emc <- make_emc(
    data, des, type = "single", n_chains = 1, compress = compress,
    verbose = FALSE, rt_resolution = NULL
  )[[1]]
  list(design = des, emc = emc)
}

set_rdmswtn_tt_values <- function(
    p, rho = NULL, v = c(1.2, .8), B = 1, A = .2, t0 = .1,
    sv = .35, tau = 1.4) {
  p[grep("^v_", names(p))] <- log(v)
  p["B"] <- if (B == 0) 0 else log(B)
  p["A"] <- if (A == 0) 0 else log(A)
  p["t0"] <- log(t0)
  p["sv"] <- if (sv == 0) 0 else log(sv)
  p["tau"] <- log(tau)
  if (!is.null(rho)) {
    p[grep("^rho", names(p))] <- qnorm((rho + 1) / 2)
  }
  p
}

rdmswtn_tt_args <- function(ctx, p) {
  dadm <- ctx$emc$data[[1]]
  model <- ctx$emc$model()
  designs <- lapply(names(model$p_types), function(nm) {
    dm <- attr(dadm, "designs")[[nm]]
    dm[attr(dm, "expand"), , drop = FALSE]
  })
  names(designs) <- names(model$p_types)
  list(
    particle_matrix = matrix(
      p, nrow = 1, dimnames = list(NULL, names(p))
    ),
    data = dadm, constants = attr(dadm, "constants"),
    designs = designs, type = model$c_name, bounds = model$bound,
    transforms = model$transform, pretransforms = model$pre_transform,
    p_types = names(model$p_types), min_ll = log(1e-10),
    trend = model$trend, marginalise = NULL
  )
}

rdmswtn_tt_ll <- function(ctx, p) {
  do.call(EMC2:::calc_ll_oo, rdmswtn_tt_args(ctx, p))
}

test_that("RDMSWTN_TT constructors expose the dedicated clock contract", {
  model <- RDMSWTN_TT()
  io <- RDMSWTN_TT(posdrift = FALSE)
  corr <- RDMSWTN_TTcorr()

  expect_identical(model$c_name, "RDMSWTN_TT")
  expect_identical(io$c_name, "RDMSWTN_TT_IO")
  expect_identical(corr$c_name, "RDMSWTN_TT_CORR")
  expect_identical(
    model$p_types_canonical, c("v", "B", "A", "t0", "s", "sv", "tau")
  )
  expect_identical(
    names(model$p_types),
    c("v", "B", "A", "t0", "s", "sv", "tau", "pContaminant", "pGuess")
  )
  expect_false(any(c(
    "mG", "mK", "lambda_g", "lambda_k", "omega"
  ) %in% names(model$p_types)))
  expect_identical(model$transform$func[["tau"]], "exp")
  expect_equal(unname(model$p_types[["tau"]]), log(1))
  expect_gt(model$bound$minmax[1, "tau"], 0)
  expect_identical(tail(names(corr$p_types), 1L), "rho")
  expect_error(
    model$log_likelihood(NULL, NULL, model),
    "only in the compiled race path"
  )
})

test_that("RDMSWTN_TT density and CDF obey the operational-time identities", {
  cases <- expand.grid(
    A = c(0, .25), sv = c(0, .4), posdrift = c(FALSE, TRUE),
    KEEP.OUT.ATTRS = FALSE
  )
  for (i in seq_len(nrow(cases))) {
    z <- cases[i, ]
    pars <- cbind(
      v = 1.15, B = 1, A = z$A, t0 = .2, s = .9,
      sv = z$sv, tau = 1.6
    )
    pars <- RDMSWTN_TT(posdrift = z$posdrift)$Ttransform(pars, NULL)
    x <- c(.03, .2, .7, 1.3)
    t <- pars[1, "t0"] + x
    q <- x - x^2 / (2 * pars[1, "tau"])
    ordinary <- pars[rep(1, length(t)), , drop = FALSE]
    ordinary[, "t0"] <- 0
    ordinary <- cbind(ordinary, lambda_g = 0, lambda_k = 0)

    expect_equal(
      EMC2:::dRDMSWTN_TT(t, pars, posdrift = z$posdrift),
      EMC2:::dRDMSWTN(q, ordinary, posdrift = z$posdrift) *
        (1 - x / pars[1, "tau"]),
      tolerance = 2e-12
    )
    expect_equal(
      EMC2:::pRDMSWTN_TT(t, pars, posdrift = z$posdrift),
      EMC2:::pRDMSWTN(q, ordinary, posdrift = z$posdrift),
      tolerance = 2e-12
    )
  }
})

test_that("the exhaustion endpoint creates a frozen defective upper tail", {
  pars <- RDMSWTN_TT()$Ttransform(cbind(
    v = 1.1, B = 1, A = .2, t0 = .15, s = 1, sv = .3, tau = 1.25
  ), NULL)
  endpoint <- pars[1, "t0"] + pars[1, "tau"]
  plateau <- EMC2:::pRDMSWTN_TT(Inf, pars)

  expect_equal(
    EMC2:::dRDMSWTN_TT(c(-Inf, 0, pars[1, "t0"], endpoint, Inf), pars),
    rep(0, 5)
  )
  expect_equal(
    EMC2:::pRDMSWTN_TT(c(-Inf, pars[1, "t0"]), pars), c(0, 0)
  )
  expect_equal(
    EMC2:::pRDMSWTN_TT(c(endpoint, endpoint + 1, Inf), pars),
    rep(plateau, 3)
  )
  expect_lt(plateau, 1)

  integrated <- integrate(
    function(t) EMC2:::dRDMSWTN_TT(t, pars),
    pars[1, "t0"], endpoint, rel.tol = 2e-7
  )$value
  expect_equal(integrated, plateau, tolerance = 3e-6)

  t <- pars[1, "t0"] + .55
  h <- 1e-5
  derivative <- (
    EMC2:::pRDMSWTN_TT(t + h, pars) -
      EMC2:::pRDMSWTN_TT(t - h, pars)
  ) / (2 * h)
  expect_equal(derivative, EMC2:::dRDMSWTN_TT(t, pars), tolerance = 2e-5)
})

test_that("clock inversion is stable and the large-tau limit is ordinary RDMSWTN", {
  tau <- 3
  u <- c(0, 1e-16, 1e-10, tau / 2 * (1 - 1e-14), tau / 2)
  x <- vapply(u, EMC2:::rdmswtn_tt_qinv, numeric(1), tau = tau)
  expect_equal(x - x^2 / (2 * tau), u, tolerance = 2e-15)

  pars <- RDMSWTN_TT()$Ttransform(cbind(
    v = 1.3, B = 1, A = .2, t0 = .1, s = 1, sv = .25, tau = 1e9
  ), NULL)
  t <- c(.25, .5, 1)
  ordinary <- cbind(pars[rep(1, length(t)), , drop = FALSE],
                    lambda_g = 0, lambda_k = 0)
  expect_equal(
    EMC2:::dRDMSWTN_TT(t, pars),
    EMC2:::dRDMSWTN(t, ordinary),
    tolerance = 2e-8
  )
  expect_equal(
    EMC2:::pRDMSWTN_TT(t, pars),
    EMC2:::pRDMSWTN(t, ordinary),
    tolerance = 2e-8
  )
})

test_that("tau = Inf is the identity-clock limit in analytic and compiled paths", {
  pars <- RDMSWTN_TT()$Ttransform(cbind(
    v = 1.25, B = 1, A = .2, t0 = .1, s = 1, sv = .3, tau = Inf
  ), NULL)
  ordinary <- cbind(pars, lambda_g = 0, lambda_k = 0)
  t <- c(.25, .5, 1.2, Inf)

  expect_equal(
    EMC2:::dRDMSWTN_TT(t, pars),
    EMC2:::dRDMSWTN(t, ordinary),
    tolerance = 2e-12
  )
  expect_equal(
    EMC2:::pRDMSWTN_TT(t, pars),
    EMC2:::pRDMSWTN(t, ordinary),
    tolerance = 2e-12
  )

  row <- c(pars[1, ], pContaminant = 0)
  sim_pars <- matrix(
    rep(row, 2L * 300L), nrow = 2L * 300L, byrow = TRUE,
    dimnames = list(NULL, names(row))
  )
  set.seed(611)
  sim <- EMC2:::rrdmswtn_tt_cpp(
    sim_pars, c("a", "b"), rep(TRUE, nrow(sim_pars)), TRUE
  )
  expect_true(all(is.finite(sim$rt)))

  corr_row <- c(row, rho = .45)
  corr_pars <- matrix(
    rep(corr_row, 2L * 300L), nrow = 2L * 300L, byrow = TRUE,
    dimnames = list(NULL, names(corr_row))
  )
  set.seed(612)
  corr <- EMC2:::rrdmswtn_tt_corr_cpp(
    corr_pars, c("a", "b"), rep(TRUE, nrow(corr_pars)), TRUE
  )
  expect_true(all(is.finite(corr$rt)))
})

test_that("compiled and reference simulators reproduce omissions and support", {
  n <- 3000L
  row <- c(
    v = 1.1, B = 1, A = .2, t0 = .1, s = 1, sv = .35,
    tau = 1.2, pContaminant = 0
  )
  row <- RDMSWTN_TT()$Ttransform(matrix(
    row, nrow = 1, dimnames = list(NULL, names(row))
  ), NULL)[1, ]
  pars <- matrix(
    rep(row, 2L * n), nrow = 2L * n, byrow = TRUE,
    dimnames = list(NULL, names(row))
  )
  lR <- factor(rep(c("a", "b"), n), levels = c("a", "b"))

  set.seed(940)
  cpp <- EMC2:::rrdmswtn_tt_cpp(pars, levels(lR), rep(TRUE, 2L * n), TRUE)
  set.seed(941)
  ref <- EMC2:::rRDMSWTN_TT(lR, pars)
  expect_lt(
    abs(mean(is.infinite(cpp$rt)) - mean(is.infinite(ref$rt))), .035
  )
  expect_true(all(cpp$rt[is.finite(cpp$rt)] <= .1 + 1.2))
  expect_true(all(ref$rt[is.finite(ref$rt)] <= .1 + 1.2))
  expect_equal(
    median(cpp$rt[is.finite(cpp$rt)]),
    median(ref$rt[is.finite(ref$rt)]), tolerance = .04
  )
})

test_that("make_data dispatches through the compiled RDMSWTN_TT simulator", {
  des <- design(
    factors = list(S = "target", subjects = 1, difficulty = c("easy", "hard")),
    Rlevels = c("a", "b"),
    formula = list(
      v ~ 0 + lR, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1, tau ~ 1
    ),
    constants = c(s = 0), model = RDMSWTN_TT(),
    report_p_vector = FALSE
  )
  p <- sampled_pars(des, doMap = FALSE)
  p[] <- 0
  p[grep("^v_", names(p))] <- log(c(1.4, 1))
  p["B"] <- log(1)
  p["A"] <- log(.2)
  p["t0"] <- log(.1)
  p["sv"] <- log(.3)
  p["tau"] <- log(1.2)
  set.seed(92)
  sim <- make_data(p, design = des, n_trials = 30, rt_resolution = NULL)
  expect_s3_class(sim, "data.frame")
  expect_true(any(is.infinite(sim$rt)))
  expect_true(all(sim$rt[is.finite(sim$rt)] <= 1.3))
})

test_that("independent compiled likelihood handles finite, omitted, and truncated trials", {
  dat <- data.frame(
    subjects = factor(1),
    R = factor(c("a", NA), levels = c("a", "b")),
    rt = c(.65, Inf),
    LT = c(.2, 0), UT = c(1.1, Inf)
  )
  ctx <- make_rdmswtn_tt_context(dat)
  p <- set_rdmswtn_tt_values(
    sampled_pars(ctx$design, doMap = FALSE), tau = 1.25
  )
  ll <- rdmswtn_tt_ll(ctx, p)
  pw_args <- rdmswtn_tt_args(ctx, p)
  pw_args$marginalise <- NULL
  pw <- do.call(EMC2:::calc_ll_oo_pw, pw_args)
  dadm <- ctx$emc$data[[1]]
  pars <- EMC2:::get_pars_matrix_oo(p, dadm, ctx$emc$model())

  rows1 <- 1:2
  win <- rows1[which(dadm$winner[rows1])]
  lose <- setdiff(rows1, win)
  rt <- dadm$rt[win]
  log_num <- log(EMC2:::dRDMSWTN_TT(rt, pars[win, , drop = FALSE])) +
    log1p(-EMC2:::pRDMSWTN_TT(rt, pars[lose, , drop = FALSE]))
  race_surv <- function(t, rows) {
    prod(1 - EMC2:::pRDMSWTN_TT(t, pars[rows, , drop = FALSE]))
  }
  log_z <- log(race_surv(.2, rows1) - race_surv(1.1, rows1))
  omitted <- log(race_surv(Inf, 3:4))
  expect_equal(as.numeric(ll), log_num - log_z + omitted, tolerance = 2e-8)
  expect_equal(sum(pw), as.numeric(ll), tolerance = 2e-8)
})

test_that("compiled no-go likelihood retains finite omission mass", {
  dat <- data.frame(
    subjects = factor(1),
    R = factor(NA_character_, levels = c("a", "b", "nogo")),
    rt = Inf, UC = .9
  )
  ctx <- make_rdmswtn_tt_context(dat)
  p <- set_rdmswtn_tt_values(
    sampled_pars(ctx$design, doMap = FALSE),
    v = c(1.2, .8, 1), tau = 1.2
  )
  ll <- rdmswtn_tt_ll(ctx, p)
  pw_args <- rdmswtn_tt_args(ctx, p)
  pw_args$marginalise <- NULL
  pw <- do.call(EMC2:::calc_ll_oo_pw, pw_args)
  expect_true(is.finite(as.numeric(ll)))
  expect_gt(as.numeric(ll), log(1e-10))
  expect_equal(sum(pw), as.numeric(ll), tolerance = 2e-8)
})

test_that("correlated likelihood uses finite plateaus and nests rho zero", {
  dat <- data.frame(
    subjects = factor(1),
    R = factor(c("a", NA), levels = c("a", "b")),
    rt = c(.7, Inf),
    LT = c(.2, 0), UT = c(1.1, Inf)
  )
  corr <- make_rdmswtn_tt_context(dat, RDMSWTN_TTcorr())
  plain <- make_rdmswtn_tt_context(dat, RDMSWTN_TT())
  pc0 <- set_rdmswtn_tt_values(
    sampled_pars(corr$design, doMap = FALSE), rho = 0
  )
  pp <- set_rdmswtn_tt_values(sampled_pars(plain$design, doMap = FALSE))
  expect_equal(rdmswtn_tt_ll(corr, pc0), rdmswtn_tt_ll(plain, pp),
               tolerance = 0)

  rho <- .55
  pc <- set_rdmswtn_tt_values(
    sampled_pars(corr$design, doMap = FALSE), rho = rho
  )
  ll <- rdmswtn_tt_ll(corr, pc)
  dadm <- corr$emc$data[[1]]
  pars <- EMC2:::get_pars_matrix_oo(pc, dadm, corr$emc$model())

  rows <- 1:2
  win <- rows[which(dadm$winner[rows])]
  lose <- setdiff(rows, win)
  t <- dadm$rt[win]
  fw <- EMC2:::dRDMSWTN_TT(t, pars[win, , drop = FALSE])
  Fw <- EMC2:::pRDMSWTN_TT(t, pars[win, , drop = FALSE])
  Fl <- EMC2:::pRDMSWTN_TT(t, pars[lose, , drop = FALSE])
  finite <- log(fw) + pnorm(
    (rho * qnorm(Fw) - qnorm(Fl)) / sqrt(1 - rho^2),
    log.p = TRUE
  )
  pair_survival <- function(t, rows) {
    F <- EMC2:::pRDMSWTN_TT(t, pars[rows, , drop = FALSE])
    EMC2:::pbvn_tvpack(-qnorm(F[1]), -qnorm(F[2]), rho)
  }
  finite <- finite - log(
    pair_survival(.2, 1:2) - pair_survival(1.1, 1:2)
  )
  Fq <- EMC2:::pRDMSWTN_TT(Inf, pars[3:4, , drop = FALSE])
  omitted <- log(EMC2:::pbvn_tvpack(
    -qnorm(Fq[1]), -qnorm(Fq[2]), rho
  ))
  expect_equal(as.numeric(ll), finite + omitted, tolerance = 3e-8)
})

test_that("correlated simulator couples omission events and respects endpoints", {
  n <- 2500L
  row <- c(
    v = 1.05, B = 1, A = .2, t0 = .1, s = 1, sv = .3,
    tau = 1.1, pContaminant = 0, rho = .7
  )
  row <- RDMSWTN_TTcorr()$Ttransform(matrix(
    row, nrow = 1, dimnames = list(NULL, names(row))
  ), NULL)[1, ]
  pars <- matrix(
    rep(row, 2L * n), nrow = 2L * n, byrow = TRUE,
    dimnames = list(NULL, names(row))
  )
  set.seed(300)
  sim <- EMC2:::rrdmswtn_tt_corr_cpp(
    pars, c("a", "b"), rep(TRUE, 2L * n), TRUE
  )
  expect_true(any(is.infinite(sim$rt)))
  expect_true(all(sim$rt[is.finite(sim$rt)] <= 1.2))
  expect_equal(
    as.numeric(prop.table(table(sim$R))), c(.5, .5), tolerance = .07
  )
})

test_that("malformed correlated designs and unrestricted active rho are rejected", {
  dat3 <- data.frame(
    subjects = factor(1),
    R = factor("a", levels = c("a", "b", "c")),
    rt = .7
  )
  expect_error(
    make_rdmswtn_tt_context(dat3, RDMSWTN_TTcorr(), rho_formula = rho ~ 1),
    "at most two accumulator rows"
  )
  model <- RDMSWTN_TTcorr(posdrift = FALSE)
  pars <- matrix(
    c(1, 1, .2, .1, 1, .3, 1, 0, 0, .5), nrow = 1,
    dimnames = list(NULL, names(model$p_types))
  )
  expect_error(model$Ttransform(pars, NULL), "posdrift = FALSE")
})
