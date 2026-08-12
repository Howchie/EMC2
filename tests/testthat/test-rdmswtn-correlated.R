library(EMC2)

rdmswtn_corr_formula <- function(correlated = TRUE, rho_formula = rho ~ 1) {
  out <- list(v ~ 0 + lR, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1,
              mG ~ 1, mK ~ 1)
  if (correlated) out <- c(out, list(rho_formula))
  out
}

make_rdmswtn_corr_context <- function(data, model = RDMSWTNcorr(),
                                      rho_formula = rho ~ 1,
                                      functions = NULL, constants = NULL,
                                      compress = FALSE) {
  correlated <- isTRUE(model$correlated)
  constants <- c(s = 0, mG = 0, mK = 0, constants)
  des <- design(
    data = data, Rlevels = levels(data$R),
    formula = rdmswtn_corr_formula(correlated, rho_formula),
    functions = functions, constants = constants, model = model,
    report_p_vector = FALSE
  )
  emc <- make_emc(data, des, type = "single", n_chains = 1,
                  compress = compress, verbose = FALSE,
                  rt_resolution = NULL)[[1]]
  list(design = des, emc = emc)
}

set_rdmswtn_corr_values <- function(p, rho = NULL, v = c(1.2, .8),
                                    B = 1, A = .2, t0 = .1, sv = .4) {
  p[grep("^v_", names(p))] <- log(v)
  p["B"] <- log(B)
  p["A"] <- if (A == 0) 0 else log(A)
  p["t0"] <- log(t0)
  # sv is sampled on the log scale, so sv = 0 is log(0) = -Inf; a plain 0
  # here would silently mean exp(0) = 1.
  p["sv"] <- log(sv)
  if (!is.null(rho)) {
    p[grep("^rho", names(p))] <- qnorm((rho + 1) / 2)
  }
  p
}

rdmswtn_corr_args <- function(ctx, p) {
  dadm <- ctx$emc$data[[1]]
  model <- ctx$emc$model()
  designs <- lapply(names(model$p_types), function(nm) {
    dm <- attr(dadm, "designs")[[nm]]
    dm[attr(dm, "expand"), , drop = FALSE]
  })
  names(designs) <- names(model$p_types)
  list(
    particle_matrix = matrix(p, nrow = 1,
                             dimnames = list(NULL, names(p))),
    data = dadm, constants = attr(dadm, "constants"),
    designs = designs, type = model$c_name, bounds = model$bound,
    transforms = model$transform, pretransforms = model$pre_transform,
    p_types = names(model$p_types), min_ll = log(1e-10),
    trend = model$trend, marginalise = NULL
  )
}

rdmswtn_corr_ll <- function(ctx, p) {
  do.call(EMC2:::calc_ll_oo, rdmswtn_corr_args(ctx, p))
}

test_that("RDMSWTN correlated constructors expose the copula contract", {
  ordinary <- RDMSWTN()
  corr <- RDMSWTN(correlated = TRUE)
  convenience <- RDMSWTNcorr()

  expect_false(ordinary$correlated)
  expect_identical(corr$c_name, "RDMSWTN_CORR")
  expect_identical(convenience$c_name, corr$c_name)
  expect_identical(corr$correlation_type, "rdmswtn_gaussian_copula")
  expect_identical(tail(names(corr$p_types), 1), "rho")
  expect_equal(corr$transform$lower[["rho"]], -1)
  expect_equal(corr$transform$upper[["rho"]], 1)
  expect_equal(corr$bound$minmax[, "rho"], c(-.99, .99))
  expect_equal(corr$bound$exception[["rho"]], 0)
  expect_error(RDMSWTNcorr(erlang_type = "local_guess"),
               "does not support guess or kill clocks")
  expect_error(RDMSWTN(correlated = TRUE, erlang_type = "local_kill"),
               "does not support guess or kill clocks")
})

test_that("binary copula likelihood matches the analytical cause formula", {
  dat <- data.frame(
    subjects = factor(1),
    R = factor(c("a", "b"), levels = c("a", "b")),
    rt = c(.7, .8)
  )
  for (rho in c(-.7, 0, .65)) {
    for (shape in list(c(sv = 0, A = 0), c(sv = .45, A = .25))) {
      ctx <- make_rdmswtn_corr_context(dat)
      p <- set_rdmswtn_corr_values(
        sampled_pars(ctx$design, doMap = FALSE), rho = rho,
        sv = shape[["sv"]], A = shape[["A"]])
      ll <- rdmswtn_corr_ll(ctx, p)
      dadm <- ctx$emc$data[[1]]
      pars <- EMC2:::get_pars_matrix_oo(p, dadm, ctx$emc$model())

      reference <- 0
      for (j in 0:1) {
        rows <- j * 2 + 1:2
        win <- rows[which(dadm$winner[rows])]
        lose <- setdiff(rows, win)
        t <- dadm$rt[win]
        fw <- EMC2:::dRDMSWTN(t, pars[win, , drop = FALSE])
        Fw <- EMC2:::pRDMSWTN(t, pars[win, , drop = FALSE])
        Fl <- EMC2:::pRDMSWTN(t, pars[lose, , drop = FALSE])
        conditional_survivor <- pnorm(
          (rho * qnorm(Fw) - qnorm(Fl)) / sqrt(1 - rho^2))
        reference <- reference + log(fw) + log(conditional_survivor)
      }
      expect_equal(as.numeric(ll), reference, tolerance = 2e-8,
                   info = paste("rho", rho, "sv", shape[["sv"]],
                                "A", shape[["A"]]))
    }
  }
})

test_that("rho zero nests ordinary RDMSWTN for finite and omitted trials", {
  dat <- data.frame(
    subjects = factor(1),
    R = factor(c("a", NA, "b", "a"), levels = c("a", "b")),
    rt = c(.7, Inf, -Inf, .8),
    LT = c(0, 0, 0, .2), UT = c(Inf, Inf, Inf, 1.5),
    LC = c(0, 0, .4, 0), UC = c(Inf, Inf, Inf, Inf)
  )
  corr <- make_rdmswtn_corr_context(dat, RDMSWTNcorr())
  plain <- make_rdmswtn_corr_context(dat, RDMSWTN())
  pc <- set_rdmswtn_corr_values(
    sampled_pars(corr$design, doMap = FALSE), rho = 0)
  pp <- set_rdmswtn_corr_values(
    sampled_pars(plain$design, doMap = FALSE))
  expect_equal(rdmswtn_corr_ll(corr, pc), rdmswtn_corr_ll(plain, pp),
               tolerance = 0)
})

test_that("two-plus-independent races use pair survival for opted-out winners", {
  dat <- data.frame(
    subjects = factor(1),
    R = factor("c", levels = c("a", "b", "c")),
    rt = .75
  )
  coupled <- function(d) factor(
    d$lR != "c", levels = c(FALSE, TRUE), labels = c("no", "yes"))
  ctx <- make_rdmswtn_corr_context(
    dat, rho_formula = rho ~ 0 + coupled,
    functions = list(coupled = coupled),
    constants = c(rho_coupledno = 0))
  p <- set_rdmswtn_corr_values(
    sampled_pars(ctx$design, doMap = FALSE), rho = .6,
    v = c(1.3, .9, 1.05))
  ll <- rdmswtn_corr_ll(ctx, p)
  dadm <- ctx$emc$data[[1]]
  pars <- EMC2:::get_pars_matrix_oo(p, dadm, ctx$emc$model())
  t <- dadm$rt[1]
  F <- EMC2:::pRDMSWTN(rep(t, 3), pars)
  f3 <- EMC2:::dRDMSWTN(t, pars[3, , drop = FALSE])
  pair_survival <- EMC2:::pbvn_tvpack(-qnorm(F[1]), -qnorm(F[2]), .6)
  expect_equal(as.numeric(ll), log(f3) + log(pair_survival),
               tolerance = 2e-8)
})

test_that("direct-pair validation rejects malformed participation", {
  dat3 <- data.frame(
    subjects = factor(1),
    R = factor("a", levels = c("a", "b", "c")),
    rt = .7
  )
  bad <- make_rdmswtn_corr_context
  expect_error(
    bad(dat3, rho_formula = rho ~ 1),
    "at most two accumulator rows")

  model <- RDMSWTNcorr(posdrift = FALSE)
  pars <- matrix(
    c(1, 1, .2, .1, 1, .4, 1, 1, 0, 0, .5),
    nrow = 1,
    dimnames = list(NULL, names(model$p_types)))
  expect_error(model$Ttransform(pars, NULL),
               "requires sv = 0 on the correlated rows")
  # sv = 0 is the supported unrestricted-drift copula; it must pass.
  pars[, "sv"] <- 0
  expect_silent(model$Ttransform(pars, NULL))
  # rho = 0 rows are never subject to the sv rule.
  pars[, c("sv", "rho")] <- c(.4, 0)
  expect_silent(model$Ttransform(pars, NULL))
})

test_that("posdrift selects the sampling scale of the RDMSWTN mean rate", {
  for (ctor in list(RDMSWTN, RDMSWTN_TT, RDMSWTNcorr, RDMSWTN_TTcorr,
                    LogicalRulesRDMSWTN)) {
    pos <- ctor(posdrift = TRUE)
    io <- ctor(posdrift = FALSE)
    expect_identical(pos$transform$func[["v"]], "exp")
    expect_equal(pos$p_types[["v"]], log(1))
    expect_equal(pos$bound$minmax[, "v"], c(1e-3, Inf))
    # Unrestricted drifts must be able to reach negative mean rates.
    expect_identical(io$transform$func[["v"]], "identity")
    expect_equal(io$p_types[["v"]], 1)
    expect_equal(io$bound$minmax[, "v"], c(-Inf, Inf))
  }
})

test_that("unrestricted-drift copula matches the defective cause formula", {
  dat <- data.frame(
    subjects = factor(1),
    R = factor(c("a", "b"), levels = c("a", "b")),
    rt = c(.7, .8)
  )
  for (rho in c(-.75, .55)) {
    ctx <- make_rdmswtn_corr_context(dat, RDMSWTNcorr(posdrift = FALSE),
                                     constants = c(sv = log(0)))
    p <- sampled_pars(ctx$design, doMap = FALSE)
    # v is on the natural scale here: one accumulator is defective.
    p[grep("^v_", names(p))] <- c(-0.35, 0.9)
    p["B"] <- log(1); p["A"] <- log(.2); p["t0"] <- log(.1)
    p[grep("^rho", names(p))] <- qnorm((rho + 1) / 2)
    ll <- rdmswtn_corr_ll(ctx, p)

    dadm <- ctx$emc$data[[1]]
    pars <- EMC2:::get_pars_matrix_oo(p, dadm, ctx$emc$model())
    reference <- 0
    for (j in 0:1) {
      rows <- j * 2 + 1:2
      win <- rows[which(dadm$winner[rows])]
      lose <- setdiff(rows, win)
      t <- dadm$rt[win]
      fw <- EMC2:::dRDMSWTN(t, pars[win, , drop = FALSE], posdrift = FALSE)
      Fw <- EMC2:::pRDMSWTN(t, pars[win, , drop = FALSE], posdrift = FALSE)
      Fl <- EMC2:::pRDMSWTN(t, pars[lose, , drop = FALSE], posdrift = FALSE)
      conditional_survivor <- pnorm(
        (rho * qnorm(Fw) - qnorm(Fl)) / sqrt(1 - rho^2))
      reference <- reference + log(fw) + log(conditional_survivor)
    }
    expect_equal(as.numeric(ll), reference, tolerance = 2e-8,
                 info = paste("rho", rho))
  }
})

test_that("the copula couples intrinsic omissions of defective marginals", {
  # Both accumulators have negative mean rates, so the pair has a joint atom
  # at +Inf.  Its probability is the bivariate normal upper orthant at the
  # marginal plateaus, NOT the product of the marginal never-finish masses.
  dat <- data.frame(
    subjects = factor(1),
    R = factor(NA, levels = c("a", "b")),
    rt = Inf
  )
  for (rho in c(-.7, .7)) {
    ctx <- make_rdmswtn_corr_context(dat, RDMSWTNcorr(posdrift = FALSE),
                                     constants = c(sv = log(0)))
    p <- sampled_pars(ctx$design, doMap = FALSE)
    p[grep("^v_", names(p))] <- c(-0.3, -0.6)
    p["B"] <- log(1); p["A"] <- log(.2); p["t0"] <- log(.1)
    p[grep("^rho", names(p))] <- qnorm((rho + 1) / 2)
    ll <- rdmswtn_corr_ll(ctx, p)

    pars <- EMC2:::get_pars_matrix_oo(p, ctx$emc$data[[1]],
                                      ctx$emc$model())
    cap <- EMC2:::pRDMSWTN(rep(Inf, 2), pars, posdrift = FALSE)
    expect_true(all(cap < 1))
    joint_atom <- EMC2:::pbvn_tvpack(-qnorm(cap[1]), -qnorm(cap[2]), rho)
    # The compiled route takes the Drezner branch of norm_cdf_2d_hybrid for
    # these moderate arguments, which carries ~1e-6 relative error against
    # tvpack; that tradeoff is shared with every other correlated route.
    expect_equal(as.numeric(ll), log(joint_atom), tolerance = 2e-5,
                 info = paste("rho", rho))
    expect_false(isTRUE(all.equal(joint_atom,
                                  prod(1 - cap), tolerance = 1e-3)))
  }
})

test_that("unrestricted-drift copula simulator reproduces the likelihood", {
  skip_on_cran()
  model <- RDMSWTNcorr(posdrift = FALSE)
  lR <- factor(rep(c("a", "b"), 2e4), levels = c("a", "b"))
  pars <- cbind(
    v = rep(c(-0.3, -0.6), 2e4), b = 1.2, A = .2, t0 = .1, s = 1,
    sv = 0, lambda_g = 0, lambda_k = 0, rho = .7)
  set.seed(11)
  sim <- EMC2:::rRDMSWTN_corr(lR, pars, posdrift = FALSE)

  cap <- EMC2:::pRDMSWTN(c(Inf, Inf), pars[1:2, , drop = FALSE],
                         posdrift = FALSE)
  joint_atom <- EMC2:::pbvn_tvpack(-qnorm(cap[1]), -qnorm(cap[2]), .7)
  cause <- function(t, w) {
    if (!is.finite(t)) return(0)  # the atom at +Inf is not part of the density
    l <- 3L - w
    fw <- EMC2:::dRDMSWTN(t, pars[w, , drop = FALSE], posdrift = FALSE)
    if (!(fw > 0)) return(0)
    Fw <- EMC2:::pRDMSWTN(t, pars[w, , drop = FALSE], posdrift = FALSE)
    Fl <- EMC2:::pRDMSWTN(t, pars[l, , drop = FALSE], posdrift = FALSE)
    fw * pnorm((.7 * qnorm(Fw) - qnorm(Fl)) / sqrt(1 - .7^2))
  }
  p_resp <- vapply(1:2, function(w) {
    integrate(function(t) vapply(t, cause, 0, w = w), .1, Inf,
              rel.tol = 1e-9)$value
  }, 0)
  # The simulator, the copula likelihood and the joint atom must agree, and
  # together they must exhaust the probability.
  expect_equal(sum(p_resp) + joint_atom, 1, tolerance = 1e-6)
  expect_equal(mean(!is.finite(sim$rt)), joint_atom, tolerance = .01)
  expect_equal(mean(sim$R == "a", na.rm = TRUE) *
                 mean(is.finite(sim$rt)), p_resp[1], tolerance = .03)

  # The compiled simulator inverts the same defective marginal: uniforms at or
  # above the plateau must become never-finish draws rather than an error.
  set.seed(12)
  cpp <- EMC2:::rrdmswtn_corr_cpp(pars, levels(lR),
                                  rep(TRUE, nrow(pars)), FALSE)
  expect_equal(mean(!is.finite(cpp$rt)), joint_atom, tolerance = .01)
  expect_equal(mean(cpp$R[is.finite(cpp$rt)] == 1L) *
                 mean(is.finite(cpp$rt)), p_resp[1], tolerance = .03)
  expect_equal(as.numeric(quantile(cpp$rt[is.finite(cpp$rt)], c(.25, .5, .75))),
               as.numeric(quantile(sim$rt[is.finite(sim$rt)], c(.25, .5, .75))),
               tolerance = .03)
})

test_that("compiled and R simulators reject sv > 0 under posdrift = FALSE", {
  lR <- factor(rep(c("a", "b"), 5), levels = c("a", "b"))
  pars <- cbind(
    v = rep(c(-0.3, 0.6), 5), b = 1.2, A = .2, t0 = .1, s = 1,
    sv = .3, lambda_g = 0, lambda_k = 0, rho = .5)
  expect_error(EMC2:::rRDMSWTN_corr(lR, pars, posdrift = FALSE),
               "requires sv = 0 on the correlated rows")
  expect_error(
    EMC2:::rrdmswtn_corr_cpp(pars, levels(lR), rep(TRUE, nrow(pars)), FALSE),
    "requires sv = 0 on the correlated rows")
  # sv > 0 stays legal once the rows are uncorrelated, and under posdrift.
  pars[, "rho"] <- 0
  expect_silent(EMC2:::rRDMSWTN_corr(lR, pars, posdrift = FALSE))
  pars[, "rho"] <- .5
  expect_silent(EMC2:::rRDMSWTN_corr(lR, pars, posdrift = TRUE))
})

test_that("compiled dispatch rejects correlated logical-rule races", {
  dat <- data.frame(
    subjects = factor(1),
    R = factor("a", levels = c("a", "b")),
    rt = .7
  )
  ctx <- make_rdmswtn_corr_context(dat)
  p <- set_rdmswtn_corr_values(
    sampled_pars(ctx$design, doMap = FALSE), rho = .5)
  args <- rdmswtn_corr_args(ctx, p)
  args$type <- "LogicalRulesRDMSWTN_CORR"
  expect_error(
    do.call(EMC2:::calc_ll_oo, args),
    "correlated RDMSWTN logical-rule races are not supported"
  )
})

test_that("extreme times and correlations remain finite and continuous", {
  dat <- data.frame(
    subjects = factor(1),
    R = factor(c("a", "b", "a", "b"), levels = c("a", "b")),
    rt = c(.1001, .2, 2, 12)
  )
  ctx <- make_rdmswtn_corr_context(dat)
  base <- sampled_pars(ctx$design, doMap = FALSE)
  values <- vapply(
    c(-.99, -.989, .989, .99),
    function(rho) {
      p <- set_rdmswtn_corr_values(base, rho = rho, A = .25, sv = .5)
      as.numeric(rdmswtn_corr_ll(ctx, p))
    },
    numeric(1)
  )
  expect_true(all(is.finite(values)))
  expect_lt(abs(values[1] - values[2]), 5)
  expect_lt(abs(values[3] - values[4]), 5)
})

test_that("correlated simulation retains symmetric race marginals", {
  n <- 1500L
  row <- c(v = 1.1, b = 1.2, A = .2, t0 = .1, s = 1, sv = .5,
           lambda_g = 0, lambda_k = 0, pContaminant = 0, rho = .75)
  pars <- matrix(rep(row, 2L * n), nrow = 2L * n, byrow = TRUE,
                 dimnames = list(NULL, names(row)))
  set.seed(841)
  sim <- EMC2:::rrdmswtn_corr_cpp(
    pars, c("a", "b"), rep(TRUE, 2L * n), TRUE)
  expect_equal(as.numeric(prop.table(table(sim$R))), c(.5, .5),
               tolerance = .04)
  expect_true(all(is.finite(sim$rt)))
  expect_true(all(sim$rt > .1))
})
