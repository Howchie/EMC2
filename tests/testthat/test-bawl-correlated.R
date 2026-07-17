library(EMC2)

log_sum_exp_test <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

bawl_formula <- function(model, rho_formula = NULL) {
  out <- list(v ~ 1, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1, k ~ 1,
              mG ~ 1, mK ~ 1)
  model_spec <- if (is.function(model)) model() else model
  if (isTRUE(model_spec$correlated)) out$rho <- rho_formula %||% (rho ~ 1)
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x

make_bawl_context <- function(data, model, rho_formula = NULL,
                               functions = NULL, constants = NULL,
                               matchfun = function(d) as.character(d$S) == as.character(d$lR),
                               formula = NULL) {
  formula <- formula %||% bawl_formula(model, rho_formula)
  des <- design(
    data = data,
    Rlevels = levels(data$R),
    matchfun = matchfun,
    functions = functions,
    formula = formula,
    constants = constants,
    model = model,
    report_p_vector = FALSE
  )
  emc <- make_emc(data, des, type = "single", n_chains = 1,
                  compress = FALSE, verbose = FALSE, rt_resolution = NULL)
  list(emc = emc[[1]], design = des)
}

context_args <- function(ctx, p) {
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
    data = dadm,
    constants = attr(dadm, "constants"),
    designs = designs,
    type = model$c_name,
    bounds = model$bound,
    transforms = model$transform,
    pretransforms = model$pre_transform,
    p_types = names(model$p_types),
    min_ll = log(1e-10),
    trend = model$trend
  )
}

set_bawl_values <- function(p, rho = NULL, v = 1, sv = .8, B = 1.4,
                            A = .2, t0 = .1, k = 0, mG = 1, mK = 1) {
  natural <- c(v = v, sv = sv, B = B, A = A, t0 = t0, k = k,
               mG = mG, mK = mK)
  for (nm in names(natural))
    p[nm] <- if (nm == "v") natural[[nm]] else log(natural[[nm]])
  if (!is.null(rho)) {
    rho_cols <- grep("^rho", names(p), value = TRUE)
    p[rho_cols] <- qnorm((rho + 1) / 2)
  }
  p
}

mapped_pars_for <- function(ctx, p) {
  EMC2:::get_pars_matrix_oo(p, ctx$emc$data[[1]], ctx$emc$model())
}

recovery_formula <- list(v ~ 0 + lR, sv ~ 1, B ~ 1, A ~ 1,
                         t0 ~ 1, k ~ 1, mG ~ 1, mK ~ 1, rho ~ 1)

set_recovery_v <- function(p) {
  v_cols <- grep("^v_", names(p), value = TRUE)
  p[v_cols] <- c(1, .5)
  p
}

reference_corr_ll <- function(pars, rho, rt, winner = 1L, posdrift = TRUE,
                              n = 80L, LT = 0, UT = Inf) {
  rule <- statmod::gauss.quad(n, kind = "hermite")
  log_num <- numeric(n)
  log_pos <- numeric(n)
  log_z <- numeric(n)
  for (i in seq_len(n)) {
    z <- sqrt(2) * rule$nodes[i]
    w <- log(rule$weights[i] / sqrt(pi))
    node <- pars
    for (r in seq_len(nrow(node))) {
      rr <- node[r, "rho"]
      node[r, "v"] <- node[r, "v"] + sign(rr) * node[r, "sv"] * sqrt(abs(rr)) * z
      node[r, "sv"] <- node[r, "sv"] * max(sqrt(1 - abs(rr)), 1e-12)
    }
    log_pos[i] <- if (posdrift)
      sum(pnorm(node[, "v"] / node[, "sv"], log.p = TRUE)) else 0
    d <- EMC2:::dBAwL(rt, node[winner, , drop = FALSE], posdrift = posdrift,
               erlang = 1L, guess = FALSE)
    survivor <- EMC2:::pBAwL(rt, node[-winner, , drop = FALSE], posdrift = posdrift,
                      erlang = 1L, guess = FALSE)
    log_num[i] <- w + log_pos[i] + log(d) + log1p(-survivor)
    if (LT == 0 && is.infinite(UT)) {
      log_z[i] <- 0
    } else {
      # pBAwL's vectorised wrapper recycles parameter rows over a vector of
      # RTs, so a multi-row matrix with one scalar RT must be evaluated
      # row-wise for the survivor product to include every racer.
      cdf_at <- function(t) vapply(seq_len(nrow(node)), function(r)
        EMC2:::pBAwL(t, node[r, , drop = FALSE], posdrift = posdrift,
                     erlang = 1L, guess = FALSE), numeric(1))
      s_lt <- prod(1 - cdf_at(LT))
      s_ut <- if (is.infinite(UT)) 0 else prod(1 - cdf_at(UT))
      log_z[i] <- log(max(s_lt - s_ut, .Machine$double.xmin))
    }
  }
  log_sum_exp_test(log_num) - if (posdrift)
    log_sum_exp_test(log(rule$weights / sqrt(pi)) + log_pos + log_z) else
    log_sum_exp_test(log(rule$weights / sqrt(pi)) + log_z)
}

test_that("BAwLcorr simulators use a jointly positive drift vector", {
  skip_on_cran()
  n <- 1500
  lR <- factor(rep(c("correct", "error", "pm"), n),
               levels = c("correct", "error", "pm"))
  pars <- cbind(
    v = rep(c(1.0, .4, .2), n), sv = 1, b = 1.5, A = .3,
    t0 = .2, k = .2, lambda_g = 0, lambda_k = 0,
    rho = rep(c(.6, -.6, 0), n)
  )
  ok <- rep(TRUE, nrow(pars))

  set.seed(101)
  cpp <- EMC2:::rbawl_corr_cpp(pars, levels(lR), ok, TRUE, 1L, FALSE, FALSE)
  set.seed(102)
  ref <- EMC2:::rBAwL_corr(lR, pars, ok = ok, posdrift = TRUE,
                           erlang = 1L, guess = FALSE, global = FALSE)

  expect_length(cpp$R, n)
  expect_equal(nrow(ref), n)
  expect_lt(abs(mean(cpp$R == 1, na.rm = TRUE) -
                mean(ref$R == "correct", na.rm = TRUE)), .06)
  expect_lt(abs(mean(cpp$R == 3, na.rm = TRUE) -
                mean(ref$R == "pm", na.rm = TRUE)), .06)
})

test_that("BAwLcorr requires lM and rejects row-varying rho", {
  skip_on_cran()
  dat <- data.frame(
    subjects = factor(1),
    S = factor("correct", levels = c("correct", "error", "pm")),
    R = factor("correct", levels = c("correct", "error", "pm")),
    rt = .7
  )
  expect_error(
    design(data = dat, model = BAwLcorr(),
           formula = bawl_formula(BAwLcorr(), rho ~ 1)),
    "requires matchfun"
  )

  des <- design(data = dat, model = BAwLcorr(),
                matchfun = function(d) as.character(d$S) == as.character(d$lR),
                formula = bawl_formula(BAwLcorr(), rho ~ 0 + lR),
                report_p_vector = FALSE)
  expect_error(
    make_emc(dat, des, type = "single", n_chains = 1,
             compress = FALSE, verbose = FALSE),
    "rho must be shared within each trial"
  )
})

test_that("cell-level rho gets role-specific signed loadings", {
  skip_on_cran()
  dat <- data.frame(
    subjects = factor(1),
    S = factor("correct", levels = c("correct", "error", "pm")),
    R = factor("correct", levels = c("correct", "error", "pm")),
    rt = .7
  )
  coupled <- function(d) factor(d$lR != "pm", levels = c(FALSE, TRUE),
                                labels = c("no", "yes"))
  ctx <- make_bawl_context(
    dat, BAwLcorr(), rho_formula = rho ~ 0 + coupled,
    functions = list(coupled = coupled), constants = c(rho_coupledno = 0)
  )
  p <- set_bawl_values(sampled_pars(ctx$design, doMap = FALSE), rho = .6)
  pars <- mapped_pars_for(ctx, p)
  expect_equal(as.numeric(pars[, "rho"]), c(.6, .6, 0), tolerance = 1e-12)

  p["rho_coupledyes"] <- qnorm(.2)
  pars <- mapped_pars_for(ctx, p)
  expect_equal(as.numeric(pars[, "rho"]), c(.6, -.6, 0), tolerance = 1e-12)
})

test_that("zero rho nests BAwL exactly across race data paths", {
  skip_on_cran()
  dat <- data.frame(
    subjects = factor(1),
    S = factor(c("correct", "error", "correct", "error"),
               levels = c("correct", "error")),
    R = factor(c("correct", NA, "error", "correct"),
               levels = c("correct", "error")),
    rt = c(.7, Inf, -Inf, .8),
    LT = c(0, 0, 0, .2), UT = c(Inf, Inf, Inf, 1.5),
    LC = c(0, 0, .4, 0), UC = c(Inf, Inf, Inf, Inf)
  )
  corr <- make_bawl_context(dat, BAwLcorr(posdrift = TRUE), rho_formula = rho ~ 1)
  ordinary <- make_bawl_context(dat, BAwL(posdrift = TRUE))
  p_corr <- set_bawl_values(sampled_pars(corr$design, doMap = FALSE), rho = 0)
  p_plain <- set_bawl_values(sampled_pars(ordinary$design, doMap = FALSE))
  ll_corr <- do.call(EMC2:::calc_ll_oo, context_args(corr, p_corr))
  ll_plain <- do.call(EMC2:::calc_ll_oo, context_args(ordinary, p_plain))
  expect_identical(as.numeric(ll_corr), as.numeric(ll_plain))
  expect_identical(
    as.numeric(do.call(EMC2:::calc_ll_oo_pw, context_args(corr, p_corr))),
    as.numeric(do.call(EMC2:::calc_ll_oo_pw, context_args(ordinary, p_plain)))
  )

  for (pos in c(FALSE, TRUE)) {
    for (kind in c("none", "local_kill_guess")) {
      corr_k <- make_bawl_context(dat, BAwLcorr(posdrift = pos,
                                                 erlang_type = kind),
                                  rho_formula = rho ~ 1)
      plain_k <- make_bawl_context(dat, BAwL(posdrift = pos,
                                             erlang_type = kind))
      pc <- set_bawl_values(sampled_pars(corr_k$design, doMap = FALSE), rho = 0,
                            mG = .8, mK = 1.2)
      pp <- set_bawl_values(sampled_pars(plain_k$design, doMap = FALSE),
                            mG = .8, mK = 1.2)
      expect_identical(
        as.numeric(do.call(EMC2:::calc_ll_oo, context_args(corr_k, pc))),
        as.numeric(do.call(EMC2:::calc_ll_oo, context_args(plain_k, pp)))
      )
    }
  }
})

test_that("an independent racer factors out of the correlated quadrature", {
  skip_on_cran()
  dat <- data.frame(
    subjects = factor(1),
    S = factor("correct", levels = c("correct", "error", "pm")),
    R = factor("correct", levels = c("correct", "error", "pm")),
    rt = .7
  )
  coupled <- function(d) factor(d$lR != "pm",
                                levels = c(FALSE, TRUE),
                                labels = c("no", "yes"))
  for (pos in c(FALSE, TRUE)) {
    ctx <- make_bawl_context(
      dat, BAwLcorr(posdrift = pos),
      rho_formula = rho ~ 0 + coupled,
      functions = list(coupled = coupled),
      constants = c(rho_coupledno = 0)
    )
    p <- set_bawl_values(sampled_pars(ctx$design, doMap = FALSE), rho = .6)
    fast <- as.numeric(do.call(EMC2:::calc_ll_oo, context_args(ctx, p)))

    generic_data <- ctx$emc$data[[1]]
    attr(generic_data, "emc2_all_finite_trials") <- FALSE
    generic_args <- context_args(ctx, p)
    generic_args$data <- generic_data
    generic <- as.numeric(do.call(EMC2:::calc_ll_oo, generic_args))
    expect_equal(fast, generic, tolerance = 1e-12,
                 info = paste("posdrift", pos))
  }
})

test_that("correlated likelihood agrees with a dense GH reference", {
  skip_on_cran()
  dat <- data.frame(subjects = factor(1),
                    S = factor("correct", levels = c("correct", "error")),
                    R = factor("correct", levels = c("correct", "error")),
                    rt = .7)
  for (pos in c(FALSE, TRUE)) {
    ctx_pos <- make_bawl_context(dat, BAwLcorr(posdrift = pos), rho_formula = rho ~ 1)
    for (rho in c(-.7, -.3, .3, .7)) {
      p <- set_bawl_values(sampled_pars(ctx_pos$design, doMap = FALSE), rho = rho)
      mapped <- mapped_pars_for(ctx_pos, p)
      got <- as.numeric(do.call(EMC2:::calc_ll_oo, context_args(ctx_pos, p)))
      got_pw <- do.call(EMC2:::calc_ll_oo_pw, context_args(ctx_pos, p))
      ref <- reference_corr_ll(mapped, rho = rho, rt = .7,
                               winner = 1L, posdrift = pos, n = 400L)
      expect_equal(got, ref, tolerance = 1e-5,
                   info = paste("rho", rho, "posdrift", pos))
      expect_equal(sum(got_pw), got, tolerance = 1e-10,
                   info = paste("pw rho", rho, "posdrift", pos))
    }
  }
})

test_that("truncated correlated likelihood uses an independent denominator", {
  skip_on_cran()
  dat <- data.frame(subjects = factor(1),
                    S = factor("correct", levels = c("correct", "error")),
                    R = factor("correct", levels = c("correct", "error")),
                    rt = .7, LT = .2, UT = 1.5)
  ctx <- make_bawl_context(dat, BAwLcorr(posdrift = TRUE), rho_formula = rho ~ 1)
  p <- set_bawl_values(sampled_pars(ctx$design, doMap = FALSE), rho = -.6)
  mapped <- mapped_pars_for(ctx, p)

  vals <- local({
    unconditional <- dat
    unconditional$LT <- 0
    unconditional$UT <- Inf
    ctx_u <- make_bawl_context(unconditional, BAwLcorr(posdrift = TRUE),
                               rho_formula = rho ~ 1)
    ll_u <- as.numeric(do.call(EMC2:::calc_ll_oo, context_args(ctx_u, p)))
    ll_t <- as.numeric(do.call(EMC2:::calc_ll_oo, context_args(ctx, p)))
    # The reference computes P(LT < T < UT) under the already positive-
    # conditioned model, independently of either C++ denominator.
    rule <- statmod::gauss.quad(200L, kind = "hermite")
    log_pos <- log_z <- numeric(200L)
    for (i in seq_len(200L)) {
      z <- sqrt(2) * rule$nodes[i]
      node <- mapped
      for (r in seq_len(nrow(node))) {
        rr <- node[r, "rho"]
        node[r, "v"] <- node[r, "v"] + sign(rr) * node[r, "sv"] * sqrt(abs(rr)) * z
        node[r, "sv"] <- node[r, "sv"] * max(sqrt(1 - abs(rr)), 1e-12)
      }
      log_pos[i] <- sum(pnorm(node[, "v"] / node[, "sv"], log.p = TRUE))
      # pBAwL's vectorised C++ wrapper recycles parameter vectors over a
      # vector of RTs, not one scalar RT over multiple parameter rows.  Call
      # it row-wise here so the independent denominator includes every racer.
      cdf_at <- function(t) vapply(seq_len(nrow(node)), function(r)
        EMC2:::pBAwL(t, node[r, , drop = FALSE], posdrift = TRUE,
                     erlang = 1L), numeric(1))
      s_lt <- prod(1 - cdf_at(.2))
      s_ut <- prod(1 - cdf_at(1.5))
      log_z[i] <- log(max(s_lt - s_ut, .Machine$double.xmin))
    }
    log_window <- log_sum_exp_test(log(rule$weights / sqrt(pi)) + log_pos + log_z) -
      log_sum_exp_test(log(rule$weights / sqrt(pi)) + log_pos)
    c(unconditional = ll_u, truncated = ll_t, window = log_window)
  })
  expect_equal(vals[["truncated"]], vals[["unconditional"]] - vals[["window"]],
               tolerance = 3e-5)
})

test_that("the direct correlation sign changes the likelihood and role sign", {
  skip_on_cran()
  dat <- data.frame(subjects = factor(1),
                    S = factor("correct", levels = c("correct", "error")),
                    R = factor("correct", levels = c("correct", "error")),
                    rt = .7)
  asym_formula <- list(v ~ 0 + lR, sv ~ 1, B ~ 1, A ~ 1,
                       t0 ~ 1, k ~ 1, mG ~ 1, mK ~ 1, rho ~ 1)
  ctx <- make_bawl_context(dat, BAwLcorr(posdrift = FALSE),
                           formula = asym_formula)
  set_asymmetric_v <- function(p) {
    v_cols <- grep("^v_", names(p), value = TRUE)
    p[v_cols] <- c(1, .45)
    p
  }
  ctx0 <- ctx
  ll_at <- function(rho, ctx = ctx0) {
    p <- set_asymmetric_v(set_bawl_values(
      sampled_pars(ctx$design, doMap = FALSE), rho = rho))
    as.numeric(do.call(EMC2:::calc_ll_oo, context_args(ctx, p)))
  }
  expect_gt(abs(ll_at(.6) - ll_at(-.6)), 1e-5)

  flipped <- make_bawl_context(
    dat, BAwLcorr(posdrift = FALSE), rho_formula = rho ~ 1,
    formula = asym_formula,
    matchfun = function(d) as.character(d$S) != as.character(d$lR)
  )
  p <- set_asymmetric_v(set_bawl_values(
    sampled_pars(flipped$design, doMap = FALSE), rho = -.6))
  pars <- mapped_pars_for(flipped, p)
  expect_equal(as.numeric(pars[, "rho"]), c(-.6, .6), tolerance = 1e-12)
  ll_flip <- as.numeric(do.call(EMC2:::calc_ll_oo, context_args(flipped, p)))
  # Reversing the match reference flips every factor loading.  Since the
  # shared latent factor is standard normal, that global sign reversal is an
  # equivalent parameterisation; the row-level mapping above is the
  # substantive assertion for this branch.
  expect_equal(ll_flip, ll_at(-.6), tolerance = 1e-12)
})

test_that("correlated likelihood quadrature covers the reviewer grid", {
  skip_on_cran()
  kinds <- c("regular", "truncated", "omission", "near_zero_large_sv")
  rhos <- c(.2, .5, .8, .95)
  make_grid_data <- function(kind) {
    R <- factor(if (kind == "omission") NA else "correct",
                levels = c("correct", "error"))
    out <- data.frame(
      subjects = factor(1),
      S = factor("correct", levels = c("correct", "error")),
      R = R,
      rt = if (kind == "omission") Inf else .7
    )
    if (kind == "truncated") {
      out$LT <- .2
      out$UT <- 1.5
    }
    out
  }
  values_for <- function(kind, rho) {
    near_zero <- kind == "near_zero_large_sv"
    ctx <- make_bawl_context(make_grid_data(kind), BAwLcorr(),
                             rho_formula = rho ~ 1)
    p <- set_bawl_values(sampled_pars(ctx$design, doMap = FALSE),
                          rho = rho, v = if (near_zero) .05 else .2,
                          sv = if (near_zero) 2 else .35)
    got <- as.numeric(do.call(EMC2:::calc_ll_oo, context_args(ctx, p)))
    # The single-trial GH reference covers the response kinds; omissions
    # exercise the generic evaluator and are asserted finite below.
    ref <- if (kind == "omission") NA_real_ else
      reference_corr_ll(mapped_pars_for(ctx, p), rho = rho, rt = .7,
                        winner = 1L, posdrift = TRUE, n = 400L,
                        LT = if (kind == "truncated") .2 else 0,
                        UT = if (kind == "truncated") 1.5 else Inf)
    c(rho = rho, got = got, ref = ref)
  }
  grid <- as.data.frame(do.call(rbind, lapply(kinds, function(kind)
    do.call(rbind, lapply(rhos, function(rho)
      c(kind = kind, values_for(kind, rho)))))),
                        stringsAsFactors = FALSE)
  numeric_cols <- c("rho", "got", "ref")
  grid[numeric_cols] <- lapply(grid[numeric_cols], as.numeric)
  expect_true(all(is.finite(grid$got)))
  with_ref <- !is.na(grid$ref)
  expect_lt(max(abs(grid$got[with_ref] - grid$ref[with_ref])), 1e-4)

  sweep <- seq(.1, .95, by = .05)
  sweep_ll <- vapply(sweep, function(rho) {
    ctx <- make_bawl_context(make_grid_data("near_zero_large_sv"),
                             BAwLcorr(), rho_formula = rho ~ 1)
    p <- set_bawl_values(sampled_pars(ctx$design, doMap = FALSE),
                          rho = rho, v = .05, sv = 2)
    as.numeric(do.call(EMC2:::calc_ll_oo, context_args(ctx, p)))
  }, numeric(1))
  expect_true(all(is.finite(sweep_ll)))
  expect_lt(max(abs(diff(sweep_ll, differences = 2))), .5)
})

test_that("coarse likelihood recovers simulator correlation", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("EMC2_RUN_LONG_BAWLCORR"), "true"))
  n <- 800L
  pars <- cbind(v = rep(c(1, .5), n), sv = .8, b = 1.6, A = .2,
                t0 = .1, k = 0, lambda_g = 0, lambda_k = 0,
                rho = rep(c(.6, .6), n))
  set.seed(12)
  sim <- EMC2:::rbawl_corr_cpp(pars, c("correct", "error"),
                               rep(TRUE, nrow(pars)), FALSE, 1L, FALSE, FALSE)
  dat <- data.frame(subjects = factor(1),
                    S = factor(rep("correct", n), levels = c("correct", "error")),
                    R = factor(ifelse(sim$R == 1, "correct", "error"),
                               levels = c("correct", "error")),
                    rt = sim$rt)
  ctx <- make_bawl_context(dat, BAwLcorr(posdrift = FALSE),
                           formula = recovery_formula)
  grid <- c(-.6, 0, .3, .6, .9)
  lls <- vapply(grid, function(rho) {
    p <- set_recovery_v(set_bawl_values(
      sampled_pars(ctx$design, doMap = FALSE), rho = rho))
    as.numeric(do.call(EMC2:::calc_ll_oo, context_args(ctx, p)))
  }, numeric(1))
  expect_equal(grid[which.max(lls)], .6)
})

test_that("long-run dependence recovery covers both correlation signs", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("EMC2_RUN_LONG_BAWLCORR"), "true"))
  # The long-runner is intentionally opt-in: it evaluates the likelihood on
  # a coarse direct-correlation grid for each generating value rather than
  # spending the ordinary package test budget on MCMC.
  truth_grid <- c(-.7, -.3, 0, .3, .7)
  recovered <- vapply(truth_grid, function(truth) {
    n <- 2000L
    cell_rho <- if (truth < 0) c(abs(truth), truth) else c(truth, truth)
    pars <- cbind(v = rep(c(1, .5), n), sv = .8, b = 1.6, A = .2,
                   t0 = .1, k = 0, lambda_g = 0, lambda_k = 0,
                   rho = rep(cell_rho, n))
    set.seed(1000 + as.integer(round(100 * (truth + 1))))
    sim <- EMC2:::rbawl_corr_cpp(pars, c("correct", "error"),
                                 rep(TRUE, nrow(pars)), FALSE, 1L, FALSE, FALSE)
    dat <- data.frame(
      subjects = factor(1),
      S = factor(rep("correct", n), levels = c("correct", "error")),
      R = factor(ifelse(sim$R == 1, "correct", "error"),
                 levels = c("correct", "error")),
      rt = sim$rt
    )
    ctx <- make_bawl_context(dat, BAwLcorr(posdrift = FALSE),
                             formula = recovery_formula)
    lls <- vapply(truth_grid, function(candidate) {
      p <- set_recovery_v(set_bawl_values(
        sampled_pars(ctx$design, doMap = FALSE), rho = candidate))
      as.numeric(do.call(EMC2:::calc_ll_oo, context_args(ctx, p)))
    }, numeric(1))
    truth_grid[which.max(lls)]
  }, numeric(1))
  expect_equal(recovered, truth_grid)
})
