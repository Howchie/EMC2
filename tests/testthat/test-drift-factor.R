library(EMC2)

# The shared correlated-drift factor (src/drift_factor.h) is used by BAwLcorr
# and by RDMSWTNcorr/RDMSWTN_TTcorr with correlate = "drifts".  These tests
# cover the RDMSWTN side and the pieces both models now share; the BAwL exact
# and fused kernels stay in test-bawl-correlated.R.

df_data <- function(n = 24, seed = 3, levels = c("a", "b")) {
  set.seed(seed)
  data.frame(
    subjects = factor(rep(1, n)),
    S = factor(rep(levels, length.out = n), levels = levels),
    R = factor(rep(levels, each = n / length(levels)), levels = levels),
    rt = round(runif(n, .3, 1.2), 3)
  )
}

df_matchfun <- function(d) as.character(d$S) == as.character(d$lR)

df_formula <- function(correlated = TRUE) {
  out <- list(v ~ 0 + lR, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1, sv ~ 1,
              mG ~ 1, mK ~ 1)
  if (correlated) out$rho <- rho ~ 1
  out
}

df_design <- function(dat, model, correlated = TRUE) {
  design(data = dat, Rlevels = levels(dat$R), matchfun = df_matchfun,
         formula = df_formula(correlated), constants = c(s = 0, mG = 0, mK = 0),
         model = model, report_p_vector = FALSE)
}

# v is on the log scale under posdrift = TRUE and the natural scale otherwise;
# see .rdmswtn_v_scale().
df_pars <- function(des, v, sv = .7, rho = NULL, posdrift = TRUE,
                    B = 1, A = .2, t0 = .15) {
  p <- sampled_pars(des, doMap = FALSE)
  p[grep("^v_", names(p))] <- if (posdrift) log(v) else v
  p["B"] <- log(B); p["A"] <- log(A); p["t0"] <- log(t0); p["sv"] <- log(sv)
  if (!is.null(rho)) p["rho"] <- qnorm((rho + 1) / 2)
  p
}

df_ll <- function(dat, des, p, min_ll = -Inf) {
  emc <- make_emc(dat, des, type = "single", n_chains = 1, compress = FALSE,
                  verbose = FALSE, rt_resolution = NULL)[[1]]
  dadm <- emc$data[[1]]
  m <- emc$model()
  designs <- lapply(names(m$p_types), function(nm) {
    dm <- attr(dadm, "designs")[[nm]]
    dm[attr(dm, "expand"), , drop = FALSE]
  })
  names(designs) <- names(m$p_types)
  EMC2:::calc_ll_oo(
    particle_matrix = matrix(p, nrow = 1, dimnames = list(NULL, names(p))),
    data = dadm, constants = attr(dadm, "constants"), designs = designs,
    type = m$c_name, bounds = m$bound, transforms = m$transform,
    pretransforms = m$pre_transform, p_types = names(m$p_types),
    min_ll = min_ll, trend = m$trend, marginalise = NULL)
}

test_that("the factor map reproduces the requested marginals and correlation", {
  # v_q = v + sign(rho) sv sqrt(|rho|) z with sv_q = sv sqrt(1 - |rho|) leaves
  # each marginal SD at sv and puts the product of the two loadings at rho.
  for (rho in c(.8, -.55, 0)) {
    sv <- c(0.7, 1.3)
    mag <- abs(rho)
    slope <- c(1, sign(rho)) * sv * sqrt(mag)   # correct row is the reference
    resid <- sv * sqrt(1 - mag)
    expect_equal(sqrt(slope^2 + resid^2), sv, tolerance = 1e-12)
    expect_equal(slope[1] * slope[2] / prod(sv), rho, tolerance = 1e-12)
  }
})

test_that("correlate is validated and selects the compiled model name", {
  expect_error(RDMSWTNcorr(correlate = "ranks"), "should be one of")
  expect_identical(RDMSWTNcorr()$c_name, "RDMSWTN_CORR")
  expect_identical(RDMSWTNcorr(correlate = "drifts")$c_name, "RDMSWTN_CORRD")
  expect_identical(RDMSWTNcorr(posdrift = FALSE, correlate = "drifts")$c_name,
                   "RDMSWTN_IO_CORRD")
  expect_identical(RDMSWTN_TTcorr(correlate = "drifts")$c_name,
                   "RDMSWTN_TT_CORRD")
  expect_identical(RDMSWTNcorr(correlate = "drifts")$correlation_type,
                   "rdmswtn_drift_factor")
  expect_identical(RDMSWTNcorr()$correlation_type,
                   "rdmswtn_gaussian_copula")
})

test_that("the drifts route requires sv > 0 on the correlated rows", {
  model <- RDMSWTNcorr(correlate = "drifts")
  # One trial's worth of natural-scale rows, in p_types order, as the other
  # RDMSWTN Ttransform tests build them.
  pars <- matrix(rep(c(1, 1, .2, .1, 1, .4, 1, 1, 0, 0, .5), each = 2),
                 nrow = 2, dimnames = list(NULL, names(model$p_types)))
  dadm <- data.frame(lM = factor(c("TRUE", "FALSE"),
                                 levels = c("FALSE", "TRUE")))
  expect_silent(model$Ttransform(pars, dadm))

  pars0 <- pars
  pars0[, "sv"] <- 0
  expect_error(model$Ttransform(pars0, dadm), "requires sv > 0")

  # rho = 0 rows are uncoupled, so sv = 0 is still legal there.
  pars_indep <- pars0
  pars_indep[, "rho"] <- 0
  expect_silent(model$Ttransform(pars_indep, dadm))

  # The correct racer is the positive reference; the incorrect racer carries
  # the cell sign, so their loadings multiply back to rho.
  neg <- pars
  neg[, "rho"] <- -.5
  expect_equal(unname(model$Ttransform(neg, dadm)[, "rho"]), c(.5, -.5))
})

test_that("the drifts route needs lM from a matchfun", {
  model <- RDMSWTNcorr(correlate = "drifts")
  pars <- matrix(c(1, 1, .2, .1, 1, .4, 1, 1, 0, 0, .5), nrow = 1,
                 dimnames = list(NULL, names(model$p_types)))
  expect_error(model$Ttransform(pars, NULL), "lM role indicator")
})

test_that("rho = 0 on the drifts route is the independent RDMSWTN likelihood", {
  dat <- df_data()
  for (posdrift in c(TRUE, FALSE)) {
    v <- if (posdrift) c(2, 1.2) else c(.8, -.3)
    des_i <- df_design(dat, RDMSWTN(posdrift = posdrift), correlated = FALSE)
    des_c <- df_design(dat, RDMSWTNcorr(posdrift = posdrift,
                                        correlate = "drifts"))
    indep <- df_ll(dat, des_i, df_pars(des_i, v, posdrift = posdrift))
    zero <- df_ll(dat, des_c,
                  df_pars(des_c, v, rho = 0, posdrift = posdrift))
    expect_equal(zero, indep, tolerance = 1e-12)
  }
})

test_that("the drifts likelihood matches an independent factor quadrature", {
  skip_if_not_installed("statmod")
  dat <- df_data()
  b <- 1.2; A <- .2; t0 <- .15; sv <- .7
  gh <- statmod::gauss.quad(160, kind = "hermite")
  zs <- sqrt(2) * gh$nodes
  ws <- gh$weights / sqrt(pi)
  kernel_pars <- function(vq, sv_res) {
    matrix(c(vq, b, A, t0, 1, sv_res, 0, 0), nrow = 1,
           dimnames = list(NULL, c("v", "b", "A", "t0", "s", "sv",
                                    "lambda_g", "lambda_k")))
  }

  # Conditional on a node the two accumulators are independent RDMSWTN
  # kernels at (v_q, sv_res).  Under posdrift the scalar kernel returns the
  # truncation-normalised quantity, so each row's Phi(v_q / sv_res) is
  # restored and the trial is divided by the joint orthant probability.
  reference <- function(v, rho, posdrift) {
    mag <- abs(rho)
    sv_res <- sv * sqrt(1 - mag)
    logZ <- if (!posdrift) 0 else log(mvtnorm::pmvnorm(
      lower = c(0, 0), mean = v,
      sigma = matrix(c(sv^2, rho * sv^2, rho * sv^2, sv^2), 2))[1])
    total <- 0
    for (k in seq_len(nrow(dat))) {
      win <- if (dat$R[k] == "a") 1L else 2L
      lose <- 3L - win
      correct <- if (dat$S[k] == "a") 1L else 2L
      sgn <- ifelse(seq_len(2) == correct, 1, sign(rho))
      slope <- sgn * sv * sqrt(mag)
      f <- vapply(zs, function(z) {
        vq <- v + slope * z
        d <- EMC2:::dRDMSWTN(dat$rt[k], kernel_pars(vq[win], sv_res),
                             posdrift = posdrift)
        s <- 1 - EMC2:::pRDMSWTN(dat$rt[k], kernel_pars(vq[lose], sv_res),
                                 posdrift = posdrift)
        if (posdrift) d <- d * pnorm(vq[win] / sv_res)
        if (posdrift) s <- s * pnorm(vq[lose] / sv_res)
        d * s
      }, 0)
      total <- total + log(sum(ws * f)) - logZ
    }
    total
  }

  # The shared-factor integral is a finite Gauss-Hermite rule, so agreement is
  # limited by its node count rather than by machine precision; see the
  # scan_default/fine_default tiers in c_log_likelihood_corr_drift().
  cases <- list(list(TRUE, c(2, 1.2), .5), list(TRUE, c(2, 1.2), -.6),
                list(TRUE, c(.9, .5), .5), list(FALSE, c(.8, -.3), .5),
                list(FALSE, c(.8, -.3), -.6))
  for (case in cases) {
    posdrift <- case[[1]]; v <- case[[2]]; rho <- case[[3]]
    des <- df_design(dat, RDMSWTNcorr(posdrift = posdrift,
                                      correlate = "drifts"))
    got <- df_ll(dat, des, df_pars(des, v, sv = sv, rho = rho,
                                   posdrift = posdrift))
    expect_equal(got, reference(v, rho, posdrift), tolerance = 2e-3,
                 info = sprintf("posdrift=%s rho=%+.2f", posdrift, rho))
  }
})

test_that("the drifts route reaches the shared generic node evaluator", {
  skip_if(!nzchar(Sys.getenv("EMC2_BAWLCORR_COUNTERS")),
          "route counters are only collected when EMC2_BAWLCORR_COUNTERS is set")
  dat <- df_data()
  des <- df_design(dat, RDMSWTNcorr(correlate = "drifts"))
  before <- EMC2:::bawl_corr_counter_values()$gh_generic_clock_trials
  df_ll(dat, des, df_pars(des, c(2, 1.2), rho = .5))
  after <- EMC2:::bawl_corr_counter_values()$gh_generic_clock_trials
  # A Wald kernel has no affine-in-drift survivor, so it must never take the
  # exact rectangle or fused no-clock routes even on a two-accumulator trial.
  expect_gt(after - before, 0)
})

test_that("the correlated-draw simulator reproduces its own likelihood", {
  dat0 <- data.frame(
    subjects = factor(rep(1, 40)),
    S = factor(rep(c("a", "b"), length.out = 40), levels = c("a", "b")),
    R = factor(rep(c("a", "b"), each = 20), levels = c("a", "b")),
    rt = .5)
  for (posdrift in c(TRUE, FALSE)) {
    for (rho_true in c(.7, -.6)) {
      des <- df_design(dat0, RDMSWTNcorr(posdrift = posdrift,
                                         correlate = "drifts"))
      v <- if (posdrift) c(.9, .5) else c(1, .2)
      p <- df_pars(des, v, sv = .7, rho = rho_true, posdrift = posdrift)
      set.seed(99)
      sim <- make_data(p, design = des, n_trials = 2500)
      grid <- seq(-.9, .9, by = .15)
      lls <- vapply(grid, function(r) {
        q <- p
        q["rho"] <- qnorm((r + 1) / 2)
        df_ll(sim, des, q)
      }, 0)
      # One grid step of slack: the profile peak is estimated from a finite
      # sample, and the correlation is the weakest-identified parameter here.
      expect_lt(abs(grid[which.max(lls)] - rho_true), .16,
                label = sprintf("posdrift=%s rho=%+.2f argmax",
                                posdrift, rho_true))
    }
  }
})

test_that("the R and compiled correlated-draw simulators agree in distribution", {
  dat0 <- data.frame(
    subjects = factor(rep(1, 2), levels = "1"),
    S = factor(c("a", "b"), levels = c("a", "b")),
    R = factor(c("a", "b"), levels = c("a", "b")),
    rt = .5)
  des <- df_design(dat0, RDMSWTNcorr(correlate = "drifts"))
  p <- df_pars(des, c(1.5, .9), sv = .7, rho = .6)
  draw <- function(use_cpp) {
    withr::with_options(list(emc2.cpp_rfun = use_cpp), {
      set.seed(4)
      make_data(p, design = des, n_trials = 6000)
    })
  }
  a <- draw(TRUE)
  b <- draw(FALSE)
  expect_equal(mean(a$R == "a"), mean(b$R == "a"), tolerance = .02)
  expect_equal(mean(a$rt, na.rm = TRUE), mean(b$rt, na.rm = TRUE),
               tolerance = .02)
})
