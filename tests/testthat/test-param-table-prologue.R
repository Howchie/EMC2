# The per-particle parameter prologue -- refill from the particle row, map the
# designs, apply the transforms, check the bounds -- has no direct test: it is
# covered only transitively, through likelihood values.  It is also the code
# most exposed to the row-constant / planned / invariant short cuts, whose whole
# point is to produce identical numbers by a cheaper route.
#
# EMC2:::pt_prologue_oo() runs exactly the lane calc_ll_oo uses.
# EMC2:::get_pars_c_batch_wrapper_oo() is an independent implementation of the
# same mapping that takes none of those short cuts.  Every test below pins the
# first against the second, or against a hand-computed value.
#
# When changing that code, also run the assertion build, which re-scans every
# column whose row-constant flag is about to be acted on:
#   EMC2_EXTRA_CPPFLAGS=-DEMC2_PT_CONST_CHECK R CMD INSTALL --library=<lib> .

RNGkind("L'Ecuyer-CMRG")

# --- fixtures ---------------------------------------------------------------

prologue_fixture <- function(formula, n_trials = 60, n_particles = 7,
                             covariates = NULL, cov_data = NULL,
                             constants = c(sv = log(1)), trend = NULL,
                             pre_transform_terms = NULL, seed = 11) {
  set.seed(seed)
  dat <- forstmann[forstmann$subjects == levels(forstmann$subjects)[1L], ]
  dat <- dat[rep(seq_len(nrow(dat)), length.out = n_trials), , drop = FALSE]
  dat$subjects <- factor("s1")
  dat$rt <- dat$rt + runif(nrow(dat), 0, 0.3)
  if (!is.null(cov_data)) for (nm in names(cov_data)) dat[[nm]] <- cov_data[[nm]]
  rownames(dat) <- NULL

  args <- list(data = dat, model = LBA, matchfun = function(d) d$S == d$lR,
               formula = formula, constants = constants)
  if (!is.null(covariates)) args$covariates <- covariates
  if (!is.null(trend)) args$trend <- trend
  if (!is.null(pre_transform_terms)) args$pre_transform_terms <- pre_transform_terms
  des <- suppressMessages(do.call(design, args))

  emc <- suppressMessages(make_emc(dat, des, n_chains = 1, compress = TRUE,
                                   rt_resolution = NULL, type = "single"))
  dadm <- emc[[1]]$data[[1]]
  model <- des$model()
  pars <- sampled_pars(des)
  # Realistic values: at the sampled_pars() zeros t0 = 1 s floors every row and
  # the fixture stops discriminating anything.
  centre <- pars; centre[] <- 0
  for (nm in names(centre)) {
    centre[nm] <- switch(sub("_.*$", "", nm),
                         v = 1.5, B = log(1), A = log(0.3), t0 = log(0.15), 0.1)
  }
  centre[grepl("_", names(centre))] <- 0.15
  prop <- matrix(rnorm(n_particles * length(pars), 0, 0.05), n_particles,
                 length(pars), dimnames = list(NULL, names(pars)))
  prop <- sweep(prop, 2, as.numeric(centre), "+")

  list(dadm = dadm, model = model, prop = prop, des = des,
       p_types = names(model$p_types),
       constants = { cst <- attr(dadm, "constants"); if (is.null(cst)) NA else cst })
}

# The prologue as calc_ll_oo runs it.
run_prologue <- function(fx, designs = NULL) {
  if (is.null(designs)) designs <- EMC2:::.oo_expanded_designs(fx$dadm, expand = FALSE)
  EMC2:::pt_prologue_oo(
    fx$prop, fx$dadm, constants = fx$constants, designs = designs,
    bounds = fx$model$bound, transforms = fx$model$transform,
    pretransforms = fx$model$pre_transform, p_types = fx$p_types,
    trend = fx$model$trend, return_pars = TRUE)
}

# The same mapping through the general wrapper, which uses neither the planned
# lane, nor the invariant lane, nor row-constant columns.
run_reference <- function(fx, designs = NULL) {
  if (is.null(designs)) designs <- EMC2:::.oo_expanded_designs(fx$dadm, expand = FALSE)
  ref <- EMC2:::get_pars_c_batch_wrapper_oo(
    fx$prop, fx$dadm, constants = fx$constants, designs = designs,
    bounds = fx$model$bound, transforms = fx$model$transform,
    pretransforms = fx$model$pre_transform, trend = fx$model$trend)
  ref[, fx$p_types, , drop = FALSE]
}

expect_prologue_matches_reference <- function(fx, designs = NULL, label = "") {
  got <- run_prologue(fx, designs)$pars
  want <- run_reference(fx, designs)
  expect_equal(dim(got), dim(want), info = label)
  expect_equal(as.numeric(got), as.numeric(want), tolerance = 0, info = label)
}

# --- the design shapes ------------------------------------------------------

test_that("mapping, transforms and bounds match the reference lane on every design shape", {
  shapes <- list(
    "intercept only"   = list(v ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
    "two-level factor" = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1),
    "three-level factor" = list(v ~ 1, B ~ E, A ~ 1, t0 ~ 1),
    "interaction"      = list(v ~ lM * E, B ~ E * lR, A ~ 1, t0 ~ 1),
    "no intercept"     = list(v ~ 0 + lM, B ~ 0 + lR, A ~ 1, t0 ~ 1)
  )
  for (nm in names(shapes)) {
    fx <- prologue_fixture(shapes[[nm]])
    expect_prologue_matches_reference(fx, label = nm)
    # And again with the design matrices expanded to trial resolution, which
    # takes the identity-expand branch of the mapping instead.
    expanded <- EMC2:::.oo_expanded_designs(fx$dadm, expand = TRUE)
    expect_prologue_matches_reference(fx, expanded, label = paste(nm, "(expanded)"))
  }
})

test_that("deferred scalar markers are cleared by varying mapped writers", {
  old_defer <- EMC2:::emc_pt_defer_scalars()
  on.exit(EMC2:::emc_pt_defer_scalars(old_defer), add = TRUE)
  EMC2:::emc_pt_defer_scalars(TRUE)

  # B is first mapped across lM.  v then reads the already-varying B column
  # through a one-cell, non-self design (the low-level shape of v ~ B).
  # The initial particle refill defers v because the structural v design is
  # single-cell.  Mapping v takes the row route because B is varying, so the
  # writer must retire v's deferred scalar marker before materialize().
  fx <- prologue_fixture(list(v ~ 1, B ~ lM, A ~ 1, t0 ~ 1),
                         n_trials = 37, n_particles = 2)
  designs <- EMC2:::.oo_expanded_designs(fx$dadm, expand = FALSE)
  v_design <- matrix(1, nrow = 1, ncol = 1,
                     dimnames = list(NULL, "B"))
  attr(v_design, "expand") <- rep.int(1L, nrow(fx$dadm))
  designs$v <- v_design
  designs <- designs[c("B", "v", setdiff(names(designs), c("B", "v")))]

  got <- run_prologue(fx, designs)$pars
  want <- run_reference(fx, designs)
  expect_equal(as.numeric(got[, "v", 2]), as.numeric(want[, "v", 2]),
               tolerance = 0)
  expect_gt(length(unique(got[, "v", 2])), 1)
})

test_that("a continuous covariate is not treated as a single design cell", {
  # This is the one input shape that would silently break the row-constant
  # short cut: the design has one row per distinct covariate value, so the
  # output column varies by trial even though the coefficient does not.
  set.seed(5)
  cov1 <- rnorm(60)
  fx <- prologue_fixture(list(v ~ cov1, B ~ 1, A ~ 1, t0 ~ 1),
                         covariates = "cov1", cov_data = list(cov1 = cov1))
  expect_prologue_matches_reference(fx, label = "covariate")
  got <- run_prologue(fx)$pars
  # v really does vary across trials -- otherwise this test proves nothing.
  expect_gt(length(unique(round(got[, "v", 1], 10))), 2)
})

test_that("split (pre-sum) transforms map identically", {
  fx <- prologue_fixture(list(v ~ E, B ~ 1, A ~ 1, t0 ~ 1),
                         pre_transform_terms = list(v = c("v", "v_Eneutral")))
  expect_prologue_matches_reference(fx, label = "split transform")
})

test_that("constants reach the parameter table unchanged", {
  fx <- prologue_fixture(list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1),
                         constants = c(sv = log(2)))
  got <- run_prologue(fx)$pars
  expect_true(all(got[, "sv", ] == 2))
})

# --- bounds -----------------------------------------------------------------

test_that("bound verdicts are per row, not per column", {
  # B ~ E with one level pushed below its lower bound: the trials in that cell
  # must be rejected and the others kept.  A row-constant short cut that
  # collapsed the column would fail this in whichever direction it got wrong.
  fx <- prologue_fixture(list(v ~ 1, B ~ E, A ~ 1, t0 ~ 1), n_particles = 1)
  ok_all <- run_prologue(fx)$ok
  expect_true(all(ok_all))

  bad <- fx
  bad$prop[1, "t0"] <- log(1e-6)         # t0 below its lower bound (0.05) everywhere
  expect_equal(sum(run_prologue(bad)$ok), 0L)

  part <- fx
  # B_Eneutral large and negative makes B negative in the neutral cells only.
  part$prop[1, "B_Eneutral"] <- -1e6
  ok_part <- run_prologue(part)$ok
  expect_gt(sum(ok_part), 0L)
  expect_lt(sum(ok_part), length(ok_part))
  neutral <- fx$dadm$E == "neutral"
  expect_true(all(!ok_part[neutral]))
  expect_true(all(ok_part[!neutral]))
})

test_that("an invariant bound is still enforced for every particle", {
  # sv is a constant, so its bound check is hoisted out of the particle loop.
  fx <- prologue_fixture(list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1),
                         constants = c(sv = log(1)))
  expect_true(all(run_prologue(fx)$ok))
})

# --- the trend lane ---------------------------------------------------------

test_that("an identity trend gives the same likelihood as no trend", {
  # The row-constant short cut is disabled whenever a TrendRuntime is present,
  # because the trend writes into the parameter table outside the pipeline.
  # A trend whose kernel is the identity therefore checks the gated and
  # un-gated lanes against each other.
  set.seed(3)
  n_trials <- 40
  dat <- forstmann[forstmann$subjects == levels(forstmann$subjects)[1L], ]
  dat <- dat[rep(seq_len(nrow(dat)), length.out = n_trials), , drop = FALSE]
  dat$subjects <- factor("s1")
  dat$rt <- dat$rt + runif(nrow(dat), 0, 0.3)
  dat$flat <- 0                      # covariate with no variation
  rownames(dat) <- NULL

  base_des <- suppressMessages(design(
    data = dat, model = LBA, matchfun = function(d) d$S == d$lR,
    formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1), constants = c(sv = log(1))))

  trend <- make_trend(par_names = "B", cov_names = "flat", kernels = "lin_incr")
  trend_des <- suppressMessages(design(
    data = dat, model = LBA, matchfun = function(d) d$S == d$lR,
    covariates = "flat", trend = trend,
    formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1), constants = c(sv = log(1))))

  base_emc <- suppressMessages(make_emc(dat, base_des, n_chains = 1,
                                        rt_resolution = NULL, type = "single"))
  trend_emc <- suppressMessages(make_emc(dat, trend_des, n_chains = 1,
                                         rt_resolution = NULL, type = "single"))

  bp <- sampled_pars(base_des); bp[] <- 0
  bp["v"] <- 1.5; bp["v_lMTRUE"] <- 0.5; bp["B"] <- log(1)
  bp["A"] <- log(0.3); bp["t0"] <- log(0.15)
  tp <- sampled_pars(trend_des); tp[] <- 0
  tp[names(bp)] <- bp
  # A linear trend on a covariate that is identically zero adds nothing.
  tp[setdiff(names(tp), names(bp))] <- 0

  base_ll <- EMC2:::calc_ll_manager(matrix(bp, 1, dimnames = list(NULL, names(bp))),
                                    base_emc[[1]]$data[[1]], base_des$model, r_cores = 1)
  trend_ll <- EMC2:::calc_ll_manager(matrix(tp, 1, dimnames = list(NULL, names(tp))),
                                     trend_emc[[1]]$data[[1]], trend_des$model, r_cores = 1)
  expect_equal(trend_ll, base_ll, tolerance = 1e-10)
})

# --- particle independence --------------------------------------------------

test_that("a particle's parameters do not depend on how many particles precede it", {
  # The planned / invariant / row-constant lanes all carry state from the
  # template particle forward.  Mapping one particle alone must give the same
  # answer as mapping it inside a batch.
  fx <- prologue_fixture(list(v ~ lM * E, B ~ E, A ~ 1, t0 ~ 1), n_particles = 6)
  batch <- run_prologue(fx)$pars
  for (i in seq_len(nrow(fx$prop))) {
    one <- fx
    one$prop <- fx$prop[i, , drop = FALSE]
    expect_equal(as.numeric(run_prologue(one)$pars),
                 as.numeric(batch[, , i]), tolerance = 0,
                 info = paste("particle", i))
  }
})

# --- route diagnostics ------------------------------------------------------

test_that("mapper route counters distinguish cell and row designs", {
  # Diagnostics are process-local and off by default.  Restore the switch even
  # when an assertion fails so later test files cannot inherit instrumentation.
  old <- EMC2:::emc_pt_mapper_stats()
  on.exit(EMC2:::emc_pt_mapper_stats(enabled = old$enabled, reset = TRUE), add = TRUE)
  old_reuse <- EMC2:::emc_pt_cell_min_reuse()
  on.exit(EMC2:::emc_pt_cell_min_reuse(old_reuse), add = TRUE)


  EMC2:::emc_pt_mapper_stats(enabled = TRUE, reset = TRUE)
  cell_fx <- prologue_fixture(list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1),
                              n_trials = 60, n_particles = 2)
  run_prologue(cell_fx)
  cell_stats <- EMC2:::emc_pt_mapper_stats()
  expect_gt(cell_stats$map_cell, 0)
  expect_gt(cell_stats$cell_admit, 0)

  EMC2:::emc_pt_mapper_stats(reset = TRUE)
  # Force the general row route for the continuous design: it is a deliberate
  # route-cliff probe, not an assertion that continuous cells are invalid.
  EMC2:::emc_pt_cell_min_reuse(1)
  cov1 <- seq(-1, 1, length.out = 60)
  row_fx <- prologue_fixture(list(v ~ cov1, B ~ 1, A ~ 1, t0 ~ 1),
                             n_trials = 60, n_particles = 2,
                             covariates = "cov1",
                             cov_data = list(cov1 = cov1))
  run_prologue(row_fx)
  row_stats <- EMC2:::emc_pt_mapper_stats()
  expect_gt(row_stats$map_row, 0)
  expect_equal(row_stats$cell_admit, 0)
})

test_that("mapper counters stay inert while disabled", {
  EMC2:::emc_pt_mapper_stats(enabled = FALSE, reset = TRUE)
  fx <- prologue_fixture(list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1),
                         n_trials = 20, n_particles = 1)
  run_prologue(fx)
  stats <- EMC2:::emc_pt_mapper_stats()
  expect_false(stats$enabled)
  expect_equal(unname(unlist(stats[setdiff(names(stats), "enabled")])), rep(0, 13))
})
