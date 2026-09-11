# Shared fixtures for the architecture-efficiency audit's acceptance gates.
#
# The audit's "Numerical and portability acceptance gates" matrix is demanded by
# nearly every commit in its sequence.  Building it inside whichever commit
# lands first would mean the second commit inherits a harness shaped around the
# first one's change, so it lives here, before any of them.
#
# Three groups of callers:
#
#   test-audit-mapper-equivalence.R  mapping/transforms/bounds against the
#                                    independent get_pars_c_batch_wrapper_oo
#                                    route, across design shapes and cell counts
#   test-audit-differential.R        likelihood and stream equivalence under
#                                    re-ordering, re-batching and pool failure
#   test-audit-numerical-grid.R      the awkward corners of the support
#   test-audit-faults.R              worker lifecycle and protocol faults
#
# Everything here is a fixture or an assertion; nothing measures time.  Timing
# belongs in WorkingTests/, which reports through R/profile_schema.R.

# --- data -------------------------------------------------------------------

# One participant's empirical data, repeated to the requested length and
# jittered so RT compression does not silently collapse the fixture to a
# handful of unique rows.  `cells` gives a numeric covariate exactly that many
# distinct values, which is how a design is driven to a chosen number of design
# rows without widening the parameter vector.  `width` instead widens the
# parameter vector through a factor.
audit_data <- function(n_trials = 400L, cells = NULL, width = 1L,
                       n_subjects = 1L, seed = 91026L, jitter = 0.1) {
  set.seed(seed)
  base <- forstmann[forstmann$subjects == levels(forstmann$subjects)[1L], ]
  base <- base[rep(seq_len(nrow(base)), length.out = n_trials), , drop = FALSE]
  out <- do.call(rbind, lapply(seq_len(n_subjects), function(s) {
    d <- base
    d$subjects <- sprintf("s%02d", s)
    d
  }))
  out$subjects <- factor(out$subjects)
  if (jitter > 0) out$rt <- out$rt + runif(nrow(out), 0, jitter)
  out$condition <- factor(rep_len(seq_len(width), nrow(out)))
  out$x <- rep_len(seq_len(if (is.null(cells)) 10L else cells), nrow(out)) / 1000
  rownames(out) <- NULL
  out
}

# --- fixtures ---------------------------------------------------------------

# Per-model formulas and the constants each one needs to be identified, so a
# test can name a model rather than restating its whole specification.
.audit_model_spec <- function(model, cells = NULL, width = 1L) {
  vform <- if (!is.null(cells)) v ~ x else if (width > 1L) v ~ 0 + condition else v ~ lM
  switch(model,
    DDM = list(
      formula = list(if (!is.null(cells) || width > 1L) vform else v ~ S,
                     a ~ E, t0 ~ 1),
      constants = c(s = log(1)),
      centre = c(v = 2, a = log(1.2), t0 = log(0.15))),
    RDM = list(
      formula = list(vform, B ~ E, A ~ 1, t0 ~ 1),
      constants = c(s = 0),
      centre = c(v = log(2), B = log(1), A = log(0.3), t0 = log(0.15))),
    LBA = list(
      formula = list(vform, B ~ E, A ~ 1, t0 ~ 1),
      constants = c(sv = log(1)),
      centre = c(v = 2, B = log(1), A = log(0.3), t0 = log(0.15))),
    BAwL = list(
      formula = list(vform, B ~ E, A ~ 1, t0 ~ 1, k ~ 1),
      constants = c(sv = 0),
      centre = c(v = 2, B = log(1), A = log(0.3), t0 = log(0.15), k = log(0.2))),
    stop("no audit spec for model ", model))
}

# A complete likelihood fixture: compressed data, the model closure, and a
# particle cloud centred somewhere the likelihood is actually finite.  At
# sampled_pars()' zeros t0 is 1 s, which floors every row and leaves the
# fixture unable to discriminate anything at all.
audit_fixture <- function(model = "RDM", n_trials = 400L, n_particles = 5L,
                          cells = NULL, width = 1L, formula = NULL,
                          constants = NULL, data = NULL, seed = 91026L,
                          spread = 0.04, rt_resolution = NULL, ...) {
  spec <- .audit_model_spec(model, cells, width)
  dat <- if (is.null(data)) audit_data(n_trials, cells, width, seed = seed) else data
  model_fun <- get(model, envir = asNamespace("EMC2"))
  des <- suppressMessages(design(
    data = dat, model = model_fun,
    formula = if (is.null(formula)) spec$formula else formula,
    matchfun = function(z) z$S == z$lR, report_p_vector = FALSE,
    constants = if (is.null(constants)) spec$constants else constants, ...))
  emc <- suppressMessages(make_emc(dat, des, n_chains = 1L, compress = TRUE,
                                   rt_resolution = rt_resolution, type = "single"))
  pars <- sampled_pars(des)
  centre <- pars
  centre[] <- 0
  for (nm in names(spec$centre)) {
    hit <- names(centre) == nm | startsWith(names(centre), paste0(nm, "_"))
    centre[hit] <- spec$centre[[nm]]
  }
  # Contrast coefficients start small rather than at the intercept's value, or
  # a wide design walks every cell far from where the density is finite.
  centre[grepl("_", names(centre), fixed = TRUE)] <- 0.1
  set.seed(seed + 1L)
  prop <- matrix(rnorm(n_particles * length(pars), 0, spread), n_particles,
                 length(pars), dimnames = list(NULL, names(pars)))
  prop <- sweep(prop, 2L, as.numeric(centre), "+")
  dadm <- emc[[1L]]$data[[1L]]
  cst <- attr(dadm, "constants")
  list(model_name = model, dadm = dadm, model = des$model, des = des,
       data = dat, emc = emc, prop = prop, centre = centre,
       p_types = names(des$model()$p_types),
       constants = if (is.null(cst)) NA else cst)
}

# --- the two mapping routes -------------------------------------------------

audit_designs <- function(fx, expand = FALSE) {
  EMC2:::.oo_expanded_designs(fx$dadm, expand = expand)
}

audit_ll_args <- function(fx, designs = NULL) {
  m <- fx$model()
  list(particle_matrix = fx$prop, data = fx$dadm, constants = fx$constants,
       designs = if (is.null(designs)) audit_designs(fx) else designs,
       bounds = m$bound, transforms = m$transform,
       pretransforms = m$pre_transform, p_types = fx$p_types, trend = m$trend)
}

# The prologue as calc_ll_oo runs it: planned lane, invariant columns,
# row-constant and design-cell short cuts all live.
audit_prologue <- function(fx, designs = NULL) {
  do.call(EMC2:::pt_prologue_oo, c(audit_ll_args(fx, designs),
                                   list(return_pars = TRUE)))
}

# The same mapping through the general wrapper, which takes none of them.
#
# The wrapper returns the model's kernel parameters and no more.  A trend adds
# weight columns (`B.w` and friends) to the prologue's `p_types` that the
# wrapper never emits, so the comparison is over the columns the two routes have
# in common -- which is still every parameter the likelihood consumes, including
# the ones the trend has already modified.
audit_reference <- function(fx, designs = NULL) {
  m <- fx$model()
  ref <- EMC2:::get_pars_c_batch_wrapper_oo(
    fx$prop, fx$dadm, constants = fx$constants,
    designs = if (is.null(designs)) audit_designs(fx) else designs,
    bounds = m$bound, transforms = m$transform,
    pretransforms = m$pre_transform, trend = m$trend)
  ref[, audit_shared_pars(fx, ref), , drop = FALSE]
}

audit_shared_pars <- function(fx, ref) {
  shared <- intersect(fx$p_types, dimnames(ref)[[2L]])
  if (!length(shared)) stop("the two mapping routes share no parameter columns")
  shared
}

audit_ll_direct <- function(fx, designs = NULL, min_ll = log(1e-10)) {
  m <- fx$model()
  do.call(EMC2:::calc_ll_oo,
          c(audit_ll_args(fx, designs), list(type = m$c_name, min_ll = min_ll)))
}

audit_ll_managed <- function(fx, r_cores = 1L) {
  EMC2:::calc_ll_manager(fx$prop, fx$dadm, fx$model, r_cores = r_cores)
}

# --- the cell-count cutoff probe -------------------------------------------

# Append `n_extra` copies of an existing design row to one parameter's design
# matrix.  The extra rows are unreferenced by `expand`, so every mapped value,
# every trial and every particle is unchanged -- only the *number of cells*
# moves.  That is what makes a cell-count cutoff measurable and testable
# without confounding it with a change in the model.
audit_pad_designs <- function(args, par, n_extra = 1L) {
  dm <- args$designs[[par]]
  keep <- attributes(dm)
  padded <- rbind(dm, dm[rep(1L, n_extra), , drop = FALSE])
  for (nm in setdiff(names(keep), c("dim", "dimnames"))) {
    attr(padded, nm) <- keep[[nm]]
  }
  args$designs[[par]] <- padded
  args
}

audit_call_ll <- function(fx, args, min_ll = log(1e-10)) {
  do.call(EMC2:::calc_ll_oo, c(args, list(type = fx$model()$c_name, min_ll = min_ll)))
}

# --- assertions -------------------------------------------------------------

# Bit-identical, not "close".  The audit's gate for scheduling, storage and
# cache changes is exact reproduction: a tolerance here would accept exactly the
# kind of silent reassociation the gate exists to catch.
expect_bit_identical <- function(got, want, label = "") {
  testthat::expect_equal(length(got), length(want), info = label)
  same <- identical(as.numeric(got), as.numeric(want))
  if (!same) {
    d <- max(abs(as.numeric(got) - as.numeric(want)))
    testthat::fail(sprintf("%s: not bit-identical, max |delta| = %.17g", label, d))
  } else {
    testthat::succeed()
  }
  invisible(same)
}

# A fixture whose likelihood is only the invalid-parameter floor proves nothing:
# every route agrees on a constant.  The bench script rejects those and so must
# the tests.
expect_informative_ll <- function(fx, ll, label = "") {
  n <- length(attr(fx$dadm, "expand"))
  if (!n) n <- nrow(fx$dadm)
  floor_ll <- n * log(1e-10)
  testthat::expect_true(all(is.finite(ll)), info = paste(label, "- all finite"))
  testthat::expect_true(any(abs(ll - floor_ll) > 1),
                        info = paste(label, "- not merely the invalid floor"))
}

# --- worker-pool fakes ------------------------------------------------------

# A pool object the lifecycle helpers accept, with no real processes behind it.
# The audit is explicit that the worker-ID boundary must be tested with a small
# simulated ID space rather than by starting 256 processes: the protocol
# question is how a width is encoded, and forking 256 R sessions to ask it
# would make the test unrunnable everywhere it matters.
audit_fake_pool <- function(n = 4L, alive = TRUE, dir = NULL) {
  list(n = as.integer(n), dir = dir, jobs = vector("list", n),
       wcs = vector("list", n), rcs = vector("list", n),
       done = NULL, alive = alive, backend = "fork")
}

# Round-trip a value through the pool's own framing over a real file, which is
# where short reads and partial payloads actually happen.
audit_frame_roundtrip <- function(obj) {
  f <- tempfile()
  on.exit(unlink(f), add = TRUE)
  con <- file(f, "wb")
  ok <- EMC2:::.emc_wpool_send(con, obj)
  close(con)
  if (!isTRUE(ok)) return(list(ok = FALSE, value = NULL))
  con <- file(f, "rb")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  list(ok = TRUE, value = EMC2:::.emc_wpool_recv(con), bytes = file.size(f))
}

# Truncate a framed message to `keep` bytes: a payload that starts arriving and
# then stops, which is the shape of a worker dying mid-reply.
audit_truncated_frame <- function(obj, keep) {
  f <- tempfile()
  con <- file(f, "wb")
  EMC2:::.emc_wpool_send(con, obj)
  close(con)
  raw_all <- readBin(f, "raw", file.size(f))
  unlink(f)
  g <- tempfile()
  writeBin(raw_all[seq_len(min(keep, length(raw_all)))], g)
  g
}

# --- joint cell partitions --------------------------------------------------

# Build a ParamTable from a fixture the way the prologue does, map its designs,
# and ask for the common refinement of a set of parameter columns.  The table is
# built here rather than reached into, because the partition depends only on the
# designs and the data and is meant to be answerable without running a particle.
audit_joint_cells <- function(fx, params, designs = NULL) {
  des <- if (is.null(designs)) audit_designs(fx) else designs
  p_vec <- fx$prop[1L, , drop = TRUE]
  cst <- fx$constants
  if (!identical(cst, NA) && length(cst)) p_vec <- c(p_vec, cst)
  pt <- EMC2:::ParamTable_create_from_pvector_designs(p_vec, des, nrow(fx$dadm))
  EMC2:::ParamTable_map_designs(pt, des,
                                stats::setNames(rep(TRUE, length(des)), names(des)))
  EMC2:::ParamTable_joint_cells(pt, params)
}

# The planned mapping route against the independent reference, which takes no
# cell short cuts at all.  Lives here rather than in one test file so every file
# that touches the mapper can make the same comparison.
expect_mapper_agrees <- function(fx, designs = NULL, label = "") {
  got <- audit_prologue(fx, designs)$pars
  want <- audit_reference(fx, designs)
  # Compare column by column rather than as one flattened block, so a failure
  # names the parameter that moved instead of only a maximum difference.
  shared <- dimnames(want)[[2L]]
  testthat::expect_true(all(shared %in% dimnames(got)[[2L]]), info = label)
  testthat::expect_equal(dim(got)[c(1L, 3L)], dim(want)[c(1L, 3L)], info = label)
  for (nm in shared) {
    expect_bit_identical(got[, nm, , drop = FALSE], want[, nm, , drop = FALSE],
                         paste(label, "-", nm))
  }
}
