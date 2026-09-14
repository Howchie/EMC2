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
# Everything here is a fixture or an assertion.  The timing helper below is
# intentionally shared with WorkingTests/; testthat callers do not make timing
# assertions, so noisy clocks cannot weaken a numerical gate.
# Timing belongs in WorkingTests, which reports the direct-call measurements.

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


# Per-model formulas and the constants each one needs to be identified, so a
# test can name a model rather than restating its whole specification.  The
# PCOUNTER branches are deliberately explicit: a benchmark case says which
# kernel path it is exercising rather than relying on an accidental zero in a
# sampled vector.
.audit_pcounter_branch <- function(branch = NULL) {
  if (is.null(branch) || !length(branch)) return("default")
  if (length(branch) != 1L || is.na(branch)) stop("branch must be one string")
  key <- tolower(trimws(as.character(branch)))
  key <- gsub("[^a-z0-9]+", "_", key)
  key <- gsub("^_|_$", "", key)
  if (!nzchar(key) || key %in% c("default", "baseline")) return("default")
  if (key %in% c("constant", "constants", "all_constant", "const")) return("constant")
  if (key %in% c("varying", "nu_varying", "nu_cells", "cells")) return("varying")
  if (key %in% c("k_varying", "k_cells", "threshold", "threshold_varying")) return("k_varying")
  if (key %in% c("continuous", "continuous_nu_t0", "continuous_covariate")) return("continuous")
  if (key %in% c("degenerate", "all_degenerate", "all_zero", "zero")) return("degenerate")
  if (key %in% c("sv_zero", "sv0", "zero_sv")) return("sv_zero")
  if (key %in% c("gamma_zero", "gamma0", "zero_gamma")) return("gamma_zero")
  if (key %in% c("omega_zero", "omega0", "zero_omega")) return("omega_zero")
  if (key %in% c("k_zero", "k0", "zero_k")) return("k_zero")
  if (key %in% c("k_1023", "k1023", "large_k", "large_threshold")) return("k_1023")
  # A combined spelling is useful at the command line and keeps the matrix
  # labels readable without making callers know an internal enum.
  if (grepl("mixed", key) ||
      (grepl("trunc", key) && grepl("omit", key))) return("mixed")
  if (grepl("trunc", key)) return("truncation")
  if (grepl("omis|omit", key)) return("omission")
  stop("unknown PCOUNTER audit branch: ", branch)
}

.audit_model_spec <- function(model, cells = NULL, width = 1L, branch = NULL) {
  vform <- if (!is.null(cells)) v ~ x else if (width > 1L) v ~ 0 + condition else v ~ lM
  if (identical(model, "PCOUNTER")) {
    b <- .audit_pcounter_branch(branch)
    shape_rhs <- if (!is.null(cells)) ~ x else if (width > 1L) ~ 0 + condition else ~ lM
    constant_rhs <- ~ 1
    f <- function(lhs, rhs) {
      rhs_txt <- paste(deparse(rhs[[length(rhs)]]), collapse = " ")
      stats::as.formula(paste(lhs, "~", rhs_txt))
    }
    forms <- list(
      nu = f("nu", constant_rhs), sv = f("sv", constant_rhs),
      gamma = f("gamma", constant_rhs), k = f("k", constant_rhs),
      omega = f("omega", constant_rhs), t0 = f("t0", constant_rhs))
    # Keep the normal fixture on the fixed-rate member.  Explicit branches
    # below then make every departure from it visible in a case label.
    forms$nu <- f("nu", shape_rhs)
    constants <- c(sv = log(0), gamma = log(0), omega = log(0),
                   pContaminant = qnorm(0), pGuess = qnorm(0))
    centre <- c(nu = log(10), sv = log(2), gamma = log(0.5),
                k = log(2), omega = log(1), t0 = log(0.1))
    if (b == "constant") {
      forms <- stats::setNames(lapply(names(forms), f, rhs = constant_rhs),
                               names(forms))
      constants <- c(pContaminant = qnorm(0), pGuess = qnorm(0))
    } else if (b == "varying") {
      forms$nu <- f("nu", shape_rhs)
    } else if (b == "k_varying") {
      forms$nu <- f("nu", shape_rhs)
      forms$k <- f("k", shape_rhs)
    } else if (b == "continuous") {
      forms$nu <- f("nu", shape_rhs)
      forms$t0 <- f("t0", shape_rhs)
    } else if (b == "degenerate") {
      forms$nu <- f("nu", constant_rhs)
      constants <- c(sv = log(0), gamma = log(0), omega = log(0), k = log(0),
                     pContaminant = qnorm(0), pGuess = qnorm(0))
    } else if (b == "sv_zero") {
      constants <- c(sv = log(0), pContaminant = qnorm(0), pGuess = qnorm(0))
    } else if (b == "gamma_zero") {
      constants <- c(gamma = log(0), pContaminant = qnorm(0), pGuess = qnorm(0))
    } else if (b == "omega_zero") {
      constants <- c(omega = log(0), pContaminant = qnorm(0), pGuess = qnorm(0))
    } else if (b == "k_zero") {
      constants <- c(k = log(0), pContaminant = qnorm(0), pGuess = qnorm(0))
    } else if (b == "k_1023") {
      constants <- c(k = log(1023), pContaminant = qnorm(0), pGuess = qnorm(0))
      centre[["k"]] <- log(1023)
    } else if (b %in% c("mixed", "omission")) {
      # Omission rows need an explicit finite omission mass; otherwise their
      # correct likelihood is zero and the benchmark would only measure the
      # invalid floor.
      constants <- c(pContaminant = qnorm(0.10), pGuess = qnorm(0))
    } else if (b == "truncation") {
      constants <- c(pContaminant = qnorm(0), pGuess = qnorm(0))
    }
    return(list(formula = unname(forms), constants = constants,
                centre = centre, branch = b))
  }
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


# Apply the data-side part of an explicit PCOUNTER branch.  The ordinary
# fixture is untouched; these deterministic edits are only for short benchmark
# cases that need omission, censoring or truncation rows.
.audit_pcounter_branch_data <- function(data, branch) {
  if (!branch %in% c("mixed", "truncation", "omission")) return(data)
  out <- data
  n <- nrow(out)
  if (branch %in% c("mixed", "omission") && n) {
    omit <- if (n >= 10L) seq.int(10L, n, by = 10L) else integer(0L)
    if (length(omit)) {
      out$rt[omit] <- Inf
      out$R[omit] <- NA
    }
  }
  if (branch %in% c("mixed", "truncation")) {
    # Keep missing rows while applying the finite observation window.  The
    # lower edge is below the fixture's fastest response; the upper edge
    # removes only a small, deterministic tail.
    lo <- 0.20
    hi <- 1.10
    keep <- !is.finite(out$rt) | (out$rt >= lo & out$rt <= hi)
    if (any(keep)) out <- out[keep, , drop = FALSE]
    out$LT <- rep(lo, nrow(out))
    out$UT <- rep(hi, nrow(out))
  }
  if (identical(branch, "mixed")) {
    out$UC <- rep(0.75, nrow(out))
    upper <- is.finite(out$rt) & out$rt > out$UC
    out$rt[upper] <- Inf
    out$R[upper] <- NA
  }
  rownames(out) <- NULL
  out
}

# A complete likelihood fixture: compressed data, the model closure, and a
# particle cloud centred somewhere the likelihood is actually finite.  At
# sampled_pars()' zeros t0 is 1 s, which floors every row and leaves the
# fixture unable to discriminate anything at all.
audit_fixture <- function(model = "RDM", n_trials = 400L, n_particles = 5L,
                          cells = NULL, width = 1L, formula = NULL,
                          constants = NULL, data = NULL, seed = 91026L,
                          spread = 0.04, rt_resolution = NULL, ...,
                          branch = NULL) {
  spec <- .audit_model_spec(model, cells, width, branch = branch)
  branch_name <- if (identical(model, "PCOUNTER")) spec$branch else "default"
  dat <- if (is.null(data)) audit_data(n_trials, cells, width, seed = seed) else data
  if (identical(model, "PCOUNTER"))
    dat <- .audit_pcounter_branch_data(dat, branch_name)
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
  list(model_name = model, branch = branch_name, dadm = dadm,
       model = des$model, des = des, data = dat, emc = emc, prop = prop,
       centre = centre, p_types = names(des$model()$p_types),
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

# Time only the direct native route.  The three warm calls make lazy loading
# and one-time native setup part of neither clock; at least seven timed calls
# leave a useful median/MAD even for a very short fixture.  Profiling is an
# explicit explanatory pass and is disabled for the speed numbers.
audit_time_direct <- function(fx, reps = 7L, profile = FALSE,
                              min_ll = log(1e-10)) {
  reps <- suppressWarnings(as.integer(reps)[1L])
  if (!is.finite(reps) || reps < 1L) stop("reps must be a positive integer")
  reps <- max(7L, reps)
  if (length(profile) != 1L || is.na(profile))
    stop("profile must be TRUE or FALSE")
  profile <- isTRUE(profile)

  has_stats <- exists("emc_kernel_stats", envir = asNamespace("EMC2"),
                      inherits = FALSE)
  previous_stats <- FALSE
  if (has_stats) {
    previous_stats <- isTRUE(tryCatch(EMC2:::emc_kernel_stats(),
                                      error = function(e) FALSE))
    # Never let a caller's profiling state contaminate a speed measurement.
    if (previous_stats)
      try(EMC2:::emc_kernel_stats(FALSE), silent = TRUE)
  }
  on.exit(if (has_stats) try(EMC2:::emc_kernel_stats(previous_stats),
                             silent = TRUE), add = TRUE)

  for (i in seq_len(3L)) invisible(audit_ll_direct(fx, min_ll = min_ll))

  if (profile && has_stats) {
    try(EMC2:::emc_kernel_stats_reset(), silent = TRUE)
    try(EMC2:::emc_kernel_stats(TRUE), silent = TRUE)
  }
  elapsed <- numeric(reps)
  cpu <- numeric(reps)
  lls <- vector("list", reps)
  for (i in seq_len(reps)) {
    tm <- system.time(lls[[i]] <- audit_ll_direct(fx, min_ll = min_ll))
    elapsed[[i]] <- unname(tm[["elapsed"]])
    cpu[[i]] <- unname(tm[["user.self"]] + tm[["sys.self"]])
  }

  st <- NULL
  if (profile && has_stats) {
    try(EMC2:::emc_kernel_stats(FALSE), silent = TRUE)
    st <- tryCatch(EMC2:::emc_kernel_stats_read(),
                   error = function(e) NULL)
  }
  kernel <- c(rows = NA_real_, cells = NA_real_, seconds = NA_real_)
  kernel_calls <- NA_real_
  pcounter <- NULL
  if (is.data.frame(st) && nrow(st)) {
    kernel <- c(rows = sum(st$rows), cells = sum(st$cells),
                seconds = sum(st$seconds))
    if ("calls" %in% names(st)) kernel_calls <- sum(st$calls)
    if (is.data.frame(st) && any(st$model == "PCOUNTER")) {
      pcounter <- st[which(st$model == "PCOUNTER")[1L], , drop = FALSE]
    }
  }

  checksums <- vapply(lls, function(x) sum(as.numeric(x)), numeric(1L))
  last_ll <- lls[[reps]]
  n_likelihood <- length(attr(fx$dadm, "expand"))
  if (!n_likelihood) n_likelihood <- nrow(fx$dadm)
  floor_ll <- n_likelihood * min_ll
  at_floor <- vapply(lls, function(x) {
    any(is.finite(x) & abs(as.numeric(x) - floor_ll) <=
          max(1e-12, abs(floor_ll) * 1e-12))
  }, logical(1L))
  med <- stats::median(elapsed)
  list(
    warm_calls = 3L, measured_calls = reps,
    elapsed = elapsed, cpu = cpu,
    median_seconds = unname(med),
    mad_seconds = unname(stats::mad(elapsed, center = med, constant = 1)),
    cpu_seconds = unname(stats::median(cpu)),
    cpu_total_seconds = unname(sum(cpu)),
    likelihood_checksum = unname(checksums[[reps]]),
    checksum_reproducible = all(vapply(checksums, identical, logical(1L),
                                       checksums[[reps]])),
    native_calls = as.integer(reps),
    native_call_count = as.integer(reps),
    kernel_calls = kernel_calls,
    kernel_rows = unname(kernel[["rows"]]),
    kernel_cells = unname(kernel[["cells"]]),
    kernel_seconds = unname(kernel[["seconds"]]),
    pcounter = pcounter,
    invalid_floor = any(at_floor),
    invalid_floor_calls = sum(at_floor),
    floor_log_likelihood = unname(floor_ll),
    all_finite = all(vapply(lls, function(x) all(is.finite(x)), logical(1L))),
    last_likelihood = last_ll
  )
}

# Descriptive alias used by a few one-off audit scripts.
audit_direct_timing <- audit_time_direct

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
