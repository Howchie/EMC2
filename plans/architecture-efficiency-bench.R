# Reproduce the measurements in architecture-efficiency-audit.md.
# Run from the repository root against a matching installed package:
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript plans/architecture-efficiency-bench.R
# No package functions or fit results are changed. Files from Rprof are temporary.
suppressPackageStartupMessages(library(EMC2))
ns <- asNamespace("EMC2")
cat("Namespace:", getNamespaceInfo(ns, "path"), "\n")
cat("DLL:", getLoadedDLLs()[["EMC2"]][["path"]], "\n")
print(tools::md5sum(c("src/EMC2.so", getLoadedDLLs()[["EMC2"]][["path"]])))
cat(R.version.string, "\n")
print(Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")))
reps <- as.integer(Sys.getenv("AUDIT_REPS", "12"))
np <- as.integer(Sys.getenv("AUDIT_PARTICLES", "100"))
nt <- as.integer(Sys.getenv("AUDIT_TRIALS", "2000"))
mode <- Sys.getenv("AUDIT_MODE", "likelihood")

fixture <- function(kind, width = 1L, cells = NULL, trials = nt, particles = np) {
  set.seed(91026)
  d <- forstmann[forstmann$subjects == levels(forstmann$subjects)[1L], ]
  d <- d[rep(seq_len(nrow(d)), length.out = trials), , drop = FALSE]
  d$subjects <- factor("s1")
  d$rt <- d$rt + runif(nrow(d), 0, 0.1)
  d$condition <- factor(rep(seq_len(width), length.out = trials))
  d$x <- rep(seq_len(if (is.null(cells)) 10L else cells), length.out = trials) / 1000
  rownames(d) <- NULL
  model <- get(kind, ns)
  vform <- if (!is.null(cells)) v ~ x else if (width > 1L) v ~ 0 + condition else v ~ lM
  forms <- if (kind == "DDM") {
    list(if (width > 1L) vform else v ~ S, a ~ E, t0 ~ 1)
  } else list(vform, B ~ E, A ~ 1, t0 ~ 1)
  if (kind == "BAwL") forms <- c(forms, list(k ~ 1))
  des <- suppressMessages(design(data = d, model = model,
    formula = forms, matchfun = function(x) x$S == x$lR, report_p_vector = FALSE,
    constants = if (kind == "BAwL") c(sv = 0) else c(s = 0)))
  emc <- suppressMessages(make_emc(d, des, n_chains = 1, compress = TRUE,
    rt_resolution = NULL, type = "single"))
  centre <- sampled_pars(des)
  centre[] <- 0
  centre[grepl("^v($|_condition)", names(centre))] <- if (kind == "RDM") log(2) else 2
  centre["t0"] <- log(0.15)
  if ("A" %in% names(centre)) centre["A"] <- log(0.3)
  if ("k" %in% names(centre)) centre["k"] <- log(0.2)
  p <- sweep(matrix(rnorm(particles * length(centre), 0, 0.04), particles,
    dimnames = list(NULL, names(centre))), 2, centre, "+")
  list(dadm = emc[[1]]$data[[1]], model = des$model, p = p, emc = emc)
}

elapsed <- function(fun, n = reps) {
  fun()
  median(replicate(3L, unname(system.time(for (i in seq_len(n)) fun())[["elapsed"]]) / n))
}
ll_args <- function(f) {
  m <- f$model()
  cc <- attr(f$dadm, "constants")
  list(particle_matrix = f$p, data = f$dadm, constants = if (is.null(cc)) NA else cc,
    designs = EMC2:::.oo_expanded_designs(f$dadm, expand = FALSE),
    bounds = m$bound, transforms = m$transform, pretransforms = m$pre_transform,
    p_types = names(m$p_types), trend = m$trend)
}
report <- function(kind, width = 1L, cells = NULL) {
  f <- fixture(kind, width, cells)
  a <- ll_args(f)
  m <- f$model()
  whole <- function() EMC2:::calc_ll_manager(f$p, f$dadm, f$model, r_cores = 1)
  pro <- function() do.call(EMC2:::pt_prologue_oo, a)
  direct <- function() do.call(EMC2:::calc_ll_oo, c(a, list(type = m$c_name, min_ll = log(1e-10))))
  stopifnot(identical(whole(), direct()), all(is.finite(whole())))
  # Exclude benchmarks that merely exercise an invalid-parameter floor.
  stopifnot(any(abs(whole() - length(attr(f$dadm, "expand")) * log(1e-10)) > 1))
  tw <- elapsed(whole); tp <- elapsed(pro); td <- elapsed(direct)
  data.frame(model = kind, width = width, cells = if (is.null(cells)) NA else cells,
    parameters = ncol(f$p), rows = nrow(f$dadm), particles = nrow(f$p),
    whole_ms = 1000 * tw, prologue_ms = 1000 * tp, direct_ms = 1000 * td,
    prologue_pct = 100 * tp / tw, first_ll = whole()[1])
}

if (mode == "likelihood") {
  for (kind in c("BAwL", "RDM", "DDM")) {
    for (width in c(1L, 64L)) print(report(kind, width), row.names = FALSE)
  }
  for (cells in c(255L, 256L, 257L, 512L)) print(report("RDM", cells = cells), row.names = FALSE)
  cat("Small-batch costs (RDM, 2000 trials, direct vs manager):\n")
  f <- fixture("RDM")
  for (n in c(1L, 2L, 5L, 25L, 100L)) {
    if (n > np) next
    ff <- f; ff$p <- f$p[seq_len(n), , drop = FALSE]
    a <- ll_args(ff); m <- ff$model()
    print(data.frame(particles = n,
      manager_ms = 1000 * elapsed(function() EMC2:::calc_ll_manager(ff$p, ff$dadm, ff$model)),
      direct_ms = 1000 * elapsed(function() do.call(EMC2:::calc_ll_oo,
        c(a, list(type = m$c_name, min_ll = log(1e-10)))))))
  }
} else if (mode == "wide") {
  for (width in c(256L, 257L)) print(report("RDM", width = width), row.names = FALSE)
  # Same particles/data/design mapping on both arms: append one unused design
  # row to cross the hard cutoff without changing any fitted quantity.
  f <- fixture("RDM", width = 256L)
  a <- ll_args(f)
  b <- a
  dm <- a$designs$v
  stopifnot(nrow(dm) == 256L)
  extra <- rbind(dm, dm[1L, , drop = FALSE])
  attr(extra, "expand") <- attr(dm, "expand")
  b$designs$v <- extra
  call <- function(z) do.call(EMC2:::calc_ll_oo, c(z, list(type = "RDM", min_ll = log(1e-10))))
  la <- call(a); lb <- call(b)
  stopifnot(identical(la, lb))
  cat("Identical-result cutoff A/B (unused 257th row):\n")
  print(data.frame(cell_ms = 1000 * elapsed(function() call(a)),
    fallback_ms = 1000 * elapsed(function() call(b)), max_ll_delta = max(abs(la - lb))))
} else if (mode == "particle") {
  # Controlled sample-stage particle step with fixed, well-conditioned proposals.
  # This exercises all four mixture components; it is not a converged fit or ESS estimate.
  for (kind in c("BAwL", "RDM", "DDM")) {
    f <- fixture(kind, trials = as.integer(Sys.getenv("AUDIT_TRIALS", "200")))
    p <- ncol(f$p); centre <- setNames(colMeans(f$p), colnames(f$p)); V <- diag(0.01, p)
    tune <- EMC2:::check_tune_settings(list(search_width = 1, components = rep(1L, p), shared_ll_idx = rep(1L, p)), p, "sample", np)
    # A real sample-stage run inherits three epsilon values from earlier stages.
    pm <- EMC2:::check_sampling_settings(list(list(epsilon = rep(0.7, 3))), "sample", p, np)
    idx <- list(rep(TRUE, p))
    sc <- EMC2:::build_subject_chol_cache(V, V, idx)
    gc <- EMC2:::build_group_chol_cache(V, idx)
    ll <- EMC2:::calc_ll_manager(matrix(centre, 1, dimnames = list(NULL, names(centre))), f$dadm, f$model)
    fun <- function() EMC2:::new_particle(1L, f$dadm, pm, centre, V, centre, V, ll,
      parameters = NULL, model = f$model, stage = "sample", type = "standard", tune = tune,
      chol_cache = sc, group_chol = gc, current_alpha = centre,
      population_mu = centre, population_var = V)
    set.seed(91026)
    fun()
    prof <- tempfile(fileext = ".Rprof")
    Rprof(prof, interval = 0.01)
    tm <- system.time(for (i in seq_len(500L)) fun())
    Rprof(NULL)
    cat("\nParticle step:", kind, "p=", p, "trials=", nrow(f$dadm), "elapsed=", tm[["elapsed"]], "\n")
    print(head(summaryRprof(prof)$by.total, 14L))
    unlink(prof)
  }
} else stop("AUDIT_MODE must be likelihood, wide, or particle")
