#!/usr/bin/env Rscript

# Baseline oracle for the EMC2 refactor.  Capture and compare are deliberately
# separate command-line operations: never load two EMC2 shared libraries in
# one R process.

oracle_usage <- function() {
  cat(paste0(
    "Usage:\n",
    "  Rscript tools/refactor_oracle.R capture --output FILE [--library DIR] \\\n",
    "      [--source-revision REV] [--build-profile PROFILE]\n",
    "  Rscript tools/refactor_oracle.R compare --old OLD.rds --new NEW.rds \\\n",
    "      [--mode exact|tolerance] [--abs-tol X] [--rel-tol X]\n",
    "\n",
    "Modes:\n",
    "  capture    Load EMC2 once from one library, run fixed direct-kernel and\n",
    "             likelihood fixtures, and write a full-precision RDS.\n",
    "  compare    Compare two independently captured RDS files.  exact is the\n",
    "             default; tolerance applies abs-tol + rel-tol * max(|old|,|new|).\n",
    "\n",
    "Options:\n",
    "  --library DIR          Package library to prepend to .libPaths/R_LIBS.\n",
    "  --output FILE          Capture output RDS path.\n",
    "  --old FILE, --new FILE Compare input RDS paths.\n",
    "  --mode MODE            exact or tolerance (default exact).\n",
    "  --abs-tol X            Absolute tolerance (default 0).\n",
    "  --rel-tol X            Relative tolerance (default 0).\n",
    "  --source-revision REV  Override source revision metadata.\n",
    "  --build-profile NAME   Override compiler/build-profile metadata.\n",
    "  --help                 Show this help.\n",
    "\n",
    "Examples:\n",
    "  Rscript tools/refactor_oracle.R capture --library /tmp/emc2-refactor-lib \\\n",
    "      --output baseline.rds\n",
    "  Rscript tools/refactor_oracle.R compare --old baseline.rds --new candidate.rds \\\n",
    "      --mode tolerance --abs-tol 1e-12 --rel-tol 1e-10\n",
    sep = ""
  ))
}

oracle_fail <- function(message, status = 2L) {
  cat(sprintf("ERROR: %s\n", message), file = stderr())
  quit(save = "no", status = status, runLast = FALSE)
}

parse_oracle_args <- function(args) {
  if (!length(args) || any(args %in% c("--help", "-h"))) {
    oracle_usage()
    quit(save = "no", status = 0L, runLast = FALSE)
  }
  mode <- args[[1L]]
  if (!mode %in% c("capture", "compare")) {
    oracle_fail(sprintf("unknown mode '%s'; expected capture or compare (use --help)", mode))
  }
  rest <- args[-1L]
  opts <- list()
  i <- 1L
  while (i <= length(rest)) {
    token <- rest[[i]]
    if (!startsWith(token, "--")) {
      oracle_fail(sprintf("unexpected positional argument '%s'", token))
    }
    eq <- regexpr("=", token, fixed = TRUE)[[1L]]
    if (eq > 0L) {
      key <- substring(token, 3L, eq - 1L)
      value <- substring(token, eq + 1L)
    } else {
      key <- substring(token, 3L)
      if (key %in% c("help")) {
        oracle_usage()
        quit(save = "no", status = 0L, runLast = FALSE)
      }
      if (i == length(rest) || startsWith(rest[[i + 1L]], "--")) {
        oracle_fail(sprintf("option '--%s' requires a value", key))
      }
      i <- i + 1L
      value <- rest[[i]]
    }
    if (!key %in% c("library", "output", "old", "new", "mode",
                    "abs-tol", "rel-tol", "source-revision", "build-profile")) {
      oracle_fail(sprintf("unknown option '--%s'", key))
    }
    opts[[key]] <- value
    i <- i + 1L
  }
  opts$command <- mode
  opts
}

as_nonnegative_number <- function(value, option) {
  out <- suppressWarnings(as.numeric(value))
  if (length(out) != 1L || is.na(out) || !is.finite(out) || out < 0) {
    oracle_fail(sprintf("--%s must be one finite non-negative number, got '%s'", option, value))
  }
  out
}

safe_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) NA_character_ else value
}

oracle_git_revision <- function() {
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)
  script_root <- character()
  if (length(file_arg)) {
    script_file <- sub("^--file=", "", file_arg[[1L]])
    script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
    script_root <- dirname(script_dir)
  }
  candidates <- unique(c(getwd(), script_root))
  for (root in candidates[nzchar(candidates)]) {
    value <- tryCatch(
      suppressWarnings(system2("git", c("-C", root, "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)),
      error = function(e) character()
    )
    value <- trimws(paste(value, collapse = ""))
    if (grepl("^[0-9a-fA-F]{7,64}$", value)) return(value)
  }
  NULL
}

package_source_revision <- function(description, override = NULL) {
  if (!is.null(override) && length(override) && !is.na(override) && nzchar(override)) {
    return(list(value = override, source = "--source-revision"))
  }
  for (name in c("EMC2_SOURCE_REVISION", "GIT_COMMIT", "SOURCE_REVISION")) {
    value <- safe_env(name)
    if (length(value) && !is.na(value) && nzchar(value)) {
      return(list(value = value, source = paste0("environment:", name)))
    }
  }
  for (name in c("RemoteSha", "GithubSHA1", "GitSHA1")) {
    value <- description[[name]]
    if (!is.null(value) && length(value) && !is.na(value) && nzchar(value)) {
      return(list(value = unname(value), source = paste0("DESCRIPTION:", name)))
    }
  }
  revision <- oracle_git_revision()
  if (!is.null(revision)) return(list(value = revision, source = "git:repository"))
  list(value = NA_character_, source = "unavailable")
}

capture_metadata <- function(description, library_dir, opts) {
  profile <- opts[["build-profile"]]
  if (is.null(profile) || !nzchar(profile)) {
    for (name in c("EMC2_BUILD_PROFILE", "EMC2_PROFILE", "R_BUILD_PROFILE")) {
      profile <- safe_env(name)
      if (!is.na(profile)) break
    }
  }
  if (is.null(profile) || !length(profile) || is.na(profile) || !nzchar(profile)) {
    profile <- "unspecified"
  }
  compiler <- list(
    CC = safe_env("CC"),
    CXX = safe_env("CXX"),
    CXXFLAGS = safe_env("CXXFLAGS"),
    CPPFLAGS = safe_env("CPPFLAGS"),
    PKG_CPPFLAGS = safe_env("PKG_CPPFLAGS"),
    PKG_LIBS = safe_env("PKG_LIBS"),
    MAKEFLAGS = safe_env("MAKEFLAGS"),
    R_MAKEVARS_USER = safe_env("R_MAKEVARS_USER"),
    package_built = if (!is.null(description$Built)) unname(description$Built) else NA_character_
  )
  list(
    package = list(
      name = unname(description$Package),
      version = unname(description$Version),
      path = normalizePath(find.package("EMC2"), winslash = "/", mustWork = FALSE),
      library_path = normalizePath(library_dir, winslash = "/", mustWork = FALSE)
    ),
    source_revision = package_source_revision(description, opts[["source-revision"]]),
    build_profile = profile,
    compiler = compiler,
    R = list(
      version = R.version$version.string,
      platform = R.version$platform,
      arch = R.version$arch,
      major = R.version$major,
      minor = R.version$minor,
      R_LIBS = safe_env("R_LIBS"),
      R_LIBS_USER = safe_env("R_LIBS_USER")
    ),
    platform = as.list(Sys.info()),
    fixture_version = "refactor-oracle-fixture-1"
  )
}

snapshot_policy <- function() {
  list(
    policy_name = "baseline-rds-not-text-snapshot",
    regenerate_snapshots = FALSE,
    snapshots_regenerated = FALSE,
    precision = "full-precision-RDS",
    comparison = "exact-or-explicit-absolute-relative-tolerance",
    missingness_column = FALSE,
    note = "Capture records an oracle artifact only; repository snapshots are not regenerated."
  )
}

capture_output <- function(value) {
  list(
    value = value,
    dimensions = dim(value),
    dim = dim(value),
    dimnames = dimnames(value),
    names = names(value),
    attributes = attributes(value)
  )
}

make_case <- function(label, model, variant, operation, input, value, provenance) {
  list(
    label = label,
    model = model,
    variant = variant,
    operation = operation,
    input = input,
    output = capture_output(value),
    provenance = provenance
  )
}


run_capture <- function(opts) {
  output_path <- opts$output
  if (is.null(output_path) || !nzchar(output_path)) {
    oracle_fail("capture requires --output FILE")
  }
  library_dir <- opts$library
  if (is.null(library_dir) || !nzchar(library_dir)) {
    library_dir <- NA_character_
  } else {
    if (!dir.exists(library_dir)) oracle_fail(sprintf("package library does not exist: %s", library_dir))
    library_dir <- normalizePath(library_dir, winslash = "/", mustWork = TRUE)
    existing_libs <- Sys.getenv("R_LIBS", unset = "")
    lib_values <- c(library_dir, if (nzchar(existing_libs)) existing_libs else character())
    Sys.setenv(R_LIBS = paste(unique(lib_values), collapse = .Platform$path.sep))
    .libPaths(unique(c(library_dir, .libPaths())))
  }
  if (!requireNamespace("EMC2", quietly = TRUE)) {
    oracle_fail("EMC2 is not available; provide an installed package library with --library DIR")
  }
  suppressPackageStartupMessages(library("EMC2", character.only = TRUE))
  description <- packageDescription("EMC2")
  if (is.na(library_dir)) library_dir <- dirname(find.package("EMC2"))

  # All direct kernels use this fixed vector.  The package's kernel APIs accept
  # non-finite and missing RT values and return their explicit boundary values.
  rt <- c(0.25, 0.75, -Inf, Inf, NA_real_)
  rt_names <- c("finite_fast", "finite_slow", "negative_infinity", "positive_infinity", "missing")
  names(rt) <- rt_names

  lba <- EMC2::LBA(posdrift = TRUE)
  lba_pars <- cbind(
    v = c(1.10, 1.10, 1.10, 1.10, 1.10),
    sv = c(.80, .80, .80, .80, .80),
    b = c(1.10, 1.10, 1.10, 1.10, 1.10),
    A = c(.30, .30, .30, .30, .30),
    t0 = c(.10, .10, .10, .10, .10)
  )
  rownames(lba_pars) <- rt_names
  lba_density <- lba$dfun(rt, lba_pars)
  lba_cdf <- lba$pfun(rt, lba_pars)

  rdm <- EMC2::RDM()
  rdm_pars <- cbind(v = rep(.80, length(rt)), B = rep(.90, length(rt)),
                    A = rep(.20, length(rt)), t0 = rep(.10, length(rt)),
                    s = rep(.70, length(rt)))
  rownames(rdm_pars) <- rt_names
  rdm_density <- rdm$dfun(rt, rdm_pars)
  rdm_cdf <- rdm$pfun(rt, rdm_pars)

  rexg <- EMC2::REXG()
  rexg_pars <- cbind(mu = rep(.40, length(rt)), sigma = rep(.08, length(rt)),
                     tau = rep(.12, length(rt)), t0 = rep(.10, length(rt)))
  rownames(rexg_pars) <- rt_names
  rexg_density <- rexg$dfun(rt, rexg_pars)
  rexg_cdf <- rexg$pfun(rt, rexg_pars)

  pcounter <- EMC2::PCOUNTER()
  pcounter_pars <- cbind(nu = rep(3, length(rt)), sv = rep(.40, length(rt)),
                         gamma = rep(.20, length(rt)), k = rep(3, length(rt)),
                         omega = rep(.50, length(rt)), t0 = rep(.10, length(rt)))
  rownames(pcounter_pars) <- rt_names
  pcounter_density <- pcounter$dfun(rt, pcounter_pars)
  pcounter_cdf <- pcounter$pfun(rt, pcounter_pars)

  # Stop helper uses natural-scale parameters and explicitly exercises a finite
  # response window and an uncensored (+Inf) upper bound.
  ssexg <- EMC2::SSEXG(stop_method = "gl", stop_n_nodes = 64L)
  ssexg_pars <- matrix(rep(c(.40, .05, .10, .30, .025, .05, 0, 0, .05, .05), 2L),
                       nrow = 2L, byrow = TRUE,
                       dimnames = list(c("finite_window", "uncensored"), names(ssexg$p_types)))
  ssexg_pars <- cbind(ssexg_pars, SSD = c(.20, .30))
  ssexg_upper <- c(.70, Inf)
  ssexg_stop_probability <- ssexg$sfun(ssexg_pars, n_acc = 1L, upper = ssexg_upper)

  # The likelihood fixture intentionally contains an NA response and +Inf RT;
  # LT and UC are censor/truncation bounds, not a synthetic missingness column.
  raw_data <- data.frame(
    subjects = factor(c("s1", "s1", "s1")),
    R = factor(c("left", "right", NA_character_), levels = c("left", "right")),
    rt = c(.45, .62, Inf),
    LT = c(.10, .10, .10),
    UC = c(Inf, Inf, Inf),
    stringsAsFactors = TRUE
  )
  likelihood_design <- NULL
  invisible(capture.output(
    likelihood_design <- EMC2::design(
      data = raw_data,
      model = EMC2::LBA,
      formula = list(v ~ 1, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1)
    )
  ))
  likelihood_emc <- NULL
  invisible(capture.output(
    likelihood_emc <- suppressWarnings(suppressMessages(EMC2::make_emc(
      raw_data, likelihood_design, type = "single", n_chains = 1L,
      compress = FALSE, verbose = FALSE
    )))
  ))
  dadm <- likelihood_emc[[1L]]$data[[1L]]
  model <- likelihood_emc[[1L]]$model()
  p_types <- attr(dadm, "sampled_p_names")
  designs <- list()
  for (p_name in p_types) {
    design_matrix <- attr(dadm, "designs")[[p_name]]
    designs[[p_name]] <- design_matrix[attr(design_matrix, "expand"), , drop = FALSE]
  }
  particle_matrix <- matrix(
    c(.50, log(.80), log(1.10), log(.30), log(.10)),
    nrow = 1L, dimnames = list("particle-1", p_types)
  )
  likelihood_value <- EMC2:::calc_ll_oo(
    particle_matrix = particle_matrix,
    data = dadm,
    constants = attr(dadm, "constants"),
    designs = designs,
    type = model$c_name,
    bounds = model$bound,
    transforms = model$transform,
    pretransforms = model$pre_transform,
    p_types = p_types,
    min_ll = log(1e-10),
    trend = model$trend
  )

  cases <- list(
    lba_density = make_case(
      "lba_density", "LBA", "posdrift=TRUE", "density", list(rt = rt, pars = lba_pars),
      lba_density, list(api = "LBA()$dfun", compiled_kernel = "dlba")
    ),
    lba_cdf = make_case(
      "lba_cdf", "LBA", "posdrift=TRUE", "cdf", list(rt = rt, pars = lba_pars),
      lba_cdf, list(api = "LBA()$pfun", compiled_kernel = "plba")
    ),
    rdm_wald_density = make_case(
      "rdm_wald_density", "RDM", "Wald/RDM", "density", list(rt = rt, pars = rdm_pars),
      rdm_density, list(api = "RDM()$dfun", compiled_kernel = "dWald")
    ),
    rdm_wald_cdf = make_case(
      "rdm_wald_cdf", "RDM", "Wald/RDM", "cdf", list(rt = rt, pars = rdm_pars),
      rdm_cdf, list(api = "RDM()$pfun", compiled_kernel = "pWald")
    ),
    rexg_density = make_case(
      "rexg_density", "REXG", "ex-Gaussian race", "density", list(rt = rt, pars = rexg_pars),
      rexg_density, list(api = "REXG()$dfun", compiled_kernel = "dtexGaussian")
    ),
    rexg_cdf = make_case(
      "rexg_cdf", "REXG", "ex-Gaussian race", "cdf", list(rt = rt, pars = rexg_pars),
      rexg_cdf, list(api = "REXG()$pfun", compiled_kernel = "ptexGaussian")
    ),
    ssexg_stop_probability = make_case(
      "ssexg_stop_probability", "SSEXG", "ex-Gaussian stop-signal helper; GL-64", "stop_probability",
      list(pars = ssexg_pars, n_acc = 1L, upper = ssexg_upper), ssexg_stop_probability,
      list(api = "SSEXG()$sfun", helper = "pstopTEXG", method = "gl", nodes = 64L)
    ),
    pcounter_density = make_case(
      "pcounter_density", "PCOUNTER", "gamma-mixed self-exciting counter", "density",
      list(rt = rt, pars = pcounter_pars), pcounter_density,
      list(api = "PCOUNTER()$dfun", compiled_kernel = "dpcounter")
    ),
    pcounter_cdf = make_case(
      "pcounter_cdf", "PCOUNTER", "gamma-mixed self-exciting counter", "cdf",
      list(rt = rt, pars = pcounter_pars), pcounter_cdf,
      list(api = "PCOUNTER()$pfun", compiled_kernel = "ppcounter")
    ),
    calc_ll_oo_lba = make_case(
      "calc_ll_oo_lba", "LBA", "likelihood dispatcher with LT/UC and +Inf/NA response", "log_likelihood",
      list(raw_data = raw_data, particle_matrix = particle_matrix,
           constants = attr(dadm, "constants"), designs = designs, type = model$c_name,
           p_types = p_types, min_ll = log(1e-10)), likelihood_value,
      list(api = "EMC2:::calc_ll_oo", compiled_kernel = model$c_name)
    )
  )

  artifact <- list(
    schema = "EMC2-refactor-oracle",
    schema_version = 1L,
    metadata = capture_metadata(description, library_dir, opts),
    snapshot_policy = snapshot_policy(),
    cases = cases
  )
  parent <- dirname(normalizePath(output_path, winslash = "/", mustWork = FALSE))
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    oracle_fail(sprintf("cannot create output directory: %s", parent))
  }
  saveRDS(artifact, output_path, version = 3L, compress = FALSE)
  cat(sprintf("Captured %d cases to %s\n", length(cases), normalizePath(output_path, winslash = "/", mustWork = FALSE)))
  invisible(artifact)
}

format_scalar <- function(value) {
  if (length(value) != 1L) return(sprintf("<length %d>", length(value)))
  if (is.nan(value)) return("NaN")
  if (is.na(value)) return("NA")
  if (is.infinite(value)) return(if (value > 0) "+Inf" else "-Inf")
  if (is.numeric(value)) return(sprintf("%.17g", value))
  if (is.character(value)) return(sprintf("%s", encodeString(value, quote = '"')))
  as.character(value)
}

path_index <- function(path, index) {
  if (!length(path)) return(sprintf("[%d]", index))
  paste0(path, "[", index, "]")
}

append_diff <- function(diffs, case_label, path, old, new, message) {
  diffs[[length(diffs) + 1L]] <- list(
    case = case_label,
    path = path,
    old = old,
    new = new,
    message = message
  )
  diffs
}

compare_nodes <- function(old, new, case_label, path, mode, abs_tol, rel_tol, diffs) {
  old_is_atomic <- is.atomic(old) && !is.null(old)
  new_is_atomic <- is.atomic(new) && !is.null(new)
  if (old_is_atomic || new_is_atomic) {
    if (!(old_is_atomic && new_is_atomic)) {
      return(append_diff(diffs, case_label, path, old, new, "type differs"))
    }
    if (!identical(typeof(old), typeof(new))) {
      diffs <- append_diff(diffs, case_label, path, old, new,
                           sprintf("type differs (%s vs %s)", typeof(old), typeof(new)))
    }
    if (length(old) != length(new)) {
      return(append_diff(diffs, case_label, path, old, new,
                         sprintf("length differs (%d vs %d)", length(old), length(new))))
    }
    if (!identical(dim(old), dim(new))) {
      diffs <- append_diff(diffs, case_label, paste0(path, ".dim"), dim(old), dim(new), "dimensions differ")
    }
    if (!identical(dimnames(old), dimnames(new))) {
      diffs <- append_diff(diffs, case_label, paste0(path, ".dimnames"), dimnames(old), dimnames(new), "dimnames differ")
    }
    if (length(old)) {
      for (i in seq_along(old)) {
        x <- old[[i]]
        y <- new[[i]]
        equal <- if (is.numeric(old) && is.numeric(new)) {
          if (is.nan(x) || is.nan(y)) is.nan(x) && is.nan(y) else if (is.na(x) || is.na(y)) is.na(x) && is.na(y) else if (is.infinite(x) || is.infinite(y)) identical(x, y) else if (mode == "exact") identical(x, y) else abs(x - y) <= abs_tol + rel_tol * max(abs(x), abs(y))
        } else {
          identical(x, y)
        }
        if (!equal) {
          if (is.numeric(old) && is.numeric(new) && is.finite(x) && is.finite(y)) {
            delta <- abs(x - y)
            allowed <- if (mode == "exact") 0 else abs_tol + rel_tol * max(abs(x), abs(y))
            msg <- sprintf("value differs: old=%s new=%s abs_diff=%.17g allowed=%.17g",
                           format_scalar(x), format_scalar(y), delta, allowed)
          } else {
            msg <- sprintf("value differs: old=%s new=%s", format_scalar(x), format_scalar(y))
          }
          diffs <- append_diff(diffs, case_label, path_index(path, i), x, y, msg)
        }
      }
    }
    old_attributes <- attributes(old)
    new_attributes <- attributes(new)
    if (!identical(old_attributes, new_attributes)) {
      diffs <- compare_nodes(old_attributes, new_attributes, case_label,
                             paste0(path, ".attributes"), mode, abs_tol, rel_tol, diffs)
    }
    return(diffs)
  }
  if (is.null(old) || is.null(new)) {
    if (!identical(old, new)) diffs <- append_diff(diffs, case_label, path, old, new, "one value is NULL")
    return(diffs)
  }
  if (is.environment(old) || is.environment(new) || is.function(old) || is.function(new)) {
    if (!identical(old, new)) diffs <- append_diff(diffs, case_label, path, old, new, "unsupported reference differs")
    return(diffs)
  }
  if (is.list(old) && is.list(new)) {
    if (length(old) != length(new)) {
      diffs <- append_diff(diffs, case_label, path, old, new,
                           sprintf("list length differs (%d vs %d)", length(old), length(new)))
    }
    common <- seq_len(min(length(old), length(new)))
    old_names <- names(old); new_names <- names(new)
    if (!identical(old_names, new_names)) {
      diffs <- append_diff(diffs, case_label, paste0(path, ".names"), old_names, new_names, "names differ")
    }
    for (i in common) {
      child <- if (!is.null(old_names) && nzchar(old_names[[i]])) paste0(path, "$", old_names[[i]]) else path_index(path, i)
      diffs <- compare_nodes(old[[i]], new[[i]], case_label, child, mode, abs_tol, rel_tol, diffs)
    }
    if (!identical(attributes(old), attributes(new))) {
      diffs <- compare_nodes(attributes(old), attributes(new), case_label,
                             paste0(path, ".attributes"), mode, abs_tol, rel_tol, diffs)
    }
    return(diffs)
  }
  if (!identical(old, new)) append_diff(diffs, case_label, path, old, new, "values differ") else diffs
}

metadata_for_compare <- function(metadata) {
  # Paths and library search paths are expected to differ between isolated
  # old/new subprocesses; all compiler/profile/platform/source fields remain
  # part of the compatibility contract.
  out <- metadata
  if (!is.null(out$package)) {
    out$package$path <- NULL
    out$package$library_path <- NULL
  }
  if (!is.null(out$R)) {
    out$R$R_LIBS <- NULL
    out$R$R_LIBS_USER <- NULL
  }
  out
}

validate_artifact <- function(x, path) {
  if (!is.list(x) || !identical(x$schema, "EMC2-refactor-oracle") ||
      !identical(x$schema_version, 1L) || !is.list(x$metadata) ||
      !is.list(x$snapshot_policy) || !is.list(x$cases) || !length(x$cases)) {
    oracle_fail(sprintf("%s is not a valid EMC2-refactor-oracle RDS (schema version 1 expected)", path))
  }
  invisible(TRUE)
}

run_compare <- function(opts) {
  old_path <- opts$old
  new_path <- opts$new
  if (is.null(old_path) || !nzchar(old_path) || is.null(new_path) || !nzchar(new_path)) {
    oracle_fail("compare requires both --old OLD.rds and --new NEW.rds")
  }
  if (!file.exists(old_path)) oracle_fail(sprintf("old RDS does not exist: %s", old_path))
  if (!file.exists(new_path)) oracle_fail(sprintf("new RDS does not exist: %s", new_path))
  mode <- opts$mode
  if (is.null(mode)) mode <- "exact"
  if (!mode %in% c("exact", "tolerance")) oracle_fail("--mode must be exact or tolerance")
  abs_tol <- if (is.null(opts[["abs-tol"]])) 0 else as_nonnegative_number(opts[["abs-tol"]], "abs-tol")
  rel_tol <- if (is.null(opts[["rel-tol"]])) 0 else as_nonnegative_number(opts[["rel-tol"]], "rel-tol")
  if (mode == "exact" && (abs_tol != 0 || rel_tol != 0)) {
    oracle_fail("exact mode cannot use non-zero --abs-tol/--rel-tol; choose --mode tolerance")
  }
  old <- tryCatch(readRDS(old_path), error = function(e) oracle_fail(sprintf("cannot read old RDS: %s", conditionMessage(e))))
  new <- tryCatch(readRDS(new_path), error = function(e) oracle_fail(sprintf("cannot read new RDS: %s", conditionMessage(e))))
  validate_artifact(old, old_path)
  validate_artifact(new, new_path)

  diffs <- list()
  diffs <- compare_nodes(metadata_for_compare(old$metadata), metadata_for_compare(new$metadata), "<metadata>", "metadata", "exact", 0, 0, diffs)
  diffs <- compare_nodes(old$snapshot_policy, new$snapshot_policy, "<snapshot-policy>", "snapshot_policy", "exact", 0, 0, diffs)
  old_labels <- names(old$cases)
  new_labels <- names(new$cases)
  if (!identical(old_labels, new_labels)) {
    diffs <- append_diff(diffs, "<inventory>", "cases.names", old_labels, new_labels,
                         "case inventories differ; capture fixtures must be identical")
  }
  common_labels <- intersect(old_labels, new_labels)
  for (label in common_labels) {
    old_case <- old$cases[[label]]
    new_case <- new$cases[[label]]
    # Fixture identity, labels, provenance, and metadata are always exact.
    for (field in c("label", "model", "variant", "operation", "input", "provenance")) {
      diffs <- compare_nodes(old_case[[field]], new_case[[field]], label,
                             paste0("cases$", label, "$", field), "exact", 0, 0, diffs)
    }
    # Only recorded numeric outputs are subject to requested tolerance.
    for (field in c("dimensions", "dim", "dimnames", "names", "attributes")) {
      diffs <- compare_nodes(old_case$output[[field]], new_case$output[[field]], label,
                             paste0("cases$", label, "$output$", field), "exact", 0, 0, diffs)
    }
    diffs <- compare_nodes(old_case$output$value, new_case$output$value, label,
                           paste0("cases$", label, "$output$value"), mode, abs_tol, rel_tol, diffs)
  }
  if (length(diffs)) {
    cat(sprintf("Comparison FAILED: %d difference(s) [mode=%s abs_tol=%.17g rel_tol=%.17g]\n",
                length(diffs), mode, abs_tol, rel_tol))
    for (d in diffs) cat(sprintf("[DIFF] case=%s path=%s %s\n", d$case, d$path, d$message))
    quit(save = "no", status = 1L, runLast = FALSE)
  }
  cat(sprintf("Comparison OK: %d cases, mode=%s abs_tol=%.17g rel_tol=%.17g\n",
              length(old$cases), mode, abs_tol, rel_tol))
  invisible(TRUE)
}

main <- function() {
  opts <- parse_oracle_args(commandArgs(trailingOnly = TRUE))
  if (opts$command == "capture") run_capture(opts) else run_compare(opts)
}

tryCatch(main(), error = function(e) {
  oracle_fail(conditionMessage(e), status = 2L)
})
