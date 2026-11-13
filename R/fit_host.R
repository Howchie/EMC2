# A fit runs in its own Rscript process, the fit host.
#
# Chains, rebuild jobs and fork-backend workers are forks of the process that
# runs the fit.  A fork shares its parent's memory only until R's garbage
# collector writes to it, and a full collection writes to every object, so
# each fork ends up holding a private copy of whatever its parent held.  Run
# from an interactive session that is the user's whole workspace, once per
# chain and once per rebuild job, every block.  The host holds the package,
# the fit's arguments and nothing else, so that is all its forks can copy.
#
# The host must end with the session that asked for it, however that session
# ends.  It leads its own process group and watches its stdin, the read end of
# a pipe whose write end only the session holds (emc_host_detach()): when the
# session dies the kernel closes the pipe and the host kills its group.  An
# interrupt in the session (Ctrl-C, RStudio's Stop button) ends the wait below,
# whose exit handler kills the group directly.  The group is killed after a
# successful fit too, so nothing the fit started can outlive it.

.emc_host_wanted <- function(n_cores) {
  if (!isTRUE(getOption("emc2.fit_host", TRUE))) return(FALSE)
  if (nzchar(Sys.getenv("EMC2_IN_HOST"))) return(FALSE)
  if (.Platform$OS.type != "unix") return(FALSE)
  if (!is.numeric(n_cores) || !isTRUE(n_cores > 1)) return(FALSE)
  ns_path <- tryCatch(getNamespaceInfo(asNamespace("EMC2"), "path"),
                      error = function(e) "")
  nzchar(ns_path) && file.exists(file.path(ns_path, "Meta", "package.rds")) &&
    file.exists(file.path(R.home("bin"), "Rscript"))
}

# Arguments the caller supplied, by name, evaluated in `env`; defaults are left
# for the host to evaluate, since some depend on arguments the fit replaces
# (`cores_for_chains = length(emc)` when `emc` is a file name).
.emc_host_args <- function(fun, env, dots) {
  nm <- setdiff(names(formals(fun)), "...")
  given <- nm[!vapply(nm, function(a) eval(call("missing", as.name(a)), env),
                      logical(1))]
  c(mget(given, envir = env), dots)
}

# Values a closure in `x` looks up in the global environment, which the host's
# global environment does not have, and the attached packages they resolve to.
.emc_host_globals <- function(x) {
  globals <- list()
  packages <- character(0)
  seen <- new.env(parent = emptyenv())
  stop_env <- function(e) {
    identical(e, globalenv()) || identical(e, baseenv()) ||
      identical(e, emptyenv()) || isNamespace(e) ||
      !is.null(attr(e, "name"))
  }
  visit_env <- function(e) {
    if (stop_env(e)) return()
    key <- format(e)
    if (exists(key, envir = seen, inherits = FALSE)) return()
    assign(key, TRUE, envir = seen)
    for (nm in ls(e, all.names = TRUE)) {
      if (bindingIsActive(nm, e)) next
      # A promise that cannot be forced (e.g. a lazy-load from another
      # machine's library) refers to no global; it is saved unforced.
      value <- tryCatch(get(nm, envir = e, inherits = FALSE),
                        error = function(err) NULL)
      visit(value)
    }
    visit_env(parent.env(e))
  }
  resolve <- function(name, e) {
    while (!identical(e, emptyenv())) {
      if (exists(name, envir = e, inherits = FALSE)) {
        if (identical(e, globalenv())) {
          if (is.null(globals[[name]])) {
            value <- get(name, envir = e)
            globals[name] <<- list(value)
            visit(value)
          }
        } else if (!is.null(attr(e, "name")) &&
                   startsWith(attr(e, "name"), "package:")) {
          packages <<- union(packages, sub("^package:", "", attr(e, "name")))
        }
        return()
      }
      e <- parent.env(e)
    }
  }
  visit <- function(o) {
    if (typeof(o) == "closure") {
      e <- environment(o)
      if (isNamespace(e)) return()
      for (g in codetools::findGlobals(o)) resolve(g, e)
      visit_env(e)
    } else if (is.environment(o)) {
      visit_env(o)
    } else if (is.list(o) || is.pairlist(o)) {
      for (el in o) visit(el)
    }
    a <- attributes(o)
    a$names <- NULL; a$dim <- NULL; a$dimnames <- NULL; a$row.names <- NULL
    for (el in a) visit(el)
  }
  visit(x)
  list(globals = globals, packages = packages)
}

.EMC_HOST_OPTIONS <- c("width", "digits", "scipen", "warn", "OutDec",
                       "nwarnings", "mc.cores", "rcpp.cache.dir")

.emc_host_payload <- function(fun, args) {
  found <- .emc_host_globals(args)
  opts <- options()
  opts <- opts[names(opts) %in% .EMC_HOST_OPTIONS |
                 grepl("^(emc2|EMC2)\\.", names(opts))]
  seed <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv()) else NULL
  list(fun = fun, args = args, globals = found$globals,
       packages = found$packages, options = opts, wd = getwd(),
       lib = .libPaths(),
       ns_lib = dirname(getNamespaceInfo(asNamespace("EMC2"), "path")),
       rng_kind = RNGkind(), seed = seed)
}

# Copy whatever the host has written to `con` since the last call, up to the
# last complete line, to the session's stdout (`err = FALSE`) or as messages.
.emc_host_relay <- function(state, err) {
  con <- state$con
  if (is.null(con)) {
    if (!file.exists(state$path)) return(state)
    con <- file(state$path, "rb")
    state$con <- con
  }
  repeat {
    chunk <- readBin(con, "raw", 65536L)
    if (!length(chunk)) break
    state$pending <- c(state$pending, chunk)
  }
  if (!length(state$pending)) return(state)
  ends <- which(state$pending == as.raw(10L) | state$pending == as.raw(13L))
  if (!length(ends) && !isTRUE(state$flush)) return(state)
  last <- if (isTRUE(state$flush)) length(state$pending) else max(ends)
  text <- rawToChar(state$pending[seq_len(last)])
  state$pending <- state$pending[-seq_len(last)]
  if (err) message(text, appendLF = FALSE) else cat(text)
  state
}

.emc_host_call <- function(fun, args) {
  dir <- tempfile("emc_fit_")
  dir.create(dir, recursive = TRUE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  payload <- file.path(dir, "payload.rds")
  result <- file.path(dir, "result.rds")
  pidfile <- file.path(dir, "host.pid")
  out <- list(path = file.path(dir, "stdout.log"), con = NULL, pending = raw())
  err <- list(path = file.path(dir, "stderr.log"), con = NULL, pending = raw())
  saveRDS(.emc_host_payload(fun, args), payload, compress = FALSE)
  # Read before the payload, whose closures would load the namespace from
  # wherever the default library paths find it.
  writeLines(dirname(getNamespaceInfo(asNamespace("EMC2"), "path")),
             file.path(dir, "lib.txt"))

  expr <- paste0(
    "b<-commandArgs(TRUE)[1L];",
    ".libPaths(unique(c(readLines(file.path(b,'lib.txt')),.libPaths())));",
    "get('.emc_host_main',loadNamespace('EMC2'))(b)"
  )
  # The shell writes its pid, which exec keeps for the host, before R starts,
  # so a host that dies while booting is noticed too.
  cmd <- paste(
    "echo $$ >", shQuote(pidfile), "&& exec env EMC2_IN_HOST=1",
    paste0("EMC2_CPP_CACHE=", shQuote(.emc_cpp_cache_dir())),
    shQuote(file.path(R.home("bin"), "Rscript")), "--vanilla",
    "-e", shQuote(expr), shQuote(dir),
    ">", shQuote(out$path), "2>", shQuote(err$path)
  )
  pid <- NA_integer_
  read_pid <- function() {
    if (!file.exists(pidfile)) return(NA_integer_)
    p <- suppressWarnings(as.integer(readLines(pidfile, n = 1L, warn = FALSE)))
    if (length(p) == 1L) p else NA_integer_
  }
  pipe_con <- pipe(cmd, open = "w")
  # Runs before the directory is removed.  A host that has not yet made its
  # group is killed alone; it has started nothing.
  on.exit({
    if (is.na(pid)) pid <- read_pid()
    if (!is.na(pid) && !emc_kill_group(pid))
      try(tools::pskill(pid, tools::SIGKILL), silent = TRUE)
    try(close(pipe_con), silent = TRUE)
    for (s in list(out, err)) if (!is.null(s$con)) try(close(s$con), silent = TRUE)
  }, add = TRUE, after = FALSE)

  wait <- 0.02
  repeat {
    out <- .emc_host_relay(out, err = FALSE)
    err <- .emc_host_relay(err, err = TRUE)
    if (file.exists(result)) break
    if (is.na(pid)) pid <- read_pid()
    if (!is.na(pid) && !.emc_wpool_pid_alive(pid)) {
      if (file.exists(result)) break
      out$flush <- err$flush <- TRUE
      out <- .emc_host_relay(out, err = FALSE)
      err <- .emc_host_relay(err, err = TRUE)
      stop("the fit process ended without a result (killed, out of memory, ",
           "or crashed)", call. = FALSE)
    }
    Sys.sleep(wait)
    wait <- min(wait * 1.5, 0.05)
  }
  # The result is renamed into place after the host's last output, so this
  # relays everything it wrote.
  out$flush <- err$flush <- TRUE
  out <- .emc_host_relay(out, err = FALSE)
  err <- .emc_host_relay(err, err = TRUE)
  res <- readRDS(result)
  if (!is.null(res$seed)) assign(".Random.seed", res$seed, envir = globalenv())
  for (w in res$warnings) warning(w, call. = FALSE)
  if (!is.null(res$error)) stop(res$error, call. = FALSE)
  res$value
}

.emc_host_main <- function(dir) {
  pgid <- emc_host_detach()
  if (identical(pgid, 0L)) {
    # Without a group of its own the host could not end its children with the
    # session; it does not run the fit.
    res <- list(error = "the fit process could not start its process group")
  } else {
    p <- readRDS(file.path(dir, "payload.rds"))
    file.remove(file.path(dir, "payload.rds"))
    here <- dirname(getNamespaceInfo(asNamespace("EMC2"), "path"))
    if (!identical(normalizePath(here), normalizePath(p$ns_lib)))
      stop("the fit process loaded EMC2 from ", here, ", not ", p$ns_lib)
    .libPaths(unique(c(p$ns_lib, p$lib)))
    setwd(p$wd)
    options(p$options)
    for (pkg in rev(p$packages))
      try(suppressPackageStartupMessages(
        library(pkg, character.only = TRUE)), silent = TRUE)
    for (nm in names(p$globals)) assign(nm, p$globals[[nm]], envir = globalenv())
    do.call(RNGkind, as.list(p$rng_kind))
    if (!is.null(p$seed)) assign(".Random.seed", p$seed, envir = globalenv())
    warnings <- character(0)
    fun <- get(p$fun, envir = asNamespace("EMC2"))
    args <- p$args
    rm(p)
    res <- tryCatch(
      list(value = withCallingHandlers(
        do.call(fun, args),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        })),
      error = function(e) list(error = conditionMessage(e)))
    res$warnings <- warnings
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      res$seed <- get(".Random.seed", envir = globalenv())
  }
  tmp <- file.path(dir, "result.rds.tmp")
  saveRDS(res, tmp, compress = FALSE)
  file.rename(tmp, file.path(dir, "result.rds"))
  invisible(NULL)
}

# Where Rcpp keeps compiled custom kernels.  Shared with the processes a fit
# starts, so each loads the library the session already built instead of
# compiling its own.
.emc_cpp_cache_dir <- function() {
  d <- Sys.getenv("EMC2_CPP_CACHE")
  if (nzchar(d)) return(d)
  getOption("rcpp.cache.dir", tempdir())
}
