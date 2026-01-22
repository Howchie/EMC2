#remotes::install_github("https://github.com/Howchie/EMC2/",ref="censoring_truncation_ah")

rm(list=ls())
library(EMC2)
library(microbenchmark)

# helper function to check system memory during fit (not sure but it looks accurate to manually monitoring)
profile_fit_rss_live <- function(expr, label = "fit", interval_sec = 0.2, keep_log = FALSE, log_path = NULL,
                                 include_children = TRUE) {
  status_path <- sprintf("/proc/%d/status", Sys.getpid())
  if (!file.exists(status_path)) {
    before <- gc(full = TRUE)[, "used"]
    out <- NULL
    err <- NULL
    t <- system.time({
      out <- tryCatch(
        eval.parent(substitute(expr)),
        error = function(e) {
          err <<- e
          NULL
        }
      )
    })
    after <- gc(full = TRUE)[, "used"]
    delta_ncells <- unname(after["Ncells"] - before["Ncells"])
    delta_vcells <- unname(after["Vcells"] - before["Vcells"])
    cat(sprintf(
      "[%s] elapsed=%.2fs delta_used(Ncells)=%s delta_used(Vcells)=%s (/proc unavailable)\n",
      label,
      unname(t[["elapsed"]]),
      format(delta_ncells, scientific = FALSE, big.mark = ","),
      format(delta_vcells, scientific = FALSE, big.mark = ",")
    ))
    if (!is.null(err)) stop(err)
    return(invisible(out))
  }

  read_kb <- function(key, lines) {
    hit <- grep(paste0("^", key, ":"), lines, value = TRUE)
    if (!length(hit)) return(NA_real_)
    as.numeric(gsub("[^0-9.]", "", hit[1]))
  }

  read_ppid <- function(lines) {
    hit <- grep("^PPid:", lines, value = TRUE)
    if (!length(hit)) return(NA_integer_)
    as.integer(gsub("[^0-9]", "", hit[1]))
  }

  snapshot_proc_tree_kb <- function(root_pid) {
    proc_dirs <- list.files("/proc", pattern = "^[0-9]+$")
    if (!length(proc_dirs)) return(list(n = 0L, rss_kb = NA_real_, hwm_kb = NA_real_))

    pids <- as.integer(proc_dirs)
    ppid <- rep(NA_integer_, length(pids))
    rss <- rep(NA_real_, length(pids))
    hwm <- rep(NA_real_, length(pids))

    for (i in seq_along(pids)) {
      lines <- tryCatch(
        readLines(sprintf("/proc/%d/status", pids[i]), warn = FALSE),
        error = function(e) character()
      )
      if (!length(lines)) next
      ppid[i] <- read_ppid(lines)
      rss[i] <- read_kb("VmRSS", lines)
      hwm[i] <- read_kb("VmHWM", lines)
    }

    if (!isTRUE(include_children)) {
      idx <- which(pids == root_pid)
      if (!length(idx)) return(list(n = 0L, rss_kb = NA_real_, hwm_kb = NA_real_))
      return(list(n = 1L, rss_kb = rss[idx[1]], hwm_kb = hwm[idx[1]]))
    }

    seen <- integer(0)
    queue <- as.integer(root_pid)

    while (length(queue)) {
      cur <- queue[1]
      queue <- queue[-1]
      if (cur %in% seen) next
      seen <- c(seen, cur)
      kids <- pids[which(ppid == cur)]
      if (length(kids)) queue <- c(queue, kids)
    }

    in_tree <- which(pids %in% seen)
    list(
      n = length(in_tree),
      rss_kb = sum(rss[in_tree], na.rm = TRUE),
      hwm_kb = sum(hwm[in_tree], na.rm = TRUE)
    )
  }

  tmp <- if (is.null(log_path)) tempfile(fileext = ".rss.tsv") else log_path
  stopflag <- tempfile(fileext = ".stop")
  root_pid <- Sys.getpid()

  monitor <- parallel::mcparallel({
    start <- Sys.time()
    writeLines("t_sec\tn_proc\ttree_VmRSS_kB\ttree_VmHWM_kB", tmp)
    repeat {
      if (file.exists(stopflag)) break
      snap <- snapshot_proc_tree_kb(root_pid)
      tsec <- as.numeric(difftime(Sys.time(), start, units = "secs"))
      cat(sprintf("%.3f\t%d\t%s\t%s\n", tsec, snap$n, snap$rss_kb, snap$hwm_kb), file = tmp, append = TRUE)
      Sys.sleep(interval_sec)
    }
    invisible(TRUE)
  })

  on.exit({
    file.create(stopflag)
    parallel::mccollect(monitor, wait = TRUE)
    unlink(stopflag)
    if (!isTRUE(keep_log) && is.null(log_path)) unlink(tmp)
  }, add = TRUE)

  out <- NULL
  err <- NULL
  t <- system.time({
    out <- tryCatch(
      eval.parent(substitute(expr)),
      error = function(e) {
        err <<- e
        NULL
      }
    )
  })

  rss_df <- tryCatch(read.table(tmp, header = TRUE, sep = "\t"), error = function(e) NULL)
  max_rss <- if (!is.null(rss_df)) max(rss_df$tree_VmRSS_kB, na.rm = TRUE) else NA_real_
  max_hwm <- if (!is.null(rss_df)) max(rss_df$tree_VmHWM_kB, na.rm = TRUE) else NA_real_

  log_msg <- if (isTRUE(keep_log) || !is.null(log_path)) tmp else "(deleted)"
  cat(sprintf(
    "[%s] elapsed=%.2fs peak_tree_VmRSS=%.0f kB peak_tree_VmHWM=%.0f kB (log=%s)\n",
    label,
    unname(t[["elapsed"]]),
    max_rss,
    max_hwm,
    log_msg
  ))

  if (!is.null(err)) stop(err)
  invisible(out)
}

## This is a helper function to directly check likelihoods
lfun <- function(p_vector, dadm, use_c=T) {
  if (use_c) {
    p_matrix <- matrix(p_vector,nrow=1)
    colnames(p_matrix) <- names(p_vector)
    model <- attr(dadm, "model")()
    p_types=names(model$p_types)
    designs <- list()
    for (p in p_types) {
      designs[[p]] <- attr(dadm,"designs")[[p]][attr(attr(dadm,"designs")[[p]],"expand"),,drop=FALSE]
    }
    constants <- attr(dadm,"constants")
    if (is.null(constants)) constants <- NA
    EMC2:::calc_ll(p_matrix, dadm, constants,designs,model$c_name,
                   model$bound,model$transform,model$pre_transform,p_types,log(1e-10),model$trend)
  } else {
    EMC2:::calc_ll_R(p_vector, attr(dadm, "model")(), dadm)
  }
}


# Here is a simple 5 parameter LBA model with contamination that uses Zach's C
designLBA <- design(
  factors=list(subjects=1,S=c("left","right")),Rlevels=c("left","right"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  model=LBA,
  formula=list(v~lM,B~1,t0~1,pContaminant~1)
)

# Old C likelihood implementation (for memory/speed comparisons)
designOLBA <- design(
  factors=list(subjects=1,S=c("left","right")),Rlevels=c("left","right"),
  matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
  model=OLBA,
  formula=list(v~lM,B~1,t0~1,pContaminant~1)
)

# This parameter vector produces reasonable RTs and accuracies with no contamination (old_c can't do that)
p_vector <- sampled_pars(designLBA,doMap = FALSE)
p_vector[] <- c(1, 1.25, log(2), log(0.2),qnorm(.0))

TC=list(LT=0,UT=Inf,LC=0,UC=Inf,verbose=TRUE)
datCT <- make_data(p_vector,designLBA,n_trials=10000,TC=TC)

# Compare likelihoods
dadmC <- EMC2:::design_model(datCT,designLBA)
dadmO <- EMC2:::design_model(datCT,designOLBA)
lfun(p_vector,dadmC)
lfun(p_vector,dadmO)

# Fitting -- check memory usage and speed equivalence
emc <- make_emc(datCT,designLBA,type="single", compress=TRUE)
emccCT <- profile_fit_rss_live(fit(emc), label = "fit: LBA (C)")
# Time difference of 1.165462 mins
# [fit: LBA (C)] elapsed=69.93s peak_tree_VmRSS=4061392 kB peak_tree_VmHWM=4098636 kB (log=(deleted))

# fit old C (OLBA) with the same data, to compare memory behaviour.
emc <- make_emc(datCT,designOLBA,type="single")
emcoCT <- profile_fit_rss_live(fit(emc), label = "fit: OLBA (Old_C)")
# Time difference of 1.135657 mins
# [fit: OLBA (Old_C)] elapsed=68.14s peak_tree_VmRSS=4033860 kB peak_tree_VmHWM=4035188 kB (log=(deleted))

#### Microbench: LBA vs OLBA fit() timing (10 iters) to average out noise in sampling ----
fit_LBA_C <- function() {
  emc <- make_emc(datCT, designLBA, type = "single")
  fit(emc)
}

fit_OLBA_OldC <- function() {
  emc <- make_emc(datCT, designOLBA, type = "single")
  fit(emc)
}

bench_fit_10 <- microbenchmark::microbenchmark(
    LBA_C = fit_LBA_C(),
    OLBA_OldC = fit_OLBA_OldC(),
    times = 10L
  )
print(bench_fit_10)

# Unit: seconds
# expr      min       lq     mean   median       uq      max neval
# LBA_C 58.87708 71.12796 76.88685 75.57930 80.36391 95.86790    10
# OLBA_OldC 66.35529 76.76941 79.94931 79.64413 85.73642 87.16758    10
