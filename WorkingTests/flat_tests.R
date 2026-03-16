# ----------------------------
# Flatness tests + EMC2 helpers
# ----------------------------
# The basic idea is is:
# 1.	Create an object interleaving the three chains together. So for example ch1=a1,a2, ch2=b1,b2, ch3=c1,c2, then the interleaved result is a1,b1,c1,a2,b2,c2. This guarantees the final number of samples in this object always divisible by 3
# 2.	Check the ratio of absolute difference in median of first and last third of this interleaved chains, relative to its IQR. This is the stationarity test. By default checks every mu, alpha and subj LL for stationarity.
# 3.	Loop fit() until every value of stat above less than some cutoff (cut.location), meaning quite flat
# 4.	While doing so, try to keep samples under maxiter. Dump some initial samples if fit is projected to break maxiter. Having a maxiter seems to help because dumping the initial tail can help to g et flat, if early samples are not yet stationary. 
# Untested (but included) option is to also test theta var and covar (upper triangle) with same tests. Also, there is untested code to compare the variability across first and third chains, rather than just median. 

# ---- Core: location drift stats per column ----
flat_stats_per_col <- function(mat, p1 = 1/3, p2 = 1/3) {
  stopifnot(is.matrix(mat), nrow(mat) > 10)
  
  N <- nrow(mat)
  xlen <- max(1L, round(N * p1))
  ylen <- max(1L, round(N * p2))
  
  apply(mat, 2, function(x) {
    denom <- IQR(x)
    if (denom == 0) return(0)
    m1 <- median(x[1:xlen])
    m2 <- median(x[(N - ylen + 1):N])
    abs(m1 - m2) / denom
  })
}

# ---- Core: scale drift stats per column (optional) ----
scale_stats_per_col <- function(mat, p1 = 1/3, p2 = 1/3) {
  stopifnot(is.matrix(mat), nrow(mat) > 10)
  
  N <- nrow(mat)
  xlen <- max(1L, round(N * p1))
  ylen <- max(1L, round(N * p2))
  
  apply(mat, 2, function(x) {
    denom <- IQR(x)
    if (denom == 0) return(0)
    s1 <- IQR(x[1:xlen])
    s2 <- IQR(x[(N - ylen + 1):N])
    abs(s1 - s2) / denom
  })
}

# ---- Wrapper: Flat check with names + reporting ----
flat_check_mat <- function(mat,
                           p1 = 1/3,
                           p2 = 1/3,
                           cut.location = 0.25,
                           cut.scale = Inf,
                           verbose = FALSE) {
  
  stopifnot(is.matrix(mat), nrow(mat) > 10)
  
  nm <- colnames(mat)
  if (is.null(nm)) nm <- as.character(seq_len(ncol(mat)))
  
  m.zs <- flat_stats_per_col(mat, p1 = p1, p2 = p2)
  names(m.zs) <- nm
  fail_loc <- any(m.zs > cut.location)
  
  s.zs <- NULL
  fail_scale <- FALSE
  if (is.finite(cut.scale)) {
    s.zs <- scale_stats_per_col(mat, p1 = p1, p2 = p2)
    names(s.zs) <- nm
    fail_scale <- any(s.zs > cut.scale)
  }
  
  fail <- fail_loc || fail_scale
  
  if (verbose) {
    cat("Flat check\n")
    print(round(m.zs, 3))
    if (!is.null(s.zs)) print(round(s.zs, 3))
    
    if (!fail) {
      cat(": OK\n")
    } else {
      out <- character(0)
      if (fail_loc) {
        w <- which.max(m.zs)
        out <- c(out, paste0(names(m.zs)[w], "=", round(m.zs[w], 2)))
      }
      if (fail_scale) {
        w <- which.max(s.zs)
        out <- c(out, paste0(names(s.zs)[w], "=", round(s.zs[w], 2)))
      }
      cat(":", paste(out, collapse = ", "), "\n")
    }
  }
  
  invisible(list(
    fail = fail,
    m.zs = m.zs,
    s.zs = s.zs,
    xlen = max(1L, round(nrow(mat) * p1)),
    ylen = max(1L, round(nrow(mat) * p2))
  ))
}

# ----------------------------
# EMC2 pooled (interleaved) matrices
# ----------------------------

# subj_ll is [subj, iter] per chain -> pooled mat is [chain*iter interleaved, subj]
subj_ll_flatmat <- function(emc) {
  C <- length(emc)
  stopifnot(C >= 1)
  
  sll1 <- emc[[1]]$samples$subj_ll
  stopifnot(is.matrix(sll1))
  
  S <- nrow(sll1)
  I <- min(vapply(emc, function(ch) ncol(ch$samples$subj_ll), integer(1)))
  
  subj_ids <- rownames(sll1)
  if (is.null(subj_ids)) subj_ids <- as.character(seq_len(S))
  
  # chain × iter × subj
  arr <- array(NA_real_, dim = c(C, I, S))
  for (c in seq_len(C)) {
    sll <- emc[[c]]$samples$subj_ll[, seq_len(I), drop = FALSE]
    arr[c, , ] <- t(sll)  # iter × subj
  }
  
  mat <- matrix(arr, ncol = S)       # rows = chain×iter (interleaved), cols = subj
  colnames(mat) <- subj_ids
  mat
}

# alpha is [param, subj, iter] per chain -> pooled mat is [chain*iter interleaved, param×subj]
alpha_pooled_mat <- function(emc) {
  C <- length(emc)
  stopifnot(C >= 1)
  
  a1 <- emc[[1]]$samples$alpha
  stopifnot(length(dim(a1)) == 3)
  
  P <- dim(a1)[1]
  S <- dim(a1)[2]
  I <- min(vapply(emc, function(ch) dim(ch$samples$alpha)[3], integer(1)))
  
  param_names <- dimnames(a1)[[1]]
  subj_ids    <- dimnames(a1)[[2]]
  if (is.null(param_names)) param_names <- paste0("par", seq_len(P))
  if (is.null(subj_ids))    subj_ids    <- paste0("subj", seq_len(S))
  
  # chain × iter × (param×subj)
  arr <- array(NA_real_, dim = c(C, I, P * S))
  
  for (c in seq_len(C)) {
    alpha_c <- emc[[c]]$samples$alpha[, , seq_len(I), drop = FALSE]  # P × S × I
    alpha_c <- aperm(alpha_c, c(3, 1, 2))                            # I × P × S
    dim(alpha_c) <- c(I, P * S)                                      # I × (P*S)
    arr[c, , ] <- alpha_c
  }
  
  mat <- matrix(arr, ncol = P * S)  # rows = chain×iter (interleaved)
  colnames(mat) <- as.vector(outer(param_names, subj_ids, paste, sep = "_"))
  mat
}

# ----------------------------
# Convenience: group LL pooled vector (sum over subjects), if you want it
# ----------------------------
group_ll_pooled_vec <- function(emc) {
  C <- length(emc)
  stopifnot(C >= 1)
  
  I <- min(vapply(emc, function(ch) ncol(ch$samples$subj_ll), integer(1)))
  arr <- matrix(NA_real_, nrow = C, ncol = I)
  
  for (c in seq_len(C)) {
    arr[c, ] <- colSums(emc[[c]]$samples$subj_ll[, seq_len(I), drop = FALSE])
  }
  
  as.vector(arr)  # c1t1,c2t1,... interleaved by iter
}



theta_mu_pooled_mat <- function(emc) {
  C <- length(emc)
  stopifnot(C >= 1)
  
  mu1 <- emc[[1]]$samples$theta_mu
  stopifnot(is.matrix(mu1))
  
  P <- nrow(mu1)
  I <- min(vapply(emc, function(ch) ncol(ch$samples$theta_mu), integer(1)))
  
  param_names <- rownames(mu1)
  if (is.null(param_names)) param_names <- paste0("par", seq_len(P))
  
  # chain × iter × param
  arr <- array(NA_real_, dim = c(C, I, P))
  
  for (c in seq_len(C)) {
    mu_c <- emc[[c]]$samples$theta_mu[, seq_len(I), drop = FALSE]  # P × I
    arr[c, , ] <- t(mu_c)  # iter × param
  }
  
  mat <- matrix(arr, ncol = P)  # interleaved chain×iter
  colnames(mat) <- param_names
  mat
}



theta_var_upper_pooled_mat <- function(emc) {
  C <- length(emc)
  stopifnot(C >= 1)
  
  var1 <- emc[[1]]$samples$theta_var
  P <- dim(var1)[1]
  I <- min(vapply(emc, function(ch) dim(ch$samples$theta_var)[3], integer(1)))
  
  param_names <- dimnames(var1)[[1]]
  if (is.null(param_names)) param_names <- paste0("par", seq_len(P))
  
  idx <- which(upper.tri(matrix(NA, P, P), diag = TRUE), arr.ind = TRUE)
  K <- nrow(idx)
  
  arr <- array(NA_real_, dim = c(C, I, K))
  
  for (c in seq_len(C)) {
    var_c <- emc[[c]]$samples$theta_var[, , seq_len(I), drop = FALSE]
    for (k in seq_len(K)) {
      i <- idx[k, 1]
      j <- idx[k, 2]
      arr[c, , k] <- var_c[i, j, ]
    }
  }
  
  mat <- matrix(arr, ncol = K)
  
  colnames(mat) <- paste0(
    "cov_",
    param_names[idx[,1]],
    "_",
    param_names[idx[,2]]
  )
  
  mat
}


theta_var_diag_pooled_mat <- function(emc) {
  C <- length(emc)
  stopifnot(C >= 1)
  
  var1 <- emc[[1]]$samples$theta_var
  stopifnot(length(dim(var1)) == 3)
  
  P <- dim(var1)[1]
  I <- min(vapply(emc, function(ch) dim(ch$samples$theta_var)[3], integer(1)))
  
  param_names <- dimnames(var1)[[1]]
  if (is.null(param_names)) param_names <- paste0("par", seq_len(P))
  
  arr <- array(NA_real_, dim = c(C, I, P))
  
  for (c in seq_len(C)) {
    var_c <- emc[[c]]$samples$theta_var[, , seq_len(I), drop = FALSE]
    for (p in seq_len(P)) {
      arr[c, , p] <- var_c[p, p, ]
    }
  }
  
  mat <- matrix(arr, ncol = P)
  colnames(mat) <- paste0("var_", param_names)
  mat
}


flat_check_emc <- function(emc,
                           p1 = 1/3,
                           p2 = 1/3,
                           cut.location = 0.5,
                           verbose = TRUE,
                           include_theta_var = FALSE) {
  
  run_test <- function(mat) {
    out <- flat_check_mat(
      mat,
      p1 = p1, p2 = p2,
      cut.location = cut.location,
      cut.scale = Inf,
      verbose = FALSE
    )
    
    worst_idx <- which.max(out$m.zs)
    list(
      fail      = out$fail,
      m.zs      = out$m.zs,
      n_total   = length(out$m.zs),
      n_fail    = sum(out$m.zs > cut.location),
      prop_fail = mean(out$m.zs > cut.location),
      worst     = names(out$m.zs)[worst_idx],
      worst_val = out$m.zs[worst_idx]
    )
  }
  
  results <- list(
    alpha    = run_test(alpha_pooled_mat(emc)),
    subj_ll  = run_test(subj_ll_flatmat(emc)),
    theta_mu = run_test(theta_mu_pooled_mat(emc))
  )
  
  # optional monitoring only
  if (include_theta_var) {
    results$theta_var <- run_test(theta_var_upper_pooled_mat(emc))
  }
  
  results$overall_fail <- any(vapply(results[c("alpha","subj_ll","theta_mu")],
                                     function(x) isTRUE(x$fail), logical(1)))
  
  if (verbose) {
    cat("\n========================================\n")
    cat("Stationarity summary (cut =", cut.location, ")\n")
    cat("========================================\n")
    for (nm in c("alpha","subj_ll","theta_mu")) {
      r <- results[[nm]]
      cat(sprintf(
        "%-8s : %-4s | %3d/%3d fail (%.1f%%) | worst = %s (%.3f)\n",
        nm,
        ifelse(r$fail, "FAIL", "OK"),
        r$n_fail, r$n_total,
        100*r$prop_fail,
        r$worst, r$worst_val
      ))
    }
    if (include_theta_var) {
      r <- results$theta_var
      cat(sprintf(
        "%-8s : %-4s | %3d/%3d fail (%.1f%%) | worst = %s (%.3f)  [monitor]\n",
        "thetaVar",
        ifelse(r$fail, "FAIL", "OK"),
        r$n_fail, r$n_total,
        100*r$prop_fail,
        r$worst, r$worst_val
      ))
    }
    cat("----------------------------------------\n")
    cat("Overall:", ifelse(results$overall_fail, "FAIL", "OK"), "\n")
    cat("========================================\n")
  }
  
  invisible(results)
}


get_emc_iters <- function(emc) {
  # Use subj_ll if present (usually the safest “iter count”)
  if (!is.null(emc[[1]]$samples$subj_ll)) {
    return(min(vapply(emc, function(ch) ncol(ch$samples$subj_ll), integer(1))))
  }
  # fallback to theta_mu
  if (!is.null(emc[[1]]$samples$theta_mu)) {
    return(min(vapply(emc, function(ch) ncol(ch$samples$theta_mu), integer(1))))
  }
  stop("Can't infer iteration count: no subj_ll or theta_mu in samples.")
}

get_planned_iter <- function(...) {
  dots <- list(...)
  if (!is.null(dots$iter)) return(as.integer(dots$iter))
  1000L  # EMC2 fit() default (per your note)
}


fit_until_stationary <- function(fileName,
                                 max_steps = 20,
                                 cut.location = 0.5,
                                 maxiter = 5000,
                                 log_file = NULL,
                                 verbose = TRUE,
                                 ...) {
  
  if (!is.character(fileName) || length(fileName) != 1L) {
    stop("fileName must be a single character path to the .RData file containing `emc`.")
  }
  
  if (is.null(log_file)) {
    log_file <- sub("\\.RData$", "_stationarity_log.RData", fileName)
  }
  
  history <- list()
  
  for (k in seq_len(max_steps)) {
    
    # ---- Load current emc ----
    load(fileName)
    if (!exists("emc")) stop("Expected object named 'emc' in ", fileName)
    
    cur_iter <- get_emc_iters(emc)
    add_iter <- get_planned_iter(...)
    
    if (verbose) {
      cat("\n\n========================================\n")
      cat("Iteration step", k, "\n")
      cat("Current iters:", cur_iter, " | Planned add:", add_iter, " | Max iters:", maxiter, "\n")
      cat("========================================\n")
    }
    
    # ---- Trim BEFORE fitting (so we don't exceed maxiter after adding draws) ----
    drop <- as.integer((cur_iter + add_iter) - maxiter)
    trimmed <- FALSE
    
    if (drop > 0L) {
      if (drop >= cur_iter) {
        stop("Requested to drop ", drop, " draws but only have ", cur_iter,
             ". Increase maxiter or reduce iter.")
      }
      if (verbose) cat("Trimming early draws: dropping first ", drop, "\n", sep = "")
      emc <- subset(emc, filter = drop)
      save(emc, file = fileName)
      trimmed <- TRUE
    }
    
    # ---- Fit FIRST ----
    if (verbose) cat("Running fit()...\n")
    fit(
      emc,
      fileName = fileName,
      ...
    )
    
    rm(emc)
    gc()
    
    # ---- Reload the updated object and THEN test ----
    load(fileName)
    if (!exists("emc")) stop("Expected object named 'emc' in ", fileName)
    
    diag <- flat_check_emc(
      emc,
      cut.location = cut.location,
      verbose = verbose
    )
    
    diag$meta <- list(
      step = k,
      iters_before = cur_iter,
      planned_add_iter = add_iter,
      maxiter = maxiter,
      trimmed = trimmed,
      dropped = if (trimmed) drop else 0L,
      iters_after = get_emc_iters(emc)
    )
    
    history[[k]] <- diag
    save(history, file = log_file)
    
    if (!diag$overall_fail) {
      if (verbose) cat("\n✅ Stationary achieved at step", k, "\n")
      return(invisible(list(ok = TRUE, history = history)))
    }
    
    rm(emc)
    gc()
  }
  
  if (verbose) cat("\n⚠️ Reached max_steps without full stationarity.\n")
  invisible(list(ok = FALSE, history = history))
}