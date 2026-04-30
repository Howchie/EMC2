library(dplyr)
library(ggplot2)
.to_correct_indicator <- function(cr, n) {
  if (is.null(cr)) {
    return(rep(TRUE, n))
  }
  out <- rep(FALSE, n)
  m <- min(length(cr), n)
  if (m > 0) {
    vals <- cr[seq_len(m)]
    out[seq_len(m)] <- !is.na(vals) & as.logical(vals)
  }
  out
}
row_to_pvec <- function(df_row) {
  v <- as.numeric(df_row[1, , drop = TRUE])
  names(v) <- names(df_row)
  v
}
.safe_ecdf <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(function(t) rep(0, length(t)))
  }
  ecdf(x)
}


## Helpers
build_ct_df <- function(times, series_specs, probs = c(.05, .95)) {
  series_df <- lapply(series_specs, function(spec) {
    ct_vals <- spec$fn(times)
    if (!is.null(spec$rt)) {
      bounds <- quantile(spec$rt, probs = probs, names = FALSE)
      ct_vals[times < bounds[1] | times > bounds[2]] <- NA_real_
    }
    data.frame(
      Time = times,
      Ct = ct_vals,
      Series = rep(spec$label, length(times))
    )
  })
  dplyr::bind_rows(series_df)
}

smooth_one_cdf <- function(t, y, smooth_cdf="none", smooth_spar = 0.65) {
  if (smooth_cdf == "none") return(y)
  ord <- order(t)
  t_ord <- t[ord]
  y_ord <- y[ord]
  t_u <- unique(t_ord)
  y_u <- y_ord[match(t_u, t_ord)]
  y_u <- pmin(1, pmax(0, cummax(y_u)))
  if (length(t_u) < 3) {
    y_s <- y_ord
  } else if (smooth_cdf == "mono") {
    sf <- splinefun(t_u, y_u, method = "monoH.FC")
    y_s <- sf(t_ord)
  } else {
    sp <- smooth_spar
    if (!is.finite(sp)) sp <- 0.65
    sp <- max(0.4, min(1.0, sp))
    fit <- tryCatch(stats::smooth.spline(t_ord, y_ord, spar = sp), error = function(e) NULL)
    if (is.null(fit)) {
      sf <- splinefun(t_u, y_u, method = "monoH.FC")
      y_s <- sf(t_ord)
    } else {
      y_s <- as.numeric(stats::predict(fit, x = t_ord)$y)
    }
  }
  y_s <- pmin(1, pmax(0, y_s))
  y_s <- cummax(y_s)
  out <- y
  out[ord] <- y_s
  out
}
## Helpers
build_sic_df <- function(times, series_specs,
                         trim_cdf_tails = TRUE, cdf_tail_cut = 5e-4,
                         smooth_cdf = c("none", "mono", "spline"),
                         smooth_spar = 0.65) {
  smooth_cdf <- match.arg(smooth_cdf)
  series_df <- lapply(series_specs, function(spec) {
    hh_cdf_raw <- spec$fn$HH(times)
    hl_cdf_raw <- spec$fn$HL(times)
    lh_cdf_raw <- spec$fn$LH(times)
    ll_cdf_raw <- spec$fn$LL(times)
    hh_cdf <- smooth_one_cdf(times, hh_cdf_raw, smooth_cdf, smooth_spar)
    hl_cdf <- smooth_one_cdf(times, hl_cdf_raw, smooth_cdf, smooth_spar)
    lh_cdf <- smooth_one_cdf(times, lh_cdf_raw, smooth_cdf, smooth_spar)
    ll_cdf <- smooth_one_cdf(times, ll_cdf_raw, smooth_cdf, smooth_spar)

    sic_vals_raw <- lh_cdf_raw + hl_cdf_raw - hh_cdf_raw - ll_cdf_raw
    sic_vals_smooth <- lh_cdf + hl_cdf - hh_cdf - ll_cdf
    out <- data.frame(
      Time = times,
      SICt = if (smooth_cdf == "none") sic_vals_raw else sic_vals_smooth,
      SICt_raw = sic_vals_raw,
      SICt_smooth = sic_vals_smooth,
      HH_s = 1 - hh_cdf,
      HL_s = 1 - hl_cdf,
      LH_s = 1 - lh_cdf,
      LL_s = 1 - ll_cdf,
      Series = rep(spec$label, length(times))
    )

    if (isTRUE(trim_cdf_tails)) {
      cut <- cdf_tail_cut
      if (!is.finite(cut) || cut <= 0 || cut >= 0.5) cut <- 5e-4
      left_tail <- (hh_cdf < cut) & (hl_cdf < cut) & (lh_cdf < cut) & (ll_cdf < cut)
      right_tail <- (hh_cdf > (1 - cut)) & (hl_cdf > (1 - cut)) &
        (lh_cdf > (1 - cut)) & (ll_cdf > (1 - cut))
      out <- out[!(left_tail | right_tail), , drop = FALSE]
    }
    out
  })
  dplyr::bind_rows(series_df)
}

build_cdf_df <- function(times, series_specs, probs = c(.05, .95),
                         smooth_cdf = c("none", "mono", "spline"),
                         smooth_spar = 0.65) {
  smooth_cdf <- match.arg(smooth_cdf)
  series_df <- lapply(series_specs, function(spec) {
    if (is.null(spec$rt)) stop("Each series spec must include an `rt` vector.")
    keep <- if (is.null(spec$cr)) rep(TRUE, length(spec$rt)) else .to_correct_indicator(spec$cr, length(spec$rt))
    rt_use <- spec$rt[keep]
    cdf_raw <- .safe_ecdf(rt_use)(times)
    cdf_smooth <- smooth_one_cdf(times, cdf_raw, smooth_cdf, smooth_spar)
    cdf_vals <- if (smooth_cdf == "none") cdf_raw else cdf_smooth
    if (length(rt_use) > 0) {
      bounds <- quantile(rt_use, probs = probs, names = FALSE, na.rm = TRUE)
      cdf_vals[times < bounds[1] | times > bounds[2]] <- NA_real_
      cdf_raw[times < bounds[1] | times > bounds[2]] <- NA_real_
      cdf_smooth[times < bounds[1] | times > bounds[2]] <- NA_real_
    } else {
      cdf_vals[] <- NA_real_
      cdf_raw[] <- NA_real_
      cdf_smooth[] <- NA_real_
    }
    data.frame(
      Time = times,
      CDF = cdf_vals,
      CDF_raw = cdf_raw,
      CDF_smooth = cdf_smooth,
      Series = rep(spec$label, length(times))
    )
  })
  dplyr::bind_rows(series_df)
}

get_time_range <- function(rt_list, probs = c(.05, .95), by = .001) {
  if (!is.list(rt_list)) {
    rt_list <- list(rt_list)
  }
  quantiles <- lapply(rt_list, quantile, probs = probs, names = FALSE)
  lower <- max(vapply(quantiles, `[`, numeric(1), 1), na.rm = TRUE)
  upper <- min(vapply(quantiles, `[`, numeric(1), 2), na.rm = TRUE)
  if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
    stop("No usable quantile window for the provided RTs.")
  }
  seq(lower, upper, by = by)
}



plot_altieri = function(n,p_vec, smooth_cdf = c("none", "mono", "spline"), smooth_spar = 0.65) {
  smooth_cdf <- match.arg(smooth_cdf)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_altieri() requires the 'ggplot2' package.")
  }
  # NB -- to correctly illustrate UCIP predictions in AND tasks we set vE very low.
  sims <- simulate_sft(
    model = "lba",
    n = n,
    p_vec = p_vec,
    design = c("DT", "A", "B", "NT"),
    logical_rules = "AND"
  )
  
  ## Construct Plots
  
  # 1. Demonstrate non-equivalence of Ct_AND and Ct_Altieri and that the three existing functions all show C(t)~1 for their respective design
  and_data <- sims$by_rule$AND
  cell_labels <- c("AB", "AN", "NB", "NN")
  cell_rts <- lapply(cell_labels, function(cell) {
    idx <- and_data$S == cell & and_data$Correct
    and_data$rt[idx]
  })
  if (any(vapply(cell_rts, length, integer(1)) == 0L)) {
    cell_rts <- lapply(cell_labels, function(cell) and_data$rt[and_data$S == cell])
  }
  times <- get_time_range(cell_rts)
  plot_df <- build_ct_df(
    times,
    list(
      list(label = "C_AND(t)", fn = sims$metrics$AND_Ct$Ct, rt = sims$by_rule$AND$rt),
      list(label = "C_Absence(t)", fn = sims$metrics$AND_Absence_Ct$Ct, rt = sims$by_rule$AND$rt),
      list(label = "C_Altieri(t)", fn = sims$metrics$Altieri_Ct$Ct_a1, rt = sims$by_rule$AND$rt),
      list(label = "C_Bound(t)", fn = sims$metrics$Altieri_Ct$Ct_a2, rt = sims$by_rule$AND$rt)
    )
  )
  plot <- ggplot(plot_df, aes(x = Time, y = Ct, color = Series, linetype=Series)) +
    geom_line(linewidth = 1, alpha=.7) +
    scale_color_manual(values = c("C_AND(t)" = "darkgreen", "C_Absence(t)" = "blue", "C_Altieri(t)" = "black", "C_Bound(t)" = "red")) +
    scale_linetype_manual(values = c("C_AND(t)" = "solid", "C_Absence(t)" = "dashed", "C_Altieri(t)" = "solid", "C_Bound(t)" = "dashed")) +
    coord_cartesian(ylim = c(0, 5)) +
    labs(x = "Time", y = "C(t)") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme_classic() + theme(axis.ticks.x=element_blank(),axis.text.x = element_blank()) +
    ggtitle("Comparison of Altieri forms and absence coefficient")
  
  cdf_df <- build_cdf_df(
      times,
      list(
      list(label = "AB", rt = sims$by_rule$AND$rt[sims$by_rule$AND$S=="AB"], cr = sims$by_rule$AND$Correct[sims$by_rule$AND$S=="AB"]),
      list(label = "NN", rt = sims$by_rule$AND$rt[sims$by_rule$AND$S=="NN"], cr = sims$by_rule$AND$Correct[sims$by_rule$AND$S=="NN"])
      ),
    smooth_cdf = smooth_cdf,
    smooth_spar = smooth_spar
  )
  cdf_plot <- ggplot(cdf_df, aes(x = Time, y = CDF, color = Series, linetype = Series)) +
    geom_line(linewidth = 1, alpha = .7) +
    scale_color_manual(values = c("AB" = "blue", "NN" = "red")) +
    scale_linetype_manual(values = c("AB" = "solid", "NN" = "solid")) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Time", y = "CDF") +
    theme_classic() + theme(axis.ticks.x=element_blank(),axis.text.x = element_blank()) +
    ggtitle("AND CDF Comparison: AB vs NN")
  return(list(sims=sims,plot_df=plot_df,cdf_df=cdf_df,plot=plot,cdf_plot=cdf_plot,p_vec=p_vec))
}

## Short helper to minimise the distance between a min and a max-time distribution
obj_max_vs_min <- function(vy, vno, t, A, b, t0, sv) {
  Fy <- rtdists::plba_norm(t, A=A, b=b, t0=t0, mean_v=vy, sd_v=sv, posdrift=TRUE)
  Fn <- rtdists::plba_norm(t, A=A, b=b, t0=t0, mean_v=vno, sd_v=sv, posdrift=TRUE)
  Fmax <- Fy^2
  Fmin <- 1 - (1 - Fn)^2
  mean((Fmax - Fmin)^2)
}

## Find vc_yes that makes AB(max) ~ NN(min) for a fixed NN drift
find_equal_vc_yes <- function(vno, t, A, b, t0, sv, interval = c(1.0, 4.5)) {
  optimize(
    obj_max_vs_min,
    interval = interval,
    vno = vno,
    t = t,
    A = A,
    b = b,
    t0 = t0,
    sv = sv
  )$minimum
}

.normalize_logical_rules <- function(logical_rules) {
  if (is.null(logical_rules)) {
    logical_rules <- c("OR", "AND", "ID")
  }
  logical_rules <- unique(logical_rules)
  bad <- setdiff(logical_rules, c("OR", "AND", "XOR", "ID"))
  if (length(bad) > 0) {
    stop(sprintf("Unsupported logical rules: %s", paste(bad, collapse = ", ")))
  }
  logical_rules
}

.default_design <- function() {
  c("AB", "AN", "NB", "NN")
}

.char_to_level <- function(ch) {
  if (ch %in% c("A", "B", "H")) return(2L)
  if (ch == "L") return(1L)
  if (ch %in% c("N", "_")) return(0L)
  stop(sprintf("Unknown design character '%s'.", ch))
}

.annotate_trials <- function(tmp, logical_rule="OR") {
  meta_row = .expand_design()
  idx <- match(as.character(tmp$S), as.character(meta_row$S))
  tmp$Condition <- logical_rule
  tmp$Channel1 <- meta_row$Channel1[idx]
  tmp$Channel2 <- meta_row$Channel2[idx]
  tmp
}

.expand_design <- function(design=NULL) {
  
  if (is.null(design)) {
    design <- c("AB", "AN", "NB", "NN")
  }
  
  out <- lapply(seq_along(design), function(i) {
    cell <- design[[i]]
    if (nchar(cell) != 2) {
      stop(sprintf("Design cell '%s' is invalid. Cells must be two characters.", cell))
    }
    c1 <- substr(cell, 1, 1)
    c2 <- substr(cell, 2, 2)
    l1 <- .char_to_level(c1)
    l2 <- .char_to_level(c2)
    s <- if (l1 > 0L && l2 > 0L) {
      "AB"
    } else if (l1 > 0L && l2 == 0L) {
      "AN"
    } else if (l1 == 0L && l2 > 0L) {
      "NB"
    } else {
      "NN"
    }
    data.frame(
      cell_raw = cell,
      S = s,
      Channel1 = l1,
      Channel2 = l2,
      A_present = l1 > 0L,
      B_present = l2 > 0L,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(out)
}

.split_by_rule <- function(data, logical_rules) {
  # Pre-split data by Condition for O(N) performance instead of O(N*S)
  split_data <- split(data, data$Condition)
  out <- list()
  for (r in logical_rules) {
    res <- split_data[[r]]
    if (is.null(res)) {
      out[[r]] <- data[0, , drop = FALSE]
    } else {
      out[[r]] <- res
    }
  }
  out
}

