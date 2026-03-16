lfun <- function(i, x, p_vector, pname, dadm, use_c) {
  p_vector[pname] <- x[i]
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


profile_plot_test <- function (data, design, p_vector, range = 0.5, layout = NA, p_min = NULL,
                               p_max = NULL, use_par = NULL, n_point = 100, n_cores = 1,use_c = FALSE,
                               round = 3, true_args = list(), figure_title = NULL,
                               figure_title_args = list(), ...)
{
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  dots <- list(...)
  figure_title_settings <- NULL
  if (!is.null(figure_title)) {
    figure_title_settings <- modifyList(list(side = 3, outer = TRUE,
                                             line = 1, cex = 1.2),
                                        figure_title_args)
    line_value <- figure_title_settings$line
    if (!is.numeric(line_value) || length(line_value) != 1 || !is.finite(line_value)) {
      line_value <- 1
    }
    figure_title_settings$line <- line_value
  }
  if (!identical(names(p_min), names(p_max)))
    stop("p_min and p_max should be specified for the same parameters")
  if (!is.null(names(p_min)) & length(p_min) == length(use_par))
    names(p_min) <- use_par
  if (!is.null(names(p_max)) & length(p_max) == length(use_par))
    names(p_max) <- use_par
  if (is.null(use_par))
    use_par <- names(p_vector)
  if (any(is.na(layout))) {
    par(mfrow = coda_setmfrow(Nchains = 1, Nparms = length(use_par),
                              nplots = 1))
  }
  else {
    par(mfrow = layout)
  }
  if (!is.null(figure_title_settings)) {
    requested_top <- figure_title_settings$line + 1
    current_oma <- par("oma")
    if (length(current_oma) < 4) {
      current_oma <- c(current_oma, rep(0, 4 - length(current_oma)))
    }
    if (current_oma[3] < requested_top) {
      par(oma = c(current_oma[1], current_oma[2], requested_top, current_oma[4]))
    }
  }
  # Make margins compact so multi-panel plots fit typical devices
  par(mar = c(3.5, 3.5, 2, 1) + 0.1, mgp = c(2, 0.7, 0), tcl = -0.2)
  if (is.null(dots$dadm)) {
    dadm <- EMC2:::design_model(data, design, verbose = FALSE)
  }
  else {
    dadm <- dots$dadm
  }
  # Helper: transform x values for a single parameter to natural scale using model's transform spec
  model_for_axis <- attr(dadm, "model")()
  transform_to_natural <- function(vals, par_name) {
    base_name <- sub("_.*$", "", par_name)  # strip factor suffixes
    fn <- model_for_axis$transform$func[[base_name]]
    if (is.null(fn)) return(vals)
    switch(as.character(fn),
           identity = vals,
           exp = exp(vals),
           pnorm = pnorm(vals),
           # default to identity if unknown
           vals)
  }
  out <- data.frame(true = rep(NA, length(use_par)), max = rep(NA,
                                                               length(use_par)), miss = rep(NA, length(use_par)))
  rownames(out) <- use_par
  for (p in 1:length(p_vector)) {
    cur_name <- names(p_vector)[p]
    if (cur_name %in% use_par) {
      cur_par <- p_vector[p]
      # Use a relative range when possible, but fall back to an absolute range
      # when the current value is ~0 (otherwise pmin==pmax and the panel collapses).
      scale_val <- abs(cur_par)
      if (!is.finite(scale_val) || scale_val < 1e-6) scale_val <- 1
      pmax_cur <- cur_par + ((range / 2) * scale_val)
      pmin_cur <- cur_par - ((range / 2) * scale_val)
      if (!is.null(p_min)) {
        if (!is.na(p_min[cur_name])) {
          pmin_cur <- p_min[cur_name]
        }
      }
      if (!is.null(p_max)) {
        if (!is.na(p_max[cur_name])) {
          pmax_cur <- p_max[cur_name]
        }
      }
      if (isTRUE(all.equal(pmin_cur, pmax_cur))) {
        pmin_cur <- cur_par - (range / 2)
        pmax_cur <- cur_par + (range / 2)
      }
      x <- seq(pmin_cur, pmax_cur, length.out = n_point)
      x <- c(x, cur_par)
      x <- unique(sort(x))
      ll <- unlist(mclapply(1:length(x), lfun, dadm = dadm, use_c = use_c,
                            x = x, p_vector = p_vector, pname = cur_name,
                            mc.cores = n_cores))
      # Plot on natural scale for the x-axis
      x_nat <- transform_to_natural(x, cur_name)
      cur_nat <- transform_to_natural(cur_par, cur_name)
      do.call(plot, c(list(x_nat, ll), EMC2:::fix_dots_plot(EMC2:::add_defaults(dots,
                                                                                type = "l", xlab = paste0(cur_name, " (natural)"), ylab = "LL"))))
      do.call(abline, c(list(v = cur_nat), EMC2:::fix_dots_plot(EMC2:::add_defaults(true_args,
                                                                                    lty = 2))))
      out[cur_name, ] <- c(p_vector[cur_name], x[which.max(ll)],
                           p_vector[cur_name] - x[which.max(ll)])
    }
  }
  if (!is.null(figure_title_settings)) {
    do.call(mtext, c(list(text = figure_title), figure_title_settings))
  }
  return(round(out, 3))
}
