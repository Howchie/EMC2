#' Two-state accumulator mixture
#'
#' Wraps a race model constructor so that each accumulator runs, with
#' probability `pmix`, in an alternative state in which the parameters named in
#' `pars` change. Accumulators draw their states independently, and an
#' accumulator with `pmix = 0` (the default) is unchanged. The likelihood is
#' exact: an independent race is linear in each accumulator's density and
#' survivor, so mixing those per accumulator is the trial-level mixture.
#'
#' `alt` sets how the alternative state is specified, and names the added
#' parameter:
#'
#' | **alt**   | **Added parameter** | **Alternative value** | **Transform**   |
#' |-----------|---------------------|-----------------------|-----------------|
#' | `"value"` | *`<par>`_mix*       | `<par>_mix`           | as `<par>`      |
#' | `"below"` | *`<par>`_below*     | `<par> - <par>_below` | log, default 0  |
#' | `"above"` | *`<par>`_above*     | `<par> + <par>_above` | log, default 0  |
#'
#' plus *pmix* (probit, natural scale \[0, 1\], default 0), the probability of
#' the alternative state. The gaps are on the natural scale of `<par>`.
#'
#' Not available for LogicalRules, correlated races, stop-signal or DDM-family
#' models.
#'
#' @param model A race model constructor, e.g. `LBA`.
#' @param pars Character vector of parameters that change in the alternative
#'   state.
#' @param alt One of `"value"`, `"below"` or `"above"`.
#' @return A model constructor that can be passed to [design()].
#' @examples
#' LBA2 <- mix_states(LBA, pars = "v", alt = "below")
#' @export
mix_states <- function(model, pars = "v", alt = c("value", "below", "above")) {
  force(model)
  force(pars)
  alt_mode <- match.arg(alt)
  function(...) {
    m <- model(...)
    if (!identical(m$type, "RACE") || grepl("LogicalRules", m$c_name) ||
        grepl("_CORR|corr", m$c_name))
      stop("mix_states() needs an independent race model.")
    bad <- setdiff(pars, names(m$p_types))
    if (length(bad))
      stop("mix_states(): not parameters of this model: ",
           paste(bad, collapse = ", "))
    suffix <- switch(alt_mode, value = "_mix", below = "_below", above = "_above")
    alt <- paste0(pars, suffix)
    sign <- switch(alt_mode, value = 0, below = -1, above = 1)

    if (sign == 0) {
      m$p_types <- c(m$p_types, stats::setNames(m$p_types[pars], alt),
                     pmix = stats::qnorm(0))
      m$transform <- lapply(m$transform, function(x) {
        if (is.null(names(x)) || !all(pars %in% names(x))) return(x)
        c(x, stats::setNames(x[pars], alt))
      })
      mm <- m$bound$minmax
      have <- pars[pars %in% colnames(mm)]
      if (length(have)) {
        add <- mm[, have, drop = FALSE]
        colnames(add) <- paste0(have, suffix)
        mm <- cbind(mm, add)
      }
      ex <- m$bound$exception
      have <- pars[pars %in% names(ex)]
      if (length(have))
        ex <- c(ex, stats::setNames(ex[have], paste0(have, suffix)))
    } else {
      m$p_types <- c(m$p_types, stats::setNames(rep(log(0), length(alt)), alt),
                     pmix = stats::qnorm(0))
      m$transform$func <- c(m$transform$func,
                            stats::setNames(rep("exp", length(alt)), alt))
      mm <- cbind(m$bound$minmax,
                  matrix(c(0, Inf), 2, length(alt), dimnames = list(NULL, alt)))
      ex <- c(m$bound$exception, stats::setNames(rep(0, length(alt)), alt))
    }
    m$transform$func <- c(m$transform$func, pmix = "pnorm")
    m$bound$minmax <- cbind(mm, pmix = c(0, 1))
    m$bound$exception <- c(ex, pmix = 0)

    # Simulate: draw each accumulator's state, then hand the substituted rows
    # to the parent simulator.  `b` is the one derived column the race
    # Ttransforms add from a mixable parameter, so it is rebuilt here.
    base_rfun <- m$rfun
    m$rfun <- function(data, pars_mat) {
      p <- pars_mat[, "pmix"]
      state <- !is.na(p) & stats::runif(length(p)) < p
      if (any(state)) {
        for (k in seq_along(pars))
          pars_mat[state, pars[k]] <- if (sign == 0) pars_mat[state, alt[k]] else
            pars_mat[state, pars[k]] + sign * pars_mat[state, alt[k]]
        if ("b" %in% colnames(pars_mat) && any(c("B", "A") %in% pars))
          pars_mat[, "b"] <- pars_mat[, "B"] + pars_mat[, "A"]
      }
      base_rfun(data, pars_mat)
    }
    m$state_mixture <- pars
    m
  }
}
