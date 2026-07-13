.oo_model_list <- function(model) {
  if (is.function(model)) {
    return(model())
  }
  model
}

.oo_expanded_designs <- function(dadm, row_idx = NULL, expand = TRUE) {
  designs <- attr(dadm, "designs")
  if (is.null(designs)) {
    stop("dadm must have a 'designs' attribute")
  }

  out <- lapply(designs, function(design_mat) {
    expand_idx <- attr(design_mat, "expand")
    if (!expand) {
      if (!is.null(row_idx)) {
        if (is.null(expand_idx)) {
          return(design_mat[row_idx, , drop = FALSE])
        }
        # Keep the compressed matrix and restrict only its row map.  The C++
        # mapper consumes this attribute directly, avoiding a large R matrix
        # allocation for every posterior-predictive draw.
        out_mat <- design_mat
        attr(out_mat, "expand") <- expand_idx[row_idx]
        return(out_mat)
      }
      return(design_mat)
    }
    if (is.null(expand_idx)) {
      expand_idx <- seq_len(nrow(design_mat))
    }

    expanded <- design_mat[expand_idx, , drop = FALSE]
    if (!is.null(row_idx)) {
      expanded <- expanded[row_idx, , drop = FALSE]
    }
    expanded
  })

  names(out) <- names(designs)
  out
}

.oo_identity_transform <- function(param_names) {
  list(
    func = stats::setNames(rep("identity", length(param_names)), param_names),
    lower = stats::setNames(rep(-Inf, length(param_names)), param_names),
    upper = stats::setNames(rep(Inf, length(param_names)), param_names)
  )
}

.oo_reorder_public_pars <- function(pars, model_list) {
  base_names <- names(model_list$p_types)
  ord <- c(intersect(base_names, colnames(pars)), setdiff(colnames(pars), base_names))
  if (identical(ord, colnames(pars))) {
    return(pars)
  }

  ok <- attr(pars, "ok")
  pars <- pars[, ord, drop = FALSE]
  attr(pars, "ok") <- ok
  pars
}

.oo_particle_matrix <- function(p, dadm, keep_all_columns = FALSE) {
  sampled_p_names <- attr(dadm, "sampled_p_names")
  p_names <- attr(dadm, "p_names")

  if (is.null(sampled_p_names)) {
    stop("dadm must have a 'sampled_p_names' attribute")
  }

  if (is.data.frame(p)) {
    p <- as.matrix(p)
  }

  target_names <- if (keep_all_columns && !is.null(p_names)) {
    p_names
  } else {
    sampled_p_names
  }

  # Data-aware sampled_pars() deliberately drops coefficients whose design
  # column is identically zero.  Such a column is not a missing parameter in
  # the likelihood: its contribution is exactly zero.  Constants are a
  # second, distinct case; restore their recorded value rather than silently
  # replacing (for example) a fixed log-SD with zero.
  constants <- attr(dadm, "constants")
  designs <- attr(dadm, "designs")
  is_zero_design <- function(nm) {
    if (is.null(designs)) return(FALSE)
    any(vapply(designs, function(design_mat) {
      if (is.null(colnames(design_mat)) || !(nm %in% colnames(design_mat))) {
        return(FALSE)
      }
      values <- design_mat[, nm]
      length(values) > 0L && all(is.finite(values) & values == 0)
    }, logical(1)))
  }

  missing_columns <- function(p_matrix) {
    missing <- setdiff(target_names, colnames(p_matrix))
    if (!length(missing)) {
      return(p_matrix[, target_names, drop = FALSE])
    }

    additions <- matrix(NA_real_, nrow = nrow(p_matrix), ncol = length(missing),
                        dimnames = list(rownames(p_matrix), missing))
    unresolved <- character()
    for (i in seq_along(missing)) {
      nm <- missing[[i]]
      if (!is.null(constants) && !is.null(names(constants)) && nm %in% names(constants) &&
          length(constants[[nm]]) == 1L && !is.na(constants[[nm]])) {
        additions[, i] <- constants[[nm]]
      } else if (is_zero_design(nm)) {
        additions[, i] <- 0
      } else {
        unresolved <- c(unresolved, nm)
      }
    }

    if (length(unresolved)) {
      stop("p matrix columns must include: ", paste(unresolved, collapse = ", "))
    }
    cbind(p_matrix, additions)[, target_names, drop = FALSE]
  }

  if (is.null(dim(p))) {
    if (is.null(names(p))) {
      if (length(p) == length(target_names)) {
        names(p) <- target_names
      } else if (!is.null(p_names) && length(p) == length(sampled_p_names)) {
        names(p) <- sampled_p_names
      } else if (!is.null(p_names) && length(p) == length(p_names)) {
        names(p) <- p_names
      } else {
        stop("Unnamed p vector must have length ", length(target_names),
             " (or ", length(sampled_p_names), " sampled columns)")
      }
    }

    out <- matrix(p, nrow = 1)
    colnames(out) <- names(p)
    return(missing_columns(out))
  }

  if (is.null(colnames(p))) {
    if (ncol(p) == length(target_names)) {
      colnames(p) <- target_names
    } else if (ncol(p) == length(sampled_p_names)) {
      colnames(p) <- sampled_p_names
    } else if (!is.null(p_names) && ncol(p) == length(p_names)) {
      colnames(p) <- p_names
    } else {
      stop("p matrix must have columns for: ", paste(target_names, collapse = ", "))
    }
  }

  missing_columns(p)
}

get_pars_oo <- function(p, dadm, model,
                        pretransformed = FALSE,
                        constants_included = FALSE,
                        return_kernel_matrix = FALSE,
                        return_all_pars = FALSE,
                        kernel_output_codes = 1L) {
  model_list <- .oo_model_list(model)
  particle_matrix <- .oo_particle_matrix(p, dadm, keep_all_columns = constants_included)
  constants <- attr(dadm, "constants")
  if (constants_included || is.null(constants)) {
    constants <- NA
  }
  pretransforms <- if (pretransformed) {
    .oo_identity_transform(colnames(particle_matrix))
  } else {
    model_list$pre_transform
  }

  call_one <- function(cur_particles, cur_dadm, row_idx = NULL) {
    get_pars_c_wrapper_oo(
      particle_matrix = cur_particles,
      data = cur_dadm,
      constants = constants,
      designs = .oo_expanded_designs(dadm, row_idx, expand = FALSE),
      bounds = model_list$bound,
      transforms = model_list$transform,
      pretransforms = pretransforms,
      trend = model_list$trend,
      return_kernel_matrix = return_kernel_matrix,
      return_all_pars = return_all_pars,
      kernel_output_codes = kernel_output_codes
    )
  }

  if (!("subjects" %in% names(dadm)) || nrow(particle_matrix) == 1L) {
    return(call_one(particle_matrix, dadm))
  }

  if (is.null(rownames(particle_matrix))) {
    stop("Multi-row p matrix must have rownames for every subject in dadm")
  }

  subj_levels <- levels(dadm$subjects)
  used_subjects <- subj_levels[subj_levels %in% as.character(dadm$subjects)]

  if (!all(used_subjects %in% rownames(particle_matrix))) {
    stop("p matrix must have rows named for every subject in dadm")
  }

  particle_matrix <- particle_matrix[used_subjects, , drop = FALSE]

  pieces <- vector("list", length(used_subjects))
  row_ids <- vector("list", length(used_subjects))

  for (i in seq_along(used_subjects)) {
    row_idx <- dadm$subjects == used_subjects[i]
    row_ids[[i]] <- which(row_idx)
    pieces[[i]] <- call_one(
      cur_particles = particle_matrix[i, , drop = FALSE],
      cur_dadm = dadm[row_idx, , drop = FALSE],
      row_idx = row_idx
    )
  }

  out <- matrix(NA_real_, nrow = nrow(dadm), ncol = ncol(pieces[[1]]))
  colnames(out) <- colnames(pieces[[1]])
  for (i in seq_along(pieces)) {
    out[row_ids[[i]], ] <- pieces[[i]]
  }
  out
}

# Batched counterpart to get_pars_oo().  It deliberately maps one subject's
# draws at a time: the model mapping is subject-specific, while the expensive
# design/transform work for all draws stays within one C++ call.
get_pars_batch_oo <- function(p, dadm, model, row_idx = NULL,
                              pretransformed = FALSE,
                              constants_included = FALSE,
                              return_kernel_matrix = FALSE,
                              return_all_pars = FALSE,
                              kernel_output_codes = 1L) {
  model_list <- .oo_model_list(model)
  if (is.null(row_idx)) row_idx <- seq_len(nrow(dadm))
  cur_dadm <- dadm[row_idx, , drop = FALSE]
  particle_matrix <- .oo_particle_matrix(p, dadm,
                                         keep_all_columns = constants_included)
  constants <- attr(dadm, "constants")
  if (constants_included || is.null(constants)) constants <- NA
  pretransforms <- if (pretransformed) {
    .oo_identity_transform(colnames(particle_matrix))
  } else {
    model_list$pre_transform
  }

  get_pars_c_batch_wrapper_oo(
    particle_matrix = particle_matrix,
    data = cur_dadm,
    constants = constants,
    designs = .oo_expanded_designs(dadm, row_idx, expand = FALSE),
    bounds = model_list$bound,
    transforms = model_list$transform,
    pretransforms = pretransforms,
    trend = model_list$trend,
    return_kernel_matrix = return_kernel_matrix,
    return_all_pars = return_all_pars,
    kernel_output_codes = kernel_output_codes
  )
}

get_pars_matrix_oo <- function(p_vector, dadm, model) {
  model_list <- .oo_model_list(model)
  pars <- get_pars_oo(p_vector, dadm, model_list)
  pars <- model_list$Ttransform(pars, dadm)
  pars <- add_bound(pars, model_list$bound, dadm$lR)
  .oo_reorder_public_pars(pars, model_list)
}
