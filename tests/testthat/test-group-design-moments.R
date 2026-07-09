# group_design_moments assembles the data-based precision/mean for the joint
# group-level Gibbs update blockwise. Pin it to the direct per-subject
# construction (the pre-optimization implementation) on random inputs.

test_that("group_design_moments matches the per-subject M_i loop", {
  loop_moments <- function(group_designs, tvinv, alpha, M) {
    p <- nrow(alpha); n <- ncol(alpha)
    prec_data <- matrix(0, M, M)
    mean_data <- numeric(M)
    for (i in seq_len(n)) {
      par_idx <- 0
      M_i <- matrix(0, nrow = p, ncol = M)
      for (k in seq_len(p)) {
        M_i[k, par_idx + 1:ncol(group_designs[[k]])] <- group_designs[[k]][i, , drop = FALSE]
        par_idx <- par_idx + ncol(group_designs[[k]])
      }
      prec_data <- prec_data + crossprod(M_i, tvinv %*% M_i)
      mean_data <- mean_data + crossprod(M_i, tvinv %*% alpha[, i, drop = FALSE])
    }
    list(prec_data = prec_data, mean_data = as.numeric(mean_data))
  }

  for (seed in 1:3) {
    set.seed(seed)
    p <- 7; n <- 23
    group_designs <- lapply(seq_len(p), function(k) {
      m_k <- sample(1:3, 1)
      if (m_k == 1) matrix(1, n, 1) else cbind(1, matrix(rnorm(n * (m_k - 1)), n))
    })
    A <- matrix(rnorm(p * p), p)
    tvinv <- crossprod(A) + diag(p)
    alpha <- matrix(rnorm(p * n), p, n)
    M <- sum(vapply(group_designs, ncol, integer(1)))

    ref <- loop_moments(group_designs, tvinv, alpha, M)
    new <- EMC2:::group_design_moments(group_designs, tvinv, alpha, M)
    expect_equal(new$prec_data, ref$prec_data, tolerance = 1e-10)
    expect_equal(new$mean_data, ref$mean_data, tolerance = 1e-10)
  }
})
