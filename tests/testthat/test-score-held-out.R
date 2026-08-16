# score_held_out(): elpd of data the fit never saw.
#
# The anchor test is an equivalence: scoring the *fitted* data must reproduce,
# exactly, the pointwise log-likelihood that WAIC/LOO already build from the
# posterior (.ll_matrix_pooled).  Same data, same draws, same likelihood -- so
# any disagreement is the held-out dadm being built differently from the fit's,
# which is the whole risk in this code path.

test_that("scoring the fitted data reproduces the in-sample pointwise elpd", {
  e   <- samples_LNR
  dat <- get_data(e)

  # samples_LNR was fitted without rt binning; score at the same resolution.
  sc <- score_held_out(e, dat, type = "trial", rt_resolution = NULL)

  ll  <- EMC2:::.ll_matrix_pooled(e, stage = "sample", filter = 0)
  ref <- apply(ll, 2, function(x) matrixStats::logSumExp(x) - log(length(x)))

  expect_equal(sc$n_points, length(ref))
  expect_equal(unname(sc$pointwise), unname(ref))   # order must match too
  expect_equal(sc$elpd, sum(ref))
  expect_equal(sc$ic, -2 * sc$elpd)
  expect_equal(sc$subject, as.character(dat$subjects[order(dat$subjects)]))
})

test_that("held-out trials from fitted subjects score finitely", {
  e   <- samples_LNR
  dat <- get_data(e)
  test_dat <- do.call(rbind, lapply(split(dat, dat$subjects), function(d)
    d[(floor(nrow(d) / 2) + 1):nrow(d), ]))

  sc <- score_held_out(e, test_dat, type = "trial", rt_resolution = NULL)
  expect_s3_class(sc, "emc.score")
  expect_equal(sc$n_points, nrow(test_dat))
  expect_true(all(is.finite(sc$pointwise)))
  expect_equal(sc$n_subjects, length(unique(test_dat$subjects)))
})

test_that("new subjects are scored marginally and cost more than known ones", {
  e   <- samples_LNR
  dat <- get_data(e)
  subs <- levels(factor(dat$subjects))[1:2]
  test_dat <- dat[dat$subjects %in% subs, ]
  test_dat$subjects <- factor(as.character(test_dat$subjects))

  set.seed(1)
  sc_marg <- score_held_out(e, test_dat, type = "subject", rt_resolution = NULL,
                            K = 25)
  expect_equal(sc_marg$n_points, length(subs))     # one point per subject
  expect_true(all(is.finite(sc_marg$pointwise)))

  # The conditional score knows each subject's alpha, so it must be at least as
  # high as the score that integrates alpha out of the group level.
  sc_cond <- score_held_out(e, test_dat, type = "trial", rt_resolution = NULL)
  cond_by_sub <- tapply(sc_cond$pointwise, sc_cond$subject, sum)
  expect_true(all(sc_marg$pointwise <= cond_by_sub[names(sc_marg$pointwise)]))
})

test_that("unseen subjects are rejected by the conditional path", {
  e   <- samples_LNR
  dat <- get_data(e)
  dat$subjects <- factor(paste0("new_", dat$subjects))
  expect_error(score_held_out(e, dat, type = "trial", rt_resolution = NULL),
               "not in the fit")
})

test_that("malformed held-out data is rejected", {
  e   <- samples_LNR
  dat <- get_data(e)
  expect_error(score_held_out(e, dat[, setdiff(names(dat), "subjects")]),
               "subjects")
  expect_error(score_held_out(e, dat[0, ], rt_resolution = NULL), "no rows")
  expect_error(score_held_out(e, list(dat, dat), rt_resolution = NULL),
               "data frame")
})

test_that("combine_scores pools folds into the whole-data score", {
  e   <- samples_LNR
  dat <- get_data(e)
  whole <- score_held_out(e, dat, type = "trial", rt_resolution = NULL)

  # Split into disjoint folds of the same data and score each separately: the
  # pooled elpd must equal the score of the undivided data, and the pooled SE
  # must be the SE of the whole pointwise vector (not any per-fold SE).
  fold_id <- rep(1:3, length.out = nrow(dat))
  folds   <- split(dat, fold_id)
  scores  <- lapply(folds, function(d) score_held_out(e, d, type = "trial",
                                                      rt_resolution = NULL))
  cv <- combine_scores(scores)

  expect_s3_class(cv, "emc.score")
  expect_equal(cv$n_folds, 3L)
  expect_equal(cv$n_points, whole$n_points)
  expect_equal(cv$elpd, whole$elpd)
  expect_equal(cv$ic, whole$ic)
  expect_equal(sort(unname(cv$pointwise)), sort(unname(whole$pointwise)))
  expect_equal(cv$se_elpd, whole$se_elpd)
  expect_equal(sum(cv$folds$elpd), cv$elpd)
  # a combined score can be combined again
  expect_equal(combine_scores(list(cv))$elpd, cv$elpd)
})

test_that("combine_scores refuses overlapping folds and mixed types", {
  e   <- samples_LNR
  dat <- get_data(e)
  s1 <- score_held_out(e, dat[1:100, ], type = "trial", rt_resolution = NULL)
  s2 <- score_held_out(e, dat[50:150, ], type = "trial", rt_resolution = NULL)
  expect_error(combine_scores(list(s1, s2)), "overlap")
  expect_equal(combine_scores(list(s1, s2), check_disjoint = FALSE)$n_points,
               s1$n_points + s2$n_points)

  subs <- levels(factor(dat$subjects))
  m1 <- score_held_out(e, dat[dat$subjects == subs[1], ], type = "subject",
                       rt_resolution = NULL, K = 5)
  expect_error(combine_scores(list(s1, m1)), "different types")

  expect_error(combine_scores(list(s1, "not a score")), "not emc.score")
  expect_error(combine_scores(list()), "non-empty")
})

test_that("combine_scores pools held-out subjects across folds", {
  e    <- samples_LNR
  dat  <- get_data(e)
  subs <- levels(factor(dat$subjects))
  set.seed(2)
  scores <- lapply(subs, function(s) {
    d <- dat[dat$subjects == s, ]
    d$subjects <- factor(as.character(d$subjects))
    score_held_out(e, d, type = "subject", rt_resolution = NULL, K = 10)
  })
  cv <- combine_scores(scores)
  expect_equal(cv$n_points, length(subs))       # one point per held-out subject
  expect_equal(cv$n_subjects, length(subs))
  expect_equal(cv$elpd, sum(vapply(scores, function(x) x$elpd, numeric(1))))
  # the same subject held out in two folds is double counting
  expect_error(combine_scores(c(scores, scores[1])), "overlap")
})

test_that("score_folds scores each fold against its own fit and pools", {
  e   <- samples_LNR
  dat <- get_data(e)
  fold_id <- rep(1:2, length.out = nrow(dat))
  folds   <- split(dat, fold_id)

  # Same emc for both folds here (this is a plumbing test, not a CV study);
  # in real use each fold has its own fit that excluded that fold's data.
  cv <- score_folds(list(e, e), folds, type = "trial", rt_resolution = NULL)
  ref <- combine_scores(lapply(folds, function(d)
    score_held_out(e, d, type = "trial", rt_resolution = NULL)))
  expect_equal(cv$elpd, ref$elpd)
  expect_equal(cv$n_points, nrow(dat))

  expect_error(score_folds(list(e, e), folds[1]), "same length")
  expect_error(score_folds(list(e), list(dat[, setdiff(names(dat), "subjects")])),
               "Fold 1")
})

test_that("an rt_resolution that disagrees with the fit warns", {
  e   <- samples_LNR
  dat <- get_data(e)
  # samples_LNR's rts are unbinned, so scoring at 1/60 is a mismatch.
  expect_warning(score_held_out(e, dat[dat$subjects == levels(dat$subjects)[1], ],
                                type = "trial", rt_resolution = 1/60),
                 "not on that grid")
})
