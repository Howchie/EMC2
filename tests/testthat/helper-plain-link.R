# A model whose exp/pnorm sides are all declared at the link's default (lower
# 0, upper 1), so design() does not move them to the bound. For tests that
# write natural values on the plain sampled scale (t0 = log(0.25) meaning
# t0 = 0.25) and compare against natural-scale references.
plain_link <- function(model) {
  force(model)
  function(...) {
    ml <- if (is.function(model)) model(...) else model
    tf <- ml$transform
    has_lo <- names(tf$func)[tf$func %in% c("exp", "pnorm")]
    has_up <- names(tf$func)[tf$func == "pnorm"]
    lo <- stats::setNames(rep(0, length(has_lo)), has_lo)
    up <- stats::setNames(rep(1, length(has_up)), has_up)
    lo[names(tf$lower)] <- tf$lower
    up[names(tf$upper)] <- tf$upper
    ml$transform$lower <- lo
    ml$transform$upper <- up
    ml
  }
}
