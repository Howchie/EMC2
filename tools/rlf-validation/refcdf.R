# Shared reference cdf.  Do NOT build this from a single rlf_fht_pdf_cdf_vec
# call at large t_max: rlf_lower_extent sizes the domain from the horizon, so a
# 22 s grid at nx = 768 leaves ~12 cells between the start point and the
# boundary -- coarser at the boundary than the nx = 128 default is at 2 s.  Go
# through pRLF instead, which buckets each time by its own horizon and so keeps
# every point on a domain sized for it.
ref_cdf <- function(a, v, tg, B, A, nx = 384L) {
  o <- options(emc2.rlf_nx = nx, emc2.rlf_richardson = TRUE,
               emc2.rlf_richardson_ratio = 1.5, emc2.rlf_horizon_split = TRUE,
               emc2.rlf_dt = 4e-3, emc2.rlf_sparse_output = TRUE)
  on.exit(options(o))
  pars <- cbind(v = rep(v, length(tg)), B = B, A = A, t0 = 0, s = 1, alpha = a)
  C <- EMC2:::pRLF(tg, pars)
  cummax(C)   # bucket boundaries can leave a sub-1e-4 dip; keep it monotone
}
