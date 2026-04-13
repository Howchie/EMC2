if (requireNamespace("EMC2", quietly = TRUE)) {
  library(EMC2)
} else {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("EMC2 is not installed and pkgload is unavailable to load local sources.")
  }
  pkgload::load_all(".", quiet = TRUE)
  library(EMC2)
}
library(microbenchmark)

ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
matchfun <- function(d) d$S == d$lR
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2], ]
dat$subjects <- droplevels(dat$subjects)

make_ctx <- function(model, formula) {
  dsg <- design(
    data = dat, model = model, matchfun = matchfun,
    formula = formula,
    contrasts = list(v = list(lM = ADmat))
  )
  emc <- make_emc(dat, dsg, compress = FALSE)
  m <- emc[[1]]$model()
  dadm <- emc[[1]]$data[[1]]
  p_types <- names(m$p_types)
  p_sampled <- sampled_pars(emc)
  designs <- list()
  for (p in p_types) {
    dsgp <- attr(dadm, "designs")[[p]]
    designs[[p]] <- dsgp[attr(dsgp, "expand"), , drop = FALSE]
  }
  constants <- attr(dadm, "constants")
  if (is.null(constants)) constants <- NA
  list(model = m, dadm = dadm, p_types = p_types, p_sampled = p_sampled, designs = designs, constants = constants)
}

ll_oo <- function(p_mat, ctx) {
  EMC2:::calc_ll_oo(
    p_mat, ctx$dadm,
    constants = ctx$constants,
    designs = ctx$designs,
    type = ctx$model$c_name,
    bounds = ctx$model$bound,
    transforms = ctx$model$transform,
    pretransforms = ctx$model$pre_transform,
    p_types = ctx$p_types,
    min_ll = log(1e-10),
    trend = ctx$model$trend
  )
}

ctx_rdm <- make_ctx(RDM, list(v ~ lM, s ~ 1, t0 ~ 1, A ~ 1, B ~ 1))
ctx_swtn <- make_ctx(RDMSWTN, list(v ~ lM, s ~ 1, t0 ~ 1, A ~ 1, B ~ 1, sv ~ 1))
ctx_lba <- make_ctx(LBA, list(v ~ lM, sv ~ 1, t0 ~ 1, A ~ 1, B ~ 1))
ctx_bawl <- make_ctx(BAwL, list(v ~ lM, sv ~ 1, t0 ~ 1, A ~ 1, B ~ 1, k ~ 1))

# Single-particle reference params in sampled-parameter space
p_rdm <- ctx_rdm$p_sampled
p_rdm[] <- 0
p_rdm[c("v", "v_lMd", "s", "A", "B", "t0")] <- c(1, 0.5, 0, -1, 0, -1.2)
p_rdm["pContaminant"] <- -5

p_swtn <- ctx_swtn$p_sampled
p_swtn[] <- 0
for (nm in names(p_rdm)) if (nm %in% names(p_swtn)) p_swtn[nm] <- p_rdm[nm]
p_swtn["sv"] <- -Inf

p_lba <- ctx_lba$p_sampled
p_lba[] <- 0
p_lba[c("v", "v_lMd", "sv", "A", "B", "t0")] <- c(1, 0.5, 0, -1, 0, -1.2)
p_lba["pContaminant"] <- -5

p_bawl <- ctx_bawl$p_sampled
p_bawl[] <- 0
for (nm in names(p_lba)) if (nm %in% names(p_bawl)) p_bawl[nm] <- p_lba[nm]
p_bawl["k"] <- -Inf

cat("Numerical check RDM vs RDMSWTN (sv=0):\n")
ll_rdm <- ll_oo(matrix(p_rdm, nrow = 1, dimnames = list(NULL, names(p_rdm))), ctx_rdm)
ll_swtn <- ll_oo(matrix(p_swtn, nrow = 1, dimnames = list(NULL, names(p_swtn))), ctx_swtn)
cat("RDM LL:", ll_rdm, "\n")
cat("SWTN LL:", ll_swtn, "\n")
cat("Difference:", ll_rdm - ll_swtn, "\n\n")

cat("Numerical check LBA vs BAwL (k=0):\n")
ll_lba <- ll_oo(matrix(p_lba, nrow = 1, dimnames = list(NULL, names(p_lba))), ctx_lba)
ll_bawl <- ll_oo(matrix(p_bawl, nrow = 1, dimnames = list(NULL, names(p_bawl))), ctx_bawl)
cat("LBA LL:", ll_lba, "\n")
cat("BAwL LL:", ll_bawl, "\n")
cat("Difference:", ll_lba - ll_bawl, "\n\n")

mk_pmat <- function(p, n = 100L) {
  # Repeat the full parameter vector per row (not each element in blocks),
  # otherwise columns get scrambled across particles.
  out <- matrix(rep(p, times = n), nrow = n, byrow = TRUE)
  colnames(out) <- names(p)
  out
}

pmat_rdm <- mk_pmat(p_rdm)
pmat_swtn <- mk_pmat(p_swtn)
pmat_lba <- mk_pmat(p_lba)
pmat_bawl <- mk_pmat(p_bawl)

cat("Benchmarking calc_ll_oo (100 particles)...\n")
res <- microbenchmark(
  RDM = ll_oo(pmat_rdm, ctx_rdm),
  RDMSWTN_sv0 = ll_oo(pmat_swtn, ctx_swtn),
  LBA = ll_oo(pmat_lba, ctx_lba),
  BAwL_k0 = ll_oo(pmat_bawl, ctx_bawl),
  times = 50
)
print(res)

pmat_swtn_act <- pmat_swtn
pmat_swtn_act[, "sv"] <- log(0.1)
pmat_bawl_act <- pmat_bawl
pmat_bawl_act[, "k"] <- log(0.1)

cat("\nBenchmarking active variability (sv=0.1, k=0.1)...\n")
res_act <- microbenchmark(
  RDMSWTN_sv_active = ll_oo(pmat_swtn_act, ctx_swtn),
  BAwL_k_active = ll_oo(pmat_bawl_act, ctx_bawl),
  times = 50
)
print(res_act)
