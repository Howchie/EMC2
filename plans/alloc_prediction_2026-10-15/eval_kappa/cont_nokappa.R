# cont.R with kappa held at one (the point-estimate allocator alone).
suppressPackageStartupMessages(library(EMC2, lib.loc = Sys.getenv("EMC_LIB")))
assignInNamespace(".emc_alloc_kappa", function(X) rep(1, nrow(X[[1L]])), "EMC2")
source(file.path(Sys.getenv("SP"), "cont.R"))
