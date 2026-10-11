pkgload::load_all(".", quiet = TRUE, compile = FALSE)
z <- readRDS("/tmp/collapsed_transport_followup/setup_saturation_start.rds")
for (multiplier in c(0.25, 0.5, 1, 2, 4)) {
  rates <- vapply(z$captured, function(item) {
    u <- item$units
    u$joint$cost <- u$joint$cost * multiplier
    u$joint$fixed_cost <- u$joint$fixed_cost * multiplier
    EMC2:::.emc_alloc_solve_units(u)$joint
  }, numeric(1))
  cat("cost multiplier", multiplier, "rates", paste(rates, collapse = ","), "\n")
}
