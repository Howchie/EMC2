# Per group coordinate at a saved state: rest-of-iteration rho, the plan's predicted lag-1 rho, and the whole-iteration lag-1 / lag-25 rho of the window's stored draws (pooled window variance); then the same on the NEXT block if present.
suppressMessages(library(EMC2, lib.loc = Sys.getenv("LIB"))); ns <- asNamespace("EMC2")
emc <- readRDS(commandArgs(TRUE)[1]); j <- 1
for (i in 2:length(emc)) for (nm in setdiff(names(emc[[1]]), names(emc[[i]]))) emc[[i]][[nm]] <- emc[[1]][[nm]]
p <- emc[[1]]$n_pars; pn <- emc[[1]]$par_names
ch <- emc[[j]]; w <- ch$alloc$window; idx <- ch$samples$idx; Dinv <- solve(ch$shear$plan$dirs)
T_g <- rbind(cbind(Dinv, matrix(0, p, p)), cbind(matrix(0, p, p), diag(p)))
coord <- function(x, k) rbind(Dinv %*% x$samples$theta_mu[, k], ns$.emc_alloc_lsd_draws(x$samples$theta_var, k))
draws <- (idx - w$n + 1):idx
X <- lapply(emc, coord, k = draws); V <- apply(do.call(cbind, X), 1, var)
jump <- function(L) rowMeans(do.call(cbind, lapply(X, function(x) (x[, -(1:L), drop = FALSE] - x[, 1:(ncol(x) - L), drop = FALSE])^2))) / L
rho_at <- function(t, L) pmax(0, 1 - t * L / (2 * V))^(1 / L); tau1 <- function(r) (1 + pmin(r, 1 - 1e-9)) / (1 - pmin(r, 1 - 1e-9))
rest <- diag(T_g %*% w$G %*% t(T_g)) / w$nG
last <- ch$alloc$last
d <- data.frame(coord = c(paste0("dir", 1:p), paste0("lsd.", pn)), rho_rest = 1 - rest / (2 * V), g = last$group_g, a = last$group_a, rounds = round(last$rounds, 2),
  pred = last$pred_group, J1 = rho_at(jump(1), 1), J10 = rho_at(jump(10), 10), J25 = rho_at(jump(25), 25))
d$tau_pred <- tau1(d$pred); d$tau_J1 <- tau1(d$J1); d$tau_J25 <- tau1(d$J25)
print(format(d, digits = 3), row.names = FALSE)
