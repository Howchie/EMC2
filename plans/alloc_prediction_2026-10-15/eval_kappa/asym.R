SP <- commandArgs(TRUE)[1]
d <- readRDS(file.path(SP, "eval.rds")); P <- readRDS(file.path(SP, "pairs.rds"))
get_F <- function(plans) do.call(rbind, lapply(seq_along(plans), function(j) do.call(rbind, lapply(plans[[j]], function(h)
  if (!is.null(h$pred_group)) data.frame(ch = j, R = h$iteration, coord = seq_along(h$pred_group), rounds = h$rounds)))))
Fs <- rbind(cbind(set = "PM burn", get_F(P$pm_burn$plans)), do.call(rbind, lapply(names(P$arms), function(nm) cbind(set = paste("arm", nm), get_F(P$arms[[nm]]$plans)))))
Fo <- Fs; Fo$R <- Fo$R + 100; names(Fo)[names(Fo) == "rounds"] <- "rounds_old"
x <- merge(merge(d, Fs), Fo); x$dn <- x$rounds - x$rounds_old
x$slowq <- x$tau_meas >= quantile(x$tau_meas, 0.75)
for (m in c("tB", "tK50m", "tA")) {
  z <- log(x[[m]] / x$tR)
  f <- function(ok) sprintf("%+.2f (n %d, se %.2f)", median(z[ok]), sum(ok), 1.25 * sd(z[ok]) / sqrt(sum(ok)))
  cat(sprintf("%-6s median log(pred/real): rounds cut >0.1 %s | about same %s | rounds added >0.1 %s | slow quartile, cut: %s\n", m,
    f(x$dn < -0.1), f(abs(x$dn) <= 0.1), f(x$dn > 0.1), f(x$dn < -0.1 & x$slowq)))
}
