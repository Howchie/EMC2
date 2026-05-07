require(dplyr)
## per-subject model-vs-data accuracy and RT quantiles (correct/incorrect), by sa.
quant_long <- function(df, value_name) {
  bind_rows(
    df %>% transmute(subjects, sa, metric, qlabel = "Q10", !!value_name := q10),
    df %>% transmute(subjects, sa, metric, qlabel = "Q50", !!value_name := q50),
    df %>% transmute(subjects, sa, metric, qlabel = "Q90", !!value_name := q90)
  )
}

plot_subj_estimates <- function(data,post_predict,plot_out=TRUE,factor = c("S"),Cfun=NULL) {
  if(!("Correct"%in%names(data)&"Correct"%in%names(post_predict))) {
    if(is.null(Cfun)) Cfun <- function(d) as.numeric(d$S)==as.numeric(d$R)
    data = data %>% mutate(Correct=Cfun(.))
    post_predict = post_predict %>% mutate(Correct=Cfun(.))
  }
  factors = paste(c(factor,"subjects"),collapse="|")
  cols <- grep(factors, colnames(data), value = TRUE)
  obs_df <- data
  mod_df <- post_predict
  
  obs_acc <- obs_df %>%
    group_by(across(all_of(cols))) %>%
    summarise(metric = "Accuracy", qlabel = "Mean", data_value = mean(Correct), .groups = "drop")
  obs_q <- obs_df %>%
    group_by(across(all_of(cols)),Correct) %>% # extend to factors
    summarise(
      q10 = as.numeric(quantile(rt, probs = 0.10, na.rm = TRUE, type = 7)),
      q50 = as.numeric(quantile(rt, probs = 0.50, na.rm = TRUE, type = 7)),
      q90 = as.numeric(quantile(rt, probs = 0.90, na.rm = TRUE, type = 7)),
      .groups = "drop"
    ) %>%
    mutate(metric = ifelse(Correct, "Correct RT", "Error RT")) %>%
    select(subjects,all_of(factor), metric, q10, q50, q90)
  obs_rt <- quant_long(obs_q, "data_value")
  obs_stats <- bind_rows(obs_acc, obs_rt)
  
  mod_acc <- mod_df %>%
    group_by(across(all_of(cols))) %>%
    summarise(metric = "Accuracy", qlabel = "Mean", model_value = mean(Correct), .groups = "drop")
  mod_q <- mod_df %>%
    group_by(across(all_of(cols)),Correct,postn) %>%
    summarise(
      q10 = as.numeric(quantile(rt, probs = 0.10, na.rm = TRUE, type = 7)),
      q50 = as.numeric(quantile(rt, probs = 0.50, na.rm = TRUE, type = 7)),
      q90 = as.numeric(quantile(rt, probs = 0.90, na.rm = TRUE, type = 7)),
      .groups = "drop"
    ) %>%
    group_by(group_by(across(all_of(cols))),Correct) %>%
    summarise(
      q10 = mean(q10, na.rm = TRUE),
      q50 = mean(q50, na.rm = TRUE),
      q90 = mean(q90, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(metric = ifelse(Correct, "Correct RT", "Error RT")) %>%
    select(subjects,all_of(factor), metric, q10, q50, q90)
  mod_rt <- quant_long(mod_q, "model_value")
  mod_stats <- bind_rows(mod_acc, mod_rt)
  
  df <- inner_join(obs_stats, mod_stats, by = c(cols, "metric", "qlabel"))
  
  if (nrow(df) > 0 & plot_out) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    par(mfrow = c(1, 3), mar = c(4, 4, 2.5, 1))
    fac_cols <- viridis::cividis(length(levels(df[[factor]])))
    names(fac_cols) = levels(df[[factor]])
    q_pch <- c("Q10" = 1, "Q50" = 16, "Q90" = 2, "Mean" = 16)
    
    for (m in c("Accuracy", "Correct RT", "Error RT")) {
      d <- df %>% filter(metric == m)
      if (nrow(d) == 0) {
        plot.new()
        title(main = m)
        next
      }
      x <- d$data_value
      y <- d$model_value
      lim <- range(c(x, y), finite = TRUE)
      plot(x, y, type = "n",
           xlab = "Data", ylab = "Model", main = m,
           xlim = lim, ylim = lim)
      abline(0, 1, lty = 2, col = "gray40")
      rmse <- sqrt(mean((y - x)^2, na.rm = TRUE))
      usr <- par("usr")
      text(usr[1], usr[4], labels = sprintf("RMSE = %.4f", rmse),
           adj = c(0, 1), xpd = NA)
      cols <- fac_cols[as.character(d[[factor]])]
      pchv <- q_pch[as.character(d$qlabel)]
      points(x, y, col = cols, pch = pchv, cex = 0.9)
    }
    legend("bottomright", legend = names(fac_cols), col = fac_cols, pch = 16, bty = "n", cex = 0.9)
  }
  invisible(list(plot_df=df))
}

# Random Wald generation function, borrowed with permission from:
# Heathcote, A. (2004). Fitting Wald and ex-Wald distributions to response time data: An example using functions for the S-PLUS package. Behavior Research Methods, Instruments, & Computers, 36(4), 678-694.
rwald <- function(n,m,a,s=0,s_random=FALSE) {
  if(length(n)>1) n <- length(n);
  if(length(m)>1 && length(m)!=n) m <- rep(m,length=n)
  if(length(a)>1 && length(a)!=n) lambda <- rep(a,length=n)
  if(s_random) s=runif(n,0,s)
  y2 <- rchisq(n,1); y2onm <- y2/m; u <- runif(n)
  r1 <- (2*a + y2onm - sqrt(y2onm*(4*a+y2onm)))/(2*m)
  r2 <- (a/m)^2/r1; ifelse(u < a/(a+m*r1), s+r1, s+r2)
}
