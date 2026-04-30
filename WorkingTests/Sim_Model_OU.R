.ou_init_state <- function(p_vec) {
  if ("uc_yes" %in% names(p_vec) && !("uc_H" %in% names(p_vec))) {
    p_vec["uc_H"] <- .param_get(p_vec, "uc_yes")
  }
  if (!("uc_L" %in% names(p_vec))) {
    p_vec["uc_L"] <- .param_get(p_vec, "uc_H")
  }
  list(
    p_vec = p_vec,
    cap_target = if ("capacity_target" %in% names(p_vec)) .param_get(p_vec, "capacity_target") else 0,
    cap_absence = if ("capacity_absence" %in% names(p_vec)) .param_get(p_vec, "capacity_absence") else 0
  )
}

.ou_cell_pars <- function(p_vec, meta_row, state) {
  pv <- state$p_vec
  pars <- matrix(NA_real_, nrow = 4, ncol = 5)
  rownames(pars) <- c("A", "n_A", "B", "n_B")
  colnames(pars) <- c("a", "u", "sv", "b", "t0")

  pars[, "b"] <- .param_get(pv, "b")
  pars[, "t0"] <- .param_get(pv, "t0")
  pars[, "a"] <- -1
  pars[, "sv"] <- .param_get(pv, "sv")

  l1 <- as.integer(meta_row$Channel1)
  l2 <- as.integer(meta_row$Channel2)

  if (l1 == 2L) {
    pars["A", "u"] <- .param_get(pv, "uc_H")
    pars["n_A", "u"] <- .param_get(pv, "ue_no")
  } else if (l1 == 1L) {
    pars["A", "u"] <- .param_get(pv, "uc_L")
    pars["n_A", "u"] <- .param_get(pv, "ue_no")
  } else {
    pars["A", "u"] <- .param_get(pv, "ue_yes")
    pars["n_A", "u"] <- .param_get(pv, "uc_no")
  }

  if (l2 == 2L) {
    pars["B", "u"] <- .param_get(pv, "uc_H")
    pars["n_B", "u"] <- .param_get(pv, "ue_no")
  } else if (l2 == 1L) {
    pars["B", "u"] <- .param_get(pv, "uc_L")
    pars["n_B", "u"] <- .param_get(pv, "ue_no")
  } else {
    pars["B", "u"] <- .param_get(pv, "ue_yes")
    pars["n_B", "u"] <- .param_get(pv, "uc_no")
  }

  is_target_pair <- (l1 > 0L) && (l2 > 0L)
  is_absence_pair <- (l1 == 0L) && (l2 == 0L)

  u_base <- pars[, "u"]
  B_in <- diag(4)
  rownames(B_in) <- rownames(pars)
  colnames(B_in) <- rownames(pars)

  if (is_target_pair) {
    B_in[c("A", "B"), c("A", "B")] <- matrix(c(1, state$cap_target, state$cap_target, 1), nrow = 2, byrow = TRUE)
  } else if (is_absence_pair) {
    B_in[c("n_A", "n_B"), c("n_A", "n_B")] <- matrix(c(1, state$cap_absence, state$cap_absence, 1), nrow = 2, byrow = TRUE)
  }

  pars[, "u"] <- as.numeric(B_in %*% u_base)
  pars
}

.ou_simulate_rule_cell <- function(n, pars, logical_rule, p_vec,
                                   dt = 0.001, max_time = 2,
                                   bridge = TRUE, bridge_exp_cutoff = -20,
                                   backend = c("auto", "Rcpp", "R"),
                                   adaptive = FALSE, adaptive_mode = c("distance", "ptol"),
                                   adapt_factor = 32.0, adapt_eps = 1.5,
                                   p_tol = 1e-9, curv_tol = 0.1) {
  LogicalRules_OU_rfun(
    n = n,
    pars = pars,
    LogicalRule = logical_rule,
    cross_talk = c(
      pos = if ("cross_talk_pos" %in% names(p_vec)) .param_get(p_vec, "cross_talk_pos") else 0,
      neg = if ("cross_talk_neg" %in% names(p_vec)) .param_get(p_vec, "cross_talk_neg") else 0
    ),
    dt = dt,
    max_time = max_time,
    bridge = bridge,
    bridge_exp_cutoff = bridge_exp_cutoff,
    backend = backend,
    adaptive = adaptive,
    adaptive_mode = adaptive_mode,
    adapt_factor = adapt_factor,
    adapt_eps = adapt_eps,
    p_tol = p_tol,
    curv_tol = curv_tol
  )
}
LogicalRules_OU_rfun <- function(n, pars, LogicalRule, cross_talk = c("pos"=0.0,"neg"=0.0), dt = 0.001, max_time = 3.0,
                                 bridge = TRUE, bridge_exp_cutoff = -20, backend = c("auto", "Rcpp", "R"),
                                 adaptive = FALSE, adaptive_mode = c("distance", "ptol"),
                                 adapt_factor = 32.0, adapt_eps = 1.5,
                                 p_tol = 1e-9, curv_tol = 0.1) {
  backend <- match.arg(backend)
  adaptive_mode <- match.arg(adaptive_mode)
  use_rcpp <- switch(
    backend,
    Rcpp = TRUE,
    R = FALSE,
    auto = .ensure_ou_rcpp()
  )
  if (use_rcpp) {
    LogicalRules_OU_rfun_Rcpp(
      n = n, pars = pars, LogicalRule = LogicalRule, cross_talk = cross_talk,
      dt = dt, max_time = max_time, bridge = bridge, bridge_exp_cutoff = bridge_exp_cutoff,
      adaptive = adaptive, adaptive_mode = adaptive_mode,
      adapt_factor = adapt_factor, adapt_eps = adapt_eps,
      p_tol = p_tol, curv_tol = curv_tol
    )
  } else {
    if (isTRUE(adaptive)) {
      warning("adaptive=TRUE is currently only supported in backend='Rcpp'; using fixed dt for backend='R'.")
    }
    LogicalRules_OU_rfun_R(
      n = n, pars = pars, LogicalRule = LogicalRule, cross_talk = cross_talk,
      dt = dt, max_time = max_time, bridge = bridge, bridge_exp_cutoff = bridge_exp_cutoff
    )
  }
}
.ensure_ou_rcpp <- local({
  compiled <- FALSE
  function() {
    if (compiled && exists("LogicalRules_OU_finish_rcpp", mode = "function")) {
      return(TRUE)
    }
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
      return(FALSE)
    }
    ok <- tryCatch({
      Rcpp::sourceCpp(
        code = '
        #include <Rcpp.h>
        #include <vector>
        #include <cmath>
        using namespace Rcpp;

        // [[Rcpp::export]]
        NumericMatrix LogicalRules_OU_finish_rcpp(
            int n,
            NumericVector diag_a,
            NumericVector u,
            NumericVector sv,
            NumericVector th,
            NumericVector t0,
            double ct_pos,
            double ct_neg,
            double dt,
            double max_time,
            bool bridge,
            double bridge_exp_cutoff,
            bool adaptive,
            int adaptive_mode,
            double adapt_factor,
            double adapt_eps,
            double p_tol,
            double curv_tol
        ) {
          NumericMatrix finish(n, 4);
          std::fill(finish.begin(), finish.end(), R_PosInf);

          if (adaptive_mode != 2) adaptive_mode = 1;
          if (!R_finite(adapt_factor) || adapt_factor < 1.0) adapt_factor = 1.0;
          if (!R_finite(adapt_eps) || adapt_eps <= 0.0) adapt_eps = 1.5;
          if (!R_finite(p_tol) || p_tol <= 0.0 || p_tol >= 1.0) p_tol = 1e-9;
          if (!R_finite(curv_tol) || curv_tol < 0.0) curv_tol = 0.1;
          const double dt_min = std::max(dt / 32.0, 1e-8);
          const double dt_max = std::max(dt, dt * adapt_factor);

          const double thresh0 = th[0];
          const double thresh1 = -th[1];
          const double thresh2 = th[2];
          const double thresh3 = -th[3];
          const double sv_max = std::max(std::max(sv[0], sv[1]), std::max(sv[2], sv[3]));
          const double dist_ref = std::max(adapt_eps * sv_max, 1e-8);

          for (int tr = 0; tr < n; ++tr) {
            double x0 = 0.0;
            double x1 = 0.0;
            double x2 = 0.0;
            double x3 = 0.0;
            double t_curr = 0.0;

            bool done0 = false;
            bool done1 = false;
            bool done2 = false;
            bool done3 = false;
            bool A_resolved = false;
            bool B_resolved = false;

            while (t_curr < max_time && !(A_resolved && B_resolved)) {
              double dt_step = std::min(dt, max_time - t_curr);

              if (adaptive) {
                if (adaptive_mode == 1) {
                  double min_dist = R_PosInf;
                  if (!done0) min_dist = std::min(min_dist, std::fabs(thresh0 - x0));
                  if (!done1) min_dist = std::min(min_dist, std::fabs(thresh1 - x1));
                  if (!done2) min_dist = std::min(min_dist, std::fabs(thresh2 - x2));
                  if (!done3) min_dist = std::min(min_dist, std::fabs(thresh3 - x3));
                  if (!R_finite(min_dist)) min_dist = 0.0;

                  const double growth = std::pow(std::max(min_dist / dist_ref, 1.0), 2.0);
                  dt_step = std::min(dt_max, dt * growth);
                  dt_step = std::max(dt_min, dt_step);
                  dt_step = std::min(dt_step, max_time - t_curr);
                } else {
                  const double log_inv_p = std::max(std::log(1.0 / p_tol), 1e-12);
                  double dt_prob = dt_max;
                  double min_dist = R_PosInf;

                  if (!done0) {
                    const double d = std::fabs(thresh0 - x0);
                    min_dist = std::min(min_dist, d);
                    if (d <= 1e-12) {
                      dt_prob = dt_min;
                    } else {
                      dt_prob = std::min(dt_prob, 2.0 * d * d / std::max(sv[0] * sv[0] * log_inv_p, 1e-12));
                    }
                  }
                  if (!done1) {
                    const double d = std::fabs(thresh1 - x1);
                    min_dist = std::min(min_dist, d);
                    if (d <= 1e-12) {
                      dt_prob = dt_min;
                    } else {
                      dt_prob = std::min(dt_prob, 2.0 * d * d / std::max(sv[1] * sv[1] * log_inv_p, 1e-12));
                    }
                  }
                  if (!done2) {
                    const double d = std::fabs(thresh2 - x2);
                    min_dist = std::min(min_dist, d);
                    if (d <= 1e-12) {
                      dt_prob = dt_min;
                    } else {
                      dt_prob = std::min(dt_prob, 2.0 * d * d / std::max(sv[2] * sv[2] * log_inv_p, 1e-12));
                    }
                  }
                  if (!done3) {
                    const double d = std::fabs(thresh3 - x3);
                    min_dist = std::min(min_dist, d);
                    if (d <= 1e-12) {
                      dt_prob = dt_min;
                    } else {
                      dt_prob = std::min(dt_prob, 2.0 * d * d / std::max(sv[3] * sv[3] * log_inv_p, 1e-12));
                    }
                  }

                  dt_step = std::max(dt_min, std::min(dt_max, dt_prob));
                  dt_step = std::min(dt_step, max_time - t_curr);

                  // Deterministic curvature guard in interacting state-space.
                  for (int shrink = 0; shrink < 10; ++shrink) {
                    const double h = dt_step;
                    const double h2 = 0.5 * h;

                    const double f0_0 = diag_a[0] * x0 + ct_pos * x2 + u[0];
                    const double f0_1 = diag_a[1] * x1 + ct_neg * x3 + u[1];
                    const double f0_2 = ct_pos * x0 + diag_a[2] * x2 + u[2];
                    const double f0_3 = ct_neg * x1 + diag_a[3] * x3 + u[3];

                    const double e1_0 = x0 + h * f0_0;
                    const double e1_1 = x1 + h * f0_1;
                    const double e1_2 = x2 + h * f0_2;
                    const double e1_3 = x3 + h * f0_3;

                    const double m_0 = x0 + h2 * f0_0;
                    const double m_1 = x1 + h2 * f0_1;
                    const double m_2 = x2 + h2 * f0_2;
                    const double m_3 = x3 + h2 * f0_3;

                    const double fm_0 = diag_a[0] * m_0 + ct_pos * m_2 + u[0];
                    const double fm_1 = diag_a[1] * m_1 + ct_neg * m_3 + u[1];
                    const double fm_2 = ct_pos * m_0 + diag_a[2] * m_2 + u[2];
                    const double fm_3 = ct_neg * m_1 + diag_a[3] * m_3 + u[3];

                    const double e2_0 = m_0 + h2 * fm_0;
                    const double e2_1 = m_1 + h2 * fm_1;
                    const double e2_2 = m_2 + h2 * fm_2;
                    const double e2_3 = m_3 + h2 * fm_3;

                    double curv = std::fabs(e2_0 - e1_0);
                    curv = std::max(curv, std::fabs(e2_1 - e1_1));
                    curv = std::max(curv, std::fabs(e2_2 - e1_2));
                    curv = std::max(curv, std::fabs(e2_3 - e1_3));

                    if (!R_finite(min_dist)) min_dist = 0.0;
                    const double scale = std::max(min_dist, sv_max * std::sqrt(std::max(h, 1e-12)));
                    if (curv <= curv_tol * std::max(scale, 1e-12)) break;

                    dt_step *= 0.5;
                    if (dt_step <= dt_min) {
                      dt_step = dt_min;
                      break;
                    }
                  }
                }
              }

              const double sq_dt = std::sqrt(dt_step);
              const double prev0 = x0;
              const double prev1 = x1;
              const double prev2 = x2;
              const double prev3 = x3;

              x0 = prev0 + (diag_a[0] * prev0 + ct_pos * prev2 + u[0]) * dt_step + sv[0] * sq_dt * R::rnorm(0.0, 1.0);
              x1 = prev1 + (diag_a[1] * prev1 + ct_neg * prev3 + u[1]) * dt_step + sv[1] * sq_dt * R::rnorm(0.0, 1.0);
              x2 = prev2 + (ct_pos * prev0 + diag_a[2] * prev2 + u[2]) * dt_step + sv[2] * sq_dt * R::rnorm(0.0, 1.0);
              x3 = prev3 + (ct_neg * prev1 + diag_a[3] * prev3 + u[3]) * dt_step + sv[3] * sq_dt * R::rnorm(0.0, 1.0);

              for (int ch = 0; ch < 4; ++ch) {
                bool already_done = (ch == 0) ? done0 : (ch == 1) ? done1 : (ch == 2) ? done2 : done3;
                if (already_done) continue;

                const bool is_neg = (ch == 1 || ch == 3);
                const double thresh = (ch == 0) ? thresh0 : (ch == 1) ? thresh1 : (ch == 2) ? thresh2 : thresh3;
                const double x_prev = (ch == 0) ? prev0 : (ch == 1) ? prev1 : (ch == 2) ? prev2 : prev3;
                const double x_next = (ch == 0) ? x0 : (ch == 1) ? x1 : (ch == 2) ? x2 : x3;

                const bool hit_straddle = is_neg ? (x_next <= thresh) : (x_next >= thresh);
                bool hit_bridge = false;
                bool hit_now = hit_straddle;

                if (bridge && !hit_straddle) {
                  const double d0 = thresh - x_prev;
                  const double d1 = thresh - x_next;
                  const double denom = std::max(sv[ch] * sv[ch] * dt_step, 1e-12);
                  const double exponent = -2.0 * d0 * d1 / denom;
                  if (exponent > bridge_exp_cutoff) {
                    const double p_cross = std::exp(exponent);
                    if (R::runif(0.0, 1.0) < p_cross) {
                      hit_bridge = true;
                      hit_now = true;
                    }
                  }
                }

                if (hit_now) {
                  double hit_t = t_curr + dt_step;
                  if (hit_straddle) {
                    const double denom_cross = x_next - x_prev;
                    double frac = (thresh - x_prev) / denom_cross;
                    if (!R_finite(frac)) frac = 1.0;
                    if (frac < 0.0) frac = 0.0;
                    if (frac > 1.0) frac = 1.0;
                    hit_t = t_curr + frac * dt_step;
                  } else if (hit_bridge) {
                    hit_t = t_curr + R::runif(0.0, 1.0) * dt_step;
                  }

                  finish(tr, ch) = hit_t + t0[ch];
                  if (ch == 0) done0 = true;
                  if (ch == 1) done1 = true;
                  if (ch == 2) done2 = true;
                  if (ch == 3) done3 = true;
                  if (ch < 2) {
                    A_resolved = true;
                  } else {
                    B_resolved = true;
                  }
                }
              }

              t_curr += dt_step;
            }
          }

          return finish;
        }'
      )
      TRUE
    }, error = function(e) {
      FALSE
    })
    compiled <<- isTRUE(ok) && exists("LogicalRules_OU_finish_rcpp", mode = "function")
    compiled
  }
})

LogicalRules_OU_rfun_R <- function(n, pars, LogicalRule, cross_talk = c("pos"=0.0,"neg"=0.0), dt = 0.001, max_time = 3.0,
                                   bridge = TRUE, bridge_exp_cutoff = -20) {
  # Simulates T&W 2004 4-channel OU model with CROSS-TALK
  # 
  # cross_talk: The interaction parameter a_ij (e.g., 0.75 or -0.75 in paper)
  #             Applied between the two POSITIVE channels (A+ and B+).
  #             Extended here to facilitate equivalent interaction on the NEGATIVE channels (A-, B-)
  #             to accommodate the "correct absence" pair developed herein
  
  racers <- c("A", "n_A", "B", "n_B")
  # pars columns: c("a", "u", "sv", "b", "t0")
  # 'a' in pars is the diagonal decay (e.g., -1.0)
  # 'u' in pars is the constant input per accumulator
  # sv, b and t0 are the Weiner diffusive noise, threshold (absolute) and non-decision times which we keep fixed across accumulators
  #
  # Performance notes:
  # - Only simulate trials that are not yet resolved at the (A+/A-) and (B+/B-) race level.
  # - Avoid full-matrix copies (`X_prev <- X`) each step; keep a per-step snapshot only for active trials.
  # - Exploit the sparse interaction structure (A+<->B+, A-<->B-) to avoid a general matrix multiply.
  
  diag_a <- pars[racers, "a"]
  p_u  <- pars[racers, "u"]
  p_sv <- pars[racers, "sv"]
  p_th <- pars[racers, "b"]
  p_t0 <- pars[racers, "t0"]
  
  ct_pos <- unname(cross_talk["pos"])
  ct_neg <- unname(cross_talk["neg"])
  
  is_neg <- c(FALSE, TRUE, FALSE, TRUE)
  thresh_vec <- p_th
  thresh_vec[is_neg] <- -thresh_vec[is_neg]
  bridge_denom <- (p_sv ^ 2) * dt
  
  X <- matrix(0, nrow = n, ncol = 4)
  FinishTimes <- matrix(Inf, nrow = n, ncol = 4)
  colnames(FinishTimes) <- racers
  channel_done <- matrix(FALSE, nrow = n, ncol = 4)
  A_resolved <- rep(FALSE, n)
  B_resolved <- rep(FALSE, n)
  
  sq_dt <- sqrt(dt)
  n_steps <- ceiling(max_time / dt)
  
  active_trials <- seq_len(n)
  
  for (step in seq_len(n_steps)) {
    if (length(active_trials) == 0) break
    t_curr <- step * dt
    
    # Remove trials whose A-race and B-race are both resolved.
    # (For resolved trials we only need their recorded FinishTimes.)
    still_unresolved <- !(A_resolved[active_trials] & B_resolved[active_trials])
    if (!all(still_unresolved)) {
      active_trials <- active_trials[still_unresolved]
      if (length(active_trials) == 0) break
    }
    
    X_prev <- X[active_trials, , drop = FALSE]
    
    dW <- matrix(rnorm(length(active_trials) * 4), nrow = length(active_trials), ncol = 4)
    
    drift_interaction <- X_prev
    drift_interaction[, 1] <- diag_a[1] * X_prev[, 1] + ct_pos * X_prev[, 3]
    drift_interaction[, 3] <- ct_pos * X_prev[, 1] + diag_a[3] * X_prev[, 3]
    drift_interaction[, 2] <- diag_a[2] * X_prev[, 2] + ct_neg * X_prev[, 4]
    drift_interaction[, 4] <- ct_neg * X_prev[, 2] + diag_a[4] * X_prev[, 4]
    
    drift_total <- drift_interaction
    drift_total[, 1] <- (drift_interaction[, 1] + p_u[1]) * dt
    drift_total[, 2] <- (drift_interaction[, 2] + p_u[2]) * dt
    drift_total[, 3] <- (drift_interaction[, 3] + p_u[3]) * dt
    drift_total[, 4] <- (drift_interaction[, 4] + p_u[4]) * dt
    
    noise <- dW
    noise[, 1] <- dW[, 1] * p_sv[1] * sq_dt
    noise[, 2] <- dW[, 2] * p_sv[2] * sq_dt
    noise[, 3] <- dW[, 3] * p_sv[3] * sq_dt
    noise[, 4] <- dW[, 4] * p_sv[4] * sq_dt
    
    X_new <- X_prev + drift_total + noise
    X[active_trials, ] <- X_new
    
    for (i in 1:4) {
      active_pos <- which(!channel_done[active_trials, i])
      if (length(active_pos) == 0) next
      
      thresh <- thresh_vec[i]
      x_prev <- X_prev[active_pos, i]
      x_curr <- X_new[active_pos, i]
      
      if (is_neg[i]) {
        hit_straddle <- x_curr <= thresh
      } else {
        hit_straddle <- x_curr >= thresh
      }
      hit_now <- hit_straddle
      hit_bridge <- rep(FALSE, length(hit_now))
      
      if (bridge) {
        still_pos <- which(!hit_straddle)
        if (length(still_pos) > 0) {
          d0 <- thresh - x_prev[still_pos]
          d1 <- thresh - x_curr[still_pos]
          exponent <- -2.0 * d0 * d1 / bridge_denom[i]
          
          # If we're far from threshold, bridge crossing is effectively 0; skip exp/runif.
          near <- exponent > bridge_exp_cutoff
          if (any(near)) {
            near_pos <- still_pos[near]
            p_cross <- exp(exponent[near])
            bridge_hits <- runif(length(near_pos)) < p_cross
            if (any(bridge_hits)) {
              bridge_idx <- near_pos[bridge_hits]
              hit_now[bridge_idx] <- TRUE
              hit_bridge[bridge_idx] <- TRUE
            }
          }
        }
      }
      
      if (any(hit_now)) {
        # Sub-step hit times:
        # - Straddle hits: linear interpolation between (x_prev, x_curr).
        # - Bridge-only hits: assume a uniform crossing time within the step.
        hit_t <- rep(t_curr, length(hit_now))
        
        if (any(hit_straddle)) {
          idx <- which(hit_straddle)
          denom <- x_curr[idx] - x_prev[idx]
          frac <- (thresh - x_prev[idx]) / denom
          frac[!is.finite(frac)] <- 1.0
          frac <- pmin(1.0, pmax(0.0, frac))
          hit_t[idx] <- (t_curr - dt) + frac * dt
        }
        
        if (any(hit_bridge)) {
          idx <- which(hit_bridge)
          hit_t[idx] <- (t_curr - dt) + runif(length(idx)) * dt
        }
        
        new_hits <- active_trials[active_pos[hit_now]]
        FinishTimes[new_hits, i] <- hit_t[hit_now] + p_t0[i]
        channel_done[new_hits, i] <- TRUE
        if (i <= 2) {
          A_resolved[new_hits] <- TRUE
        } else {
          B_resolved[new_hits] <- TRUE
        }
      }
    }
  }
  
  apply_logic_gates(
    LogicalRule = LogicalRule,
    A_t  = FinishTimes[, "A"],
    nA_t = FinishTimes[, "n_A"],
    B_t  = FinishTimes[, "B"],
    nB_t = FinishTimes[, "n_B"]
  )
}

LogicalRules_OU_rfun_Rcpp <- function(n, pars, LogicalRule, cross_talk = c("pos"=0.0,"neg"=0.0), dt = 0.001, max_time = 3.0,
                                      bridge = TRUE, bridge_exp_cutoff = -20,
                                      adaptive = FALSE, adaptive_mode = c("distance", "ptol"),
                                      adapt_factor = 32.0, adapt_eps = 1.5,
                                      p_tol = 1e-9, curv_tol = 0.1) {
  if (!.ensure_ou_rcpp()) {
    stop("Rcpp backend not available; set backend='R' or install Rcpp.")
  }
  adaptive_mode <- match.arg(adaptive_mode)
  adaptive_mode_i <- if (adaptive_mode == "ptol") 2L else 1L
  racers <- c("A", "n_A", "B", "n_B")
  diag_a <- pars[racers, "a"]
  p_u  <- pars[racers, "u"]
  p_sv <- pars[racers, "sv"]
  p_th <- pars[racers, "b"]
  p_t0 <- pars[racers, "t0"]
  
  ct_pos <- unname(cross_talk["pos"])
  ct_neg <- unname(cross_talk["neg"])
  
  FinishTimes <- LogicalRules_OU_finish_rcpp(
    n = as.integer(n),
    diag_a = as.numeric(diag_a),
    u = as.numeric(p_u),
    sv = as.numeric(p_sv),
    th = as.numeric(p_th),
    t0 = as.numeric(p_t0),
    ct_pos = as.numeric(ct_pos),
    ct_neg = as.numeric(ct_neg),
    dt = as.numeric(dt),
    max_time = as.numeric(max_time),
    bridge = isTRUE(bridge),
    bridge_exp_cutoff = as.numeric(bridge_exp_cutoff),
    adaptive = isTRUE(adaptive),
    adaptive_mode = as.integer(adaptive_mode_i),
    adapt_factor = as.numeric(adapt_factor),
    adapt_eps = as.numeric(adapt_eps),
    p_tol = as.numeric(p_tol),
    curv_tol = as.numeric(curv_tol)
  )
  colnames(FinishTimes) <- racers
  
  apply_logic_gates(
    LogicalRule = LogicalRule,
    A_t  = FinishTimes[, "A"],
    nA_t = FinishTimes[, "n_A"],
    B_t  = FinishTimes[, "B"],
    nB_t = FinishTimes[, "n_B"]
  )
}


