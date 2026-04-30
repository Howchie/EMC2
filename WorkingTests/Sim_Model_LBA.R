LogicalRules_rfun_lba <- function(n, pars, LogicalRule) {
  racers <- c("A", "n_A", "B", "n_B")
  if (!all(racers %in% rownames(pars))) stop("pars must have rows: A, n_A, B, n_B")
  req_cols <- c("A", "b", "t0", "v", "sv")
  if (!all(req_cols %in% colnames(pars))) stop("pars must have cols: A, b, t0, v, sv")
  if (!LogicalRule %in% c("OR", "AND", "XOR", "ID")) stop("LogicalRule must be one of: OR, AND, XOR, ID")

  # Vectorized across racers to reduce function call overhead
  n_racers <- length(racers)
  all_rt <- rtdists::rlba_norm(
    n = n * n_racers,
    A = rep(pars[racers, "A"], each = n),
    b = rep(pars[racers, "b"], each = n),
    t0 = rep(pars[racers, "t0"], each = n),
    mean_v = rep(pars[racers, "v"], each = n),
    sd_v = rep(pars[racers, "sv"], each = n),
    st0 = 0,
    posdrift = TRUE
  )[, "rt"]

  Rrti <- matrix(all_rt, nrow = n, ncol = n_racers, dimnames = list(NULL, racers))

  apply_logic_gates(
    LogicalRule,
    A_t = Rrti[, "A"],
    nA_t = Rrti[, "n_A"],
    B_t = Rrti[, "B"],
    nB_t = Rrti[, "n_B"]
  )
}

.lba_init_state <- function(p_vec) {
  vc_yes <- if ("vc_yes" %in% names(p_vec)) .param_get(p_vec, "vc_yes") else NA_real_
  b <- if ("b" %in% names(p_vec)) .param_get(p_vec, "b") else NA_real_
  vc_yesA_H <- if ("vc_yesA_H" %in% names(p_vec)) {
    .param_get(p_vec, "vc_yesA_H")
  } else if ("vc_yesA" %in% names(p_vec)) {
    .param_get(p_vec, "vc_yesA")
  } else if ("vc_yesH" %in% names(p_vec)) {
    .param_get(p_vec, "vc_yesH")
  } else if (is.finite(vc_yes)) {
    vc_yes
  } else {
    stop("p_vec must include one of: vc_yesA_H, vc_yesA, or vc_yes.")
  }
  vc_yesA_L <- if ("vc_yesA_L" %in% names(p_vec)) {
    .param_get(p_vec, "vc_yesA_L")
  } else if ("vc_yesL" %in% names(p_vec)) {
    .param_get(p_vec, "vc_yesL")
  } else if ("vc_yesA" %in% names(p_vec)) {
    .param_get(p_vec, "vc_yesA")
  } else {
    vc_yesA_H
  }

  vc_yesB_H <- if ("vc_yesB_H" %in% names(p_vec)) {
    .param_get(p_vec, "vc_yesB_H")
  } else if ("vc_yesB" %in% names(p_vec)) {
    .param_get(p_vec, "vc_yesB")
  } else if ("vc_yesH" %in% names(p_vec)) {
    .param_get(p_vec, "vc_yesH")
  } else if (is.finite(vc_yes)) {
    vc_yes
  } else {
    stop("p_vec must include one of: vc_yesB_H, vc_yesHm, vc_yesB, or vc_yes.")
  }
  vc_yesB_L <- if ("vc_yesB_L" %in% names(p_vec)) {
    .param_get(p_vec, "vc_yesB_L")
  } else if ("vc_yesL" %in% names(p_vec)) {
    .param_get(p_vec, "vc_yesL")
  } else if ("vc_yesB" %in% names(p_vec)) {
    .param_get(p_vec, "vc_yesB")
  } else {
    vc_yesB_H
  }
  b_yes <- if ("b_yes" %in% names(p_vec)) {
    .param_get(p_vec, "b_yes")
  } else {
    b
  }
  b_no <- if ("b_no" %in% names(p_vec)) {
    .param_get(p_vec, "b_no")
  } else {
    b
  }

  list(
    vc_yesA_H = vc_yesA_H,
    vc_yesA_L = vc_yesA_L,
    vc_yesB_H = vc_yesB_H,
    vc_yesB_L = vc_yesB_L,
    b_yes = b_yes,
    b_no = b_no,
    cap_target = if ("capacity_target" %in% names(p_vec)) .param_get(p_vec, "capacity_target") else if ("capacity" %in% names(p_vec)) .param_get(p_vec, "capacity") else 0,
    cap_absence = if ("capacity_absence" %in% names(p_vec)) .param_get(p_vec, "capacity_absence") else 0
  )
}

.lba_cell_pars <- function(p_vec, meta_row, state) {
  pars <- matrix(NA_real_, nrow = 4, ncol = 5)
  rownames(pars) <- c("A", "n_A", "B", "n_B")
  colnames(pars) <- c("A", "b", "t0", "v", "sv")

  l1 <- as.integer(meta_row$Channel1)
  l2 <- as.integer(meta_row$Channel2)
  A_present <- l1 > 0L
  B_present <- l2 > 0L

  pars[, "A"] <- .param_get(p_vec, "A")
  pars[, "t0"] <- .param_get(p_vec, "t0")
  pars[c("A", "B"), "b"] <- state$b_yes
  pars[c("n_A", "n_B"), "b"] <- state$b_no

  pars["A", "v"] <- if (l1 == 2L) {
    state$vc_yesA_H
  } else if (l1 == 1L) {
    state$vc_yesA_L
  } else {
    .param_get(p_vec, "ve")
  }
  pars["n_A", "v"] <- if (A_present) .param_get(p_vec, "ve") else .param_get(p_vec, "vc_no")
  pars["B", "v"] <- if (l2 == 2L) {
    state$vc_yesB_H
  } else if (l2 == 1L) {
    state$vc_yesB_L
  } else {
    .param_get(p_vec, "ve")
  }
  pars["n_B", "v"] <- if (B_present) .param_get(p_vec, "ve") else .param_get(p_vec, "vc_no")

  if (A_present && B_present) {
    pars[c("A", "B"), "v"] <- pars[c("A", "B"), "v"] + state$cap_target
  } else if (!A_present && !B_present) {
    pars[c("n_A", "n_B"), "v"] <- pars[c("n_A", "n_B"), "v"] + state$cap_absence
  }

  pars["A", "sv"] <- if (A_present) .param_get(p_vec, "sv_c") else .param_get(p_vec, "sv_e")
  pars["n_A", "sv"] <- if (A_present) .param_get(p_vec, "sv_e") else .param_get(p_vec, "sv_c")
  pars["B", "sv"] <- if (B_present) .param_get(p_vec, "sv_c") else .param_get(p_vec, "sv_e")
  pars["n_B", "sv"] <- if (B_present) .param_get(p_vec, "sv_e") else .param_get(p_vec, "sv_c")

  pars
}

.lba_simulate_rule_cell <- function(n, pars, logical_rule) {
  LogicalRules_rfun_lba(n = n, pars = pars, LogicalRule = logical_rule)
}
