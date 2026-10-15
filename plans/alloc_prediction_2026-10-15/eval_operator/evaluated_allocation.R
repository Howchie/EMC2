.emc_alloc_act_part <-
function (Xf, add, fixed) 
{
    if (fixed <= 0 || !nrow(Xf)) 
        return(NULL)
    fixed * emc_alloc_act_term_cpp(Xf, rep_len(as.numeric(add), 
        nrow(Xf)))
}
.emc_alloc_act_rows <-
function (units) 
list(shear = which(vapply(units, function(u) identical(u$kind, 
    "shear"), logical(1))), hmc = which(vapply(units, function(u) identical(u$kind, 
    "subject") && !is.null(u$h), logical(1))))
.emc_alloc_act_sum <-
function (shear, hmc, n_t) 
{
    out <- numeric(n_t)
    if (!is.null(shear)) 
        out <- out + shear
    if (!is.null(hmc)) 
        out <- out + hmc
    out
}
.emc_alloc_activation <-
function (units, X, r_min, shear_fixed, hmc_fixed) 
{
    rows <- .emc_alloc_act_rows(units)
    .emc_alloc_act_sum(.emc_alloc_act_part(X[rows$shear, , drop = FALSE], 
        0, shear_fixed), .emc_alloc_act_part(X[rows$hmc, , drop = FALSE], 
        r_min[rows$hmc], hmc_fixed), ncol(X))
}
.emc_alloc_as_subject <-
function (u) 
{
    if (!identical(u$kind, "hmc")) 
        return(u)
    list(kind = "subject", rho = u$rho, p = NULL, c_p = 0, h = u$h, 
        c_h = u$cost)
}
.emc_alloc_block_iters <-
function (chain) 
{
    n <- chain$alloc$window$n
    if (is.null(n) || !is.finite(n)) 
        0
    else as.numeric(n)
}
.emc_alloc_cov <-
function (x) 
{
    if (ncol(x) < 2L) 
        return(matrix(NA_real_, nrow(x), nrow(x)))
    xc <- x - rowMeans(x)
    tcrossprod(xc)/(ncol(x) - 1L)
}
.emc_alloc_credit <-
function (units, g_floor, ratio) 
{
    for (k in seq_along(g_floor)) {
        x <- units[[k]]
        s <- x$rho - g_floor[k]
        units[[k]]$rho <- if (ratio <= 1 || x$rho >= 1) 
            min(1, g_floor[k] + s * ratio)
        else 1 - (1 - x$rho) * exp(-s * (ratio - 1)/(1 - x$rho))
    }
    units
}
.emc_alloc_data_rows <-
function (d) 
{
    if (is.data.frame(d)) 
        return(as.numeric(nrow(d)))
    if (is.list(d)) 
        return(sum(vapply(d, .emc_alloc_data_rows, numeric(1))))
    1
}
.emc_alloc_dev_factor <-
function (q, r, P, H) 
{
    both <- max(0, q + r - 1)
    (1 - q - r + both) + (q - both) * P + (r - both) * H + both * 
        P * H
}
.emc_alloc_dev_ratio <-
function (units, q, r) 
{
    num <- den <- 0
    for (i in seq_along(units)) {
        x <- units[[i]]
        if (!identical(x$kind, "subject") || is.null(x$d_rho) || 
            !all(is.finite(x$d_rho))) 
            next
        k <- length(x$d_rho)
        q0 <- if (is.null(x$q_old)) 
            1
        else x$q_old
        r0 <- if (is.null(x$r_old)) 
            0
        else x$r_old
        H <- if (is.null(x$h_dev)) 
            rep(1, k)
        else x$h_dev
        P <- x$p_dev
        qi <- q[i]
        if (is.null(P)) {
            P <- rep(1, k)
            qi <- q0
        }
        ri <- if (is.null(x$h_dev)) 
            0
        else r[i]
        base <- x$d_rho/pmax(1 - q0 * (1 - P), 9.9999999999999995e-07)
        num <- num + sum(base * .emc_alloc_dev_factor(qi, ri, 
            P, H))
        den <- den + sum(base * .emc_alloc_dev_factor(q0, r0, 
            P, H))
    }
    if (den > 0) 
        num/den
    else NA_real_
}
.emc_alloc_draw_stats <-
function (emc, vdraws, dev, cache = NULL) 
{
    st <- if (!is.null(cache) && identical(cache$vdraws, vdraws)) 
        cache$stats
    else NULL
    if (is.null(st)) {
        p <- emc[[1L]]$n_pars
        mu <- do.call(cbind, lapply(seq_along(emc), function(i) emc[[i]]$samples$theta_mu[, 
            vdraws[[i]], drop = FALSE]))
        lsd <- do.call(cbind, lapply(seq_along(emc), function(i) .emc_alloc_lsd_draws(emc[[i]]$samples$theta_var, 
            vdraws[[i]])))
        al <- lapply(seq_along(emc), function(i) aperm(emc[[i]]$samples$alpha[, 
            , vdraws[[i]], drop = FALSE], c(1L, 3L, 2L)))
        st <- list(mu = mu, mu_cov = .emc_alloc_cov(mu), lsd_var = .emc_alloc_row_var(lsd), 
            al = al, a_cov = lapply(seq_len(dim(al[[1L]])[3L]), 
                function(s) .emc_alloc_cov(do.call(cbind, lapply(al, 
                  function(a) matrix(a[, , s], p))))))
    }
    if (dev && is.null(st$d_cov)) {
        p <- nrow(st$mu)
        st$d_cov <- lapply(seq_along(st$a_cov), function(s) .emc_alloc_cov(do.call(cbind, 
            lapply(st$al, function(a) matrix(a[, , s], p))) - 
            st$mu))
    }
    if (!is.null(cache)) {
        cache$vdraws <- vdraws
        cache$stats <- st
    }
    st
}
.emc_alloc_draws <-
function (chain, m) 
{
    i <- chain$samples$idx
    seq.int(max(1L, i - max(2L, m) + 1L), i)
}
.emc_alloc_eff <-
function (rho, kappa) 
{
    if (is.null(kappa)) 
        return(rho)
    .emc_alloc_rho_of_tau(kappa * .emc_alloc_tau(rho))
}
.emc_alloc_eff_inv <-
function (rs, kappa) 
{
    if (is.null(kappa)) 
        return(rs)
    target <- .emc_alloc_rho_of_tau(.emc_alloc_tau(rs)/kappa)
    ifelse(rs < .emc_alloc_eff(0, kappa), -1, target)
}
.emc_alloc_empty_window <-
function (n, p) 
{
    list(n = 0L, nG = 0L, G = matrix(0, 2L * p, 2L * p), nA = 0L, 
        A = replicate(n, matrix(0, p, p), simplify = FALSE), 
        D = replicate(n, matrix(0, p, p), simplify = FALSE), 
        cpu_rest = 0, fam = .emc_alloc_fam_empty(), cpu_hmc = numeric(n), 
        hmc_att = integer(n), rows_hmc = numeric(n), P = replicate(n, 
            matrix(0, p, p), simplify = FALSE), nP = integer(n), 
        cpu_part = numeric(n))
}
.EMC_ALLOC_EVIDENCE_DECAY <-
0.5
.emc_alloc_factor <-
function (g, x, a = g) 
{
    m <- floor(x)
    f <- x - m
    (1 - f) * .emc_alloc_factor_k(g, a, m) + f * .emc_alloc_factor_k(g, 
        a, m + 1)
}
.emc_alloc_factor_inverse <-
function (g, a, t) 
{
    out <- rep(Inf, length(t))
    if (g >= 1 - 1.0000000000000001e-09) 
        return(out)
    geo <- a >= g - 9.9999999999999998e-13
    if (geo) {
        ok <- t > 0
        m <- ceiling(log(t[ok])/log(g) - 1.0000000000000001e-09)
    }
    else {
        c <- (1 - g)/(1 - a)
        s <- 1 - (1 - t)/c
        ok <- s > 0
        m <- if (a <= 0) 
            rep(1, sum(ok))
        else ceiling(log(s[ok])/log(a) - 1.0000000000000001e-09)
    }
    m <- pmax(1, m)
    f0 <- .emc_alloc_factor_k(g, a, m - 1)
    f1 <- .emc_alloc_factor_k(g, a, m)
    out[ok] <- (m - 1) + pmin(1, pmax(0, (f0 - t[ok])/(f0 - f1)))
    out
}
.emc_alloc_factor_k <-
function (g, a, k) 
{
    n <- max(length(g), length(a), length(k))
    g <- rep_len(g, n)
    a <- rep_len(a, n)
    k <- rep_len(k, n)
    out <- g^k
    sat <- a < g - 9.9999999999999998e-13
    out[sat] <- 1 - (1 - g[sat])/(1 - a[sat]) * (1 - a[sat]^k[sat])
    out
}
.emc_alloc_fam_empty <-
function () 
{
    stats::setNames(lapply(.EMC_ALLOC_FAMILIES, function(f) list(cpu = 0, 
        own = 0, iters = 0L, moves = 0L)), .EMC_ALLOC_FAMILIES)
}
.EMC_ALLOC_FAMILIES <-
c("shear", "collapsed", "hmc")
.emc_alloc_geyer_tau <-
function (X) 
{
    n <- min(vapply(X, nrow, integer(1)))
    m <- length(X)
    d <- ncol(X[[1L]])
    if (n < 4L) 
        return(rep(NA_real_, d))
    X <- lapply(X, function(x) x[nrow(x) - n + seq_len(n), , 
        drop = FALSE])
    N <- 2L^ceiling(log2(2L * n))
    acov <- Reduce(`+`, lapply(X, function(x) {
        xc <- sweep(x, 2L, colMeans(x))
        f <- stats::mvfft(rbind(xc, matrix(0, N - n, d)))
        Re(stats::mvfft(Mod(f)^2, inverse = TRUE))[seq_len(n), 
            , drop = FALSE]/(N * n)
    }))/m
    w <- acov[1L, ] * n/(n - 1)
    b <- if (m > 1L) 
        .emc_alloc_row_var(matrix(vapply(X, colMeans, numeric(d)), 
            d))
    else numeric(d)
    vplus <- w * (n - 1)/n + b
    vapply(seq_len(d), function(i) {
        if (!is.finite(vplus[i]) || vplus[i] <= 0) 
            return(NA_real_)
        r <- 1 - (w[i] - acov[, i])/vplus[i]
        k <- floor(n/2)
        P <- r[2L * seq_len(k) - 1L] + r[2L * seq_len(k)]
        stop_at <- which(P <= 0)
        if (length(stop_at)) 
            P <- P[seq_len(stop_at[1L] - 1L)]
        if (!length(P)) 
            return(1)
        max(1, -1 + 2 * sum(cummin(P)))
    }, numeric(1))
}
.emc_alloc_hmc_factor <-
function (ev, W, var, min_n) 
{
    if (is.null(ev) || ev$n < min_n) 
        return(NULL)
    esjd <- rowSums((W %*% ev$M) * W)/ev$n
    pmin(1, pmax(0, 1 - esjd/(2 * pmax(var, .Machine$double.xmin))))
}
.EMC_ALLOC_HMC_REGIMES <-
c(0.5, 1)
.emc_alloc_joint_grid <-
function (cap) 
{
    x <- seq(0.125, 4, by = 0.125)
    if (cap >= 4.5) 
        x <- c(x, seq(4.5, cap, by = 0.5))
    sort(unique(c(x[x <= cap], cap)))
}
.emc_alloc_kappa <-
function (tau_stage, rho_whole) 
{
    k <- tau_stage/.emc_alloc_tau(rho_whole)
    ifelse(is.finite(k), pmax(1, k), 1)
}
.emc_alloc_khist_add <-
function (K, rows, k) 
{
    if (!length(rows) || k < 1L) 
        return(K)
    if (ncol(K) < k) 
        K <- cbind(K, matrix(0L, nrow(K), k - ncol(K)))
    K[rows, k] <- K[rows, k] + 1L
    K
}
.emc_alloc_khist_sum <-
function (A, B) 
{
    if (is.null(A)) 
        return(B)
    if (is.null(B)) 
        return(A)
    w <- max(ncol(A), ncol(B))
    pad <- function(M) cbind(M, matrix(0L, nrow(M), w - ncol(M)))
    pad(A) + pad(B)
}
.emc_alloc_level <-
function (u) 
{
    if (identical(u$kind, "shear") && !is.null(u$floor) && u$floor > 
        0) 
        return(u$rho * .emc_alloc_factor(u$g, u$floor, if (is.null(u$a)) u$g else u$a))
    if (!identical(u$kind, "subject")) 
        return(u$rho)
    base <- .emc_alloc_subject_base(u)
    if (!is.null(u$r_min) && !is.null(u$h)) 
        base * (1 - u$r_min * (1 - u$h))
    else base
}
.emc_alloc_lsd <-
function (V) 
0.5 * log(pmax(diag(as.matrix(V)), .Machine$double.xmin))
.emc_alloc_lsd_draws <-
function (V, k) 
{
    p <- dim(V)[1L]
    idx <- c(outer(seq_len(p) + (seq_len(p) - 1L) * p, (k - 1) * 
        p * p, "+"))
    matrix(0.5 * log(pmax(V[idx], .Machine$double.xmin)), p, 
        dimnames = list(dimnames(V)[[1L]], NULL))
}
.emc_alloc_marginal <-
function (w, fam, last = NULL) 
{
    f <- if (is.null(w$fam)) 
        NULL
    else w$fam[[fam]]
    if (is.null(f)) 
        f <- list(cpu = 0, own = 0, iters = 0L, moves = 0L)
    carried <- if (is.list(last) && isTRUE(is.finite(last$fixed[[fam]]))) 
        last$fixed[[fam]]
    else 0
    list(c = if (f$moves >= 5L) f$own/f$moves else NA_real_, 
        fixed = if (f$iters >= .EMC_ALLOC_MIN_PAIRS) max(0, f$cpu - 
            f$own)/f$iters else carried)
}
.EMC_ALLOC_MAX_ROUNDS <-
4
.emc_alloc_measure <-
function (emc, j, hmc_state, Vbar, stage = NULL, cache = NULL) 
{
    ch <- emc[[j]]
    p <- ch$n_pars
    n <- ch$n_subjects
    w <- ch$alloc$window
    stored <- is.null(w) || w$n < .EMC_ALLOC_MIN_PAIRS
    m <- if (stored) 
        200L
    else w$n
    draws <- lapply(emc, .emc_alloc_draws, m = m)
    production <- identical(stage, "sample")
    pool <- NULL
    if (production && identical(ch$alloc$stage, "sample")) {
        pool <- ch$alloc$pool
        if (!stored) {
            if (is.list(pool)) {
                w$G <- w$G + pool$G
                w$nG <- w$nG + pool$nG
            }
            pool <- list(G = w$G, nG = w$nG)
        }
    }
    vdraws <- draws
    if (production) {
        prod <- lapply(emc, function(x) {
            st <- x$samples$stage[seq_len(x$samples$idx)]
            which(st == "sample")
        })
        if (all(lengths(prod) >= lengths(draws))) 
            vdraws <- prod
    }
    if (stored) {
        k <- draws[[j]]
        w <- .emc_alloc_empty_window(n, p)
        if (length(k) >= 2L) {
            s <- ch$samples
            for (t in k[-1L]) {
                d <- c(s$theta_mu[, t], .emc_alloc_lsd(s$theta_var[, 
                  , t])) - c(s$theta_mu[, t - 1L], .emc_alloc_lsd(s$theta_var[, 
                  , t - 1L]))
                if (all(is.finite(d))) {
                  w$G <- w$G + tcrossprod(d)
                  w$nG <- w$nG + 1L
                }
                da <- matrix(s$alpha[, , t] - s$alpha[, , t - 
                  1L], p)
                if (all(is.finite(da))) {
                  for (ss in seq_len(n)) w$A[[ss]] <- w$A[[ss]] + 
                    tcrossprod(da[, ss])
                  w$nA <- w$nA + 1L
                }
            }
        }
    }
    shear <- .emc_shear_active(ch) && is.list(ch$shear$plan) && 
        is.matrix(ch$shear$plan$dirs)
    dirs <- if (shear) 
        ch$shear$plan$dirs
    else diag(p)
    Dinv <- tryCatch(solve(dirs), error = function(e) NULL)
    if (is.null(Dinv)) {
        shear <- FALSE
        dirs <- diag(p)
        Dinv <- diag(p)
    }
    T_g <- rbind(cbind(Dinv, matrix(0, p, p)), cbind(matrix(0, 
        p, p), diag(p)))
    has_D <- !stored && !is.null(w$D) && w$nA > 0
    ds <- .emc_alloc_draw_stats(emc, vdraws, has_D, cache)
    g_var <- c(.emc_alloc_quad_diag(Dinv, ds$mu_cov), ds$lsd_var)
    g_esjd <- if (w$nG > 0) 
        rowSums((T_g %*% w$G) * T_g)/w$nG
    else rep(0, 2L * p)
    g_rho <- if (w$nG > 0) 
        .emc_alloc_rho(g_esjd, g_var)
    else rep(0, 2L * p)
    g_tau_stage <- if (!stored) 
        .emc_alloc_stage_tau(emc, Dinv, stage)
    else rep(NA_real_, 2L * p)
    g_kappa <- .emc_alloc_kappa(g_tau_stage, .emc_alloc_whole_rho(emc, 
        draws, Dinv, g_var))
    g_fac <- g_sat <- rep(NA_real_, 2L * p)
    g_measured <- logical(2L * p)
    if (shear) {
        last <- .emc_shear_record(ch$shear$plan$last, p)
        if (production && identical(ch$shear$plan$stage, "sample")) {
            ks <- p + seq_len(p)
            for (i in seq_along(emc)) {
                if (i == j) 
                  next
                o <- emc[[i]]$shear$plan
                if (!is.list(o) || !identical(o$stage, "sample") || 
                  !is.list(o$last)) 
                  next
                ol <- .emc_shear_record(o$last, p)
                last$E[ks] <- last$E[ks] + ol$E[ks]
                last$n[ks] <- last$n[ks] + ol$n[ks]
                last$moved[ks] <- last$moved[ks] + ol$moved[ks]
                last$N[ks] <- last$N[ks] + ol$N[ks]
                oK <- ol$K
                oK[seq_len(p), ] <- 0
                last$K <- .emc_alloc_khist_sum(last$K, oK)
            }
        }
        nk <- last$n
        last_n <- nk
        last_iters <- last$iters
        g_measured <- nk > 0
        g_obs <- ifelse(g_measured, pmin(1, pmax(0, 1 - (last$E/pmax(nk, 
            9.9999999999999998e-13))/(2 * g_var))), NA_real_)
        g_fac <- .emc_alloc_shrink_g(g_obs, nk, rep(1:2, each = p))
        g_sat <- .emc_alloc_saturation(g_fac, last$N/pmax(1, 
            rowSums(last$K))/(2 * g_var), last$K)
    }
    coll <- .emc_collapsed_active(ch) && is.list(ch$collapsed$plan) && 
        is.matrix(ch$collapsed$plan$C)
    c_fac <- c_sat <- rep(.EMC_ALLOC_SHEAR_PRIOR_G, 2L * p)
    c_measured <- FALSE
    if (coll) {
        last <- ch$collapsed$plan$last
        c_measured <- is.list(last) && last$n > 0
        if (c_measured) {
            last <- .emc_collapsed_record(last, p)
            c_obs <- pmin(1, pmax(0, 1 - (rowSums((T_g %*% last$E) * 
                T_g)/last$n)/(2 * g_var)))
            n0 <- .EMC_ALLOC_ROUND_EVIDENCE
            c_fac <- (last$n * c_obs + n0 * mean(c_obs))/(last$n + 
                n0)
            c_sat <- .emc_alloc_saturation(c_fac, rowSums((T_g %*% 
                last$N) * T_g)/max(1, sum(last$K))/(2 * g_var), 
                last$K)
        }
    }
    lik <- ch$lik_laplace
    W <- lapply(seq_len(n), function(s) {
        if (!is.null(hmc_state$subjects[[s]]$W)) 
            return(hmc_state$subjects[[s]]$W)
        prec <- if (length(lik) >= s && is.list(lik[[s]])) 
            lik[[s]]$prec
        else NULL
        Ws <- if (is.matrix(prec) && all(dim(prec) == p) && all(is.finite(prec))) 
            .emc_hmc_ref_whiten(prec, Vbar)
        if (is.null(Ws)) 
            diag(p)
        else Ws
    })
    s_rho <- matrix(0, n, p)
    s_var <- matrix(NA_real_, n, p)
    d_var <- matrix(NA_real_, n, p)
    d_rho <- matrix(NA_real_, n, p)
    for (s in seq_len(n)) {
        s_var[s, ] <- .emc_alloc_quad_diag(W[[s]], ds$a_cov[[s]])
        if (w$nA > 0) 
            s_rho[s, ] <- .emc_alloc_rho(rowSums((W[[s]] %*% 
                w$A[[s]]) * W[[s]])/w$nA, s_var[s, ])
        if (has_D) {
            d_var[s, ] <- .emc_alloc_quad_diag(W[[s]], ds$d_cov[[s]])
            d_rho[s, ] <- .emc_alloc_rho(rowSums((W[[s]] %*% 
                w$D[[s]]) * W[[s]])/w$nA, d_var[s, ])
        }
    }
    r_bar <- if (has_D && any(is.finite(d_rho))) 
        mean(d_rho[is.finite(d_rho)])
    else NA_real_
    prev <- ch$alloc$particle
    nP <- if (stored || is.null(w$nP)) 
        integer(n)
    else w$nP
    p_fac <- matrix(NA_real_, n, p)
    p_dev <- matrix(NA_real_, n, p)
    c_part <- rep(NA_real_, n)
    for (s in seq_len(n)) {
        if (nP[s] >= .EMC_ALLOC_MIN_PAIRS) {
            esjd <- rowSums((W[[s]] %*% w$P[[s]]) * W[[s]])/nP[s]
            p_fac[s, ] <- pmin(1, pmax(0, 1 - esjd/(2 * pmax(s_var[s, 
                ], .Machine$double.xmin))))
            if (all(is.finite(d_var[s, ]))) 
                p_dev[s, ] <- pmin(1, pmax(0, 1 - esjd/(2 * pmax(d_var[s, 
                  ], .Machine$double.xmin))))
            c_part[s] <- w$cpu_part[s]/nP[s]
        }
        else if (is.matrix(prev$fac) && all(dim(prev$fac) == 
            c(n, p)) && all(is.finite(prev$fac[s, ])) && is.finite(prev$cost[s])) {
            p_fac[s, ] <- prev$fac[s, ]
            if (is.matrix(prev$fac_dev) && all(dim(prev$fac_dev) == 
                c(n, p))) 
                p_dev[s, ] <- prev$fac_dev[s, ]
            c_part[s] <- prev$cost[s]
        }
    }
    q_old <- if (is.null(prev$rate)) 
        rep(1, n)
    else prev$rate
    r_old <- ch$alloc$last$rate
    if (length(r_old) != n || !all(is.finite(r_old))) 
        r_old <- numeric(n)
    t_rho <- s_rho
    for (s in which(is.finite(c_part))) {
        den <- 1 - q_old[s] * (1 - p_fac[s, ])
        t_rho[s, ] <- ifelse(den > 9.9999999999999995e-07, pmin(1, 
            s_rho[s, ]/pmax(den, 9.9999999999999995e-07)), 1)
    }
    part_in_rest <- if (stored || is.null(w$cpu_part)) 
        0
    else sum(w$cpu_part[is.finite(c_part)])
    C_rest <- if (!stored && w$n > 0 && w$cpu_rest > part_in_rest) 
        (w$cpu_rest - part_in_rest)/w$n
    else NA_real_
    data_rows <- vapply(ch$data, .emc_alloc_data_rows, numeric(1))
    hmc_rows <- if (!stored && length(w$rows_hmc) == n) 
        sum(w$rows_hmc * data_rows)
    else 0
    hmc_row_cpu <- if (hmc_rows > 0) 
        sum(w$cpu_hmc)/hmc_rows
    else {
        carried <- ch$alloc$last$hmc_row_cpu
        if (isTRUE(is.finite(carried))) 
            carried
        else NA_real_
    }
    marg <- stats::setNames(lapply(.EMC_ALLOC_FAMILIES, function(f) .emc_alloc_marginal(if (stored) 
        NULL
    else w, f, ch$alloc$last)), .EMC_ALLOC_FAMILIES)
    list(stored = stored, pairs = c(group = w$nG, subjects = w$nA), 
        pool = pool, shear = shear, g_rho = g_rho, g_var = g_var, 
        g_fac = g_fac, g_kappa = g_kappa, g_tau_stage = g_tau_stage, 
        g_sat = g_sat, g_measured = g_measured, g_n = if (shear) last_n else numeric(2L * 
            p), g_iters = if (shear) last_iters else numeric(2L * 
            p), c_n = if (coll && c_measured) last$n else 0, 
        c_iters = if (coll && c_measured) last$iters else 0, 
        coll = coll, c_fac = c_fac, c_sat = c_sat, c_measured = c_measured, 
        c_coll = if (stored) NA_real_ else marg$collapsed$c, 
        coll_fixed = marg$collapsed$fixed, c_rate = if (!stored && 
            w$n > 0 && !is.null(w$fam)) w$fam$collapsed$moves/w$n else 0, 
        W = W, s_rho = s_rho, t_rho = t_rho, p_fac = p_fac, p_dev = p_dev, 
        c_part = c_part, q_old = q_old, r_old = r_old, s_var = s_var, 
        C_rest = C_rest, d_var = d_var, d_rho = d_rho, r_bar = r_bar, 
        g_floor = if (is.finite(r_bar) && r_bar < 1) pmin(g_rho, 
            pmax(0, (g_rho - r_bar)/(1 - r_bar))) else g_rho, 
        c_round = if (stored) NA_real_ else marg$shear$c, shear_fixed = marg$shear$fixed, 
        c_hmc = if (!stored) ifelse(w$hmc_att > 0, w$cpu_hmc/pmax(1L, 
            w$hmc_att), NA_real_) else rep(NA_real_, n), data_rows = data_rows, 
        hmc_row_cpu = hmc_row_cpu, hmc_fixed = marg$hmc$fixed, 
        fixed = list(shear = marg$shear$fixed, collapsed = marg$collapsed$fixed, 
            hmc = marg$hmc$fixed))
}
.EMC_ALLOC_MIN_GAIN <-
0.90000000000000002
.EMC_ALLOC_MIN_MOVE <-
0.001
.EMC_ALLOC_MIN_PAIRS <-
10L
.emc_alloc_movable <-
function (units, max_rounds, joint = NULL, parts = NULL) 
{
    if (is.null(parts)) 
        parts <- .emc_alloc_range_parts(units, max_rounds)
    movable <- .emc_alloc_range_of(parts)$movable
    if (!is.null(joint)) {
        idx <- as.integer(joint$idx)
        g <- rep_len(joint$g, length(idx))
        movable[idx] <- movable[idx] | (is.finite(g) & g < 1 - 
            1.0000000000000001e-09)
    }
    movable
}
.emc_alloc_need_layout <-
function (us, max_rounds) 
{
    kind <- vapply(us, function(u) u$kind, character(1))
    shear <- kind == "shear" & lengths(lapply(us, `[[`, "rho")) == 
        1L
    sh <- us[shear]
    field <- function(f, d) vapply(sh, function(u) if (is.null(u[[f]])) 
        d
    else as.numeric(u[[f]]), numeric(1))
    list(kind = kind, shear = shear, g = field("g", NA_real_), 
        a = vapply(sh, function(u) as.numeric(if (is.null(u$a)) u$g else u$a), 
            numeric(1)), cost = field("cost", NA_real_), cap = field("cap", 
            max_rounds), floor = field("floor", 0), kappa = field("kappa", 
            1))
}
.emc_alloc_particle_held <-
function (units) 
{
    any(vapply(units, function(x) identical(x$kind, "subject") && 
        !is.null(x$q_fixed) && is.finite(x$q_fixed) && x$q_fixed < 
        1, logical(1)))
}
.EMC_ALLOC_PARTICLE_RATES <-
c(1, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0)
.emc_alloc_plan_rates <-
function (units, sol) 
{
    n_u <- length(units)
    q <- if (length(sol$q) == n_u) 
        sol$q
    else rep(1, n_u)
    r <- vapply(seq_len(n_u), function(i) {
        x <- units[[i]]
        if (!identical(x$kind, "subject") || is.null(x$h)) 
            return(0)
        min(1, sol$x[i] + (if (is.null(x$r_min)) 0 else x$r_min))
    }, numeric(1))
    list(q = q, r = r)
}
.emc_alloc_predict <-
function (units, sol, joint = NULL, max_rounds = .EMC_ALLOC_MAX_ROUNDS) 
{
    xc <- if (is.null(sol$joint)) 
        0
    else sol$joint
    idx <- if (is.null(joint)) 
        integer(0)
    else as.integer(joint$idx)
    fac_c <- if (xc > 0 && length(idx)) 
        rep_len(.emc_alloc_factor(joint$g, xc, if (is.null(joint$a)) 
            joint$g
        else joint$a), length(idx))
    lapply(seq_along(units), function(i) {
        u <- .emc_alloc_as_subject(units[[i]])
        rho <- u$rho
        k <- match(i, idx)
        if (!is.null(fac_c) && !is.na(k)) 
            rho <- rho * fac_c[k]
        x <- if (length(sol$x) >= i && is.finite(sol$x[i])) 
            sol$x[i]
        else 0
        switch(u$kind, fixed = .emc_alloc_eff(rho, u$kappa), 
            shear = .emc_alloc_eff(rho * .emc_alloc_factor(u$g, 
                x, if (is.null(u$a)) u$g else u$a), u$kappa), 
            subject = {
                q <- if (is.null(u$p) || length(sol$q) < i) 1 else sol$q[i]
                r <- if (is.null(u$h)) 0 else min(1, x + (if (is.null(u$r_min)) 0 else u$r_min))
                rho * .emc_alloc_dev_factor(q, r, if (is.null(u$p)) 1 else u$p, 
                  if (is.null(u$h)) 1 else u$h)
            }, rho)
    })
}
.emc_alloc_quad_diag <-
function (A, S) 
rowSums((A %*% S) * A)
.emc_alloc_range <-
function (units, max_rounds, movable = NULL) 
.emc_alloc_range_of(.emc_alloc_range_parts(units, max_rounds), 
    movable)
.emc_alloc_range_of <-
function (parts, movable = NULL) 
{
    levels <- parts$levels
    best_rho <- parts$best
    if (is.null(movable) || length(movable) != length(levels)) 
        movable <- best_rho < levels - 1.0000000000000001e-09
    list(levels = levels, best = best_rho, movable = movable, 
        hi = if (any(movable)) max(levels[movable], 0) else 0, 
        lo = if (any(movable)) max(best_rho[movable], 0) else 0, 
        hi_all = max(levels, 0))
}
.emc_alloc_range_parts <-
function (units, max_rounds) 
list(levels = vapply(units, function(u) .emc_alloc_eff(max(.emc_alloc_level(u)), 
    u$kappa), numeric(1)), best = vapply(units, function(u) switch(u$kind, 
    fixed = .emc_alloc_eff(max(u$rho), u$kappa), shear = .emc_alloc_eff(u$rho * 
        .emc_alloc_factor(u$g, if (is.null(u$cap)) max_rounds else u$cap, 
            if (is.null(u$a)) u$g else u$a), u$kappa), subject = max(if (is.null(u$h)) .emc_alloc_level(u) else .emc_alloc_subject_base(u) * 
        u$h)), numeric(1)))
.emc_alloc_rebuild <-
function (emc, stage, n_cores, B) 
{
    if (identical(stage, "preburn") || !.emc_alloc_wanted(emc[[1L]])) 
        return(emc)
    emc <- .emc_prepare_chain_rng(emc)
    nch <- length(emc)
    Vbar <- lapply(emc, function(ch) {
        k <- .emc_alloc_draws(ch, if (is.null(ch$alloc$window)) 
            200L
        else max(2L, ch$alloc$window$n))
        apply(ch$samples$theta_var[, , k, drop = FALSE], 1:2, 
            mean)
    })
    hst <- lapply(seq_len(nch), function(j) .emc_hmc_plan_chain(emc[[j]], 
        stage, Vbar[[j]]))
    stats_cache <- new.env()
    ms <- lapply(seq_len(nch), function(j) .emc_alloc_measure(emc, 
        j, hst[[j]], Vbar[[j]], stage, cache = stats_cache))
    shear_on <- vapply(seq_len(nch), function(j) ms[[j]]$shear, 
        logical(1))
    jobs <- NULL
    for (j in seq_len(nch)) {
        if (!isTRUE(hst[[j]]$supported)) 
            next
        u <- .emc_alloc_units(ms[[j]], hst[[j]], shear_on[j], 
            optimistic = TRUE, B = B)
        sol <- .emc_alloc_solve_units(u)
        x <- sol$x
        regime <- if (is.null(sol$regime)) 
            0
        else sol$regime
        n_g <- length(ms[[j]]$g_rho)
        for (s in seq_len(emc[[j]]$n_subjects)) {
            sj <- hst[[j]]$subjects[[s]]
            if (!isTRUE(sj$eligible) || isTRUE(sj$warmup)) 
                next
            if (regime + x[n_g + s] <= 0) 
                next
            ev <- hst[[j]]$evidence[[s]]
            lacking <- is.null(ev) || ev$n < .EMC_HMC_MIN_EVIDENCE
            bracket <- is.na(sj$eps)
            pilots <- if (lacking) 
                .EMC_HMC_PILOTS
            else 0L
            if (bracket || pilots > 0L) 
                jobs <- rbind(jobs, data.frame(chain = j, subject = s, 
                  bracket = bracket, pilots = pilots))
        }
    }
    if (!is.null(jobs)) {
        res <- .emc_with_preserved_rng(auto_mclapply(seq_len(nrow(jobs)), 
            function(k) {
                j <- jobs$chain[k]
                .emc_quiet_fork(.emc_hmc_calibrate(emc[[j]], 
                  hst[[j]], jobs$subject[k], bracket = jobs$bracket[k], 
                  pilots = jobs$pilots[k]))
            }, mc.cores = min(nrow(jobs), max(1L, n_cores)), 
            mc.set.seed = FALSE))
        for (k in seq_len(nrow(jobs))) {
            j <- jobs$chain[k]
            s <- jobs$subject[k]
            r <- res[[k]]
            if (!is.list(r) || is.null(r$seed)) 
                next
            .emc_reject_add(r$rejects, "pilot")
            emc[[j]]$rng$subjects[[s]] <- r$seed
            hst[[j]] <- .emc_hmc_apply_calibration(hst[[j]], 
                s, r)
        }
    }
    plan_chain <- function(j, choice = NULL, q_fixed = NULL) {
        ch <- emc[[j]]
        m <- ms[[j]]
        st <- hst[[j]]
        if (!is.null(st) && isTRUE(st$supported)) 
            for (s in seq_len(ch$n_subjects)) {
                if (isTRUE(st$subjects[[s]]$eligible) && !is.finite(st$subjects[[s]]$eps)) {
                  st$subjects[[s]]$eligible <- FALSE
                  st$subjects[[s]]$reason <- "no usable step size"
                }
            }
        u <- .emc_alloc_units(m, st, shear_on[j], optimistic = FALSE, 
            B = B)
        n_g <- length(m$g_rho)
        if (!is.null(q_fixed)) 
            for (s in seq_len(ch$n_subjects)) {
                k <- n_g + s
                if (!is.null(u$units[[k]]$p) && !is.null(u$units[[k]]$h)) 
                  u$units[[k]]$q_fixed <- q_fixed[s]
            }
        u$regime_fixed <- choice
        sol <- .emc_alloc_solve_units(u)
        if (shear_on[j]) 
            ch$shear$plan$rounds <- sol$x[seq_len(n_g)]
        if (isTRUE(m$coll)) 
            ch$collapsed$plan$rounds <- sol$joint
        regime <- if (is.null(sol$regime)) 
            0
        else sol$regime
        rate <- pmin(1, regime + sol$x[n_g + seq_len(ch$n_subjects)])
        q <- rep(1, ch$n_subjects)
        if (!is.null(st) && isTRUE(st$supported)) {
            for (s in seq_len(ch$n_subjects)) {
                sj <- st$subjects[[s]]
                sj$rate <- if (!isTRUE(sj$eligible)) 
                  0
                else if (isTRUE(sj$warmup)) 
                  .EMC_HMC_WARMUP_RATE
                else rate[s]
                if (isTRUE(sj$eligible) && !isTRUE(sj$warmup) && 
                  length(sol$q) >= n_g + s) 
                  q[s] <- sol$q[n_g + s]
                st$subjects[[s]] <- sj
            }
        }
        rec <- list(stage = stage, B = B, iteration = ch$samples$idx, 
            transitions = if (is.null(ch$sample_transitions)) NA_integer_ else as.integer(ch$sample_transitions), 
            stored = m$stored, pairs = m$pairs, C_rest = u$C_rest, 
            c_round = u$c_round, c_hmc = m$c_hmc, hmc_row_cpu = m$hmc_row_cpu, 
            c_part = m$c_part, hmc_fixed = u$hmc_fixed, shear_fixed = u$shear_fixed, 
            fixed = m$fixed, tau = sol$tau, tau0 = sol$tau0, 
            tau0_all = sol$tau0_all, cost = sol$cost, regime = regime, 
            choice = sol$choice, options = sol$options, r_bar = m$r_bar, 
            objective = sol$objective, objective0 = sol$objective0, 
            pred_group = if (length(sol$pred) >= n_g) vapply(sol$pred[seq_len(n_g)], 
                max, numeric(1)), pred_subject = if (length(sol$pred) >= 
                n_g + ch$n_subjects) do.call(rbind, sol$pred[n_g + 
                seq_len(ch$n_subjects)]), pred0_group = m$g_rho, 
            pred0_subject = do.call(rbind, lapply(u$units[n_g + 
                seq_len(ch$n_subjects)], .emc_alloc_subject_base)), 
            group_rho = m$g_rho, group_g = m$g_fac, group_a = m$g_sat, 
            group_kappa = m$g_kappa, group_tau_stage = m$g_tau_stage, 
            rounds = if (shear_on[j]) sol$x[seq_len(n_g)] else numeric(n_g), 
            collapsed_rounds = sol$joint, collapsed_g = if (isTRUE(m$coll)) m$c_fac, 
            collapsed_a = if (isTRUE(m$coll)) m$c_sat, collapsed_cost = if (is.null(u$joint)) NA_real_ else u$joint$cost, 
            collapsed_eps = if (isTRUE(m$coll)) ch$collapsed$plan$eps else NA_real_, 
            subject_rho = apply(m$s_rho, 1L, max), hmc_h = vapply(u$units[n_g + 
                seq_len(ch$n_subjects)], function(x) if (!is.null(x$h)) max(x$h) else NA_real_, 
                numeric(1)), rate = if (is.null(st$subjects)) numeric(ch$n_subjects) else vapply(st$subjects, 
                function(x) x$rate, numeric(1)), particle_rate = q, 
            particle_fixed = !is.null(q_fixed), choice_fixed = !is.null(choice), 
            particle_p = apply(m$p_fac, 1L, function(x) if (all(is.finite(x))) max(x) else NA_real_), 
            eps = if (is.null(st$subjects)) rep(NA_real_, ch$n_subjects) else vapply(st$subjects, 
                function(x) x$eps, numeric(1)), calib_cpu = if (is.null(st$subjects)) 0 else sum(vapply(st$subjects, 
                function(x) x$calib_cpu, numeric(1))))
        ch$alloc <- list(window = .emc_alloc_empty_window(ch$n_subjects, 
            ch$n_pars), stage = stage, pool = m$pool, last = rec, 
            history = c(ch$alloc$history, list(rec)), particle = list(rate = q, 
                fac = m$p_fac, fac_dev = m$p_dev, cost = m$c_part))
        if (!is.null(st)) {
            st$alloc <- rec
            ch$subject_hmc <- st
        }
        ch
    }
    for (j in seq_len(nch)) emc[[j]] <- plan_chain(j)
    opts <- lapply(emc, function(ch) ch$alloc$last$options)
    choice <- NULL
    if (all(lengths(opts) > 0L)) {
        keys <- Reduce(intersect, lapply(opts, names))
        if (length(keys)) {
            tot <- Reduce(`+`, lapply(opts, function(o) o[keys]))
            choice <- keys[which.min(tot)]
            for (j in seq_len(nch)) if (!identical(emc[[j]]$alloc$last$choice, 
                choice)) 
                emc[[j]] <- plan_chain(j, choice)
        }
    }
    q <- do.call(pmax, lapply(emc, function(ch) ch$alloc$particle$rate))
    for (j in seq_len(nch)) if (any(emc[[j]]$alloc$particle$rate < 
        q - 1.0000000000000001e-09)) 
        emc[[j]] <- plan_chain(j, choice, q_fixed = q)
    for (j in seq_len(nch)) {
        emc[[j]]$alloc$particle$rate <- q
        emc[[j]]$alloc$last$particle_shared <- q
        h <- emc[[j]]$alloc$history
        h[[length(h)]]$particle_shared <- q
        emc[[j]]$alloc$history <- h
    }
    emc
}
.emc_alloc_record_cost <-
function (sampler, rest, tasks, moves = list()) 
{
    w <- .emc_alloc_window(sampler)
    w$n <- w$n + 1L
    w$cpu_rest <- w$cpu_rest + rest
    if (is.null(w$fam)) 
        w$fam <- .emc_alloc_fam_empty()
    work <- 0
    attempts <- 0L
    for (s in seq_along(tasks)) {
        tk <- tasks[[s]]
        if (is.null(tk)) 
            next
        c_s <- if (!is.null(tk$work)) 
            tk$work
        else if (!is.null(tk$cost)) 
            tk$cost
        else tk$cpu
        w$cpu_hmc[s] <- w$cpu_hmc[s] + c_s
        w$hmc_att[s] <- w$hmc_att[s] + 1L
        if (is.null(w$rows_hmc)) 
            w$rows_hmc <- numeric(length(tasks))
        w$rows_hmc[s] <- w$rows_hmc[s] + if (is.null(tk$rows)) 
            0
        else tk$rows
        work <- work + c_s
        attempts <- attempts + 1L
    }
    if (attempts > 0L) 
        moves$hmc <- c(cpu = if (is.null(moves$hmc)) work else as.numeric(moves$hmc)[1], 
            own = work, n = attempts)
    for (fam in intersect(names(moves), .EMC_ALLOC_FAMILIES)) {
        m <- moves[[fam]]
        if (is.null(m) || length(m) < 3L || !is.finite(m[[3]]) || 
            m[[3]] <= 0) 
            next
        f <- w$fam[[fam]]
        f$cpu <- f$cpu + as.numeric(m[[1]])
        f$own <- f$own + as.numeric(m[[2]])
        f$iters <- f$iters + 1L
        f$moves <- f$moves + as.integer(round(m[[3]]))
        w$fam[[fam]] <- f
    }
    sampler$alloc$window <- w
    sampler
}
.emc_alloc_record_group <-
function (sampler, pars) 
{
    idx <- sampler$samples$idx
    w <- .emc_alloc_window(sampler)
    d <- c(as.numeric(pars$tmu), .emc_alloc_lsd(pars$tvar)) - 
        c(sampler$samples$theta_mu[, idx], .emc_alloc_lsd(sampler$samples$theta_var[, 
            , idx]))
    if (all(is.finite(d))) {
        w$G <- w$G + tcrossprod(d)
        w$nG <- w$nG + 1L
    }
    sampler$alloc$window <- w
    sampler
}
.emc_alloc_record_particle <-
function (sampler, before, after, ran, cost, times) 
{
    w <- .emc_alloc_window(sampler)
    if (is.null(w$nP)) {
        e <- .emc_alloc_empty_window(sampler$n_subjects, sampler$n_pars)
        w$P <- e$P
        w$nP <- e$nP
        w$cpu_part <- e$cpu_part
    }
    idx <- which(ran)
    if (!length(idx)) 
        return(sampler)
    t_ran <- pmax(0, times[idx])
    share <- if (sum(t_ran) > 0) 
        t_ran/sum(t_ran)
    else rep(1/length(idx), length(idx))
    for (k in seq_along(idx)) {
        s <- idx[k]
        d <- after[, s] - before[, s]
        if (!all(is.finite(d))) 
            next
        w$P[[s]] <- w$P[[s]] + tcrossprod(d)
        w$nP[s] <- w$nP[s] + 1L
        w$cpu_part[s] <- w$cpu_part[s] + if (sum(t_ran) > 0) 
            t_ran[k]
        else cost * share[k]
    }
    sampler$alloc$window <- w
    sampler
}
.emc_alloc_record_subjects <-
function (sampler, alpha, mu = NULL) 
{
    idx <- sampler$samples$idx
    w <- .emc_alloc_window(sampler)
    if (is.null(w$D)) 
        w$D <- replicate(ncol(alpha), matrix(0, nrow(alpha), 
            nrow(alpha)), simplify = FALSE)
    old <- matrix(sampler$samples$alpha[, , idx], nrow(alpha))
    d <- alpha - old
    if (all(is.finite(d))) {
        dd <- if (is.null(mu)) 
            d
        else d - (as.numeric(mu) - as.numeric(sampler$samples$theta_mu[, 
            idx]))
        for (s in seq_len(ncol(d))) {
            w$A[[s]] <- w$A[[s]] + tcrossprod(d[, s])
            w$D[[s]] <- w$D[[s]] + tcrossprod(dd[, s])
        }
        w$nA <- w$nA + 1L
    }
    sampler$alloc$window <- w
    sampler
}
.emc_alloc_regime <-
function (u, r) 
{
    units <- u$units
    cost <- 0
    covered <- FALSE
    q <- r_new <- numeric(length(units))
    for (i in seq_along(units)) {
        x <- units[[i]]
        if (!identical(x$kind, "subject")) 
            next
        q[i] <- if (is.null(x$q_old)) 
            1
        else x$q_old
        if (is.null(x$d_rho) || !all(is.finite(x$d_rho)) || is.null(x$h_dev) || 
            is.null(x$h)) 
            next
        covered <- TRUE
        units[[i]]$r_min <- r
        r_new[i] <- r
        cost <- cost + r * x$c_h
    }
    if (!covered) 
        return(NULL)
    ratio <- .emc_alloc_dev_ratio(units, q, r_new)
    if (!is.finite(ratio)) 
        return(NULL)
    list(units = .emc_alloc_credit(units, u$g_floor, ratio), 
        cost = cost, credit = ratio)
}
.emc_alloc_rho <-
function (esjd, var) 
{
    pmin(1, pmax(0, 1 - esjd/(2 * pmax(var, .Machine$double.xmin))))
}
.emc_alloc_rho_of_tau <-
function (tau) 
{
    tau <- pmax(tau, 1)
    (tau - 1)/(tau + 1)
}
.EMC_ALLOC_ROUND_EVIDENCE <-
20
.emc_alloc_row_var <-
function (x) 
{
    if (ncol(x) < 2L) 
        return(rep(NA_real_, nrow(x)))
    rowSums((x - rowMeans(x))^2)/(ncol(x) - 1L)
}
.emc_alloc_saturation <-
function (g, D, K) 
{
    a <- g
    if (is.null(K) || !length(D)) 
        return(a)
    K <- as.matrix(K)
    for (i in seq_along(g)) {
        w <- K[if (nrow(K) == 1L) 
            1L
        else i, ]
        if (sum(w) <= 0 || sum(w[-1L]) <= 0 || !is.finite(D[i]) || 
            g[i] >= 1 - 1.0000000000000001e-09) 
            next
        w <- w/sum(w)
        k <- seq_along(w)
        h <- function(x) (1 - g[i]) * sum(w * (1 - x^k))/(1 - 
            x)
        if (D[i] <= h(0)) {
            a[i] <- 0
            next
        }
        if (D[i] >= h(g[i])) 
            next
        lo <- 0
        hi <- g[i]
        for (it in 1:60) {
            mid <- (lo + hi)/2
            if (h(mid) < D[i]) 
                lo <- mid
            else hi <- mid
        }
        a[i] <- (lo + hi)/2
    }
    a
}
.emc_alloc_scaled_parts <-
function (units, max_rounds) 
{
    kind <- vapply(units, function(u) u$kind, character(1))
    if (!all(kind %in% c("shear", "fixed")) || any(lengths(lapply(units, 
        `[[`, "rho")) != 1L)) 
        return(function(us) .emc_alloc_range_parts(us, max_rounds))
    sh <- kind == "shear"
    a_of <- function(u) if (is.null(u$a)) 
        u$g
    else u$a
    f_floor <- vapply(units, function(u) if (identical(u$kind, 
        "shear") && !is.null(u$floor) && u$floor > 0) 
        .emc_alloc_factor(u$g, u$floor, a_of(u))
    else NA_real_, numeric(1))
    f_cap <- vapply(units, function(u) if (identical(u$kind, 
        "shear")) 
        .emc_alloc_factor(u$g, if (is.null(u$cap)) 
            max_rounds
        else u$cap, a_of(u))
    else NA_real_, numeric(1))
    floored <- !is.na(f_floor)
    kappa <- vapply(units, function(u) if (is.null(u$kappa)) 
        1
    else u$kappa, numeric(1))
    function(us) {
        rho <- as.numeric(unlist(lapply(us, `[[`, "rho")))
        levels <- rho
        levels[floored] <- rho[floored] * f_floor[floored]
        best <- rho
        best[sh] <- rho[sh] * f_cap[sh]
        list(levels = .emc_alloc_eff(levels, kappa), best = .emc_alloc_eff(best, 
            kappa))
    }
}
.emc_alloc_settle <-
function (u, sol, units, credit, C_rest, movable, shear_act, 
    hmc_fixed) 
{
    rates <- .emc_alloc_plan_rates(units, sol)
    ratio <- .emc_alloc_dev_ratio(units, rates$q, rates$r)
    if (!is.finite(ratio) || ratio <= credit * (1 + 9.9999999999999995e-07)) 
        return(sol)
    n_g <- length(u$g_floor)
    held <- .emc_alloc_particle_held(units)
    raw <- u$units
    hold <- function(us, q, r) {
        for (i in seq_along(us)) {
            if (!identical(us[[i]]$kind, "subject") || is.null(us[[i]]$h)) 
                next
            if (!is.null(us[[i]]$p)) 
                us[[i]]$q_fixed <- q[i]
            us[[i]]$r_fixed <- r[i]
        }
        us
    }
    solve <- function(us) .emc_alloc_solve(us, C_rest, joint = u$joint, 
        movable = movable, shear_fixed = shear_act, hmc_fixed = hmc_fixed)
    own <- units
    own[seq_len(n_g)] <- .emc_alloc_credit(raw[seq_len(n_g)], 
        u$g_floor, ratio)
    own <- solve(hold(own, rates$q, rates$r))
    if (held) {
        own$regime <- sol$regime
        return(own)
    }
    all_q <- rep(1, length(units))
    kept <- units
    ratio_k <- .emc_alloc_dev_ratio(units, all_q, rates$r)
    if (is.finite(ratio_k)) 
        kept[seq_len(n_g)] <- .emc_alloc_credit(raw[seq_len(n_g)], 
            u$g_floor, ratio_k)
    kept <- solve(hold(kept, all_q, rates$r))
    best <- if (isTRUE(own$objective < kept$objective)) 
        own
    else kept
    best$regime <- sol$regime
    best
}
.EMC_ALLOC_SHEAR_PRIOR_G <-
0.5
.emc_alloc_shrink_g <-
function (g_obs, n, kind) 
{
    out <- rep(.EMC_ALLOC_SHEAR_PRIOR_G, length(n))
    meas <- n > 0 & is.finite(g_obs)
    if (!any(meas)) 
        return(out)
    l <- log(pmin(1, pmax(.EMC_ALLOC_MIN_MOVE, 1 - g_obs)))
    v <- 2/pmax(n, 9.9999999999999998e-13)
    t2_min <- 2/.EMC_ALLOC_ROUND_EVIDENCE
    spread <- function(sel) {
        if (sum(sel) < 3L) 
            return(NA_real_)
        w <- 1/v[sel]
        lbar <- sum(w * l[sel])/sum(w)
        sum(w * (l[sel] - lbar)^2)/sum(w) - sum(w * v[sel])/sum(w)
    }
    t2_all <- spread(meas)
    centre <- function(sel, t2) {
        w <- 1/(t2 + v[sel])
        sum(w * l[sel])/sum(w)
    }
    for (kk in unique(kind)) {
        sel <- meas & kind == kk
        t2 <- spread(sel)
        if (!is.finite(t2)) 
            t2 <- t2_all
        t2 <- max(t2_min, if (is.finite(t2)) t2 else t2_min)
        ref <- if (any(sel)) 
            centre(sel, t2)
        else centre(meas, t2)
        lk <- ifelse(sel, (l/v + ref/t2)/(1/v + 1/t2), ref)
        m <- exp(lk[kind == kk])
        out[kind == kk] <- ifelse(m <= .EMC_ALLOC_MIN_MOVE * 
            (1 + 1.0000000000000001e-09), 1, 1 - m)
    }
    out
}
.emc_alloc_solve <-
function (units, C_rest, max_rounds = .EMC_ALLOC_MAX_ROUNDS, 
    joint = NULL, movable = NULL, shear_fixed = 0, hmc_fixed = 0) 
{
    units <- lapply(units, .emc_alloc_as_subject)
    n_u <- length(units)
    moved_idx <- if (is.null(joint)) 
        integer(0)
    else as.integer(joint$idx)
    j_floor <- if (is.null(joint) || is.null(joint$floor)) 
        0
    else min(joint$floor, joint$cap)
    xcs <- if (is.null(joint)) 
        0
    else {
        grid <- .emc_alloc_joint_grid(joint$cap)
        unique(c(j_floor, grid[grid > j_floor]))
    }
    cand <- lapply(xcs, function(xc) {
        if (xc == 0) 
            return(units)
        fac <- .emc_alloc_factor(joint$g, xc, if (is.null(joint$a)) 
            joint$g
        else joint$a)
        moved <- units
        for (i in seq_along(moved_idx)) moved[[moved_idx[i]]]$rho <- moved[[moved_idx[i]]]$rho * 
            fac[i]
        moved
    })
    parts0 <- .emc_alloc_range_parts(units, max_rounds)
    if (is.null(movable)) 
        movable <- .emc_alloc_movable(units, max_rounds, joint, 
            parts0)
    moved_parts <- .emc_alloc_scaled_parts(units[moved_idx], 
        max_rounds)
    ranges <- lapply(seq_along(cand), function(c) {
        if (xcs[c] == 0 || !length(moved_idx)) 
            return(.emc_alloc_range_of(parts0, movable))
        pm <- moved_parts(cand[[c]][moved_idx])
        pc <- parts0
        pc$levels[moved_idx] <- pm$levels
        pc$best[moved_idx] <- pm$best
        .emc_alloc_range_of(pc, movable)
    })
    r0 <- ranges[[1L]]
    c_p0 <- sum(vapply(units, function(u) if (identical(u$kind, 
        "subject") && !is.null(u$p) && !is.null(u$c_p)) u$c_p else 0, 
        numeric(1)))
    x0 <- vapply(units, function(u) if (identical(u$kind, "shear") && 
        !is.null(u$floor)) 
        u$floor
    else 0, numeric(1))
    r_min <- vapply(units, function(u) if (identical(u$kind, 
        "subject") && !is.null(u$h) && !is.null(u$r_min)) 
        u$r_min
    else 0, numeric(1))
    activation <- function(X) .emc_alloc_activation(units, X, 
        r_min, shear_fixed, hmc_fixed)
    j_fixed <- function(xc) if (xc > 0 && !is.null(joint$fixed)) 
        joint$fixed * min(1, xc)
    else 0
    c_f0 <- sum(x0 * vapply(units, function(u) if (is.null(u$cost)) 0 else u$cost, 
        numeric(1))) + (if (j_floor > 0) 
        j_floor * joint$cost
    else 0) + j_fixed(j_floor) + activation(matrix(x0, ncol = 1L))
    r00 <- if (j_floor > 0) 
        ranges[[1L]]
    else r0
    zero <- list(x = x0, q = rep(1, n_u), rho = r00$hi, tau = .emc_alloc_tau(r00$hi), 
        cost = C_rest + c_p0 + c_f0, objective = .emc_alloc_tau(r00$hi) * 
            (C_rest + c_p0 + c_f0), joint = j_floor)
    lo <- min(vapply(ranges, `[[`, numeric(1), "lo"))
    held <- .emc_alloc_particle_held(units) || any(vapply(units, 
        function(x) identical(x$kind, "subject") && !is.null(x$h) && 
            isTRUE(x$r_fixed > 0), logical(1)))
    best <- if (held) 
        list(objective = Inf)
    else zero
    if (r0$hi > lo + 9.9999999999999998e-13 || c_p0 > 0 || held) {
        tau_grid <- exp(seq(log(.emc_alloc_tau(lo)), log(.emc_alloc_tau(r0$hi)), 
            length.out = 400L))
        lv <- unlist(lapply(ranges, function(r) r$levels[r$levels >= 
            r$lo & r$levels <= r$hi]))
        targets <- sort(unique(c((tau_grid - 1)/(tau_grid + 1), 
            lv)))
        tau_t <- .emc_alloc_tau(targets)
        still <- setdiff(seq_len(n_u), moved_idx)
        need_of <- function(us, movable, layout = NULL) .emc_alloc_units_need(us, 
            movable, targets, max_rounds, layout)
        moved_layout <- .emc_alloc_need_layout(units[moved_idx], 
            max_rounds)
        need <- vector("list", n_u)
        need[still] <- need_of(units[still], r0$movable[still])
        base <- C_rest + Reduce(`+`, lapply(need[still], `[[`, 
            "cost"), numeric(length(targets)))
        act_rows <- .emc_alloc_act_rows(units)
        act_add <- list(shear = 0, hmc = r_min[act_rows$hmc])
        act_fixed <- list(shear = shear_fixed, hmc = hmc_fixed)
        act_part <- function(nd, f) {
            rows <- act_rows[[f]]
            if (act_fixed[[f]] <= 0 || !length(rows)) 
                return(NULL)
            act_fixed[[f]] * emc_alloc_act_term_list_cpp(lapply(nd[rows], 
                `[[`, "x"), rep_len(as.numeric(act_add[[f]]), 
                length(rows)))
        }
        act_held <- lapply(c(shear = "shear", hmc = "hmc"), function(f) if (any(act_rows[[f]] %in% 
            moved_idx)) 
            NULL
        else list(act_part(need, f)))
        best_k <- NULL
        for (c in seq_along(xcs)) {
            nm <- need_of(cand[[c]][moved_idx], ranges[[c]]$movable[moved_idx], 
                moved_layout)
            need_c <- need
            need_c[moved_idx] <- nm
            act <- lapply(c(shear = "shear", hmc = "hmc"), function(f) if (is.null(act_held[[f]])) 
                act_part(need_c, f)
            else act_held[[f]][[1L]])
            cost <- base + xcs[c] * (if (xcs[c] > 0) 
                joint$cost
            else 0) + j_fixed(xcs[c]) + Reduce(`+`, lapply(nm, 
                `[[`, "cost"), numeric(length(targets))) + .emc_alloc_act_sum(act$shear, 
                act$hmc, length(targets))
            obj <- tau_t * cost
            k <- which.min(obj)
            if (!length(k) || !is.finite(obj[k]) || obj[k] >= 
                best$objective * (1 - 1.0000000000000001e-09)) 
                next
            best_nm <- nm
            best_k <- k
            best <- list(x = NULL, q = NULL, rho = targets[k], 
                tau = tau_t[k], cost = cost[k], objective = obj[k], 
                joint = xcs[c])
        }
        if (!is.null(best_k)) {
            need[moved_idx] <- best_nm
            best$x <- vapply(need, function(z) z$x[best_k], numeric(1))
            best$q <- vapply(need, function(z) z$q[best_k], numeric(1))
        }
    }
    if (is.null(best$x) || !held && best$objective > .EMC_ALLOC_MIN_GAIN * 
        zero$objective) 
        best <- zero
    best$objective0 <- zero$objective
    best$tau0 <- zero$tau
    best$tau0_all <- .emc_alloc_tau(r0$hi_all)
    best$pred <- .emc_alloc_predict(units, best, joint, max_rounds)
    best
}
.emc_alloc_solve_units <-
function (u) 
{
    fixed_shear <- if (is.null(u$shear_fixed)) 
        0
    else u$shear_fixed
    if (fixed_shear > 0 && any(vapply(u$units, function(x) identical(x$kind, 
        "shear"), logical(1)))) {
        movable <- .emc_alloc_movable(u$units, .EMC_ALLOC_MAX_ROUNDS, 
            u$joint)
        no_shear <- u
        explored <- FALSE
        no_shear$units <- lapply(u$units, function(x) {
            if (!identical(x$kind, "shear")) 
                return(x)
            if (!is.null(x$floor) && x$floor > 0) {
                explored <<- TRUE
                x$cap <- x$floor
                return(x)
            }
            list(kind = "fixed", rho = x$rho, cost = 0, kappa = x$kappa)
        })
        no_shear$shear_fixed <- 0
        no_shear$shear_act <- if (explored) 
            fixed_shear
        else 0
        no_shear$movable <- movable
        without <- .emc_alloc_solve_units(no_shear)
        on <- u
        on$shear_fixed <- 0
        on$shear_act <- fixed_shear
        on$movable <- movable
        with <- .emc_alloc_solve_units(on)
        if (isTRUE(with$objective < without$objective && with$objective <= 
            .EMC_ALLOC_MIN_GAIN * without$objective0)) {
            with$objective0 <- without$objective0
            with$tau0 <- without$tau0
            return(with)
        }
        return(without)
    }
    off <- lapply(u$units, function(x) {
        if (identical(x$kind, "hmc")) 
            x <- .emc_alloc_as_subject(x)
        if (identical(x$kind, "subject")) 
            x["h"] <- list(NULL)
        x
    })
    movable <- if (!is.null(u$movable)) 
        u$movable
    else .emc_alloc_movable(u$units, .EMC_ALLOC_MAX_ROUNDS, u$joint)
    has_h <- any(vapply(u$units, function(x) !is.null(x$h) || 
        identical(x$kind, "hmc"), logical(1)))
    regimes <- list()
    if (has_h) {
        levels0 <- .emc_alloc_range(u$units, .EMC_ALLOC_MAX_ROUNDS)$levels
        for (r in .EMC_ALLOC_HMC_REGIMES) {
            reg <- .emc_alloc_regime(u, r)
            if (is.null(reg)) 
                break
            reg$r <- r
            regimes[[length(regimes) + 1L]] <- reg
            movable <- movable | (.emc_alloc_range(reg$units, 
                .EMC_ALLOC_MAX_ROUNDS)$levels < levels0 - 1.0000000000000001e-09)
        }
    }
    shear_act <- if (is.null(u$shear_act)) 
        0
    else u$shear_act
    hmc_fixed <- if (is.null(u$hmc_fixed)) 
        0
    else u$hmc_fixed
    if (has_h && !is.null(u$g_floor)) {
        ratio <- .emc_alloc_dev_ratio(u$units, rep(1, length(u$units)), 
            numeric(length(u$units)))
        if (is.finite(ratio)) {
            n_g <- length(u$g_floor)
            off[seq_len(n_g)] <- .emc_alloc_credit(off[seq_len(n_g)], 
                u$g_floor, ratio)
        }
    }
    without <- .emc_alloc_solve(off, u$C_rest, joint = u$joint, 
        movable = movable, shear_fixed = shear_act)
    without$regime <- 0
    if (!has_h) 
        return(without)
    settle <- function(sol, used) {
        if (is.null(u$g_floor)) 
            return(sol)
        .emc_alloc_settle(u, sol, used$units, used$credit, used$C_rest, 
            movable, shear_act, hmc_fixed)
    }
    plans <- list()
    sol <- .emc_alloc_solve(u$units, u$C_rest, joint = u$joint, 
        movable = movable, shear_fixed = shear_act, hmc_fixed = hmc_fixed)
    sol$regime <- 0
    plans[["0"]] <- settle(sol, list(units = u$units, credit = 1, 
        C_rest = u$C_rest))
    for (reg in regimes) {
        sol <- .emc_alloc_solve(reg$units, u$C_rest + reg$cost, 
            joint = u$joint, movable = movable, shear_fixed = shear_act, 
            hmc_fixed = hmc_fixed)
        sol$regime <- reg$r
        plans[[as.character(reg$r)]] <- settle(sol, list(units = reg$units, 
            credit = reg$credit, C_rest = u$C_rest + reg$cost))
    }
    held <- .emc_alloc_particle_held(u$units)
    gain_ok <- function(x) isTRUE(x$objective <= .EMC_ALLOC_MIN_GAIN * 
        without$objective0)
    options <- c(none = if (held) Inf else without$objective, 
        vapply(plans, function(x) if (held || gain_ok(x)) x$objective else Inf, 
            numeric(1)))
    choice <- u$regime_fixed
    if (is.null(choice) || !(choice %in% names(options))) 
        choice <- names(options)[which.min(options)]
    if (identical(choice, "none") || !is.finite(options[[choice]]) && 
        !held) {
        without$options <- options
        without$choice <- "none"
        return(without)
    }
    out <- plans[[choice]]
    out$objective0 <- without$objective0
    out$tau0 <- without$tau0
    out$options <- options
    out$choice <- choice
    out
}
.EMC_ALLOC_STAGE_DRAWS <-
1000L
.emc_alloc_stage_tau <-
function (emc, Dinv, stage) 
{
    p <- emc[[1L]]$n_pars
    X <- lapply(emc, function(x) {
        i <- x$samples$idx
        lab <- if (is.null(stage)) 
            x$samples$stage[i]
        else stage
        k <- which(x$samples$stage[seq_len(i)] == lab)
        k <- utils::tail(k, .EMC_ALLOC_STAGE_DRAWS)
        if (length(k) < .EMC_ALLOC_MIN_PAIRS) 
            return(NULL)
        cbind(t(Dinv %*% x$samples$theta_mu[, k, drop = FALSE]), 
            t(.emc_alloc_lsd_draws(x$samples$theta_var, k)))
    })
    if (any(vapply(X, is.null, logical(1)))) 
        return(rep(NA_real_, 2L * p))
    X <- lapply(X, function(x) x[apply(is.finite(x), 1L, all), 
        , drop = FALSE])
    .emc_alloc_geyer_tau(X)
}
.emc_alloc_subject_base <-
function (u) 
if (is.null(u$p)) u$rho else u$rho * u$p
.emc_alloc_tau <-
function (rho) 
{
    rho <- pmin(rho, 1 - 9.9999999999999998e-13)
    (1 + rho)/(1 - rho)
}
.emc_alloc_unit_need <-
function (u, rs, max_rounds) 
{
    n_t <- length(rs)
    if (identical(u$kind, "fixed")) 
        return(list(x = numeric(n_t), q = rep(1, n_t), cost = ifelse(.emc_alloc_eff(max(u$rho), 
            u$kappa) > rs, Inf, 0)))
    if (identical(u$kind, "shear")) {
        rs <- .emc_alloc_eff_inv(rs, u$kappa)
        r <- u$rho
        g <- u$g
        cap <- if (is.null(u$cap)) 
            max_rounds
        else u$cap
        floor <- if (is.null(u$floor)) 
            0
        else u$floor
        x <- rep(floor, n_t)
        over <- r > rs
        if (any(over)) {
            xo <- .emc_alloc_factor_inverse(g, if (is.null(u$a)) 
                g
            else u$a, rs[over]/r)
            xo <- pmax(xo, floor)
            xo[!is.finite(xo) | xo > cap + 1.0000000000000001e-09] <- Inf
            x[over] <- xo
        }
        return(list(x = x, q = rep(1, n_t), cost = ifelse(is.finite(x), 
            x * u$cost, Inf)))
    }
    rho <- u$rho
    level <- .emc_alloc_subject_base(u)
    qs <- if (is.null(u$p) || is.null(u$h)) 
        1
    else if (!is.null(u$q_fixed)) 
        u$q_fixed
    else .EMC_ALLOC_PARTICLE_RATES
    r_min <- if (is.null(u$h) || is.null(u$r_min)) 
        0
    else u$r_min
    r_need <- emc_alloc_subject_need_cpp(as.numeric(rho), if (is.null(u$p)) 
        rep(1, length(rho))
    else as.numeric(u$p), if (is.null(u$h)) 
        NULL
    else as.numeric(u$h), as.numeric(level), as.numeric(qs), 
        as.numeric(rs), as.numeric(r_min))
    r_need[r_need > 1 + 1.0000000000000001e-09] <- Inf
    fin <- is.finite(r_need)
    r_need[fin] <- pmin(r_need[fin], 1)
    if (!is.null(u$h) && !is.null(u$r_fixed)) 
        r_need[] <- ifelse(r_need <= u$r_fixed + 1.0000000000000001e-09, 
            u$r_fixed, Inf)
    fin <- is.finite(r_need)
    c_h <- if (is.null(u$c_h)) 
        0
    else u$c_h
    c_p <- if (is.null(u$c_p)) 
        0
    else u$c_p
    cost <- qs * c_p + (r_need - r_min) * c_h
    cost[!fin] <- Inf
    k <- max.col(t(-cost), ties.method = "first")
    pick <- cbind(k, seq_len(n_t))
    list(x = r_need[pick] - r_min, q = qs[k], cost = cost[pick])
}
.emc_alloc_unit_still <-
function (u, n_t) 
{
    c_p <- if (identical(u$kind, "subject") && !is.null(u$p) && 
        !is.null(u$c_p)) 
        u$c_p
    else 0
    fl <- if (identical(u$kind, "shear") && !is.null(u$floor)) 
        u$floor
    else 0
    list(x = rep(fl, n_t), q = rep(1, n_t), cost = rep(c_p + 
        fl * (if (is.null(u$cost)) 0 else u$cost), n_t))
}
.emc_alloc_units <-
function (ms, hmc_state, shear_on, optimistic, B = NULL) 
{
    n <- nrow(ms$s_rho)
    bounds <- function(n_ev, iters) {
        if (is.null(B) || !is.finite(B) || B <= 0) 
            return(list(floor = 0, cap = Inf))
        n0 <- .EMC_ALLOC_ROUND_EVIDENCE
        floor <- max(0, n0 - n_ev)/B
        span <- if (length(iters) && is.finite(iters) && iters > 
            0) 
            iters
        else B
        list(floor = floor, cap = max(floor, 2 * n_ev/span))
    }
    C_rest <- if (is.finite(ms$C_rest)) 
        ms$C_rest
    else 1
    C_iter <- C_rest + sum(ms$c_part[is.finite(ms$c_part)])
    c_round <- if (is.finite(ms$c_round)) 
        ms$c_round
    else 0.10000000000000001 * C_iter
    kappa <- if (is.null(ms$g_kappa)) 
        rep(1, length(ms$g_rho))
    else ms$g_kappa
    units <- lapply(seq_along(ms$g_rho), function(k) if (shear_on && 
        is.finite(ms$g_fac[k])) {
        b <- bounds(if (is.null(ms$g_n)) 
            0
        else ms$g_n[k], if (is.null(ms$g_iters)) 
            0
        else ms$g_iters[k])
        list(kind = "shear", rho = ms$g_rho[k], g = ms$g_fac[k], 
            a = if (is.null(ms$g_sat)) ms$g_fac[k] else ms$g_sat[k], 
            cost = c_round, floor = min(b$floor, .EMC_ALLOC_MAX_ROUNDS), 
            cap = min(b$cap, .EMC_ALLOC_MAX_ROUNDS), kappa = kappa[k])
    }
    else list(kind = "fixed", rho = ms$g_rho[k], cost = 0, kappa = kappa[k]))
    for (s in seq_len(n)) {
        sj <- if (is.null(hmc_state)) 
            NULL
        else hmc_state$subjects[[s]]
        priced <- !is.null(ms$c_part) && is.finite(ms$c_part[s])
        u <- list(kind = "subject", rho = ms$s_rho[s, ], p = NULL, 
            c_p = 0, h = NULL, c_h = 0)
        if (isTRUE(sj$eligible)) {
            ev <- hmc_state$evidence[[s]]
            h <- .emc_alloc_hmc_factor(ev, ms$W[[s]], ms$s_var[s, 
                ], .EMC_HMC_MIN_EVIDENCE)
            if (is.null(h) && optimistic) 
                h <- numeric(length(ms$s_rho[s, ]))
            else if (!is.null(h)) 
                u$h_n <- ev$n
            if (!is.null(h)) {
                c_s <- if (is.finite(ms$c_hmc[s])) 
                  ms$c_hmc[s]
                else if (!is.null(ev) && length(ev$arows) && 
                  isTRUE(is.finite(ms$hmc_row_cpu)) && !is.null(ms$data_rows)) 
                  mean(ev$arows) * ms$data_rows[s] * ms$hmc_row_cpu
                else if (!is.null(ev) && ev$n > 0 && ev$cpu > 
                  0) 
                  ev$cpu/ev$n
                else 3 * C_iter/n
                u$h <- h
                u$c_h <- c_s
                if (!is.null(ms$d_var) && all(is.finite(ms$d_var[s, 
                  ]))) {
                  u$h_dev <- .emc_alloc_hmc_factor(ev, ms$W[[s]], 
                    ms$d_var[s, ], .EMC_HMC_MIN_EVIDENCE)
                  if (is.null(u$h_dev)) 
                    u$h_dev <- h
                }
            }
        }
        if (!is.null(ms$d_rho) && all(is.finite(ms$d_rho[s, ]))) {
            u$d_rho <- ms$d_rho[s, ]
            if (!is.null(ms$p_dev) && all(is.finite(ms$p_dev[s, 
                ]))) 
                u$p_dev <- ms$p_dev[s, ]
            u$q_old <- if (is.null(ms$q_old)) 
                1
            else ms$q_old[s]
            u$r_old <- if (is.null(ms$r_old)) 
                0
            else ms$r_old[s]
        }
        if (priced) {
            u$rho <- ms$t_rho[s, ]
            u$p <- ms$p_fac[s, ]
            u$c_p <- ms$c_part[s]
        }
        units[[length(units) + 1L]] <- u
    }
    joint <- if (isTRUE(ms$coll)) {
        b <- bounds(if (is.null(ms$c_n)) 
            0
        else ms$c_n, if (is.null(ms$c_iters)) 
            0
        else ms$c_iters)
        list(idx = seq_along(ms$g_rho), g = ms$c_fac, a = if (is.null(ms$c_sat)) ms$c_fac else ms$c_sat, 
            cost = if (is.finite(ms$c_coll)) ms$c_coll else 0.10000000000000001 * 
                C_iter, fixed = if (isTRUE(is.finite(ms$coll_fixed))) ms$coll_fixed else 0, 
            floor = b$floor, cap = min(b$cap, max(.EMC_ALLOC_MAX_ROUNDS, 
                2 * ms$c_rate)))
    }
    n0 <- .EMC_ALLOC_ROUND_EVIDENCE
    hs <- lapply(units[-seq_along(ms$g_rho)], function(x) x$h)
    ns <- vapply(units[-seq_along(ms$g_rho)], function(x) if (is.null(x$h_n)) 
        0
    else x$h_n, numeric(1))
    with_h <- which(ns > 0 & !vapply(hs, is.null, logical(1)))
    if (length(with_h)) {
        h_ref <- Reduce(`+`, lapply(with_h, function(i) ns[i] * 
            hs[[i]]))/sum(ns[with_h])
        hd <- lapply(units[-seq_along(ms$g_rho)], function(x) x$h_dev)
        with_d <- with_h[!vapply(hd[with_h], is.null, logical(1))]
        hd_ref <- if (length(with_d)) 
            Reduce(`+`, lapply(with_d, function(i) ns[i] * hd[[i]]))/sum(ns[with_d])
        for (i in with_h) {
            k <- length(ms$g_rho) + i
            units[[k]]$h <- (ns[i] * hs[[i]] + n0 * h_ref)/(ns[i] + 
                n0)
            if (!is.null(units[[k]]$h_dev) && !is.null(hd_ref)) 
                units[[k]]$h_dev <- (ns[i] * units[[k]]$h_dev + 
                  n0 * hd_ref)/(ns[i] + n0)
        }
    }
    list(units = units, C_rest = C_rest, c_round = c_round, joint = joint, 
        hmc_fixed = if (isTRUE(is.finite(ms$hmc_fixed))) ms$hmc_fixed else 0, 
        shear_fixed = if (isTRUE(is.finite(ms$shear_fixed))) ms$shear_fixed else 0, 
        r_bar = ms$r_bar, g_floor = ms$g_floor)
}
.emc_alloc_units_need <-
function (us, movable, rs, max_rounds, layout = NULL) 
{
    if (is.null(layout)) 
        layout <- .emc_alloc_need_layout(us, max_rounds)
    out <- vector("list", length(us))
    shear <- movable & layout$shear
    subj <- movable & layout$kind == "subject"
    fixed <- movable & layout$kind == "fixed"
    if (any(fixed)) 
        out[fixed] <- emc_alloc_fixed_need_cpp(vapply(us[fixed], 
            function(u) .emc_alloc_eff(max(u$rho), u$kappa), 
            numeric(1)), as.numeric(rs))
    for (i in which(!subj & !shear & !fixed)) out[[i]] <- if (movable[i]) 
        .emc_alloc_unit_need(us[[i]], rs = rs, max_rounds = max_rounds)
    else .emc_alloc_unit_still(us[[i]], length(rs))
    sh <- which(shear)
    if (length(sh)) {
        kap <- layout$kappa[match(sh, which(layout$shear))]
        plain <- sh[kap == 1]
        price <- function(i, r) {
            k <- match(i, which(layout$shear))
            emc_alloc_shears_need_cpp(as.numeric(unlist(lapply(us[i], 
                `[[`, "rho"))), layout$g[k], layout$a[k], layout$cost[k], 
                layout$cap[k], layout$floor[k], as.numeric(r))
        }
        if (length(plain)) 
            out[plain] <- price(plain, rs)
        for (i in sh[kap != 1]) out[i] <- price(i, .emc_alloc_eff_inv(rs, 
            us[[i]]$kappa))
    }
    idx <- which(subj)
    if (!length(idx)) 
        return(out)
    num <- function(x, d) if (is.null(x)) 
        d
    else as.numeric(x)
    rho <- lapply(us[idx], function(u) as.numeric(u$rho))
    k <- lengths(rho)
    has_p <- vapply(us[idx], function(u) !is.null(u$p), logical(1))
    has_h <- vapply(us[idx], function(u) !is.null(u$h), logical(1))
    got <- emc_alloc_subjects_need_cpp(unlist(rho), unlist(lapply(seq_along(idx), 
        function(j) num(us[[idx[j]]]$p, rep(1, k[j])))), unlist(lapply(seq_along(idx), 
        function(j) num(us[[idx[j]]]$h, rep(1, k[j])))), unlist(lapply(us[idx], 
        function(u) as.numeric(.emc_alloc_subject_base(u)))), 
        c(0L, cumsum(k)), has_p, has_h, vapply(us[idx], function(u) num(u$q_fixed, 
            NA_real_), numeric(1)), vapply(us[idx], function(u) num(u$r_fixed, 
            NA_real_), numeric(1)), vapply(us[idx], function(u) num(u$r_min, 
            0), numeric(1)), vapply(us[idx], function(u) num(u$c_p, 
            0), numeric(1)), vapply(us[idx], function(u) num(u$c_h, 
            0), numeric(1)), .EMC_ALLOC_PARTICLE_RATES, as.numeric(rs))
    for (j in seq_along(idx)) out[[idx[j]]] <- list(x = got$x[j, 
        ], q = got$q[j, ], cost = got$cost[j, ])
    out
}
.emc_alloc_wanted <-
function (sampler) 
{
    .emc_shear_active(sampler) || .emc_collapsed_active(sampler) || 
        isTRUE(.emc_hmc_support(sampler)$ok)
}
.emc_alloc_whole_rho <-
function (emc, draws, Dinv, var) 
{
    p <- emc[[1L]]$n_pars
    sq <- numeric(2L * p)
    n <- 0L
    for (i in seq_along(emc)) {
        k <- draws[[i]]
        if (length(k) < 2L) 
            next
        s <- emc[[i]]$samples
        x <- rbind(Dinv %*% s$theta_mu[, k, drop = FALSE], .emc_alloc_lsd_draws(s$theta_var, 
            k))
        d <- x[, -1L, drop = FALSE] - x[, -ncol(x), drop = FALSE]
        d <- d[, apply(is.finite(d), 2L, all), drop = FALSE]
        sq <- sq + rowSums(d^2)
        n <- n + ncol(d)
    }
    if (n == 0L) 
        return(rep(NA_real_, 2L * p))
    .emc_alloc_rho(sq/n, var)
}
.emc_alloc_window <-
function (sampler) 
{
    w <- sampler$alloc$window
    if (is.null(w)) 
        w <- .emc_alloc_empty_window(sampler$n_subjects, sampler$n_pars)
    w
}
