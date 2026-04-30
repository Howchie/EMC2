source("WorkingTests/Helper_Functions.R")

.capacity_bounds_or <- function(times, rt, CR, denom_H, ratio = TRUE) {
  F1 <- .safe_ecdf(rt[[2]][.to_correct_indicator(CR[[2]], length(rt[[2]]))])
  F2 <- .safe_ecdf(rt[[3]][.to_correct_indicator(CR[[3]], length(rt[[3]]))])
  pmax01 <- function(x) pmin(1, pmax(0, x))
  
  F_upper <- pmax01(F1(times) + F2(times))
  F_lower <- pmax01(pmax(F1(times), F2(times)))
  H_upper <- -log1p(-F_upper)
  H_lower <- -log1p(-F_lower)
  
  if (ratio) {
    Ct_upper <- H_upper / denom_H(times)
    Ct_lower <- H_lower / denom_H(times)
  } else {
    Ct_upper <- H_upper - denom_H(times)
    Ct_lower <- H_lower - denom_H(times)
  }
  
  Ct_upper[is.nan(Ct_upper) | is.infinite(Ct_upper)] <- NA_real_
  Ct_lower[is.nan(Ct_lower) | is.infinite(Ct_lower)] <- NA_real_
  
  list(
    Ct_upper = approxfun(times, Ct_upper),
    Ct_lower = approxfun(times, Ct_lower),
    H_upper = approxfun(times, H_upper),
    H_lower = approxfun(times, H_lower),
    F_upper = approxfun(times, F_upper),
    F_lower = approxfun(times, F_lower)
  )
}

.capacity_bounds_and <- function(times, rt, CR, denom_K, ratio = TRUE) {
  F1 <- .safe_ecdf(rt[[2]][.to_correct_indicator(CR[[2]], length(rt[[2]]))])
  F2 <- .safe_ecdf(rt[[3]][.to_correct_indicator(CR[[3]], length(rt[[3]]))])
  pmax01 <- function(x) pmin(1, pmax(0, x))
  
  F_upper <- pmax01(pmin(F1(times), F2(times)))       # Miller upper bound in CDF space
  F_lower <- pmax01(F1(times) + F2(times) - 1)        # Grice lower bound in CDF space
  # In this codebase, reverse cumulative hazard is represented as log(F(t)) (<= 0).
  K_from_upperF <- log(pmax(F_upper, .Machine$double.eps))
  K_from_lowerF <- log(pmax(F_lower, .Machine$double.eps))
  
  if (ratio) {
    # Handle denominator sign/near-zero robustly by ordering numeric candidates.
    ratio_a <- K_from_lowerF / denom_K(times)
    ratio_b <- K_from_upperF / denom_K(times)
    Ct_upper <- pmax(ratio_a, ratio_b)
    Ct_lower <- pmin(ratio_a, ratio_b)
  } else {
    # C(t) = K_UCIP(t) - K_AB(t), so upper C uses smaller admissible K_AB.
    Ct_upper <- denom_K(times) - K_from_lowerF
    Ct_lower <- denom_K(times) - K_from_upperF
  }
  
  Ct_upper[is.nan(Ct_upper) | is.infinite(Ct_upper)] <- NA_real_
  Ct_lower[is.nan(Ct_lower) | is.infinite(Ct_lower)] <- NA_real_
  
  list(
    Ct_upper = approxfun(times, Ct_upper),
    Ct_lower = approxfun(times, Ct_lower),
    K_upper = approxfun(times, K_from_upperF),
    K_lower = approxfun(times, K_from_lowerF),
    F_upper = approxfun(times, F_upper),
    F_lower = approxfun(times, F_lower)
  )
}

capacity.or <- function(rt, CR=NULL, ratio=TRUE) {
  if ( is.null(CR) | (length(CR) != length(rt)) ) {
    CR <- lapply(rt, function(x) rep(1, length(x)))
  } 
  # Bolt: Use unlist(..., use.names=FALSE) instead of c(..., recursive=TRUE) for faster array flattening
  times <- sort(unique(unlist(rt, use.names=FALSE)))
  ncond <- length(rt) - 1 
  
  # Find Nelson-Aalen Cumulative Hazard Estimates
  numer <- estimateNAH(rt=rt[[1]], CR=CR[[1]])
  denom <- estimateUCIPor(rt=rt[1+(1:ncond)], CR=CR[1+(1:ncond)])
  
  
  if (ratio) {
    C.or <- numer$H(times) / denom$H(times)
    
    C.or[is.nan(C.or)] <- NA
    C.or[is.infinite(C.or)] <- NA
    C.or <- approxfun(times, C.or)
    bounds <- .capacity_bounds_or(
      times = times,
      rt = rt,
      CR = CR,
      denom_H = denom$H,
      ratio = TRUE
    )
    return(list(Ct = C.or, Ct_upper = bounds$Ct_upper, Ct_lower = bounds$Ct_lower, times=times))
  } else {
    C.or <- numer$H(times) - denom$H(times)
    C.or <- approxfun(c(0,times), c(0,C.or))
    Var.or <- numer$Var(times) + denom$Var(times)
    Var.or <- approxfun(c(0,times), c(0,Var.or))
    bounds <- .capacity_bounds_or(
      times = times,
      rt = rt,
      CR = CR,
      denom_H = denom$H,
      ratio = FALSE
    )
    return(list(Ct = C.or, Var = Var.or, Ct_upper = bounds$Ct_upper, Ct_lower = bounds$Ct_lower, times=times))
  }
}

capacity.and <- function(rt, CR=NULL, ratio=TRUE) {
  if ( is.null(CR) | (length(CR) != length(rt)) ) {
    CR <- lapply(rt, function(x) rep(1, length(x)))
  } 
  # Bolt: Use unlist(..., use.names=FALSE) instead of c(..., recursive=TRUE) for faster array flattening
  times <- sort(unique(unlist(rt, use.names=FALSE)))
  
  ncond <- length(rt) - 1 
  
  
  # Find Nelson-Aalen Reverse Cumulative Hazard Estimates
  denom <- estimateNAK(rt[[1]], CR[[1]])
  numer <- estimateUCIPand(rt=rt[1+(1:ncond)], CR=CR[1+(1:ncond)])
  
  # Calculate the and capacity coefficient
  if (ratio) {
    C.and <- numer$K(times) / denom$K(times)
    C.and[is.nan(C.and)] <- NA
    C.and[is.infinite(C.and)] <- NA
    C.and <- approxfun(times, C.and)
    bounds <- .capacity_bounds_and(
      times = times,
      rt = rt,
      CR = CR,
      denom_K = denom$K,
      ratio = TRUE
    )
    return(list(Ct = C.and, Ct_upper = bounds$Ct_upper, Ct_lower = bounds$Ct_lower, times=times))
  } else {
    C.and <- denom$K(times) - numer$K(times) 
    C.and <- approxfun(c(times,Inf), c(C.and,0))
    Var.and <- numer$Var(times) + denom$Var(times)
    Var.and <- approxfun(c(times,Inf), c(Var.and,0))
    bounds <- .capacity_bounds_and(
      times = times,
      rt = rt,
      CR = CR,
      denom_K = denom$K,
      ratio = FALSE
    )
    return(list(Ct = C.and, Var = Var.and, Ct_upper = bounds$Ct_upper, Ct_lower = bounds$Ct_lower, times=times))
  }
}

capacity.id <- function(dt.rt, nt.rt, st.rts, dt.cr=NULL, nt.cr=NULL, st.crs=NULL, ratio=TRUE) {
  n_single <- length(st.rts)
  
  if ( is.null(dt.cr) ) { 
    dt.cr <- rep(1, length(dt.rt))
  } 
  if ( is.null(nt.cr) ) { 
    nt.cr <- rep(1, length(nt.rt))
  } 
  if ( is.null(st.crs) | (length(st.crs) != n_single) ) {
    st.crs <- lapply(st.rts, function(x) rep(1, length(x)))
  } 
  # Bolt: Use unlist(..., use.names=FALSE) instead of c(..., recursive=TRUE) for faster array flattening
  times <- sort(unique(c(dt.rt, nt.rt, unlist(st.rts, use.names=FALSE))))
  
  
  # Find Nelson-Aalen Reverse Cumulative Hazard Estimates
  NAK.dual <- estimateNAK(dt.rt, dt.cr) 
  NAK.no <- estimateNAK(nt.rt, nt.cr)
  NAK.single <- vector("list", n_single)
  for ( i in 1:n_single) { 
    NAK.single[[i]] <- estimateNAK(st.rts[[i]], st.crs[[i]])
  }
  
  # Calculate the and capacity coefficient
  denom <- NAK.dual$K(times) + NAK.no$K(times)
  if (!ratio) {
    Var.id <- NAK.dual$Var(times) + NAK.no$Var(times)
  }
  
  # Bolt: Use rowSums and matrix to push accumulation into C code for better performance
  nt_id <- length(times)
  numer <- rowSums(matrix(vapply(NAK.single, function(nak) nak$K(times), numeric(nt_id)), nrow = nt_id))
  if (!ratio) {
    Var.id <- Var.id + rowSums(matrix(vapply(NAK.single, function(nak) nak$Var(times), numeric(nt_id)), nrow = nt_id))
  }
  
  if (ratio) {
    C.id <- numer / denom
    C.id[is.nan(C.id)] <- NA
    C.id[is.infinite(C.id)] <- NA
    C.id <- approxfun(times, C.id)
    return( list(Ct=C.id, times=times) )
  } else {
    C.id <- denom - numer
    C.id <- approxfun(c(times,Inf), c(C.id,0))
    Var.id <- approxfun(c(times,Inf), c(Var.id,0))
    return( list(Ct=C.id, Var=Var.id, times=times) )
  }
}

capacity.altieri <- function(rt, CR=NULL, ratio=TRUE) {
  if ( is.null(CR) | (length(CR) != length(rt)) ) {
    CR <- lapply(rt, function(x) rep(1, length(x)))
  } 
  # Bolt: Use unlist(..., use.names=FALSE) instead of c(..., recursive=TRUE) for faster array flattening
  times <- sort(unique(unlist(rt, use.names=FALSE)))
  ncond <- length(rt) - 2 
  
  # Find Nelson-Aalen Cumulative Hazard Estimates
  numer <- estimateNAH(rt=rt[[1]], CR=CR[[1]])
  denom1 <- estimateUCIPor(rt=rt[1+(1:ncond)], CR=CR[1+(1:ncond)])
  denom2 <- estimateNAH(rt=rt[[4]], CR=CR[[4]])
  
  if (ratio) {
    C.altieri1 <- numer$H(times) / denom1$H(times)
    
    C.altieri1[is.nan(C.altieri1)] <- NA
    C.altieri1[is.infinite(C.altieri1)] <- NA
    C.altieri1 <- approxfun(times, C.altieri1)
    
    C.altieri2 <- numer$H(times) / denom2$H(times)
    
    C.altieri2[is.nan(C.altieri2)] <- NA
    C.altieri2[is.infinite(C.altieri2)] <- NA
    C.altieri2 <- approxfun(times, C.altieri2)

    return( list(Ct_a1=C.altieri1,Ct_a2=C.altieri2) )
  } else {
    C.altieri1 <- numer$H(times) - denom1$H(times)
    C.altieri1 <- approxfun(c(0,times), c(0,C.altieri1))
    C.altieri2 <- numer$H(times) - denom2$H(times)
    C.altieri2 <- approxfun(c(0,times), c(0,C.altieri2))
    return( list(Ct_a1=C.altieri1,Ct_a2=C.altieri2, times=times) )
  }
}

.ensure_hazard_rcpp <- local({
  compiled <- FALSE
  function() {
    if (compiled &&
        exists("nah_nak_eval_rcpp", mode = "function") &&
        exists("ucip_eval_rcpp", mode = "function")) {
      return(TRUE)
    }
    if (!requireNamespace("Rcpp", quietly = TRUE)) {
      return(FALSE)
    }
    ok <- tryCatch({
      Rcpp::sourceCpp(code = '
        #include <Rcpp.h>
        #include <vector>
        #include <algorithm>
        #include <cmath>
        using namespace Rcpp;

        struct RTCR {
          double rt;
          bool cr;
        };

        static inline NumericMatrix eval_one_hazard(
            NumericVector rt,
            LogicalVector CR,
            NumericVector query,
            bool reverse_hazard
        ) {
          const int n = rt.size();
          const int nq = query.size();
          NumericMatrix out(nq, 2);
          if (n == 0 || nq == 0) return out;

          std::vector<RTCR> v;
          v.reserve(n);
          for (int i = 0; i < n; ++i) {
            if (NumericVector::is_na(rt[i]) || !R_finite(rt[i])) continue;
            bool cr_i = true;
            if (i < CR.size()) {
              if (CR[i] == NA_LOGICAL) cr_i = false;
              else cr_i = (CR[i] == TRUE);
            }
            RTCR x;
            x.rt = rt[i];
            x.cr = cr_i;
            v.push_back(x);
          }

          const int m = (int)v.size();
          if (m == 0) return out;

          std::sort(v.begin(), v.end(), [](const RTCR& a, const RTCR& b) {
            return a.rt < b.rt;
          });

          std::vector<double> event_rt;
          std::vector<double> cum_h;
          std::vector<double> cum_v;
          event_rt.reserve(m);
          cum_h.reserve(m);
          cum_v.reserve(m);

          double h = 0.0;
          double vv = 0.0;
          for (int i = 0; i < m; ++i) {
            if (!v[i].cr) continue;
            const double risk = reverse_hazard ? (double)(i + 1) : (double)(m - i);
            const double inc = 1.0 / risk;
            h += inc;
            vv += inc * inc;
            event_rt.push_back(v[i].rt);
            cum_h.push_back(h);
            cum_v.push_back(vv);
          }

          const int ne = (int)event_rt.size();
          if (ne == 0) return out;

          if (!reverse_hazard) {
            // NAH: H(t) = sum_{x_i <= t} 1/Y_i
            for (int j = 0; j < nq; ++j) {
              const double q = query[j];
              if (!R_finite(q)) continue;
              auto it = std::upper_bound(event_rt.begin(), event_rt.end(), q);
              int idx = (int)(it - event_rt.begin()) - 1;
              if (idx >= 0) {
                out(j, 0) = cum_h[idx];
                out(j, 1) = cum_v[idx];
              }
            }
          } else {
            // NAK: K(t) = -sum_{x_i >= t} 1/G_i ; Var uses +sum_{x_i >= t} 1/G_i^2
            std::vector<double> tail_h(ne, 0.0), tail_v(ne, 0.0);
            double th = 0.0, tv = 0.0;
            for (int i = ne - 1; i >= 0; --i) {
              const double inc_h = (i == 0 ? cum_h[0] : cum_h[i] - cum_h[i - 1]);
              const double inc_v = (i == 0 ? cum_v[0] : cum_v[i] - cum_v[i - 1]);
              th += inc_h;
              tv += inc_v;
              tail_h[i] = -th;
              tail_v[i] = tv;
            }
            for (int j = 0; j < nq; ++j) {
              const double q = query[j];
              if (!R_finite(q)) continue;
              auto it = std::lower_bound(event_rt.begin(), event_rt.end(), q);
              int idx = (int)(it - event_rt.begin());
              if (idx < ne) {
                out(j, 0) = tail_h[idx];
                out(j, 1) = tail_v[idx];
              }
            }
          }
          return out;
        }

        // [[Rcpp::export]]
        NumericMatrix nah_nak_eval_rcpp(
            NumericVector rt,
            LogicalVector CR,
            NumericVector query,
            bool reverse_hazard
        ) {
          return eval_one_hazard(rt, CR, query, reverse_hazard);
        }

        // [[Rcpp::export]]
        NumericMatrix ucip_eval_rcpp(
            List RT_list,
            List CR_list,
            NumericVector query,
            bool reverse_hazard
        ) {
          const int nq = query.size();
          NumericMatrix out(nq, 2);
          const int ncond = RT_list.size();
          for (int i = 0; i < ncond; ++i) {
            NumericVector RTi = as<NumericVector>(RT_list[i]);
            LogicalVector CRi;
            if (i < CR_list.size() && !Rf_isNull(CR_list[i])) {
              CRi = as<LogicalVector>(CR_list[i]);
            } else {
              CRi = LogicalVector(RTi.size(), true);
            }
            NumericMatrix one = eval_one_hazard(RTi, CRi, query, reverse_hazard);
            for (int j = 0; j < nq; ++j) {
              out(j, 0) += one(j, 0);
              out(j, 1) += one(j, 1);
            }
          }
          return out;
        }')
      TRUE
    }, error = function(e) FALSE)
    compiled <<- isTRUE(ok) &&
      exists("nah_nak_eval_rcpp", mode = "function") &&
      exists("ucip_eval_rcpp", mode = "function")
    compiled
  }
})

estimateNAH <- function(rt, CR=NULL) {
  nt <- length(rt)
  zero_fun <- function(t) rep(0, length(t))
  if (nt == 0) {
    return(list(H = zero_fun, Var = zero_fun))
  }
  
  if ( is.null(CR) || length(CR) != nt ) {
    CR <- rep(1, nt)
  }
  if (.ensure_hazard_rcpp()) {
    RTv <- as.numeric(rt)
    CRv <- as.logical(CR)
    H_fun <- function(t) nah_nak_eval_rcpp(RTv, CRv, as.numeric(t), FALSE)[, 1]
    V_fun <- function(t) nah_nak_eval_rcpp(RTv, CRv, as.numeric(t), FALSE)[, 2]
    return(list(H = H_fun, Var = V_fun))
  }
  
  RTx <- sort(rt,index.return=TRUE)
  rt <- RTx$x
  CR <- as.logical(CR)[RTx$ix]
  hit_idx <- which(CR)
  if (length(hit_idx) == 0) {
    return(list(H = zero_fun, Var = zero_fun))
  }
  
  Y <- nt:1
  
  H <- stepfun(rt[hit_idx], c(0, cumsum(1 / Y[hit_idx])))
  H.v <- stepfun(rt[hit_idx], c(0, cumsum(1 / (Y[hit_idx] ^ 2))))
  return(list(H=H, Var=H.v))
}

estimateNAK <- function(rt, CR=NULL) {
  nt <- length(rt)
  zero_fun <- function(t) rep(0, length(t))
  if (nt == 0) {
    return(list(K = zero_fun, Var = zero_fun))
  }
  
  if ( is.null(CR) || length(CR) != nt ) {
    CR <- rep(1, nt)
  }
  if (.ensure_hazard_rcpp()) {
    RTv <- as.numeric(rt)
    CRv <- as.logical(CR)
    K_fun <- function(t) nah_nak_eval_rcpp(RTv, CRv, as.numeric(t), TRUE)[, 1]
    V_fun <- function(t) nah_nak_eval_rcpp(RTv, CRv, as.numeric(t), TRUE)[, 2]
    return(list(K = K_fun, Var = V_fun))
  }
  
  RTx <- sort(rt,index.return=TRUE)
  rt <- RTx$x
  CR <- as.logical(CR)[RTx$ix]
  hit_idx <- which(CR)
  if (length(hit_idx) == 0) {
    return(list(K = zero_fun, Var = zero_fun))
  }
  
  G <- seq_len(nt)
  
  K <- stepfun(rt[hit_idx], c(rev(-cumsum(rev(1 / G[hit_idx]))), 0), right = TRUE)
  K.v <- stepfun(rt[hit_idx], c(rev(cumsum(rev(1 / (G[hit_idx] ^ 2)))), 0), right = TRUE)
  return(list(K=K, Var=K.v))
}

estimateUCIPor <- function(rt, CR=NULL) {
  # Bolt: Use unlist(..., use.names=FALSE) instead of c(..., recursive=TRUE) for faster array flattening
  allRT <- sort(unique(unlist(rt, use.names=FALSE)))
  ncond <- length(rt)
  nt <- length(allRT)
  zero_fun <- function(t) rep(0, length(t))
  if (nt == 0) {
    return(list(H = zero_fun, Var = zero_fun))
  }
  
  if ( is.null(CR) || length(CR) != length(rt) ){
    CR <- vector("list", length(rt))
  }
  if (.ensure_hazard_rcpp()) {
    RT_list <- lapply(rt, as.numeric)
    CR_list <- vector("list", length(RT_list))
    for (i in seq_along(RT_list)) {
      if (i <= length(CR) && !is.null(CR[[i]])) {
        CR_list[[i]] <- as.logical(CR[[i]])
      } else {
        CR_list[[i]] <- rep(TRUE, length(RT_list[[i]]))
      }
    }
    H_fun <- function(t) ucip_eval_rcpp(RT_list, CR_list, as.numeric(t), FALSE)[, 1]
    V_fun <- function(t) ucip_eval_rcpp(RT_list, CR_list, as.numeric(t), FALSE)[, 2]
    return(list(H = H_fun, Var = V_fun))
  }
  
  # Bolt: Use lapply and rowSums to push vector accumulation into C code for better performance
  Hi_list <- lapply(seq_len(ncond), function(i) estimateNAH(rt[[i]], CR[[i]]))
  Hucip <- rowSums(matrix(vapply(Hi_list, function(h) h$H(allRT), numeric(nt)), nrow = nt))
  Hucip.v <- rowSums(matrix(vapply(Hi_list, function(h) h$Var(allRT), numeric(nt)), nrow = nt))
  
  Hucip <- stepfun(allRT, c(0, Hucip))
  Hucip.v <- stepfun(allRT, c(0, Hucip.v))
  
  return(list(H=Hucip, Var=Hucip.v))
}

estimateUCIPand <- function(rt, CR=NULL) {
  # Bolt: Use unlist(..., use.names=FALSE) instead of c(..., recursive=TRUE) for faster array flattening
  allRT <- sort(unique(unlist(rt, use.names=FALSE)))
  ncond <- length(rt)
  nt <- length(allRT)
  zero_fun <- function(t) rep(0, length(t))
  if (nt == 0) {
    return(list(K = zero_fun, Var = zero_fun))
  }
  
  if ( is.null(CR) || length(CR) != length(rt) ) {
    CR <- vector("list", length(rt))
  }
  if (.ensure_hazard_rcpp()) {
    RT_list <- lapply(rt, as.numeric)
    CR_list <- vector("list", length(RT_list))
    for (i in seq_along(RT_list)) {
      if (i <= length(CR) && !is.null(CR[[i]])) {
        CR_list[[i]] <- as.logical(CR[[i]])
      } else {
        CR_list[[i]] <- rep(TRUE, length(RT_list[[i]]))
      }
    }
    K_fun <- function(t) ucip_eval_rcpp(RT_list, CR_list, as.numeric(t), TRUE)[, 1]
    V_fun <- function(t) ucip_eval_rcpp(RT_list, CR_list, as.numeric(t), TRUE)[, 2]
    return(list(K = K_fun, Var = V_fun))
  }
  
  # Bolt: Use lapply and rowSums to push vector accumulation into C code for better performance
  Ki_list <- lapply(seq_len(ncond), function(i) estimateNAK(rt[[i]], CR[[i]]))
  Kucip <- rowSums(matrix(vapply(Ki_list, function(k) k$K(allRT), numeric(nt)), nrow = nt))
  Kucip.v <- rowSums(matrix(vapply(Ki_list, function(k) k$Var(allRT), numeric(nt)), nrow = nt))
  
  Kucip <- stepfun(allRT, c(Kucip,0),right=TRUE)
  Kucip.v <- stepfun(allRT, c(Kucip.v,0),right=TRUE)
  
  return(list(K=Kucip, Var=Kucip.v))
}

sic <- function(HH, HL, LH, LL, interpolate = FALSE) {
  
  all_rt <- unique(c(HH, HL, LH, LL))
  if (!is.finite(min(all_rt, na.rm = TRUE)) || !is.finite(max(all_rt, na.rm = TRUE))) {
    stop("sic(): rt inputs must contain finite values.")
  }
  
  RTall <- sort(all_rt)
  if (length(RTall) < 2) {
    stop("sic(): Need at least 2 time points to compute SIC.")
  }
  HH.ecdf <- ecdf(HH)
  HL.ecdf <- ecdf(HL)
  LH.ecdf <- ecdf(LH)
  LL.ecdf <- ecdf(LL)
  
  sicall <- LH.ecdf(RTall) + HL.ecdf(RTall) - HH.ecdf(RTall) - LL.ecdf(RTall)
  if (interpolate) {
    SIC <- approxfun(RTall, sicall, yleft = 0, yright = 0, rule = 2)
  } else {
    SIC <- stepfun(RTall, c(0, sicall))
  }
  MIC <- (mean(LL,na.rm=TRUE) - mean(LH,na.rm=TRUE))  - (mean(HL,na.rm=TRUE) - mean(HH,na.rm=TRUE))
  
  return(list(SIC=SIC, MIC=MIC, HH=HH.ecdf, HL=HL.ecdf, LH=LH.ecdf, LL=LL.ecdf))
}

resilience <- function(rt, CR=NULL, ratio=TRUE, rho=0) {
  # Bolt: Use unlist(..., use.names=FALSE) instead of c(..., recursive=TRUE) for faster array flattening
  times <- sort(unique(unlist(rt, use.names=FALSE)))
  
  Hab <- estimateNAH(rt[[1]], CR[[1]])
  Hay <- estimateNAH(rt[[2]], CR[[2]])
  Hxb <- estimateNAH(rt[[3]], CR[[3]])
  
  if(ratio) {
    R <- Hab$H(times)/ (Hay$H(times) + Hxb$H(times))
    R[is.nan(R)] <- NA
    R[is.infinite(R)] <- NA
    R <- approxfun(times, R)
    return(list(rt=R, Rtest=rtest))
  }
  else {
    R <- Hab$H(times) - (Hay$H(times) + Hxb$H(times))
    R[is.nan(R)] <- NA
    R[is.infinite(R)] <- NA
    R <- approxfun(times, R)
    Var.R <- Hab$Var(times) + Hay$Var(times) + Hxb$Var(times)
    Var.R <- approxfun(times, Var.R)
    return(list(rt=R, Var=Var.R))
  }
}

.safe_sic <- function(df_rule) {
  if (!is.data.frame(df_rule) || nrow(df_rule) == 0) return(list())
  hh <- df_rule$rt[df_rule$Channel1 == 2 & df_rule$Channel2 == 2 & df_rule$Correct == 1]
  hl <- df_rule$rt[df_rule$Channel1 == 2 & df_rule$Channel2 == 1 & df_rule$Correct == 1]
  lh <- df_rule$rt[df_rule$Channel1 == 1 & df_rule$Channel2 == 2 & df_rule$Correct == 1]
  ll <- df_rule$rt[df_rule$Channel1 == 1 & df_rule$Channel2 == 1 & df_rule$Correct == 1]
  if (length(hh) == 0 || length(hl) == 0 || length(lh) == 0 || length(ll) == 0) {
    return(list())
  }
  sic(hh, hl, lh, ll)
}

.compute_capacity <- function(data, logical_rules) {
  metrics <- list()
  if ("OR" %in% logical_rules) {
    OR <- data[data$Condition == "OR", , drop = FALSE]
    metrics$OR_Ct <- capacity.or(
      rt = list(OR$rt[OR$S == "AB"], OR$rt[OR$S == "AN"], OR$rt[OR$S == "NB"]),
      CR = list(OR$Correct[OR$S == "AB"], OR$Correct[OR$S == "AN"], OR$Correct[OR$S == "NB"])
    )
  }
  
  if ("AND" %in% logical_rules) {
    AND <- data[data$Condition == "AND", , drop = FALSE]
    metrics$AND_Ct <- capacity.and(
      rt = list(AND$rt[AND$S == "AB"], AND$rt[AND$S == "AN"], AND$rt[AND$S == "NB"]),
      CR = list(AND$Correct[AND$S == "AB"], AND$Correct[AND$S == "AN"], AND$Correct[AND$S == "NB"])
    )
    metrics$AND_Absence_Ct <- capacity.or(
      rt = list(AND$rt[AND$S == "NN"], AND$rt[AND$S == "AN"], AND$rt[AND$S == "NB"]),
      CR = list(AND$Correct[AND$S == "NN"], AND$Correct[AND$S == "AN"], AND$Correct[AND$S == "NB"])
    )
  }
  
  if ("ID" %in% logical_rules) {
    ID <- data[data$Condition == "ID", , drop = FALSE]
    metrics$ID_Ct <- capacity.id(
      dt.rt = ID$rt[ID$S == "AB"],
      nt.rt = ID$rt[ID$S == "NN"],
      st.rts = list(ID$rt[ID$S == "AN"], ID$rt[ID$S == "NB"]),
      dt.cr = ID$Correct[ID$S == "AB"],
      nt.cr = ID$Correct[ID$S == "NN"],
      st.crs = list(ID$Correct[ID$S == "AN"], ID$Correct[ID$S == "NB"])
    )
  }
  
  metrics
}

.compute_sics <- function(data, logical_rules) {
  out <- list()
  # Split the data once by Condition to avoid repeated O(N) subsetting inside the loop.
  # This improves complexity from O(N * R) to O(N + R) where N=nrow(data) and R=length(logical_rules).
  data_by_rule <- split(data, data$Condition)
  for (r in logical_rules) {
    df_rule <- data_by_rule[[r]]
    if (is.null(df_rule)) next
    sic_obj <- .safe_sic(df_rule)
    if (length(sic_obj) > 0) {
      out[[paste0(r, "_SIC")]] <- sic_obj
    }
  }
  out
}

.compute_cdfs <- function(data, logical_rules) {
  out <- list()
  cell_order <- c("AB", "AN", "NB", "NN")
  # Split the data once by Condition to avoid repeated O(N) subsetting inside the loop.
  data_by_rule <- split(data, data$Condition)
  for (r in logical_rules) {
    df_rule <- data_by_rule[[r]]
    if (is.null(df_rule)) next
    cdf_by_cell <- list()
    # Nested splitting by S further optimizes inner cell-level rt extraction.
    df_rule_by_s <- split(df_rule, df_rule$S)
    cells <- intersect(cell_order, names(df_rule_by_s))
    for (s in cells) {
      df_s <- df_rule_by_s[[s]]
      rt <- df_s$rt[df_s$Correct == 1]
      cdf_by_cell[[s]] <- .safe_ecdf(rt)
    }
    out[[r]] <- cdf_by_cell
  }
  out
}
