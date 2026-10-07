# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8365969 -0.8843513
      m_lMd -0.3481828 -0.5614779
      s     -0.9692018 -0.7138364
      t0    -2.0846798 -1.8364538

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.9070235 -0.4278811 -0.8494327 -1.9082165 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                      m      m_lMd          s         t0
      m     0.001196942 0.00000000 0.00000000 0.00000000
      m_lMd 0.000000000 0.02880734 0.00000000 0.00000000
      s     0.000000000 0.00000000 0.04242218 0.00000000
      t0    0.000000000 0.00000000 0.00000000 0.01292785

---

    Code
      compare(list(diag = LNR_diag), stage = "preburn", cores_for_props = 1)
    Output
             MD wMD  WAIC wWAIC  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      diag -344   1 35515     1 2990    1 4808     1       1818  1172  -532 -646

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -1.0044262 -0.9758783
      m_lMd -0.4541240 -0.6225409
      s     -0.6811102 -0.5186594
      t0    -1.6398392 -1.5839028

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
                m       m_lMd           s          t0 
      -0.61141142 -0.72782154  0.01278949 -1.52050320 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                    m      m_lMd         s         t0
      m     0.8570679 0.00000000 0.0000000 0.00000000
      m_lMd 0.0000000 0.02948191 0.0000000 0.00000000
      s     0.0000000 0.00000000 0.8770423 0.10574899
      t0    0.0000000 0.00000000 0.1057490 0.01390065

---

    Code
      compare(list(blocked = LNR_blocked), BayesFactor = FALSE, stage = "preburn",
      cores_for_props = 1)
    Output
               WAIC wWAIC  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      blocked 15354     1 2219    1 3660     1       1441   779  -309 -662

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -1.3012780 -1.1078428
      m_lMd -0.4937164 -0.7099383
      s     -0.6118380 -0.3283319
      t0    -1.4657967 -1.4620786

---

    Code
      compare(list(single = LNR_single), stage = "preburn", cores_for_props = 1)
    Output
               MD wMD WAIC wWAIC DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single -250   1 3249     1 833    1 1498     1        665   168  -294 -497

