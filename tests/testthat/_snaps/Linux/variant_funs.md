# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.7325784 -0.7327157
      m_lMd -0.2677144 -0.5187610
      s     -1.2416886 -0.7869572
      t0    -6.5464200 -2.6152976

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
                m       m_lMd           s          t0 
      -0.76414951 -0.40270515 -0.95096128  0.02110497 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                      m      m_lMd          s       t0
      m     0.001543597 0.00000000 0.00000000  0.00000
      m_lMd 0.000000000 0.02049045 0.00000000  0.00000
      s     0.000000000 0.00000000 0.06423483  0.00000
      t0    0.000000000 0.00000000 0.00000000 47.36962

---

    Code
      compare(list(diag = LNR_diag), stage = "preburn", cores_for_props = 1)
    Output
             MD wMD  WAIC wWAIC  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      diag -247   1 35537     1 3012    1 4839     1       1827  1185   568 -642

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                 as1t       bd6t
      m     -1.064849 -1.1168726
      m_lMd -0.383237 -0.7583110
      s     -0.678852 -0.3574443
      t0    -1.887388 -1.6166929

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.7029649 -0.8303607  0.4317363 -1.3768481 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                    m      m_lMd         s        t0
      m     0.8815767 0.00000000 0.0000000 0.0000000
      m_lMd 0.0000000 0.06774806 0.0000000 0.0000000
      s     0.0000000 0.00000000 1.8318692 0.5891878
      t0    0.0000000 0.00000000 0.5891878 0.1996984

---

    Code
      compare(list(blocked = LNR_blocked), BayesFactor = FALSE, stage = "preburn",
      cores_for_props = 1)
    Output
               WAIC wWAIC  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      blocked 17054     1 2709    1 4402     1       1693  1017  -220 -676

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -1.2767377 -1.1170794
      m_lMd -0.4025261 -0.6014827
      s     -0.6403422 -0.3110400
      t0    -1.7936388 -1.7910356

---

    Code
      compare(list(single = LNR_single), stage = "preburn", cores_for_props = 1)
    Output
               MD wMD WAIC wWAIC DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single -181   1 1219     1 423    1  891     1        468   -45  -446 -514

