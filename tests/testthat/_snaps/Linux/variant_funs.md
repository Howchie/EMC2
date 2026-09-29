# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                 as1t       bd6t
      m     -1.067415 -0.9219158
      m_lMd -0.454594 -0.6025142
      s     -0.653663 -0.5549081
      t0    -1.711956 -1.7550733

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -1.0421700 -0.5143877 -0.6175581 -1.6922719 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                      m      m_lMd          s         t0
      m     0.002189247 0.00000000 0.00000000 0.00000000
      m_lMd 0.000000000 0.01107427 0.00000000 0.00000000
      s     0.000000000 0.00000000 0.01312722 0.00000000
      t0    0.000000000 0.00000000 0.00000000 0.00352035

---

    Code
      compare(list(diag = LNR_diag), stage = "preburn", cores_for_props = 1)
    Output
             MD wMD  WAIC wWAIC  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      diag -325   1 35762     1 3109    1 4990     1       1881  1227  -517 -654

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -1.0997939 -1.0151308
      m_lMd -0.3956896 -0.6092596
      s     -0.6551528 -0.5078630
      t0    -1.6895472 -1.6479004

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
                  m         m_lMd             s            t0 
      -0.6667633370 -0.7062561056  0.0007264556 -1.5971139997 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                   m      m_lMd          s          t0
      m     0.896541 0.00000000 0.00000000 0.000000000
      m_lMd 0.000000 0.03590755 0.00000000 0.000000000
      s     0.000000 0.00000000 0.80516389 0.080072263
      t0    0.000000 0.00000000 0.08007226 0.008721289

---

    Code
      compare(list(blocked = LNR_blocked), BayesFactor = FALSE, stage = "preburn",
      cores_for_props = 1)
    Output
               WAIC wWAIC  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      blocked 16653     1 2298    1 3787     1       1489   809  -415 -680

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8757191 -0.6840740
      m_lMd -0.8752876 -0.5154932
      s     -0.3120741 -0.8428329
      t0    -1.8689038 -2.3061279

---

    Code
      compare(list(single = LNR_single), stage = "preburn", cores_for_props = 1)
    Output
               MD wMD WAIC wWAIC DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single -148   1 1021     1 566    1 1119     1        553    13  -312 -540

