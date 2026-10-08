# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8538151 -0.8634263
      m_lMd -0.3858405 -0.5331963
      s     -0.9058225 -0.6449558
      t0    -1.9716145 -1.7915876

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.9337886 -0.4474248 -0.8324346 -1.7778773 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                      m      m_lMd     s         t0
      m     0.002543296 0.00000000 0.000 0.00000000
      m_lMd 0.000000000 0.04675125 0.000 0.00000000
      s     0.000000000 0.00000000 0.111 0.00000000
      t0    0.000000000 0.00000000 0.000 0.04459593

---

    Code
      compare(list(diag = LNR_diag), stage = "preburn", cores_for_props = 1)
    Output
             MD wMD  WAIC wWAIC  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      diag -377   1 35500     1 2980    1 4802     1       1822  1158  -528 -664

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -1.0967399 -1.0655185
      m_lMd -0.4086515 -0.6763077
      s     -0.7068633 -0.3945553
      t0    -1.6100372 -1.4279430

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.7016637 -0.7645540  0.2702122 -1.2878767 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                    m     m_lMd         s         t0
      m     0.8525755 0.0000000 0.0000000 0.00000000
      m_lMd 0.0000000 0.0451548 0.0000000 0.00000000
      s     0.0000000 0.0000000 1.4599027 0.33595238
      t0    0.0000000 0.0000000 0.3359524 0.08245812

---

    Code
      compare(list(blocked = LNR_blocked), BayesFactor = FALSE, stage = "preburn",
      cores_for_props = 1)
    Output
               WAIC wWAIC  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      blocked 15219     1 2095    1 3483     1       1388   707  -355 -681

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
      single -382   1 2493     1 487    1  988     1        501   -15  -378 -516

