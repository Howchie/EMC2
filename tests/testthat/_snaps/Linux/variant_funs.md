# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -1.0615914 -1.0443479
      m_lMd -0.5174121 -0.5738478
      s     -0.6466361 -0.5030446
      t0    -1.6892925 -1.6701134

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -1.1021507 -0.5083796 -0.5498154 -1.6548610 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                      m      m_lMd          s          t0
      m     0.001196309 0.00000000 0.00000000 0.000000000
      m_lMd 0.000000000 0.01020018 0.00000000 0.000000000
      s     0.000000000 0.00000000 0.02559477 0.000000000
      t0    0.000000000 0.00000000 0.00000000 0.002010657

---

    Code
      compare(list(diag = LNR_diag), stage = "preburn", cores_for_props = 1)
    Output
             MD wMD  WAIC wWAIC  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      diag -252   1 52732     1 4780    1 7513     1       2733  2048  -477 -685

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -1.0259319 -1.0785466
      m_lMd -0.4415315 -0.7958626
      s     -0.6836014 -0.3666746
      t0    -1.7751772 -1.5652330

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.6645021 -0.8620663  0.3361565 -1.4065058 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                    m      m_lMd         s        t0
      m     0.8841082 0.00000000 0.0000000 0.0000000
      m_lMd 0.0000000 0.05999124 0.0000000 0.0000000
      s     0.0000000 0.00000000 1.5803127 0.3962279
      t0    0.0000000 0.00000000 0.3962279 0.1056520

---

    Code
      compare(list(blocked = LNR_blocked), stage = "preburn", cores_for_props = 1)
    Output
               MD wMD  WAIC wWAIC  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      blocked -29   1 19593     1 2582    1 4209     1       1627   955   -76 -673

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9864667 -0.9569580
      m_lMd -0.5497707 -0.7095726
      s     -0.6295003 -0.2639983
      t0    -1.8805018 -1.6472347

---

    Code
      compare(list(single = LNR_single), stage = "preburn", cores_for_props = 1)
    Output
               MD wMD WAIC wWAIC DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single -380   1  298     1 118    1  450     1        333  -215  -544 -548

