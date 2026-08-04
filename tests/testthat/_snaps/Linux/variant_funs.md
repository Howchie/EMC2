# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9667867 -0.8904360
      m_lMd -0.3652931 -0.5496101
      s     -0.7973441 -0.6263701
      t0    -1.8995477 -1.8653071

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.8480291 -0.5006714 -0.6496058 -1.8798856 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                     m      m_lMd          s           t0
      m     0.03864755 0.00000000 0.00000000 0.0000000000
      m_lMd 0.00000000 0.01579569 0.00000000 0.0000000000
      s     0.00000000 0.00000000 0.05342463 0.0000000000
      t0    0.00000000 0.00000000 0.00000000 0.0005324957

---

    Code
      compare(list(diag = LNR_diag), stage = "preburn", cores_for_props = 1)
    Output
             MD wMD  WAIC wWAIC  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      diag -477   1 12178     1 1340    1 2349     1       1008   332  -568 -677

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8411010 -0.9929560
      m_lMd -0.3674576 -0.5308419
      s     -0.7712796 -0.5672621
      t0    -2.1668809 -1.7320239

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.8828587 -0.6090358 -0.5186590 -1.2020880 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                     m      m_lMd           s          t0
      m     0.09996624 0.00000000  0.00000000  0.00000000
      m_lMd 0.00000000 0.01558194  0.00000000  0.00000000
      s     0.00000000 0.00000000  0.03091205 -0.01111874
      t0    0.00000000 0.00000000 -0.01111874  0.28305044

---

    Code
      compare(list(blocked = LNR_blocked), stage = "preburn", cores_for_props = 1)
    Output
                MD wMD WAIC wWAIC DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      blocked -460   1 2799     1 826    1 1581     1        754    72  -499 -682

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8382160 -0.7452781
      m_lMd -0.3282055 -0.7591053
      s     -1.0620691 -0.4747831
      t0    -2.3725880 -1.9081209

---

    Code
      compare(list(single = LNR_single), stage = "preburn", cores_for_props = 1)
    Output
               MD wMD WAIC wWAIC DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single -436   1  645     1 230    1  665     1        436  -206  -582 -641

