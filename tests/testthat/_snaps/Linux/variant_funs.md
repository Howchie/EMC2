# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9839639 -0.9891628
      m_lMd -0.3737356 -0.6324346
      s     -0.7318134 -0.5006487
      t0    -1.8540617 -1.6935394

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -1.0409894 -0.4781946 -0.6267372 -1.7319364 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                      m      m_lMd         s         t0
      m     0.001527041 0.00000000 0.0000000 0.00000000
      m_lMd 0.000000000 0.03132476 0.0000000 0.00000000
      s     0.000000000 0.00000000 0.0391503 0.00000000
      t0    0.000000000 0.00000000 0.0000000 0.00761244

---

    Code
      compare(list(diag = LNR_diag), stage = "preburn", cores_for_props = 1)
    Output
             MD wMD  WAIC wWAIC  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      diag -534   1 35518     1 2964    1 4783     1       1819  1146  -544 -673

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -1.0122111 -1.0269133
      m_lMd -0.3439454 -0.6090662
      s     -0.7502023 -0.4850585
      t0    -1.8568173 -1.6852100

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.6964647 -0.6973460  0.2749363 -1.4936724 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                    m      m_lMd         s        t0
      m     0.6694671 0.00000000 0.0000000 0.0000000
      m_lMd 0.0000000 0.04453015 0.0000000 0.0000000
      s     0.0000000 0.00000000 1.6444167 0.4118455
      t0    0.0000000 0.00000000 0.4118455 0.1091788

---

    Code
      compare(list(blocked = LNR_blocked), stage = "preburn", cores_for_props = 1)
    Output
               MD wMD  WAIC wWAIC  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      blocked 213   1 17257     1 2542    1 4150     1       1608   933  -193 -675

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8913198 -0.6035807
      m_lMd -0.2686419 -0.6859228
      s     -0.9715802 -0.8603207
      t0    -2.6292223 -2.6116865

---

    Code
      compare(list(single = LNR_single), stage = "preburn", cores_for_props = 1)
    Output
               MD wMD WAIC wWAIC DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single -122   1  555     1 543    1 1066     1        523    19  -306 -504

