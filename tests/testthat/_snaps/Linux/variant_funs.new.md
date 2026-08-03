# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9772701 -1.0213333
      m_lMd -0.4197685 -0.5570124
      s     -0.6955529 -0.4367439
      t0    -1.7993578 -1.6277934

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -1.0909752 -0.6990911 -0.5416060 -1.8165074 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                     m     m_lMd          s         t0
      m     0.02060934 0.0000000 0.00000000 0.00000000
      m_lMd 0.00000000 0.0833289 0.00000000 0.00000000
      s     0.00000000 0.0000000 0.02655852 0.00000000
      t0    0.00000000 0.0000000 0.00000000 0.02349171

---

    Code
      compare(list(diag = LNR_diag), stage = "preburn", cores_for_props = 1)
    Output
             MD wMD WAIC wWAIC DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      diag -521   1 8493     1 770    1 1494     1        724    46  -545 -678

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9446515 -0.9762585
      m_lMd -0.3853960 -0.6073520
      s     -0.8171087 -0.4904382
      t0    -1.9353535 -1.6775655

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.9919638 -1.5683960 -0.7658148 -1.7901868 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                      m     m_lMd            s           t0
      m     0.004998066 0.0000000  0.000000000  0.000000000
      m_lMd 0.000000000 0.8127034  0.000000000  0.000000000
      s     0.000000000 0.0000000  0.009510822 -0.002298335
      t0    0.000000000 0.0000000 -0.002298335  0.070525428

---

    Code
      compare(list(blocked = LNR_blocked), stage = "preburn", cores_for_props = 1)
    Output
                MD wMD  WAIC wWAIC  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      blocked -495   1 64562     1 5841    1 9099     1       3258  2583  -614 -675

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                   as1t       bd6t
      m     -0.86010555 -0.9356648
      m_lMd -0.03729243 -0.8312491
      s     -0.74467589 -0.5590139
      t0    -2.21601123 -1.6965391

---

    Code
      compare(list(single = LNR_single), stage = "preburn", cores_for_props = 1)
    Output
               MD wMD WAIC wWAIC DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single -314   1 1697     1 471    1  987     1        516   -45  -434 -561

