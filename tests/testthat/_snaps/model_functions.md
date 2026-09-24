# LNR

    Code
      init_chains(LNR_s, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $theta_mu
                         [,1]
      m            0.10295498
      m_lMd       -0.17485347
      m_Eneutral   0.15362521
      m_Eaccuracy  0.04786557
      s           -0.11250747
      t0           0.47721481
      
      $theta_var
      , , 1
      
                              m        m_lMd  m_Eneutral   m_Eaccuracy           s
      m            6.587165e-02  0.007385361 -0.01742800  7.080807e-05 -0.02937366
      m_lMd        7.385361e-03  0.070488318  0.00962596  7.109282e-03 -0.01161209
      m_Eneutral  -1.742800e-02  0.009625960  0.10169241 -3.978216e-02  0.01122457
      m_Eaccuracy  7.080807e-05  0.007109282 -0.03978216  1.042876e-01  0.02904821
      s           -2.937366e-02 -0.011612087  0.01122457  2.904821e-02  0.07505536
      t0           1.996821e-02  0.007166038 -0.03560956  2.469992e-02 -0.01637705
                            t0
      m            0.019968211
      m_lMd        0.007166038
      m_Eneutral  -0.035609562
      m_Eaccuracy  0.024699915
      s           -0.016377047
      t0           0.088132661
      
      
      $a_half
                       [,1]
      m           1.1049396
      m_lMd       0.8879827
      m_Eneutral  0.3506475
      m_Eaccuracy 1.1383694
      s           2.6196376
      t0          0.4212093
      
      $alpha
      , , 1
      
                         as1t        bd6t
      m           -0.03213819 -0.13610410
      m_lMd       -0.04820415 -0.34999024
      m_Eneutral  -0.05186968  0.01024248
      m_Eaccuracy  0.37267253  0.27619277
      s            0.05738331  0.41990689
      t0           0.39541015  0.30842879
      
      
      $stage
      [1] "init"
      
      $subj_ll
                [,1]
      as1t -18650.94
      bd6t -19526.53
      
      $last_theta_var_inv
                 [,1]       [,2]      [,3]       [,4]      [,5]       [,6]
      [1,] 19.5197554 -0.7711759  1.536017 -0.8501984  7.145064 -2.1732755
      [2,] -0.7711759 15.8200358 -3.675052 -3.2959066  3.758356 -0.9743898
      [3,]  1.5360172 -3.6750520 14.370256  5.9479347 -3.675075  3.4071615
      [4,] -0.8501984 -3.2959066  5.947935 15.0417806 -8.180328 -2.8718235
      [5,]  7.1450642  3.7583563 -3.675075 -8.1803284 21.025694  2.7903130
      [6,] -2.1732755 -0.9743898  3.407161 -2.8718235  2.790313 14.6181601
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_LNR, design_LNR, n_trials = 10)
    Output
          trials subjects        E     S     R        rt LT  UT LC  UC
      1        1     as1t    speed  left  left 0.3224880  0 Inf  0 Inf
      3        2     as1t    speed  left right 0.2272760  0 Inf  0 Inf
      5        3     as1t    speed  left right 0.2354363  0 Inf  0 Inf
      7        4     as1t    speed  left right 0.2201006  0 Inf  0 Inf
      9        5     as1t    speed  left  left 0.3530253  0 Inf  0 Inf
      11       6     as1t    speed  left right 0.6494876  0 Inf  0 Inf
      13       7     as1t    speed  left right 0.2152369  0 Inf  0 Inf
      15       8     as1t    speed  left right 0.2585436  0 Inf  0 Inf
      17       9     as1t    speed  left right 0.2026706  0 Inf  0 Inf
      19      10     as1t    speed  left right 0.4304292  0 Inf  0 Inf
      21      11     as1t  neutral  left right 0.4319677  0 Inf  0 Inf
      23      12     as1t  neutral  left right 0.9394224  0 Inf  0 Inf
      25      13     as1t  neutral  left right 0.2243063  0 Inf  0 Inf
      27      14     as1t  neutral  left  left 0.2467995  0 Inf  0 Inf
      29      15     as1t  neutral  left right 0.3424219  0 Inf  0 Inf
      31      16     as1t  neutral  left right 0.2583192  0 Inf  0 Inf
      33      17     as1t  neutral  left right 0.6128815  0 Inf  0 Inf
      35      18     as1t  neutral  left right 0.2976982  0 Inf  0 Inf
      37      19     as1t  neutral  left right 0.2067883  0 Inf  0 Inf
      39      20     as1t  neutral  left right 0.3028746  0 Inf  0 Inf
      41      21     as1t accuracy  left right 0.4863136  0 Inf  0 Inf
      43      22     as1t accuracy  left right 1.0110174  0 Inf  0 Inf
      45      23     as1t accuracy  left right 0.2288424  0 Inf  0 Inf
      47      24     as1t accuracy  left right 0.2120292  0 Inf  0 Inf
      49      25     as1t accuracy  left right 0.3197284  0 Inf  0 Inf
      51      26     as1t accuracy  left right 0.2548833  0 Inf  0 Inf
      53      27     as1t accuracy  left right 0.2974396  0 Inf  0 Inf
      55      28     as1t accuracy  left right 0.2021762  0 Inf  0 Inf
      57      29     as1t accuracy  left  left 0.2312319  0 Inf  0 Inf
      59      30     as1t accuracy  left right 0.2014580  0 Inf  0 Inf
      61      31     as1t    speed right  left 0.2269831  0 Inf  0 Inf
      63      32     as1t    speed right  left 0.4616264  0 Inf  0 Inf
      65      33     as1t    speed right  left 0.5132418  0 Inf  0 Inf
      67      34     as1t    speed right right 0.4121504  0 Inf  0 Inf
      69      35     as1t    speed right  left 0.2313245  0 Inf  0 Inf
      71      36     as1t    speed right right 0.2719558  0 Inf  0 Inf
      73      37     as1t    speed right  left 0.2434625  0 Inf  0 Inf
      75      38     as1t    speed right  left 0.5606815  0 Inf  0 Inf
      77      39     as1t    speed right  left 0.2716943  0 Inf  0 Inf
      79      40     as1t    speed right  left 0.2462127  0 Inf  0 Inf
      81      41     as1t  neutral right  left 0.3365330  0 Inf  0 Inf
      83      42     as1t  neutral right  left 0.4470369  0 Inf  0 Inf
      85      43     as1t  neutral right  left 1.1962713  0 Inf  0 Inf
      87      44     as1t  neutral right  left 0.3131830  0 Inf  0 Inf
      89      45     as1t  neutral right  left 0.3994020  0 Inf  0 Inf
      91      46     as1t  neutral right  left 1.1008538  0 Inf  0 Inf
      93      47     as1t  neutral right  left 0.2045263  0 Inf  0 Inf
      95      48     as1t  neutral right right 0.4579290  0 Inf  0 Inf
      97      49     as1t  neutral right  left 0.4289439  0 Inf  0 Inf
      99      50     as1t  neutral right  left 0.3078732  0 Inf  0 Inf
      101     51     as1t accuracy right  left 0.3484332  0 Inf  0 Inf
      103     52     as1t accuracy right  left 0.4319588  0 Inf  0 Inf
      105     53     as1t accuracy right right 0.3255497  0 Inf  0 Inf
      107     54     as1t accuracy right right 0.3777698  0 Inf  0 Inf
      109     55     as1t accuracy right  left 0.5117870  0 Inf  0 Inf
      111     56     as1t accuracy right  left 0.2840452  0 Inf  0 Inf
      113     57     as1t accuracy right  left 0.2077781  0 Inf  0 Inf
      115     58     as1t accuracy right  left 0.4775184  0 Inf  0 Inf
      117     59     as1t accuracy right  left 0.3134105  0 Inf  0 Inf
      119     60     as1t accuracy right  left 0.2164368  0 Inf  0 Inf
      121      1     bd6t    speed  left right 0.2500167  0 Inf  0 Inf
      123      2     bd6t    speed  left  left 0.5153744  0 Inf  0 Inf
      125      3     bd6t    speed  left right 0.2136713  0 Inf  0 Inf
      127      4     bd6t    speed  left right 0.2067060  0 Inf  0 Inf
      129      5     bd6t    speed  left right 0.4384221  0 Inf  0 Inf
      131      6     bd6t    speed  left  left 0.4027332  0 Inf  0 Inf
      133      7     bd6t    speed  left right 0.2998384  0 Inf  0 Inf
      135      8     bd6t    speed  left right 0.2280791  0 Inf  0 Inf
      137      9     bd6t    speed  left right 0.3031266  0 Inf  0 Inf
      139     10     bd6t    speed  left right 0.4160953  0 Inf  0 Inf
      141     11     bd6t  neutral  left  left 0.2150449  0 Inf  0 Inf
      143     12     bd6t  neutral  left right 0.3923806  0 Inf  0 Inf
      145     13     bd6t  neutral  left right 0.2521765  0 Inf  0 Inf
      147     14     bd6t  neutral  left right 0.8420605  0 Inf  0 Inf
      149     15     bd6t  neutral  left right 0.5083832  0 Inf  0 Inf
      151     16     bd6t  neutral  left right 0.2130774  0 Inf  0 Inf
      153     17     bd6t  neutral  left right 0.2278625  0 Inf  0 Inf
      155     18     bd6t  neutral  left right 0.2352806  0 Inf  0 Inf
      157     19     bd6t  neutral  left right 0.2522588  0 Inf  0 Inf
      159     20     bd6t  neutral  left right 0.2330901  0 Inf  0 Inf
      161     21     bd6t accuracy  left right 0.2146149  0 Inf  0 Inf
      163     22     bd6t accuracy  left right 0.3949964  0 Inf  0 Inf
      165     23     bd6t accuracy  left right 0.2933916  0 Inf  0 Inf
      167     24     bd6t accuracy  left right 0.6367135  0 Inf  0 Inf
      169     25     bd6t accuracy  left right 0.2121097  0 Inf  0 Inf
      171     26     bd6t accuracy  left right 0.2100779  0 Inf  0 Inf
      173     27     bd6t accuracy  left right 0.2590189  0 Inf  0 Inf
      175     28     bd6t accuracy  left right 0.2225345  0 Inf  0 Inf
      177     29     bd6t accuracy  left right 0.4171558  0 Inf  0 Inf
      179     30     bd6t accuracy  left  left 0.3514937  0 Inf  0 Inf
      181     31     bd6t    speed right  left 4.6347640  0 Inf  0 Inf
      183     32     bd6t    speed right  left 0.2294992  0 Inf  0 Inf
      185     33     bd6t    speed right  left 0.2856515  0 Inf  0 Inf
      187     34     bd6t    speed right  left 0.2041104  0 Inf  0 Inf
      189     35     bd6t    speed right  left 0.2435389  0 Inf  0 Inf
      191     36     bd6t    speed right  left 0.3757083  0 Inf  0 Inf
      193     37     bd6t    speed right  left 0.2938369  0 Inf  0 Inf
      195     38     bd6t    speed right right 0.2505960  0 Inf  0 Inf
      197     39     bd6t    speed right  left 0.7438186  0 Inf  0 Inf
      199     40     bd6t    speed right right 0.2862129  0 Inf  0 Inf
      201     41     bd6t  neutral right  left 0.2076989  0 Inf  0 Inf
      203     42     bd6t  neutral right  left 0.2546633  0 Inf  0 Inf
      205     43     bd6t  neutral right  left 0.3357014  0 Inf  0 Inf
      207     44     bd6t  neutral right  left 0.2485620  0 Inf  0 Inf
      209     45     bd6t  neutral right  left 0.2284170  0 Inf  0 Inf
      211     46     bd6t  neutral right  left 0.7839503  0 Inf  0 Inf
      213     47     bd6t  neutral right  left 0.2234772  0 Inf  0 Inf
      215     48     bd6t  neutral right  left 0.7135456  0 Inf  0 Inf
      217     49     bd6t  neutral right  left 0.3465151  0 Inf  0 Inf
      219     50     bd6t  neutral right  left 0.3611702  0 Inf  0 Inf
      221     51     bd6t accuracy right  left 0.2770322  0 Inf  0 Inf
      223     52     bd6t accuracy right  left 0.4660723  0 Inf  0 Inf
      225     53     bd6t accuracy right right 0.5665229  0 Inf  0 Inf
      227     54     bd6t accuracy right right 0.2634043  0 Inf  0 Inf
      229     55     bd6t accuracy right  left 0.2907393  0 Inf  0 Inf
      231     56     bd6t accuracy right  left 0.5179346  0 Inf  0 Inf
      233     57     bd6t accuracy right  left 0.2532600  0 Inf  0 Inf
      235     58     bd6t accuracy right  left 0.2306803  0 Inf  0 Inf
      237     59     bd6t accuracy right  left 0.2222017  0 Inf  0 Inf
      239     60     bd6t accuracy right  left 0.9323968  0 Inf  0 Inf

# LBA

    Code
      init_chains(LBA_s, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $theta_mu
                         [,1]
      v            0.10295498
      v_lMd       -0.17485347
      sv_lMTRUE    0.15362521
      B            0.04786557
      B_Eneutral  -0.11250747
      B_Eaccuracy  0.47721481
      B_lRright    0.41720870
      A           -0.03732249
      t0          -0.77359023
      
      $theta_var
      , , 1
      
                              v        v_lMd     sv_lMTRUE            B   B_Eneutral
      v            0.0259490148  0.005010764  0.0008039212  0.001039186  0.004671743
      v_lMd        0.0050107643  0.054198668 -0.0063395659  0.016853096  0.011241011
      sv_lMTRUE    0.0008039212 -0.006339566  0.0660491934 -0.021845123  0.024543796
      B            0.0010391859  0.016853096 -0.0218451231  0.058738576  0.005181317
      B_Eneutral   0.0046717429  0.011241011  0.0245437963  0.005181317  0.057368805
      B_Eaccuracy  0.0062510813 -0.002114676 -0.0040508422  0.004327501 -0.002720295
      B_lRright    0.0031684940  0.003367429  0.0059128710  0.014675459  0.019979670
      A           -0.0131009964 -0.029069087  0.0167338768 -0.022679508 -0.005251711
      t0           0.0040141996  0.016844431 -0.0079536906  0.011913763  0.008338580
                   B_Eaccuracy    B_lRright            A           t0
      v            0.006251081  0.003168494 -0.013100996  0.004014200
      v_lMd       -0.002114676  0.003367429 -0.029069087  0.016844431
      sv_lMTRUE   -0.004050842  0.005912871  0.016733877 -0.007953691
      B            0.004327501  0.014675459 -0.022679508  0.011913763
      B_Eneutral  -0.002720295  0.019979670 -0.005251711  0.008338580
      B_Eaccuracy  0.029195495  0.005568241 -0.009362256 -0.002326190
      B_lRright    0.005568241  0.037848420 -0.001496246  0.005725596
      A           -0.009362256 -0.001496246  0.078727638 -0.016972357
      t0          -0.002326190  0.005725596 -0.016972357  0.044669825
      
      
      $a_half
                       [,1]
      v           1.3638374
      v_lMd       0.3851234
      sv_lMTRUE   0.4683564
      B           0.8473333
      B_Eneutral  0.1829583
      B_Eaccuracy 0.4256923
      B_lRright   0.2080372
      A           0.8331091
      t0          1.0314865
      
      $alpha
      , , 1
      
                         as1t       bd6t
      v           -0.18365864  0.2351246
      v_lMd        0.02070972 -0.3027253
      sv_lMTRUE    0.15264692  0.2118511
      B           -0.08907052 -0.1551697
      B_Eneutral  -0.36914892 -0.2107319
      B_Eaccuracy  0.36955936  0.3663802
      B_lRright    0.38985219  0.4408931
      A           -0.16545291  0.4656651
      t0          -0.98002284 -0.9267676
      
      
      $stage
      [1] "init"
      
      $subj_ll
                [,1]
      as1t -15911.36
      bd6t -16519.59
      
      $last_theta_var_inv
                  [,1]       [,2]        [,3]        [,4]       [,5]       [,6]
       [1,] 44.9040894 -0.8014289  -0.9812410   3.3620715  -2.191829 -8.1337372
       [2,] -0.8014289 26.2380707   0.2153557  -3.8519366  -3.769284  4.0043411
       [3,] -0.9812410  0.2153557  23.0950016   8.2391363 -10.580642  1.0256725
       [4,]  3.3620715 -3.8519366   8.2391363  25.6470100  -1.349849 -0.6209953
       [5,] -2.1918285 -3.7692843 -10.5806416  -1.3498493  27.998252  4.1803956
       [6,] -8.1337372  4.0043411   1.0256725  -0.6209953   4.180396 40.5675700
       [7,] -1.9141836  1.6771437  -1.8938780 -10.0302584 -12.222883 -8.1641299
       [8,]  6.7725554  7.5357230  -2.5434395   4.1758907   1.649528  5.5506770
       [9,] -2.0002672 -5.1959610   3.2265984  -1.1309205  -2.720966  4.0568270
                  [,7]      [,8]      [,9]
       [1,]  -1.914184  6.772555 -2.000267
       [2,]   1.677144  7.535723 -5.195961
       [3,]  -1.893878 -2.543439  3.226598
       [4,] -10.030258  4.175891 -1.130920
       [5,] -12.222883  1.649528 -2.720966
       [6,]  -8.164130  5.550677  4.056827
       [7,]  38.521104 -3.812009 -2.651828
       [8,]  -3.812009 19.684888  2.932193
       [9,]  -2.651828  2.932193 27.574874
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_LBA, design_LBA, n_trials = 10)
    Output
          trials subjects        E     S     R        rt LT  UT LC  UC
      1        1     as1t    speed  left  left 1.2464691  0 Inf  0 Inf
      3        2     as1t    speed  left  left 1.3033527  0 Inf  0 Inf
      5        3     as1t    speed  left  left 1.1804751  0 Inf  0 Inf
      7        4     as1t    speed  left  left 1.3417653  0 Inf  0 Inf
      9        5     as1t    speed  left  left 0.9860655  0 Inf  0 Inf
      11       6     as1t    speed  left  left 1.4131841  0 Inf  0 Inf
      13       7     as1t    speed  left  left 1.4020763  0 Inf  0 Inf
      15       8     as1t    speed  left  left 1.1330933  0 Inf  0 Inf
      17       9     as1t    speed  left  left 1.2275937  0 Inf  0 Inf
      19      10     as1t    speed  left right 1.3591356  0 Inf  0 Inf
      21      11     as1t  neutral  left  left 1.3229579  0 Inf  0 Inf
      23      12     as1t  neutral  left  left 2.5118615  0 Inf  0 Inf
      25      13     as1t  neutral  left  left 1.7405644  0 Inf  0 Inf
      27      14     as1t  neutral  left right 1.5416069  0 Inf  0 Inf
      29      15     as1t  neutral  left  left 1.7208276  0 Inf  0 Inf
      31      16     as1t  neutral  left  left 2.4579638  0 Inf  0 Inf
      33      17     as1t  neutral  left  left 1.6632346  0 Inf  0 Inf
      35      18     as1t  neutral  left  left 2.2951800  0 Inf  0 Inf
      37      19     as1t  neutral  left  left 1.7296549  0 Inf  0 Inf
      39      20     as1t  neutral  left  left 2.2140710  0 Inf  0 Inf
      41      21     as1t accuracy  left  left 1.9294449  0 Inf  0 Inf
      43      22     as1t accuracy  left  left 1.7144434  0 Inf  0 Inf
      45      23     as1t accuracy  left right 1.1862126  0 Inf  0 Inf
      47      24     as1t accuracy  left  left 1.3224174  0 Inf  0 Inf
      49      25     as1t accuracy  left right 1.4690119  0 Inf  0 Inf
      51      26     as1t accuracy  left  left 1.5071558  0 Inf  0 Inf
      53      27     as1t accuracy  left  left 1.3793004  0 Inf  0 Inf
      55      28     as1t accuracy  left right 1.9042128  0 Inf  0 Inf
      57      29     as1t accuracy  left  left 1.2244016  0 Inf  0 Inf
      59      30     as1t accuracy  left  left 1.6414049  0 Inf  0 Inf
      61      31     as1t    speed right right 1.6145551  0 Inf  0 Inf
      63      32     as1t    speed right right 0.9085772  0 Inf  0 Inf
      65      33     as1t    speed right right 1.2890678  0 Inf  0 Inf
      67      34     as1t    speed right  left 1.2889023  0 Inf  0 Inf
      69      35     as1t    speed right  left 1.1192877  0 Inf  0 Inf
      71      36     as1t    speed right right 1.0806063  0 Inf  0 Inf
      73      37     as1t    speed right right 1.8217263  0 Inf  0 Inf
      75      38     as1t    speed right right 0.9863785  0 Inf  0 Inf
      77      39     as1t    speed right right 1.3769000  0 Inf  0 Inf
      79      40     as1t    speed right right 1.4303812  0 Inf  0 Inf
      81      41     as1t  neutral right right 2.6435485  0 Inf  0 Inf
      83      42     as1t  neutral right right 2.3567246  0 Inf  0 Inf
      85      43     as1t  neutral right right 1.5802799  0 Inf  0 Inf
      87      44     as1t  neutral right  left 1.7137110  0 Inf  0 Inf
      89      45     as1t  neutral right  left 2.1747273  0 Inf  0 Inf
      91      46     as1t  neutral right right 1.8815617  0 Inf  0 Inf
      93      47     as1t  neutral right right 1.5957860  0 Inf  0 Inf
      95      48     as1t  neutral right  left 1.6588285  0 Inf  0 Inf
      97      49     as1t  neutral right right 2.4957353  0 Inf  0 Inf
      99      50     as1t  neutral right right 2.0791869  0 Inf  0 Inf
      101     51     as1t accuracy right right 1.5873774  0 Inf  0 Inf
      103     52     as1t accuracy right right 1.4924681  0 Inf  0 Inf
      105     53     as1t accuracy right right 1.5113153  0 Inf  0 Inf
      107     54     as1t accuracy right right 1.8164823  0 Inf  0 Inf
      109     55     as1t accuracy right right 1.3439070  0 Inf  0 Inf
      111     56     as1t accuracy right  left 1.0432531  0 Inf  0 Inf
      113     57     as1t accuracy right  left 1.5194730  0 Inf  0 Inf
      115     58     as1t accuracy right right 2.0471897  0 Inf  0 Inf
      117     59     as1t accuracy right right 1.9632373  0 Inf  0 Inf
      119     60     as1t accuracy right right 2.0072728  0 Inf  0 Inf
      121      1     bd6t    speed  left  left 1.6003743  0 Inf  0 Inf
      123      2     bd6t    speed  left  left 1.2668039  0 Inf  0 Inf
      125      3     bd6t    speed  left  left 1.1422959  0 Inf  0 Inf
      127      4     bd6t    speed  left  left 1.1051867  0 Inf  0 Inf
      129      5     bd6t    speed  left  left 1.3786669  0 Inf  0 Inf
      131      6     bd6t    speed  left  left 1.5071506  0 Inf  0 Inf
      133      7     bd6t    speed  left  left 1.3780718  0 Inf  0 Inf
      135      8     bd6t    speed  left  left 1.2390632  0 Inf  0 Inf
      137      9     bd6t    speed  left  left 1.4961801  0 Inf  0 Inf
      139     10     bd6t    speed  left  left 1.0043952  0 Inf  0 Inf
      141     11     bd6t  neutral  left right 2.3222674  0 Inf  0 Inf
      143     12     bd6t  neutral  left  left 1.3377450  0 Inf  0 Inf
      145     13     bd6t  neutral  left right 1.5519956  0 Inf  0 Inf
      147     14     bd6t  neutral  left right 1.6327317  0 Inf  0 Inf
      149     15     bd6t  neutral  left right 1.8546351  0 Inf  0 Inf
      151     16     bd6t  neutral  left right 1.7403799  0 Inf  0 Inf
      153     17     bd6t  neutral  left right 1.5687975  0 Inf  0 Inf
      155     18     bd6t  neutral  left  left 1.9929754  0 Inf  0 Inf
      157     19     bd6t  neutral  left  left 2.0723467  0 Inf  0 Inf
      159     20     bd6t  neutral  left  left 1.6356848  0 Inf  0 Inf
      161     21     bd6t accuracy  left  left 1.2839819  0 Inf  0 Inf
      163     22     bd6t accuracy  left right 2.0638401  0 Inf  0 Inf
      165     23     bd6t accuracy  left right 1.3536935  0 Inf  0 Inf
      167     24     bd6t accuracy  left  left 1.6755419  0 Inf  0 Inf
      169     25     bd6t accuracy  left  left 1.8506097  0 Inf  0 Inf
      171     26     bd6t accuracy  left  left 2.1221108  0 Inf  0 Inf
      173     27     bd6t accuracy  left right 1.2298956  0 Inf  0 Inf
      175     28     bd6t accuracy  left  left 2.2236751  0 Inf  0 Inf
      177     29     bd6t accuracy  left  left 1.3812292  0 Inf  0 Inf
      179     30     bd6t accuracy  left  left 1.6827496  0 Inf  0 Inf
      181     31     bd6t    speed right right 1.3194418  0 Inf  0 Inf
      183     32     bd6t    speed right right 1.2071600  0 Inf  0 Inf
      185     33     bd6t    speed right right 1.0587408  0 Inf  0 Inf
      187     34     bd6t    speed right  left 1.1781135  0 Inf  0 Inf
      189     35     bd6t    speed right  left 1.4797612  0 Inf  0 Inf
      191     36     bd6t    speed right right 1.1431422  0 Inf  0 Inf
      193     37     bd6t    speed right right 1.2590853  0 Inf  0 Inf
      195     38     bd6t    speed right right 1.1011947  0 Inf  0 Inf
      197     39     bd6t    speed right right 1.2614584  0 Inf  0 Inf
      199     40     bd6t    speed right right 1.1919023  0 Inf  0 Inf
      201     41     bd6t  neutral right right 1.3749543  0 Inf  0 Inf
      203     42     bd6t  neutral right right 2.2333520  0 Inf  0 Inf
      205     43     bd6t  neutral right  left 3.8657458  0 Inf  0 Inf
      207     44     bd6t  neutral right right 1.6643480  0 Inf  0 Inf
      209     45     bd6t  neutral right right 2.0649714  0 Inf  0 Inf
      211     46     bd6t  neutral right right 3.0167270  0 Inf  0 Inf
      213     47     bd6t  neutral right right 1.9764564  0 Inf  0 Inf
      215     48     bd6t  neutral right right 2.7169333  0 Inf  0 Inf
      217     49     bd6t  neutral right  left 1.6947490  0 Inf  0 Inf
      219     50     bd6t  neutral right right 1.6200159  0 Inf  0 Inf
      221     51     bd6t accuracy right right 1.2978131  0 Inf  0 Inf
      223     52     bd6t accuracy right  left 1.3096785  0 Inf  0 Inf
      225     53     bd6t accuracy right right 2.4621440  0 Inf  0 Inf
      227     54     bd6t accuracy right  left 1.6115293  0 Inf  0 Inf
      229     55     bd6t accuracy right  left 1.1581419  0 Inf  0 Inf
      231     56     bd6t accuracy right right 1.5756386  0 Inf  0 Inf
      233     57     bd6t accuracy right  left 1.3103896  0 Inf  0 Inf
      235     58     bd6t accuracy right right 1.5694869  0 Inf  0 Inf
      237     59     bd6t accuracy right right 2.2055039  0 Inf  0 Inf
      239     60     bd6t accuracy right right 1.6128635  0 Inf  0 Inf

# RDM

    Code
      init_chains(RDM_s, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $theta_mu
                         [,1]
      v            0.10295498
      v_lMd       -0.17485347
      s_lMTRUE     0.15362521
      B            0.04786557
      B_Eneutral  -0.11250747
      B_Eaccuracy  0.47721481
      B_lRright    0.41720870
      A           -0.03732249
      t0          -0.77359023
      
      $theta_var
      , , 1
      
                              v        v_lMd      s_lMTRUE            B   B_Eneutral
      v            0.0259490148  0.005010764  0.0008039212  0.001039186  0.004671743
      v_lMd        0.0050107643  0.054198668 -0.0063395659  0.016853096  0.011241011
      s_lMTRUE     0.0008039212 -0.006339566  0.0660491934 -0.021845123  0.024543796
      B            0.0010391859  0.016853096 -0.0218451231  0.058738576  0.005181317
      B_Eneutral   0.0046717429  0.011241011  0.0245437963  0.005181317  0.057368805
      B_Eaccuracy  0.0062510813 -0.002114676 -0.0040508422  0.004327501 -0.002720295
      B_lRright    0.0031684940  0.003367429  0.0059128710  0.014675459  0.019979670
      A           -0.0131009964 -0.029069087  0.0167338768 -0.022679508 -0.005251711
      t0           0.0040141996  0.016844431 -0.0079536906  0.011913763  0.008338580
                   B_Eaccuracy    B_lRright            A           t0
      v            0.006251081  0.003168494 -0.013100996  0.004014200
      v_lMd       -0.002114676  0.003367429 -0.029069087  0.016844431
      s_lMTRUE    -0.004050842  0.005912871  0.016733877 -0.007953691
      B            0.004327501  0.014675459 -0.022679508  0.011913763
      B_Eneutral  -0.002720295  0.019979670 -0.005251711  0.008338580
      B_Eaccuracy  0.029195495  0.005568241 -0.009362256 -0.002326190
      B_lRright    0.005568241  0.037848420 -0.001496246  0.005725596
      A           -0.009362256 -0.001496246  0.078727638 -0.016972357
      t0          -0.002326190  0.005725596 -0.016972357  0.044669825
      
      
      $a_half
                       [,1]
      v           1.3638374
      v_lMd       0.3851234
      s_lMTRUE    0.4683564
      B           0.8473333
      B_Eneutral  0.1829583
      B_Eaccuracy 0.4256923
      B_lRright   0.2080372
      A           0.8331091
      t0          1.0314865
      
      $alpha
      , , 1
      
                         as1t       bd6t
      v            0.04862164  0.2351246
      v_lMd       -0.59169531 -0.3027253
      s_lMTRUE     0.19585407  0.2118511
      B           -0.03412426 -0.1551697
      B_Eneutral  -0.20813009 -0.2107319
      B_Eaccuracy  0.40342201  0.3663802
      B_lRright    0.09966816  0.4408931
      A            0.58638382  0.4656651
      t0          -0.97303424 -0.9267676
      
      
      $stage
      [1] "init"
      
      $subj_ll
                 [,1]
      as1t  -7991.728
      bd6t -10390.021
      
      $last_theta_var_inv
                  [,1]       [,2]        [,3]        [,4]       [,5]       [,6]
       [1,] 44.9040894 -0.8014289  -0.9812410   3.3620715  -2.191829 -8.1337372
       [2,] -0.8014289 26.2380707   0.2153557  -3.8519366  -3.769284  4.0043411
       [3,] -0.9812410  0.2153557  23.0950016   8.2391363 -10.580642  1.0256725
       [4,]  3.3620715 -3.8519366   8.2391363  25.6470100  -1.349849 -0.6209953
       [5,] -2.1918285 -3.7692843 -10.5806416  -1.3498493  27.998252  4.1803956
       [6,] -8.1337372  4.0043411   1.0256725  -0.6209953   4.180396 40.5675700
       [7,] -1.9141836  1.6771437  -1.8938780 -10.0302584 -12.222883 -8.1641299
       [8,]  6.7725554  7.5357230  -2.5434395   4.1758907   1.649528  5.5506770
       [9,] -2.0002672 -5.1959610   3.2265984  -1.1309205  -2.720966  4.0568270
                  [,7]      [,8]      [,9]
       [1,]  -1.914184  6.772555 -2.000267
       [2,]   1.677144  7.535723 -5.195961
       [3,]  -1.893878 -2.543439  3.226598
       [4,] -10.030258  4.175891 -1.130920
       [5,] -12.222883  1.649528 -2.720966
       [6,]  -8.164130  5.550677  4.056827
       [7,]  38.521104 -3.812009 -2.651828
       [8,]  -3.812009 19.684888  2.932193
       [9,]  -2.651828  2.932193 27.574874
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_RDM, design_RDM, n_trials = 10)
    Output
          trials subjects        E     S     R        rt LT  UT LC  UC
      1        1     as1t    speed  left  left 1.3918623  0 Inf  0 Inf
      3        2     as1t    speed  left  left 1.0688901  0 Inf  0 Inf
      5        3     as1t    speed  left  left 0.9512526  0 Inf  0 Inf
      7        4     as1t    speed  left  left 0.8309385  0 Inf  0 Inf
      9        5     as1t    speed  left  left 0.8944488  0 Inf  0 Inf
      11       6     as1t    speed  left  left 0.8511627  0 Inf  0 Inf
      13       7     as1t    speed  left  left 1.1662794  0 Inf  0 Inf
      15       8     as1t    speed  left  left 1.1627099  0 Inf  0 Inf
      17       9     as1t    speed  left  left 0.9653640  0 Inf  0 Inf
      19      10     as1t    speed  left  left 1.0248635  0 Inf  0 Inf
      21      11     as1t  neutral  left  left 1.4453645  0 Inf  0 Inf
      23      12     as1t  neutral  left  left 1.6061440  0 Inf  0 Inf
      25      13     as1t  neutral  left  left 1.7985809  0 Inf  0 Inf
      27      14     as1t  neutral  left  left 1.4021645  0 Inf  0 Inf
      29      15     as1t  neutral  left  left 1.7208506  0 Inf  0 Inf
      31      16     as1t  neutral  left  left 1.4408238  0 Inf  0 Inf
      33      17     as1t  neutral  left  left 1.3571863  0 Inf  0 Inf
      35      18     as1t  neutral  left  left 1.5092219  0 Inf  0 Inf
      37      19     as1t  neutral  left  left 1.4253424  0 Inf  0 Inf
      39      20     as1t  neutral  left  left 1.8877953  0 Inf  0 Inf
      41      21     as1t accuracy  left  left 1.1156122  0 Inf  0 Inf
      43      22     as1t accuracy  left  left 1.4540218  0 Inf  0 Inf
      45      23     as1t accuracy  left  left 1.4829104  0 Inf  0 Inf
      47      24     as1t accuracy  left  left 1.2635398  0 Inf  0 Inf
      49      25     as1t accuracy  left  left 1.0630212  0 Inf  0 Inf
      51      26     as1t accuracy  left  left 1.3051940  0 Inf  0 Inf
      53      27     as1t accuracy  left  left 1.3294166  0 Inf  0 Inf
      55      28     as1t accuracy  left right 1.3184059  0 Inf  0 Inf
      57      29     as1t accuracy  left  left 1.4759064  0 Inf  0 Inf
      59      30     as1t accuracy  left  left 1.2388894  0 Inf  0 Inf
      61      31     as1t    speed right right 0.8283962  0 Inf  0 Inf
      63      32     as1t    speed right right 0.9386218  0 Inf  0 Inf
      65      33     as1t    speed right right 0.9937264  0 Inf  0 Inf
      67      34     as1t    speed right right 0.8586215  0 Inf  0 Inf
      69      35     as1t    speed right right 1.2849061  0 Inf  0 Inf
      71      36     as1t    speed right right 1.1000486  0 Inf  0 Inf
      73      37     as1t    speed right right 0.7752526  0 Inf  0 Inf
      75      38     as1t    speed right right 0.8087958  0 Inf  0 Inf
      77      39     as1t    speed right right 1.1313398  0 Inf  0 Inf
      79      40     as1t    speed right right 1.0757195  0 Inf  0 Inf
      81      41     as1t  neutral right  left 1.3512642  0 Inf  0 Inf
      83      42     as1t  neutral right right 1.4729514  0 Inf  0 Inf
      85      43     as1t  neutral right right 1.4667028  0 Inf  0 Inf
      87      44     as1t  neutral right right 1.4469030  0 Inf  0 Inf
      89      45     as1t  neutral right right 1.5636693  0 Inf  0 Inf
      91      46     as1t  neutral right right 1.7366819  0 Inf  0 Inf
      93      47     as1t  neutral right right 1.4481670  0 Inf  0 Inf
      95      48     as1t  neutral right right 1.4824821  0 Inf  0 Inf
      97      49     as1t  neutral right right 1.6676795  0 Inf  0 Inf
      99      50     as1t  neutral right right 1.2922862  0 Inf  0 Inf
      101     51     as1t accuracy right right 1.2565651  0 Inf  0 Inf
      103     52     as1t accuracy right right 1.3922900  0 Inf  0 Inf
      105     53     as1t accuracy right right 0.9328332  0 Inf  0 Inf
      107     54     as1t accuracy right right 1.5982685  0 Inf  0 Inf
      109     55     as1t accuracy right right 1.2804435  0 Inf  0 Inf
      111     56     as1t accuracy right right 1.1191248  0 Inf  0 Inf
      113     57     as1t accuracy right right 1.6224066  0 Inf  0 Inf
      115     58     as1t accuracy right right 1.1666073  0 Inf  0 Inf
      117     59     as1t accuracy right right 1.3193807  0 Inf  0 Inf
      119     60     as1t accuracy right right 1.5705005  0 Inf  0 Inf
      121      1     bd6t    speed  left  left 0.8401741  0 Inf  0 Inf
      123      2     bd6t    speed  left right 1.0684374  0 Inf  0 Inf
      125      3     bd6t    speed  left  left 1.3242524  0 Inf  0 Inf
      127      4     bd6t    speed  left  left 1.2068959  0 Inf  0 Inf
      129      5     bd6t    speed  left right 1.0586767  0 Inf  0 Inf
      131      6     bd6t    speed  left  left 1.5112399  0 Inf  0 Inf
      133      7     bd6t    speed  left  left 1.3361156  0 Inf  0 Inf
      135      8     bd6t    speed  left  left 0.9845948  0 Inf  0 Inf
      137      9     bd6t    speed  left  left 1.1412832  0 Inf  0 Inf
      139     10     bd6t    speed  left  left 0.7257711  0 Inf  0 Inf
      141     11     bd6t  neutral  left  left 1.7248749  0 Inf  0 Inf
      143     12     bd6t  neutral  left  left 1.2763894  0 Inf  0 Inf
      145     13     bd6t  neutral  left  left 1.5615051  0 Inf  0 Inf
      147     14     bd6t  neutral  left right 1.5348809  0 Inf  0 Inf
      149     15     bd6t  neutral  left  left 1.5772825  0 Inf  0 Inf
      151     16     bd6t  neutral  left  left 1.2166392  0 Inf  0 Inf
      153     17     bd6t  neutral  left  left 1.5141789  0 Inf  0 Inf
      155     18     bd6t  neutral  left  left 1.6731905  0 Inf  0 Inf
      157     19     bd6t  neutral  left  left 1.3833226  0 Inf  0 Inf
      159     20     bd6t  neutral  left  left 1.2214668  0 Inf  0 Inf
      161     21     bd6t accuracy  left  left 1.6739206  0 Inf  0 Inf
      163     22     bd6t accuracy  left  left 1.2534416  0 Inf  0 Inf
      165     23     bd6t accuracy  left  left 1.0997884  0 Inf  0 Inf
      167     24     bd6t accuracy  left  left 1.2228005  0 Inf  0 Inf
      169     25     bd6t accuracy  left  left 1.1551982  0 Inf  0 Inf
      171     26     bd6t accuracy  left  left 1.1270968  0 Inf  0 Inf
      173     27     bd6t accuracy  left  left 1.3991320  0 Inf  0 Inf
      175     28     bd6t accuracy  left  left 1.1191917  0 Inf  0 Inf
      177     29     bd6t accuracy  left right 0.9213955  0 Inf  0 Inf
      179     30     bd6t accuracy  left  left 1.1606488  0 Inf  0 Inf
      181     31     bd6t    speed right right 1.1421160  0 Inf  0 Inf
      183     32     bd6t    speed right right 0.8012760  0 Inf  0 Inf
      185     33     bd6t    speed right right 1.0507990  0 Inf  0 Inf
      187     34     bd6t    speed right right 0.9727477  0 Inf  0 Inf
      189     35     bd6t    speed right right 1.4838783  0 Inf  0 Inf
      191     36     bd6t    speed right right 1.0686530  0 Inf  0 Inf
      193     37     bd6t    speed right right 0.7512505  0 Inf  0 Inf
      195     38     bd6t    speed right  left 0.8442915  0 Inf  0 Inf
      197     39     bd6t    speed right  left 0.6995642  0 Inf  0 Inf
      199     40     bd6t    speed right right 0.9171390  0 Inf  0 Inf
      201     41     bd6t  neutral right right 1.5097716  0 Inf  0 Inf
      203     42     bd6t  neutral right right 1.1157394  0 Inf  0 Inf
      205     43     bd6t  neutral right right 1.3224846  0 Inf  0 Inf
      207     44     bd6t  neutral right right 1.2412219  0 Inf  0 Inf
      209     45     bd6t  neutral right right 1.8140746  0 Inf  0 Inf
      211     46     bd6t  neutral right right 1.5120106  0 Inf  0 Inf
      213     47     bd6t  neutral right right 1.2495202  0 Inf  0 Inf
      215     48     bd6t  neutral right  left 1.4333431  0 Inf  0 Inf
      217     49     bd6t  neutral right right 1.6915711  0 Inf  0 Inf
      219     50     bd6t  neutral right right 1.6786080  0 Inf  0 Inf
      221     51     bd6t accuracy right right 1.2911611  0 Inf  0 Inf
      223     52     bd6t accuracy right  left 1.0626044  0 Inf  0 Inf
      225     53     bd6t accuracy right right 1.1362981  0 Inf  0 Inf
      227     54     bd6t accuracy right right 0.8687058  0 Inf  0 Inf
      229     55     bd6t accuracy right right 1.2517573  0 Inf  0 Inf
      231     56     bd6t accuracy right right 1.1841919  0 Inf  0 Inf
      233     57     bd6t accuracy right right 1.3515198  0 Inf  0 Inf
      235     58     bd6t accuracy right right 1.0408941  0 Inf  0 Inf
      237     59     bd6t accuracy right right 1.1624195  0 Inf  0 Inf
      239     60     bd6t accuracy right right 1.4358274  0 Inf  0 Inf

# DDM

    Code
      init_chains(DDM_s, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $theta_mu
                         [,1]
      v_Sleft      0.10295498
      v_Sright    -0.17485347
      a            0.15362521
      a_Eneutral   0.04786557
      a_Eaccuracy -0.11250747
      t0           0.47721481
      Z            0.41720870
      sv          -0.03732249
      SZ          -0.77359023
      
      $theta_var
      , , 1
      
                        v_Sleft     v_Sright             a   a_Eneutral  a_Eaccuracy
      v_Sleft      0.0259490148  0.005010764  0.0008039212  0.001039186  0.004671743
      v_Sright     0.0050107643  0.054198668 -0.0063395659  0.016853096  0.011241011
      a            0.0008039212 -0.006339566  0.0660491934 -0.021845123  0.024543796
      a_Eneutral   0.0010391859  0.016853096 -0.0218451231  0.058738576  0.005181317
      a_Eaccuracy  0.0046717429  0.011241011  0.0245437963  0.005181317  0.057368805
      t0           0.0062510813 -0.002114676 -0.0040508422  0.004327501 -0.002720295
      Z            0.0031684940  0.003367429  0.0059128710  0.014675459  0.019979670
      sv          -0.0131009964 -0.029069087  0.0167338768 -0.022679508 -0.005251711
      SZ           0.0040141996  0.016844431 -0.0079536906  0.011913763  0.008338580
                            t0            Z           sv           SZ
      v_Sleft      0.006251081  0.003168494 -0.013100996  0.004014200
      v_Sright    -0.002114676  0.003367429 -0.029069087  0.016844431
      a           -0.004050842  0.005912871  0.016733877 -0.007953691
      a_Eneutral   0.004327501  0.014675459 -0.022679508  0.011913763
      a_Eaccuracy -0.002720295  0.019979670 -0.005251711  0.008338580
      t0           0.029195495  0.005568241 -0.009362256 -0.002326190
      Z            0.005568241  0.037848420 -0.001496246  0.005725596
      sv          -0.009362256 -0.001496246  0.078727638 -0.016972357
      SZ          -0.002326190  0.005725596 -0.016972357  0.044669825
      
      
      $a_half
                       [,1]
      v_Sleft     1.3638374
      v_Sright    0.3851234
      a           0.4683564
      a_Eneutral  0.8473333
      a_Eaccuracy 0.1829583
      t0          0.4256923
      Z           0.2080372
      sv          0.8331091
      SZ          1.0314865
      
      $alpha
      , , 1
      
                          as1t        bd6t
      v_Sleft      0.001153243  0.11584132
      v_Sright    -0.739398191  0.27169948
      a           -0.287978912  0.54288903
      a_Eneutral   0.207325046  0.15432663
      a_Eaccuracy -0.719387674  0.03189456
      t0           0.717560487  0.31347144
      Z            0.338818498  0.12447050
      sv           0.206385123 -0.50800449
      SZ          -1.043411727 -0.81921406
      
      
      $stage
      [1] "init"
      
      $subj_ll
                [,1]
      as1t -18650.94
      bd6t -19536.69
      
      $last_theta_var_inv
                  [,1]       [,2]        [,3]        [,4]       [,5]       [,6]
       [1,] 44.9040894 -0.8014289  -0.9812410   3.3620715  -2.191829 -8.1337372
       [2,] -0.8014289 26.2380707   0.2153557  -3.8519366  -3.769284  4.0043411
       [3,] -0.9812410  0.2153557  23.0950016   8.2391363 -10.580642  1.0256725
       [4,]  3.3620715 -3.8519366   8.2391363  25.6470100  -1.349849 -0.6209953
       [5,] -2.1918285 -3.7692843 -10.5806416  -1.3498493  27.998252  4.1803956
       [6,] -8.1337372  4.0043411   1.0256725  -0.6209953   4.180396 40.5675700
       [7,] -1.9141836  1.6771437  -1.8938780 -10.0302584 -12.222883 -8.1641299
       [8,]  6.7725554  7.5357230  -2.5434395   4.1758907   1.649528  5.5506770
       [9,] -2.0002672 -5.1959610   3.2265984  -1.1309205  -2.720966  4.0568270
                  [,7]      [,8]      [,9]
       [1,]  -1.914184  6.772555 -2.000267
       [2,]   1.677144  7.535723 -5.195961
       [3,]  -1.893878 -2.543439  3.226598
       [4,] -10.030258  4.175891 -1.130920
       [5,] -12.222883  1.649528 -2.720966
       [6,]  -8.164130  5.550677  4.056827
       [7,]  38.521104 -3.812009 -2.651828
       [8,]  -3.812009 19.684888  2.932193
       [9,]  -2.651828  2.932193 27.574874
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_DDM, design_DDM, n_trials = 10)
    Output
          trials subjects        E     S     R        rt LT  UT LC  UC
      1        1     as1t    speed  left right 0.8292741  0 Inf  0 Inf
      2        2     as1t    speed  left right 1.1244544  0 Inf  0 Inf
      3        3     as1t    speed  left right 0.5102559  0 Inf  0 Inf
      4        4     as1t    speed  left  left 0.2624711  0 Inf  0 Inf
      5        5     as1t    speed  left  left 0.4860668  0 Inf  0 Inf
      6        6     as1t    speed  left  left 0.3080312  0 Inf  0 Inf
      7        7     as1t    speed  left  left 0.2831042  0 Inf  0 Inf
      8        8     as1t    speed  left  left 0.4433934  0 Inf  0 Inf
      9        9     as1t    speed  left right 0.3684199  0 Inf  0 Inf
      10      10     as1t    speed  left  left 0.7706748  0 Inf  0 Inf
      21      11     as1t  neutral  left  left 2.1807246  0 Inf  0 Inf
      22      12     as1t  neutral  left  left 0.3915126  0 Inf  0 Inf
      23      13     as1t  neutral  left  left 0.5574094  0 Inf  0 Inf
      24      14     as1t  neutral  left  left 0.5282457  0 Inf  0 Inf
      25      15     as1t  neutral  left  left 0.9902447  0 Inf  0 Inf
      26      16     as1t  neutral  left  left 0.3662669  0 Inf  0 Inf
      27      17     as1t  neutral  left  left 0.4342010  0 Inf  0 Inf
      28      18     as1t  neutral  left  left 0.5744900  0 Inf  0 Inf
      29      19     as1t  neutral  left right 0.9113978  0 Inf  0 Inf
      30      20     as1t  neutral  left  left 1.1391006  0 Inf  0 Inf
      41      21     as1t accuracy  left  left 1.4771946  0 Inf  0 Inf
      42      22     as1t accuracy  left  left 0.5487755  0 Inf  0 Inf
      43      23     as1t accuracy  left  left 0.9717475  0 Inf  0 Inf
      44      24     as1t accuracy  left right 1.4863021  0 Inf  0 Inf
      45      25     as1t accuracy  left  left 0.5843454  0 Inf  0 Inf
      46      26     as1t accuracy  left  left 0.5242357  0 Inf  0 Inf
      47      27     as1t accuracy  left  left 0.6869696  0 Inf  0 Inf
      48      28     as1t accuracy  left  left 1.0398874  0 Inf  0 Inf
      49      29     as1t accuracy  left  left 1.2009965  0 Inf  0 Inf
      50      30     as1t accuracy  left  left 0.6043799  0 Inf  0 Inf
      61      31     as1t    speed right right 0.2997463  0 Inf  0 Inf
      62      32     as1t    speed right right 0.2892761  0 Inf  0 Inf
      63      33     as1t    speed right  left 0.3728619  0 Inf  0 Inf
      64      34     as1t    speed right right 0.3861277  0 Inf  0 Inf
      65      35     as1t    speed right right 0.2710263  0 Inf  0 Inf
      66      36     as1t    speed right right 1.4417330  0 Inf  0 Inf
      67      37     as1t    speed right right 0.3978968  0 Inf  0 Inf
      68      38     as1t    speed right right 0.3022248  0 Inf  0 Inf
      69      39     as1t    speed right right 0.5108373  0 Inf  0 Inf
      70      40     as1t    speed right right 0.4101871  0 Inf  0 Inf
      81      41     as1t  neutral right right 0.2941752  0 Inf  0 Inf
      82      42     as1t  neutral right  left 0.9465992  0 Inf  0 Inf
      83      43     as1t  neutral right right 0.6176995  0 Inf  0 Inf
      84      44     as1t  neutral right right 0.4293582  0 Inf  0 Inf
      85      45     as1t  neutral right right 0.4141462  0 Inf  0 Inf
      86      46     as1t  neutral right right 0.3249108  0 Inf  0 Inf
      87      47     as1t  neutral right  left 0.5610629  0 Inf  0 Inf
      88      48     as1t  neutral right right 0.3906846  0 Inf  0 Inf
      89      49     as1t  neutral right  left 0.8440436  0 Inf  0 Inf
      90      50     as1t  neutral right right 0.5268619  0 Inf  0 Inf
      101     51     as1t accuracy right right 0.5435503  0 Inf  0 Inf
      102     52     as1t accuracy right right 1.5774766  0 Inf  0 Inf
      103     53     as1t accuracy right right 1.4329561  0 Inf  0 Inf
      104     54     as1t accuracy right right 0.3647529  0 Inf  0 Inf
      105     55     as1t accuracy right right 0.7885833  0 Inf  0 Inf
      106     56     as1t accuracy right right 2.5388481  0 Inf  0 Inf
      107     57     as1t accuracy right  left 3.3845968  0 Inf  0 Inf
      108     58     as1t accuracy right right 0.4596461  0 Inf  0 Inf
      109     59     as1t accuracy right right 0.5804869  0 Inf  0 Inf
      110     60     as1t accuracy right right 0.9339778  0 Inf  0 Inf
      11       1     bd6t    speed  left  left 0.5113699  0 Inf  0 Inf
      12       2     bd6t    speed  left  left 0.8324528  0 Inf  0 Inf
      13       3     bd6t    speed  left  left 0.2857536  0 Inf  0 Inf
      14       4     bd6t    speed  left  left 1.6465173  0 Inf  0 Inf
      15       5     bd6t    speed  left right 0.4440789  0 Inf  0 Inf
      16       6     bd6t    speed  left  left 0.5638909  0 Inf  0 Inf
      17       7     bd6t    speed  left right 0.3872566  0 Inf  0 Inf
      18       8     bd6t    speed  left  left 0.4102117  0 Inf  0 Inf
      19       9     bd6t    speed  left  left 0.4385320  0 Inf  0 Inf
      20      10     bd6t    speed  left  left 0.5151820  0 Inf  0 Inf
      31      11     bd6t  neutral  left  left 0.6754373  0 Inf  0 Inf
      32      12     bd6t  neutral  left  left 0.5341780  0 Inf  0 Inf
      33      13     bd6t  neutral  left  left 0.6738859  0 Inf  0 Inf
      34      14     bd6t  neutral  left right 0.4096221  0 Inf  0 Inf
      35      15     bd6t  neutral  left right 0.6590889  0 Inf  0 Inf
      36      16     bd6t  neutral  left  left 0.7598441  0 Inf  0 Inf
      37      17     bd6t  neutral  left  left 0.7733685  0 Inf  0 Inf
      38      18     bd6t  neutral  left  left 3.0105664  0 Inf  0 Inf
      39      19     bd6t  neutral  left  left 1.2588182  0 Inf  0 Inf
      40      20     bd6t  neutral  left  left 0.3977863  0 Inf  0 Inf
      51      21     bd6t accuracy  left right 0.9841399  0 Inf  0 Inf
      52      22     bd6t accuracy  left  left 0.4919777  0 Inf  0 Inf
      53      23     bd6t accuracy  left  left 0.4202342  0 Inf  0 Inf
      54      24     bd6t accuracy  left  left 0.8096352  0 Inf  0 Inf
      55      25     bd6t accuracy  left  left 0.5073642  0 Inf  0 Inf
      56      26     bd6t accuracy  left  left 0.3106668  0 Inf  0 Inf
      57      27     bd6t accuracy  left  left 0.4004328  0 Inf  0 Inf
      58      28     bd6t accuracy  left  left 0.6001172  0 Inf  0 Inf
      59      29     bd6t accuracy  left  left 1.2838355  0 Inf  0 Inf
      60      30     bd6t accuracy  left  left 0.9690398  0 Inf  0 Inf
      71      31     bd6t    speed right right 0.5134839  0 Inf  0 Inf
      72      32     bd6t    speed right right 0.3945113  0 Inf  0 Inf
      73      33     bd6t    speed right right 0.6975348  0 Inf  0 Inf
      74      34     bd6t    speed right right 1.0157935  0 Inf  0 Inf
      75      35     bd6t    speed right right 0.4431862  0 Inf  0 Inf
      76      36     bd6t    speed right right 3.1578335  0 Inf  0 Inf
      77      37     bd6t    speed right right 0.8361330  0 Inf  0 Inf
      78      38     bd6t    speed right right 0.6229737  0 Inf  0 Inf
      79      39     bd6t    speed right right 0.5663145  0 Inf  0 Inf
      80      40     bd6t    speed right  left 0.9061478  0 Inf  0 Inf
      91      41     bd6t  neutral right right 1.6592027  0 Inf  0 Inf
      92      42     bd6t  neutral right right 0.4236181  0 Inf  0 Inf
      93      43     bd6t  neutral right right 0.6476295  0 Inf  0 Inf
      94      44     bd6t  neutral right right 0.4088404  0 Inf  0 Inf
      95      45     bd6t  neutral right right 0.7285072  0 Inf  0 Inf
      96      46     bd6t  neutral right right 0.5097442  0 Inf  0 Inf
      97      47     bd6t  neutral right right 0.8411554  0 Inf  0 Inf
      98      48     bd6t  neutral right right 0.4949010  0 Inf  0 Inf
      99      49     bd6t  neutral right right 0.3862297  0 Inf  0 Inf
      100     50     bd6t  neutral right  left 1.7644152  0 Inf  0 Inf
      111     51     bd6t accuracy right  left 0.8200140  0 Inf  0 Inf
      112     52     bd6t accuracy right  left 0.6925654  0 Inf  0 Inf
      113     53     bd6t accuracy right  left 0.8251542  0 Inf  0 Inf
      114     54     bd6t accuracy right right 1.2094080  0 Inf  0 Inf
      115     55     bd6t accuracy right right 0.4710820  0 Inf  0 Inf
      116     56     bd6t accuracy right  left 1.5934566  0 Inf  0 Inf
      117     57     bd6t accuracy right  left 1.5077570  0 Inf  0 Inf
      118     58     bd6t accuracy right  left 0.6312352  0 Inf  0 Inf
      119     59     bd6t accuracy right  left 0.8091109  0 Inf  0 Inf
      120     60     bd6t accuracy right right 0.4777607  0 Inf  0 Inf

