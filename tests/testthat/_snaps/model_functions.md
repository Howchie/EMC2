# LNR

    Code
      init_chains(LNR_s, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $theta_mu
                        [,1]
      m           -0.9685927
      m_lMd        0.7061091
      m_Eneutral   1.4890213
      m_Eaccuracy -1.8150926
      s            0.3304096
      t0          -1.1421557
      
      $theta_var
      , , 1
      
                            m        m_lMd   m_Eneutral m_Eaccuracy           s
      m            0.15548530 -0.120772081  0.018078294  0.04633176 -0.04484492
      m_lMd       -0.12077208  0.184889995 -0.009471693 -0.04763271  0.02012377
      m_Eneutral   0.01807829 -0.009471693  0.058583553 -0.01291727 -0.02728142
      m_Eaccuracy  0.04633176 -0.047632707 -0.012917266  0.09477040 -0.01673469
      s           -0.04484492  0.020123770 -0.027281415 -0.01673469  0.07945096
      t0           0.03785971 -0.018083101 -0.008237056  0.01937868 -0.01285052
                            t0
      m            0.037859710
      m_lMd       -0.018083101
      m_Eneutral  -0.008237056
      m_Eaccuracy  0.019378684
      s           -0.012850519
      t0           0.076944459
      
      
      $a_half
                       [,1]
      m           0.4056687
      m_lMd       1.9319610
      m_Eneutral  0.2242516
      m_Eaccuracy 0.4969591
      s           0.3334446
      t0          0.2464756
      
      $alpha
      , , 1
      
                        as1t         bd6t
      m           -1.7437564 -0.451325041
      m_lMd        0.8894434 -0.004573125
      m_Eneutral   1.3453151  1.767460316
      m_Eaccuracy -1.6999044 -2.108929606
      s            0.8516853  0.442019143
      t0          -1.5816652 -1.285432743
      
      
      $stage
      [1] "init"
      
      $subj_ll
                [,1]
      as1t -1283.287
      bd6t -1121.183
      
      $last_theta_var_inv
                [,1]       [,2]       [,3]      [,4]      [,5]      [,6]
      [1,] 17.928942  9.9292394 -2.6234713 -2.086301  5.400555 -5.341685
      [2,]  9.929239 11.9481645  0.3450765  2.123976  2.803122 -2.107418
      [3,] -2.623471  0.3450765 23.5046084  5.314472  8.257526  3.928791
      [4,] -2.086301  2.1239763  5.3144716 14.075925  2.918336 -0.963037
      [5,]  5.400555  2.8031216  8.2575264  2.918336 18.577442  1.253107
      [6,] -5.341685 -2.1074180  3.9287913 -0.963037  1.253107 16.001843
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_LNR, design_LNR, n_trials = 10)
    Output
          trials subjects        E     S     R        rt LT  UT LC  UC
      1        1     as1t    speed  left right 0.4882775  0 Inf  0 Inf
      3        2     as1t    speed  left  left 0.3759703  0 Inf  0 Inf
      5        3     as1t    speed  left right 0.3549159  0 Inf  0 Inf
      7        4     as1t    speed  left  left 0.7546314  0 Inf  0 Inf
      9        5     as1t    speed  left right 0.4066846  0 Inf  0 Inf
      11       6     as1t    speed  left  left 0.2488817  0 Inf  0 Inf
      13       7     as1t    speed  left right 0.3245281  0 Inf  0 Inf
      15       8     as1t    speed  left right 0.2107385  0 Inf  0 Inf
      17       9     as1t    speed  left right 0.4971083  0 Inf  0 Inf
      19      10     as1t    speed  left right 0.3012542  0 Inf  0 Inf
      21      11     as1t  neutral  left right 0.2512849  0 Inf  0 Inf
      23      12     as1t  neutral  left right 0.5129122  0 Inf  0 Inf
      25      13     as1t  neutral  left right 0.3665155  0 Inf  0 Inf
      27      14     as1t  neutral  left right 0.3394392  0 Inf  0 Inf
      29      15     as1t  neutral  left right 0.4688527  0 Inf  0 Inf
      31      16     as1t  neutral  left  left 0.2423591  0 Inf  0 Inf
      33      17     as1t  neutral  left right 0.2675480  0 Inf  0 Inf
      35      18     as1t  neutral  left right 0.3261101  0 Inf  0 Inf
      37      19     as1t  neutral  left right 0.2754475  0 Inf  0 Inf
      39      20     as1t  neutral  left  left 0.2880217  0 Inf  0 Inf
      41      21     as1t accuracy  left right 0.2223220  0 Inf  0 Inf
      43      22     as1t accuracy  left right 0.2093809  0 Inf  0 Inf
      45      23     as1t accuracy  left right 0.2070884  0 Inf  0 Inf
      47      24     as1t accuracy  left  left 0.6534595  0 Inf  0 Inf
      49      25     as1t accuracy  left right 0.2326188  0 Inf  0 Inf
      51      26     as1t accuracy  left  left 0.4120461  0 Inf  0 Inf
      53      27     as1t accuracy  left  left 0.2932776  0 Inf  0 Inf
      55      28     as1t accuracy  left right 0.3896610  0 Inf  0 Inf
      57      29     as1t accuracy  left right 0.7487846  0 Inf  0 Inf
      59      30     as1t accuracy  left right 0.2187082  0 Inf  0 Inf
      61      31     as1t    speed right right 0.7730265  0 Inf  0 Inf
      63      32     as1t    speed right  left 0.2399823  0 Inf  0 Inf
      65      33     as1t    speed right  left 0.2074797  0 Inf  0 Inf
      67      34     as1t    speed right  left 0.2359580  0 Inf  0 Inf
      69      35     as1t    speed right  left 0.8071387  0 Inf  0 Inf
      71      36     as1t    speed right  left 0.2139748  0 Inf  0 Inf
      73      37     as1t    speed right  left 0.2095997  0 Inf  0 Inf
      75      38     as1t    speed right  left 0.3731836  0 Inf  0 Inf
      77      39     as1t    speed right  left 0.3146687  0 Inf  0 Inf
      79      40     as1t    speed right  left 0.2211629  0 Inf  0 Inf
      81      41     as1t  neutral right  left 0.2966333  0 Inf  0 Inf
      83      42     as1t  neutral right  left 0.6458612  0 Inf  0 Inf
      85      43     as1t  neutral right  left 0.2064647  0 Inf  0 Inf
      87      44     as1t  neutral right  left 0.2932197  0 Inf  0 Inf
      89      45     as1t  neutral right  left 0.6361902  0 Inf  0 Inf
      91      46     as1t  neutral right  left 0.9006914  0 Inf  0 Inf
      93      47     as1t  neutral right  left 0.3468783  0 Inf  0 Inf
      95      48     as1t  neutral right right 0.2191004  0 Inf  0 Inf
      97      49     as1t  neutral right  left 0.6793734  0 Inf  0 Inf
      99      50     as1t  neutral right  left 0.2025083  0 Inf  0 Inf
      101     51     as1t accuracy right  left 0.2200028  0 Inf  0 Inf
      103     52     as1t accuracy right  left 0.4584434  0 Inf  0 Inf
      105     53     as1t accuracy right  left 0.2213869  0 Inf  0 Inf
      107     54     as1t accuracy right  left 0.8649434  0 Inf  0 Inf
      109     55     as1t accuracy right  left 0.2068187  0 Inf  0 Inf
      111     56     as1t accuracy right  left 0.2936895  0 Inf  0 Inf
      113     57     as1t accuracy right  left 0.2394335  0 Inf  0 Inf
      115     58     as1t accuracy right right 1.1782178  0 Inf  0 Inf
      117     59     as1t accuracy right right 0.3448678  0 Inf  0 Inf
      119     60     as1t accuracy right right 0.3024092  0 Inf  0 Inf
      121      1     bd6t    speed  left  left 0.2538393  0 Inf  0 Inf
      123      2     bd6t    speed  left  left 0.8051629  0 Inf  0 Inf
      125      3     bd6t    speed  left right 0.3794030  0 Inf  0 Inf
      127      4     bd6t    speed  left right 0.3124322  0 Inf  0 Inf
      129      5     bd6t    speed  left  left 0.2753972  0 Inf  0 Inf
      131      6     bd6t    speed  left right 0.8570522  0 Inf  0 Inf
      133      7     bd6t    speed  left right 0.2687622  0 Inf  0 Inf
      135      8     bd6t    speed  left right 0.2160403  0 Inf  0 Inf
      137      9     bd6t    speed  left right 0.9311225  0 Inf  0 Inf
      139     10     bd6t    speed  left right 0.2149328  0 Inf  0 Inf
      141     11     bd6t  neutral  left right 0.2860992  0 Inf  0 Inf
      143     12     bd6t  neutral  left right 0.7721988  0 Inf  0 Inf
      145     13     bd6t  neutral  left right 0.2396286  0 Inf  0 Inf
      147     14     bd6t  neutral  left  left 1.1336319  0 Inf  0 Inf
      149     15     bd6t  neutral  left right 0.2736985  0 Inf  0 Inf
      151     16     bd6t  neutral  left right 0.4228318  0 Inf  0 Inf
      153     17     bd6t  neutral  left right 0.2649706  0 Inf  0 Inf
      155     18     bd6t  neutral  left  left 0.3984080  0 Inf  0 Inf
      157     19     bd6t  neutral  left right 0.8836455  0 Inf  0 Inf
      159     20     bd6t  neutral  left right 0.2865105  0 Inf  0 Inf
      161     21     bd6t accuracy  left  left 0.4545996  0 Inf  0 Inf
      163     22     bd6t accuracy  left right 0.5733221  0 Inf  0 Inf
      165     23     bd6t accuracy  left right 0.3093538  0 Inf  0 Inf
      167     24     bd6t accuracy  left right 0.3356477  0 Inf  0 Inf
      169     25     bd6t accuracy  left right 0.2678329  0 Inf  0 Inf
      171     26     bd6t accuracy  left  left 1.8014528  0 Inf  0 Inf
      173     27     bd6t accuracy  left right 0.2141539  0 Inf  0 Inf
      175     28     bd6t accuracy  left right 0.2535888  0 Inf  0 Inf
      177     29     bd6t accuracy  left right 0.2359511  0 Inf  0 Inf
      179     30     bd6t accuracy  left right 0.5717901  0 Inf  0 Inf
      181     31     bd6t    speed right  left 0.2184325  0 Inf  0 Inf
      183     32     bd6t    speed right  left 0.6617118  0 Inf  0 Inf
      185     33     bd6t    speed right  left 0.2278260  0 Inf  0 Inf
      187     34     bd6t    speed right right 0.2331189  0 Inf  0 Inf
      189     35     bd6t    speed right  left 0.3696220  0 Inf  0 Inf
      191     36     bd6t    speed right  left 0.2186500  0 Inf  0 Inf
      193     37     bd6t    speed right  left 0.2268752  0 Inf  0 Inf
      195     38     bd6t    speed right  left 0.3869918  0 Inf  0 Inf
      197     39     bd6t    speed right right 0.3559524  0 Inf  0 Inf
      199     40     bd6t    speed right  left 0.2168651  0 Inf  0 Inf
      201     41     bd6t  neutral right  left 0.2589507  0 Inf  0 Inf
      203     42     bd6t  neutral right  left 0.2379628  0 Inf  0 Inf
      205     43     bd6t  neutral right  left 0.2693013  0 Inf  0 Inf
      207     44     bd6t  neutral right  left 0.3082749  0 Inf  0 Inf
      209     45     bd6t  neutral right  left 0.2324702  0 Inf  0 Inf
      211     46     bd6t  neutral right  left 0.4127583  0 Inf  0 Inf
      213     47     bd6t  neutral right  left 0.8947093  0 Inf  0 Inf
      215     48     bd6t  neutral right  left 0.2914081  0 Inf  0 Inf
      217     49     bd6t  neutral right  left 0.5776475  0 Inf  0 Inf
      219     50     bd6t  neutral right  left 0.6227479  0 Inf  0 Inf
      221     51     bd6t accuracy right right 0.9460963  0 Inf  0 Inf
      223     52     bd6t accuracy right  left 0.2237662  0 Inf  0 Inf
      225     53     bd6t accuracy right  left 0.9027387  0 Inf  0 Inf
      227     54     bd6t accuracy right  left 0.8432729  0 Inf  0 Inf
      229     55     bd6t accuracy right right 0.4870274  0 Inf  0 Inf
      231     56     bd6t accuracy right  left 0.2965542  0 Inf  0 Inf
      233     57     bd6t accuracy right  left 0.9262431  0 Inf  0 Inf
      235     58     bd6t accuracy right  left 0.2087114  0 Inf  0 Inf
      237     59     bd6t accuracy right  left 0.2572218  0 Inf  0 Inf
      239     60     bd6t accuracy right  left 0.2066247  0 Inf  0 Inf

# LBA

    Code
      init_chains(LBA_s, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $theta_mu
                        [,1]
      v           -0.9685927
      v_lMd        0.7061091
      sv_lMTRUE    1.4890213
      B           -1.8150926
      B_Eneutral   0.3304096
      B_Eaccuracy -1.1421557
      B_lRright    0.1571934
      A           -2.0654072
      t0          -0.4405469
      
      $theta_var
      , , 1
      
                             v         v_lMd    sv_lMTRUE             B   B_Eneutral
      v            0.059329585 -0.0104497238 -0.010774234  0.0021450832 -0.018128601
      v_lMd       -0.010449724  0.0559162133 -0.002872401  0.0014127637  0.004728529
      sv_lMTRUE   -0.010774234 -0.0028724006  0.059773323  0.0087906529  0.022576320
      B            0.002145083  0.0014127637  0.008790653  0.0447136854  0.001458968
      B_Eneutral  -0.018128601  0.0047285288  0.022576320  0.0014589675  0.038522063
      B_Eaccuracy  0.018325623  0.0019110134 -0.037503336 -0.0100477306 -0.023423360
      B_lRright   -0.012850266 -0.0030964417  0.020333128  0.0098693855  0.019288722
      A           -0.009137082 -0.0181065189  0.007641907 -0.0029721272  0.008676374
      t0          -0.019762600  0.0006819411  0.001329633  0.0006671515  0.007687211
                   B_Eaccuracy    B_lRright            A            t0
      v            0.018325623 -0.012850266 -0.009137082 -0.0197625997
      v_lMd        0.001911013 -0.003096442 -0.018106519  0.0006819411
      sv_lMTRUE   -0.037503336  0.020333128  0.007641907  0.0013296327
      B           -0.010047731  0.009869385 -0.002972127  0.0006671515
      B_Eneutral  -0.023423360  0.019288722  0.008676374  0.0076872108
      B_Eaccuracy  0.078051779 -0.028702651 -0.004381484  0.0022878527
      B_lRright   -0.028702651  0.042243108  0.005224348 -0.0090916441
      A           -0.004381484  0.005224348  0.052137019  0.0067800486
      t0           0.002287853 -0.009091644  0.006780049  0.0477904171
      
      
      $a_half
                       [,1]
      v           0.4941231
      v_lMd       2.5931756
      sv_lMTRUE   1.4684537
      B           0.8437531
      B_Eneutral  0.6261293
      B_Eaccuracy 0.8155694
      B_lRright   2.2752573
      A           0.7726374
      t0          0.3149199
      
      $alpha
      , , 1
      
                        as1t       bd6t
      v           -1.2748892 -0.8494871
      v_lMd        0.5956117  0.4824282
      sv_lMTRUE    1.3952474  1.2472335
      B           -1.9260735 -2.0320307
      B_Eneutral   0.2458292  0.3543701
      B_Eaccuracy -1.3690513 -0.8626310
      B_lRright    0.3166578  0.3197739
      A           -2.1157142 -1.9422547
      t0          -0.7838674 -0.8355143
      
      
      $stage
      [1] "init"
      
      $subj_ll
                 [,1]
      as1t -10576.397
      bd6t  -9887.026
      
      $last_theta_var_inv
                 [,1]      [,2]       [,3]      [,4]       [,5]      [,6]
       [1,] 26.028515  6.224107 -1.3572503 -3.811330   3.070339 -3.657354
       [2,]  6.224107 22.834612  1.1572385 -1.918507  -6.426290 -1.390907
       [3,] -1.357250  1.157238 27.2973400 -2.546630  -9.560950  9.393468
       [4,] -3.811330 -1.918507 -2.5466296 25.241309   4.568563  1.562002
       [5,]  3.070339 -6.426290 -9.5609496  4.568563  46.584345  3.445213
       [6,] -3.657354 -1.390907  9.3934682  1.562002   3.445213 22.028224
       [7,]  8.055800  5.016142 -2.0467420 -8.084231 -16.170739  7.104309
       [8,]  3.587576  8.795323 -1.3430047  1.224921  -4.769682 -1.611160
       [9,] 11.470411  3.049327 -0.4123182 -4.351590  -8.494190 -1.804338
                    [,7]        [,8]       [,9]
       [1,]   8.05580010  3.58757606 11.4704105
       [2,]   5.01614235  8.79532339  3.0493269
       [3,]  -2.04674198 -1.34300467 -0.4123182
       [4,]  -8.08423140  1.22492096 -4.3515897
       [5,] -16.17073884 -4.76968243 -8.4941901
       [6,]   7.10430863 -1.61116037 -1.8043382
       [7,]  44.63330665 -0.03626134 14.1866918
       [8,]  -0.03626134 23.94575450 -1.1814287
       [9,]  14.18669184 -1.18142869 30.0158974
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_LBA, design_LBA, n_trials = 10)
    Output
          trials subjects        E     S     R        rt LT  UT LC  UC
      1        1     as1t    speed  left  left 1.0419108  0 Inf  0 Inf
      3        2     as1t    speed  left  left 1.1313822  0 Inf  0 Inf
      5        3     as1t    speed  left  left 1.0121655  0 Inf  0 Inf
      7        4     as1t    speed  left right 1.3372172  0 Inf  0 Inf
      9        5     as1t    speed  left  left 1.2873095  0 Inf  0 Inf
      11       6     as1t    speed  left  left 1.0781206  0 Inf  0 Inf
      13       7     as1t    speed  left  left 1.3021941  0 Inf  0 Inf
      15       8     as1t    speed  left right 1.0815568  0 Inf  0 Inf
      17       9     as1t    speed  left  left 1.0943317  0 Inf  0 Inf
      19      10     as1t    speed  left right 1.0432316  0 Inf  0 Inf
      21      11     as1t  neutral  left right 2.8778065  0 Inf  0 Inf
      23      12     as1t  neutral  left  left 1.6457073  0 Inf  0 Inf
      25      13     as1t  neutral  left  left 1.6281656  0 Inf  0 Inf
      27      14     as1t  neutral  left  left 1.8192458  0 Inf  0 Inf
      29      15     as1t  neutral  left  left 1.5957928  0 Inf  0 Inf
      31      16     as1t  neutral  left  left 2.3281691  0 Inf  0 Inf
      33      17     as1t  neutral  left right 1.4379477  0 Inf  0 Inf
      35      18     as1t  neutral  left  left 2.4051755  0 Inf  0 Inf
      37      19     as1t  neutral  left  left 1.5487176  0 Inf  0 Inf
      39      20     as1t  neutral  left  left 2.9673579  0 Inf  0 Inf
      41      21     as1t accuracy  left  left 2.6195623  0 Inf  0 Inf
      43      22     as1t accuracy  left  left 1.5626070  0 Inf  0 Inf
      45      23     as1t accuracy  left right 1.6082487  0 Inf  0 Inf
      47      24     as1t accuracy  left right 1.7777390  0 Inf  0 Inf
      49      25     as1t accuracy  left  left 1.7294351  0 Inf  0 Inf
      51      26     as1t accuracy  left  left 1.5220097  0 Inf  0 Inf
      53      27     as1t accuracy  left right 1.4265891  0 Inf  0 Inf
      55      28     as1t accuracy  left  left 1.8362953  0 Inf  0 Inf
      57      29     as1t accuracy  left  left 1.5278419  0 Inf  0 Inf
      59      30     as1t accuracy  left  left 1.3516826  0 Inf  0 Inf
      61      31     as1t    speed right right 1.7630421  0 Inf  0 Inf
      63      32     as1t    speed right  left 1.2578173  0 Inf  0 Inf
      65      33     as1t    speed right  left 1.4169701  0 Inf  0 Inf
      67      34     as1t    speed right right 1.4811086  0 Inf  0 Inf
      69      35     as1t    speed right right 1.1423523  0 Inf  0 Inf
      71      36     as1t    speed right right 1.2638843  0 Inf  0 Inf
      73      37     as1t    speed right right 1.8524553  0 Inf  0 Inf
      75      38     as1t    speed right right 1.2754474  0 Inf  0 Inf
      77      39     as1t    speed right right 0.8839711  0 Inf  0 Inf
      79      40     as1t    speed right right 1.3630025  0 Inf  0 Inf
      81      41     as1t  neutral right right 2.7865551  0 Inf  0 Inf
      83      42     as1t  neutral right  left 1.7678774  0 Inf  0 Inf
      85      43     as1t  neutral right  left 2.0785731  0 Inf  0 Inf
      87      44     as1t  neutral right right 2.8520288  0 Inf  0 Inf
      89      45     as1t  neutral right right 1.9722206  0 Inf  0 Inf
      91      46     as1t  neutral right right 1.6667824  0 Inf  0 Inf
      93      47     as1t  neutral right right 1.8442542  0 Inf  0 Inf
      95      48     as1t  neutral right  left 1.1390565  0 Inf  0 Inf
      97      49     as1t  neutral right right 1.6982394  0 Inf  0 Inf
      99      50     as1t  neutral right right 1.6889953  0 Inf  0 Inf
      101     51     as1t accuracy right right 1.7240903  0 Inf  0 Inf
      103     52     as1t accuracy right  left 2.1659172  0 Inf  0 Inf
      105     53     as1t accuracy right right 1.4324819  0 Inf  0 Inf
      107     54     as1t accuracy right  left 1.2209504  0 Inf  0 Inf
      109     55     as1t accuracy right right 1.7460738  0 Inf  0 Inf
      111     56     as1t accuracy right right 1.4488193  0 Inf  0 Inf
      113     57     as1t accuracy right right 2.0781390  0 Inf  0 Inf
      115     58     as1t accuracy right right 1.3395594  0 Inf  0 Inf
      117     59     as1t accuracy right  left 1.6660284  0 Inf  0 Inf
      119     60     as1t accuracy right  left 1.2940604  0 Inf  0 Inf
      121      1     bd6t    speed  left  left 1.0564745  0 Inf  0 Inf
      123      2     bd6t    speed  left right 1.6189798  0 Inf  0 Inf
      125      3     bd6t    speed  left right 1.2051955  0 Inf  0 Inf
      127      4     bd6t    speed  left  left 1.2701861  0 Inf  0 Inf
      129      5     bd6t    speed  left right 1.5379538  0 Inf  0 Inf
      131      6     bd6t    speed  left  left 1.0584713  0 Inf  0 Inf
      133      7     bd6t    speed  left  left 1.3707626  0 Inf  0 Inf
      135      8     bd6t    speed  left  left 1.4758537  0 Inf  0 Inf
      137      9     bd6t    speed  left  left 1.1858372  0 Inf  0 Inf
      139     10     bd6t    speed  left right 0.7726886  0 Inf  0 Inf
      141     11     bd6t  neutral  left  left 1.2424591  0 Inf  0 Inf
      143     12     bd6t  neutral  left  left 1.4939337  0 Inf  0 Inf
      145     13     bd6t  neutral  left  left 1.4808143  0 Inf  0 Inf
      147     14     bd6t  neutral  left  left 1.8477564  0 Inf  0 Inf
      149     15     bd6t  neutral  left right 1.7857550  0 Inf  0 Inf
      151     16     bd6t  neutral  left  left 1.6173031  0 Inf  0 Inf
      153     17     bd6t  neutral  left  left 2.4171455  0 Inf  0 Inf
      155     18     bd6t  neutral  left  left 1.3290809  0 Inf  0 Inf
      157     19     bd6t  neutral  left  left 1.8970796  0 Inf  0 Inf
      159     20     bd6t  neutral  left right 1.5519690  0 Inf  0 Inf
      161     21     bd6t accuracy  left  left 2.1936652  0 Inf  0 Inf
      163     22     bd6t accuracy  left  left 1.5336339  0 Inf  0 Inf
      165     23     bd6t accuracy  left  left 1.2778961  0 Inf  0 Inf
      167     24     bd6t accuracy  left right 2.2959055  0 Inf  0 Inf
      169     25     bd6t accuracy  left  left 1.8071073  0 Inf  0 Inf
      171     26     bd6t accuracy  left  left 2.0747441  0 Inf  0 Inf
      173     27     bd6t accuracy  left  left 1.8204269  0 Inf  0 Inf
      175     28     bd6t accuracy  left  left 1.7547447  0 Inf  0 Inf
      177     29     bd6t accuracy  left right 1.3534998  0 Inf  0 Inf
      179     30     bd6t accuracy  left right 1.2378976  0 Inf  0 Inf
      181     31     bd6t    speed right right 1.1159573  0 Inf  0 Inf
      183     32     bd6t    speed right right 1.2882513  0 Inf  0 Inf
      185     33     bd6t    speed right right 1.1563917  0 Inf  0 Inf
      187     34     bd6t    speed right right 1.1198214  0 Inf  0 Inf
      189     35     bd6t    speed right  left 1.0411377  0 Inf  0 Inf
      191     36     bd6t    speed right right 1.4660228  0 Inf  0 Inf
      193     37     bd6t    speed right right 1.0623836  0 Inf  0 Inf
      195     38     bd6t    speed right  left 1.2195528  0 Inf  0 Inf
      197     39     bd6t    speed right  left 1.0540358  0 Inf  0 Inf
      199     40     bd6t    speed right right 1.2974444  0 Inf  0 Inf
      201     41     bd6t  neutral right right 1.7656980  0 Inf  0 Inf
      203     42     bd6t  neutral right right 1.6497181  0 Inf  0 Inf
      205     43     bd6t  neutral right right 1.8682176  0 Inf  0 Inf
      207     44     bd6t  neutral right right 1.7987898  0 Inf  0 Inf
      209     45     bd6t  neutral right right 1.8094139  0 Inf  0 Inf
      211     46     bd6t  neutral right right 2.4099597  0 Inf  0 Inf
      213     47     bd6t  neutral right right 2.2386035  0 Inf  0 Inf
      215     48     bd6t  neutral right right 2.3083533  0 Inf  0 Inf
      217     49     bd6t  neutral right  left 1.2111311  0 Inf  0 Inf
      219     50     bd6t  neutral right right 1.9707261  0 Inf  0 Inf
      221     51     bd6t accuracy right right 1.3309782  0 Inf  0 Inf
      223     52     bd6t accuracy right right 1.1722820  0 Inf  0 Inf
      225     53     bd6t accuracy right right 2.2283015  0 Inf  0 Inf
      227     54     bd6t accuracy right right 2.5522187  0 Inf  0 Inf
      229     55     bd6t accuracy right right 2.1286917  0 Inf  0 Inf
      231     56     bd6t accuracy right right 1.6609632  0 Inf  0 Inf
      233     57     bd6t accuracy right right 1.4865106  0 Inf  0 Inf
      235     58     bd6t accuracy right right 2.7291277  0 Inf  0 Inf
      237     59     bd6t accuracy right  left 1.3126655  0 Inf  0 Inf
      239     60     bd6t accuracy right right 1.4585764  0 Inf  0 Inf

# RDM

    Code
      init_chains(RDM_s, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $theta_mu
                        [,1]
      v           -0.9685927
      v_lMd        0.7061091
      s_lMTRUE     1.4890213
      B           -1.8150926
      B_Eneutral   0.3304096
      B_Eaccuracy -1.1421557
      B_lRright    0.1571934
      A           -2.0654072
      t0          -0.4405469
      
      $theta_var
      , , 1
      
                             v         v_lMd     s_lMTRUE             B   B_Eneutral
      v            0.059329585 -0.0104497238 -0.010774234  0.0021450832 -0.018128601
      v_lMd       -0.010449724  0.0559162133 -0.002872401  0.0014127637  0.004728529
      s_lMTRUE    -0.010774234 -0.0028724006  0.059773323  0.0087906529  0.022576320
      B            0.002145083  0.0014127637  0.008790653  0.0447136854  0.001458968
      B_Eneutral  -0.018128601  0.0047285288  0.022576320  0.0014589675  0.038522063
      B_Eaccuracy  0.018325623  0.0019110134 -0.037503336 -0.0100477306 -0.023423360
      B_lRright   -0.012850266 -0.0030964417  0.020333128  0.0098693855  0.019288722
      A           -0.009137082 -0.0181065189  0.007641907 -0.0029721272  0.008676374
      t0          -0.019762600  0.0006819411  0.001329633  0.0006671515  0.007687211
                   B_Eaccuracy    B_lRright            A            t0
      v            0.018325623 -0.012850266 -0.009137082 -0.0197625997
      v_lMd        0.001911013 -0.003096442 -0.018106519  0.0006819411
      s_lMTRUE    -0.037503336  0.020333128  0.007641907  0.0013296327
      B           -0.010047731  0.009869385 -0.002972127  0.0006671515
      B_Eneutral  -0.023423360  0.019288722  0.008676374  0.0076872108
      B_Eaccuracy  0.078051779 -0.028702651 -0.004381484  0.0022878527
      B_lRright   -0.028702651  0.042243108  0.005224348 -0.0090916441
      A           -0.004381484  0.005224348  0.052137019  0.0067800486
      t0           0.002287853 -0.009091644  0.006780049  0.0477904171
      
      
      $a_half
                       [,1]
      v           0.4941231
      v_lMd       2.5931756
      s_lMTRUE    1.4684537
      B           0.8437531
      B_Eneutral  0.6261293
      B_Eaccuracy 0.8155694
      B_lRright   2.2752573
      A           0.7726374
      t0          0.3149199
      
      $alpha
      , , 1
      
                        as1t       bd6t
      v           -1.2748892 -0.7532267
      v_lMd        0.5956117  0.2242688
      s_lMTRUE     1.3952474  1.3894611
      B           -1.9260735 -1.4897176
      B_Eneutral   0.2458292  0.3957119
      B_Eaccuracy -1.3690513 -1.3320533
      B_lRright    0.3166578  0.3242525
      A           -2.1157142 -2.1016525
      t0          -0.7838674 -0.8287774
      
      
      $stage
      [1] "init"
      
      $subj_ll
                [,1]
      as1t -10097.64
      bd6t -10282.65
      
      $last_theta_var_inv
                 [,1]      [,2]       [,3]      [,4]       [,5]      [,6]
       [1,] 26.028515  6.224107 -1.3572503 -3.811330   3.070339 -3.657354
       [2,]  6.224107 22.834612  1.1572385 -1.918507  -6.426290 -1.390907
       [3,] -1.357250  1.157238 27.2973400 -2.546630  -9.560950  9.393468
       [4,] -3.811330 -1.918507 -2.5466296 25.241309   4.568563  1.562002
       [5,]  3.070339 -6.426290 -9.5609496  4.568563  46.584345  3.445213
       [6,] -3.657354 -1.390907  9.3934682  1.562002   3.445213 22.028224
       [7,]  8.055800  5.016142 -2.0467420 -8.084231 -16.170739  7.104309
       [8,]  3.587576  8.795323 -1.3430047  1.224921  -4.769682 -1.611160
       [9,] 11.470411  3.049327 -0.4123182 -4.351590  -8.494190 -1.804338
                    [,7]        [,8]       [,9]
       [1,]   8.05580010  3.58757606 11.4704105
       [2,]   5.01614235  8.79532339  3.0493269
       [3,]  -2.04674198 -1.34300467 -0.4123182
       [4,]  -8.08423140  1.22492096 -4.3515897
       [5,] -16.17073884 -4.76968243 -8.4941901
       [6,]   7.10430863 -1.61116037 -1.8043382
       [7,]  44.63330665 -0.03626134 14.1866918
       [8,]  -0.03626134 23.94575450 -1.1814287
       [9,]  14.18669184 -1.18142869 30.0158974
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_RDM, design_RDM, n_trials = 10)
    Output
          trials subjects        E     S     R        rt LT  UT LC  UC
      1        1     as1t    speed  left right 1.1034717  0 Inf  0 Inf
      3        2     as1t    speed  left  left 1.3525361  0 Inf  0 Inf
      5        3     as1t    speed  left  left 1.0570703  0 Inf  0 Inf
      7        4     as1t    speed  left  left 1.1516012  0 Inf  0 Inf
      9        5     as1t    speed  left  left 1.2978926  0 Inf  0 Inf
      11       6     as1t    speed  left  left 0.8962148  0 Inf  0 Inf
      13       7     as1t    speed  left  left 0.7928530  0 Inf  0 Inf
      15       8     as1t    speed  left  left 1.1872952  0 Inf  0 Inf
      17       9     as1t    speed  left  left 0.8909292  0 Inf  0 Inf
      19      10     as1t    speed  left right 1.2538669  0 Inf  0 Inf
      21      11     as1t  neutral  left  left 1.4233884  0 Inf  0 Inf
      23      12     as1t  neutral  left  left 1.7676882  0 Inf  0 Inf
      25      13     as1t  neutral  left  left 1.5297115  0 Inf  0 Inf
      27      14     as1t  neutral  left  left 1.7292717  0 Inf  0 Inf
      29      15     as1t  neutral  left  left 1.1853714  0 Inf  0 Inf
      31      16     as1t  neutral  left  left 1.4454578  0 Inf  0 Inf
      33      17     as1t  neutral  left  left 1.4815033  0 Inf  0 Inf
      35      18     as1t  neutral  left  left 1.4792022  0 Inf  0 Inf
      37      19     as1t  neutral  left  left 1.3016436  0 Inf  0 Inf
      39      20     as1t  neutral  left  left 1.5235295  0 Inf  0 Inf
      41      21     as1t accuracy  left  left 1.0580971  0 Inf  0 Inf
      43      22     as1t accuracy  left  left 1.2717121  0 Inf  0 Inf
      45      23     as1t accuracy  left  left 1.1332819  0 Inf  0 Inf
      47      24     as1t accuracy  left  left 1.3346238  0 Inf  0 Inf
      49      25     as1t accuracy  left  left 1.2156782  0 Inf  0 Inf
      51      26     as1t accuracy  left  left 1.2739077  0 Inf  0 Inf
      53      27     as1t accuracy  left  left 0.9279105  0 Inf  0 Inf
      55      28     as1t accuracy  left  left 1.3706267  0 Inf  0 Inf
      57      29     as1t accuracy  left  left 1.1491989  0 Inf  0 Inf
      59      30     as1t accuracy  left  left 1.3015414  0 Inf  0 Inf
      61      31     as1t    speed right right 0.8819235  0 Inf  0 Inf
      63      32     as1t    speed right right 1.0588649  0 Inf  0 Inf
      65      33     as1t    speed right right 1.0610025  0 Inf  0 Inf
      67      34     as1t    speed right right 0.9954035  0 Inf  0 Inf
      69      35     as1t    speed right  left 1.3513162  0 Inf  0 Inf
      71      36     as1t    speed right right 0.8241207  0 Inf  0 Inf
      73      37     as1t    speed right right 0.9660952  0 Inf  0 Inf
      75      38     as1t    speed right right 1.2340796  0 Inf  0 Inf
      77      39     as1t    speed right right 0.7747656  0 Inf  0 Inf
      79      40     as1t    speed right right 1.1975226  0 Inf  0 Inf
      81      41     as1t  neutral right right 1.3617426  0 Inf  0 Inf
      83      42     as1t  neutral right right 1.4225849  0 Inf  0 Inf
      85      43     as1t  neutral right right 1.3528814  0 Inf  0 Inf
      87      44     as1t  neutral right right 1.5540439  0 Inf  0 Inf
      89      45     as1t  neutral right right 1.9168751  0 Inf  0 Inf
      91      46     as1t  neutral right right 1.2206239  0 Inf  0 Inf
      93      47     as1t  neutral right right 1.3801029  0 Inf  0 Inf
      95      48     as1t  neutral right right 1.5767468  0 Inf  0 Inf
      97      49     as1t  neutral right right 1.8313306  0 Inf  0 Inf
      99      50     as1t  neutral right right 2.0640986  0 Inf  0 Inf
      101     51     as1t accuracy right right 1.1725837  0 Inf  0 Inf
      103     52     as1t accuracy right right 1.3322051  0 Inf  0 Inf
      105     53     as1t accuracy right right 1.1438255  0 Inf  0 Inf
      107     54     as1t accuracy right right 1.1675735  0 Inf  0 Inf
      109     55     as1t accuracy right right 1.2315378  0 Inf  0 Inf
      111     56     as1t accuracy right right 0.9599531  0 Inf  0 Inf
      113     57     as1t accuracy right right 1.2737596  0 Inf  0 Inf
      115     58     as1t accuracy right right 1.2413910  0 Inf  0 Inf
      117     59     as1t accuracy right right 1.4681549  0 Inf  0 Inf
      119     60     as1t accuracy right  left 1.5653743  0 Inf  0 Inf
      121      1     bd6t    speed  left  left 0.9827485  0 Inf  0 Inf
      123      2     bd6t    speed  left  left 0.9096272  0 Inf  0 Inf
      125      3     bd6t    speed  left  left 0.8387087  0 Inf  0 Inf
      127      4     bd6t    speed  left  left 0.8558674  0 Inf  0 Inf
      129      5     bd6t    speed  left  left 0.9685154  0 Inf  0 Inf
      131      6     bd6t    speed  left  left 0.9837590  0 Inf  0 Inf
      133      7     bd6t    speed  left  left 0.8939876  0 Inf  0 Inf
      135      8     bd6t    speed  left  left 1.2481454  0 Inf  0 Inf
      137      9     bd6t    speed  left  left 0.8196589  0 Inf  0 Inf
      139     10     bd6t    speed  left  left 1.0853970  0 Inf  0 Inf
      141     11     bd6t  neutral  left  left 1.6229030  0 Inf  0 Inf
      143     12     bd6t  neutral  left  left 1.5817553  0 Inf  0 Inf
      145     13     bd6t  neutral  left  left 1.2124257  0 Inf  0 Inf
      147     14     bd6t  neutral  left  left 1.1611487  0 Inf  0 Inf
      149     15     bd6t  neutral  left  left 1.6600982  0 Inf  0 Inf
      151     16     bd6t  neutral  left  left 1.2348317  0 Inf  0 Inf
      153     17     bd6t  neutral  left  left 0.9798273  0 Inf  0 Inf
      155     18     bd6t  neutral  left  left 1.2604914  0 Inf  0 Inf
      157     19     bd6t  neutral  left  left 1.6356537  0 Inf  0 Inf
      159     20     bd6t  neutral  left  left 1.6713326  0 Inf  0 Inf
      161     21     bd6t accuracy  left  left 1.3448155  0 Inf  0 Inf
      163     22     bd6t accuracy  left  left 1.5054217  0 Inf  0 Inf
      165     23     bd6t accuracy  left  left 1.3345295  0 Inf  0 Inf
      167     24     bd6t accuracy  left  left 1.2569495  0 Inf  0 Inf
      169     25     bd6t accuracy  left  left 1.4043120  0 Inf  0 Inf
      171     26     bd6t accuracy  left  left 1.5902863  0 Inf  0 Inf
      173     27     bd6t accuracy  left  left 1.0782699  0 Inf  0 Inf
      175     28     bd6t accuracy  left  left 1.5920665  0 Inf  0 Inf
      177     29     bd6t accuracy  left  left 1.2827341  0 Inf  0 Inf
      179     30     bd6t accuracy  left  left 1.1715033  0 Inf  0 Inf
      181     31     bd6t    speed right right 1.1717299  0 Inf  0 Inf
      183     32     bd6t    speed right right 0.7868848  0 Inf  0 Inf
      185     33     bd6t    speed right right 0.8653696  0 Inf  0 Inf
      187     34     bd6t    speed right  left 0.9007193  0 Inf  0 Inf
      189     35     bd6t    speed right right 0.8702976  0 Inf  0 Inf
      191     36     bd6t    speed right right 1.2817368  0 Inf  0 Inf
      193     37     bd6t    speed right right 0.8047711  0 Inf  0 Inf
      195     38     bd6t    speed right right 0.9146416  0 Inf  0 Inf
      197     39     bd6t    speed right right 1.1630303  0 Inf  0 Inf
      199     40     bd6t    speed right right 1.0911270  0 Inf  0 Inf
      201     41     bd6t  neutral right right 1.8329683  0 Inf  0 Inf
      203     42     bd6t  neutral right right 1.4730071  0 Inf  0 Inf
      205     43     bd6t  neutral right right 1.4312788  0 Inf  0 Inf
      207     44     bd6t  neutral right right 1.2220378  0 Inf  0 Inf
      209     45     bd6t  neutral right right 1.5538855  0 Inf  0 Inf
      211     46     bd6t  neutral right right 1.3733378  0 Inf  0 Inf
      213     47     bd6t  neutral right right 1.9586144  0 Inf  0 Inf
      215     48     bd6t  neutral right right 2.0643858  0 Inf  0 Inf
      217     49     bd6t  neutral right right 1.3853451  0 Inf  0 Inf
      219     50     bd6t  neutral right  left 1.4190414  0 Inf  0 Inf
      221     51     bd6t accuracy right right 1.2128551  0 Inf  0 Inf
      223     52     bd6t accuracy right right 1.2976874  0 Inf  0 Inf
      225     53     bd6t accuracy right  left 1.7298243  0 Inf  0 Inf
      227     54     bd6t accuracy right right 1.3249850  0 Inf  0 Inf
      229     55     bd6t accuracy right right 1.6483030  0 Inf  0 Inf
      231     56     bd6t accuracy right right 1.4835889  0 Inf  0 Inf
      233     57     bd6t accuracy right right 0.9261821  0 Inf  0 Inf
      235     58     bd6t accuracy right right 1.3110864  0 Inf  0 Inf
      237     59     bd6t accuracy right right 1.3092125  0 Inf  0 Inf
      239     60     bd6t accuracy right right 1.1064641  0 Inf  0 Inf

# DDM

    Code
      init_chains(DDM_s, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $theta_mu
                        [,1]
      v_Sleft     -0.9685927
      v_Sright     0.7061091
      a            1.4890213
      a_Eneutral  -1.8150926
      a_Eaccuracy  0.3304096
      t0          -1.1421557
      Z            0.1571934
      sv          -2.0654072
      SZ          -0.4405469
      
      $theta_var
      , , 1
      
                       v_Sleft      v_Sright            a    a_Eneutral  a_Eaccuracy
      v_Sleft      0.059329585 -0.0104497238 -0.010774234  0.0021450832 -0.018128601
      v_Sright    -0.010449724  0.0559162133 -0.002872401  0.0014127637  0.004728529
      a           -0.010774234 -0.0028724006  0.059773323  0.0087906529  0.022576320
      a_Eneutral   0.002145083  0.0014127637  0.008790653  0.0447136854  0.001458968
      a_Eaccuracy -0.018128601  0.0047285288  0.022576320  0.0014589675  0.038522063
      t0           0.018325623  0.0019110134 -0.037503336 -0.0100477306 -0.023423360
      Z           -0.012850266 -0.0030964417  0.020333128  0.0098693855  0.019288722
      sv          -0.009137082 -0.0181065189  0.007641907 -0.0029721272  0.008676374
      SZ          -0.019762600  0.0006819411  0.001329633  0.0006671515  0.007687211
                            t0            Z           sv            SZ
      v_Sleft      0.018325623 -0.012850266 -0.009137082 -0.0197625997
      v_Sright     0.001911013 -0.003096442 -0.018106519  0.0006819411
      a           -0.037503336  0.020333128  0.007641907  0.0013296327
      a_Eneutral  -0.010047731  0.009869385 -0.002972127  0.0006671515
      a_Eaccuracy -0.023423360  0.019288722  0.008676374  0.0076872108
      t0           0.078051779 -0.028702651 -0.004381484  0.0022878527
      Z           -0.028702651  0.042243108  0.005224348 -0.0090916441
      sv          -0.004381484  0.005224348  0.052137019  0.0067800486
      SZ           0.002287853 -0.009091644  0.006780049  0.0477904171
      
      
      $a_half
                       [,1]
      v_Sleft     0.4941231
      v_Sright    2.5931756
      a           1.4684537
      a_Eneutral  0.8437531
      a_Eaccuracy 0.6261293
      t0          0.8155694
      Z           2.2752573
      sv          0.7726374
      SZ          0.3149199
      
      $alpha
      , , 1
      
                         as1t       bd6t
      v_Sleft     -0.81572253 -1.0867281
      v_Sright     1.03870104  0.7321494
      a            1.22854033  1.3829047
      a_Eneutral  -2.14289844 -1.4698452
      a_Eaccuracy -0.09083814  0.3671716
      t0          -1.13678510 -1.6263327
      Z            0.10742647  0.2940408
      sv          -2.61001680 -2.2368029
      SZ          -0.56886020 -0.2263479
      
      
      $stage
      [1] "init"
      
      $subj_ll
                [,1]
      as1t -4820.904
      bd6t -3994.957
      
      $last_theta_var_inv
                 [,1]      [,2]       [,3]      [,4]       [,5]      [,6]
       [1,] 26.028515  6.224107 -1.3572503 -3.811330   3.070339 -3.657354
       [2,]  6.224107 22.834612  1.1572385 -1.918507  -6.426290 -1.390907
       [3,] -1.357250  1.157238 27.2973400 -2.546630  -9.560950  9.393468
       [4,] -3.811330 -1.918507 -2.5466296 25.241309   4.568563  1.562002
       [5,]  3.070339 -6.426290 -9.5609496  4.568563  46.584345  3.445213
       [6,] -3.657354 -1.390907  9.3934682  1.562002   3.445213 22.028224
       [7,]  8.055800  5.016142 -2.0467420 -8.084231 -16.170739  7.104309
       [8,]  3.587576  8.795323 -1.3430047  1.224921  -4.769682 -1.611160
       [9,] 11.470411  3.049327 -0.4123182 -4.351590  -8.494190 -1.804338
                    [,7]        [,8]       [,9]
       [1,]   8.05580010  3.58757606 11.4704105
       [2,]   5.01614235  8.79532339  3.0493269
       [3,]  -2.04674198 -1.34300467 -0.4123182
       [4,]  -8.08423140  1.22492096 -4.3515897
       [5,] -16.17073884 -4.76968243 -8.4941901
       [6,]   7.10430863 -1.61116037 -1.8043382
       [7,]  44.63330665 -0.03626134 14.1866918
       [8,]  -0.03626134 23.94575450 -1.1814287
       [9,]  14.18669184 -1.18142869 30.0158974
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_DDM, design_DDM, n_trials = 10)
    Output
          trials subjects        E     S     R        rt LT  UT LC  UC
      1        1     as1t    speed  left right 0.8475210  0 Inf  0 Inf
      2        2     as1t    speed  left right 1.3977571  0 Inf  0 Inf
      3        3     as1t    speed  left  left 0.5479478  0 Inf  0 Inf
      4        4     as1t    speed  left  left 0.5043682  0 Inf  0 Inf
      5        5     as1t    speed  left  left 0.6161943  0 Inf  0 Inf
      6        6     as1t    speed  left  left 0.3536291  0 Inf  0 Inf
      7        7     as1t    speed  left right 0.3705458  0 Inf  0 Inf
      8        8     as1t    speed  left right 0.4237434  0 Inf  0 Inf
      9        9     as1t    speed  left  left 0.4199297  0 Inf  0 Inf
      10      10     as1t    speed  left  left 0.3581555  0 Inf  0 Inf
      21      11     as1t  neutral  left right 0.6422508  0 Inf  0 Inf
      22      12     as1t  neutral  left  left 0.5621679  0 Inf  0 Inf
      23      13     as1t  neutral  left  left 0.3364405  0 Inf  0 Inf
      24      14     as1t  neutral  left  left 0.5094026  0 Inf  0 Inf
      25      15     as1t  neutral  left  left 0.8267394  0 Inf  0 Inf
      26      16     as1t  neutral  left right 0.6065282  0 Inf  0 Inf
      27      17     as1t  neutral  left  left 0.3984644  0 Inf  0 Inf
      28      18     as1t  neutral  left  left 0.4752852  0 Inf  0 Inf
      29      19     as1t  neutral  left right 2.7330970  0 Inf  0 Inf
      30      20     as1t  neutral  left  left 0.9813589  0 Inf  0 Inf
      41      21     as1t accuracy  left right 0.4583993  0 Inf  0 Inf
      42      22     as1t accuracy  left right 0.4116267  0 Inf  0 Inf
      43      23     as1t accuracy  left  left 0.5088002  0 Inf  0 Inf
      44      24     as1t accuracy  left  left 0.6617732  0 Inf  0 Inf
      45      25     as1t accuracy  left  left 0.9081772  0 Inf  0 Inf
      46      26     as1t accuracy  left right 0.5397453  0 Inf  0 Inf
      47      27     as1t accuracy  left  left 0.3820394  0 Inf  0 Inf
      48      28     as1t accuracy  left  left 2.2692272  0 Inf  0 Inf
      49      29     as1t accuracy  left  left 0.6142529  0 Inf  0 Inf
      50      30     as1t accuracy  left right 1.1070516  0 Inf  0 Inf
      61      31     as1t    speed right right 1.4224980  0 Inf  0 Inf
      62      32     as1t    speed right right 0.6981349  0 Inf  0 Inf
      63      33     as1t    speed right  left 1.9124266  0 Inf  0 Inf
      64      34     as1t    speed right right 0.6658440  0 Inf  0 Inf
      65      35     as1t    speed right  left 0.3721390  0 Inf  0 Inf
      66      36     as1t    speed right right 0.3419176  0 Inf  0 Inf
      67      37     as1t    speed right right 1.0700717  0 Inf  0 Inf
      68      38     as1t    speed right right 0.8066363  0 Inf  0 Inf
      69      39     as1t    speed right right 0.9798979  0 Inf  0 Inf
      70      40     as1t    speed right right 0.2891666  0 Inf  0 Inf
      81      41     as1t  neutral right right 1.9109804  0 Inf  0 Inf
      82      42     as1t  neutral right right 0.3832704  0 Inf  0 Inf
      83      43     as1t  neutral right right 0.3035224  0 Inf  0 Inf
      84      44     as1t  neutral right  left 0.8363085  0 Inf  0 Inf
      85      45     as1t  neutral right right 0.5065339  0 Inf  0 Inf
      86      46     as1t  neutral right right 0.6796286  0 Inf  0 Inf
      87      47     as1t  neutral right right 0.3470856  0 Inf  0 Inf
      88      48     as1t  neutral right  left 0.8067003  0 Inf  0 Inf
      89      49     as1t  neutral right right 0.4305781  0 Inf  0 Inf
      90      50     as1t  neutral right right 0.4647672  0 Inf  0 Inf
      101     51     as1t accuracy right right 0.4311426  0 Inf  0 Inf
      102     52     as1t accuracy right right 0.4152764  0 Inf  0 Inf
      103     53     as1t accuracy right right 0.3863490  0 Inf  0 Inf
      104     54     as1t accuracy right right 0.6896771  0 Inf  0 Inf
      105     55     as1t accuracy right right 0.5634301  0 Inf  0 Inf
      106     56     as1t accuracy right  left 0.4589499  0 Inf  0 Inf
      107     57     as1t accuracy right right 0.4440003  0 Inf  0 Inf
      108     58     as1t accuracy right right 0.3316251  0 Inf  0 Inf
      109     59     as1t accuracy right  left 2.2897616  0 Inf  0 Inf
      110     60     as1t accuracy right right 1.3066917  0 Inf  0 Inf
      11       1     bd6t    speed  left right 0.4831361  0 Inf  0 Inf
      12       2     bd6t    speed  left right 3.5592794  0 Inf  0 Inf
      13       3     bd6t    speed  left  left 0.3608753  0 Inf  0 Inf
      14       4     bd6t    speed  left  left 0.3871631  0 Inf  0 Inf
      15       5     bd6t    speed  left  left 0.3527741  0 Inf  0 Inf
      16       6     bd6t    speed  left right 0.9599148  0 Inf  0 Inf
      17       7     bd6t    speed  left  left 0.3550514  0 Inf  0 Inf
      18       8     bd6t    speed  left  left 0.4878590  0 Inf  0 Inf
      19       9     bd6t    speed  left  left 0.3188284  0 Inf  0 Inf
      20      10     bd6t    speed  left  left 0.6249439  0 Inf  0 Inf
      31      11     bd6t  neutral  left  left 0.9004394  0 Inf  0 Inf
      32      12     bd6t  neutral  left  left 0.5566473  0 Inf  0 Inf
      33      13     bd6t  neutral  left right 0.9992054  0 Inf  0 Inf
      34      14     bd6t  neutral  left right 0.7996416  0 Inf  0 Inf
      35      15     bd6t  neutral  left  left 1.8117929  0 Inf  0 Inf
      36      16     bd6t  neutral  left  left 0.5378516  0 Inf  0 Inf
      37      17     bd6t  neutral  left  left 0.5007347  0 Inf  0 Inf
      38      18     bd6t  neutral  left  left 1.1165995  0 Inf  0 Inf
      39      19     bd6t  neutral  left  left 0.5918927  0 Inf  0 Inf
      40      20     bd6t  neutral  left  left 2.6302307  0 Inf  0 Inf
      51      21     bd6t accuracy  left  left 0.5961940  0 Inf  0 Inf
      52      22     bd6t accuracy  left  left 1.0309690  0 Inf  0 Inf
      53      23     bd6t accuracy  left  left 0.7042805  0 Inf  0 Inf
      54      24     bd6t accuracy  left right 0.3834298  0 Inf  0 Inf
      55      25     bd6t accuracy  left  left 0.4855750  0 Inf  0 Inf
      56      26     bd6t accuracy  left  left 1.9538256  0 Inf  0 Inf
      57      27     bd6t accuracy  left  left 0.5784256  0 Inf  0 Inf
      58      28     bd6t accuracy  left right 0.4473647  0 Inf  0 Inf
      59      29     bd6t accuracy  left right 0.2976676  0 Inf  0 Inf
      60      30     bd6t accuracy  left right 0.7603309  0 Inf  0 Inf
      71      31     bd6t    speed right  left 0.3008438  0 Inf  0 Inf
      72      32     bd6t    speed right right 0.4273279  0 Inf  0 Inf
      73      33     bd6t    speed right right 0.5284606  0 Inf  0 Inf
      74      34     bd6t    speed right  left 1.7334288  0 Inf  0 Inf
      75      35     bd6t    speed right  left 0.6242915  0 Inf  0 Inf
      76      36     bd6t    speed right  left 0.8399532  0 Inf  0 Inf
      77      37     bd6t    speed right right 0.3349009  0 Inf  0 Inf
      78      38     bd6t    speed right right 1.4237554  0 Inf  0 Inf
      79      39     bd6t    speed right right 0.3445491  0 Inf  0 Inf
      80      40     bd6t    speed right right 0.4406860  0 Inf  0 Inf
      91      41     bd6t  neutral right right 2.5831702  0 Inf  0 Inf
      92      42     bd6t  neutral right  left 0.4957517  0 Inf  0 Inf
      93      43     bd6t  neutral right right 0.4600258  0 Inf  0 Inf
      94      44     bd6t  neutral right  left 0.7558276  0 Inf  0 Inf
      95      45     bd6t  neutral right right 2.1385071  0 Inf  0 Inf
      96      46     bd6t  neutral right  left 1.3245201  0 Inf  0 Inf
      97      47     bd6t  neutral right  left 1.5780413  0 Inf  0 Inf
      98      48     bd6t  neutral right right 0.5081798  0 Inf  0 Inf
      99      49     bd6t  neutral right right 0.3131300  0 Inf  0 Inf
      100     50     bd6t  neutral right right 0.4269199  0 Inf  0 Inf
      111     51     bd6t accuracy right right 0.9688415  0 Inf  0 Inf
      112     52     bd6t accuracy right right 0.5108115  0 Inf  0 Inf
      113     53     bd6t accuracy right  left 1.2113939  0 Inf  0 Inf
      114     54     bd6t accuracy right right 2.2105861  0 Inf  0 Inf
      115     55     bd6t accuracy right right 0.4278465  0 Inf  0 Inf
      116     56     bd6t accuracy right  left 0.4004046  0 Inf  0 Inf
      117     57     bd6t accuracy right right 0.4267699  0 Inf  0 Inf
      118     58     bd6t accuracy right  left 0.7489315  0 Inf  0 Inf
      119     59     bd6t accuracy right right 0.3209121  0 Inf  0 Inf
      120     60     bd6t accuracy right right 0.7879015  0 Inf  0 Inf

