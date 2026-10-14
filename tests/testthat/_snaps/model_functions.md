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
      
                         as1t       bd6t
      m           -0.03213819 -0.5094484
      m_lMd       -0.04820415 -0.0517495
      m_Eneutral  -0.05186968  0.1085556
      m_Eaccuracy  0.37267253 -0.3556700
      s            0.05738331 -0.4030297
      t0           0.39541015  0.6241218
      
      
      $stage
      [1] "init"
      
      $subj_ll
                [,1]
      as1t -18650.94
      bd6t -19548.95
      
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
      1        1     as1t    speed  left  left 0.3724880  0 Inf  0 Inf
      3        2     as1t    speed  left right 0.2772760  0 Inf  0 Inf
      5        3     as1t    speed  left right 0.2854363  0 Inf  0 Inf
      7        4     as1t    speed  left right 0.2701006  0 Inf  0 Inf
      9        5     as1t    speed  left  left 0.4030253  0 Inf  0 Inf
      11       6     as1t    speed  left right 0.6994876  0 Inf  0 Inf
      13       7     as1t    speed  left right 0.2652369  0 Inf  0 Inf
      15       8     as1t    speed  left right 0.3085436  0 Inf  0 Inf
      17       9     as1t    speed  left right 0.2526706  0 Inf  0 Inf
      19      10     as1t    speed  left right 0.4804292  0 Inf  0 Inf
      21      11     as1t  neutral  left right 0.4819677  0 Inf  0 Inf
      23      12     as1t  neutral  left right 0.9894224  0 Inf  0 Inf
      25      13     as1t  neutral  left right 0.2743063  0 Inf  0 Inf
      27      14     as1t  neutral  left  left 0.2967995  0 Inf  0 Inf
      29      15     as1t  neutral  left right 0.3924219  0 Inf  0 Inf
      31      16     as1t  neutral  left right 0.3083192  0 Inf  0 Inf
      33      17     as1t  neutral  left right 0.6628815  0 Inf  0 Inf
      35      18     as1t  neutral  left right 0.3476982  0 Inf  0 Inf
      37      19     as1t  neutral  left right 0.2567883  0 Inf  0 Inf
      39      20     as1t  neutral  left right 0.3528746  0 Inf  0 Inf
      41      21     as1t accuracy  left right 0.5363136  0 Inf  0 Inf
      43      22     as1t accuracy  left right 1.0610174  0 Inf  0 Inf
      45      23     as1t accuracy  left right 0.2788424  0 Inf  0 Inf
      47      24     as1t accuracy  left right 0.2620292  0 Inf  0 Inf
      49      25     as1t accuracy  left right 0.3697284  0 Inf  0 Inf
      51      26     as1t accuracy  left right 0.3048833  0 Inf  0 Inf
      53      27     as1t accuracy  left right 0.3474396  0 Inf  0 Inf
      55      28     as1t accuracy  left right 0.2521762  0 Inf  0 Inf
      57      29     as1t accuracy  left  left 0.2812319  0 Inf  0 Inf
      59      30     as1t accuracy  left right 0.2514580  0 Inf  0 Inf
      61      31     as1t    speed right  left 0.2769831  0 Inf  0 Inf
      63      32     as1t    speed right  left 0.5116264  0 Inf  0 Inf
      65      33     as1t    speed right  left 0.5632418  0 Inf  0 Inf
      67      34     as1t    speed right right 0.4621504  0 Inf  0 Inf
      69      35     as1t    speed right  left 0.2813245  0 Inf  0 Inf
      71      36     as1t    speed right right 0.3219558  0 Inf  0 Inf
      73      37     as1t    speed right  left 0.2934625  0 Inf  0 Inf
      75      38     as1t    speed right  left 0.6106815  0 Inf  0 Inf
      77      39     as1t    speed right  left 0.3216943  0 Inf  0 Inf
      79      40     as1t    speed right  left 0.2962127  0 Inf  0 Inf
      81      41     as1t  neutral right  left 0.3865330  0 Inf  0 Inf
      83      42     as1t  neutral right  left 0.4970369  0 Inf  0 Inf
      85      43     as1t  neutral right  left 1.2462713  0 Inf  0 Inf
      87      44     as1t  neutral right  left 0.3631830  0 Inf  0 Inf
      89      45     as1t  neutral right  left 0.4494020  0 Inf  0 Inf
      91      46     as1t  neutral right  left 1.1508538  0 Inf  0 Inf
      93      47     as1t  neutral right  left 0.2545263  0 Inf  0 Inf
      95      48     as1t  neutral right right 0.5079290  0 Inf  0 Inf
      97      49     as1t  neutral right  left 0.4789439  0 Inf  0 Inf
      99      50     as1t  neutral right  left 0.3578732  0 Inf  0 Inf
      101     51     as1t accuracy right  left 0.3984332  0 Inf  0 Inf
      103     52     as1t accuracy right  left 0.4819588  0 Inf  0 Inf
      105     53     as1t accuracy right right 0.3755497  0 Inf  0 Inf
      107     54     as1t accuracy right right 0.4277698  0 Inf  0 Inf
      109     55     as1t accuracy right  left 0.5617870  0 Inf  0 Inf
      111     56     as1t accuracy right  left 0.3340452  0 Inf  0 Inf
      113     57     as1t accuracy right  left 0.2577781  0 Inf  0 Inf
      115     58     as1t accuracy right  left 0.5275184  0 Inf  0 Inf
      117     59     as1t accuracy right  left 0.3634105  0 Inf  0 Inf
      119     60     as1t accuracy right  left 0.2664368  0 Inf  0 Inf
      121      1     bd6t    speed  left right 0.3000167  0 Inf  0 Inf
      123      2     bd6t    speed  left  left 0.5653744  0 Inf  0 Inf
      125      3     bd6t    speed  left right 0.2636713  0 Inf  0 Inf
      127      4     bd6t    speed  left right 0.2567060  0 Inf  0 Inf
      129      5     bd6t    speed  left right 0.4884221  0 Inf  0 Inf
      131      6     bd6t    speed  left  left 0.4527332  0 Inf  0 Inf
      133      7     bd6t    speed  left right 0.3498384  0 Inf  0 Inf
      135      8     bd6t    speed  left right 0.2780791  0 Inf  0 Inf
      137      9     bd6t    speed  left right 0.3531266  0 Inf  0 Inf
      139     10     bd6t    speed  left right 0.4660953  0 Inf  0 Inf
      141     11     bd6t  neutral  left  left 0.2650449  0 Inf  0 Inf
      143     12     bd6t  neutral  left right 0.4423806  0 Inf  0 Inf
      145     13     bd6t  neutral  left right 0.3021765  0 Inf  0 Inf
      147     14     bd6t  neutral  left right 0.8920605  0 Inf  0 Inf
      149     15     bd6t  neutral  left right 0.5583832  0 Inf  0 Inf
      151     16     bd6t  neutral  left right 0.2630774  0 Inf  0 Inf
      153     17     bd6t  neutral  left right 0.2778625  0 Inf  0 Inf
      155     18     bd6t  neutral  left right 0.2852806  0 Inf  0 Inf
      157     19     bd6t  neutral  left right 0.3022588  0 Inf  0 Inf
      159     20     bd6t  neutral  left right 0.2830901  0 Inf  0 Inf
      161     21     bd6t accuracy  left right 0.2646149  0 Inf  0 Inf
      163     22     bd6t accuracy  left right 0.4449964  0 Inf  0 Inf
      165     23     bd6t accuracy  left right 0.3433916  0 Inf  0 Inf
      167     24     bd6t accuracy  left right 0.6867135  0 Inf  0 Inf
      169     25     bd6t accuracy  left right 0.2621097  0 Inf  0 Inf
      171     26     bd6t accuracy  left right 0.2600779  0 Inf  0 Inf
      173     27     bd6t accuracy  left right 0.3090189  0 Inf  0 Inf
      175     28     bd6t accuracy  left right 0.2725345  0 Inf  0 Inf
      177     29     bd6t accuracy  left right 0.4671558  0 Inf  0 Inf
      179     30     bd6t accuracy  left  left 0.4014937  0 Inf  0 Inf
      181     31     bd6t    speed right  left 4.6847640  0 Inf  0 Inf
      183     32     bd6t    speed right  left 0.2794992  0 Inf  0 Inf
      185     33     bd6t    speed right  left 0.3356515  0 Inf  0 Inf
      187     34     bd6t    speed right  left 0.2541104  0 Inf  0 Inf
      189     35     bd6t    speed right  left 0.2935389  0 Inf  0 Inf
      191     36     bd6t    speed right  left 0.4257083  0 Inf  0 Inf
      193     37     bd6t    speed right  left 0.3438369  0 Inf  0 Inf
      195     38     bd6t    speed right right 0.3005960  0 Inf  0 Inf
      197     39     bd6t    speed right  left 0.7938186  0 Inf  0 Inf
      199     40     bd6t    speed right right 0.3362129  0 Inf  0 Inf
      201     41     bd6t  neutral right  left 0.2576989  0 Inf  0 Inf
      203     42     bd6t  neutral right  left 0.3046633  0 Inf  0 Inf
      205     43     bd6t  neutral right  left 0.3857014  0 Inf  0 Inf
      207     44     bd6t  neutral right  left 0.2985620  0 Inf  0 Inf
      209     45     bd6t  neutral right  left 0.2784170  0 Inf  0 Inf
      211     46     bd6t  neutral right  left 0.8339503  0 Inf  0 Inf
      213     47     bd6t  neutral right  left 0.2734772  0 Inf  0 Inf
      215     48     bd6t  neutral right  left 0.7635456  0 Inf  0 Inf
      217     49     bd6t  neutral right  left 0.3965151  0 Inf  0 Inf
      219     50     bd6t  neutral right  left 0.4111702  0 Inf  0 Inf
      221     51     bd6t accuracy right  left 0.3270322  0 Inf  0 Inf
      223     52     bd6t accuracy right  left 0.5160723  0 Inf  0 Inf
      225     53     bd6t accuracy right right 0.6165229  0 Inf  0 Inf
      227     54     bd6t accuracy right right 0.3134043  0 Inf  0 Inf
      229     55     bd6t accuracy right  left 0.3407393  0 Inf  0 Inf
      231     56     bd6t accuracy right  left 0.5679346  0 Inf  0 Inf
      233     57     bd6t accuracy right  left 0.3032600  0 Inf  0 Inf
      235     58     bd6t accuracy right  left 0.2806803  0 Inf  0 Inf
      237     59     bd6t accuracy right  left 0.2722017  0 Inf  0 Inf
      239     60     bd6t accuracy right  left 0.9823968  0 Inf  0 Inf

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
      as1t -16826.52
      bd6t -17318.42
      
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
      1        1     as1t    speed  left  left 1.2965552  0 Inf  0 Inf
      3        2     as1t    speed  left  left 1.3534176  0 Inf  0 Inf
      5        3     as1t    speed  left  left 1.2305393  0 Inf  0 Inf
      7        4     as1t    speed  left  left 1.3918411  0 Inf  0 Inf
      9        5     as1t    speed  left  left 1.0361262  0 Inf  0 Inf
      11       6     as1t    speed  left  left 1.4632735  0 Inf  0 Inf
      13       7     as1t    speed  left  left 1.4521725  0 Inf  0 Inf
      15       8     as1t    speed  left  left 1.1831674  0 Inf  0 Inf
      17       9     as1t    speed  left  left 1.2776537  0 Inf  0 Inf
      19      10     as1t    speed  left right 1.4092066  0 Inf  0 Inf
      21      11     as1t  neutral  left  left 1.3730172  0 Inf  0 Inf
      23      12     as1t  neutral  left  left 2.5619601  0 Inf  0 Inf
      25      13     as1t  neutral  left  left 1.7906202  0 Inf  0 Inf
      27      14     as1t  neutral  left right 1.5916791  0 Inf  0 Inf
      29      15     as1t  neutral  left  left 1.7709011  0 Inf  0 Inf
      31      16     as1t  neutral  left  left 2.5080629  0 Inf  0 Inf
      33      17     as1t  neutral  left  left 1.7133046  0 Inf  0 Inf
      35      18     as1t  neutral  left  left 2.3452663  0 Inf  0 Inf
      37      19     as1t  neutral  left  left 1.7797220  0 Inf  0 Inf
      39      20     as1t  neutral  left  left 2.2641737  0 Inf  0 Inf
      41      21     as1t accuracy  left  left 1.9795338  0 Inf  0 Inf
      43      22     as1t accuracy  left  left 1.7645306  0 Inf  0 Inf
      45      23     as1t accuracy  left right 1.2362568  0 Inf  0 Inf
      47      24     as1t accuracy  left  left 1.3724747  0 Inf  0 Inf
      49      25     as1t accuracy  left right 1.5190753  0 Inf  0 Inf
      51      26     as1t accuracy  left  left 1.5572071  0 Inf  0 Inf
      53      27     as1t accuracy  left  left 1.4293693  0 Inf  0 Inf
      55      28     as1t accuracy  left right 1.9542846  0 Inf  0 Inf
      57      29     as1t accuracy  left  left 1.2744618  0 Inf  0 Inf
      59      30     as1t accuracy  left  left 1.6914952  0 Inf  0 Inf
      61      31     as1t    speed right right 1.6646671  0 Inf  0 Inf
      63      32     as1t    speed right right 0.9586171  0 Inf  0 Inf
      65      33     as1t    speed right right 1.3391580  0 Inf  0 Inf
      67      34     as1t    speed right  left 1.3389733  0 Inf  0 Inf
      69      35     as1t    speed right  left 1.1693483  0 Inf  0 Inf
      71      36     as1t    speed right right 1.1306680  0 Inf  0 Inf
      73      37     as1t    speed right right 1.8718574  0 Inf  0 Inf
      75      38     as1t    speed right right 1.0364355  0 Inf  0 Inf
      77      39     as1t    speed right right 1.4269596  0 Inf  0 Inf
      79      40     as1t    speed right right 1.4804771  0 Inf  0 Inf
      81      41     as1t  neutral right right 2.6936757  0 Inf  0 Inf
      83      42     as1t  neutral right right 2.4068372  0 Inf  0 Inf
      85      43     as1t  neutral right right 1.6303361  0 Inf  0 Inf
      87      44     as1t  neutral right  left 1.7637833  0 Inf  0 Inf
      89      45     as1t  neutral right  left 2.2247968  0 Inf  0 Inf
      91      46     as1t  neutral right right 1.9316486  0 Inf  0 Inf
      93      47     as1t  neutral right right 1.6458289  0 Inf  0 Inf
      95      48     as1t  neutral right  left 1.7089033  0 Inf  0 Inf
      97      49     as1t  neutral right right 2.5458337  0 Inf  0 Inf
      99      50     as1t  neutral right right 2.1292790  0 Inf  0 Inf
      101     51     as1t accuracy right right 1.6374558  0 Inf  0 Inf
      103     52     as1t accuracy right right 1.5425385  0 Inf  0 Inf
      105     53     as1t accuracy right right 1.5613662  0 Inf  0 Inf
      107     54     as1t accuracy right right 1.8665443  0 Inf  0 Inf
      109     55     as1t accuracy right right 1.3939576  0 Inf  0 Inf
      111     56     as1t accuracy right  left 1.0932851  0 Inf  0 Inf
      113     57     as1t accuracy right  left 1.5695506  0 Inf  0 Inf
      115     58     as1t accuracy right right 2.0972658  0 Inf  0 Inf
      117     59     as1t accuracy right right 2.0133100  0 Inf  0 Inf
      119     60     as1t accuracy right right 2.0573826  0 Inf  0 Inf
      121      1     bd6t    speed  left  left 1.6504488  0 Inf  0 Inf
      123      2     bd6t    speed  left  left 1.3168622  0 Inf  0 Inf
      125      3     bd6t    speed  left  left 1.1923500  0 Inf  0 Inf
      127      4     bd6t    speed  left  left 1.1552406  0 Inf  0 Inf
      129      5     bd6t    speed  left  left 1.4287524  0 Inf  0 Inf
      131      6     bd6t    speed  left  left 1.5572416  0 Inf  0 Inf
      133      7     bd6t    speed  left  left 1.4281470  0 Inf  0 Inf
      135      8     bd6t    speed  left  left 1.2891465  0 Inf  0 Inf
      137      9     bd6t    speed  left  left 1.5462707  0 Inf  0 Inf
      139     10     bd6t    speed  left  left 1.0544393  0 Inf  0 Inf
      141     11     bd6t  neutral  left right 2.3723386  0 Inf  0 Inf
      143     12     bd6t  neutral  left  left 1.3877926  0 Inf  0 Inf
      145     13     bd6t  neutral  left right 1.6020387  0 Inf  0 Inf
      147     14     bd6t  neutral  left right 1.6827961  0 Inf  0 Inf
      149     15     bd6t  neutral  left right 1.9047115  0 Inf  0 Inf
      151     16     bd6t  neutral  left right 1.7904609  0 Inf  0 Inf
      153     17     bd6t  neutral  left right 1.6188581  0 Inf  0 Inf
      155     18     bd6t  neutral  left  left 2.0430401  0 Inf  0 Inf
      157     19     bd6t  neutral  left  left 2.1224393  0 Inf  0 Inf
      159     20     bd6t  neutral  left  left 1.6857447  0 Inf  0 Inf
      161     21     bd6t accuracy  left  left 1.3340224  0 Inf  0 Inf
      163     22     bd6t accuracy  left right 2.1139143  0 Inf  0 Inf
      165     23     bd6t accuracy  left right 1.4037548  0 Inf  0 Inf
      167     24     bd6t accuracy  left  left 1.7256191  0 Inf  0 Inf
      169     25     bd6t accuracy  left  left 1.9006757  0 Inf  0 Inf
      171     26     bd6t accuracy  left  left 2.1722071  0 Inf  0 Inf
      173     27     bd6t accuracy  left right 1.2799516  0 Inf  0 Inf
      175     28     bd6t accuracy  left  left 2.2737814  0 Inf  0 Inf
      177     29     bd6t accuracy  left  left 1.4312823  0 Inf  0 Inf
      179     30     bd6t accuracy  left  left 1.7328319  0 Inf  0 Inf
      181     31     bd6t    speed right right 1.3695219  0 Inf  0 Inf
      183     32     bd6t    speed right right 1.2572264  0 Inf  0 Inf
      185     33     bd6t    speed right right 1.1087944  0 Inf  0 Inf
      187     34     bd6t    speed right  left 1.2281858  0 Inf  0 Inf
      189     35     bd6t    speed right  left 1.5298429  0 Inf  0 Inf
      191     36     bd6t    speed right right 1.1931936  0 Inf  0 Inf
      193     37     bd6t    speed right right 1.3091507  0 Inf  0 Inf
      195     38     bd6t    speed right right 1.1512442  0 Inf  0 Inf
      197     39     bd6t    speed right right 1.3115294  0 Inf  0 Inf
      199     40     bd6t    speed right right 1.2419832  0 Inf  0 Inf
      201     41     bd6t  neutral right right 1.4249955  0 Inf  0 Inf
      203     42     bd6t  neutral right right 2.2834612  0 Inf  0 Inf
      205     43     bd6t  neutral right  left 3.9159141  0 Inf  0 Inf
      207     44     bd6t  neutral right right 1.7144087  0 Inf  0 Inf
      209     45     bd6t  neutral right right 2.1150569  0 Inf  0 Inf
      211     46     bd6t  neutral right right 3.0668672  0 Inf  0 Inf
      213     47     bd6t  neutral right right 2.0265500  0 Inf  0 Inf
      215     48     bd6t  neutral right right 2.7670114  0 Inf  0 Inf
      217     49     bd6t  neutral right  left 1.7448116  0 Inf  0 Inf
      219     50     bd6t  neutral right right 1.6700714  0 Inf  0 Inf
      221     51     bd6t accuracy right right 1.3478598  0 Inf  0 Inf
      223     52     bd6t accuracy right  left 1.3597237  0 Inf  0 Inf
      225     53     bd6t accuracy right right 2.5122839  0 Inf  0 Inf
      227     54     bd6t accuracy right  left 1.6615850  0 Inf  0 Inf
      229     55     bd6t accuracy right  left 1.2081940  0 Inf  0 Inf
      231     56     bd6t accuracy right right 1.6257200  0 Inf  0 Inf
      233     57     bd6t accuracy right  left 1.3604489  0 Inf  0 Inf
      235     58     bd6t accuracy right right 1.6195406  0 Inf  0 Inf
      237     59     bd6t accuracy right right 2.2556193  0 Inf  0 Inf
      239     60     bd6t accuracy right right 1.6629167  0 Inf  0 Inf

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
      as1t -10892.17
      bd6t -12483.34
      
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
      1        1     as1t    speed  left  left 1.4413789  0 Inf  0 Inf
      3        2     as1t    speed  left  left 1.1185888  0 Inf  0 Inf
      5        3     as1t    speed  left  left 1.0010308  0 Inf  0 Inf
      7        4     as1t    speed  left  left 0.8807496  0 Inf  0 Inf
      9        5     as1t    speed  left  left 0.9442313  0 Inf  0 Inf
      11       6     as1t    speed  left  left 0.9009904  0 Inf  0 Inf
      13       7     as1t    speed  left  left 1.2159235  0 Inf  0 Inf
      15       8     as1t    speed  left  left 1.2123634  0 Inf  0 Inf
      17       9     as1t    speed  left  left 1.0151265  0 Inf  0 Inf
      19      10     as1t    speed  left  left 1.0745897  0 Inf  0 Inf
      21      11     as1t  neutral  left  left 1.4949170  0 Inf  0 Inf
      23      12     as1t  neutral  left  left 1.6556311  0 Inf  0 Inf
      25      13     as1t  neutral  left  left 1.8479449  0 Inf  0 Inf
      27      14     as1t  neutral  left  left 1.4517764  0 Inf  0 Inf
      29      15     as1t  neutral  left  left 1.7702969  0 Inf  0 Inf
      31      16     as1t  neutral  left  left 1.4903889  0 Inf  0 Inf
      33      17     as1t  neutral  left  left 1.4068067  0 Inf  0 Inf
      35      18     as1t  neutral  left  left 1.5587795  0 Inf  0 Inf
      37      19     as1t  neutral  left  left 1.4749124  0 Inf  0 Inf
      39      20     as1t  neutral  left  left 1.9371228  0 Inf  0 Inf
      41      21     as1t accuracy  left  left 1.1652946  0 Inf  0 Inf
      43      22     as1t accuracy  left  left 1.5035249  0 Inf  0 Inf
      45      23     as1t accuracy  left  left 1.5324305  0 Inf  0 Inf
      47      24     as1t accuracy  left  left 1.3131464  0 Inf  0 Inf
      49      25     as1t accuracy  left  left 1.1127378  0 Inf  0 Inf
      51      26     as1t accuracy  left  left 1.3548187  0 Inf  0 Inf
      53      27     as1t accuracy  left  left 1.3789867  0 Inf  0 Inf
      55      28     as1t accuracy  left right 1.3678645  0 Inf  0 Inf
      57      29     as1t accuracy  left  left 1.5254417  0 Inf  0 Inf
      59      30     as1t accuracy  left  left 1.2885054  0 Inf  0 Inf
      61      31     as1t    speed right right 0.8782070  0 Inf  0 Inf
      63      32     as1t    speed right right 0.9883937  0 Inf  0 Inf
      65      33     as1t    speed right right 1.0434413  0 Inf  0 Inf
      67      34     as1t    speed right right 0.9084312  0 Inf  0 Inf
      69      35     as1t    speed right right 1.3344970  0 Inf  0 Inf
      71      36     as1t    speed right right 1.1497505  0 Inf  0 Inf
      73      37     as1t    speed right right 0.8250868  0 Inf  0 Inf
      75      38     as1t    speed right right 0.8586063  0 Inf  0 Inf
      77      39     as1t    speed right right 1.1809869  0 Inf  0 Inf
      79      40     as1t    speed right right 1.1254231  0 Inf  0 Inf
      81      41     as1t  neutral right  left 1.4007508  0 Inf  0 Inf
      83      42     as1t  neutral right right 1.5225026  0 Inf  0 Inf
      85      43     as1t  neutral right right 1.5162371  0 Inf  0 Inf
      87      44     as1t  neutral right right 1.4964466  0 Inf  0 Inf
      89      45     as1t  neutral right right 1.6131718  0 Inf  0 Inf
      91      46     as1t  neutral right right 1.7860756  0 Inf  0 Inf
      93      47     as1t  neutral right right 1.4977565  0 Inf  0 Inf
      95      48     as1t  neutral right right 1.5320221  0 Inf  0 Inf
      97      49     as1t  neutral right right 1.7171076  0 Inf  0 Inf
      99      50     as1t  neutral right right 1.3419129  0 Inf  0 Inf
      101     51     as1t accuracy right right 1.3061966  0 Inf  0 Inf
      103     52     as1t accuracy right right 1.4418711  0 Inf  0 Inf
      105     53     as1t accuracy right right 0.9826146  0 Inf  0 Inf
      107     54     as1t accuracy right right 1.6477141  0 Inf  0 Inf
      109     55     as1t accuracy right right 1.3300668  0 Inf  0 Inf
      111     56     as1t accuracy right right 1.1688516  0 Inf  0 Inf
      113     57     as1t accuracy right right 1.6718517  0 Inf  0 Inf
      115     58     as1t accuracy right right 1.2162995  0 Inf  0 Inf
      117     59     as1t accuracy right right 1.3689552  0 Inf  0 Inf
      119     60     as1t accuracy right right 1.6199666  0 Inf  0 Inf
      121      1     bd6t    speed  left  left 0.8899750  0 Inf  0 Inf
      123      2     bd6t    speed  left right 1.1179794  0 Inf  0 Inf
      125      3     bd6t    speed  left  left 1.3737672  0 Inf  0 Inf
      127      4     bd6t    speed  left  left 1.2565384  0 Inf  0 Inf
      129      5     bd6t    speed  left right 1.1082506  0 Inf  0 Inf
      131      6     bd6t    speed  left  left 1.5606688  0 Inf  0 Inf
      133      7     bd6t    speed  left  left 1.3856710  0 Inf  0 Inf
      135      8     bd6t    speed  left  left 1.0343146  0 Inf  0 Inf
      137      9     bd6t    speed  left  left 1.1909469  0 Inf  0 Inf
      139     10     bd6t    speed  left  left 0.7756392  0 Inf  0 Inf
      141     11     bd6t  neutral  left  left 1.7742525  0 Inf  0 Inf
      143     12     bd6t  neutral  left  left 1.3260621  0 Inf  0 Inf
      145     13     bd6t  neutral  left  left 1.6109828  0 Inf  0 Inf
      147     14     bd6t  neutral  left right 1.5841187  0 Inf  0 Inf
      149     15     bd6t  neutral  left  left 1.6267794  0 Inf  0 Inf
      151     16     bd6t  neutral  left  left 1.2663493  0 Inf  0 Inf
      153     17     bd6t  neutral  left  left 1.5636835  0 Inf  0 Inf
      155     18     bd6t  neutral  left  left 1.7226510  0 Inf  0 Inf
      157     19     bd6t  neutral  left  left 1.4329461  0 Inf  0 Inf
      159     20     bd6t  neutral  left  left 1.2711419  0 Inf  0 Inf
      161     21     bd6t accuracy  left  left 1.7233471  0 Inf  0 Inf
      163     22     bd6t accuracy  left  left 1.3031005  0 Inf  0 Inf
      165     23     bd6t accuracy  left  left 1.1494813  0 Inf  0 Inf
      167     24     bd6t accuracy  left  left 1.2724348  0 Inf  0 Inf
      169     25     bd6t accuracy  left  left 1.2048800  0 Inf  0 Inf
      171     26     bd6t accuracy  left  left 1.1767682  0 Inf  0 Inf
      173     27     bd6t accuracy  left  left 1.4486781  0 Inf  0 Inf
      175     28     bd6t accuracy  left  left 1.1689022  0 Inf  0 Inf
      177     29     bd6t accuracy  left right 0.9711361  0 Inf  0 Inf
      179     30     bd6t accuracy  left  left 1.2103535  0 Inf  0 Inf
      181     31     bd6t    speed right right 1.1917876  0 Inf  0 Inf
      183     32     bd6t    speed right right 0.8511150  0 Inf  0 Inf
      185     33     bd6t    speed right right 1.1005202  0 Inf  0 Inf
      187     34     bd6t    speed right right 1.0224577  0 Inf  0 Inf
      189     35     bd6t    speed right right 1.5332909  0 Inf  0 Inf
      191     36     bd6t    speed right right 1.1183719  0 Inf  0 Inf
      193     37     bd6t    speed right right 0.8011024  0 Inf  0 Inf
      195     38     bd6t    speed right  left 0.8939908  0 Inf  0 Inf
      197     39     bd6t    speed right  left 0.7494244  0 Inf  0 Inf
      199     40     bd6t    speed right right 0.9669068  0 Inf  0 Inf
      201     41     bd6t  neutral right right 1.5592809  0 Inf  0 Inf
      203     42     bd6t  neutral right right 1.1654861  0 Inf  0 Inf
      205     43     bd6t  neutral right right 1.3721378  0 Inf  0 Inf
      207     44     bd6t  neutral right right 1.2909112  0 Inf  0 Inf
      209     45     bd6t  neutral right right 1.8634184  0 Inf  0 Inf
      211     46     bd6t  neutral right right 1.5615342  0 Inf  0 Inf
      213     47     bd6t  neutral right right 1.2991772  0 Inf  0 Inf
      215     48     bd6t  neutral right  left 1.4827300  0 Inf  0 Inf
      217     49     bd6t  neutral right right 1.7409714  0 Inf  0 Inf
      219     50     bd6t  neutral right right 1.7280616  0 Inf  0 Inf
      221     51     bd6t accuracy right right 1.3408036  0 Inf  0 Inf
      223     52     bd6t accuracy right  left 1.1122038  0 Inf  0 Inf
      225     53     bd6t accuracy right right 1.1859820  0 Inf  0 Inf
      227     54     bd6t accuracy right right 0.9185411  0 Inf  0 Inf
      229     55     bd6t accuracy right right 1.3013777  0 Inf  0 Inf
      231     56     bd6t accuracy right right 1.2338889  0 Inf  0 Inf
      233     57     bd6t accuracy right right 1.4011179  0 Inf  0 Inf
      235     58     bd6t accuracy right right 1.0906175  0 Inf  0 Inf
      237     59     bd6t accuracy right right 1.2121172  0 Inf  0 Inf
      239     60     bd6t accuracy right right 1.4853299  0 Inf  0 Inf

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
      
                          as1t          bd6t
      v_Sleft      0.001153243  0.1958094525
      v_Sright    -0.739398191 -0.1689060112
      a           -0.287978912  0.2995111904
      a_Eneutral   0.207325046  0.5271980707
      a_Eaccuracy -0.719387674 -0.0751964646
      t0           0.717560487  0.4095226870
      Z            0.338818498  0.5809248993
      sv           0.206385123  0.0001835331
      SZ          -1.043411727 -1.1563846152
      
      
      $stage
      [1] "init"
      
      $subj_ll
                [,1]
      as1t -18650.94
      bd6t -19548.95
      
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
      1        1     as1t    speed  left right 0.8805799  0 Inf  0 Inf
      2        2     as1t    speed  left right 1.1747003  0 Inf  0 Inf
      3        3     as1t    speed  left right 0.5597313  0 Inf  0 Inf
      4        4     as1t    speed  left  left 0.3121141  0 Inf  0 Inf
      5        5     as1t    speed  left  left 0.5356980  0 Inf  0 Inf
      6        6     as1t    speed  left  left 0.3578257  0 Inf  0 Inf
      7        7     as1t    speed  left  left 0.3331043  0 Inf  0 Inf
      8        8     as1t    speed  left  left 0.4927992  0 Inf  0 Inf
      9        9     as1t    speed  left right 0.4184083  0 Inf  0 Inf
      10      10     as1t    speed  left  left 0.8194661  0 Inf  0 Inf
      21      11     as1t  neutral  left  left 2.2340234  0 Inf  0 Inf
      22      12     as1t  neutral  left  left 0.4405468  0 Inf  0 Inf
      23      13     as1t  neutral  left  left 0.6059885  0 Inf  0 Inf
      24      14     as1t  neutral  left  left 0.5771193  0 Inf  0 Inf
      25      15     as1t  neutral  left  left 1.0424137  0 Inf  0 Inf
      26      16     as1t  neutral  left  left 0.4153175  0 Inf  0 Inf
      27      17     as1t  neutral  left  left 0.4827289  0 Inf  0 Inf
      28      18     as1t  neutral  left  left 0.6252632  0 Inf  0 Inf
      29      19     as1t  neutral  left right 0.9538454  0 Inf  0 Inf
      30      20     as1t  neutral  left  left 1.1921494  0 Inf  0 Inf
      41      21     as1t accuracy  left  left 1.5267784  0 Inf  0 Inf
      42      22     as1t accuracy  left  left 0.5974937  0 Inf  0 Inf
      43      23     as1t accuracy  left  left 1.0257809  0 Inf  0 Inf
      44      24     as1t accuracy  left right 1.5296556  0 Inf  0 Inf
      45      25     as1t accuracy  left  left 0.6330345  0 Inf  0 Inf
      46      26     as1t accuracy  left  left 0.5747440  0 Inf  0 Inf
      47      27     as1t accuracy  left  left 0.7390723  0 Inf  0 Inf
      48      28     as1t accuracy  left  left 1.0881356  0 Inf  0 Inf
      49      29     as1t accuracy  left  left 1.2555026  0 Inf  0 Inf
      50      30     as1t accuracy  left  left 0.6523659  0 Inf  0 Inf
      61      31     as1t    speed right right 0.3493858  0 Inf  0 Inf
      62      32     as1t    speed right right 0.3391079  0 Inf  0 Inf
      63      33     as1t    speed right  left 0.4213549  0 Inf  0 Inf
      64      34     as1t    speed right right 0.4365562  0 Inf  0 Inf
      65      35     as1t    speed right right 0.3205181  0 Inf  0 Inf
      66      36     as1t    speed right right 1.4876486  0 Inf  0 Inf
      67      37     as1t    speed right right 0.4475905  0 Inf  0 Inf
      68      38     as1t    speed right right 0.3513510  0 Inf  0 Inf
      69      39     as1t    speed right right 0.5618381  0 Inf  0 Inf
      70      40     as1t    speed right right 0.4596839  0 Inf  0 Inf
      81      41     as1t  neutral right right 0.3439475  0 Inf  0 Inf
      82      42     as1t  neutral right  left 0.9939673  0 Inf  0 Inf
      83      43     as1t  neutral right right 0.6684577  0 Inf  0 Inf
      84      44     as1t  neutral right right 0.4785465  0 Inf  0 Inf
      85      45     as1t  neutral right right 0.4632035  0 Inf  0 Inf
      86      46     as1t  neutral right right 0.3738223  0 Inf  0 Inf
      87      47     as1t  neutral right  left 0.6073673  0 Inf  0 Inf
      88      48     as1t  neutral right right 0.4406135  0 Inf  0 Inf
      89      49     as1t  neutral right  left 0.8876660  0 Inf  0 Inf
      90      50     as1t  neutral right right 0.5759519  0 Inf  0 Inf
      101     51     as1t accuracy right right 0.5924398  0 Inf  0 Inf
      102     52     as1t accuracy right right 1.6281539  0 Inf  0 Inf
      103     53     as1t accuracy right right 1.4881673  0 Inf  0 Inf
      104     54     as1t accuracy right right 0.4137741  0 Inf  0 Inf
      105     55     as1t accuracy right right 0.8407528  0 Inf  0 Inf
      106     56     as1t accuracy right right 2.5914160  0 Inf  0 Inf
      107     57     as1t accuracy right  left 3.4075450  0 Inf  0 Inf
      108     58     as1t accuracy right right 0.5086852  0 Inf  0 Inf
      109     59     as1t accuracy right right 0.6296841  0 Inf  0 Inf
      110     60     as1t accuracy right right 0.9865453  0 Inf  0 Inf
      11       1     bd6t    speed  left  left 0.5616933  0 Inf  0 Inf
      12       2     bd6t    speed  left  left 0.8821916  0 Inf  0 Inf
      13       3     bd6t    speed  left  left 0.3356000  0 Inf  0 Inf
      14       4     bd6t    speed  left  left 1.6979495  0 Inf  0 Inf
      15       5     bd6t    speed  left right 0.4947166  0 Inf  0 Inf
      16       6     bd6t    speed  left  left 0.6144373  0 Inf  0 Inf
      17       7     bd6t    speed  left right 0.4366723  0 Inf  0 Inf
      18       8     bd6t    speed  left  left 0.4605378  0 Inf  0 Inf
      19       9     bd6t    speed  left  left 0.4874733  0 Inf  0 Inf
      20      10     bd6t    speed  left  left 0.5647492  0 Inf  0 Inf
      31      11     bd6t  neutral  left  left 0.7253865  0 Inf  0 Inf
      32      12     bd6t  neutral  left  left 0.5831586  0 Inf  0 Inf
      33      13     bd6t  neutral  left  left 0.7239770  0 Inf  0 Inf
      34      14     bd6t  neutral  left right 0.4588492  0 Inf  0 Inf
      35      15     bd6t  neutral  left right 0.7045694  0 Inf  0 Inf
      36      16     bd6t  neutral  left  left 0.8104721  0 Inf  0 Inf
      37      17     bd6t  neutral  left  left 0.8233182  0 Inf  0 Inf
      38      18     bd6t  neutral  left  left 3.0701048  0 Inf  0 Inf
      39      19     bd6t  neutral  left  left 1.3077321  0 Inf  0 Inf
      40      20     bd6t  neutral  left  left 0.4470535  0 Inf  0 Inf
      51      21     bd6t accuracy  left right 1.0250510  0 Inf  0 Inf
      52      22     bd6t accuracy  left  left 0.5399017  0 Inf  0 Inf
      53      23     bd6t accuracy  left  left 0.4700714  0 Inf  0 Inf
      54      24     bd6t accuracy  left  left 0.8607142  0 Inf  0 Inf
      55      25     bd6t accuracy  left  left 0.5554308  0 Inf  0 Inf
      56      26     bd6t accuracy  left  left 0.3598135  0 Inf  0 Inf
      57      27     bd6t accuracy  left  left 0.4506607  0 Inf  0 Inf
      58      28     bd6t accuracy  left  left 0.6494674  0 Inf  0 Inf
      59      29     bd6t accuracy  left  left 1.3346860  0 Inf  0 Inf
      60      30     bd6t accuracy  left  left 1.0101512  0 Inf  0 Inf
      71      31     bd6t    speed right right 0.5608489  0 Inf  0 Inf
      72      32     bd6t    speed right right 0.4445678  0 Inf  0 Inf
      73      33     bd6t    speed right right 0.7460522  0 Inf  0 Inf
      74      34     bd6t    speed right right 1.0679428  0 Inf  0 Inf
      75      35     bd6t    speed right right 0.4936465  0 Inf  0 Inf
      76      36     bd6t    speed right right 3.2095989  0 Inf  0 Inf
      77      37     bd6t    speed right right 0.8861243  0 Inf  0 Inf
      78      38     bd6t    speed right right 0.6711227  0 Inf  0 Inf
      79      39     bd6t    speed right right 0.6157516  0 Inf  0 Inf
      80      40     bd6t    speed right  left 0.9536515  0 Inf  0 Inf
      91      41     bd6t  neutral right right 1.7096533  0 Inf  0 Inf
      92      42     bd6t  neutral right right 0.4720140  0 Inf  0 Inf
      93      43     bd6t  neutral right right 0.6970336  0 Inf  0 Inf
      94      44     bd6t  neutral right right 0.4573710  0 Inf  0 Inf
      95      45     bd6t  neutral right right 0.7757981  0 Inf  0 Inf
      96      46     bd6t  neutral right right 0.5576210  0 Inf  0 Inf
      97      47     bd6t  neutral right right 0.8889030  0 Inf  0 Inf
      98      48     bd6t  neutral right right 0.5441407  0 Inf  0 Inf
      99      49     bd6t  neutral right right 0.4358213  0 Inf  0 Inf
      100     50     bd6t  neutral right  left 1.8191597  0 Inf  0 Inf
      111     51     bd6t accuracy right  left 0.8719539  0 Inf  0 Inf
      112     52     bd6t accuracy right  left 0.7368114  0 Inf  0 Inf
      113     53     bd6t accuracy right  left 0.8720341  0 Inf  0 Inf
      114     54     bd6t accuracy right right 1.2560636  0 Inf  0 Inf
      115     55     bd6t accuracy right right 0.5217418  0 Inf  0 Inf
      116     56     bd6t accuracy right  left 1.6370467  0 Inf  0 Inf
      117     57     bd6t accuracy right  left 1.5575868  0 Inf  0 Inf
      118     58     bd6t accuracy right  left 0.6817668  0 Inf  0 Inf
      119     59     bd6t accuracy right  left 0.8550562  0 Inf  0 Inf
      120     60     bd6t accuracy right right 0.5270226  0 Inf  0 Inf

