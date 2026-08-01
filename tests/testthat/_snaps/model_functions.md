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
      
                        as1t       bd6t
      m           -0.6599300 -0.5836701
      m_lMd        0.3772913  0.1802469
      m_Eneutral   1.6147295  1.7008123
      m_Eaccuracy -2.1286744 -1.3967109
      s            0.5289465  0.4253892
      t0          -1.3065202 -1.3823828
      
      
      $stage
      [1] "init"
      
      $subj_ll
                 [,1]
      as1t -1119.9841
      bd6t  -849.3553
      
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
      1        1     as1t    speed  left right 0.6621568  0 Inf  0 Inf
      3        2     as1t    speed  left right 0.2057630  0 Inf  0 Inf
      5        3     as1t    speed  left right 0.2148357  0 Inf  0 Inf
      7        4     as1t    speed  left  left 0.3104419  0 Inf  0 Inf
      9        5     as1t    speed  left  left 0.3391626  0 Inf  0 Inf
      11       6     as1t    speed  left right 0.3029197  0 Inf  0 Inf
      13       7     as1t    speed  left right 1.3993717  0 Inf  0 Inf
      15       8     as1t    speed  left right 0.2306459  0 Inf  0 Inf
      17       9     as1t    speed  left  left 1.1688779  0 Inf  0 Inf
      19      10     as1t    speed  left right 0.4420536  0 Inf  0 Inf
      21      11     as1t  neutral  left right 0.6309831  0 Inf  0 Inf
      23      12     as1t  neutral  left  left 0.2536847  0 Inf  0 Inf
      25      13     as1t  neutral  left right 0.2319419  0 Inf  0 Inf
      27      14     as1t  neutral  left right 0.2159625  0 Inf  0 Inf
      29      15     as1t  neutral  left right 0.3556382  0 Inf  0 Inf
      31      16     as1t  neutral  left  left 0.2119769  0 Inf  0 Inf
      33      17     as1t  neutral  left right 0.2055206  0 Inf  0 Inf
      35      18     as1t  neutral  left right 0.2933270  0 Inf  0 Inf
      37      19     as1t  neutral  left right 0.9621727  0 Inf  0 Inf
      39      20     as1t  neutral  left right 0.2841160  0 Inf  0 Inf
      41      21     as1t accuracy  left right 0.2082286  0 Inf  0 Inf
      43      22     as1t accuracy  left right 1.0291990  0 Inf  0 Inf
      45      23     as1t accuracy  left right 0.2387766  0 Inf  0 Inf
      47      24     as1t accuracy  left right 0.3454905  0 Inf  0 Inf
      49      25     as1t accuracy  left  left 0.3402073  0 Inf  0 Inf
      51      26     as1t accuracy  left right 0.2931404  0 Inf  0 Inf
      53      27     as1t accuracy  left right 0.5631042  0 Inf  0 Inf
      55      28     as1t accuracy  left right 0.2650557  0 Inf  0 Inf
      57      29     as1t accuracy  left right 0.2311833  0 Inf  0 Inf
      59      30     as1t accuracy  left right 0.3249946  0 Inf  0 Inf
      61      31     as1t    speed right  left 0.3348222  0 Inf  0 Inf
      63      32     as1t    speed right  left 0.3124411  0 Inf  0 Inf
      65      33     as1t    speed right right 0.3009970  0 Inf  0 Inf
      67      34     as1t    speed right right 0.3324145  0 Inf  0 Inf
      69      35     as1t    speed right  left 0.2808354  0 Inf  0 Inf
      71      36     as1t    speed right  left 0.5841607  0 Inf  0 Inf
      73      37     as1t    speed right  left 0.2735157  0 Inf  0 Inf
      75      38     as1t    speed right  left 0.2484093  0 Inf  0 Inf
      77      39     as1t    speed right right 0.2299376  0 Inf  0 Inf
      79      40     as1t    speed right  left 0.2606403  0 Inf  0 Inf
      81      41     as1t  neutral right right 0.2842899  0 Inf  0 Inf
      83      42     as1t  neutral right  left 0.2101295  0 Inf  0 Inf
      85      43     as1t  neutral right  left 1.4312545  0 Inf  0 Inf
      87      44     as1t  neutral right  left 0.2745290  0 Inf  0 Inf
      89      45     as1t  neutral right right 0.4214576  0 Inf  0 Inf
      91      46     as1t  neutral right  left 0.2168920  0 Inf  0 Inf
      93      47     as1t  neutral right right 0.2591430  0 Inf  0 Inf
      95      48     as1t  neutral right  left 1.0453442  0 Inf  0 Inf
      97      49     as1t  neutral right  left 0.3268598  0 Inf  0 Inf
      99      50     as1t  neutral right right 0.2522106  0 Inf  0 Inf
      101     51     as1t accuracy right right 2.6889300  0 Inf  0 Inf
      103     52     as1t accuracy right  left 0.2239633  0 Inf  0 Inf
      105     53     as1t accuracy right right 0.3311753  0 Inf  0 Inf
      107     54     as1t accuracy right right 0.3193368  0 Inf  0 Inf
      109     55     as1t accuracy right  left 0.4696293  0 Inf  0 Inf
      111     56     as1t accuracy right  left 0.3360226  0 Inf  0 Inf
      113     57     as1t accuracy right  left 0.4199177  0 Inf  0 Inf
      115     58     as1t accuracy right right 0.7826471  0 Inf  0 Inf
      117     59     as1t accuracy right  left 0.7467092  0 Inf  0 Inf
      119     60     as1t accuracy right right 0.2285965  0 Inf  0 Inf
      121      1     bd6t    speed  left right 0.3521819  0 Inf  0 Inf
      123      2     bd6t    speed  left  left 2.5742499  0 Inf  0 Inf
      125      3     bd6t    speed  left  left 0.2447084  0 Inf  0 Inf
      127      4     bd6t    speed  left right 0.3353470  0 Inf  0 Inf
      129      5     bd6t    speed  left right 0.3184436  0 Inf  0 Inf
      131      6     bd6t    speed  left right 0.4884841  0 Inf  0 Inf
      133      7     bd6t    speed  left right 0.4399945  0 Inf  0 Inf
      135      8     bd6t    speed  left right 0.2215287  0 Inf  0 Inf
      137      9     bd6t    speed  left right 0.2100897  0 Inf  0 Inf
      139     10     bd6t    speed  left right 0.2133570  0 Inf  0 Inf
      141     11     bd6t  neutral  left right 0.2446073  0 Inf  0 Inf
      143     12     bd6t  neutral  left  left 0.3726230  0 Inf  0 Inf
      145     13     bd6t  neutral  left right 0.4092510  0 Inf  0 Inf
      147     14     bd6t  neutral  left  left 0.5294351  0 Inf  0 Inf
      149     15     bd6t  neutral  left right 0.2143891  0 Inf  0 Inf
      151     16     bd6t  neutral  left right 0.9475786  0 Inf  0 Inf
      153     17     bd6t  neutral  left right 0.2133931  0 Inf  0 Inf
      155     18     bd6t  neutral  left right 1.2623103  0 Inf  0 Inf
      157     19     bd6t  neutral  left right 0.4351044  0 Inf  0 Inf
      159     20     bd6t  neutral  left right 0.4247098  0 Inf  0 Inf
      161     21     bd6t accuracy  left  left 0.3207819  0 Inf  0 Inf
      163     22     bd6t accuracy  left right 0.2389074  0 Inf  0 Inf
      165     23     bd6t accuracy  left right 0.2820200  0 Inf  0 Inf
      167     24     bd6t accuracy  left right 0.2198767  0 Inf  0 Inf
      169     25     bd6t accuracy  left right 0.4524300  0 Inf  0 Inf
      171     26     bd6t accuracy  left right 0.2300322  0 Inf  0 Inf
      173     27     bd6t accuracy  left right 0.2416458  0 Inf  0 Inf
      175     28     bd6t accuracy  left right 0.2191954  0 Inf  0 Inf
      177     29     bd6t accuracy  left right 0.6286309  0 Inf  0 Inf
      179     30     bd6t accuracy  left right 0.2628880  0 Inf  0 Inf
      181     31     bd6t    speed right  left 0.6651145  0 Inf  0 Inf
      183     32     bd6t    speed right  left 0.2299034  0 Inf  0 Inf
      185     33     bd6t    speed right right 2.4026591  0 Inf  0 Inf
      187     34     bd6t    speed right right 0.3498691  0 Inf  0 Inf
      189     35     bd6t    speed right right 0.3024430  0 Inf  0 Inf
      191     36     bd6t    speed right  left 0.3085880  0 Inf  0 Inf
      193     37     bd6t    speed right  left 0.2907286  0 Inf  0 Inf
      195     38     bd6t    speed right  left 0.2056963  0 Inf  0 Inf
      197     39     bd6t    speed right  left 0.2142630  0 Inf  0 Inf
      199     40     bd6t    speed right  left 0.4062331  0 Inf  0 Inf
      201     41     bd6t  neutral right  left 0.2521459  0 Inf  0 Inf
      203     42     bd6t  neutral right right 0.5737047  0 Inf  0 Inf
      205     43     bd6t  neutral right  left 0.2679879  0 Inf  0 Inf
      207     44     bd6t  neutral right  left 0.6547885  0 Inf  0 Inf
      209     45     bd6t  neutral right  left 0.7484957  0 Inf  0 Inf
      211     46     bd6t  neutral right  left 0.7427598  0 Inf  0 Inf
      213     47     bd6t  neutral right  left 0.2639333  0 Inf  0 Inf
      215     48     bd6t  neutral right  left 0.2970921  0 Inf  0 Inf
      217     49     bd6t  neutral right  left 0.4422466  0 Inf  0 Inf
      219     50     bd6t  neutral right  left 0.9803098  0 Inf  0 Inf
      221     51     bd6t accuracy right  left 0.7469548  0 Inf  0 Inf
      223     52     bd6t accuracy right right 0.2676699  0 Inf  0 Inf
      225     53     bd6t accuracy right right 1.3892407  0 Inf  0 Inf
      227     54     bd6t accuracy right  left 0.2842021  0 Inf  0 Inf
      229     55     bd6t accuracy right  left 0.3571977  0 Inf  0 Inf
      231     56     bd6t accuracy right  left 0.8165039  0 Inf  0 Inf
      233     57     bd6t accuracy right  left 0.2736779  0 Inf  0 Inf
      235     58     bd6t accuracy right  left 0.4250408  0 Inf  0 Inf
      237     59     bd6t accuracy right  left 0.3621805  0 Inf  0 Inf
      239     60     bd6t accuracy right  left 0.3439226  0 Inf  0 Inf

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
      v           -0.8659366 -0.8274847
      v_lMd        0.8828966  0.5761732
      sv_lMTRUE    1.4190220  1.7831332
      B           -1.5914170 -1.6746325
      B_Eneutral   0.2736373  0.5774843
      B_Eaccuracy -0.7775480 -1.8565900
      B_lRright    0.1379297  0.5908926
      A           -1.9607101 -2.4159174
      t0          -0.6633743 -0.9088031
      
      
      $stage
      [1] "init"
      
      $subj_ll
                 [,1]
      as1t -14281.188
      bd6t  -8083.852
      
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
      1        1     as1t    speed  left  left 1.2896914  0 Inf  0 Inf
      3        2     as1t    speed  left  left 1.0285371  0 Inf  0 Inf
      5        3     as1t    speed  left  left 1.2373603  0 Inf  0 Inf
      7        4     as1t    speed  left  left 1.1524328  0 Inf  0 Inf
      9        5     as1t    speed  left  left 0.8333630  0 Inf  0 Inf
      11       6     as1t    speed  left  left 1.7243277  0 Inf  0 Inf
      13       7     as1t    speed  left  left 1.3136383  0 Inf  0 Inf
      15       8     as1t    speed  left  left 1.1226859  0 Inf  0 Inf
      17       9     as1t    speed  left  left 1.2613008  0 Inf  0 Inf
      19      10     as1t    speed  left  left 1.3511772  0 Inf  0 Inf
      21      11     as1t  neutral  left right 1.4992827  0 Inf  0 Inf
      23      12     as1t  neutral  left  left 1.5886886  0 Inf  0 Inf
      25      13     as1t  neutral  left  left 1.6352497  0 Inf  0 Inf
      27      14     as1t  neutral  left  left 2.2574440  0 Inf  0 Inf
      29      15     as1t  neutral  left  left 1.6082515  0 Inf  0 Inf
      31      16     as1t  neutral  left right 1.3587948  0 Inf  0 Inf
      33      17     as1t  neutral  left right 3.1745879  0 Inf  0 Inf
      35      18     as1t  neutral  left  left 1.7081291  0 Inf  0 Inf
      37      19     as1t  neutral  left  left 1.9135045  0 Inf  0 Inf
      39      20     as1t  neutral  left  left 1.7989982  0 Inf  0 Inf
      41      21     as1t accuracy  left  left 1.4079495  0 Inf  0 Inf
      43      22     as1t accuracy  left  left 1.7598214  0 Inf  0 Inf
      45      23     as1t accuracy  left  left 1.8853903  0 Inf  0 Inf
      47      24     as1t accuracy  left  left 0.9389688  0 Inf  0 Inf
      49      25     as1t accuracy  left  left 1.5189026  0 Inf  0 Inf
      51      26     as1t accuracy  left right 2.0054210  0 Inf  0 Inf
      53      27     as1t accuracy  left  left 1.2857265  0 Inf  0 Inf
      55      28     as1t accuracy  left right 1.4371268  0 Inf  0 Inf
      57      29     as1t accuracy  left  left 1.2420155  0 Inf  0 Inf
      59      30     as1t accuracy  left  left 1.6205307  0 Inf  0 Inf
      61      31     as1t    speed right  left 1.0359667  0 Inf  0 Inf
      63      32     as1t    speed right right 1.0424088  0 Inf  0 Inf
      65      33     as1t    speed right right 1.2043272  0 Inf  0 Inf
      67      34     as1t    speed right right 1.2190427  0 Inf  0 Inf
      69      35     as1t    speed right right 1.0134876  0 Inf  0 Inf
      71      36     as1t    speed right  left 1.2197044  0 Inf  0 Inf
      73      37     as1t    speed right right 1.3743929  0 Inf  0 Inf
      75      38     as1t    speed right right 1.7841597  0 Inf  0 Inf
      77      39     as1t    speed right right 1.2256365  0 Inf  0 Inf
      79      40     as1t    speed right right 1.7222360  0 Inf  0 Inf
      81      41     as1t  neutral right  left 2.2956287  0 Inf  0 Inf
      83      42     as1t  neutral right right 2.7010827  0 Inf  0 Inf
      85      43     as1t  neutral right right 1.6281372  0 Inf  0 Inf
      87      44     as1t  neutral right right 2.1973495  0 Inf  0 Inf
      89      45     as1t  neutral right right 1.4796700  0 Inf  0 Inf
      91      46     as1t  neutral right right 2.9906534  0 Inf  0 Inf
      93      47     as1t  neutral right  left 1.3987644  0 Inf  0 Inf
      95      48     as1t  neutral right  left 2.1625135  0 Inf  0 Inf
      97      49     as1t  neutral right  left 2.4094435  0 Inf  0 Inf
      99      50     as1t  neutral right right 1.3683235  0 Inf  0 Inf
      101     51     as1t accuracy right right 1.5178946  0 Inf  0 Inf
      103     52     as1t accuracy right right 1.7793677  0 Inf  0 Inf
      105     53     as1t accuracy right right 1.2988182  0 Inf  0 Inf
      107     54     as1t accuracy right right 1.6421192  0 Inf  0 Inf
      109     55     as1t accuracy right right 1.4054138  0 Inf  0 Inf
      111     56     as1t accuracy right  left 1.6433920  0 Inf  0 Inf
      113     57     as1t accuracy right right 1.7162056  0 Inf  0 Inf
      115     58     as1t accuracy right right 1.2697505  0 Inf  0 Inf
      117     59     as1t accuracy right right 1.4206307  0 Inf  0 Inf
      119     60     as1t accuracy right  left 1.7459946  0 Inf  0 Inf
      121      1     bd6t    speed  left  left 1.3714492  0 Inf  0 Inf
      123      2     bd6t    speed  left  left 1.4364864  0 Inf  0 Inf
      125      3     bd6t    speed  left  left 1.1778969  0 Inf  0 Inf
      127      4     bd6t    speed  left right 0.9819309  0 Inf  0 Inf
      129      5     bd6t    speed  left  left 1.1987254  0 Inf  0 Inf
      131      6     bd6t    speed  left  left 1.0798519  0 Inf  0 Inf
      133      7     bd6t    speed  left  left 1.0386992  0 Inf  0 Inf
      135      8     bd6t    speed  left  left 1.3614354  0 Inf  0 Inf
      137      9     bd6t    speed  left  left 1.2041459  0 Inf  0 Inf
      139     10     bd6t    speed  left  left 1.0689792  0 Inf  0 Inf
      141     11     bd6t  neutral  left right 1.0109397  0 Inf  0 Inf
      143     12     bd6t  neutral  left right 1.6625619  0 Inf  0 Inf
      145     13     bd6t  neutral  left right 1.8192494  0 Inf  0 Inf
      147     14     bd6t  neutral  left right 1.8034078  0 Inf  0 Inf
      149     15     bd6t  neutral  left  left 2.4293131  0 Inf  0 Inf
      151     16     bd6t  neutral  left  left 2.7545452  0 Inf  0 Inf
      153     17     bd6t  neutral  left  left 1.9841839  0 Inf  0 Inf
      155     18     bd6t  neutral  left  left 1.5321501  0 Inf  0 Inf
      157     19     bd6t  neutral  left right 1.4199504  0 Inf  0 Inf
      159     20     bd6t  neutral  left  left 1.9732700  0 Inf  0 Inf
      161     21     bd6t accuracy  left  left 2.1056345  0 Inf  0 Inf
      163     22     bd6t accuracy  left  left 1.7901854  0 Inf  0 Inf
      165     23     bd6t accuracy  left right 1.5089911  0 Inf  0 Inf
      167     24     bd6t accuracy  left right 1.2822764  0 Inf  0 Inf
      169     25     bd6t accuracy  left right 1.1000061  0 Inf  0 Inf
      171     26     bd6t accuracy  left  left 2.2118107  0 Inf  0 Inf
      173     27     bd6t accuracy  left  left 1.9370737  0 Inf  0 Inf
      175     28     bd6t accuracy  left  left 1.3774867  0 Inf  0 Inf
      177     29     bd6t accuracy  left  left 1.2905306  0 Inf  0 Inf
      179     30     bd6t accuracy  left right 1.1681654  0 Inf  0 Inf
      181     31     bd6t    speed right right 1.3152703  0 Inf  0 Inf
      183     32     bd6t    speed right right 1.4990331  0 Inf  0 Inf
      185     33     bd6t    speed right  left 1.0154286  0 Inf  0 Inf
      187     34     bd6t    speed right right 0.8989288  0 Inf  0 Inf
      189     35     bd6t    speed right right 1.1673860  0 Inf  0 Inf
      191     36     bd6t    speed right right 1.0976647  0 Inf  0 Inf
      193     37     bd6t    speed right right 1.1752087  0 Inf  0 Inf
      195     38     bd6t    speed right right 1.3387905  0 Inf  0 Inf
      197     39     bd6t    speed right right 1.0504423  0 Inf  0 Inf
      199     40     bd6t    speed right right 1.0435325  0 Inf  0 Inf
      201     41     bd6t  neutral right right 2.0063859  0 Inf  0 Inf
      203     42     bd6t  neutral right right 1.3282215  0 Inf  0 Inf
      205     43     bd6t  neutral right right 2.3469504  0 Inf  0 Inf
      207     44     bd6t  neutral right right 1.2746677  0 Inf  0 Inf
      209     45     bd6t  neutral right right 2.4306889  0 Inf  0 Inf
      211     46     bd6t  neutral right right 2.3053142  0 Inf  0 Inf
      213     47     bd6t  neutral right right 2.0149578  0 Inf  0 Inf
      215     48     bd6t  neutral right right 2.2116807  0 Inf  0 Inf
      217     49     bd6t  neutral right right 1.8316487  0 Inf  0 Inf
      219     50     bd6t  neutral right right 1.4542518  0 Inf  0 Inf
      221     51     bd6t accuracy right right 1.4135664  0 Inf  0 Inf
      223     52     bd6t accuracy right  left 1.3147865  0 Inf  0 Inf
      225     53     bd6t accuracy right  left 0.8385127  0 Inf  0 Inf
      227     54     bd6t accuracy right right 2.0408400  0 Inf  0 Inf
      229     55     bd6t accuracy right right 1.4847450  0 Inf  0 Inf
      231     56     bd6t accuracy right  left 2.1690447  0 Inf  0 Inf
      233     57     bd6t accuracy right  left 1.6631792  0 Inf  0 Inf
      235     58     bd6t accuracy right right 1.4560110  0 Inf  0 Inf
      237     59     bd6t accuracy right right 1.1012568  0 Inf  0 Inf
      239     60     bd6t accuracy right right 1.4341022  0 Inf  0 Inf

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
      v           -0.8659366 -0.8274847
      v_lMd        0.8828966  0.5761732
      s_lMTRUE     1.4190220  1.7831332
      B           -1.5914170 -1.6746325
      B_Eneutral   0.2736373  0.5774843
      B_Eaccuracy -0.7775480 -1.8565900
      B_lRright    0.1379297  0.5908926
      A           -1.9607101 -2.4159174
      t0          -0.6633743 -0.9088031
      
      
      $stage
      [1] "init"
      
      $subj_ll
                 [,1]
      as1t -13319.641
      bd6t  -8344.067
      
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
      1        1     as1t    speed  left  left 1.0347592  0 Inf  0 Inf
      3        2     as1t    speed  left  left 0.7821024  0 Inf  0 Inf
      5        3     as1t    speed  left  left 0.9624820  0 Inf  0 Inf
      7        4     as1t    speed  left  left 0.7665585  0 Inf  0 Inf
      9        5     as1t    speed  left  left 0.9098368  0 Inf  0 Inf
      11       6     as1t    speed  left  left 1.2086700  0 Inf  0 Inf
      13       7     as1t    speed  left right 1.0375593  0 Inf  0 Inf
      15       8     as1t    speed  left  left 1.2600122  0 Inf  0 Inf
      17       9     as1t    speed  left  left 1.4212036  0 Inf  0 Inf
      19      10     as1t    speed  left right 0.9174280  0 Inf  0 Inf
      21      11     as1t  neutral  left  left 1.5491453  0 Inf  0 Inf
      23      12     as1t  neutral  left  left 1.0641637  0 Inf  0 Inf
      25      13     as1t  neutral  left  left 1.6852806  0 Inf  0 Inf
      27      14     as1t  neutral  left  left 1.4563528  0 Inf  0 Inf
      29      15     as1t  neutral  left  left 1.6797563  0 Inf  0 Inf
      31      16     as1t  neutral  left right 1.3651094  0 Inf  0 Inf
      33      17     as1t  neutral  left  left 1.5790399  0 Inf  0 Inf
      35      18     as1t  neutral  left  left 1.5117983  0 Inf  0 Inf
      37      19     as1t  neutral  left right 1.0361019  0 Inf  0 Inf
      39      20     as1t  neutral  left  left 1.4389256  0 Inf  0 Inf
      41      21     as1t accuracy  left  left 1.3294067  0 Inf  0 Inf
      43      22     as1t accuracy  left  left 1.3568084  0 Inf  0 Inf
      45      23     as1t accuracy  left  left 1.3457983  0 Inf  0 Inf
      47      24     as1t accuracy  left  left 1.3768192  0 Inf  0 Inf
      49      25     as1t accuracy  left  left 1.1842639  0 Inf  0 Inf
      51      26     as1t accuracy  left  left 1.1972452  0 Inf  0 Inf
      53      27     as1t accuracy  left  left 1.2529548  0 Inf  0 Inf
      55      28     as1t accuracy  left  left 0.7966244  0 Inf  0 Inf
      57      29     as1t accuracy  left  left 1.0104944  0 Inf  0 Inf
      59      30     as1t accuracy  left  left 1.1564672  0 Inf  0 Inf
      61      31     as1t    speed right  left 1.1112721  0 Inf  0 Inf
      63      32     as1t    speed right right 1.1805067  0 Inf  0 Inf
      65      33     as1t    speed right  left 0.8850711  0 Inf  0 Inf
      67      34     as1t    speed right right 1.0790560  0 Inf  0 Inf
      69      35     as1t    speed right right 1.1603104  0 Inf  0 Inf
      71      36     as1t    speed right right 0.9123440  0 Inf  0 Inf
      73      37     as1t    speed right right 1.0184145  0 Inf  0 Inf
      75      38     as1t    speed right right 0.9051514  0 Inf  0 Inf
      77      39     as1t    speed right right 0.9930121  0 Inf  0 Inf
      79      40     as1t    speed right right 1.0800397  0 Inf  0 Inf
      81      41     as1t  neutral right right 1.3587212  0 Inf  0 Inf
      83      42     as1t  neutral right right 1.5467528  0 Inf  0 Inf
      85      43     as1t  neutral right right 1.4929551  0 Inf  0 Inf
      87      44     as1t  neutral right right 1.3161518  0 Inf  0 Inf
      89      45     as1t  neutral right right 1.5583743  0 Inf  0 Inf
      91      46     as1t  neutral right right 1.4229619  0 Inf  0 Inf
      93      47     as1t  neutral right right 1.5496307  0 Inf  0 Inf
      95      48     as1t  neutral right right 1.5800015  0 Inf  0 Inf
      97      49     as1t  neutral right right 1.3601328  0 Inf  0 Inf
      99      50     as1t  neutral right right 1.8111763  0 Inf  0 Inf
      101     51     as1t accuracy right right 1.0285247  0 Inf  0 Inf
      103     52     as1t accuracy right right 1.3594421  0 Inf  0 Inf
      105     53     as1t accuracy right right 1.1771948  0 Inf  0 Inf
      107     54     as1t accuracy right right 1.3919244  0 Inf  0 Inf
      109     55     as1t accuracy right right 1.0275638  0 Inf  0 Inf
      111     56     as1t accuracy right right 1.1989792  0 Inf  0 Inf
      113     57     as1t accuracy right  left 1.2949567  0 Inf  0 Inf
      115     58     as1t accuracy right right 1.1265176  0 Inf  0 Inf
      117     59     as1t accuracy right right 1.1076321  0 Inf  0 Inf
      119     60     as1t accuracy right right 1.4518215  0 Inf  0 Inf
      121      1     bd6t    speed  left  left 1.3922462  0 Inf  0 Inf
      123      2     bd6t    speed  left  left 0.6925170  0 Inf  0 Inf
      125      3     bd6t    speed  left  left 1.3164219  0 Inf  0 Inf
      127      4     bd6t    speed  left  left 0.9731763  0 Inf  0 Inf
      129      5     bd6t    speed  left  left 0.7860885  0 Inf  0 Inf
      131      6     bd6t    speed  left  left 0.9035583  0 Inf  0 Inf
      133      7     bd6t    speed  left  left 0.9123796  0 Inf  0 Inf
      135      8     bd6t    speed  left  left 1.2288798  0 Inf  0 Inf
      137      9     bd6t    speed  left  left 1.2575851  0 Inf  0 Inf
      139     10     bd6t    speed  left  left 1.2119507  0 Inf  0 Inf
      141     11     bd6t  neutral  left  left 1.5226392  0 Inf  0 Inf
      143     12     bd6t  neutral  left  left 1.5671165  0 Inf  0 Inf
      145     13     bd6t  neutral  left  left 1.2726183  0 Inf  0 Inf
      147     14     bd6t  neutral  left  left 1.3531255  0 Inf  0 Inf
      149     15     bd6t  neutral  left  left 1.3484313  0 Inf  0 Inf
      151     16     bd6t  neutral  left  left 1.5242674  0 Inf  0 Inf
      153     17     bd6t  neutral  left  left 1.5889658  0 Inf  0 Inf
      155     18     bd6t  neutral  left  left 1.5116453  0 Inf  0 Inf
      157     19     bd6t  neutral  left  left 1.7800438  0 Inf  0 Inf
      159     20     bd6t  neutral  left right 1.4899301  0 Inf  0 Inf
      161     21     bd6t accuracy  left  left 1.5127973  0 Inf  0 Inf
      163     22     bd6t accuracy  left  left 1.2771573  0 Inf  0 Inf
      165     23     bd6t accuracy  left  left 1.1142639  0 Inf  0 Inf
      167     24     bd6t accuracy  left  left 1.2217839  0 Inf  0 Inf
      169     25     bd6t accuracy  left  left 1.1251490  0 Inf  0 Inf
      171     26     bd6t accuracy  left  left 1.1405746  0 Inf  0 Inf
      173     27     bd6t accuracy  left  left 1.7402099  0 Inf  0 Inf
      175     28     bd6t accuracy  left  left 1.4027512  0 Inf  0 Inf
      177     29     bd6t accuracy  left  left 1.6235649  0 Inf  0 Inf
      179     30     bd6t accuracy  left  left 1.1659195  0 Inf  0 Inf
      181     31     bd6t    speed right right 0.9975004  0 Inf  0 Inf
      183     32     bd6t    speed right right 1.0920639  0 Inf  0 Inf
      185     33     bd6t    speed right right 0.9999527  0 Inf  0 Inf
      187     34     bd6t    speed right right 0.9085117  0 Inf  0 Inf
      189     35     bd6t    speed right  left 0.6663475  0 Inf  0 Inf
      191     36     bd6t    speed right right 0.8962133  0 Inf  0 Inf
      193     37     bd6t    speed right right 0.9046192  0 Inf  0 Inf
      195     38     bd6t    speed right right 0.9463857  0 Inf  0 Inf
      197     39     bd6t    speed right right 1.2566465  0 Inf  0 Inf
      199     40     bd6t    speed right right 0.9136549  0 Inf  0 Inf
      201     41     bd6t  neutral right right 1.6188236  0 Inf  0 Inf
      203     42     bd6t  neutral right right 1.3315123  0 Inf  0 Inf
      205     43     bd6t  neutral right right 1.5011992  0 Inf  0 Inf
      207     44     bd6t  neutral right right 1.4922071  0 Inf  0 Inf
      209     45     bd6t  neutral right right 1.7380602  0 Inf  0 Inf
      211     46     bd6t  neutral right right 1.4534910  0 Inf  0 Inf
      213     47     bd6t  neutral right right 1.6021723  0 Inf  0 Inf
      215     48     bd6t  neutral right right 1.7835590  0 Inf  0 Inf
      217     49     bd6t  neutral right right 1.2510612  0 Inf  0 Inf
      219     50     bd6t  neutral right right 1.4445173  0 Inf  0 Inf
      221     51     bd6t accuracy right right 1.4527156  0 Inf  0 Inf
      223     52     bd6t accuracy right right 1.2398796  0 Inf  0 Inf
      225     53     bd6t accuracy right right 1.3300044  0 Inf  0 Inf
      227     54     bd6t accuracy right right 1.2901084  0 Inf  0 Inf
      229     55     bd6t accuracy right right 1.1579839  0 Inf  0 Inf
      231     56     bd6t accuracy right right 1.5016884  0 Inf  0 Inf
      233     57     bd6t accuracy right right 1.1262655  0 Inf  0 Inf
      235     58     bd6t accuracy right right 1.3265131  0 Inf  0 Inf
      237     59     bd6t accuracy right right 1.0309981  0 Inf  0 Inf
      239     60     bd6t accuracy right right 1.3712989  0 Inf  0 Inf

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
      
                          as1t        bd6t
      v_Sleft     -0.795745564 -0.61298014
      v_Sright     0.831296796  0.53982132
      a            1.424407847  1.19546709
      a_Eneutral  -1.922987294 -1.62502932
      a_Eaccuracy  0.049542401  0.03054263
      t0          -1.396100362 -1.38269638
      Z           -0.005054203 -0.21436007
      sv          -2.197044667 -2.08693786
      SZ          -0.436668005 -0.46579217
      
      
      $stage
      [1] "init"
      
      $subj_ll
                [,1]
      as1t -3601.088
      bd6t -3112.817
      
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
      1        1     as1t    speed  left right 0.4918896  0 Inf  0 Inf
      2        2     as1t    speed  left  left 2.1635739  0 Inf  0 Inf
      3        3     as1t    speed  left  left 0.2698151  0 Inf  0 Inf
      4        4     as1t    speed  left  left 0.6766773  0 Inf  0 Inf
      5        5     as1t    speed  left  left 0.3295339  0 Inf  0 Inf
      6        6     as1t    speed  left  left 0.9199389  0 Inf  0 Inf
      7        7     as1t    speed  left  left 0.5543087  0 Inf  0 Inf
      8        8     as1t    speed  left  left 0.3507439  0 Inf  0 Inf
      9        9     as1t    speed  left right 0.6624050  0 Inf  0 Inf
      10      10     as1t    speed  left  left 0.4215988  0 Inf  0 Inf
      21      11     as1t  neutral  left  left 0.9370477  0 Inf  0 Inf
      22      12     as1t  neutral  left  left 0.5277394  0 Inf  0 Inf
      23      13     as1t  neutral  left  left 0.6777917  0 Inf  0 Inf
      24      14     as1t  neutral  left  left 0.7699879  0 Inf  0 Inf
      25      15     as1t  neutral  left  left 1.3854253  0 Inf  0 Inf
      26      16     as1t  neutral  left  left 0.6911855  0 Inf  0 Inf
      27      17     as1t  neutral  left  left 0.6287586  0 Inf  0 Inf
      28      18     as1t  neutral  left  left 0.3620831  0 Inf  0 Inf
      29      19     as1t  neutral  left  left 0.8407541  0 Inf  0 Inf
      30      20     as1t  neutral  left right 1.0489656  0 Inf  0 Inf
      41      21     as1t accuracy  left  left 0.4813305  0 Inf  0 Inf
      42      22     as1t accuracy  left  left 1.6113097  0 Inf  0 Inf
      43      23     as1t accuracy  left  left 0.5893647  0 Inf  0 Inf
      44      24     as1t accuracy  left  left 0.3988592  0 Inf  0 Inf
      45      25     as1t accuracy  left  left 0.8051490  0 Inf  0 Inf
      46      26     as1t accuracy  left  left 1.0891681  0 Inf  0 Inf
      47      27     as1t accuracy  left  left 1.0257608  0 Inf  0 Inf
      48      28     as1t accuracy  left right 1.5787319  0 Inf  0 Inf
      49      29     as1t accuracy  left  left 0.5356605  0 Inf  0 Inf
      50      30     as1t accuracy  left right 0.7410944  0 Inf  0 Inf
      61      31     as1t    speed right right 1.5454972  0 Inf  0 Inf
      62      32     as1t    speed right right 2.9916319  0 Inf  0 Inf
      63      33     as1t    speed right right 0.5117038  0 Inf  0 Inf
      64      34     as1t    speed right right 0.3290761  0 Inf  0 Inf
      65      35     as1t    speed right right 0.7918841  0 Inf  0 Inf
      66      36     as1t    speed right right 0.4269244  0 Inf  0 Inf
      67      37     as1t    speed right right 0.4079138  0 Inf  0 Inf
      68      38     as1t    speed right right 0.3534925  0 Inf  0 Inf
      69      39     as1t    speed right  left 0.8016991  0 Inf  0 Inf
      70      40     as1t    speed right right 1.2769941  0 Inf  0 Inf
      81      41     as1t  neutral right right 0.5308851  0 Inf  0 Inf
      82      42     as1t  neutral right right 0.5480776  0 Inf  0 Inf
      83      43     as1t  neutral right right 0.8859953  0 Inf  0 Inf
      84      44     as1t  neutral right right 1.0519749  0 Inf  0 Inf
      85      45     as1t  neutral right right 0.3584496  0 Inf  0 Inf
      86      46     as1t  neutral right right 0.4053901  0 Inf  0 Inf
      87      47     as1t  neutral right right 0.5509626  0 Inf  0 Inf
      88      48     as1t  neutral right right 1.0218690  0 Inf  0 Inf
      89      49     as1t  neutral right right 0.4517298  0 Inf  0 Inf
      90      50     as1t  neutral right right 1.3101657  0 Inf  0 Inf
      101     51     as1t accuracy right right 0.3339661  0 Inf  0 Inf
      102     52     as1t accuracy right right 0.6139724  0 Inf  0 Inf
      103     53     as1t accuracy right right 0.3721116  0 Inf  0 Inf
      104     54     as1t accuracy right right 0.7174872  0 Inf  0 Inf
      105     55     as1t accuracy right right 0.8513820  0 Inf  0 Inf
      106     56     as1t accuracy right right 5.9493359  0 Inf  0 Inf
      107     57     as1t accuracy right right 0.6835362  0 Inf  0 Inf
      108     58     as1t accuracy right right 0.3788342  0 Inf  0 Inf
      109     59     as1t accuracy right  left 0.4092784  0 Inf  0 Inf
      110     60     as1t accuracy right right 0.4775454  0 Inf  0 Inf
      11       1     bd6t    speed  left right 0.4325714  0 Inf  0 Inf
      12       2     bd6t    speed  left  left 0.5118116  0 Inf  0 Inf
      13       3     bd6t    speed  left  left 0.3219160  0 Inf  0 Inf
      14       4     bd6t    speed  left right 0.4170058  0 Inf  0 Inf
      15       5     bd6t    speed  left  left 1.2714175  0 Inf  0 Inf
      16       6     bd6t    speed  left  left 0.4190258  0 Inf  0 Inf
      17       7     bd6t    speed  left right 0.8166139  0 Inf  0 Inf
      18       8     bd6t    speed  left  left 0.8160875  0 Inf  0 Inf
      19       9     bd6t    speed  left right 1.0603326  0 Inf  0 Inf
      20      10     bd6t    speed  left  left 0.3661472  0 Inf  0 Inf
      31      11     bd6t  neutral  left  left 0.4651464  0 Inf  0 Inf
      32      12     bd6t  neutral  left  left 0.4881943  0 Inf  0 Inf
      33      13     bd6t  neutral  left  left 0.6469333  0 Inf  0 Inf
      34      14     bd6t  neutral  left right 1.3517407  0 Inf  0 Inf
      35      15     bd6t  neutral  left  left 0.4617880  0 Inf  0 Inf
      36      16     bd6t  neutral  left right 5.3458456  0 Inf  0 Inf
      37      17     bd6t  neutral  left  left 0.5566515  0 Inf  0 Inf
      38      18     bd6t  neutral  left  left 2.6954209  0 Inf  0 Inf
      39      19     bd6t  neutral  left  left 0.9965662  0 Inf  0 Inf
      40      20     bd6t  neutral  left  left 0.3362758  0 Inf  0 Inf
      51      21     bd6t accuracy  left  left 0.5661682  0 Inf  0 Inf
      52      22     bd6t accuracy  left right 2.0680112  0 Inf  0 Inf
      53      23     bd6t accuracy  left  left 1.7238075  0 Inf  0 Inf
      54      24     bd6t accuracy  left  left 0.4823011  0 Inf  0 Inf
      55      25     bd6t accuracy  left  left 1.0530585  0 Inf  0 Inf
      56      26     bd6t accuracy  left right 1.6529589  0 Inf  0 Inf
      57      27     bd6t accuracy  left  left 0.4848049  0 Inf  0 Inf
      58      28     bd6t accuracy  left  left 0.7076725  0 Inf  0 Inf
      59      29     bd6t accuracy  left  left 0.4436648  0 Inf  0 Inf
      60      30     bd6t accuracy  left right 0.9370205  0 Inf  0 Inf
      71      31     bd6t    speed right  left 0.7293802  0 Inf  0 Inf
      72      32     bd6t    speed right right 0.3304065  0 Inf  0 Inf
      73      33     bd6t    speed right right 0.6348638  0 Inf  0 Inf
      74      34     bd6t    speed right right 0.4922937  0 Inf  0 Inf
      75      35     bd6t    speed right right 0.3681049  0 Inf  0 Inf
      76      36     bd6t    speed right right 0.7115350  0 Inf  0 Inf
      77      37     bd6t    speed right right 0.4972506  0 Inf  0 Inf
      78      38     bd6t    speed right right 1.8552103  0 Inf  0 Inf
      79      39     bd6t    speed right  left 2.1091180  0 Inf  0 Inf
      80      40     bd6t    speed right right 0.5811276  0 Inf  0 Inf
      91      41     bd6t  neutral right right 0.5548770  0 Inf  0 Inf
      92      42     bd6t  neutral right right 0.4551236  0 Inf  0 Inf
      93      43     bd6t  neutral right right 0.7856112  0 Inf  0 Inf
      94      44     bd6t  neutral right right 0.5094467  0 Inf  0 Inf
      95      45     bd6t  neutral right right 0.4823116  0 Inf  0 Inf
      96      46     bd6t  neutral right  left 0.6135672  0 Inf  0 Inf
      97      47     bd6t  neutral right right 2.3302518  0 Inf  0 Inf
      98      48     bd6t  neutral right right 1.5269199  0 Inf  0 Inf
      99      49     bd6t  neutral right right 0.8128263  0 Inf  0 Inf
      100     50     bd6t  neutral right  left 3.7695008  0 Inf  0 Inf
      111     51     bd6t accuracy right right 0.3649574  0 Inf  0 Inf
      112     52     bd6t accuracy right  left 1.3181878  0 Inf  0 Inf
      113     53     bd6t accuracy right right 0.9815775  0 Inf  0 Inf
      114     54     bd6t accuracy right right 0.5260662  0 Inf  0 Inf
      115     55     bd6t accuracy right right 0.4114870  0 Inf  0 Inf
      116     56     bd6t accuracy right  left 0.6621542  0 Inf  0 Inf
      117     57     bd6t accuracy right  left 1.3903490  0 Inf  0 Inf
      118     58     bd6t accuracy right right 0.7186431  0 Inf  0 Inf
      119     59     bd6t accuracy right right 0.8569282  0 Inf  0 Inf
      120     60     bd6t accuracy right  left 0.6147701  0 Inf  0 Inf

