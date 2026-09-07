# mapped_pars

    Code
      mapped_pars(des)
    Output
      $v 
        covariate           E         
       -0.513716908811104  speed     : v - 0.514 * v_covariate 
       0.601813121611908   neutral   : v + 0.602 * v_covariate + v_Eneutral + 0.602 * v_covariate:Eneutral 
       -1.54982186331112   accuracy  : v - 1.550 * v_covariate + v_Eaccuracy - 1.550 * v_covariate:Eaccuracy
       -1.70962282511647   speed     : v - 1.710 * v_covariate 
       0.7701488026301     neutral   : v + 0.770 * v_covariate + v_Eneutral + 0.770 * v_covariate:Eneutral 
       -0.716873043602934  accuracy  : v - 0.717 * v_covariate + v_Eaccuracy - 0.717 * v_covariate:Eaccuracy
      
      $B 
        E         
       speed     : exp(B )
       neutral   : exp(B + B_Eneutral )
       accuracy  : exp(B + B_Eaccuracy)
      
      $t0 
        S      
       left   : exp(t0 )
       right  : exp(t0 + t0_Sright)
      

---

    Code
      mapped_pars(des, p_vector = rnorm(length(sampled_pars(des))))
    Output
                 E     S covariate    lR          v sv     B A    t0 eta pContaminant
      1      speed  left         3  left     -2.851  1 3.048 0 0.292   0            0
      2      speed  left         3 right     -2.851  1 3.048 0 0.292   0            0
      3      speed  left      1853  left  -2520.117  1 3.048 0 0.292   0            0
      4      speed  left      1853 right  -2520.117  1 3.048 0 0.292   0            0
      5      speed  left      3463  left  -4710.819  1 3.048 0 0.292   0            0
      6      speed  left      3463 right  -4710.819  1 3.048 0 0.292   0            0
      7      speed  left      5353  left  -7282.513  1 3.048 0 0.292   0            0
      8      speed  left      5353 right  -7282.513  1 3.048 0 0.292   0            0
      9      speed  left      7153  left  -9731.745  1 3.048 0 0.292   0            0
      10     speed  left      7153 right  -9731.745  1 3.048 0 0.292   0            0
      11     speed  left      9061  left -12327.931  1 3.048 0 0.292   0            0
      12     speed  left      9061 right -12327.931  1 3.048 0 0.292   0            0
      13     speed  left     10735  left -14605.716  1 3.048 0 0.292   0            0
      14     speed  left     10735 right -14605.716  1 3.048 0 0.292   0            0
      15     speed  left     12327  left -16771.926  1 3.048 0 0.292   0            0
      16     speed  left     12327 right -16771.926  1 3.048 0 0.292   0            0
      17     speed  left     14033  left -19093.254  1 3.048 0 0.292   0            0
      18     speed  left     14033 right -19093.254  1 3.048 0 0.292   0            0
      19     speed  left     15809  left -21509.829  1 3.048 0 0.292   0            0
      20     speed  left     15809 right -21509.829  1 3.048 0 0.292   0            0
      21   neutral  left         6  left     -0.942  1 1.793 0 0.292   0            0
      22   neutral  left         6 right     -0.942  1 1.793 0 0.292   0            0
      23   neutral  left      1680  left   -517.645  1 1.793 0 0.292   0            0
      24   neutral  left      1680 right   -517.645  1 1.793 0 0.292   0            0
      25   neutral  left      3409  left  -1051.324  1 1.793 0 0.292   0            0
      26   neutral  left      3409 right  -1051.324  1 1.793 0 0.292   0            0
      27   neutral  left      5071  left  -1564.323  1 1.793 0 0.292   0            0
      28   neutral  left      5071 right  -1564.323  1 1.793 0 0.292   0            0
      29   neutral  left      6804  left  -2099.237  1 1.793 0 0.292   0            0
      30   neutral  left      6804 right  -2099.237  1 1.793 0 0.292   0            0
      31   neutral  left      8739  left  -2696.501  1 1.793 0 0.292   0            0
      32   neutral  left      8739 right  -2696.501  1 1.793 0 0.292   0            0
      33   neutral  left     10461  left  -3228.019  1 1.793 0 0.292   0            0
      34   neutral  left     10461 right  -3228.019  1 1.793 0 0.292   0            0
      35   neutral  left     12287  left  -3791.639  1 1.793 0 0.292   0            0
      36   neutral  left     12287 right  -3791.639  1 1.793 0 0.292   0            0
      37   neutral  left     14137  left  -4362.667  1 1.793 0 0.292   0            0
      38   neutral  left     14137 right  -4362.667  1 1.793 0 0.292   0            0
      39   neutral  left     15816  left  -4880.913  1 1.793 0 0.292   0            0
      40   neutral  left     15816 right  -4880.913  1 1.793 0 0.292   0            0
      41  accuracy  left         5  left    -11.877  1 3.051 0 0.292   0            0
      42  accuracy  left         5 right    -11.877  1 3.051 0 0.292   0            0
      43  accuracy  left      1764  left  -4228.095  1 3.051 0 0.292   0            0
      44  accuracy  left      1764 right  -4228.095  1 3.051 0 0.292   0            0
      45  accuracy  left      3543  left  -8492.252  1 3.051 0 0.292   0            0
      46  accuracy  left      3543 right  -8492.252  1 3.051 0 0.292   0            0
      47  accuracy  left      5235  left -12547.875  1 3.051 0 0.292   0            0
      48  accuracy  left      5235 right -12547.875  1 3.051 0 0.292   0            0
      49  accuracy  left      7069  left -16943.863  1 3.051 0 0.292   0            0
      50  accuracy  left      7069 right -16943.863  1 3.051 0 0.292   0            0
      51  accuracy  left      8725  left -20913.196  1 3.051 0 0.292   0            0
      52  accuracy  left      8725 right -20913.196  1 3.051 0 0.292   0            0
      53  accuracy  left     10462  left -25076.681  1 3.051 0 0.292   0            0
      54  accuracy  left     10462 right -25076.681  1 3.051 0 0.292   0            0
      55  accuracy  left     12289  left -29455.891  1 3.051 0 0.292   0            0
      56  accuracy  left     12289 right -29455.891  1 3.051 0 0.292   0            0
      57  accuracy  left     14136  left -33883.040  1 3.051 0 0.292   0            0
      58  accuracy  left     14136 right -33883.040  1 3.051 0 0.292   0            0
      59  accuracy  left     15817  left -37912.296  1 3.051 0 0.292   0            0
      60  accuracy  left     15817 right -37912.296  1 3.051 0 0.292   0            0
      61     speed right         2  left     -1.490  1 3.048 0 0.189   0            0
      62     speed right         2 right     -1.490  1 3.048 0 0.189   0            0
      63     speed right      1812  left  -2464.329  1 3.048 0 0.189   0            0
      64     speed right      1812 right  -2464.329  1 3.048 0 0.189   0            0
      65     speed right      3494  left  -4753.000  1 3.048 0 0.189   0            0
      66     speed right      3494 right  -4753.000  1 3.048 0 0.189   0            0
      67     speed right      5243  left  -7132.837  1 3.048 0 0.189   0            0
      68     speed right      5243 right  -7132.837  1 3.048 0 0.189   0            0
      69     speed right      7089  left  -9644.661  1 3.048 0 0.189   0            0
      70     speed right      7089 right  -9644.661  1 3.048 0 0.189   0            0
      71     speed right      8655  left -11775.493  1 3.048 0 0.189   0            0
      72     speed right      8655 right -11775.493  1 3.048 0 0.189   0            0
      73     speed right     10545  left -14347.186  1 3.048 0 0.189   0            0
      74     speed right     10545 right -14347.186  1 3.048 0 0.189   0            0
      75     speed right     12274  left -16699.810  1 3.048 0 0.189   0            0
      76     speed right     12274 right -16699.810  1 3.048 0 0.189   0            0
      77     speed right     13900  left -18912.283  1 3.048 0 0.189   0            0
      78     speed right     13900 right -18912.283  1 3.048 0 0.189   0            0
      79     speed right     15818  left -21522.075  1 3.048 0 0.189   0            0
      80     speed right     15818 right -21522.075  1 3.048 0 0.189   0            0
      81   neutral right         8  left     -1.559  1 1.793 0 0.189   0            0
      82   neutral right         8 right     -1.559  1 1.793 0 0.189   0            0
      83   neutral right      1733  left   -534.004  1 1.793 0 0.189   0            0
      84   neutral right      1733 right   -534.004  1 1.793 0 0.189   0            0
      85   neutral right      3619  left  -1116.143  1 1.793 0 0.189   0            0
      86   neutral right      3619 right  -1116.143  1 1.793 0 0.189   0            0
      87   neutral right      5276  left  -1627.599  1 1.793 0 0.189   0            0
      88   neutral right      5276 right  -1627.599  1 1.793 0 0.189   0            0
      89   neutral right      6853  left  -2114.361  1 1.793 0 0.189   0            0
      90   neutral right      6853 right  -2114.361  1 1.793 0 0.189   0            0
      91   neutral right      8684  left  -2679.524  1 1.793 0 0.189   0            0
      92   neutral right      8684 right  -2679.524  1 1.793 0 0.189   0            0
      93   neutral right     10391  left  -3206.413  1 1.793 0 0.189   0            0
      94   neutral right     10391 right  -3206.413  1 1.793 0 0.189   0            0
      95   neutral right     12206  left  -3766.637  1 1.793 0 0.189   0            0
      96   neutral right     12206 right  -3766.637  1 1.793 0 0.189   0            0
      97   neutral right     14038  left  -4332.109  1 1.793 0 0.189   0            0
      98   neutral right     14038 right  -4332.109  1 1.793 0 0.189   0            0
      99   neutral right     15811  left  -4879.369  1 1.793 0 0.189   0            0
      100  neutral right     15811 right  -4879.369  1 1.793 0 0.189   0            0
      101 accuracy right         1  left     -2.289  1 3.051 0 0.189   0            0
      102 accuracy right         1 right     -2.289  1 3.051 0 0.189   0            0
      103 accuracy right      1706  left  -4089.073  1 3.051 0 0.189   0            0
      104 accuracy right      1706 right  -4089.073  1 3.051 0 0.189   0            0
      105 accuracy right      3611  left  -8655.244  1 3.051 0 0.189   0            0
      106 accuracy right      3611 right  -8655.244  1 3.051 0 0.189   0            0
      107 accuracy right      5429  left -13012.881  1 3.051 0 0.189   0            0
      108 accuracy right      5429 right -13012.881  1 3.051 0 0.189   0            0
      109 accuracy right      7224  left -17315.389  1 3.051 0 0.189   0            0
      110 accuracy right      7224 right -17315.389  1 3.051 0 0.189   0            0
      111 accuracy right      8890  left -21308.691  1 3.051 0 0.189   0            0
      112 accuracy right      8890 right -21308.691  1 3.051 0 0.189   0            0
      113 accuracy right     10666  left -25565.657  1 3.051 0 0.189   0            0
      114 accuracy right     10666 right -25565.657  1 3.051 0 0.189   0            0
      115 accuracy right     12404  left -29731.539  1 3.051 0 0.189   0            0
      116 accuracy right     12404 right -29731.539  1 3.051 0 0.189   0            0
      117 accuracy right     14094  left -33782.368  1 3.051 0 0.189   0            0
      118 accuracy right     14094 right -33782.368  1 3.051 0 0.189   0            0
      119 accuracy right     15810  left -37895.518  1 3.051 0 0.189   0            0
      120 accuracy right     15810 right -37895.518  1 3.051 0 0.189   0            0
          pGuess     b
      1        0 3.048
      2        0 3.048
      3        0 3.048
      4        0 3.048
      5        0 3.048
      6        0 3.048
      7        0 3.048
      8        0 3.048
      9        0 3.048
      10       0 3.048
      11       0 3.048
      12       0 3.048
      13       0 3.048
      14       0 3.048
      15       0 3.048
      16       0 3.048
      17       0 3.048
      18       0 3.048
      19       0 3.048
      20       0 3.048
      21       0 1.793
      22       0 1.793
      23       0 1.793
      24       0 1.793
      25       0 1.793
      26       0 1.793
      27       0 1.793
      28       0 1.793
      29       0 1.793
      30       0 1.793
      31       0 1.793
      32       0 1.793
      33       0 1.793
      34       0 1.793
      35       0 1.793
      36       0 1.793
      37       0 1.793
      38       0 1.793
      39       0 1.793
      40       0 1.793
      41       0 3.051
      42       0 3.051
      43       0 3.051
      44       0 3.051
      45       0 3.051
      46       0 3.051
      47       0 3.051
      48       0 3.051
      49       0 3.051
      50       0 3.051
      51       0 3.051
      52       0 3.051
      53       0 3.051
      54       0 3.051
      55       0 3.051
      56       0 3.051
      57       0 3.051
      58       0 3.051
      59       0 3.051
      60       0 3.051
      61       0 3.048
      62       0 3.048
      63       0 3.048
      64       0 3.048
      65       0 3.048
      66       0 3.048
      67       0 3.048
      68       0 3.048
      69       0 3.048
      70       0 3.048
      71       0 3.048
      72       0 3.048
      73       0 3.048
      74       0 3.048
      75       0 3.048
      76       0 3.048
      77       0 3.048
      78       0 3.048
      79       0 3.048
      80       0 3.048
      81       0 1.793
      82       0 1.793
      83       0 1.793
      84       0 1.793
      85       0 1.793
      86       0 1.793
      87       0 1.793
      88       0 1.793
      89       0 1.793
      90       0 1.793
      91       0 1.793
      92       0 1.793
      93       0 1.793
      94       0 1.793
      95       0 1.793
      96       0 1.793
      97       0 1.793
      98       0 1.793
      99       0 1.793
      100      0 1.793
      101      0 3.051
      102      0 3.051
      103      0 3.051
      104      0 3.051
      105      0 3.051
      106      0 3.051
      107      0 3.051
      108      0 3.051
      109      0 3.051
      110      0 3.051
      111      0 3.051
      112      0 3.051
      113      0 3.051
      114      0 3.051
      115      0 3.051
      116      0 3.051
      117      0 3.051
      118      0 3.051
      119      0 3.051
      120      0 3.051

---

    Code
      mapped_pars(prior(des, mu_mean = c(v_covariate = 1)))
    Output
                 E     S covariate    lR     v sv B A t0 eta pContaminant pGuess b
      1      speed  left         3  left     3  1 1 0  1   0            0      0 1
      2      speed  left         3 right     3  1 1 0  1   0            0      0 1
      3      speed  left      1853  left  1853  1 1 0  1   0            0      0 1
      4      speed  left      1853 right  1853  1 1 0  1   0            0      0 1
      5      speed  left      3463  left  3463  1 1 0  1   0            0      0 1
      6      speed  left      3463 right  3463  1 1 0  1   0            0      0 1
      7      speed  left      5353  left  5353  1 1 0  1   0            0      0 1
      8      speed  left      5353 right  5353  1 1 0  1   0            0      0 1
      9      speed  left      7153  left  7153  1 1 0  1   0            0      0 1
      10     speed  left      7153 right  7153  1 1 0  1   0            0      0 1
      11     speed  left      9061  left  9061  1 1 0  1   0            0      0 1
      12     speed  left      9061 right  9061  1 1 0  1   0            0      0 1
      13     speed  left     10735  left 10735  1 1 0  1   0            0      0 1
      14     speed  left     10735 right 10735  1 1 0  1   0            0      0 1
      15     speed  left     12327  left 12327  1 1 0  1   0            0      0 1
      16     speed  left     12327 right 12327  1 1 0  1   0            0      0 1
      17     speed  left     14033  left 14033  1 1 0  1   0            0      0 1
      18     speed  left     14033 right 14033  1 1 0  1   0            0      0 1
      19     speed  left     15809  left 15809  1 1 0  1   0            0      0 1
      20     speed  left     15809 right 15809  1 1 0  1   0            0      0 1
      21   neutral  left         6  left     6  1 1 0  1   0            0      0 1
      22   neutral  left         6 right     6  1 1 0  1   0            0      0 1
      23   neutral  left      1680  left  1680  1 1 0  1   0            0      0 1
      24   neutral  left      1680 right  1680  1 1 0  1   0            0      0 1
      25   neutral  left      3409  left  3409  1 1 0  1   0            0      0 1
      26   neutral  left      3409 right  3409  1 1 0  1   0            0      0 1
      27   neutral  left      5071  left  5071  1 1 0  1   0            0      0 1
      28   neutral  left      5071 right  5071  1 1 0  1   0            0      0 1
      29   neutral  left      6804  left  6804  1 1 0  1   0            0      0 1
      30   neutral  left      6804 right  6804  1 1 0  1   0            0      0 1
      31   neutral  left      8739  left  8739  1 1 0  1   0            0      0 1
      32   neutral  left      8739 right  8739  1 1 0  1   0            0      0 1
      33   neutral  left     10461  left 10461  1 1 0  1   0            0      0 1
      34   neutral  left     10461 right 10461  1 1 0  1   0            0      0 1
      35   neutral  left     12287  left 12287  1 1 0  1   0            0      0 1
      36   neutral  left     12287 right 12287  1 1 0  1   0            0      0 1
      37   neutral  left     14137  left 14137  1 1 0  1   0            0      0 1
      38   neutral  left     14137 right 14137  1 1 0  1   0            0      0 1
      39   neutral  left     15816  left 15816  1 1 0  1   0            0      0 1
      40   neutral  left     15816 right 15816  1 1 0  1   0            0      0 1
      41  accuracy  left         5  left     5  1 1 0  1   0            0      0 1
      42  accuracy  left         5 right     5  1 1 0  1   0            0      0 1
      43  accuracy  left      1764  left  1764  1 1 0  1   0            0      0 1
      44  accuracy  left      1764 right  1764  1 1 0  1   0            0      0 1
      45  accuracy  left      3543  left  3543  1 1 0  1   0            0      0 1
      46  accuracy  left      3543 right  3543  1 1 0  1   0            0      0 1
      47  accuracy  left      5235  left  5235  1 1 0  1   0            0      0 1
      48  accuracy  left      5235 right  5235  1 1 0  1   0            0      0 1
      49  accuracy  left      7069  left  7069  1 1 0  1   0            0      0 1
      50  accuracy  left      7069 right  7069  1 1 0  1   0            0      0 1
      51  accuracy  left      8725  left  8725  1 1 0  1   0            0      0 1
      52  accuracy  left      8725 right  8725  1 1 0  1   0            0      0 1
      53  accuracy  left     10462  left 10462  1 1 0  1   0            0      0 1
      54  accuracy  left     10462 right 10462  1 1 0  1   0            0      0 1
      55  accuracy  left     12289  left 12289  1 1 0  1   0            0      0 1
      56  accuracy  left     12289 right 12289  1 1 0  1   0            0      0 1
      57  accuracy  left     14136  left 14136  1 1 0  1   0            0      0 1
      58  accuracy  left     14136 right 14136  1 1 0  1   0            0      0 1
      59  accuracy  left     15817  left 15817  1 1 0  1   0            0      0 1
      60  accuracy  left     15817 right 15817  1 1 0  1   0            0      0 1
      61     speed right         2  left     2  1 1 0  1   0            0      0 1
      62     speed right         2 right     2  1 1 0  1   0            0      0 1
      63     speed right      1812  left  1812  1 1 0  1   0            0      0 1
      64     speed right      1812 right  1812  1 1 0  1   0            0      0 1
      65     speed right      3494  left  3494  1 1 0  1   0            0      0 1
      66     speed right      3494 right  3494  1 1 0  1   0            0      0 1
      67     speed right      5243  left  5243  1 1 0  1   0            0      0 1
      68     speed right      5243 right  5243  1 1 0  1   0            0      0 1
      69     speed right      7089  left  7089  1 1 0  1   0            0      0 1
      70     speed right      7089 right  7089  1 1 0  1   0            0      0 1
      71     speed right      8655  left  8655  1 1 0  1   0            0      0 1
      72     speed right      8655 right  8655  1 1 0  1   0            0      0 1
      73     speed right     10545  left 10545  1 1 0  1   0            0      0 1
      74     speed right     10545 right 10545  1 1 0  1   0            0      0 1
      75     speed right     12274  left 12274  1 1 0  1   0            0      0 1
      76     speed right     12274 right 12274  1 1 0  1   0            0      0 1
      77     speed right     13900  left 13900  1 1 0  1   0            0      0 1
      78     speed right     13900 right 13900  1 1 0  1   0            0      0 1
      79     speed right     15818  left 15818  1 1 0  1   0            0      0 1
      80     speed right     15818 right 15818  1 1 0  1   0            0      0 1
      81   neutral right         8  left     8  1 1 0  1   0            0      0 1
      82   neutral right         8 right     8  1 1 0  1   0            0      0 1
      83   neutral right      1733  left  1733  1 1 0  1   0            0      0 1
      84   neutral right      1733 right  1733  1 1 0  1   0            0      0 1
      85   neutral right      3619  left  3619  1 1 0  1   0            0      0 1
      86   neutral right      3619 right  3619  1 1 0  1   0            0      0 1
      87   neutral right      5276  left  5276  1 1 0  1   0            0      0 1
      88   neutral right      5276 right  5276  1 1 0  1   0            0      0 1
      89   neutral right      6853  left  6853  1 1 0  1   0            0      0 1
      90   neutral right      6853 right  6853  1 1 0  1   0            0      0 1
      91   neutral right      8684  left  8684  1 1 0  1   0            0      0 1
      92   neutral right      8684 right  8684  1 1 0  1   0            0      0 1
      93   neutral right     10391  left 10391  1 1 0  1   0            0      0 1
      94   neutral right     10391 right 10391  1 1 0  1   0            0      0 1
      95   neutral right     12206  left 12206  1 1 0  1   0            0      0 1
      96   neutral right     12206 right 12206  1 1 0  1   0            0      0 1
      97   neutral right     14038  left 14038  1 1 0  1   0            0      0 1
      98   neutral right     14038 right 14038  1 1 0  1   0            0      0 1
      99   neutral right     15811  left 15811  1 1 0  1   0            0      0 1
      100  neutral right     15811 right 15811  1 1 0  1   0            0      0 1
      101 accuracy right         1  left     1  1 1 0  1   0            0      0 1
      102 accuracy right         1 right     1  1 1 0  1   0            0      0 1
      103 accuracy right      1706  left  1706  1 1 0  1   0            0      0 1
      104 accuracy right      1706 right  1706  1 1 0  1   0            0      0 1
      105 accuracy right      3611  left  3611  1 1 0  1   0            0      0 1
      106 accuracy right      3611 right  3611  1 1 0  1   0            0      0 1
      107 accuracy right      5429  left  5429  1 1 0  1   0            0      0 1
      108 accuracy right      5429 right  5429  1 1 0  1   0            0      0 1
      109 accuracy right      7224  left  7224  1 1 0  1   0            0      0 1
      110 accuracy right      7224 right  7224  1 1 0  1   0            0      0 1
      111 accuracy right      8890  left  8890  1 1 0  1   0            0      0 1
      112 accuracy right      8890 right  8890  1 1 0  1   0            0      0 1
      113 accuracy right     10666  left 10666  1 1 0  1   0            0      0 1
      114 accuracy right     10666 right 10666  1 1 0  1   0            0      0 1
      115 accuracy right     12404  left 12404  1 1 0  1   0            0      0 1
      116 accuracy right     12404 right 12404  1 1 0  1   0            0      0 1
      117 accuracy right     14094  left 14094  1 1 0  1   0            0      0 1
      118 accuracy right     14094 right 14094  1 1 0  1   0            0      0 1
      119 accuracy right     15810  left 15810  1 1 0  1   0            0      0 1
      120 accuracy right     15810 right 15810  1 1 0  1   0            0      0 1

---

    Code
      mapped_pars(samples_LNR)
    Output
                E     S    lR    lM      m     s    t0
      1     speed  left  left  TRUE -1.225 0.587 0.197
      2     speed  left right FALSE -0.711 0.587 0.197
      3   neutral  left  left  TRUE -1.225 0.587 0.197
      4   neutral  left right FALSE -0.711 0.587 0.197
      5  accuracy  left  left  TRUE -1.225 0.587 0.197
      6  accuracy  left right FALSE -0.711 0.587 0.197
      7     speed right  left FALSE -0.711 0.587 0.197
      8     speed right right  TRUE -1.225 0.587 0.197
      9   neutral right  left FALSE -0.711 0.587 0.197
      10  neutral right right  TRUE -1.225 0.587 0.197
      11 accuracy right  left FALSE -0.711 0.587 0.197
      12 accuracy right right  TRUE -1.225 0.587 0.197

---

    Code
      mapped_pars(get_prior(samples_LNR))
    Output
                E     S    lR    lM m s t0
      1     speed  left  left  TRUE 0 1  1
      2     speed  left right FALSE 0 1  1
      3   neutral  left  left  TRUE 0 1  1
      4   neutral  left right FALSE 0 1  1
      5  accuracy  left  left  TRUE 0 1  1
      6  accuracy  left right FALSE 0 1  1
      7     speed right  left FALSE 0 1  1
      8     speed right right  TRUE 0 1  1
      9   neutral right  left FALSE 0 1  1
      10  neutral right right  TRUE 0 1  1
      11 accuracy right  left FALSE 0 1  1
      12 accuracy right right  TRUE 0 1  1

---

    Code
      mapped_pars(get_design(samples_LNR))
    Output
      $m 
        lM     
       TRUE   : m + 0.5 * m_lMd
       FALSE  : m - 0.5 * m_lMd
      

# map

    Code
      credint(samples_LNR, selection = "mu", map = "E")
    Output
      $mu
                     2.5%    50%  97.5%
      m_Eaccuracy  -1.008 -0.964 -0.919
      m_Eneutral   -1.008 -0.964 -0.919
      m_Espeed     -1.008 -0.964 -0.919
      s_Eaccuracy   0.554  0.583  0.619
      s_Eneutral    0.554  0.583  0.619
      s_Espeed      0.554  0.583  0.619
      t0_Eaccuracy  0.177  0.193  0.206
      t0_Eneutral   0.177  0.193  0.206
      t0_Espeed     0.177  0.193  0.206
      

---

    Code
      credint(samples_LNR, selection = "mu", map = list(~ E * S))
    Output
      $mu
                            2.5%    50%  97.5%
      m_(Intercept)       -1.008 -0.964 -0.919
      m_Eneutral          -1.008 -0.964 -0.919
      m_Eaccuracy         -1.008 -0.964 -0.919
      m_Sright            -1.008 -0.964 -0.919
      m_Eneutral:Sright   -1.008 -0.964 -0.919
      m_Eaccuracy:Sright  -1.008 -0.964 -0.919
      s_(Intercept)        0.554  0.583  0.619
      s_Eneutral           0.554  0.583  0.619
      s_Eaccuracy          0.554  0.583  0.619
      s_Sright             0.554  0.583  0.619
      s_Eneutral:Sright    0.554  0.583  0.619
      s_Eaccuracy:Sright   0.554  0.583  0.619
      t0_(Intercept)       0.177  0.193  0.206
      t0_Eneutral          0.177  0.193  0.206
      t0_Eaccuracy         0.177  0.193  0.206
      t0_Sright            0.177  0.193  0.206
      t0_Eneutral:Sright   0.177  0.193  0.206
      t0_Eaccuracy:Sright  0.177  0.193  0.206
      

---

    Code
      credint(samples_LNR, selection = "mu", map = TRUE)
    Output
      $mu
                  2.5%    50%  97.5%
      m_lMFALSE -0.752 -0.700 -0.655
      m_lMTRUE  -1.271 -1.228 -1.167
      s          0.554  0.583  0.619
      t0         0.177  0.193  0.206
      

