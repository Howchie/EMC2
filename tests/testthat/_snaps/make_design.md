# design

    Code
      str(design_data, give.attr = FALSE)
    Output
      List of 13
       $ Flist      :List of 6
        ..$ :Class 'formula'  language v ~ lM
        ..$ :Class 'formula'  language sv ~ lM
        ..$ :Class 'formula'  language B ~ E + lR
        ..$ :Class 'formula'  language t0 ~ E2 + CO
        ..$ :Class 'formula'  language A ~ 1
        ..$ :Class 'formula'  language pContaminant ~ 1
       $ Ffactors   :List of 3
        ..$ subjects: chr [1:19] "as1t" "bd6t" "bl1t" "hsft" ...
        ..$ E       : chr [1:3] "speed" "neutral" "accuracy"
        ..$ S       : chr [1:2] "left" "right"
       $ Rlevels    : chr [1:2] "left" "right"
       $ Clist      :List of 1
        ..$ v:List of 1
        .. ..$ lM: num [1:2, 1] -0.5 0.5
       $ matchfun   :function (d)  
       $ constants  : Named num [1:3] 0 -Inf -Inf
       $ Fcovariates: chr "CO"
       $ Ffunctions :List of 1
        ..$ E2:function (d)  
       $ model      :function ()  
       $ LT         : num 0
       $ LC         : num 0
       $ UC         : num Inf
       $ UT         : num Inf

---

    Code
      str(design_custom, give.attr = FALSE)
    Output
      List of 13
       $ Flist      :List of 4
        ..$ :Class 'formula'  language m ~ 0 + S
        ..$ :Class 'formula'  language s ~ 1
        ..$ :Class 'formula'  language t0 ~ 1
        ..$ :Class 'formula'  language pContaminant ~ 1
       $ Ffactors   :List of 2
        ..$ S       : Factor w/ 2 levels "left","right": 1 2
        ..$ subjects: Factor w/ 3 levels "1","2","3": 1 2 3
       $ Rlevels    : chr [1:2] "left" "right"
       $ Clist      : NULL
       $ matchfun   : NULL
       $ constants  : Named num [1:2] 0 -Inf
       $ Fcovariates: NULL
       $ Ffunctions : NULL
       $ model      :function ()  
       $ LT         : NULL
       $ LC         : NULL
       $ UC         : NULL
       $ UT         : NULL

