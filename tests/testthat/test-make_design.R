design_data <- design(data = data.frame(forstmann, CO = 1:nrow(forstmann)),
            model=LBA,matchfun=function(d)d$S==d$lR,
            formula=list(v~lM,sv~lM,B~E+lR,t0~E2 + CO),
            contrasts=list(v = list(lM=matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d")))),
            constants=c(sv=log(1)),
            functions = list(
              E2 = function(d) factor(d$E!="speed",labels=c("speed","nonspeed"))
            ))


design_custom <- design(factors = list(S = c("left", "right"),
                                         subjects = 1:3),
                          Rlevels = c("left", "right"), model = LNR,
                          formula =list(m~0+S,s~1, t0~1),
                          constants=c(s=log(1)))

test_that("design", {
  expect_snapshot(
    str(design_data, give.attr = FALSE)
  )
  expect_snapshot(
    str(design_custom, give.attr = FALSE)
  )
})

test_that("auto covariate detection ignores unused numeric columns", {
  dat <- data.frame(forstmann, CO = seq_len(nrow(forstmann)),
                    UNUSED_NUM = seq_len(nrow(forstmann)))
  des <- design(
    data = dat,
    model = LBA,
    matchfun = function(d) d$S == d$lR,
    formula = list(v ~ lM, sv ~ lM, B ~ E + lR, t0 ~ E2 + CO),
    contrasts = list(v = list(lM = matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d")))),
    constants = c(sv = log(1)),
    functions = list(
      E2 = function(d) factor(d$E != "speed", labels = c("speed", "nonspeed"))
    )
  )

  expect_true("CO" %in% des$Fcovariates)
  expect_false("UNUSED_NUM" %in% des$Fcovariates)
})
