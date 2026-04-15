library(EMC2)
library(microbenchmark)

set.seed(1)
n <- 10000
t <- runif(n, 0.3, 2)
v <- runif(n, 0.5, 3)
B <- runif(n, 0.1, 1)
A <- runif(n, 0.1, 1)
t0 <- rep(0.2, n)

cat("Benchmarking dWald with n =", n, "\n")
res_d <- microbenchmark(
  out_d <- EMC2:::dWald(t + 0, v, B, A, t0),
  times = 100
)
print(res_d)
cat("Mean dWald output:", mean(out_d), "\n")

cat("\nBenchmarking pWald with n =", n, "\n")
res_p <- microbenchmark(
  out_p <- EMC2:::pWald(t + 0, v, B, A, t0),
  times = 100
)
print(res_p)
cat("Mean pWald output:", mean(out_p), "\n")
