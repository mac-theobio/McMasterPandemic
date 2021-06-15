library(testthat)
library(microbenchmark)
library(Matrix)

sweep_sparse <- function(x, margin, stats, fun = "*") {
  f <- match.fun(fun)
  if (margin == 1) {
    idx <- x@i + 1
  } else {
    idx <- x@j + 1
  }
  x@x <- f(x@x, stats[idx])
  return(x)
}

generate_dense_matrix <- function(n){
  names <- c(99:(99 + n - 1))

  M <- matrix(runif(n * n), nrow = n, ncol = n)
  rownames(M) <- names
  colnames(M) <- names
  v <- runif(n)
  names(v) <- names
  return(list(M=M, v=v))
}

generate_sparse_matrix <- function(n, sparsity=0.1){
  names <- c(99:(99 + n - 1))

  M <- rsparsematrix(nrow=n, ncol=n, density = sparsity, rand.x = runif)
  rownames(M) <- names
  colnames(M) <- names
  v <- runif(n)
  names(v) <- names
  return(list(M=M, v=v))
}

# -----------------------------------------------------------------------
L <- generate_dense_matrix(20)
M <- L$M
v <- L$v

baseline <- sweep(M, v, MARGIN = 1, FUN = "*")
other_0 <- t(t(M) * rep(v, rep.int(nrow(M), length(v))))
other_1 = M * v

stopifnot(baseline == other_0)
stopifnot(baseline == other_1)

microbenchmark(sweep(M, v, MARGIN = 1, FUN = "*"))
microbenchmark(t(t(M) * rep(v, rep.int(nrow(M), length(v)))))
microbenchmark( M * v)
# -----------------------------------------------------------------------
# quick checks to ensure equivalence of these methods
L <- generate_sparse_matrix(20)
M <- L$M
v <- L$v

baseline <- sweep(M, v, MARGIN = 1, FUN = "*")
other_0 <- M * v
other_1 <- sweep_sparse(M, margin=1, v, fun = "*")


stopifnot(all(baseline == other_0))
stopifnot(all(baseline == other_1))
microbenchmark(sweep_sparse(M, margin=1, v, fun = "*"))
microbenchmark(sweep(M, v, MARGIN = 1, FUN = "*"))
microbenchmark(M * v)






