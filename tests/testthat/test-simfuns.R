library(testthat)
library(McMasterPandemic)

n = 15
names = c(99:(99+n-1))

M = matrix(runif(n*n), nrow = n, ncol = n)
rownames(M) = names
colnames(M) = names
v = c(1:n)
names(v) = names


test_that('New sweep works correctly', {
baseline = sweep(M, v, MARGIN=1, FUN="*")
expect(identical(col_multiply(M, v), baseline), 'col_multiply not working as expected!')
})

#microbenchmark(baseline = sweep(M, v, MARGIN=1, FUN="*"))
#microbenchmark({new =t(t(M) %*% diag(v));
#n = names(v);
#rownames(new) = n})
#microbenchmark(new= t(t(M) * rep(v, rep.int(nrow(M), length(v)))))
