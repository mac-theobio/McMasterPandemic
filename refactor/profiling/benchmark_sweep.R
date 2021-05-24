library(microbenchmark)

n = 20
names = c(99:(99+n-1))

M = matrix(runif(n*n), nrow = n, ncol = n)
rownames(M) = names
colnames(M) = names
v = c(1:n)
names(v) = names


microbenchmark(baseline = sweep(M, v, MARGIN=1, FUN="*"))
microbenchmark({new =t(t(M) %*% diag(v));
               n = names(v);
               rownames(new) = n})
microbenchmark(new= t(t(M) * rep(v, rep.int(nrow(M), length(v)))))

baseline = sweep(M, v, MARGIN=1, FUN="*")
identical(t(t(M) * rep(v, rep.int(nrow(M), length(v)))), baseline)