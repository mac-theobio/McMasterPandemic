library(microbenchmark)

M <- matrix(1:81, 9, 9)

env <- list2env(list(v=M))

by_value <- function(v) {
    v[5] <- 1000
    return(v)
}
    
by_reference <- function(e) {
    e$v[5] <- 1000
}


microbenchmark({by_value(M)})
microbenchmark({by_reference(env)})
rm(list=ls())