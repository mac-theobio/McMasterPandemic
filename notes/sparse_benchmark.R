library(McMasterPandemic)
library(rbenchmark)
params <- read_params("PHAC_testify.csv")
state <- make_state(params=params)
f <- function(sparse=FALSE,testify=FALSE) {
    run_sim(params,state,start_date="1-March-2020",end_date="1-Jun-2020",
            ratemat_args=list(sparse=sparse,testify=testify))
}
benchmark(f(),f(sparse=TRUE),f(testify=TRUE),f(sparse=TRUE,testify=TRUE))
