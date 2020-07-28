library(McMasterPandemic)

source("../R/testify.R")

p <- read_params("PHAC_testify.csv")
state<- make_state(params=p)


ratemat <- make_ratemat(state,p)

ratemat_testify <- make_ratemat(state,p,testify=TRUE)

beta_vec <- make_betavec(state,p)
beta_vec_testify <- make_betavec(attr(state,"testify"),p)
