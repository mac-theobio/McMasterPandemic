library(testthat)
library(McMasterPandemic)

context("very basic simulation")

params <- read_params(system.file("params","ICU1.csv",
                                  package="McMasterPandemic"))
s1 <- run_sim(params,state2,start_date=sdate,end_date="1-Jun-2020")
