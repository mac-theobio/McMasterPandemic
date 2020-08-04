# Goal: to setup 10 simulations and fit each with natural spline, 
# using log_beta0 transformation, no stochasticity, (whatelse?)

library(McMasterPandemic)
library(tidyverse)
library(anytime)
library(bbmle)
library(parallel)

source("makestuff/makeRfuns.R")

# #######################################################
# run_params
# options(mc.cores=8)

## setup 
params <- fix_pars(read_params("../inst/params/PHAC.csv"), target=c(R0=3, Gbar=6))
summary(params) 
print(params)

# set_initial_state
state1 <- make_state(params=params)

# start and end dates
start_date <- anydate("2020-01-01")
end_date <- anydate("2020-03-31") 

res1 <- run_sim(params=params, state=state1, start_date=start_date, end_date=end_date)
summary(res1)

# saveVars(res1) #doesn't work for now?
plot(res1,log=TRUE)

# look at ./rshape.R and find out the description of Rt that fits here. 
R0 <- summary(res1)$R0
(res1$R)
res1$date


# 







