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
# setwd('\\\\wsl$\\Ubuntu\\home\\ali\\projects\\McMasterPandemic\\inst\\params')

# use_true_start <- TRUE
# use_cut <- TRUE
# nsim <- 1 # for now to check
# options(mc.cores=1)

## setup 
# params <- fix_pars(read_params("PHAC.csv"), target=c(R0=3, Gbar=6))
params <- fix_pars(read_params("../inst/params/PHAC.csv"), target=c(R0=3, Gbar=6))
# params[["beta0"]] <- 0.9   ## slightly rounded for convenience
# params[["obs_disp"]] <- 100 ## BMB: less noise
# params[["N"]] <- 1e8
summary(params) 

print(params)

# set_initial_state
state1 <- make_state(params=params)

# start and end dates
start_date <- anydate("2020-01-01")
end_date <- anydate("2020-03-31") 

res1 <- run_sim(params=params, state=state1, start_date=start_date, end_date=end_date)
summary(res1)

saveVars(res1)

plot(res1,log=TRUE)

echo "end"

