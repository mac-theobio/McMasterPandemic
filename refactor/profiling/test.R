library(profvis)
library(McMasterPandemic)
library(anytime)
library(ggplot2)
theme_set(theme_bw())
library(cowplot)
library(tidyverse)
library(zoo)
# Required to make profvis tell us which problematic lines exist in the original run_sim function


source("../../R/sim_funs.R")
source("../../R/calibrate.R")
source("../../R/afuns.R")
source("../../R/utils.R")
source("../../R/testify.R")
source("../../R/methods.R")
source("../../R/kfuns.R")



library(bbmle)

params1 <- read_params("ICU1.csv")
state1 <- make_state(params = params1)
sdate <- "2020-Feb-10"
edate <- "2020-Jun-1"

#################################################################################################
# Non-stochastic simulation

state1 <- state1 
#state1['S'] <- 0
params1['N'] <- sum(state1)
res1 <- run_sim(params = params1, state = state1, start_date = sdate, end_date = edate)

