## run this from McMasterPandemic/sandbox on the master branch

library(tidyverse)
devtools::load_all("..")
options(macpan_pfun_method = "grep")

fit_end_date <- "2021-06-18"

n_sim <- 1
n_cores <- 1

mod <- readRDS("calibration_2021-06-18_ON.RDS")

rr <- forecast_ensemble(mod$fit
  , nsim=n_sim 
  , qvec = c(0.025,0.5,0.975)
  , seed = 1
  , parallel = TRUE
  , n_cores = n_cores
)
