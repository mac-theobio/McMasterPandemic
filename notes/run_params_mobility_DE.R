library(McMasterPandemic)
library(tidyverse)
library(anytime)
library(bbmle)
library(parallel)


use_true_start <- TRUE
use_cut <- TRUE
nsim <- 20
options(mc.cores=1)

## setup 

params <- fix_pars(read_params("ICU1.csv"), target=c(R0=3, Gbar=6))
params[["beta0"]] <- 0.9   ## slightly rounded for convenience
params[["obs_disp"]] <- 100 ## BMB: less noise
params[["N"]] <- 1e8
summary(params) 

print(params)

start_date <- anydate("2020-01-01")
end_date <- anydate("2020-03-31") ## BMB: don't run as long

mob_length <- length(start_date:end_date)
mob_vals <- (1:mob_length)^(-0.3)

noresponse <- 15
mob_vals <- c(rep(1,noresponse),mob_vals[1:(mob_length-noresponse-1)])
mob_start <- anydate("2020-01-01")

cut_start <- anydate("2020-01-01")
cut_end <- anydate("2020-03-31")


true_pars <- list(params = c(
	log_E0 = log(params[["E0"]]), 
	log_beta0 = log(params[["beta0"]])
	)
	, logit_mob_power=0.5
	, log_nb_disp = log(params[["obs_disp"]])
)

opt_pars <- list(params=c(
	log_E0 = log(params[["E0"]]*2), 
   log_beta0=log(params[["beta0"]]*1.2)
   )
	, logit_mob_power=0.5
   , log_nb_disp = log(params[["obs_disp"]])
)

true_pars <- unlist(true_pars)


### DEoptim setups

## (lwr/upr not currently used in run_caltest.R?)
## note these must match param set if used ...
lwr <- c(
    params.log_E0 = 1
  , params.log_beta0=-1
  , logit_rel_beta01=-1
  , log_nb_disp = 1
)

upr <- c(
    params.log_E0 = 5
  , params.log_beta0=1
  , logit_rel_beta01 = 4
  , log_nb_disp = 7
)




use_DEoptim <- FALSE
do_plots=FALSE
