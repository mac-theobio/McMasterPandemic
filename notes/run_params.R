library(McMasterPandemic)
library(tidyverse)
library(anytime)
library(bbmle)
library(parallel)

use_true_start <- TRUE
nsim <- 100
options(mc.cores=1)

## setup 

params <- fix_pars(read_params("ICU1.csv"), target=c(R0=3, Gbar=6))
params[["beta0"]] <- 0.9   ## slightly rounded for convenience
params[["obs_disp"]] <- 100 ## BMB: less noise
params[["N"]] <- 1e8
summary(params) 


start_date <- anydate("2020-01-01")
end_date <- anydate("2020-03-31") ## BMB: don't run as long


breaks <- FALSE
break1 <- anydate("2020-01-22")
rel_break1 <- 0.01
bd <- NULL

## Start with true parameters

true_pars <- list(
    params=c(log_beta0=log(params[["beta0"]]))
	   , log_nb_disp = log(params[["obs_disp"]])
) 

##    optimization breaks if factor of true value is > about 1.5 ... ?
opt_pars <- list(
   params=c(log_beta0=log(params[["beta0"]]*1.2))
		, log_nb_disp = log(params[["obs_disp"]])
)

if(breaks){
	true_pars <- list(params = c(log_beta0 = log(params[["beta0"]]))
		, logit_rel_beta0 = qlogis(rel_break1)
		, log_nb_disp = log(params[["obs_disp"]])
	)

	opt_pars <- list(params=c(log_beta0=log(params[["beta0"]]*1.2))
		, logit_rel_beta0 = qlogis(rel_break1)
		, log_nb_disp = log(params[["obs_disp"]])
		)
	bd <- break1
}

true_pars <- unlist(true_pars)

