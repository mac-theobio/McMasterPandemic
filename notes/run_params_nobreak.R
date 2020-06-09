library(McMasterPandemic)
library(tidyverse)
library(anytime)
library(bbmle)

use_true_start <- TRUE
use_cut <- TRUE
nsim <- 49
options(mc.cores=1)

## setup 

params <- fix_pars(read_params("ICU1.csv"), target=c(R0=3, Gbar=6))
params[["E0"]] <- 10
params[["beta0"]] <- 0.9   ## slightly rounded for convenience
params[["obs_disp"]] <- 100 ## BMB: less noise
params[["N"]] <- 1e8
summary(params) 


start_date <- anydate("2020-01-01")
end_date <- anydate("2020-03-31") ## BMB: don't run as long

cut_start <- anydate("2020-01-01")
cut_end <- anydate("2020-03-01")


add_breaks <- FALSE
bd <- NULL

true_pars <- list(
	params=c(
		log_E0 = log(params[["E0"]]), 
		log_beta0=log(params[["beta0"]]))
   , log_nb_disp = log(params[["obs_disp"]])
)

opt_pars <- list(params=c(log_E0= log(params[["E0"]])
		, log_beta0=log(params[["beta0"]])
		)
	, log_nb_disp = c(report=log(params[["obs_disp"]]))
)
	
true_pars <- unlist(true_pars)

