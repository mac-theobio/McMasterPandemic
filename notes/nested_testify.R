library(McMasterPandemic)
library(tidyverse)
library(furrr)
library(future.batchtools)

source("makestuff/makeRfuns.R")
commandEnvironments()

fn <- if (interactive()) "PHAC_testify.csv" else matchFile(".csv$")
pars <- (read_params(fn)
    %>% fix_pars(target=c(R0=R0, Gbar=Gbar))
    %>% update(N=pop)
)

## different combinations
testing_intensity <- c(0.001, 0.01, 0.1)
keep_vars <- c("postest", "H/death", "postest/H/death")
opt_testify <- c(TRUE,FALSE)

comboframe <- expand.grid(testing_intensity=testing_intensity
	, keep_vars = keep_vars
	, opt_testify = opt_testify
)

sim_and_calibrate <- function(x){
	pp <- update_pars(x)
	simdat <- simulate_testify_sim(pp)
	#calib_mod <- calibrate_sim(dd=simdat, pars=pp, p=x)
	calib_mod <- NULL
	res_list <- list(fit=calib_mod,params=pp, data=simdat)
	saveRDS(res_list,file=paste0("cachestuff/simcalib.",x,".rds"))
	return(res_list)
}

batch_setup(ncpus=6)

res_list <- future_map(seq(nrow(comboframe)),function(x){sim_and_calibrate(comboframe[x,])})

saveVars(res_list)


