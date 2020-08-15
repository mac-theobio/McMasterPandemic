library(McMasterPandemic)
library(tidyverse)

source("makestuff/makeRfuns.R")
commandEnvironments()

## functions: update parameters, simulate, calibrate

update_pars <- function(pp){
	pars[["testing_intensity"]] <- pp[[1]]
	return(pars)
}

simulate_testify_sim <- function(p){
	sims <- (run_sim(params = p
		, ratemat_args = list(testify=TRUE)
		, start_date = start
		, end_date = end
		, use_ode = use_ode
		, step_args = list(testwt_scale = testwt_scale)
		, condense_args = list(keep_all = keep_all
			, add_reports = !keep_all
		)
	)
	%>% mutate(testing_intensity=p[["testing_intensity"]])
	)
	return(sims)
}

calibrate_sim <- function(dd, pars, p){
	dat <- (dd
		%>% transmute(date
			, postest
			, death
			, H
			)
		%>% gather(key="var",value="value",-date)
		%>% mutate(value=round(value))
	)
	dat2 <- dat %>% rowwise() %>% filter(grepl(var,p[[2]]))
	opt_pars <- with(as.list(pars)
		, list(params=c(log_beta0 = log(beta0)
			, log_E0 = log(E0)
			)
		)
	)
	if(p[[3]]){
		opt_pars <- with(as.list(pars)
			, list(params=c(log_beta0 = log(beta0)
				, log_E0 = log(E0)
				, log_testing_intensity = log(testing_intensity)
				)
			)
		)
	}
	sim_args <- list(ratemat_args = list(testify=TRUE, testing_time="report"))
	mod <- do.call(calibrate_comb
		, c(nlist(params = pars
			, use_DEoptim = FALSE
			, data = dat2
			, opt_pars = opt_pars
			, sim_args = sim_args
			, maxit = 1000
			)
		)
	)
	return(mod)
}

saveVars(update_pars, simulate_testify_sim, calibrate_sim)
