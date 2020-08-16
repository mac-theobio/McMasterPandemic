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

simtestify <- function(p,testing_data){
   ## create fake data for dates
   dd <- (simulate_testify_sim(p)
      %>% transmute(date
			, postest
			, death
			, H
			)
		%>% gather(key="var",value="value",-date)
		%>% mutate(value=round(value))
	)
   opt_pars <- with(as.list(p)
		, list(params=c(log_beta0 = log(beta0)
			, log_E0 = log(E0)
			)
		)
	)
	sim_args <- list(ratemat_args = list(testify=TRUE)
	   , start_date = start
		, end_date = end
		, use_ode = use_ode
		, step_args = list(testwt_scale = testwt_scale)
		, condense_args = list(keep_all = keep_all
			, add_reports = !keep_all
		)
	)
   time_args <- do.call(calibrate_comb
		, c(nlist(params = p
			, use_DEoptim = FALSE
			, use_spline = FALSE
			, data = dd
			, sim_args = sim_args
			, maxit = 1000
			, return_val = "time_args"
			, start_date = start
                        , end_date = end
			, testing_data = testing_data
			, use_testing = TRUE
			)
		)
	)
   sim_args <- c(sim_args, list(params_timevar = time_args$testing_data))
    sims <- do.call(run_sim,
                    c(list(params=p),
	   # , opt_pars = opt_pars
	   # , base_params = p
                    sim_args)
	   # , start_date = start, end_date = end
                    )   # %>% mutate(testing_intensity=p[["testing_intensity"]])
	# )
	return(sims)
}

calibrate_sim <- function(dd, pars, p,testing_data){
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
			, use_spline = FALSE
			, debug_plot = 
			, data = dat2
			, opt_pars = opt_pars
			, sim_args = sim_args
			, use_testing = TRUE
			, testing_data = testing_data
			, maxit = 1000
			)
		)
	)
	return(mod)
}

saveVars(update_pars, simulate_testify_sim, simtestify, calibrate_sim)
