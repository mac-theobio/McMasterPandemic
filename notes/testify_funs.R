library(McMasterPandemic)
library(tidyverse)

source("makestuff/makeRfuns.R")
commandEnvironments()

##' @param p a set of parameters
##' @param testing_data a data frame with dates, per capita testing intensity per day
simtestify <- function(p, testing_data){
   testing_data <- with(testing_data
      , data.frame(Date=Date
         , Symbol="testing_intensity"
         , Relative_value=c(1, intensity[-1]/intensity[1])
      )
   )
   
	sim_args <- list(ratemat_args = list(testing_time=testing_time)
	   , start_date = start
           , end_date = end
           , use_ode = use_ode
           , stoch=c(obs=stoch_obs, proc=FALSE)  
           , step_args = list(testwt_scale = "sum_smooth")
           , condense_args = list(keep_all = keep_all
                                , add_reports = !keep_all
                                  )
	        , params_timevar = testing_data
	)
   sims <- do.call(run_sim,c(list(params=p),sim_args))
	return(sims)
}

calibrate_sim <- function(dd, pars, p,testing_data,debug_plot=FALSE,
                          debug=FALSE, debug_hist=FALSE){
    ## change sim output to input format
    dat <- (dd %>% select(date
			, postest
			, death
			, H
                          )
		%>% gather(key="var",value="value",-date)
		%>% mutate(value=round(value))
    )
    dat2 <- dat %>% rowwise() %>% filter(grepl(var,p$keep_vars))
    opt_pars <- with(as.list(pars)
                   , list(params=c(log_beta0 = log(beta0)
                                 , log_E0 = log(E0)
                                   )
                          )
                     )
    if(p$opt_testify){
        opt_pars <- c(opt_pars,
                      list(log_testing_intensity = log(pars[["testing_intensity"]])))
    }
    sim_args <- list(ratemat_args = list(testing_time=testing_time))
    mod <- do.call(calibrate_comb
		, c(nlist(params = pars
			, use_DEoptim = FALSE
			, use_spline = FALSE
			, debug_plot = debug_plot
                        , debug_hist = debug_hist
                        , debug = debug
			, data = dat2
			, opt_pars = opt_pars
			, sim_args = sim_args
			, use_testing = TRUE
			, testing_data = testing_data
			, start_date_offset=0
			, maxit = 1000
			)
		)
	)
	return(mod)
}

saveVars(simtestify, calibrate_sim)
