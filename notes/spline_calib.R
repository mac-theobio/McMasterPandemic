library(McMasterPandemic)
library(tidyverse)
library(zoo)
library(parallel)
library(furrr)
library(future.batchtools)
library(splines)

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

objects()

start_date <- first_date
end_date <- start_date -1 + fitmax
obs_disp <- 50
proc_disp <- 1

params <- read_params(matchFile(".csv$"))

print(bb)

params <- fix_pars(params, target=c(R0=Rt[[1]]))

params["obs_disp"] <- obs_disp
params["proc_disp"] <- proc_disp
# params["obs_disp_report"] <- obs_disp
# params["obs_disp_death"] <- obs_disp

opt_pars <- list(params=c(log_beta0 = as.numeric(log(params["beta0"]))
								  # , log_E0=log(as.numeric(params["E0"]))
)
# , log_nb_disp = c(report=3,death=3)
)

opt_parsE0 <- list(
	params=c(log_beta0 = as.numeric(log(params["beta0"]))
		, log_E0=log(as.numeric(params["E0"]))
	)
	# , log_nb_disp = c(report=3,death=3)
)

sim_calib <- function(x){
	set.seed(x)

ddfull_sim<- (forecast_sim(p = unlist(opt_pars)
	, opt_pars = opt_pars
	, base_params = params
	, stoch = list(obs=TRUE,proc=TRUE)
	, stoch_start = c(proc=min(dd),obs=min(dd))
	, time_args = list(X_date=dd, X=X,extra_pars=list(time_beta=bb))
	, start_date = min(dd)
	, end_date = max(dd)
	, sim_fun = run_sim_loglin
	)
	%>% filter(var %in% c("report"))
#	%>% mutate(value = ifelse(is.na(value),0,value))
)


dd_sim <- (ddfull_sim
	%>% filter(date <= end_date)
	)

ff <- calibrate_comb(params = params
	, debug_plot=FALSE
	, use_DEoptim=TRUE
	, DE_cores = 6
	, opt_pars = opt_pars
	, use_spline = TRUE
	, spline_df = ndf
	, spline_type = "bs"
	, spline_int = FALSE
	, data= dd_sim
	, start_date = min(dd_sim$date)
	, start_date_offset = 0
)

ff2 <- calibrate_comb(params = params
							, debug_plot=FALSE
							, use_DEoptim=TRUE
							, DE_cores = 6
							, opt_pars = opt_pars
							, use_spline = TRUE
							, spline_df = ndf
							, spline_pen = 5
							, spline_type = "bs"
							, spline_int=FALSE
							, data= dd_sim
							, start_date = min(dd_sim$date)
							, start_date_offset = 0
)

ffE0 <- calibrate_comb(params = params
	, debug_plot=FALSE
	, use_DEoptim=TRUE
	, DE_cores = 6
	, opt_pars = opt_parsE0
	, use_spline = TRUE
	, spline_df = ndf
	, spline_type = "bs"
	, spline_int = FALSE
	, data= dd_sim
	, start_date = min(dd_sim$date)
	, start_date_offset = 0
)

ffE02 <- calibrate_comb(params = params
							  , debug_plot=FALSE
							  , use_DEoptim=TRUE
							  , DE_cores = 6
							  , opt_pars = opt_parsE0
							  , use_spline = TRUE
							  , spline_df = ndf
							  , spline_type = "bs"
							  , spline_pen = 5
							  , spline_int = FALSE
							  , data= dd_sim
							  , start_date = min(dd_sim$date)
							  , start_date_offset = 0
)

ff_list <- list(fit=ff, fitE0=ffE0, fitpen=ff2, fitE0pen=ffE02, fitdat=dd_sim, full_dat=ddfull_sim)

saveRDS(object=ff_list, file=paste0("./cachestuff/spline_calib.",x,".RDS"))
}

batch_setup()

future_map(1:10,function(x)sim_calib(x))


