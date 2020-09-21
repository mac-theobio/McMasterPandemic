library(McMasterPandemic)
library(tidyverse)
library(zoo)
library(parallel)
library(furrr)
library(future.batchtools)

# callArgs <- "spline_recalib.Rout spline_recalib.R cachestuff/spline_calib.rda"


source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()


ndf <- 5

cc <- coef(ff,"fitted")

opt_pars <- list(params=c(log_beta0=log(as.numeric(cc$params[1]))
# , log_E0=log(as.numeric(cc$params[2]))
		)
	, log_nb_disp = log(cc$nb_disp)
)

params <- ff$forecast_args$base_params

sim_recalib <- function(x){
	set.seed(x)

dd_resim <-(predict(ff,ensemble=TRUE,nsim=1) 
	%>% filter(var %in% c("report","death"))
	%>% mutate(value = round(value))
)

ff_refit <- calibrate_comb(params = params
	, debug_plot=TRUE
	, use_DEoptim=TRUE
	, DE_cores = 3
	, opt_pars = opt_pars
	, use_spline = TRUE
	, spline_df = ndf
	, spline_type = "ns"
	, data= dd_resim
	, start_date = min(dd_resim$date)
	, start_date_offset = 0
)				

ff_list <- list(fit=ff_refit,data=dd_resim)

saveRDS(object=ff_refit, file=paste0("./cachestuff/spline_recalib.noE0",x,".RDS"))
}

batch_setup()

future_map(1:10,function(x)sim_recalib(x))

