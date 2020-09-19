library(McMasterPandemic)
library(tidyverse)

# callArgs <- "spline_recalib.Rout spline_recalib.R cachestuff/spline_calib.rda"


source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

ndf <- 3

print(ff)

cc <- coef(ff,"fitted")

opt_pars <- list(params=c(log_beta0=log(as.numeric(cc$params[1]))
#	, log_E0=log(as.numeric(cc$params[2]))
		)
# , logit_c_prop = plogis(c_prop)
# , logit_phi1 = plogis(phi1)
	, log_nb_disp = log(cc$nb_disp)
#	, time_beta = cc$time_beta
)

params <- ff$forecast_args$base_params

# debug(calibrate_comb)
# debug(calibrate)
# debug(run_sim_loglin)
# x11()

dd_resim <-(predict(ff,ensembles=FALSE) 
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

saveVars(ff_refit,dd_resim)

