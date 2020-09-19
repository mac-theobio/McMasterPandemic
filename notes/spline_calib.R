library(McMasterPandemic)
library(tidyverse)

callArgs <- "spline_calib.Rout spline_calib.R spline_sim.rda spline.csv"

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()


dat <- (sims 
	%>% select(date, report, death)
	%>% gather(key = "var", value="value", -date)
	%>% mutate(value = round(value))
)
opt_pars <- with(as.list(params)
	, list(params=c(log_beta0=log(beta0)
		# , log_E0=log(E0)
		# , logit_c_prop = plogis(c_prop)
		# , logit_phi1 = plogis(phi1)
		)
	, log_nb_disp = c(report=3, death=3)
#	, time_beta = c(bb[-1])
	)
)

# debug(calibrate_comb)
# debug(calibrate)
# debug(run_sim_loglin)
# x11()

ff <- calibrate_comb(params = params
	, debug_plot=FALSE
	, use_DEoptim=FALSE
	, DE_cores = 3
	, opt_pars = opt_pars
	, use_spline = TRUE
	, spline_df = 7
	, spline_type = "ns"
	, data= dat
	, start_date_offset=0
	# , 
)

saveVars(ff, dat)
