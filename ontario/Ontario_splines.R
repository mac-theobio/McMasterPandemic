library(McMasterPandemic)
library(tidyverse)
library(splines)
source("makestuff/makeRfuns.R")
commandFiles()


trim_dat <- calibrate_data_fill %>% filter(date < as.Date("2020-07-01"))

NS_splines <- do.call(calibrate_comb
	, c(nlist(params=params
		, debug_plot=FALSE
#		, sim_args = list(use_ode = TRUE)
     	, data = trim_dat
     	, mob_data = clean_mobility
     	, opt_pars = opt_pars
     	, use_DEoptim = FALSE
		, DE_cores = 7
		, use_phenomhet = FALSE
		, use_mobility = FALSE
		, mob_breaks = "2020-04-15"
		, mob_breaks_int = TRUE
		, mob_logist_scale = 3
		, use_spline = TRUE
		, spline_type = "ns"
		, spline_df = NA
		, spline_days = 14
		)
	)
)

BS_splines <- do.call(calibrate_comb
	, c(nlist(params=params
		, debug_plot=FALSE
#		, sim_args = list(use_ode = TRUE)
     	, data = trim_dat
     	, mob_data = clean_mobility
     	, opt_pars = opt_pars
     	, use_DEoptim = TRUE
		, DE_cores = 7
		, use_phenomhet = FALSE
		, use_mobility = FALSE
		, mob_breaks = "2020-04-15"
		, mob_breaks_int = TRUE
		, mob_logist_scale = 3
		, use_spline = TRUE
		, spline_type = "bs"
		, spline_df = NA
		, spline_days = 14
)
)
)

print(plot(NS_splines, data=calibrate_data_fill) + ggtitle("NS splines"))

print(plot(BS_splines, data=calibrate_data_fill) + ggtitle("BS splines"))

saveEnvironment()


