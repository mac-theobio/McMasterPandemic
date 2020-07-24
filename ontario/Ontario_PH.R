library(McMasterPandemic)
library(tidyverse)
source("makestuff/makeRfuns.R")
commandFiles()


trim_dat <- calibrate_data_fill %>% filter(date < as.Date("2020-07-01"))

PH <- do.call(calibrate_comb
, c(nlist(params=params
	, debug_plot=FALSE
     , data=trim_dat
     , mob_data = clean_mobility
     , opt_pars = opt_pars
     , use_DEoptim = TRUE
	, DE_cores = 7
	, use_phenomhet = TRUE
	, use_mobility = FALSE
	, mob_breaks = "2020-04-15"
	, mob_breaks_int = TRUE
	, mob_logist_scale = 3
	, use_spline = FALSE
	, spline_df = NA
	, spline_days = 14
)
)
)

PH_ode <- do.call(calibrate_comb
	, c(nlist(params=params
		, debug_plot=FALSE
		, sim_args = list(use_ode = TRUE)
     	, data=trim_dat
     	, mob_data = clean_mobility
     	, opt_pars = opt_pars
     	, use_DEoptim = TRUE
		, DE_cores = 7
		, use_phenomhet = TRUE
		, use_mobility = FALSE
		, mob_breaks = "2020-04-15"
		, mob_breaks_int = TRUE
		, mob_logist_scale = 3
		, use_spline = FALSE
		, spline_df = NA
		, spline_days = 14
)
)
)

PH_ndt4 <- do.call(calibrate_comb
	, c(nlist(params=params
		, debug_plot=FALSE
		, sim_args = list(ndt = 4)
     	, data=trim_dat
     	, mob_data = clean_mobility
     	, opt_pars = opt_pars
     	, use_DEoptim = TRUE
		, DE_cores = 7
		, use_phenomhet = TRUE
		, use_mobility = FALSE
		, mob_breaks = "2020-04-15"
		, mob_breaks_int = TRUE
		, mob_logist_scale = 3
		, use_spline = FALSE
		, spline_df = NA
		, spline_days = 14
)
)
)

print(plot(PH, data=calibrate_data_fill) + ggtitle("PH"))

print(plot(PH_ode, data=calibrate_data_fill) + ggtitle("PH ode"))

print(plot(PH_ndt4, data=calibrate_data_fill) + ggtitle("PH ndt=4"))

saveEnvironment()


