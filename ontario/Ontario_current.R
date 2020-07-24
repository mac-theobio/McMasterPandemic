library(McMasterPandemic)
library(tidyverse)
source("makestuff/makeRfuns.R")
commandFiles()


trim_dat <- calibrate_data_fill %>% filter(date < as.Date("2020-07-01"))

clean_mobility_cap <- (clean_mobility
   %>% mutate(rel_activity = ifelse(rel_activity > 1, 1, rel_activity))
)


current <- do.call(calibrate_comb
	, c(nlist(params=params
		, debug_plot=FALSE
		, sim_args = list(ndt = 4)
     	, data = trim_dat
     	, mob_data = clean_mobility_cap
     	, opt_pars = opt_pars
     	, use_DEoptim = TRUE
		, DE_cores = 7
		, use_phenomhet = TRUE
		, use_mobility = TRUE
		, mob_breaks = "2020-04-15"
		, mob_breaks_int = TRUE
		, mob_logist_scale = 3
		, use_spline = FALSE
		, spline_type = "ns"
		, spline_df = NA
		, spline_days = 14
		)
	)
)

current_ode <- do.call(calibrate_comb
	, c(nlist(params=params
		, debug_plot=FALSE
		, sim_args = list(use_ode = TRUE)
     	, data = trim_dat
     	, mob_data = clean_mobility_cap
     	, opt_pars = opt_pars
     	, use_DEoptim = TRUE
		, DE_cores = 7
		, use_phenomhet = TRUE
		, use_mobility = TRUE
		, mob_breaks = "2020-04-15"
		, mob_breaks_int = TRUE
		, mob_logist_scale = 3
		, use_spline = FALSE
		, spline_type = "ns"
		, spline_df = NA
		, spline_days = 14
		)
	)
)

print(plot(current, data=calibrate_data_fill) + ggtitle("Current model: PH + mobility cap"))

print(plot(current_ode, data=calibrate_data_fill) + ggtitle("Current model: PH + mobility cap ode"))


saveEnvironment()


