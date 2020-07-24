library(McMasterPandemic)
library(tidyverse)
source("makestuff/makeRfuns.R")
commandFiles()


trim_dat <- calibrate_data_fill %>% filter(date < as.Date("2020-07-01"))

## Don't need?
clean_mobility_cap <- (clean_mobility
	%>% mutate(rel_activity = ifelse(rel_activity > 1, 1, rel_activity))
)

print(clean_mobility)

print(clean_mobility_cap)


Mobility <- do.call(calibrate_comb
	, c(nlist(params=params
		, debug_plot=FALSE
		, sim_args = list(ndt = 4)
     	, data = trim_dat
     	, mob_data = clean_mobility
     	, opt_pars = opt_pars
     	, use_DEoptim = TRUE
		, DE_cores = 7
		, use_phenomhet = FALSE
		, use_mobility = TRUE
		, mob_breaks = "2020-04-15"
		, mob_breaks_int = TRUE
		, mob_logist_scale = 3
		, use_spline = FALSE
		, spline_df = NA
		, spline_days = 14
		)
	)
)

Mobility_ode <- do.call(calibrate_comb
	, c(nlist(params=params
		, debug_plot=FALSE
		, sim_args = list(use_ode = TRUE)
     	, data = trim_dat
     	, mob_data = clean_mobility
     	, opt_pars = opt_pars
     	, use_DEoptim = TRUE
		, DE_cores = 7
		, use_phenomhet = FALSE
		, use_mobility = TRUE
		, mob_breaks = "2020-04-15"
		, mob_breaks_int = TRUE
		, mob_logist_scale = 3
		, use_spline = FALSE
		, spline_df = NA
		, spline_days = 14
)
)
)

print(plot(Mobility, data=calibrate_data_fill) + ggtitle("Mobility"))

print(plot(Mobility_ode, data=calibrate_data_fill) + ggtitle("Mobility ode"))

saveEnvironment()


