library(McMasterPandemic)

source("makestuff/makeRfuns.R")
commandEnvironments()

Ontario_fit_spline <- do.call(calibrate_comb
, c(nlist(params=params
	, debug_plot=FALSE
     , data=calibrate_data_fill
     , mob_data = clean_mobility
     , opt_pars = opt_pars
     , use_DEoptim = TRUE
	, DE_cores = 7
	, use_phenomhet = TRUE
	, use_mobility = TRUE
	, mob_breaks = "2020-04-15"
	, mob_breaks_int = TRUE
	, mob_logist_scale = 3
	, use_spline = TRUE
	, spline_df = NA
	, spline_days = 14
)
)
)

#save.image(list=c("Ontario_fit", "calibrate_data_fill", "clean_mobility") file = "Ontario_basic.rda")


saveEnvironment()


