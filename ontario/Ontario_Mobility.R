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
##### Changing to ICU1 parameters

### MacPan setup
paramsICU1 <- fix_pars(read_params("ICU1.csv")
	, target=c(Gbar=6)
	, pars_adj=list(c("sigma","gamma_s","gamma_m","gamma_a"))
)

paramsICU1[["N"]] <- 14.57e6 ## Population of Ontario (2019)

		### parameters we are trying to estimate
opt_pars <- list(params=c(log_E0=2, log_beta0=-1, logit_c_prop=-1, logit_mu = -1, logit_phi1 = -1),log_nb_disp=c(report=1, death=1,H=1))


Mobility_ICU1 <- do.call(calibrate_comb
	, c(nlist(params=paramsICU1
		, debug_plot=FALSE
		, sim_args = list(ndt = 2)
     	, data = trim_dat
     	, mob_data = clean_mobility
		, mle2_args = list(skip.hessian=TRUE)
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

Mobility_PHAC <- do.call(calibrate_comb
	, c(nlist(params=params
		, debug_plot=FALSE
		, sim_args = list(ndt = 2)
     	, data = trim_dat
		# , mle2_control = list(maxit = 10)
		, mle2_args = list(skip.hessian=TRUE)
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

print(plot(Mobility_ICU1, data=calibrate_data_fill) + ggtitle("Mobility ICU1"))

print(plot(Mobility_PHAC, data=calibrate_data_fill) + ggtitle("Mobility PHAC"))

saveEnvironment()


