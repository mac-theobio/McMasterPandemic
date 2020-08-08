library(McMasterPandemic)
library(tidyverse)

source("makestuff/makeRfuns.R")
commandEnvironments()

load("cachestuff/Ontario_calibration_setup.rda")

calibrate_data_testify <- (calibrate_data_fill
		%>% mutate(var = ifelse(var == "report", "postest", var))
	## first intensity cannot be zero
		%>% filter(date >= test_data_fill$Date[which(test_data_fill$intensity>0)[1]])
)

test_data_fill <- test_data_fill %>% filter(Date >= test_data_fill$Date[which(test_data_fill$intensity>0)[1]])


## loading parameters

params <- fix_pars(read_params("PHAC_testify.csv")
	, target=c(Gbar=6)	
   , pars_adj=list(c("sigma","gamma_s","gamma_m","gamma_a"))
)

params[["N"]] <- 14.57e6 ## Population of Ontario (2019)

### What parameters do we want to optimize?

opt_pars <- list(params=c(log_E0=2, log_beta0=-1, logit_phi1 = -1),log_nb_disp=c(report=1, death=1,H=1))

current <- do.call(calibrate_comb
	, c(nlist(params=params
		, debug_plot=TRUE
		, sim_args = list(ndt = 1, ratemat_args = list(testify=TRUE,testing_time="report"))
		# , mle2_args = list(skip.hessian=TRUE)
     	, data = calibrate_data_testify
     	, mob_data = clean_mobility
     	, opt_pars = opt_pars
     	, use_DEoptim = FALSE
		, DE_cores = 1
		, use_phenomhet = TRUE
		, use_mobility = TRUE
		, mob_breaks = "2020-04-15"
		, mob_breaks_int = TRUE
		, mob_logist_scale = 3
		, use_spline = FALSE
		, spline_type = "ns"
		, spline_df = NA
		, spline_days = 14
		, testing_data = test_data_fill
		, use_testing = TRUE
		)
	)
)

print(plot(current, data=calibrate_data_fill) + ggtitle("Current model: PH + mobility cap"))


saveEnvironment()


