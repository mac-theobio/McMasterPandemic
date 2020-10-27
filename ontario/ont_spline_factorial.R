library(McMasterPandemic)
library(tidyverse)

source("makestuff/makeRfuns.R")
commandEnvironments()

params <- fix_pars(read_params("PHAC.csv")
	, target = c(R0 = 2, Gbar=6)
)

opt_pars <- list(
	params = c(log_E0 = 2, log_beta0=-1, logit_c_prop=-1, logit_mu= -1, logit_phi1=-1)
)

params["obs_disp"] <- 50
params["N"] <- 1.45e7

start_date <- ifelse(grepl("full",targetname())
	, as.Date("2020-03-10")
	, as.Date("2020-09-01")
)

type <- ifelse(grepl("full",targetname()),"full","short")

trimdat <- filter(ont_dat,date>=start_date)

maxdf <- floor(nrow(trimdat)/14)

spline_df = maxdf - c(0, 1, 2, 4)
spline_pen = c(4,6,8,10)


splinef <- expand.grid(spline_df=pmax(spline_df,0)
	, spline_pen = spline_pen
)

print(splinef)

calibrate_factorial <- function(x){
	spline_params <- splinef[x,]
	print(spline_params)
	ff <- calibrate_comb(params = params
		, use_DEoptim=TRUE
		, DE_core = 7
		, opt_pars = opt_pars
		, use_spline = TRUE
		, spline_df = spline_params[["spline_df"]]
		, spline_type = "bs"
		, spline_int = FALSE
		, data = trimdat
		, start_date = min(trimdat)
		, start_date_offset = 0
	)
	ff_list <- list(fit=ff, fitdat=trimdat, spline_params=spline_params)
	saveRDS(object=ff_list, file=paste0("./cachestuff/ont_spline.",type,x,".RDS"))
}

batch_setup()
future_map(1:nrow(splinef),function(x)calibrate_factorial(x))





