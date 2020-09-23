library(McMasterPandemic)
library(tidyverse)
library(zoo)
library(parallel)
library(furrr)
library(future.batchtools)

## callArgs <- "spline_recalib.Rout spline_recalib.R batchtools.rda cachestuff/spline_calib.rda spline_fit.rda spline.csv"


source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

R0 <- 2
fitmax <- 75 ## Copied from dependency 
ndf <- 3

params <- read_params(matchFile(".csv$"))


X <- cbind(1,mod_ns$model[,-1])
start_date <- as.Date("2020-01-01")
end_date <- start_date -1 + nrow(mod_ns$model)
time_beta <- c(0,as.numeric(coef(mod_ns)[-1]))

params <- fix_pars(params, target=c(R0=R0))

opt_pars <- list(params=c(log_beta0 = log(params["beta0"])
# , log_E0=log(as.numeric(cc$params[2]))
		)
	, log_nb_disp = c(report=3,death=3)
)


sim_recalib <- function(x){
	set.seed(x)

dd_sim <- (run_sim_loglin(
	sim_args=list(start_date=start_date
		, end_date = end_date
		# , time_beta = time_beta
		)
	, extra_pars = list(time_beta = time_beta)
	, time_args = list(X=X
		, X_date=start_date:end_date
		)
	, params = params
	) 
	%>% filter(var %in% c("report","death"))
	%>% mutate(value = round(value))
)
	

ff_refit <- calibrate_comb(params = params
	, debug_plot=TRUE
	, use_DEoptim=FALSE
	, DE_cores = 3
	, opt_pars = opt_pars
	, use_spline = TRUE
	, spline_df = ndf
	, spline_type = "ns"
	, data= dd_resim
	, start_date = min(dd_resim$date)
	, start_date_offset = 0
)				

ff_list <- list(fit=ff_refit,data=dd_resim)

saveRDS(object=ff_refit, file=paste0("./cachestuff/spline_recalib.",x,".RDS"))
}

batch_setup()

future_map(1:10,function(x)sim_recalib(x))

