library(McMasterPandemic)
library(tidyverse)
library(zoo)
library(parallel)
library(furrr)
library(future.batchtools)
library(splines)

## callArgs <- "spline_recalib.Rout spline_recalib.R batchtools.rda spline.csv"

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

R0 <- 3
splinedf <- 3
fitmax <- 150
spline_pars <- c(0.5,-0.3,0.2)


params <- read_params(matchFile(".csv$"))

## construct intercept and spline basis
X <- ns(1:fitmax,df=splinedf)
start_date <- as.Date("2020-01-01")
end_date <- start_date -1 + fitmax
pred_days <- 30
obs_disp <- 50

## Readjust beta0 and other parameters wrt R0
params <- fix_pars(params, target=c(R0=R0))
params["obs_disp"] <- obs_disp

## I just care about the shape here right?
time_beta <- c(spline_pars)

print(plot(exp(X %*% time_beta)))

opt_pars <- list(params=c(log_beta0 = as.numeric(log(params["beta0"]))
# , log_E0=log(as.numeric(params["E0"]))
		)
	# , log_nb_disp = c(report=3,death=3)
)

opt_parsE0 <- list(params=c(log_beta0 = as.numeric(log(params["beta0"]))
	, log_E0=log(as.numeric(params["E0"]))
)
# , log_nb_disp = c(report=3,death=3)
)


sim_recalib <- function(x){
	set.seed(x)

ddfull_sim <- (run_sim_loglin(
	sim_args=list(start_date=start_date
		, end_date = end_date + pred_days
		, stoch = c(obs = TRUE, proc = FALSE)
		)
	, extra_pars = list(time_beta = time_beta)
	, time_args = list(X=X
		, X_date=start_date:end_date
		)
	, params = params
	) 
	%>% gather(key = "var", value = "value", -date)
	%>% filter(var %in% c("report","death"))
	%>% mutate(value = round(value))
)

dd_sim <- (ddfull_sim
	%>% filter(date <= end_date)
	)

dd_sim

ff <- calibrate_comb(params = params
	, debug_plot=FALSE
	, use_DEoptim=TRUE
	, DE_cores = 6
	, opt_pars = opt_pars
	, use_spline = TRUE
	, spline_df = splinedf
	, spline_type = "ns"
	, data= dd_sim
	, start_date = min(dd_sim$date)
	, start_date_offset = 0
)

ffE0 <- calibrate_comb(params = params
	, debug_plot=FALSE
	, use_DEoptim=TRUE
	, DE_cores = 6
	, opt_pars = opt_parsE0
	, use_spline = TRUE
	, spline_df = splinedf
	, spline_type = "ns"
	, data= dd_sim
	, start_date = min(dd_sim$date)
	, start_date_offset = 0
)

ff_list <- list(fit=ff, fitE0=ffE0, fitdat=dd_sim, full_dat=ddfull_sim)

saveRDS(object=ff_list, file=paste0("./cachestuff/spline_recalib.",x,".RDS"))
}

batch_setup()

future_map(1:10,function(x)sim_recalib(x))

