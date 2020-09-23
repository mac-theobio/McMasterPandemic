library(McMasterPandemic)
library(tidyverse)
library(zoo)
library(parallel)
library(furrr)
library(future.batchtools)

## callArgs <- "spline_recalib.Rout spline_recalib.R batchtools.rda spline_fit.rda spline.csv"

source("makestuff/makeRfuns.R")
commandEnvironments()
makeGraphics()

R0 <- 2
ndf <- 5

params <- read_params(matchFile(".csv$"))

## FIXME: construct X here independently from mod_ns
X <- cbind(1,mod_ns$model[,-1])
start_date <- as.Date("2020-01-01")
end_date <- start_date -1 + nrow(mod_ns$model)
pred_days <- 30
obs_disp <- 3

params <- fix_pars(params, target=c(R0=R0))
params["obs_disp"] <- obs_disp 

## I should just worry about the shape and not the intercept right?
time_beta <- c(1,as.numeric(coef(mod_ns)[-1]))

print(plot(exp(X %*% time_beta)))

opt_pars <- list(params=c(log_beta0 = log(params["beta0"])
# , log_E0=log(as.numeric(cc$params[2]))
		)
	, log_nb_disp = c(report=3,death=3)
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

ff_refit <- calibrate_comb(params = params
	, debug_plot=FALSE
	, use_DEoptim=FALSE
	, DE_cores = 6
	, opt_pars = opt_pars
	, use_spline = TRUE
	, spline_df = ndf
	, spline_type = "ns"
	, data= dd_sim
	, start_date = min(dd_sim$date)
	, start_date_offset = 0
)				

ff_list <- list(fit=ff_refit,fitdat=dd_sim, full_dat=ddfull_sim)

saveRDS(object=ff_list, file=paste0("./cachestuff/spline_recalib.",x,".RDS"))
}

batch_setup()

future_map(1:10,function(x)sim_recalib(x))

