library(tidyverse)
library(McMasterPandemic)
library(zoo)
library(parallel)

source("makestuff/makeRfuns.R")
commandEnvironments()

flist <- list.files(path="cachestuff/",pattern="ont_spline.short")

print(flist)

forecast_dat <- function(x){
	modlist <- readRDS(paste0("cachestuff/",x))
	mle2Sigma <- bbmle::vcov(modlist$fit$mle2)
	ensembles <- predict(modlist$fit
		, ensemble = TRUE
		, start_date = min(modlist$fitdat$date) - modlist$start_date_offset
		, end_date = max(modlist$fitdat$date)+ 30
		, stoch = c(proc =TRUE, obs = TRUE)
		, Sigma = mle2Sigma
		, nsim = 200
		, new_params = c(obs_disp=5, proc_disp=1)
		, stoch_start = c(proc=max(modlist$fitdat$date)+1, obs=min(modlist$fitdat$date)-modlist$start_date_offset)
		, keep_vars = c("report","Rt")
	)
	
	combodat <- (ensembles
		%>% transmute(date = date
			, var
			, lwr
			, med = value
			, upr
			, spline_df = modlist$spline_params$spline_df
			, spline_pen = modlist$spline_params$spline_pen
			, trim = ifelse(grepl("short",x),"short","full")
			, convergence_code = modlist$fit$mle2@details$convergence
		)
		%>% left_join(.,modlist$fitdat)
	)
	return(combodat)
}

ensembles_list <- mclapply(X=flist,FUN=forecast_dat,mc.cores = 3)

ensembles_dat <- bind_rows(ensembles_list)

saveVars(ensembles_dat)



