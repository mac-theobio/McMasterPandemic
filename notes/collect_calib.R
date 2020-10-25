library(McMasterPandemic)
library(tidyverse)

callArgs <- "collect_calib.Rout collect_calib.R spline_sim.rda cachestuff/simcalib.1.RDS cachestuff/spline_calib.8.RDS cachestuff/spline_calib.9.RDS cachestuff/spline_recalib.10.RDS cachestuff/simcalib.2.RDS cachestuff/spline_calib.3.RDS cachestuff/spline_recalib.7.RDS cachestuff/spline_recalib.8.RDS cachestuff/spline_calib.6.RDS cachestuff/spline_calib.4.RDS cachestuff/spline_recalib.3.RDS cachestuff/spline_calib.1.RDS cachestuff/spline_recalib.9.RDS cachestuff/spline_recalib.4.RDS cachestuff/spline_recalib.1.RDS cachestuff/spline_calib.5.RDS cachestuff/spline_recalib.6.RDS cachestuff/spline_calib.2.RDS cachestuff/spline_recalib.5.RDS cachestuff/spline_recalib.2.RDS cachestuff/spline_calib.7.RDS"


source("makestuff/makeRfuns.R") ## Will eventually be a package
commandEnvironments() ## Read in any environments specified as dependencies
## makeGraphics()

flist <- list.files(path="cachestuff/",pattern="_calib[.]")

print(flist)

## true beta0 is the same, it lives inside base_params

tempmod <- readRDS(paste0("cachestuff/",flist[1]))
base_params <- tempmod$fit$forecast_args$base_params

## We expect the ratio of R/β to stay constant
rmult <- get_R0(base_params)/base_params[["beta0"]]
print(rmult)

print(summary(tempmod$fit)$R0)

## Calculate time-varying betas
btfun <- function(cc, X){
	bt <- cc$params[["beta0"]] * exp(X %*% matrix(cc$time_beta, ncol=1))
}

collect_pars <- function(x){
	modlist <- readRDS(paste0("cachestuff/",x))
	mods <- c("fit","fitpen","fitE0","fitE0pen")
	df_funs <- function(seed,y){
		tempmod <- modlist[[y]]
		cc <- coef(tempmod, "fitted")
		parsdf <- data.frame(beta0 = cc$params[["beta0"]]
			, E0 = as.numeric(cc$params[2]) ## cannot say E0, it will give an error for mods without E0
#			, tbint = cc$time_beta[1]
			, tb1 = cc$time_beta[1]
			, tb2 = cc$time_beta[2]
			, tb3 = cc$time_beta[3]
			, tb4 = cc$time_beta[4]
			, tb5 = cc$time_beta[5]
			, tb6 = cc$time_beta[6]
			, tb7 = cc$time_beta[7]
#			, beta0sum = cc$time_beta[1] + cc$params[["beta0"]]
			, seed = seed
			, type = "sim"
			, mod = y
		)
		return(parsdf)
	}
	dflist <- lapply(mods,function(y){df_funs(seed=x,y)})
	return(bind_rows(dflist))
}

pars_df <- bind_rows(lapply(flist,collect_pars))

print(pars_df)


true_pars_df <- data.frame(beta0 = base_params["beta0"]
	, E0 = base_params["E0"]
#	, tbint = bb[1]
	, tb1 = bb[1]
	, tb2 = bb[2]
	, tb3 = bb[3]
	, tb4 = bb[4]
	, tb5 = bb[5]
	, tb6 = bb[6]
	, tb7 = bb[7]
#	, beta0sum = base_params["beta0"] + bb[1]
	, seed = NA
	, type = "true"
	, mod = "true"
	, spline_pen=NA
)


combo_pars <- (bind_rows(true_pars_df, pars_df)
	%>% gather(key = "var", value = "value", -seed, -type, -mod, -spline_pen)
)

### spline shape

collect_splines <- function(x){
	modlist <- readRDS(paste0("cachestuff/",x))
	mods <- c("fit","fitpen","fitE0","fitE0pen")
	df_funs <- function(seed,y){
		tempmod <- modlist[[y]]
		R0t <- summary(tempmod)[["R0"]]
		cc <- coef(tempmod, "fitted")
		spline_df <- data.frame(time = 1:nrow(X)
		, bt = btfun(cc=cc,X=tempmod$forecast_args$time_args$X)/base_params[["beta0"]]
		, seed = x
		, Rt = R0t[-1]
		, type = "sim"
		, mod = y
		)
		return(spline_df)
	}
	spline_dflist <- lapply(mods,function(y){df_funs(seed=x,y)})
	return(bind_rows(spline_dflist))
}


spline_df <- bind_rows(lapply(flist,collect_splines))

## copied from spline_recalib.R 

true_splines <- data.frame(time=1:nrow(X)
	, bt = Rt
	, Rt = Rt
	, seed = NA
	, type = "true"
	, mod = "true"
)

spline_df <- bind_rows(spline_df, true_splines)

print(spline_df)

saveVars(combo_pars, spline_df, rmult, base_params)

