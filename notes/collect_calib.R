library(McMasterPandemic)
library(tidyverse)

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
#	bt <- exp(X %*% matrix(cc$time_beta, ncol=1))
}

#bt <- btfun(
#	cc = coef(tempmod$fit,"fitted")
#	, X = tempmod$fit$forecast_args$time_args$X
#)

#print(bt)

collect_pars <- function(x){
	modlist <- readRDS(paste0("cachestuff/",x))
	cc <- coef(modlist$fit,"fitted")
	parsdf <- data.frame(beta0 = cc$params[1]
		, E0 = cc$params[2]
		, seed = x
		, type = "sim"
		, mod = "withoutE0"
		)
	cc2 <- coef(modlist$fitE0,"fitted")
	parsdf2 <- data.frame(beta0 = cc2$params[1]
								, E0 = cc2$params[2]
								, seed = x
								, type = "sim"
								, mod = "withE0"
	)
	return(bind_rows(parsdf,parsdf2))
}

pars_df <- bind_rows(lapply(flist,collect_pars))

print(pars_df)


true_pars_df <- data.frame(beta0 = base_params["beta0"]
	, E0 = base_params["E0"]
	, seed = NA
	, type = "true"
	, mod = "true"
)


combo_pars <- (bind_rows(true_pars_df, pars_df)
	%>% gather(key = "var", value = "value", -seed, -type, -mod)
)

### spline shape

X <- cbind(1,tempmod$fit$forecast_args$time_args$X)
print(X)

collect_splines <- function(x){
	modlist <- readRDS(paste0("cachestuff/",x))
	R0t <- summary(modlist$fit)$R0
	cc <- coef(modlist$fit,"fitted")
	spline_df <- (data.frame(time = 1:nrow(X)
		, bt = btfun(cc=coef(modlist$fit,"fitted"),X=modlist$fit$forecast_args$time_args$X)/base_params["beta0"]
		, seed = x
		, Rt = R0t[-1]
		, type = "sim"
		, mod = "withoutE0"
	))
	cc2 <- coef(modlist$fitE0,"fitted")
	R0t2 <- summary(modlist$fitE0)$R0
	spline_df2 <- (data.frame(time = 1:nrow(X)
		, Rt = R0t2[-1]
		, bt = btfun(cc=coef(modlist$fitE0,"fitted"),X=modlist$fitE0$forecast_args$time_args$X)/base_params["beta0"]
		, seed = x 
		, type = "sim"
		, mod = "withE0"
	))
	return(bind_rows(spline_df,spline_df2))
}


spline_df <- bind_rows(lapply(flist,collect_splines))

## copied from spline_recalib.R 

true_splines <- data.frame(time=1:nrow(X)
	, bt = bt
	, Rt = Rt
	, seed = NA
	, type = "true"
	, mod = "true"
)

spline_df <- bind_rows(spline_df, true_splines)

print(spline_df)

saveVars(combo_pars, spline_df, rmult, base_params)

