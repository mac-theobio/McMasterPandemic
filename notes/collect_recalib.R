library(McMasterPandemic)
library(tidyverse)

source("makestuff/makeRfuns.R") ## Will eventually be a package
commandEnvironments() ## Read in any environments specified as dependencies
## makeGraphics()


flist <- list.files(path="cachestuff/",pattern="recalib.no")

print(flist)

collect_pars <- function(x){
	modlist <- readRDS(paste0("cachestuff/",x))
	cc <- coef(modlist,"fitted")
	parsdf <- data.frame(beta0 = cc$params[1]
		, E0 = cc$params[2]
		, seed = x
		, type = "sim"
		)
}

pars_df <- bind_rows(lapply(flist,collect_pars))
true_pars <- coef(ff,"fitted")
true_pars_df <- data.frame(beta0 = true_pars$params[1]
	, E0 = true_pars$params[2]
	, seed = NA
	, type = "true"
)


combo_pars <- (bind_rows(true_pars_df, pars_df)
	%>% gather(key = "var", value = "value", -seed, -type)
)

### spline shape

print(X)

collect_splines <- function(x){
	modlist <- readRDS(paste0("cachestuff/",x))
	cc <- coef(modlist,"fitted")
	spline_df <- data.frame(time = 1:nrow(X)
	, bt = exp(X[,-1] %*% matrix(cc$time_beta, ncol=1))
	, seed = x 
	, type = "sim"
	)
}


spline_df <- bind_rows(lapply(flist,collect_splines))

tp <- coef(ff,"fitted")

true_splines <- data.frame(time=1:nrow(X)
	, bt = exp(X[,-1] %*% matrix(tp$time_beta,ncol=1))
	, seed = NA
	, type = "true"
)

spline_df <- bind_rows(spline_df, true_splines)

saveVars(combo_pars, spline_df)

