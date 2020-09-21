library(McMasterPandemic)
library(tidyverse)

source("makestuff/makeRfuns.R") ## Will eventually be a package
commandEnvironments() ## Read in any environments specified as dependencies
## makeGraphics()


flist <- list.files(path="cachestuff/",pattern="recalib")

print(flist)

collect_pars <- function(x){
	modlist <- readRDS(paste0("cachestuff/",x))
	cc <- coef(modlist,"fitted")
	parsdf <- data.frame(beta0 = cc$params[1]
		, E0 = cc$params[2]
		, seed = x
	)
	parsdf <- (parsdf
		%>% mutate(seed = gsub(".RDS",replace = "", seed)
			, seed = gsub("spline_recalib.",replace = "", seed)
		)
		%>% gather(key="param",value="value",-seed)
	)
}

pars_df <- bind_rows(lapply(flist,collect_pars))
true_pars <- coef(ff,"fitted")
true_pars$params


### spline shape

# collect_splines <- function(x){
# 	modlist <- readRDS(paste0("cachestuff/",x))
# 	cc <- coef(modlist,"fitted")
# 	spline_df <- data.frame(time = 1:nrow(X)
# 	, calib_spline = exp(X[,-1] %*% matrix(aa$time_beta, ncol=1))
# 	, true_spline = exp(X[,-1] %*% matrix(bb[-1], ncol=1))
# 	)
# 	
# }
