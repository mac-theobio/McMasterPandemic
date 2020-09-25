library(McMasterPandemic)
library(dplyr)


# callArgs <- "spline_forecast.Rout spline_forecast.R collect_recalib.rda cachestuff/spline_calib.rda"

source("makestuff/makeRfuns.R")
commandEnvironments()


flist <- list.files(path="cachestuff/",pattern="spline_recalib")

forecasting <- function(x){
	
	modlist <- readRDS(paste0("cachestuff/",x))
	truedat <- (modlist$full_dat
		%>% filter(var == "report")
		%>% mutate(seed = x
			, value = ifelse(is.na(value),NA,value)
			, type = "true"
			, mod = "true"
			, lwr = value
			, upr = value
			)
	)
	
	sim_forecast <- (predict(readRDS(paste0("cachestuff/",x))$fit,ensemble=TRUE,end_date=max(truedat$date)
		, stoch=c(proc=FALSE,obs=TRUE))
	%>% filter(var == "report")
	%>% mutate(seed = x
		, type = "sim"
		, mod = "withoutE0"
		)
	%>% select(-vtype)
)
	
	sim_forecast2 <- (predict(readRDS(paste0("cachestuff/",x))$fitE0,ensemble=TRUE,end_date=max(truedat$date)
									 , stoch=c(proc=FALSE,obs=TRUE))
						  %>% filter(var == "report")
						  %>% mutate(seed = x
						  			  , type = "sim"
						  			  , mod = "withE0"
						  )
						  %>% select(-vtype)
	)
	
	
combodat <- bind_rows(truedat,sim_forecast, sim_forecast2)

return(combodat)

}

forecastdat <- bind_rows(lapply(flist,function(x)forecasting(x)))

saveVars(forecastdat)




