library(McMasterPandemic)
library(zoo)
library(tidyverse)

## Copying example from run_sim and modifying it

##1. testing if Relative_value=1 and non-timevar run_sim are the same

params <- read_params("ICU1.csv")
paramsS <- update(params,c(proc_disp=0.1,obs_disp=100))
paramsSz <- update(paramsS, zeta=5)
state <- make_state(params=params)
startdate <- as.Date("2020-01-01")
enddate <- as.Date("2020-05-01")
time_pars <- data.frame(Date=as.Date(startdate:enddate),
								Symbol="beta0",
								Relative_value=1,
								stringsAsFactors=FALSE)

## This is checking if we can get the same thing if we don't add stoch

res1 <- run_sim(params,state,start_date=startdate,end_date=enddate)
res2 <- run_sim(params,state,start_date=startdate,end_date=enddate)
res1 == res2
all.equal(res1,res2)


## This fits a timevar dataframe where beta0 = 1 
res1_t <- update(res1, params_timevar=time_pars)


all.equal(res1,res1_t)
## Only attributes are different?
all.equal(res1$report,res1_t$report)


## Checking predict and run_sim from mre

modlist <- readRDS("ON.short.spline.RDS")

spline_lims <- range(modlist$fit$forecast_args$time_args$X_date)

pp <- predict(modlist$fit,ensembles=FALSE,
				  end_date="2020-12-01",
				  keep_vars = c("hosp", "death"))

## This is already using mli stitching code

print(gg <- ggplot(pp,aes(date,value))
		+ geom_line()
		+ facet_wrap(~var, nrow=1,scale="free")
)



X <- modlist$fit$forecast_args$time_args$X
X_date <- modlist$fit$forecast_args$time_args$X_date
time_beta <- coef(modlist$fit,"fitted")$time_beta

rel_beta <- exp(X %*% matrix(time_beta,ncol=1))

time_pars <- data.frame(Date=as.Date(min(pp$date):max(pp$date)),
								stringsAsFactors=FALSE)


rel_betaf <- data.frame(Date = X_date
	, Symbol = "beta0"
	, Relative_value=rel_beta)

freeze_rel_value <- rel_betaf$Relative_value[nrow(rel_betaf)]

## let start_date_offset have relative value = 1 (i.e. this is the time before the first data point)
## freeze last relative value in fit and extend it for forecast period

time_pars_calc <- (left_join(time_pars,rel_betaf)
	%>% mutate(Relative_value = ifelse(Date < min(rel_betaf$Date),1,Relative_value)
			, Relative_value=ifelse(Date>max(rel_betaf$Date),freeze_rel_value,Relative_value)
		
		, Symbol = "beta0")	
)

## manually creating the stitching beta(t) and use timevar to simulate it

pp2 <- run_sim(params=coef(modlist$fit,"all")
					, state = make_state(params = coef(modlist$fit,"all"))
					,start_date=min(pp$date)
					,end_date=max(pp$date)
					, ndt = 2
					, params_timevar = time_pars_calc)

pp2 <- (pp2 %>% select(date,hosp,death)
	%>% gather(var,value,-date))

print(gg %+% pp2)


## This is experimenting if we only restrict to X_date, which is the data time series (i.e. nothing for start date offset/pre data period and nothing for forecast period)
## I.e. we are only taking the middle chunk (the data series chunk)

pp3 <- run_sim(params=coef(modlist$fit,"all")
					, state = make_state(params = coef(modlist$fit,"all"))
					,start_date=min(pp$date)
					,end_date=max(pp$date)
					, ndt = 2
					, params_timevar = rel_betaf)

pp3 <- (pp3 %>% select(date,hosp,death)
		  %>% gather(var,value,-date))

print(gg %+% pp3)

