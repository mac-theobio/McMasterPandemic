library(McMasterPandemic)
library(tidyverse)

source("makestuff/makeRfuns.R")
## Don't worry about this when coding interactively/sourcing
commandFiles()

## 1 Downloading latest data (note: the two models above did not fit up to the very last datapoint)
### Ideally, we should be able to do this via source, but makeR comandFiles() is giving us problems.
### source("Ontario_calibration.R")


tsdat_url <- "https://wzmli.github.io/COVID19-Canada/git_push/clean.Rout.csv"

tsdat <- read_csv(tsdat_url)

Ontario_dat <- (tsdat
	%>% filter(Province=="ON")
   %>% select(Province,Date,Hospitalization,ICU,Ventilator,deceased,newConfirmations,newTests)
	%>% mutate(newDeaths=c(NA,diff(deceased))
   	## ON hosp includes ICU, our model compartment is just acute care
   	, Hospitalization=Hospitalization-ICU)
   %>% select(-deceased)
   %>% pivot_longer(names_to="var",-c(Date,Province))
   %>% setNames(tolower(names(.)))
	%>% ungroup()
)

keep_vars <- c("H","ICU","death","report","newTests")

newdat <- (Ontario_dat
	%>% mutate_at("var",trans_state_vars)
   %>% filter(var %in% keep_vars)
	
)

## 2 Loading refernce model (PH, mobility)

load("cachestuff/Ontario_basic.rda")

PH_mob <- Ontario_fit
	
## 3 Loading refernce model + splines (PH, mobility, splines)

load("cachestuff/Ontario_calibration_spline.rda")
PH_mob_spline <- Ontario_fit_spline

## Default plots

plot(PH_mob, data=newdat)


plot(PH_mob_spline, data=newdat)


## Forecasting

nsim = 50
keep_vars <- c("Rt","H","report","death","ICU","incidence")
end <- Sys.Date()
prediction_offset <- 15
epistart <- min(newdat$date) - prediction_offset

function(mod,importance=FALSE,modtype){
ensemble1 <- predict(object=mod, ensemble=TRUE
        , stoch=c(proc=FALSE,obs=TRUE)
        , new_params=c(obs_disp=5.0)
        , end_date=end
        # , stoch_start = c(obs=epistart)
        , nsim=nsim
        , keep_vars = keep_vars
        , qvec=NULL
        , Sigma = PH_mob$mle2@vcov
)

predf1 <- (reshape2::melt(ensemble1[keep_vars,,])
    %>% as_tibble()
    %>% mutate_at("date",as.Date)
    %>% mutate_at("var",as.character)
)


## set up a realization wts df
wts_df <- data.frame(sim = 1:nsim, wts=1)
if(importance){
        wts_df$wts <- attr(ensemble1,"imp_wts")
} 

keepdat <- (predf1
        %>% filter(var %in% c("H","ICU","Rt","incidence","report","death"))      
)
modsum <- (keepdat
        %>% left_join(.,wts_df)
        %>% group_by(var,date)
        %>% summarise(lwr = Hmisc::wtd.quantile(value, probs=q[[1]], weights=wts)
                        , med = Hmisc::wtd.quantile(value, probs=q[[2]], weights=wts)
                        , upr = Hmisc::wtd.quantile(value, probs=q[[3]], weights=wts)
        )
        %>% mutate(mod = modtype
        )
)




}

saveEnvironment()


