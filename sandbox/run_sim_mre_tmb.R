library(McMasterPandemic)
library(zoo)
library(tidyverse)
library(anytime)
library(lubridate)

## Copying example from run_sim and modifying it

##1. testing if Relative_value=1 and non-timevar run_sim are the same (done)
## MLi: This means I can always construct my own timevar for any time-varying simulation
## Move this to testhat?!?

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
time_pars = data.frame(
  Date = as.Date(startdate:enddate),
  Symbol = "beta0",
  Value = 1,
  Type = "rel_orig"
)

#modlist$fit$forecast_args$time_args$break_dates %>% dput
bd = structure(c(18538, 18589, 18607, 18622, 18641, 18666), class = "Date")
#coef(modlist$fit, 'fitted')$rel_beta0 %>% dput
vl = c(0.891575471271562, 0.833574240453344, 0.997156343645233, 0.785687213490393,
       0.624149179855572, 0.878772008882617)
# coef(modlist$fit) %>% dput
pm = structure(c(beta0 = 0.551210287937185, Ca = 0.5, Cp = 0.5, Cm = 0.5,
                 Cs = 0.1, alpha = 0.2, sigma = 0.48104483915902, gamma_a = 0.276600782516436,
                 gamma_m = 0.276600782516436, gamma_s = 0.110640313006575, gamma_p = 0.8,
                 rho = 0.0555555555555556, delta = 0.4, mu = 0.973567835677866,
                 N = 13448494, E0 = 400, nonhosp_mort = 0.3, iso_m = 0, iso_s = 0,
                 phi1 = 0.531430782575431, phi2 = 0.3, psi1 = 0.1, psi2 = 0.5,
                 psi3 = 0.2, c_prop = 0.45, c_delay_mean = 11, c_delay_cv = 0.25,
                 proc_disp = 0, zeta = 0, obs_disp_report = 40, obs_disp_H = 7,
                 obs_disp_ICU = 7, obs_disp_death = 7), class = 'params_pansim')
time_pars = data.frame(
  Date = bd,
  Symbol = "beta0",
  Value = vl,
  Type = "rel_orig"
)

options(MP_warn_bad_breaks = FALSE)
flex_model = make_base_model(
  params = params,
  start_date = "2020-10-01",
  end_date = "2021-04-01",
  params_timevar = time_pars,
  do_make_state = TRUE
)

(flex_model
  %>% simulate(do_condensation = TRUE)
  %>% ggplot
  +  facet_wrap(~variable, scales = 'free')
  +  geom_line(aes(Date, value))
)

(run_sim(params, start_date = flex_model$start_date, end_date = flex_model$end_date, params_timevar = time_pars, condense = TRUE)
  %>% tidyr::pivot_longer(-date)
  %>% filter(name %in% c("death", "foi", "H", "hosp", "ICU", "incidence", "report"))
  %>% ggplot
  +  facet_wrap(~name, scales = 'free')
  +  geom_line(aes(date, value))
)

## This is checking if we can get the same thing if we don't add stoch

res1 <- run_sim(params,state,start_date=startdate,end_date=enddate)
res2 <- run_sim(params,state,start_date=startdate,end_date=enddate)
stopifnot(identical(res1,res2))

## This fits a timevar dataframe where beta0 = 1
res1_t <- update(res1, params_timevar=time_pars)

## keep only numeric values
stopifnot(identical(c(res1),c(res1_t)))

## This is the latest ON calibration using the break_date model
modlist <- readRDS("sandbox/ON.short.breaks.RDS")

print(plot(modlist$fit,data=modlist$trimdat))

## Not plotting the data exposes the spikes in death and hosp (these are the compartments where we do diff from accumulating boxes
print(plot(modlist$fit))
## FIXME: may break with testify_eigvec branch?

## Projecting to Dec 2021 using the default settings
pp <- predict(modlist$fit,ensembles=FALSE
              , end_date = "2021-12-01"
              , keep_vars = c("hosp", "death","report","Rt")
)
## warning about testing time?

print(gg <- ggplot(pp,aes(date,value))
      + geom_line()
      + facet_wrap(~var,scale="free",nrow=2)
)

## We can see the spikes and Rt does not account for depletion of S


## Now I am going to manually construct the break_date timevar dataframe and use run_sim to simulate

break_date <- modlist$fit$forecast_args$time_args$break_dates
rel_beta <- coef(modlist$fit,"fitted")$rel_beta0

## including the first date and last date (pp is the simulation above using the default settings, I am just using it to set up the same projection window

break_date2 <- c(min(pp$date),break_date,max(pp$date))

## Before the first break_date, relative_val = 1
rel_beta2 <- c(1,rel_beta)

rel_betaf <- data.frame()
for(i in 1:(length(break_date2)-1)){
  tempdf <- data.frame(Date = seq.Date(from = break_date2[i], to = break_date2[i+1], by=1)
                       , Symbol = "beta0"
                       , Relative_value = rel_beta2[i]
  )
  rel_betaf <- rbind(rel_betaf,tempdf)
}

time_pars = data.frame(
  Date = break_date,
  Symbol = "beta0",
  Value = rel_beta,
  Type = "rel_orig"
)

## manually  and use timevar to simulate it

pp2 <- run_sim(params=coef(modlist$fit,"all")
               , state = make_state(params = coef(modlist$fit,"all"))
               ,start_date=min(pp$date)
               ,end_date=max(pp$date)
               , params_timevar = rel_betaf
)

pp3 <- run_sim(params=coef(modlist$fit,"all")
               , state = make_state(params = coef(modlist$fit,"all"))
               ,start_date=min(pp$date)
               ,end_date=max(pp$date)
               , params_timevar = time_pars
)

pp2 <- (pp2 %>% select(date,hosp,death,report)
        %>% gather(var,value,-date))
pp3 <- (pp3 %>% select(date,hosp,death,report)
        %>% gather(var,value,-date))

print(gg %+% pp2)
print(gg %+% pp3)

## compare modlist$fit params_timevar component with rel_betaf
## No more spikes!!

## 1. Allow MLi hack (construct timepars by hand) to work with ensembles
## 2. fix whatever's causing the spikes in the first place
