## This test break points where it changes beta0
## We are assuming the effect to happen immediately right after the break date to incidence
## 
library(McMasterPandemic)
library(tidyverse)
library(zoo)
library(anytime)
library(ggplot2); theme_set(theme_bw())

##  load parameters
params <- fix_pars(read_params("ICU1.csv")
   , target=c(Gbar=6)
   , pars_adj=list(c("sigma","gamma_s","gamma_m","gamma_a"))
)
                   
first_date <- as.Date("2020-03-15")
end_date <- as.Date("2020-03-25")

## breakpoints
change1 <- anydate(first_date + 19)
change2 <- anydate(first_date + 29)

dd <- data.frame(date = as.Date(first_date:(first_date+39))
   , I = c(seq(2,2000,length.out = 20)
      , seq(1800,1500,length.out=10)
      , seq(1600,5000,length.out=10)
      )
   , var = "report"
)

## convoluting incidence to get report
## we need to round because we are fitting negbinom and getting rid of the NAs in the beginning
dd$value <- round(calc_conv(dd$I,params))
dd <- dd %>% filter(!is.na(value))

## setting some population size
params[["N"]] <- 10e6  

## The break points actually don't match the reports because we want to match incidence
(ggplot(dd,aes(x=date,y=value))
    + geom_line()
    + geom_point()
    + geom_vline(xintercept=as.Date(c(change1,change2)))
)

bd <- anydate(c(change1,change2))
## print(bd)

opt_pars <- list(
    ## these params go to run_sim
    params=c(log_E0=4
             , log_beta0=-1
             ## fraction of mild (non-hosp) cases
             , log_mu=log(params[["mu"]])
             ## fraction of incidence reported
             ## logit_c_prop=qlogis(params[["c_prop"]]),
             ## fraction of hosp to acute (non-ICU)
             , logit_phi1=qlogis(params[["phi1"]])
             ## fraction of ICU cases dying
             ## logit_phi2=qlogis(params[["phi2"]])
    ),
    log_rel_beta0 = rep(-1, length(bd)),
    log_nb_disp=0)

t1 <- system.time(g1 <- calibrate(data=dd
   , base_params=params
   , optim_args=list(control=list(maxit=10000),hessian=TRUE)
   , opt_pars = opt_pars
   , break_dates = bd
   # , debug=TRUE
   # , debug_plot = TRUE
)
) ## system.time


f_args <-attr(g1,"forecast_args")
i1 <- invlink_trans(restore(g1$par,f_args$opt_pars))

fc <- (do.call(forecast_sim,c(list(p=g1$par), f_args))) 
forecast <- fc %>% filter(var %in% c("report","incidence"))
(ggplot(forecast,aes(x=date,y=value,colour=var))
    + geom_line()
    + geom_vline(xintercept = as.Date(c(change1,change2)))
    + scale_x_date(breaks="2 day")
    + scale_y_log10()
)
