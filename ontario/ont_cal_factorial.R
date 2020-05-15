library(McMasterPandemic)
library(splines)
library(dplyr)

## flags
use_mobility <- FALSE
use_spline <- FALSE
use_zeta <- TRUE
use_DEoptim <- FALSE

load("ontario_clean.RData")
load("ont_cal.RData")  ## baseline calibration: ont_cal1 (calibrated object), bd (breakpoint dates)

## select the part of the mobility data that lies within the calibration data set
comb_sub2 <- (comb_sub
    %>% right_join(ont_all %>% select(date) %>% unique(),by="date")
    %>% na.omit()
    %>% mutate_at("rel_activity",pmin,1) ## cap relative mobility at 1
)

print(unique(ont_all_sub$var))
ont_noICU <- dplyr::filter(ont_all_sub, var != "ICU")

X_date <- unique(ont_noICU$date)
comb_sub3 <- (full_join(unique(select(ont_noICU,date)),comb_sub2,by="date")
    ## assume missing data = constant at last value
    %>% mutate_at("rel_activity",zoo::na.locf)
    %>% mutate(t_vec=as.numeric(date-min(date)))
)

plot(rel_activity~date,data=comb_sub3)

opt_pars <- list(
    ## these params are part of the main parameter vector: go to run_sim()
    params=c(log_E0=4      ## initial exposed
           , log_beta0=-1  ## initial baseline transmission
             ## fraction of mild (non-hosp) cases
           , log_mu=log(params[["mu"]])
             ## fraction of incidence reported
             ## logit_c_prop=qlogis(params[["c_prop"]]),
             ## fraction of hosp to acute (non-ICU)
           , logit_phi1=qlogis(params[["phi1"]])
             )
  , log_nb_disp=NULL
)

## now modify opt_pars according to flags
if (use_zeta) opt_pars$params <- c(opt_pars$params,log_zeta=-1)
loglin_terms <- "-1"
if (use_mobility) loglin_terms <- c(loglin_terms, "log(rel_activity)")
if (use_spline) loglin_terms <- c(loglin_terms, "bs(t_vec,df=round(length(t_vec)/14))")
form <- reformulate(loglin_terms)
X <- model.matrix(form, data = comb_sub3)

matplot(comb_sub3$t_vec,X,type="l")

opt_pars$time_beta <- rep(0,ncol(X))  ## mob-power is incorporated (param 1)
time_args <- nlist(X,X_date=comb_sub3$date)

run_sim_loglin(params=params, extra_pars=opt_pars["time_beta"] ## n.b. single brackets here are on purpose
             , time_args=time_args
             , sim_args=list(start_date=min(ont_all_sub$date)-15
                           , end_date=max(ont_all_sub$date)))

## do the calibration
t_ont_cal_fac <- system.time(ont_cal_fac <-
                                 calibrate(data=ont_noICU
                                         , use_DEoptim=use_DEoptim
                                         , debug_plot = TRUE
                                         , base_params=params
                                         , opt_pars = opt_pars
                                         , time_args=time_args
                                         , sim_fun = run_sim_loglin
                                           )
                                 ) ## system.time
