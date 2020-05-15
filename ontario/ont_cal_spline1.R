library(McMasterPandemic)
library(splines)

load("ontario_clean.RData")
load("ont_cal.RData")  ## baseline calibration: ont_cal1 (calibrated object), bd (breakpoint dates)

print(unique(ont_all_sub$var))
ont_noICU <- dplyr::filter(ont_all_sub, var != "ICU")

## penalization on beta parameters ... could be working but not sure?
pen_prior <- list(~sum(dnorm(time_beta,mean=0,sd=1)))
## stronger penalization
pen_prior2 <- list(~sum(dnorm(time_beta,mean=0,sd=0.5)))

X_date <- unique(ont_noICU$date)
t_vec <- as.numeric(X_date-min(X_date))
kvec <- seq(7,length(t_vec),by=14)
## spline basis with no intercept (the intercept is the main (log_)beta0 parameter)
X <- model.matrix(~ns(t_vec,knots=kvec)-1)
## take a look at the spline basis
matplot(t_vec,X,type="l")

opt_pars_spline1 <- list(
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
  , time_beta=rep(0,ncol(X))   ## start linear params from 0
    ## NB dispersion (filled with log_nb_disp=0 for each state variable, internally)
  , log_nb_disp=NULL
)

## single test sim
run_sim_loglin(params=params, extra_pars=opt_pars_spline1["time_beta"],
               time_args=nlist(X,X_date),
               sim_args=list(start_date=min(ont_all_sub$date)-15,
                             end_date=max(ont_all_sub$date)))

## do the calibration: Nelder-Mead only
t_ont_cal_spline1_noDE <- system.time(ont_cal_spline1_noDE <-
                                     calibrate(data=ont_noICU
                                             , use_DEoptim=FALSE
                                             , debug_plot = TRUE
                                             , base_params=params
                                             , opt_pars = opt_pars_spline1
                                             , time_args=nlist(X,X_date)  ## model matrix and date vector
                                             , sim_fun = run_sim_loglin
                                        )
                                 ) ## system.time

## penalize log-lin parameters (joint ridge penalty)
t_ont_cal_spline1_noDE_pen <- system.time(ont_cal_spline1_noDE_pen <-
                                              update(ont_cal_spline1_noDE, prior=pen_prior2)
                                          )        

## now try with DEoptim (doesn't seem like a big improvement)
## switched parameter bounds for spline from [-3,3] to [-2,2] (still too loose?)
## NOTE there are some problems with update() + DE_cores>1
t_ont_cal_spline1 <- system.time(ont_cal_spline1 <-
                                     update(ont_cal_spline1_noDE
                                         , DE_cores=1
                                         , prior=pen_prior
                                         , use_DEoptim=TRUE))

