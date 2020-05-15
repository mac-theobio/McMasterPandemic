library(McMasterPandemic)

load("ontario_clean.RData")
load("ont_cal.RData")  ## baseline calibration: ont_cal1 (calibrated object), bd (breakpoint dates)

print(unique(ont_all_sub$var))
ont_noICU <- dplyr::filter(ont_all_sub, var != "ICU")

opt_pars_zeta1 <- list(
    ## these params are part of the main parameter vector: go to run_sim()
    params=c(log_E0=4      ## initial exposed
           , log_beta0=-1  ## initial baseline transmission
             ## fraction of mild (non-hosp) cases
           , log_mu=log(params[["mu"]])
             ## fraction of incidence reported
             ## logit_c_prop=qlogis(params[["c_prop"]]),
             ## fraction of hosp to acute (non-ICU)
           , logit_phi1=qlogis(params[["phi1"]])
           , log_zeta=-1
             )
    ## NB dispersion
    , log_nb_disp=NULL
)

t_ont_cal_zeta1 <- system.time(ont_cal_zeta1 <-
                                   calibrate(data=ont_noICU
                                           , use_DEoptim=TRUE
                                           , DE_cores=4
                                           ## , debug_plot = TRUE
                                           ## , debug=TRUE
                                           , base_params=params
                                           , opt_pars = opt_pars_zeta1
                                           , time_args = NULL
                                             ## default is run_sim_break
                                        )
                                 ) ## system.time



