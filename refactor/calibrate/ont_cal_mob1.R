library(McMasterPandemic)
library(dplyr)
## comb_sub comes in MacPan package
params <- fix_pars(read_params("ICU1.csv"))

## select the part of the mobility data that lies within the calibration data set
comb_sub2 <- (comb_sub
    %>% right_join(ont_all %>% select(date) %>% unique(),by="date")
    %>% na.omit()
    %>% mutate_at("rel_activity",pmin,1) ## cap relative mobility at 1
)

plot(rel_activity ~ date, comb_sub)
with(comb_sub2, lines(date, rel_activity, col=2))



opt_pars <- list(params=c(log_E0=4
                        , log_beta0=-1
                          ## fraction of mild cases (-> report vs hosp)
                        , log_mu=log(params[["mu"]])
                          ## fraction to ICU (hosp vs death)
                        , logit_phi1=qlogis(params[["phi1"]])),
                 logit_mob_power=0.5,  ## include power parameter for relative mobility
                 log_nb_disp=NULL)

dd <- (ont_all %>% trans_state_vars() %>% filter(var %in% c("report", "death", "H")))

set.seed(101)
ont_cal_mob1 <- calibrate(data=dd, base_params=params, opt_pars=opt_pars,
                          ## debug_plot=TRUE,
                          ## debug=TRUE,
                      use_DEoptim=TRUE,
                      time_args = list(mob_value=comb_sub2$rel_activity,
                                       mob_startdate=comb_sub2$date[1]),
                      sim_fun=run_sim_mobility)

save("ont_cal_mob1", file=sprintf("data/ONcalib_mob1.rda",format(Sys.time(),"%Y%b%d")))

