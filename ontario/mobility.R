library(McMasterPandemic)
library(dplyr)
load("ontario_clean.RData")
## comb_sub
## restrict to dates with data
comb_sub2 <- (comb_sub
    %>% right_join(ont_all %>% select(date),by="date")
    %>% na.omit()
    %>% mutate_at("rel_activity",pmin,1)
)

plot(rel_activity ~ date, comb_sub2)

opt_pars <- list(params=c(log_E0=4
                        , log_beta0=-1
                        , log_mu=log(params[["mu"]])
                        , logit_phi1=qlogis(params[["phi1"]])),
                 mob_value=comb_sub2$rel_activity,
                 mob_startdate=comb_sub2$date[1],
                 logit_mob_power=0.5,
                 log_nb_disp=NULL)

dd <- (ont_all %>% trans_state_vars() %>% filter(var %in% c("report", "death", "H")))

debug(calibrate)

calibrate(data=dd, base_params=params, opt_pars=opt_pars, debug_plot=TRUE,
          time_args = list(mob_value=comb_sub2$rel_activity,
                           mob_startdate=comb_sub2$date[1]),
          sim_fun=run_sim_mobility)
