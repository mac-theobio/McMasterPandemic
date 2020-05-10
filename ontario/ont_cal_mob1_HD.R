library(McMasterPandemic)
library(dplyr)

## update opt_pars; we no longer use reports, so don't want to try to calibrate reports vs hosp
opt_pars$log_mu <- NULL
dd <- dd %>% filter(var %in% c("death", "H"))

set.seed(101)
ont_cal_mob1_HD <- calibrate(data=dd, base_params=params, opt_pars=opt_pars,
                         use_DEoptim=TRUE,
                         time_args = list(mob_value=comb_sub2$rel_activity,
                                          mob_startdate=comb_sub2$date[1]),
                         sim_fun=run_sim_mobility)
