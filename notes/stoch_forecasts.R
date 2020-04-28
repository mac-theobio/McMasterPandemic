library(McMasterPandemic)
library(cowplot)
load("../ontario/ontario_calibration_noICU_2brks_prior.RData")
load("../ontario/ontario_calibration.RData")
fit <- ont_cal_noICU_2brks_prior  ## ont_cal_2brks

p0 <- predict(fit)


f_args <- fit$forecast_args
p1 <- predict(fit, ensemble=TRUE)
p1_obs <- predict(fit, ensemble=TRUE,
                  stoch=c(proc=FALSE,obs=TRUE),
                  new_params=list(obs_disp=20))
p1_proc <- predict(fit, ensemble=TRUE,
                  stoch=c(proc=TRUE,obs=FALSE),
                  new_params=list(proc_disp=0))
p1_obs_proc <- predict(fit, ensemble=TRUE,
                  stoch=c(proc=TRUE,obs=TRUE),
                  new_params=list(obs_disp=20,proc_disp=0.2))
pred_list <- list(p1,p1_obs,p1_proc,p1_obs_proc)
cowplot::plot_grid(plotlist=purrr::map(pred_list,plot))
