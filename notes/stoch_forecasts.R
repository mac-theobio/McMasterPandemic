library(McMasterPandemic)
library(cowplot)
load("../ontario/ontario_calibration_2brks_ndt.RData")
fit <- ont_cal_2brks_ndt  ## ont_cal_2brks
p0_proc <- predict(fit, 
                  stoch=c(proc=TRUE,obs=FALSE),
                  new_params=list(proc_disp=0))


f_args <- attr(fit, "forecast_args")
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
plot_grid(plotlist=purrr::map(pred_list,plot))
