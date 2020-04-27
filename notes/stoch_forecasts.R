library(McMasterPandemic)
library(cowplot)
load("../ontario/ontario_calibration_2brks.RData")
load("../ontario/ontario_calibration.RData")
load("../ontario/ontario_calibration_hosponly.RData")
fit <- ont_cal_2brks  ## ont_cal_2brks
fit <- ont_cal1
fit <- ont_cal2  ## hosp only
p0 <- predict(fit)


eigen(vcov(fit))
np <- nrow(vcov(fit))
eigen(vcov(fit))$vec[,np]
eigen(solve(fit$hessian[-np,-np]))$val
heatmap(fit$hessian,Rowv=NA,Colv=NA)
fit$par
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
plot_grid(plotlist=purrr::map(pred_list,plot))
