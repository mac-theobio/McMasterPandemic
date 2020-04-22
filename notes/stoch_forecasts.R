library(McMasterPandemic)
library(cowplot)
load("../ontario/ontario_calibration_2brks.RData")
load("../ontario/ontario_calibration.RData")
load("../ontario/ontario_calibration_hosponly.RData")
fit <- ont_cal_2brks  ## ont_cal_2brks
fit <- ont_cal1
fit <- ont_cal2  ## hosp only
p0 <- predict(fit)

## pp <- coef(fit)
## s0 <- make_state(params=pp)
## M <- make_ratemat(params=pp,state=s0)

## myEM <- function(n, size, rate, dt) {
##     S  <- sum(rate)
##     p0 <- (1-exp(-S*dt))
##     p  <- c(1-p0,p0*rate/S)
##     return(rmultinom(n,size=size, prob=p)[-1,,drop=FALSE])
## }


## ## checking hand-rolled version vs pomp version
## set.seed(101)
## pomp_res <- pomp::reulermultinom(100,size=10000,rate=c(a=1,b=2,c=3),dt=0.1)
## set.seed(101)
## my_res <- myEM(100,size=10000,rate=c(a=1,b=2,c=3),dt=0.1)
## all.equal(pomp_res,my_res)

## dump(c("s0","M"),file="dumpdata.txt")
## source("dumpdata.txt")
## i <- 1

## floor(s0[[1]]+0.5)==s0[[1]]
## (s0[[1]]+0.5) - s0[[1]]
## floor(s0[[1]]+0.5) - s0[[1]]

## pomp::reulermultinom(n=1,size=s0[[1]],rate=M[i,-i],dt=1)
## myEM(n=1,size=s0[[1]],rate=M[i,-i],dt=1)

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
