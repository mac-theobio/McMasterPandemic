library(McMasterPandemic)
## devtools::load_all("../../McMasterPandemic")
library(tidyverse)
library(furrr)

p1 <- fix_pars(read_params("ICU1.csv"))
summary(p1)
p2 <- update(p1, obs_disp=1, proc_disp=0, zeta=5)
set.seed(101)
r1 <- run_sim(p2, stoch=c(obs=TRUE, proc=TRUE), end_date="2020-05-31")
plot(r1, log=TRUE)
dd_r <- r1 %>% select(date,report) %>% pivot() %>% na.omit()
dd_rh <- r1 %>% select(date,report,hosp) %>% pivot() %>% na.omit()


c_rh <- calibrate_comb(params=p2, use_phenomhet=FALSE, debug_plot=TRUE,
                       data=dd_rh, use_DEoptim=FALSE,
                       use_spline=FALSE)
## re-run calibrate() [NOT calibrate_comb]
c_r <- update(c_rh, data=dd_r)
## re-run calibrate_comb (only change is data)
c_r2 <- calibrate_comb(params=p2, use_phenomhet=FALSE, debug_plot=TRUE,
                       data=dd_r, use_DEoptim=FALSE,
                       use_spline=FALSE)

## inspect components
t1 <- c_rh$forecast_args$time_args
t2 <- c_r2$forecast_args$time_args
cmp <- map2(t1, t2, ~all.equal(.x,.y))
cmp[!map_lgl(cmp,isTRUE)]
## start_date, time_args, opt_pars

range(t1$X_date)  ## starts 03-20
range(t2$X_date)  ## starts 03-29







