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
dd <- r1 %>% select(date,report) %>% pivot() %>% na.omit()

dd_rh <- r1 %>% select(date,report,hosp) %>% pivot() %>% na.omit()

c1 <- calibrate_comb(params=p2, use_phenomhet=FALSE, debug_plot=TRUE, data=dd, use_DEoptim=FALSE,
                    use_spline=FALSE)

c2 <- calibrate_comb(params=p2, use_phenomhet=FALSE, debug_plot=TRUE, data=dd_rh, use_DEoptim=FALSE,
                    use_spline=FALSE)

