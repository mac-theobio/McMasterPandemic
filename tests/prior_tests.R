library(McMasterPandemic)
library(dplyr)

## sim example
params <- read_params("ICU1.csv")
paramsS <- update(params,c(proc_disp=0.1,obs_disp=100))
state <- make_state(params=params)
sdate <- "2020-02-10" ## arbitrary!
set.seed(101)
res1 <- run_sim(params,state,start_date=sdate,end_date="2020-06-01")
res1_S <- update(res1, params=paramsS, stoch=c(obs=TRUE, proc=TRUE))

cdat <- (res1_S
    %>% pivot()
    %>% filter(var=="report",
               date>as.Date("2020-03-01"),
               date<as.Date("2020-04-15"))
)

priors <- list(~dlnorm(rel_beta0[1],meanlog=-1,sd=0.5))
c0 <- calibrate(data=cdat,base_params=params) ## ,debug_plot=TRUE,debug=TRUE)
c1 <- calibrate(data=cdat,base_params=params,priors=priors) ## debug_plot=TRUE)
summary(c0)
summary(c1)
                
