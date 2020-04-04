library("McMasterPandemic")
library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())

load(system.file("testdata","calib_test.RData",package="McMasterPandemic"))
### 
summary(cparams)
cparams[["obs_disp"]] <- 200

set.seed(101)
sim1S <- run_sim(cparams, cstate, start_date="1-Mar-2020",
                 end_date="31-Mar-2020",
                 stoch=c(obs = TRUE, proc = FALSE))
plot(sim1S,log=TRUE)

## aggregate/subset simulated data to a short time window (15 Mar - 29 Mar)/
simdat <- (aggregate(sim1S,pivot=TRUE)
    %>% filter(date>as.Date("2020-03-15") & date<as.Date("2020-03-29"))
)
## reasonable trajectories
## note H > ICU > D as it should be
plot(sim1S,log=TRUE) + geom_point(data=simdat)

## DRY: condense this some more

## extract H data and set t==0 at beginning
regdatS <- (simdat
  %>% mutate(t0=as.numeric(date-min(date)))
  %>% filter(var=="H")
)
##  fit to H data
g1S <- MASS::glm.nb(value~t0,data=regdatS)

## calibrate (use orig, sim params as starting point)
ccS <- calibrate(coef(g1S)[1],coef(g1S)[2],
                pop=cparams[["N"]],
                params=cparams,
                date0="1-Mar-2020",
                date1=min(regdatS$date))
summary(ccS$params)
## run simulation with calibrated data
simScal <- run_sim(ccS$params, ccS$state, start_date="1-Mar-2020",
                   end_date="1-Apr-2020")
## step_args=list(do_hazard=TRUE))

pframeS <- data.frame(date=seq(as.Date("2020-03-01"),
                          as.Date("2020-03-30"),
                          by="1 day"),
                 t0=seq(-15,14))
pframeS <- (pframeS
    %>% mutate(value=predict(g1S,newdata=pframeS,type="response"),
               var="H")
)


print(plot(simScal,log=TRUE)
      + geom_point(data=simdat)
      + geom_line(data=aggregate(sim1S,pivot=TRUE))
      + geom_line(data=pframeS,lty=2)
      )
##  dashed line is NB regression fit
## conclusion: slope of simulation is not giving the expected growth rate.
## problem with make_jac?

