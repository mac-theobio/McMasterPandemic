library("McMasterPandemic")
library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())

L <- load(system.file("testdata","calib_test.RData",package="McMasterPandemic"))
summary(cparams)
cparams[["obs_disp"]] <- 20

## FIXME: thinning interacts with ndt?
set.seed(101)
sim1S <- run_sim(cparams, cstate, start_date="1-Mar-2020",
                 end_date="31-Mar-2020",
                 ndt=10,
                 stoch=c(obs = TRUE, proc = FALSE))

## aggregate/subset simulated data to a short time window (15 Mar - 29 Mar)/
simdat <- (aggregate(sim1S,pivot=TRUE)
    %>% filter(date>as.Date("2020-03-15") & date<as.Date("2020-03-29"))
)

simdat %>% filter(date==as.Date("2020-03-17"),var=="D")
sim1S %>% filter(date==as.Date("2020-03-17")) %>% select(D)
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

## brute-force (shooting) calibration
ccS2 <- calibrate(coef(g1S)[1],coef(g1S)[2],
                pop=cparams[["N"]],
                params=cparams,
                date0="1-Mar-2020",
                date1=min(regdatS$date),
                init_target=head(regdatS$value,1),
                sim_args=list(ndt=10))

summary(ccS$params)
## run (deterministic) simulation with calibrated data
simScal <- run_sim(ccS$params, ccS$state,
                   start_date="1-Mar-2020",
                   end_date="1-Apr-2020",
                   ndt=10)

L <- load("tmp.RData")
ccS2$params[["E0"]] <- ccS2$state[["E"]]
all.equal(params,ccS2$params)
all.equal(state,ccS2$state)

## try brute-force-calibrated initial conditions
simScal_brute <- run_sim(ccS2$params, ccS2$state,
                   start_date="1-Mar-2020",
                   end_date="1-Apr-2020",
                   ndt=10)
simScal_brute %>% filter(date==date1)  ## 45 + 1(H2),not 36???
## step_args=list(do_hazard=TRUE))

## predicted values from regression
pframeS <- data.frame(date=seq(as.Date("2020-03-01"),
                               as.Date("2020-03-30"),
                               by="1 day"),
                      t0=seq(-15,14))
pframeS <- (pframeS
    %>% mutate(value=predict(g1S,newdata=pframeS,type="response"),
               var="H")
)

print(gg1 <- plot(simScal,log=TRUE)
      + geom_point(data=simdat)
      + geom_line(data=aggregate(sim1S,pivot=TRUE),lty=3)
      + geom_line(data=pframeS,lty=2)
      + geom_line(data=aggregate(simScal_brute,pivot=TRUE))
      )
print(gg1
      + geom_hline(yintercept=36,lty=2)
      + geom_vline(xintercept=as.Date("2020-03-16"),lty=2)
      )
##  dashed line is NB regression fit
## conclusion: slope of simulation is not giving the expected growth rate.
## problem with make_jac?

## what is the actual r?
simAgg <- aggregate(simScal)[,c("H","ICU","D")]
n <- nrow(simAgg)
print(log(unlist(simAgg[n,]/simAgg[n-10,]))/10)
##         H       ICU         D 
## 0.1942182 0.1946885 0.2026567 
summary(ccS$params)

## LESSON: *small* details in 


