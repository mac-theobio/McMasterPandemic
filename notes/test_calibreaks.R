library(McMasterPandemic)
library(ggplot2)
library(tidyverse)
library(anytime)

params <- fix_pars(read_params("ICU1.csv"))

start_date <- anydate("2020-01-01")
end_date <- anydate("2020-04-01")
break1 <- "2020-02-22"
bd <- anydate(c(break1))
rel_break1 <- 0.2

params[["N"]] <- 1e7

sim1break <- run_sim_break(params
   , start_date=start_date
   , end_date=end_date
   , break_dates = bd
   , rel_beta0= rel_break1
     )
plot(sim1break)

simdat <- pivot(condense(sim1break))
simdat <- simdat %>% filter(var %in% c("report","I","incidence"))
print(ggsim <- ggplot(simdat,aes(x=date,y=value,color=var))
   + geom_line()
   + geom_vline(xintercept = bd)
   + scale_y_log10()
   + geom_point()
)

## What is wrong with the break date? Why isn't it turning at the break?

## Need to round value because of negative binomial fit
dd <- (simdat 
   %>% filter(var == "report") 
   %>% filter(!is.na(value)) 
   %>% mutate(value = round(value))
)

opt_pars <- list(
   ## these params go to run_sim
   params=c(log_E0=4, log_beta0=-1)
   , log_rel_beta0 = rep(-1, length(bd))
   , log_nb_disp=0
)

g1 <- calibrate(data=dd, base_params=params
   , opt_pars = opt_pars
   , break_dates = bd
)

pp <- invlink_trans(restore(g1$mle2@coef,opt_pars))

print(pp)

print(rel_break1)

print(plot(g1))
print(ggsim)
