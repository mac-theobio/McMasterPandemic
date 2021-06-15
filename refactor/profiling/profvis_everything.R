
library(profvis)

########################################
## Profiling the base run_sim function

profvis({
library(McMasterPandemic)
library(ggplot2); theme_set(theme_bw())
library(cowplot)
library(tidyverse)
library(zoo)


params1 <- read_params("ICU1.csv")
knitr::kable(describe_params(params1))
knitr::kable(round(t(summary(params1)),2))
knitr::kable(round(t(get_R0(params1, components=TRUE)),2))


state1 <- make_state(params=params1)

########################################
## Profiling the base run_sim function

res1 <- run_sim(params=params1, state=state1, start_date=sdate, end_date=edate)


summary(res1)
plot_grid(plot(res1, log=TRUE), ## logarithmic
           plot(res1)) ## linear

############################################
## stochastic run_sim
set.seed(101)
params1obs <- update(params1, obs_disp=200)


res1obs <- run_sim(params1obs, state1, start_date=sdate, end_date=edate,
                   stoch=c(obs=TRUE, proc=FALSE))

summary(res1obs)
plot_grid(plot(res1obs, log=TRUE),
           plot(res1obs))

######################################################
## demographic stochastic run_sim

params1proc <- update(params1,E0=200,proc_disp=0) ## demog stoch only


res1proc <- run_sim(params1proc, start_date=sdate, end_date=edate,
                       stoch=c(obs=FALSE, proc=TRUE))


params1proc2 <- update(params1,E0=200, proc_disp=0.5, obs_disp=5)

res1proc2 <- run_sim(params1proc2, start_date=sdate, end_date=edate,
                      stoch=c(obs=FALSE, proc=TRUE))
plot_grid(plot(res1proc2, log=TRUE), plot(res1proc2))

########################################################################3
# Time-dependent
time_pars <- data.frame(Date=c("2020-Mar-10","2020-Mar-25"),
                         Symbol=c("beta0","beta0"),
                         Relative_value=c(0.5,0.1))



restimedep <- run_sim(params1,state1,start_date=sdate,end_date=edate,
                      params_timevar=time_pars,ndt=20, condense=FALSE)


summary(restimedep)
plot_grid(plot(restimedep, log=TRUE, condense=FALSE),
           plot(restimedep, condense=FALSE))

######################################################################
print(summary(params1))
newparams1 <- fix_pars(params1, target=c(R0=2))
print(summary(newparams1))



######################################################################
report_data <- (res1obs %>%
      mutate(value=round(report), var="report") %>%
      select(date, value, var) %>%
      na.omit() )
head(report_data)

report_death_data <- (res1obs %>%
    select(date, report, death) %>%
    pivot_longer(names_to = "var", -date) %>%
    mutate(value=round(value)) %>%
    na.omit())
head(report_death_data, n=12)


opt_pars <- list(params = c(beta0=0.1))


fitted.mod <- calibrate(
   data = report_death_data,
   start_date = sdate
   ## skip breaks that are present by default:
   , time_args = list(break_dates = NULL)
   , base_params = params1obs
   , opt_pars = opt_pars
   ##, debug_plot = TRUE # instructive plotting during optimization
     )


plot(fitted.mod, data=report_death_data)
plot(fitted.mod, data=report_death_data,
      predict_args=list(keep_vars=c("report","death")))

######################################################################
})
