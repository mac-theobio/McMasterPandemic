library(styler)
library(profvis)
library(McMasterPandemic)
library(anytime)
library(ggplot2)
theme_set(theme_bw())
library(cowplot)
library(tidyverse)
library(zoo)
# Required to make profvis tell us which problematic lines exist in the original run_sim function


source("../../R/sim_funs.R")
source("../../R/calibrate.R")
source("../../R/afuns.R")
source("../../R/utils.R")
source("../../R/testify.R")
source("../../R/methods.R")
source("../../R/kfuns.R")



library(bbmle)

params1 <- read_params("ICU1.csv")
state1 <- make_state(params = params1)
sdate <- "2020-Feb-10"
edate <- "2020-Jun-1"

#################################################################################################
# Non-stochastic simulation

classic <- profvis(
  {
    for (i in c(1:100)) {
      res1 <- run_sim(params = params1, state = state1, start_date = sdate, end_date = edate)
    }
  },
  prof_output = "classic.out"
)
classic_summary <- summaryRprof("classic.out")

# Note: in afuns, expsim is being used by default
write.csv(classic_summary$`by.self`, 'classic.csv')
htmlwidgets::saveWidget(classic, "classic_profile.html")

################################################################################################
# Demographic stochastic simulation
params1proc2 <- update(params1, E0 = 200, proc_disp = 0.5, obs_disp = 5)

demographic_stochastic <- profvis(
  {
    for (i in c(1:100)) {
      res1proc2 <- run_sim(params1proc2,
        start_date = sdate, end_date = edate,
        stoch = c(obs = FALSE, proc = TRUE)
      )
    }
  },
  prof_output = "demographic_stochastic.out"
)


summary = summaryRprof("demographic_stochastic.out")
write.csv(summary$`by.self`, 'demographic_stochastic.csv')
htmlwidgets::saveWidget(demographic_stochastic, "demographic_stochastic_profile.html")

#########################################################
# Profile calibrate

set.seed(101)
params1obs <- update(params1, obs_disp = 200)
res1obs <- run_sim(params1obs, state1,
  start_date = sdate, end_date = edate,
  stoch = c(obs = TRUE, proc = FALSE)
)

report_data <- (res1obs %>%
  mutate(value = round(report), var = "report") %>%
  select(date, value, var) %>%
  na.omit())
head(report_data)

report_death_data <- (res1obs %>%
  select(date, report, death) %>%
  pivot_longer(names_to = "var", -date) %>%
  mutate(value = round(value)) %>%
  na.omit())
head(report_death_data, n = 12)


opt_pars <- list(params = c(beta0 = 0.1))


calibrated <- profvis(
  {
    for (i in c(1:10)) {
      fitted.mod <- calibrate(
        data = report_death_data,
        start_date = sdate
        ## skip breaks that are present by default:
        , time_args = list(break_dates = NULL),
        base_params = params1obs,
        opt_pars = opt_pars
        ## , debug_plot = TRUE # instructive plotting during optimization
      )
    }
  },
  prof_output = "calibrate.out"
)


summary = summaryRprof("calibrate.out")
write.csv(summary$`by.self`, 'calibrate.csv')
htmlwidgets::saveWidget(calibrated, "calibrate.html")
